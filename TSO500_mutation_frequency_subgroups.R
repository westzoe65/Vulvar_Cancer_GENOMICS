###############################################
# TSO500 ALTERATION FREQUENCIES BY STAGE AND MOLECULAR SUBTYPE
###############################################

## ============================================================
## 0. Libraries
## ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(writexl)
})

## ============================================================
## 0b. User toggles
## ============================================================
INCLUDE_SYNONYMOUS <- FALSE

RESTRICT_TO_STATUS_SAMPLES <- TRUE

RUN_FISHER <- FALSE

STAGE_GROUPS <- c("Stage I", "Stage III")

SUBTYPE_LEVELS <- c(
  "HPV-associated",
  "HPV-independent/p53-wildtype",
  "HPV-independent/p53-mutant",
  "Unknown"
)

SUBTYPE_TEST_EXCLUDES_UNKNOWN <- TRUE

## ============================================================
## 1. Paths
## ============================================================
dir_small    <- "<SMALL_VARIANT_DIR>"
dir_combined <- "<COMBINED_VARIANT_DIR>"
dir_out      <- "<OUTPUT_DIR>"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

status_file  <- "<METASTATIC_STATUS_FILE>"
file_stage   <- "<STAGE_FILE>"
file_subtype <- "<MOLECULAR_SUBTYPE_FILE>"

for (.f in c(file_stage, file_subtype,
             if (RESTRICT_TO_STATUS_SAMPLES) status_file)) {
  if (!file.exists(.f)) stop("File does not exist: ", .f)
}
rm(.f)

## ============================================================
## 2. Helper functions
## ============================================================
normalize_sample_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^\\d+$", x), paste0("PQ", as.integer(x)), x)
}

empty_variants <- function() {
  tibble(sample = character(), gene = character(),
         consequence_term = character(), alteration_class = character())
}

read_event_annotation <- function(file, id_col_candidates, value_col_candidates) {
  df <- readxl::read_xlsx(file, .name_repair = "minimal") %>%
    janitor::clean_names()
  
  id_clean  <- janitor::make_clean_names(id_col_candidates)
  val_clean <- janitor::make_clean_names(value_col_candidates)
  id_col    <- intersect(id_clean,  names(df))[1]
  val_col   <- intersect(val_clean, names(df))[1]
  
  if (is.na(id_col) || is.na(val_col)) {
    stop("Could not find required columns in: ", file, "\n",
         "Tried id columns: ",    paste(id_clean,  collapse = ", "), "\n",
         "Tried value columns: ", paste(val_clean, collapse = ", "), "\n",
         "Columns present: ",     paste(names(df), collapse = ", "))
  }
  
  df %>%
    transmute(sample = normalize_sample_id(.data[[id_col]]),
              value  = na_if(str_squish(as.character(.data[[val_col]])), "")) %>%
    filter(!is.na(sample), sample != "") %>%
    arrange(is.na(value)) %>%
    distinct(sample, .keep_all = TRUE)
}

snpeff_to_class <- function(x) {
  x <- tolower(as.character(x))
  case_when(
    str_detect(x, "frameshift")           ~ "Frameshift",
    str_detect(x, "stop_gained|nonsense") ~ "Nonsense",
    str_detect(x, "intron")               ~ "Intron",
    str_detect(x, "splice")               ~ "Splice",
    str_detect(x, "missense")             ~ "Missense",
    str_detect(x, "inframe")              ~ "Inframe",
    str_detect(x, "synonymous")           ~ "Synonymous",
    TRUE                                  ~ "Other"
  )
}

wilson_ci <- function(x, n, conf_level = 0.95) {
  z     <- qnorm(1 - (1 - conf_level) / 2)
  lower <- upper <- rep(NA_real_, length(x))
  ok    <- !is.na(n) & n > 0
  p     <- x[ok] / n[ok]
  nn    <- n[ok]
  denom  <- 1 + z^2 / nn
  centre <- p + z^2 / (2 * nn)
  margin <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2))
  lower[ok] <- pmax(0, (centre - margin) / denom)
  upper[ok] <- pmin(1, (centre + margin) / denom)
  list(lower = lower, upper = upper)
}

fmt_pct_ci <- function(pct, lo, hi) {
  ifelse(is.na(pct), NA_character_,
         sprintf("%.1f%% (%.1f\u2013%.1f)", pct, lo, hi))
}

recode_stage_group <- function(x) {
  s <- toupper(as.character(x))
  s <- str_remove_all(s, "STAGE|FIGO|\\s")
  case_when(
    is.na(s) | s == ""          ~ "Unknown",
    str_detect(s, "^(IV|4)")    ~ "Stage IV",
    str_detect(s, "^(III|3)")   ~ "Stage III",
    str_detect(s, "^(II|2)")    ~ "Stage II",
    str_detect(s, "^(I|1)")     ~ "Stage I",
    TRUE                        ~ "Unknown"
  )
}

recode_subtype <- function(x) {
  s <- tolower(as.character(x))
  s <- str_remove_all(s, "[^a-z0-9]")
  indep <- str_detect(s, "independent|hpvi|hpvneg")
  out <- case_when(
    is.na(s) | s == "" | str_detect(s, "unknown|notrecorded|^na$") ~ "Unknown",
    str_detect(s, "hpvassoc|hpvpos|^hpva$")                         ~ "HPV-associated",
    indep & str_detect(s, "abnormal|mutant|mut|aberrant")           ~ "HPV-independent/p53-mutant",
    indep & str_detect(s, "wild|wt|normal")                         ~ "HPV-independent/p53-wildtype",
    TRUE                                                            ~ "Unknown"
  )
  unname(out)
}

## ============================================================
## 3. Discover samples
## ============================================================
small_files <- list.files(dir_small,
                          pattern    = "^PQ[0-9]+_filtered_[0-9]{8}\\.xlsx$",
                          full.names = TRUE)
if (length(small_files) == 0) stop("No small variant xlsx files found in: ", dir_small)

file_map <- tibble(
  sample     = sub("_filtered_[0-9]{8}\\.xlsx$", "", basename(small_files)),
  small_file = small_files
) %>%
  mutate(combined_file = file.path(dir_combined,
                                   paste0(sample, "_CombinedVariantOutput_filtered.txt")),
         sample_num    = as.integer(str_remove(sample, "^PQ"))) %>%
  arrange(sample_num)

if (RESTRICT_TO_STATUS_SAMPLES) {
  status_raw <- readxl::read_xlsx(status_file, .name_repair = "minimal") %>%
    janitor::clean_names()
  s_col <- intersect(names(status_raw),
                     c("sample", "randomisation_no", "randomization_no", "id", "patient"))[1]
  m_col <- intersect(names(status_raw),
                     c("status", "metastatic_status", "met_status", "metastasis", "group"))[1]
  if (is.na(s_col) || is.na(m_col)) {
    stop("Could not find sample/status columns in: ", status_file,
         "\nColumns present: ", paste(names(status_raw), collapse = ", "))
  }
  status_samples <- status_raw %>%
    transmute(sample = normalize_sample_id(.data[[s_col]]),
              status = toupper(str_remove_all(as.character(.data[[m_col]]), "[^A-Za-z]"))) %>%
    filter(status %in% c("METASTATIC", "NONMETASTATIC")) %>%
    pull(sample)
  file_map <- file_map %>% filter(sample %in% status_samples)
}

cohort <- file_map$sample
if (length(cohort) == 0) stop("No samples left in the cohort.")
message("Cohort size: ", length(cohort), " samples")

## ============================================================
## 4. Read variants (small variants, CNVs, large rearrangements)
## ============================================================
read_small_variants <- function(file, sample_id) {
  dat <- readxl::read_xlsx(file, .name_repair = "minimal") %>%
    janitor::clean_names()
  gene_col  <- grep("^gene",                             names(dat), value = TRUE)[1]
  coseq_col <- grep("^coseq$|^consequence$|consequence", names(dat), value = TRUE)[1]
  if (is.na(gene_col) || is.na(coseq_col)) return(empty_variants())
  
  tibble(sample = sample_id,
         gene   = toupper(as.character(dat[[gene_col]])),
         coseq  = as.character(dat[[coseq_col]])) %>%
    filter(!is.na(gene), gene != "", !is.na(coseq), coseq != "") %>%
    mutate(consequence_term = str_trim(tolower(coseq)),
           alteration_class = snpeff_to_class(consequence_term)) %>%
    select(sample, gene, consequence_term, alteration_class)
}

extract_section <- function(raw, header_regex) {
  idx <- grep(header_regex, raw)
  if (length(idx) == 0) return(NULL)
  all_idx <- grep("^\\[.*\\]", raw)
  nxt     <- all_idx[all_idx > idx][1]
  end_idx <- ifelse(is.na(nxt), length(raw), nxt - 1)
  lines   <- raw[(idx + 1):end_idx]
  lines   <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 2) return(NULL)
  tryCatch(
    readr::read_tsv(I(lines), show_col_types = FALSE, name_repair = "minimal") %>%
      janitor::clean_names(),
    error = function(e) NULL
  )
}

read_cnv_variants <- function(file, sample_id) {
  if (!file.exists(file)) {
    warning("Combined variant file missing for ", sample_id, ": ", file)
    return(empty_variants())
  }
  cnv_df <- extract_section(readLines(file, warn = FALSE), "^\\[Copy Number Variants\\]")
  if (is.null(cnv_df) || nrow(cnv_df) == 0) return(empty_variants())
  
  gene_col <- grep("^gene",               names(cnv_df), value = TRUE)[1]
  cnv_col  <- grep("copy_number_variant", names(cnv_df), value = TRUE)[1]
  if (is.na(gene_col) || is.na(cnv_col)) return(empty_variants())
  
  cnv_vec <- as.character(cnv_df[[cnv_col]])
  tibble(sample           = sample_id,
         gene             = toupper(as.character(cnv_df[[gene_col]])),
         consequence_term = NA_character_,
         alteration_class = case_when(
           str_detect(cnv_vec, regex("dup", TRUE)) ~ "Duplication",
           str_detect(cnv_vec, regex("del", TRUE)) ~ "Deletion",
           TRUE ~ NA_character_)) %>%
    filter(!is.na(alteration_class), !is.na(gene), gene != "")
}

read_lr_variants <- function(file, sample_id) {
  if (!file.exists(file)) return(empty_variants())
  lr_df <- extract_section(readLines(file, warn = FALSE), "^\\[Large Rearrangements\\]")
  if (is.null(lr_df) || nrow(lr_df) == 0) return(empty_variants())
  
  gene_col    <- grep("^gene",    names(lr_df), value = TRUE)[1]
  cnvtype_col <- grep("cnv_type", names(lr_df), value = TRUE)[1]
  if (is.na(gene_col) || is.na(cnvtype_col)) return(empty_variants())
  
  cnvtype_vec <- toupper(trimws(as.character(lr_df[[cnvtype_col]])))
  tibble(sample           = sample_id,
         gene             = toupper(trimws(as.character(lr_df[[gene_col]]))),
         consequence_term = NA_character_,
         alteration_class = case_when(
           str_detect(cnvtype_vec, "GAIN") ~ "Duplication",
           str_detect(cnvtype_vec, "LOSS") ~ "Deletion",
           TRUE ~ NA_character_)) %>%
    filter(!is.na(alteration_class), !is.na(gene), gene != "")
}

all_variants <- bind_rows(
  map2_dfr(file_map$small_file,    file_map$sample, read_small_variants),
  map2_dfr(file_map$combined_file, file_map$sample, read_cnv_variants),
  map2_dfr(file_map$combined_file, file_map$sample, read_lr_variants)
) %>%
  mutate(sample = normalize_sample_id(sample),
         gene   = toupper(str_squish(gene))) %>%
  filter(sample %in% cohort, !is.na(gene), gene != "")

if (!INCLUDE_SYNONYMOUS) {
  all_variants <- all_variants %>% filter(alteration_class != "Synonymous")
}

## ============================================================
## 5. Group assignments
## ============================================================
id_candidates <- c("Randomisation No.", "Randomization No.", "Sample", "ID", "Patient")

stage_samples <- tibble(sample = cohort) %>%
  left_join(read_event_annotation(file_stage, id_candidates, c("Stage")),
            by = "sample") %>%
  rename(stage_raw = value) %>%
  mutate(stage_group = recode_stage_group(stage_raw))

message("\nStage recoding (check this looks right):")
print(stage_samples %>% count(stage_raw, stage_group), n = Inf)

subtype_samples <- tibble(sample = cohort) %>%
  left_join(read_event_annotation(file_subtype, id_candidates,
                                  c("Molecular subtype", "Molecular_subtype",
                                    "Subtype", "IHC molecular subtype")),
            by = "sample") %>%
  rename(subtype_raw = value) %>%
  mutate(subtype_group = recode_subtype(subtype_raw))

message("\nMolecular subtype recoding (check this looks right):")
print(subtype_samples %>% count(subtype_raw, subtype_group), n = Inf)

## ============================================================
## 6. Frequency engine
## ============================================================
group_frequency_long <- function(variants, groups_df, group_levels) {
  groups_df <- groups_df %>% filter(group %in% group_levels)
  
  sizes <- tibble(group = group_levels) %>%
    left_join(count(groups_df, group, name = "N"), by = "group") %>%
    mutate(N = replace_na(N, 0L))
  
  categories <- list(
    Any_alteration = variants,
    SNV_indel      = filter(variants, !alteration_class %in% c("Duplication", "Deletion")),
    Duplication    = filter(variants, alteration_class == "Duplication"),
    Deletion       = filter(variants, alteration_class == "Deletion")
  )
  
  long <- imap(categories, function(v, cat_name) {
    counts <- v %>%
      distinct(sample, gene) %>%
      inner_join(groups_df, by = "sample") %>%
      count(gene, group, name = "n")
    if (nrow(counts) == 0) return(NULL)
    
    out <- expand_grid(gene = sort(unique(counts$gene)), group = group_levels) %>%
      left_join(counts, by = c("gene", "group")) %>%
      left_join(sizes,  by = "group") %>%
      mutate(n = replace_na(n, 0L))
    
    ci <- wilson_ci(out$n, out$N)
    out %>%
      mutate(category = cat_name,
             pct      = ifelse(N > 0, 100 * n / N, NA_real_),
             ci_lo    = 100 * ci$lower,
             ci_hi    = 100 * ci$upper,
             summary  = fmt_pct_ci(pct, ci_lo, ci_hi)) %>%
      select(category, gene, group, n, N, pct, ci_lo, ci_hi, summary)
  })
  
  list(long = compact(long), sizes = sizes)
}

fisher_by_gene <- function(long_cat, test_levels) {
  long_cat %>%
    filter(group %in% test_levels, N > 0) %>%
    group_by(gene) %>%
    summarise(
      `Fisher p` = {
        m <- rbind(n, N - n)
        if (ncol(m) < 2 || sum(m[1, ]) == 0) NA_real_ else
          tryCatch(fisher.test(m, workspace = 2e7)$p.value,
                   error = function(e)
                     fisher.test(m, simulate.p.value = TRUE, B = 1e5)$p.value)
      },
      .groups = "drop"
    ) %>%
    mutate(`Fisher p (BH-adjusted)` = p.adjust(`Fisher p`, method = "BH"))
}

make_wide <- function(long_cat) {
  totals <- long_cat %>%
    group_by(gene) %>%
    summarise(`Total altered (n)` = sum(n), .groups = "drop")
  
  long_cat %>%
    mutate(label                = paste0(group, " (N=", N, ")"),
           `%`                  = round(pct,   1),
           `95% CI low`         = round(ci_lo, 1),
           `95% CI high`        = round(ci_hi, 1),
           `% (95% CI)`         = summary) %>%
    select(gene, label, n, `%`, `95% CI low`, `95% CI high`, `% (95% CI)`) %>%
    pivot_wider(names_from  = label,
                values_from = c(n, `%`, `95% CI low`, `95% CI high`, `% (95% CI)`),
                names_glue  = "{label} {.value}",
                names_vary  = "slowest") %>%
    left_join(totals, by = "gene") %>%
    relocate(`Total altered (n)`, .after = gene) %>%
    rename(Gene = gene) %>%
    arrange(desc(`Total altered (n)`), Gene)
}

run_group_analysis <- function(variants, groups_df, group_levels, test_levels,
                               sample_sheet, group_size_sheet, notes, out_file) {
  res    <- group_frequency_long(variants, groups_df, group_levels)
  sheets <- list(Notes = notes, Group_sizes = group_size_sheet)
  
  for (cat_name in names(res$long)) {
    long_cat <- res$long[[cat_name]]
    wide     <- make_wide(long_cat)
    if (RUN_FISHER && length(test_levels) >= 2) {
      wide <- wide %>%
        left_join(fisher_by_gene(long_cat, test_levels) %>% rename(Gene = gene),
                  by = "Gene")
    }
    sheets[[cat_name]] <- wide
  }
  
  sheets$Long_format <- bind_rows(res$long) %>%
    mutate(across(c(pct, ci_lo, ci_hi), ~ round(.x, 1))) %>%
    rename(Category = category, Gene = gene, Group = group,
           `%` = pct, `95% CI low` = ci_lo, `95% CI high` = ci_hi,
           `% (95% CI)` = summary)
  sheets$Samples <- sample_sheet
  
  write_xlsx(sheets, out_file)
  message("Saved: ", out_file)
  invisible(res)
}

common_notes <- c(
  "Values are the number and percentage of samples in each group with >=1 alteration in the gene.",
  "Denominator (N) = number of cohort samples in that group, including samples with no alterations.",
  "95% CI = Wilson score interval.",
  paste0("Synonymous variants ", if (INCLUDE_SYNONYMOUS) "INCLUDED." else "EXCLUDED."),
  "SNV_indel = missense, inframe, splice, frameshift, nonsense, intron, other.",
  "Duplication/Deletion include CNV calls plus large-rearrangement gains/losses (BRCA1/2).",
  paste0("Cohort: ", length(cohort), " samples",
         if (RESTRICT_TO_STATUS_SAMPLES) " (restricted to Metastatic_status.xlsx)." else ".")
)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

## ============================================================
## 7. Output 1: Stage I vs Stage III
## ============================================================
stage_sizes_sheet <- stage_samples %>%
  count(stage_group, name = "N") %>%
  mutate(`Included in comparison` = stage_group %in% STAGE_GROUPS) %>%
  rename(Group = stage_group)

run_group_analysis(
  variants         = all_variants,
  groups_df        = stage_samples %>% select(sample, group = stage_group),
  group_levels     = STAGE_GROUPS,
  test_levels      = STAGE_GROUPS,
  sample_sheet     = stage_samples %>%
    rename(Sample = sample, `Stage (raw)` = stage_raw, `Stage group` = stage_group),
  group_size_sheet = stage_sizes_sheet,
  notes            = tibble(Note = c(
    "Comparison: Stage I vs Stage III (stage sub-groups such as IB or IIIA-C are pooled).",
    "Samples with other or unknown stage are excluded (see Group_sizes).",
    if (RUN_FISHER) "Fisher p = two-sided Fisher's exact test, Stage I vs Stage III; BH-adjusted within each sheet.",
    common_notes)),
  out_file         = file.path(dir_out, paste0("TSO500_Frequency_StageI_vs_StageIII_", ts, ".xlsx"))
)

## ============================================================
## 8. Output 2: Molecular subtype
## ============================================================
subtype_test_levels <- if (SUBTYPE_TEST_EXCLUDES_UNKNOWN) {
  setdiff(SUBTYPE_LEVELS, "Unknown")
} else {
  SUBTYPE_LEVELS
}

subtype_sizes_sheet <- subtype_samples %>%
  count(subtype_group, name = "N") %>%
  mutate(subtype_group = factor(subtype_group, levels = SUBTYPE_LEVELS)) %>%
  complete(subtype_group, fill = list(N = 0L)) %>%
  arrange(subtype_group) %>%
  mutate(subtype_group = as.character(subtype_group),
         `Included in Fisher test` = subtype_group %in% subtype_test_levels) %>%
  rename(Group = subtype_group)

run_group_analysis(
  variants         = all_variants,
  groups_df        = subtype_samples %>% select(sample, group = subtype_group),
  group_levels     = SUBTYPE_LEVELS,
  test_levels      = subtype_test_levels,
  sample_sheet     = subtype_samples %>%
    rename(Sample = sample, `Molecular subtype (raw)` = subtype_raw,
           `Molecular subtype group` = subtype_group),
  group_size_sheet = subtype_sizes_sheet,
  notes            = tibble(Note = c(
    "Groups: HPV-associated, HPV-independent/p53-wildtype, HPV-independent/p53-mutant, Unknown.",
    "'p53-abnormal' in the IHC file is reported as 'p53-mutant'.",
    "Samples missing from the subtype file are counted as Unknown.",
    if (RUN_FISHER) paste0("Fisher p = Fisher's exact test across: ",
                           paste(subtype_test_levels, collapse = ", "),
                           "; BH-adjusted within each sheet."),
    common_notes)),
  out_file         = file.path(dir_out, paste0("TSO500_Frequency_MolecularSubtype_", ts, ".xlsx"))
)

message("Done.")