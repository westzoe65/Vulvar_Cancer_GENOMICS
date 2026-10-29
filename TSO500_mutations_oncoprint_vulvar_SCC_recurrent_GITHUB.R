###############################################
# TSO500 ONCOPRINT PIPELINE — RECURRENT DISEASE

## ============================================================
## 0. Libraries
## ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(ComplexHeatmap)
  library(grid)
  library(stringr)
  library(readr)
  library(writexl)
  library(tibble)
  library(circlize)
})

## ============================================================
## 0b. User toggles
## ============================================================
INCLUDE_SYNONYMOUS_IN_TABLES    <- TRUE
INCLUDE_SYNONYMOUS_IN_ONCOPRINT <- FALSE

PLOT_ALL_GENES <- FALSE
TOP_N_SNV      <- 50
TOP_N_DUP      <- 5
TOP_N_DEL      <- 5

BASE_FONT_SIZE   <- 11
ROW_FONT_SIZE    <- 11
ANN_FONT_SIZE    <- 10
TITLE_FONT_SIZE  <- 13
SAMPLE_FONT_SIZE <- 9

COL_GAP <- unit(0.5, "mm")
COL_WIDTH <- unit(11.65, "mm")

## ============================================================
## 1. Helper functions
## ============================================================
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  ifelse(grepl("^\\d+$", x), paste0("PQ", as.integer(x)), x)
}

read_event_annotation <- function(file,
                                  id_col_candidates,
                                  value_col_candidates,
                                  recode_map = NULL) {
  df <- readxl::read_xlsx(file, .name_repair = "minimal") %>%
    janitor::clean_names()
  
  id_clean  <- janitor::make_clean_names(id_col_candidates)
  val_clean <- janitor::make_clean_names(value_col_candidates)
  
  id_col  <- intersect(id_clean,  names(df))[1]
  val_col <- intersect(val_clean, names(df))[1]
  
  if (is.na(id_col) || is.na(val_col)) {
    stop(
      "Could not find required columns in: ", file, "\n",
      "Tried id columns: ",    paste(id_clean,  collapse = ", "), "\n",
      "Tried value columns: ", paste(val_clean, collapse = ", "), "\n",
      "Columns present:\n",    paste(names(df), collapse = ", ")
    )
  }
  
  out <- df %>%
    dplyr::select(sample = all_of(id_col), value = all_of(val_col)) %>%
    dplyr::mutate(
      sample = normalize_sample_id(sample),
      value  = as.character(value),
      value  = stringr::str_squish(value),
      value  = dplyr::na_if(value, "")
    )
  
  if (!is.null(recode_map)) {
    out <- out %>% dplyr::mutate(value = dplyr::recode(value, !!!recode_map, .default = value))
  }
  out
}

align_to_onco <- function(df, sample_order, label = "annotation") {
  df2 <- df %>%
    mutate(
      sample = normalize_sample_id(sample),
      value  = as.character(value),
      value  = str_squish(value),
      value  = na_if(value, "")
    ) %>%
    group_by(sample) %>%
    summarise(
      value = {
        vv <- value[!is.na(value)]
        if (length(vv) == 0) NA_character_ else vv[length(vv)]
      },
      .groups = "drop"
    )
  
  v   <- setNames(df2$value, df2$sample)
  out <- unname(v[normalize_sample_id(sample_order)])
  
  if (any(is.na(out))) {
    message(label, ": missing values for ", sum(is.na(out)), " samples (left blank).")
  }
  out
}

replace_na_label <- function(x, label = "Not recorded") {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[is.na(x) | x == ""] <- label
  factor(x)
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

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "N/A", "Not recorded", "not recorded")] <- NA_character_
  suppressWarnings(as.numeric(str_extract(x, "[-+]?[0-9]*\\.?[0-9]+")))
}

get_type <- function(x) {
  x <- as.character(x)
  if (is.na(x) || x == "") character(0) else strsplit(x, ";")[[1]]
}

wilson_ci <- function(x, n, conf_level = 0.95) {
  z      <- qnorm(1 - (1 - conf_level) / 2)
  p_hat  <- x / n
  denom  <- 1 + (z^2) / n
  centre <- p_hat + (z^2) / (2 * n)
  margin <- z * sqrt((p_hat * (1 - p_hat) / n) + (z^2 / (4 * n^2)))
  lower  <- pmax(0, (centre - margin) / denom)
  upper  <- pmin(1, (centre + margin) / denom)
  list(lower = lower, upper = upper)
}

recode_tumour_size <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    x == "0" ~ "<2 cm",
    x == "1" ~ "2-4 cm",
    x == "2" ~ ">4 cm",
    TRUE     ~ "Not recorded"
  )
}

recode_nodes <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    x == "0" ~ "0",
    x == "1" ~ "1",
    x == "2" ~ "\u22652",
    TRUE     ~ "Not recorded"
  )
}

recode_lvsi <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    x %in% c("0", "No",  "no",  "NO",  "Absent",  "absent")  ~ "Absent",
    x %in% c("1", "Yes", "yes", "YES", "Present", "present") ~ "Present",
    TRUE ~ "Not recorded"
  )
}

recode_doi <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    x == "0" ~ "<2 mm",
    x == "1" ~ "\u22652 mm",
    TRUE     ~ "Not recorded"
  )
}

recode_hpv <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(
    x %in% c("positive", "pos", "1", "yes") ~ "Positive",
    x %in% c("negative", "neg", "0", "no")  ~ "Negative",
    TRUE ~ "Not recorded"
  )
}

recode_p53 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(
    x %in% c("abnormal") ~ "Abnormal",
    x %in% c("wild-type", "wildtype", "wt") ~ "Wild-type",
    TRUE ~ "Not recorded"
  )
}

pull_char_col <- function(df, candidates) {
  col <- intersect(candidates, names(df))[1]
  if (is.na(col)) {
    message("Patient characteristics: column not found among: ",
            paste(candidates, collapse = ", "), " — filling with NA")
    return(rep(NA_character_, nrow(df)))
  }
  as.character(df[[col]])
}

## ============================================================
## 2. Input directories + status file
## ============================================================
dir_small    <- "<SMALL_VARIANT_DIR>"
dir_combined <- "<COMBINED_VARIANT_DIR>"

status_file <- "<METASTATIC_STATUS_FILE>"

dir_out <- "<OUTPUT_DIR>"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
## ============================================================
## 3. Discover samples + numeric order
## ============================================================
small_files <- list.files(
  dir_small,
  pattern    = "^PQ[0-9]+_filtered_[0-9]{8}\\.xlsx$",
  full.names = TRUE
)

if (length(small_files) == 0) stop("No small variant xlsx files found in dir_small: ", dir_small)

sample_ids <- sub("_filtered_[0-9]{8}\\.xlsx$", "", basename(small_files))

all_combined_in_dir <- list.files(dir_combined, full.names = TRUE)
message("All files in dir_combined (", length(all_combined_in_dir), " total):")
message(paste(basename(all_combined_in_dir), collapse = "\n"))

combined_files_all <- list.files(
  dir_combined,
  pattern     = "CombinedVariantOutput",
  full.names  = TRUE,
  ignore.case = TRUE
)

message("\nMatched CombinedVariantOutput files (", length(combined_files_all), "):")
message(paste(basename(combined_files_all), collapse = "\n"))

combined_ids <- stringr::str_extract(basename(combined_files_all), "^PQ[0-9]+")
combined_map <- setNames(combined_files_all, toupper(combined_ids))

file_map <- tibble(
  sample        = sample_ids,
  small_file    = small_files,
  combined_file = unname(combined_map[toupper(sample_ids)])
)

message("\nfile_map combined_file column:")
print(file_map %>% select(sample, combined_file))

missing_combined <- file_map$sample[is.na(file_map$combined_file)]
if (length(missing_combined) > 0) {
  message("WARNING: No combined variant file found for: ",
          paste(missing_combined, collapse = ", "),
          "\nDUP/DEL events will be missing for these samples.")
} else {
  message("Combined variant files found for all ", nrow(file_map), " samples.")
}

sample_numbers <- suppressWarnings(as.integer(str_remove(file_map$sample, "^PQ")))
if (any(is.na(sample_numbers))) {
  stop("Sample names not in format PQ<number>: ",
       paste(file_map$sample[is.na(sample_numbers)], collapse = ", "))
}

file_map <- file_map %>%
  mutate(sample_num = sample_numbers) %>%
  arrange(sample_num)

sample_order_all <- file_map$sample

## ============================================================
## 3B. Read metastatic status and order samples
## ============================================================
status_raw <- readxl::read_xlsx(status_file, .name_repair = "minimal") %>%
  janitor::clean_names()

sample_col <- intersect(
  names(status_raw),
  c("sample", "randomisation_no", "randomization_no", "id", "patient")
)[1]
status_col <- intersect(
  names(status_raw),
  c("status", "metastatic_status", "met_status", "metastasis", "group")
)[1]

if (is.na(sample_col) || is.na(status_col)) {
  stop("Could not find sample/status columns in: ", status_file,
       "\nColumns present: ", paste(names(status_raw), collapse = ", "))
}

sample_status <- status_raw %>%
  transmute(sample = .data[[sample_col]], status = .data[[status_col]]) %>%
  mutate(
    sample = normalize_sample_id(sample),
    status = as.character(status) %>% str_squish() %>% na_if(""),
    status = case_when(
      status %in% c("Non metastatic", "Non-metastatic", "Non metastatic ",
                    "Nonmetastatic", "NON METASTATIC", "NON-METASTATIC",
                    "NONMETASTATIC") ~ "Non-met",
      status %in% c("Metastatic", "METASTATIC")                            ~ "Met",
      TRUE ~ status
    )
  ) %>%
  filter(status %in% c("Non-met", "Met")) %>%
  distinct(sample, status)

nonmet_samples <- sample_status %>% filter(status == "Non-met") %>%
  pull(sample) %>% normalize_sample_id()
met_samples    <- sample_status %>% filter(status == "Met") %>%
  pull(sample) %>% normalize_sample_id()

nonmet_samples <- sample_order_all[sample_order_all %in% nonmet_samples]
met_samples    <- sample_order_all[sample_order_all %in% met_samples]

unassigned_samples <- sample_order_all[
  !sample_order_all %in% c(nonmet_samples, met_samples)
]
if (length(unassigned_samples) > 0) {
  message("NOTE: ", length(unassigned_samples),
          " recurrence sample(s) not found in metastatic status file — appended at end: ",
          paste(unassigned_samples, collapse = ", "))
}

sample_order <- c(nonmet_samples, met_samples, unassigned_samples)

if (length(sample_order) == 0) stop("No recurrence samples found in: ", dir_small)

file_map <- file_map %>%
  filter(sample %in% sample_order) %>%
  arrange(match(sample, sample_order))

N_TOTAL_SAMPLES <- length(sample_order)

## ============================================================
## 4. Read variants
## ============================================================
read_small_variants <- function(file, sample_id) {
  dat <- readxl::read_xlsx(file, .name_repair = "minimal") %>%
    janitor::clean_names()
  
  gene_col  <- grep("^gene",                             names(dat), value = TRUE)[1]
  coseq_col <- grep("^coseq$|^consequence$|consequence", names(dat), value = TRUE)[1]
  
  if (is.na(gene_col) || is.na(coseq_col)) {
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  tibble(
    sample = sample_id,
    gene   = toupper(as.character(dat[[gene_col]])),
    coseq  = as.character(dat[[coseq_col]])
  ) %>%
    filter(!is.na(gene), gene != "", !is.na(coseq), coseq != "") %>%
    mutate(
      consequence_term = str_trim(tolower(coseq)),
      alteration_class = snpeff_to_class(consequence_term)
    ) %>%
    { if (!INCLUDE_SYNONYMOUS_IN_TABLES) filter(., alteration_class != "Synonymous") else . } %>%
    select(sample, gene, consequence_term, alteration_class)
}

read_cnv_variants <- function(file, sample_id) {
  if (is.na(file) || !file.exists(file)) {
    if (!is.na(file)) warning("CNV file missing for ", sample_id, ": ", file)
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  raw     <- readLines(file, warn = FALSE)
  cnv_idx <- grep("^\\[Copy Number Variants\\]", raw)
  if (length(cnv_idx) == 0) {
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  next_idx <- grep("^\\[.*\\]", raw)
  next_idx <- next_idx[next_idx > cnv_idx][1]
  end_idx  <- ifelse(is.na(next_idx), length(raw), next_idx - 1)
  
  cnv_df <- readr::read_tsv(I(raw[(cnv_idx + 1):end_idx]), show_col_types = FALSE) %>%
    janitor::clean_names()
  
  gene_col <- grep("^gene",               names(cnv_df), value = TRUE)[1]
  cnv_col  <- grep("copy_number_variant", names(cnv_df), value = TRUE)[1]
  
  if (is.na(gene_col) || is.na(cnv_col)) {
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  gene_vec <- toupper(as.character(cnv_df[[gene_col]]))
  cnv_vec  <- as.character(cnv_df[[cnv_col]])
  
  tibble(
    sample           = rep(sample_id, length(gene_vec)),
    gene             = gene_vec,
    consequence_term = NA_character_,
    alteration_class = case_when(
      str_detect(cnv_vec, regex("dup", TRUE)) ~ "Duplication",
      str_detect(cnv_vec, regex("del", TRUE)) ~ "Deletion",
      TRUE ~ NA_character_
    )
  ) %>%
    filter(!is.na(alteration_class), !is.na(gene), gene != "")
}

read_lr_variants <- function(file, sample_id) {
  if (is.na(file) || !file.exists(file)) {
    if (!is.na(file)) warning("Combined variant file missing for ", sample_id, ": ", file)
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  raw    <- readLines(file, warn = FALSE)
  lr_idx <- grep("^\\[Large Rearrangements\\]", raw)
  if (length(lr_idx) == 0) {
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  all_section_idx <- grep("^\\[.*\\]", raw)
  next_section    <- all_section_idx[all_section_idx > lr_idx][1]
  end_idx         <- ifelse(is.na(next_section), length(raw), next_section - 1)
  
  block_lines <- raw[(lr_idx + 1):end_idx]
  block_lines <- block_lines[nchar(trimws(block_lines)) > 0]
  
  if (length(block_lines) < 2) {
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  lr_df <- tryCatch(
    readr::read_tsv(I(block_lines), show_col_types = FALSE, name_repair = "minimal") %>%
      janitor::clean_names(),
    error = function(e) {
      warning("Could not parse [Large Rearrangements] for ", sample_id, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(lr_df) || nrow(lr_df) == 0) {
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  gene_col    <- grep("^gene",    names(lr_df), value = TRUE)[1]
  cnvtype_col <- grep("cnv_type", names(lr_df), value = TRUE)[1]
  
  if (is.na(gene_col)) {
    warning("No gene column in [Large Rearrangements] for ", sample_id)
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  if (is.na(cnvtype_col)) {
    warning("No CNV Type column in [Large Rearrangements] for ", sample_id)
    return(tibble(sample = character(), gene = character(),
                  consequence_term = character(), alteration_class = character()))
  }
  
  gene_vec    <- toupper(trimws(as.character(lr_df[[gene_col]])))
  cnvtype_vec <- trimws(as.character(lr_df[[cnvtype_col]]))
  
  tibble(
    sample           = rep(sample_id, length(gene_vec)),
    gene             = gene_vec,
    consequence_term = NA_character_,
    alteration_class = case_when(
      str_detect(toupper(cnvtype_vec), "GAIN") ~ "Duplication",
      str_detect(toupper(cnvtype_vec), "LOSS") ~ "Deletion",
      TRUE ~ NA_character_
    )
  ) %>%
    filter(!is.na(alteration_class), !is.na(gene), gene != "")
}

small_variants <- map2_dfr(file_map$small_file,    file_map$sample, read_small_variants)
cnv_variants   <- map2_dfr(file_map$combined_file, file_map$sample, read_cnv_variants)
lr_variants    <- map2_dfr(file_map$combined_file, file_map$sample, read_lr_variants)

all_variants <- bind_rows(small_variants, cnv_variants, lr_variants) %>%
  mutate(sample = normalize_sample_id(sample)) %>%
  filter(sample %in% sample_order)

message("\nVariant counts by alteration class:")
print(table(all_variants$alteration_class))

## ============================================================
## 5. Split alteration types
## ============================================================
snv_variants_all <- all_variants %>%
  mutate(gene = toupper(gene), alteration_class = as.character(alteration_class)) %>%
  filter(!is.na(gene), gene != "",
         !(alteration_class %in% c("Duplication", "Deletion")))

cnv_dup_all <- all_variants %>%
  mutate(gene = toupper(gene), alteration_class = as.character(alteration_class)) %>%
  filter(!is.na(gene), gene != "", alteration_class == "Duplication")

cnv_del_all <- all_variants %>%
  mutate(gene = toupper(gene), alteration_class = as.character(alteration_class)) %>%
  filter(!is.na(gene), gene != "", alteration_class == "Deletion")

snv_rank_table <- snv_variants_all %>%
  filter(alteration_class != "Synonymous") %>%
  distinct(sample, gene) %>%
  count(gene, name = "n_samples_snv") %>%
  mutate(freq_snv = n_samples_snv / length(sample_order)) %>%
  arrange(desc(n_samples_snv), gene)

top_snv_genes <- if (PLOT_ALL_GENES) snv_rank_table %>% pull(gene) else
  snv_rank_table %>% slice_head(n = TOP_N_SNV) %>% pull(gene)

plot_snv_variants <- snv_variants_all %>%
  filter(gene %in% top_snv_genes) %>%
  { if (!INCLUDE_SYNONYMOUS_IN_ONCOPRINT) filter(., alteration_class != "Synonymous") else . }

snv_gene_order <- snv_rank_table %>% filter(gene %in% top_snv_genes) %>% pull(gene)

dup_rank_table <- cnv_dup_all %>%
  distinct(sample, gene) %>%
  count(gene, name = "n_samples_dup") %>%
  mutate(freq_dup = n_samples_dup / length(sample_order)) %>%
  arrange(desc(n_samples_dup), gene)

dup_gene_order <- if (PLOT_ALL_GENES) dup_rank_table %>% pull(gene) else
  dup_rank_table %>% slice_head(n = TOP_N_DUP) %>% pull(gene)

del_rank_table <- cnv_del_all %>%
  distinct(sample, gene) %>%
  count(gene, name = "n_samples_del") %>%
  mutate(freq_del = n_samples_del / length(sample_order)) %>%
  arrange(desc(n_samples_del), gene)

del_gene_order <- if (PLOT_ALL_GENES) del_rank_table %>% pull(gene) else
  del_rank_table %>% slice_head(n = TOP_N_DEL) %>% pull(gene)

## ============================================================
## 6. Build matrices
## ============================================================
build_onco_mat <- function(long_df, gene_order, sample_order,
                           alteration_col   = "alteration_class",
                           fixed_alteration = NULL) {
  if (length(gene_order) == 0) {
    m <- matrix("", nrow = 0, ncol = length(sample_order))
    colnames(m) <- sample_order
    return(m)
  }
  
  df <- long_df %>% filter(gene %in% gene_order)
  if (!is.null(fixed_alteration)) df <- df %>% mutate(alteration_class = fixed_alteration)
  
  wide <- df %>%
    group_by(gene, sample) %>%
    summarise(alteration = paste(unique(.data[[alteration_col]]), collapse = ";"),
              .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = alteration, values_fill = "") %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  miss <- setdiff(sample_order, colnames(wide))
  if (length(miss) > 0) {
    wide <- cbind(wide,
                  matrix("", nrow = nrow(wide), ncol = length(miss),
                         dimnames = list(rownames(wide), miss)))
  }
  
  gene_use <- intersect(gene_order, rownames(wide))
  wide[gene_use, sample_order, drop = FALSE]
}

snv_mat <- build_onco_mat(plot_snv_variants, snv_gene_order, sample_order)
dup_mat <- build_onco_mat(cnv_dup_all, dup_gene_order, sample_order,
                          fixed_alteration = "Duplication")
del_mat <- build_onco_mat(cnv_del_all, del_gene_order, sample_order,
                          fixed_alteration = "Deletion")

samps <- colnames(snv_mat)

patient_labels <- paste0("Patient ", str_remove(samps, "^PQ"))

## ============================================================
## 7. Oncoprint colours
## ============================================================
alt_col <- c(
  Missense    = "#1E88E5",
  Inframe     = "#00ACC1",
  Splice      = "#8E24AA",
  Frameshift  = "#D81B60",
  Nonsense    = "#F48FB1",
  Duplication = "#26A69A",
  Deletion    = "#5E35B1",
  Intron      = "#90CAF9",
  Other       = "#B0BEC5",
  Synonymous  = "#ECEFF1"
)

alter_fun <- list(
  background = function(x, y, w, h)
    grid.rect(x, y, w * 0.92, h * 0.55,
              gp = gpar(fill = "#EEEEEE", col = "#FAFAFA", lwd = 0.8)),
  Missense    = function(x, y, w, h)
    grid.rect(x, y, w * 0.90, h * 0.35, gp = gpar(fill = alt_col["Missense"],    col = NA)),
  Inframe     = function(x, y, w, h)
    grid.rect(x, y, w * 0.90, h * 0.35, gp = gpar(fill = alt_col["Inframe"],     col = NA)),
  Splice      = function(x, y, w, h)
    grid.rect(x, y, w * 0.90, h * 0.35, gp = gpar(fill = alt_col["Splice"],      col = NA)),
  Frameshift  = function(x, y, w, h)
    grid.rect(x, y, w * 0.45, h * 0.35, gp = gpar(fill = alt_col["Frameshift"],  col = NA)),
  Nonsense    = function(x, y, w, h)
    grid.rect(x, y, w * 0.45, h * 0.35, gp = gpar(fill = alt_col["Nonsense"],    col = NA)),
  Intron      = function(x, y, w, h)
    grid.rect(x, y, w * 0.45, h * 0.35, gp = gpar(fill = alt_col["Intron"],      col = NA)),
  Other       = function(x, y, w, h)
    grid.rect(x, y, w * 0.90, h * 0.35, gp = gpar(fill = alt_col["Other"],       col = NA)),
  Duplication = function(x, y, w, h)
    grid.rect(x, y, w * 0.90, h * 0.35, gp = gpar(fill = alt_col["Duplication"], col = NA)),
  Deletion    = function(x, y, w, h)
    grid.rect(x, y, w * 0.90, h * 0.35, gp = gpar(fill = alt_col["Deletion"],    col = NA))
)

col_map_all <- alt_col
if (!INCLUDE_SYNONYMOUS_IN_ONCOPRINT) {
  alter_fun$Synonymous <- NULL
  col_map_all <- col_map_all[names(col_map_all) != "Synonymous"]
}

col_map_dup <- c(Duplication = alt_col["Duplication"])
col_map_del <- c(Deletion    = alt_col["Deletion"])

## ============================================================
## 8. Top annotation: mutation count barplot
## ============================================================
snv_long_counts <- snv_variants_all %>%
  filter(sample %in% samps, alteration_class != "Synonymous") %>%
  distinct(sample, gene, alteration_class) %>%
  count(sample, alteration_class, name = "n")

dup_long_counts <- cnv_dup_all %>%
  filter(sample %in% samps) %>% distinct(sample, gene) %>%
  count(sample, name = "n") %>% mutate(alteration_class = "Duplication")

del_long_counts <- cnv_del_all %>%
  filter(sample %in% samps) %>% distinct(sample, gene) %>%
  count(sample, name = "n") %>% mutate(alteration_class = "Deletion")

bar_long <- bind_rows(snv_long_counts, dup_long_counts, del_long_counts)

cat_order <- c(
  "Missense", "Inframe", "Splice", "Frameshift", "Nonsense",
  "Intron", "Other", "Duplication", "Deletion"
)
bar_long$alteration_class <- factor(bar_long$alteration_class, levels = cat_order)

bar_wide <- bar_long %>%
  pivot_wider(names_from = alteration_class, values_from = n, values_fill = 0)

missing_samps <- setdiff(samps, bar_wide$sample)
if (length(missing_samps) > 0) bar_wide <- bind_rows(bar_wide, tibble(sample = missing_samps))

for (cat in cat_order) {
  if (!cat %in% names(bar_wide)) bar_wide[[cat]] <- 0L
}

bar_wide <- bar_wide %>%
  mutate(across(all_of(cat_order), ~ replace_na(.x, 0L))) %>%
  filter(sample %in% samps) %>%
  arrange(match(sample, samps))

stopifnot(
  "bar_wide rows don't match samps — check sample IDs" =
    identical(bar_wide$sample, samps)
)

bar_mat           <- t(as.matrix(bar_wide[, cat_order]))
rownames(bar_mat) <- cat_order
colnames(bar_mat) <- samps

message("snv_mat ncol : ", ncol(snv_mat))
message("bar_mat ncol : ", ncol(bar_mat))
message("Columns match: ", identical(colnames(bar_mat), colnames(snv_mat)))

bar_fill_col <- c(
  Missense    = alt_col["Missense"],
  Inframe     = alt_col["Inframe"],
  Splice      = alt_col["Splice"],
  Frameshift  = alt_col["Frameshift"],
  Nonsense    = alt_col["Nonsense"],
  Intron      = alt_col["Intron"],
  Other       = alt_col["Other"],
  Duplication = alt_col["Duplication"],
  Deletion    = alt_col["Deletion"]
)

y_max <- max(colSums(bar_mat)) * 1.05

top_bar_anno <- HeatmapAnnotation(
  `Mutation count` = anno_barplot(
    t(bar_mat),
    gp        = gpar(fill = unname(bar_fill_col), col = NA),
    border    = TRUE,
    bar_width = 0.7,
    height    = unit(2.2, "cm"),
    ylim      = c(0, y_max),
    axis_param = list(side = "left", gp = gpar(fontsize = ANN_FONT_SIZE))
  ),
  which                = "column",
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = ANN_FONT_SIZE, fontface = "plain"),
  annotation_height    = unit(2.2, "cm")
)

top_bar_heatmap <- Heatmap(
  matrix("x", nrow = 1, ncol = length(samps), dimnames = list("x", samps)),
  col                 = c("x" = "#FFFFFF"),
  rect_gp             = gpar(col = NA, fill = NA),
  show_heatmap_legend = FALSE,
  show_row_names      = FALSE,
  show_column_names   = FALSE,
  height              = unit(0.01, "mm"),
  width               = COL_WIDTH * length(samps),
  border              = FALSE,
  column_gap          = COL_GAP,
  top_annotation      = top_bar_anno
)

## ============================================================
## 9. Clinical annotations: MSI, TMB
## ============================================================
file_msi   <- "<MSI_FILE>"
file_tmb   <- "<TMB_FILE>"

msi_df <- read_event_annotation(
  file_msi,
  id_col_candidates    = c("Randomisation No.", "Randomization No.", "Sample", "ID", "Patient_no.", "Patient No."),
  value_col_candidates = c("MSI", "MSI score", "MSI_Score", "MSI_Event"))

tmb_df <- read_event_annotation(
  file_tmb,
  id_col_candidates    = c("Randomisation No.", "Randomization No.", "Sample", "ID", "Patient_no.", "Patient No."),
  value_col_candidates = c("TMB", "Total TMB", "Total_TMB", "TMB_Event"))

## ============================================================
## 9b. Patient characteristics
## ============================================================
file_char <- "<PATIENT_CHARACTERISTICS_FILE>"

char_raw <- readxl::read_xlsx(file_char, .name_repair = "minimal") %>%
  janitor::clean_names()

char_pt_col <- intersect(names(char_raw), c("patient_no", "patient_number", "patient",
                                            "id", "randomisation_no", "randomization_no",
                                            "sample"))[1]
if (is.na(char_pt_col)) {
  stop("Could not find patient ID column in: ", file_char,
       "\nColumns present: ", paste(names(char_raw), collapse = ", "))
}

char_df_full <- tibble(
  sample              = normalize_sample_id(as.character(char_raw[[char_pt_col]])),
  `Tumour size`       = recode_tumour_size(pull_char_col(char_raw, c("tumour_size", "tumor_size", "size"))),
  `Positive nodes`    = recode_nodes(pull_char_col(char_raw, c("positive_nodes", "node_count", "nodes"))),
  LVSI                = recode_lvsi(pull_char_col(char_raw, c("lvsi", "lvi"))),
  `Depth of invasion` = recode_doi(pull_char_col(char_raw, c("depth_of_invasion", "doi", "depth"))),
  HPV                 = recode_hpv(pull_char_col(char_raw, c("hpv", "hpv_status"))),
  p53                 = recode_p53(pull_char_col(char_raw, c("p53", "tp53")))
) %>%
  distinct(sample, .keep_all = TRUE)

align_char <- function(col_name, levels_order = NULL) {
  v   <- setNames(char_df_full[[col_name]], char_df_full$sample)
  out <- unname(v[normalize_sample_id(samps)])
  out[is.na(out)] <- "Not recorded"
  if (!is.null(levels_order)) factor(out, levels = levels_order) else factor(out)
}

tumoursize_vec <- align_char("Tumour size",
                             levels_order = c("<2 cm", "2-4 cm", ">4 cm", "Not recorded"))
nodes_vec      <- align_char("Positive nodes",
                             levels_order = c("0", "1", "\u22652", "Not recorded"))
lvsi_vec       <- align_char("LVSI",
                             levels_order = c("Absent", "Present", "Not recorded"))
doi_vec        <- align_char("Depth of invasion",
                             levels_order = c("<2 mm", "\u22652 mm", "Not recorded"))
hpv_vec        <- align_char("HPV",
                             levels_order = c("Negative", "Positive", "Not recorded"))
p53_vec        <- align_char("p53",
                             levels_order = c("Wild-type", "Abnormal", "Not recorded"))

## ============================================================
## 9c. Age
## ============================================================
file_age <- "<AGE_FILE>"

age_raw <- readxl::read_xlsx(file_age, .name_repair = "minimal") %>%
  janitor::clean_names()

age_id_col  <- intersect(names(age_raw),
                         c("patient_no", "patient_number", "patient", "id",
                           "randomisation_no", "randomization_no", "sample"))[1]
age_val_col <- intersect(names(age_raw),
                         c("age", "age_at_diagnosis", "age_event",
                           "age_at_diagnosis_2", "age_diagnosis"))[1]

if (is.na(age_id_col) || is.na(age_val_col)) {
  stop(
    "Could not find required columns in age file: ", file_age, "\n",
    "Columns present: ", paste(names(age_raw), collapse = ", ")
  )
}

age_lookup <- age_raw %>%
  transmute(
    sample = normalize_sample_id(.data[[age_id_col]]),
    age    = suppressWarnings(as.numeric(as.character(.data[[age_val_col]])))
  ) %>%
  distinct(sample, .keep_all = TRUE)

age_vec <- age_lookup$age[match(normalize_sample_id(samps), age_lookup$sample)]
if (any(is.na(age_vec))) {
  message("Age: missing values for ", sum(is.na(age_vec)), " sample(s) — plotted as 0.")
}

age_plot <- age_vec
age_plot[is.na(age_plot) | !is.finite(age_plot)] <- 0

age_data_max <- max(age_vec, na.rm = TRUE)
age_max      <- ceiling(age_data_max / 10) * 10
age_axis_at  <- seq(0, age_max, by = 20)
age_thresh   <- mean(age_vec, na.rm = TRUE)

## ============================================================
## 10. MSI / TMB numeric vectors
## ============================================================
msi_max <- 35; msi_thresh <- 30
tmb_max <- 45; tmb_thresh <- 10

msi_num  <- to_numeric_safe(align_to_onco(msi_df, samps, "MSI"))
tmb_num  <- to_numeric_safe(align_to_onco(tmb_df, samps, "TMB"))
msi_plot <- msi_num; msi_plot[is.na(msi_plot) | !is.finite(msi_plot)] <- 0
tmb_plot <- tmb_num; tmb_plot[is.na(tmb_plot) | !is.finite(tmb_plot)] <- 0

## ============================================================
## 11. Bottom annotations
## ============================================================
na_col <- "#ECEFF1"

col_tumour_size <- c(
  "<2 cm"        = "#BBDEFB",
  "2-4 cm"       = "#42A5F5",
  ">4 cm"        = "#1565C0",
  "Not recorded" = na_col
)

col_doi <- c(
  "<2 mm"        = "#EDE7F6",
  "\u22652 mm"   = "#6A1B9A",
  "Not recorded" = na_col
)

col_lvsi <- c(
  "Absent"       = "#FFEBEE",
  "Present"      = "#B71C1C",
  "Not recorded" = na_col
)

col_nodes <- c(
  "0"            = "#E8F5E9",
  "1"            = "#66BB6A",
  "\u22652"      = "#1B5E20",
  "Not recorded" = na_col
)

col_hpv <- c(
  "Negative"     = "#E0F7FA",
  "Positive"     = "#00695C",
  "Not recorded" = na_col
)

col_p53 <- c(
  "Wild-type"    = "#FFF3E0",
  "Abnormal"     = "#E65100",
  "Not recorded" = na_col
)

bottom_anno_only <- HeatmapAnnotation(
  MSI = anno_barplot(
    msi_plot,
    ylim      = c(0, msi_max),
    gp        = gpar(fill = "#26A69A", col = NA),
    border    = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    height    = unit(2.2, "cm"),
    axis_param = list(
      at     = c(0, 10, 20, 30),
      labels = c("0", "10", "20", "30"),
      side   = "left",
      gp     = gpar(fontsize = ANN_FONT_SIZE)
    )
  ),
  TMB = anno_barplot(
    tmb_plot,
    ylim      = c(0, tmb_max),
    gp        = gpar(fill = "#D81B60", col = NA),
    border    = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    height    = unit(2.2, "cm"),
    axis_param = list(
      at     = c(0, 10, 20, 30, 40),
      labels = c("0", "10", "20", "30", "40"),
      side   = "left",
      gp     = gpar(fontsize = ANN_FONT_SIZE)
    )
  ),
  `Tumour size`       = tumoursize_vec,
  `Depth of invasion` = doi_vec,
  HPV                 = hpv_vec,
  p53                 = p53_vec,
  LVSI                = lvsi_vec,
  `Positive nodes`    = nodes_vec,
  Age = anno_points(
    age_plot,
    ylim      = c(0, age_max),
    pch       = 16,
    size      = unit(2, "mm"),
    gp        = gpar(col = "black", fill = "black"),
    border    = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    height    = unit(2.2, "cm"),
    axis_param = list(
      at     = age_axis_at,
      labels = as.character(age_axis_at),
      side   = "left",
      gp     = gpar(fontsize = ANN_FONT_SIZE)
    )
  ),
  col = list(
    `Tumour size`       = col_tumour_size,
    `Depth of invasion` = col_doi,
    HPV                 = col_hpv,
    p53                 = col_p53,
    LVSI                = col_lvsi,
    `Positive nodes`    = col_nodes
  ),
  which                = "column",
  gap                  = unit(rep(2, 8), "mm"),
  gp                   = gpar(col = "white", lwd = 0.8),
  border               = FALSE,
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = ANN_FONT_SIZE, fontface = "plain"),
  show_annotation_name = c(
    MSI = TRUE, TMB = TRUE,
    `Tumour size` = TRUE, `Depth of invasion` = TRUE,
    HPV = TRUE, p53 = TRUE, LVSI = TRUE,
    `Positive nodes` = TRUE,
    Age = TRUE
  ),
  simple_anno_size     = unit(0.32, "cm"),
  annotation_height    = unit.c(
    unit(2.2, "cm"),
    unit(2.2, "cm"),
    unit(0.32, "cm"),
    unit(0.32, "cm"),
    unit(0.32, "cm"),
    unit(0.32, "cm"),
    unit(0.32, "cm"),
    unit(0.32, "cm"),
    unit(2.2, "cm")
  )
)

patient_label_anno_only <- HeatmapAnnotation(
  Patient = anno_text(
    patient_labels,
    rot      = 90,
    just     = "right",
    location = 1,
    gp       = gpar(fontsize = SAMPLE_FONT_SIZE, fontface = "plain", col = "black")
  ),
  which                = "column",
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = ANN_FONT_SIZE, fontface = "plain"),
  show_annotation_name = TRUE,
  annotation_height    = unit(2.5, "cm")
)

bottom_heatmap <- Heatmap(
  matrix("x", nrow = 1, ncol = length(samps), dimnames = list("x", samps)),
  col                 = c("x" = "#FFFFFF"),
  rect_gp             = gpar(col = NA, fill = NA),
  show_heatmap_legend = FALSE,
  show_row_names      = FALSE,
  show_column_names   = FALSE,
  height              = unit(0.01, "mm"),
  width               = COL_WIDTH * length(samps),
  border              = FALSE,
  column_gap          = COL_GAP,
  bottom_annotation   = bottom_anno_only
)

patient_heatmap <- Heatmap(
  matrix("x", nrow = 1, ncol = length(samps), dimnames = list("x", samps)),
  col                 = c("x" = "#FFFFFF"),
  rect_gp             = gpar(col = NA, fill = NA),
  show_heatmap_legend = FALSE,
  show_row_names      = FALSE,
  show_column_names   = FALSE,
  height              = unit(0.01, "mm"),
  width               = COL_WIDTH * length(samps),
  border              = FALSE,
  column_gap          = COL_GAP,
  bottom_annotation   = patient_label_anno_only
)

## ============================================================
## 12. Legends
## ============================================================
alteration_legend <- Legend(
  title     = "Alterations",
  labels    = names(col_map_all),
  legend_gp = gpar(fill = unname(col_map_all), col = NA),
  ncol      = 1
)

## ============================================================
## 13. Build oncoprint heatmaps
## ============================================================
ROW_HEIGHT <- unit(5, "mm")

main_onco <- oncoPrint(
  snv_mat,
  alter_fun               = alter_fun,
  alter_fun_is_vectorized = FALSE,
  get_type                = get_type,
  col                     = col_map_all,
  top_annotation          = NULL,
  remove_empty_rows       = TRUE,
  remove_empty_columns    = FALSE,
  show_column_names       = FALSE,
  show_heatmap_legend     = FALSE,
  show_pct                = FALSE,
  row_names_side          = "right",
  row_names_gp            = gpar(fontsize = ROW_FONT_SIZE),
  column_title            = NULL,
  width                   = COL_WIDTH * ncol(snv_mat),
  border                  = FALSE,
  row_gap                 = unit(0.2, "mm"),
  column_gap              = COL_GAP,
  height                  = ROW_HEIGHT * nrow(snv_mat),
  row_order               = seq_len(nrow(snv_mat)),
  column_order            = seq_len(ncol(snv_mat))
)

make_sub_onco <- function(mat, col_map) {
  if (nrow(mat) == 0) return(NULL)
  oncoPrint(
    mat,
    alter_fun               = alter_fun,
    alter_fun_is_vectorized = FALSE,
    get_type                = get_type,
    col                     = col_map,
    remove_empty_rows       = TRUE,
    remove_empty_columns    = FALSE,
    show_column_names       = FALSE,
    row_names_side          = "right",
    row_names_gp            = gpar(fontsize = ROW_FONT_SIZE),
    column_title            = NULL,
    width                   = COL_WIDTH * ncol(mat),
    show_pct                = FALSE,
    border                  = FALSE,
    show_heatmap_legend     = FALSE,
    row_gap                 = unit(0.2, "mm"),
    column_gap              = COL_GAP,
    height                  = ROW_HEIGHT * nrow(mat),
    row_order               = seq_len(nrow(mat)),
    column_order            = seq_len(ncol(mat))
  )
}

dup_onco <- make_sub_onco(dup_mat, col_map_dup)
del_onco <- make_sub_onco(del_mat, col_map_del)

## ============================================================
## 14. Assemble heatmap list
## ============================================================
ht_list <- top_bar_heatmap %v% main_onco
if (!is.null(dup_onco)) ht_list <- ht_list %v% dup_onco
if (!is.null(del_onco)) ht_list <- ht_list %v% del_onco
ht_list <- ht_list %v% bottom_heatmap %v% patient_heatmap

n_gaps <- length(ht_list@ht_list) - 1

draw_args <- list(
  ht_list,
  merge_legends          = FALSE,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  heatmap_legend_list    = list(alteration_legend),
  ht_gap                 = unit(rep(2, n_gaps), "mm"),
  padding                = unit(c(5, 5, 5, 5), "mm")
)

decorate_plots <- function() {
  decorate_annotation("MSI", {
    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(msi_thresh / msi_max, "npc"),
               gp = gpar(col = "black", lty = 2, lwd = 1.2))
  })
  decorate_annotation("TMB", {
    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(tmb_thresh / tmb_max, "npc"),
               gp = gpar(col = "black", lty = 2, lwd = 1.2))
  })
  decorate_annotation("Age", {
    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(age_thresh / age_max, "npc"),
               gp = gpar(col = "black", lty = 2, lwd = 1.2))
  })
}

## ============================================================
## 14a. PDF
## ============================================================
n_snv_rows <- sum(rowSums(snv_mat != "") > 0)
n_dup_rows <- if (!is.null(dup_onco)) nrow(dup_mat) else 0
n_del_rows <- if (!is.null(del_onco)) nrow(del_mat) else 0
total_rows <- n_snv_rows + n_dup_rows + n_del_rows
pdf_height <- max(14, total_rows * 0.20 + 8) + 1

LEGEND_WIDTH_IN    <- 4.0
ROW_LABEL_WIDTH_IN <- 1.5
n_samps      <- length(samps)
col_width_in <- as.numeric(COL_WIDTH) / 25.4
pdf_width    <- col_width_in * n_samps + LEGEND_WIDTH_IN + ROW_LABEL_WIDTH_IN

message("PDF dimensions: ", round(pdf_width, 2), " x ", round(pdf_height, 1),
        " inches  (", n_samps, " samples at ", round(col_width_in, 3),
        " in/col — same as 39-sample master)")

ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
pdf_file <- file.path(dir_out, paste0("TSO500_Oncoprint_Recurrence_", ts, ".pdf"))

pdf(pdf_file, width = pdf_width, height = pdf_height)
do.call(ComplexHeatmap::draw, draw_args)
decorate_plots()
dev.off()
message("Saved PDF: ", pdf_file)

## ============================================================
## 14b. TIFF
## ============================================================
TIFF_DPI <- 300
tif_file <- sub("\\.pdf$", ".tif", pdf_file)

tiff(filename = tif_file, width = pdf_width, height = pdf_height,
     units = "in", res = TIFF_DPI, compression = "lzw")
do.call(ComplexHeatmap::draw, draw_args)
decorate_plots()
dev.off()
message("Saved TIFF: ", tif_file)

## ============================================================
## 14c. SVG
## ============================================================
svg_file <- sub("\\.pdf$", ".svg", pdf_file)

svg(filename = svg_file, width = pdf_width, height = pdf_height)
do.call(ComplexHeatmap::draw, draw_args)
decorate_plots()
dev.off()
message("Saved SVG: ", svg_file)

## ============================================================
## 15. Ranking tables (with Wilson score 95% CI)
## ============================================================

rank_df2 <- all_variants %>%
  mutate(sample           = normalize_sample_id(sample),
         gene             = toupper(as.character(gene)),
         alteration_class = as.character(alteration_class)) %>%
  filter(!is.na(gene), gene != "")

rank_any <- rank_df2 %>% distinct(sample, gene) %>%
  count(gene, name = "n_samples_any") %>%
  arrange(desc(n_samples_any), gene)
.ci_any <- wilson_ci(rank_any$n_samples_any, N_TOTAL_SAMPLES)
rank_any <- rank_any %>%
  mutate(
    freq_any  = n_samples_any / N_TOTAL_SAMPLES,
    ci_any_lo = .ci_any$lower,
    ci_any_hi = .ci_any$upper
  )

rank_snv <- rank_df2 %>%
  filter(!(alteration_class %in% c("Duplication", "Deletion")),
         alteration_class != "Synonymous") %>%
  distinct(sample, gene) %>%
  count(gene, name = "n_samples_snv") %>%
  arrange(desc(n_samples_snv), gene)
.ci_snv <- wilson_ci(rank_snv$n_samples_snv, N_TOTAL_SAMPLES)
rank_snv <- rank_snv %>%
  mutate(
    freq_snv  = n_samples_snv / N_TOTAL_SAMPLES,
    ci_snv_lo = .ci_snv$lower,
    ci_snv_hi = .ci_snv$upper
  )

rank_dup <- rank_df2 %>% filter(alteration_class == "Duplication") %>%
  distinct(sample, gene) %>%
  count(gene, name = "n_samples_dup") %>%
  arrange(desc(n_samples_dup), gene)
.ci_dup <- wilson_ci(rank_dup$n_samples_dup, N_TOTAL_SAMPLES)
rank_dup <- rank_dup %>%
  mutate(
    freq_dup  = n_samples_dup / N_TOTAL_SAMPLES,
    ci_dup_lo = .ci_dup$lower,
    ci_dup_hi = .ci_dup$upper
  )

rank_del <- rank_df2 %>% filter(alteration_class == "Deletion") %>%
  distinct(sample, gene) %>%
  count(gene, name = "n_samples_del") %>%
  arrange(desc(n_samples_del), gene)
.ci_del <- wilson_ci(rank_del$n_samples_del, N_TOTAL_SAMPLES)
rank_del <- rank_del %>%
  mutate(
    freq_del  = n_samples_del / N_TOTAL_SAMPLES,
    ci_del_lo = .ci_del$lower,
    ci_del_hi = .ci_del$upper
  )

rank_merged <- rank_any %>%
  full_join(rank_snv,    by = "gene") %>%
  full_join(rank_dup,    by = "gene") %>%
  full_join(rank_del,    by = "gene") %>%
  mutate(across(starts_with("n_samples"), ~ replace_na(.x, 0L))) %>%
  arrange(desc(n_samples_snv), desc(n_samples_dup),
          desc(n_samples_del), gene)

ts2      <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_xlsx <- file.path(dir_out, paste0("TSO500_Oncoprint_Ranking_Recurrence_", ts2, ".xlsx"))

write_xlsx(
  list(
    Ranking_All               = rank_merged,
    Ranking_SNV_Top25         = rank_snv %>% slice_head(n = TOP_N_SNV),
    Duplications              = rank_dup,
    Duplications_Top5         = rank_dup %>% slice_head(n = TOP_N_DUP),
    Deletions                 = rank_del,
    Deletions_Top5            = rank_del %>% slice_head(n = TOP_N_DEL)
  ),
  out_xlsx
)

message("Saved ranking tables: ", out_xlsx)
message("Outputs: ", pdf_file, " | ", tif_file, " | ", svg_file, " | ", out_xlsx)