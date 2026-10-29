###############################################
# TSO500 CO-OCCURRENCE ANALYSIS (FREQUENCY ONLY)
###############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(stringr)
  library(readr)
  library(writexl)
})

## ============================================================
## 1) Paths
## ============================================================
dir_small <- "<SMALL_VARIANT_DIR>"
dir_cnv   <- "<COMBINED_VARIANT_DIR>"

out_root  <- "<OUTPUT_DIR>"
out_dir   <- file.path(out_root, "CoOccurrence_ExcelOnly")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

## ============================================================
## 2) Global settings
## ============================================================
DROP_SYNONYMOUS <- TRUE

N_TOTAL_SAMPLES <- 39

MIN_SNV_SAMPLES_PER_GENE <- 2
MIN_CNV_SAMPLES_PER_GENE <- 2
MIN_PAIR_N_BOTH          <- 2

## ============================================================
## 3) Helpers
## ============================================================
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  ifelse(grepl("^\\d+$", x), paste0("PQ", as.integer(x)), x)
}

snpeff_to_class <- function(term) {
  term <- str_trim(tolower(term))
  dplyr::case_when(
    str_detect(term, "frameshift")                     ~ "Frameshift",
    str_detect(term, "stop_gained|nonsense")           ~ "Nonsense",
    str_detect(term, "intron")                         ~ "Intron",
    str_detect(term, "splice|polypyrimidine_tract")    ~ "Splice",
    str_detect(term, "inframe")                        ~ "Inframe",
    str_detect(term, "missense")                       ~ "Missense",
    str_detect(term, "synonymous")                     ~ "Synonymous",
    TRUE                                               ~ "Other"
  )
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

pad_sample_cols <- function(mat, samples_all) {
  if (is.null(mat) || nrow(mat) == 0) {
    return(matrix(0L, nrow = 0, ncol = length(samples_all),
                  dimnames = list(NULL, samples_all)))
  }
  missing <- setdiff(samples_all, colnames(mat))
  if (length(missing) > 0) {
    add <- matrix(0L, nrow(mat), length(missing),
                  dimnames = list(rownames(mat), missing))
    mat <- cbind(mat, add)
  }
  mat[, samples_all, drop = FALSE]
}

make_gene_sample_mat <- function(df) {
  if (nrow(df) == 0) return(matrix(0L, 0, 0))
  df %>%
    distinct(sample, gene) %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = sample, values_from = val, values_fill = 0L) %>%
    column_to_rownames("gene") %>%
    as.matrix()
}

pairwise_square_counts <- function(bin_mat, min_pair_n = MIN_PAIR_N_BOTH,
                                   n_total = N_TOTAL_SAMPLES) {
  genes_use <- rownames(bin_mat)
  n         <- length(genes_use)
  if (n < 2) stop("Need >=2 genes for pairwise square analysis.")
  
  pair_list <- vector("list", 0)
  
  for (i in seq_len(n - 1)) {
    xi <- bin_mat[genes_use[i], ]
    for (j in (i + 1):n) {
      xj   <- bin_mat[genes_use[j], ]
      both <- sum(xi == 1 & xj == 1)
      if (both < min_pair_n) next
      n1 <- sum(xi == 1)
      n2 <- sum(xj == 1)
      ci_both <- wilson_ci(both, n_total)
      pair_list[[length(pair_list) + 1]] <- tibble(
        gene1          = genes_use[i],
        gene2          = genes_use[j],
        n_gene1        = n1,
        n_gene2        = n2,
        n_both         = both,
        n_either       = n1 + n2 - both,
        n_total_cohort = n_total,
        pct_cohort     = round(100 * both / n_total, 1),
        ci95_lower     = round(100 * ci_both$lower, 1),
        ci95_upper     = round(100 * ci_both$upper, 1),
        pct_of_gene1   = round(100 * both / n1, 1),
        pct_of_gene2   = round(100 * both / n2, 1),
        jaccard        = round(both / (n1 + n2 - both), 3),
        samples_with_both = paste(sort(colnames(bin_mat)[xi == 1 & xj == 1]), collapse = ", ")
      )
    }
  }
  
  if (length(pair_list) == 0) return(tibble())
  bind_rows(pair_list) %>%
    arrange(desc(n_both), desc(jaccard), gene1, gene2)
}

pairwise_rect_tbl <- function(row_mat, col_mat, col_a, col_b,
                              min_pair_n = MIN_PAIR_N_BOTH,
                              n_total = N_TOTAL_SAMPLES) {
  g1 <- rownames(row_mat); g2 <- rownames(col_mat)
  
  pair_list <- vector("list", 0)
  
  for (i in seq_along(g1)) {
    x  <- row_mat[g1[i], ]
    n1 <- sum(x == 1)
    for (j in seq_along(g2)) {
      y    <- col_mat[g2[j], ]
      n2   <- sum(y == 1)
      both <- sum(x == 1 & y == 1)
      if (both < min_pair_n) next
      ci_both <- wilson_ci(both, n_total)
      pair_list[[length(pair_list) + 1]] <- tibble(
        !!col_a := g1[i],
        !!col_b := g2[j],
        n_a            = n1,
        n_b            = n2,
        n_both         = both,
        n_total_cohort = n_total,
        pct_cohort     = round(100 * both / n_total, 1),
        ci95_lower     = round(100 * ci_both$lower, 1),
        ci95_upper     = round(100 * ci_both$upper, 1),
        pct_of_a       = round(100 * both / n1, 1),
        pct_of_b       = round(100 * both / n2, 1),
        jaccard        = round(both / (n1 + n2 - both), 3),
        samples_with_both = paste(sort(colnames(row_mat)[x == 1 & y == 1]), collapse = ", ")
      )
    }
  }
  
  if (length(pair_list) == 0) return(tibble())
  bind_rows(pair_list) %>%
    arrange(desc(n_both), desc(jaccard), .data[[col_a]], .data[[col_b]])
}

## ============================================================
## 4) Read SNVs
## ============================================================
read_small_variants <- function(file) {
  sample_id <- sub("_filtered_[0-9]{8}\\.xlsx$", "", basename(file))
  sample_id <- normalize_sample_id(sample_id)
  
  dat <- readxl::read_xlsx(file, .name_repair = "minimal") %>%
    janitor::clean_names()
  
  gene_col  <- grep("^gene(\\.|$|_)", names(dat), value = TRUE)[1]
  coseq_col <- grep("^coseq$|^consequence$|^variant_consequence$|consequence",
                    names(dat), value = TRUE)[1]
  pass_col  <- grep("^pass$", names(dat), value = TRUE)[1]
  
  if (is.na(gene_col) || is.na(coseq_col)) {
    message("Skipping (missing gene/coseq): ", basename(file))
    return(tibble(sample = character(), gene = character()))
  }
  if (!is.na(pass_col))
    dat <- dat %>% filter(.data[[pass_col]] == TRUE | .data[[pass_col]] == "PASS")
  
  tibble(
    sample = sample_id,
    gene   = toupper(str_squish(as.character(dat[[gene_col]]))),
    coseq  = as.character(dat[[coseq_col]])
  ) %>%
    filter(!is.na(gene), gene != "", !is.na(coseq), coseq != "") %>%
    mutate(alteration = snpeff_to_class(coseq)) %>%
    { if (DROP_SYNONYMOUS) filter(., alteration != "Synonymous") else . } %>%
    distinct(sample, gene)
}

small_files <- list.files(
  dir_small, pattern = "^PQ[0-9]+_filtered_[0-9]{8}\\.xlsx$", full.names = TRUE
)
stopifnot(length(small_files) > 0)
snv_long <- purrr::map_dfr(small_files, read_small_variants)
message("SNV calls: rows=", nrow(snv_long),
        " samples=", n_distinct(snv_long$sample),
        " genes=", n_distinct(snv_long$gene))

## ============================================================
## 5) Read AMP / DEL from the combined variant file
## ============================================================
split_row_trim_trailing <- function(l) {
  parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
  while (length(parts) > 0 && trimws(utils::tail(parts, 1)) == "") {
    parts <- parts[-length(parts)]
  }
  parts
}

get_section_table <- function(lines, section_name) {
  section_idx   <- which(str_detect(lines, "^\\[.*\\]\\s*$"))
  section_names <- str_match(lines[section_idx], "^\\[(.*)\\]\\s*$")[, 2]
  hit <- which(section_names == section_name)
  if (length(hit) == 0) return(NULL)
  
  start <- section_idx[hit[1]] + 1
  later <- section_idx[section_idx > section_idx[hit[1]]]
  end   <- if (length(later) > 0) min(later) - 1 else length(lines)
  if (start > end) return(NULL)
  
  seg <- lines[start:end]
  seg <- seg[trimws(seg) != ""]
  if (length(seg) < 2) return(NULL)
  
  rows   <- lapply(seg, split_row_trim_trailing)
  header <- rows[[1]]
  ncol_h <- length(header)
  body   <- lapply(rows[-1], function(r) {
    if (length(r) < ncol_h) r <- c(r, rep(NA_character_, ncol_h - length(r)))
    r[seq_len(ncol_h)]
  })
  
  df <- as.data.frame(do.call(rbind, body), stringsAsFactors = FALSE)
  names(df) <- header
  janitor::clean_names(df)
}

classify_cnv <- function(x) {
  x <- toupper(str_remove_all(str_squish(as.character(x)), "[<>]"))
  dplyr::case_when(
    x %in% c("DUP", "AMP", "AMPLIFICATION") ~ "AMP",
    x %in% c("DEL", "DELETION")             ~ "DEL",
    TRUE                                    ~ NA_character_
  )
}

classify_lr <- function(x) {
  x <- toupper(str_squish(as.character(x)))
  dplyr::case_when(
    x == "GAIN" ~ "AMP",
    x == "LOSS" ~ "DEL",
    TRUE        ~ NA_character_
  )
}

empty_cnv <- function() tibble(sample = character(), gene = character(),
                               cnv_type = character(), source = character(),
                               raw_call = character())

read_amp_del <- function(file) {
  sample_id <- normalize_sample_id(str_extract(basename(file), "^PQ[0-9]+"))
  lines     <- readr::read_lines(file)
  out       <- list()
  
  cnv_df <- get_section_table(lines, "Copy Number Variants")
  if (!is.null(cnv_df) && all(c("gene", "copy_number_variant") %in% names(cnv_df))) {
    out$cnv <- tibble(
      sample   = sample_id,
      gene     = toupper(str_squish(cnv_df$gene)),
      raw_call = cnv_df$copy_number_variant,
      source   = "Copy Number Variants"
    ) %>% mutate(cnv_type = classify_cnv(raw_call))
  } else if (!is.null(cnv_df)) {
    warning(sample_id, ": [Copy Number Variants] has no gene/copy_number_variant column. ",
            "Columns: ", paste(names(cnv_df), collapse = ", "))
  }
  
  lr_df <- get_section_table(lines, "Large Rearrangements")
  if (!is.null(lr_df) && all(c("gene", "cnv_type") %in% names(lr_df))) {
    out$lr <- tibble(
      sample   = sample_id,
      gene     = toupper(str_squish(lr_df$gene)),
      raw_call = lr_df$cnv_type,
      source   = "Large Rearrangements"
    ) %>%
      filter(!is.na(raw_call), trimws(raw_call) != "", toupper(trimws(raw_call)) != "NA") %>%
      mutate(cnv_type = classify_lr(raw_call))
  } else if (!is.null(lr_df)) {
    warning(sample_id, ": [Large Rearrangements] has no gene/cnv_type column. ",
            "Columns: ", paste(names(lr_df), collapse = ", "))
  }
  
  if (length(out) == 0) return(empty_cnv())
  bind_rows(out) %>%
    filter(!is.na(gene), gene != "") %>%
    select(sample, gene, cnv_type, source, raw_call)
}

cnv_files <- list.files(
  dir_cnv, pattern = "^PQ[0-9]+_CombinedVariantOutput_filtered\\.txt$", full.names = TRUE
)
stopifnot(length(cnv_files) > 0)
cnv_raw <- purrr::map_dfr(cnv_files, read_amp_del)

unmapped <- cnv_raw %>% filter(is.na(cnv_type)) %>% count(source, raw_call)
if (nrow(unmapped) > 0) {
  message("\nNOTE: these calls were NOT classified as AMP or DEL and were excluded:")
  print(unmapped)
}

cnv_long <- cnv_raw %>%
  filter(!is.na(cnv_type)) %>%
  distinct(sample, gene, cnv_type, .keep_all = TRUE)

message("\nAMP/DEL calls by source:")
print(cnv_long %>% count(source, cnv_type))

## ============================================================
## 6) Harmonise samples
## ============================================================
snv_file_samples <- normalize_sample_id(sub("_filtered_[0-9]{8}\\.xlsx$", "", basename(small_files)))
cnv_file_samples <- normalize_sample_id(str_extract(basename(cnv_files), "^PQ[0-9]+"))

samples_all <- sort(intersect(snv_file_samples, cnv_file_samples))

no_cnv_file <- setdiff(snv_file_samples, cnv_file_samples)
no_snv_file <- setdiff(cnv_file_samples, snv_file_samples)
if (length(no_cnv_file) > 0) message("No combined variant file for: ", paste(no_cnv_file, collapse = ", "))
if (length(no_snv_file) > 0) message("No small-variant file for: ",   paste(no_snv_file, collapse = ", "))

snv_long <- snv_long %>% filter(sample %in% samples_all)
cnv_long <- cnv_long %>% filter(sample %in% samples_all)

message("Samples analysed: ", length(samples_all))
if (length(samples_all) == 0) stop("No samples with both an SNV file and a combined variant file.")

if (length(samples_all) != N_TOTAL_SAMPLES) {
  warning(
    "Samples analysed (", length(samples_all),
    ") does not match the fixed cohort denominator N_TOTAL_SAMPLES (",
    N_TOTAL_SAMPLES, "). Percentages/CIs are still calculated against ",
    "N_TOTAL_SAMPLES = ", N_TOTAL_SAMPLES, "."
  )
}

## ============================================================
## 7) Build binary gene x sample matrices
## ============================================================
snv_mat <- snv_long %>% make_gene_sample_mat()
amp_mat <- cnv_long %>% filter(cnv_type == "AMP") %>% select(sample, gene) %>% make_gene_sample_mat()
del_mat <- cnv_long %>% filter(cnv_type == "DEL") %>% select(sample, gene) %>% make_gene_sample_mat()

snv_mat <- pad_sample_cols(snv_mat, samples_all)
amp_mat <- pad_sample_cols(amp_mat, samples_all)
del_mat <- pad_sample_cols(del_mat, samples_all)

snv_mat <- snv_mat[rowSums(snv_mat) >= MIN_SNV_SAMPLES_PER_GENE, , drop = FALSE]
amp_mat <- amp_mat[rowSums(amp_mat) >= MIN_CNV_SAMPLES_PER_GENE, , drop = FALSE]
del_mat <- del_mat[rowSums(del_mat) >= MIN_CNV_SAMPLES_PER_GENE, , drop = FALSE]

message(sprintf("Matrix dims: SNV=%dx%d  AMP=%dx%d  DEL=%dx%d",
                nrow(snv_mat), ncol(snv_mat),
                nrow(amp_mat), ncol(amp_mat),
                nrow(del_mat), ncol(del_mat)))

## ============================================================
## 7b) Per-gene alteration frequency table (fixed denominator, Wilson CI)
## ============================================================
gene_freq_tbl <- function(mat, label, n_total = N_TOTAL_SAMPLES) {
  if (nrow(mat) == 0) return(tibble())
  n_samples <- as.integer(rowSums(mat))
  ci        <- wilson_ci(n_samples, n_total)
  tibble(
    gene           = rownames(mat),
    alteration     = label,
    n_samples      = n_samples,
    n_total_cohort = n_total,
    pct_cohort     = round(100 * n_samples / n_total, 1),
    ci95_lower     = round(100 * ci$lower, 1),
    ci95_upper     = round(100 * ci$upper, 1)
  ) %>%
    arrange(desc(n_samples), gene)
}

gene_frequency <- bind_rows(
  gene_freq_tbl(snv_mat, "SNV"),
  gene_freq_tbl(amp_mat, "AMP"),
  gene_freq_tbl(del_mat, "DEL")
)

## ============================================================
## 8) Run all six pairwise co-alteration comparisons
## ============================================================
empty_note <- function(label) tibble(note = paste0("No pairs met MIN_PAIR_N_BOTH for ", label, "."))
or_note    <- function(t, label) if (nrow(t) == 0) empty_note(label) else t

snv_snv <- if (nrow(snv_mat) >= 2) or_note(pairwise_square_counts(snv_mat), "SNV x SNV") else empty_note("SNV x SNV")
amp_amp <- if (nrow(amp_mat) >= 2) or_note(pairwise_square_counts(amp_mat), "AMP x AMP") else empty_note("AMP x AMP")
del_del <- if (nrow(del_mat) >= 2) or_note(pairwise_square_counts(del_mat), "DEL x DEL") else empty_note("DEL x DEL")

snv_amp <- if (nrow(snv_mat) >= 1 && nrow(amp_mat) >= 1)
  or_note(pairwise_rect_tbl(snv_mat, amp_mat, "snv_gene", "amp_gene"), "SNV x AMP") else empty_note("SNV x AMP")
snv_del <- if (nrow(snv_mat) >= 1 && nrow(del_mat) >= 1)
  or_note(pairwise_rect_tbl(snv_mat, del_mat, "snv_gene", "del_gene"), "SNV x DEL") else empty_note("SNV x DEL")
amp_del <- if (nrow(amp_mat) >= 1 && nrow(del_mat) >= 1)
  or_note(pairwise_rect_tbl(amp_mat, del_mat, "amp_gene", "del_gene"), "AMP x DEL") else empty_note("AMP x DEL")

## ============================================================
## 9) Write single Excel file
## ============================================================
out_xlsx <- file.path(out_dir, paste0("CoMutation_AllPairs_", ts, ".xlsx"))

writexl::write_xlsx(
  list(
    Gene_Frequency = gene_frequency,
    SNV_SNV        = snv_snv,
    SNV_AMP        = snv_amp,
    SNV_DEL        = snv_del,
    AMP_AMP        = amp_amp,
    DEL_DEL        = del_del,
    AMP_DEL        = amp_del,
    CNV_Calls      = cnv_long %>% arrange(sample, cnv_type, gene)
  ),
  out_xlsx
)

message("\nDONE. Excel file saved: ", out_xlsx)
message("All percentages/CIs calculated against fixed cohort N = ", N_TOTAL_SAMPLES)