###############################################
# CLINVAR ID COVERAGE AND ONCOGENICITY CALL COUNTS
###############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(readr)
  library(writexl)
})

## ============================================================
## 0. User settings
## ============================================================
denovo_dir    <- "<SMALL_VARIANT_DIR>"
recurrent_dir <- "<RECURRENT_SMALL_VARIANT_DIR>"

data_dir     <- "<OUTPUT_DIR>"
results_file <- file.path(data_dir, "somatic_clinvar_results.csv")
out_xlsx     <- file.path(data_dir, "clinvar_id_coverage_and_calls.xlsx")

recurrent_samples <- c("PQ7", "PQ14", "PQ21", "PQ33", "PQ46", "PQ47")

## ============================================================
## 1. Helpers
## ============================================================
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  ifelse(grepl("^\\d+$", x), paste0("PQ", as.integer(x)), x)
}

sample_id_from_filename <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  m <- str_match(stem, "^([A-Za-z0-9]+)_")[, 2]
  ifelse(is.na(m), stem, m)
}

scan_folder <- function(d) {
  if (!dir.exists(d)) {
    warning("Folder does not exist -- skipped: ", d)
    return(character(0))
  }
  found <- list.files(d, pattern = "\\.xlsx?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  found[!str_starts(basename(found), fixed("~$"))]
}

## ============================================================
## 2. Count variant rows with a ClinVar VariationID
## ============================================================
denovo_files    <- scan_folder(denovo_dir)
recurrent_files <- scan_folder(recurrent_dir)

message("[De novo] ", length(denovo_files), " file(s) found in ", denovo_dir)
message("[Recurrence] ", length(recurrent_files), " file(s) found in ", recurrent_dir)

all_files <- tibble(
  path   = c(denovo_files, recurrent_files),
  cohort = c(rep("De novo", length(denovo_files)), rep("Recurrence", length(recurrent_files)))
)

if (nrow(all_files) == 0) stop("No .xlsx/.xls files found in either folder -- check the paths above.")

read_one_file <- function(path, cohort) {
  df <- tryCatch(read_excel(path), error = function(e) {
    warning("Could not read ", basename(path), ": ", conditionMessage(e))
    NULL
  })
  if (is.null(df)) return(NULL)
  
  if (!"clinvar.sig" %in% names(df)) {
    warning(basename(path), " has no 'clinvar.sig' column -- treating all rows as having no ID.")
    df$clinvar.sig <- NA
  }
  
  tibble(
    sample    = sample_id_from_filename(path),
    cohort    = cohort,
    n_rows    = nrow(df),
    n_with_id = sum(!is.na(df$clinvar.sig))
  )
}

per_file <- map2_dfr(all_files$path, all_files$cohort, read_one_file)

per_file <- per_file %>%
  mutate(
    sample = normalize_sample_id(sample),
    cohort = ifelse(sample %in% recurrent_samples, "Recurrence", "De novo")
  )

## ============================================================
## 3. Summarise ID coverage, overall and per cohort
## ============================================================
summarise_coverage <- function(df, label) {
  tibble(
    Cohort                    = label,
    `Total variant rows`      = sum(df$n_rows),
    `Rows with ClinVar ID`    = sum(df$n_with_id),
    `Rows without ClinVar ID` = sum(df$n_rows) - sum(df$n_with_id),
    `% with ClinVar ID`       = round(100 * sum(df$n_with_id) / sum(df$n_rows), 1)
  )
}

id_coverage <- bind_rows(
  summarise_coverage(per_file, "All samples"),
  per_file %>% group_split(cohort) %>% map_dfr(~ summarise_coverage(.x, unique(.x$cohort)))
)

message(
  "\nClinVar ID coverage:\n",
  paste(capture.output(print(id_coverage, width = Inf)), collapse = "\n")
)

## ============================================================
## 4. Oncogenicity call counts
## ============================================================
if (!file.exists(results_file)) {
  stop(
    results_file, " not found -- run tso500_clinvar_somatic.py first to ",
    "generate it before running this script."
  )
}

results <- read_csv(results_file, show_col_types = FALSE) %>%
  mutate(
    sample = normalize_sample_id(sample),
    cohort = ifelse(sample %in% recurrent_samples, "Recurrence", "De novo")
  )

required_cols <- c("sample", "gene", "somatic_calling", "classification_type")
missing_cols <- setdiff(required_cols, names(results))
if (length(missing_cols) > 0) {
  stop(
    results_file, " is missing required column(s): ", paste(missing_cols, collapse = ", "),
    ". Re-run tso500_clinvar_somatic.py to regenerate it."
  )
}

oncogenicity_rows <- results %>% filter(classification_type == "oncogenicity")

benign_terms     <- c("Benign", "Likely benign")
uncertain_terms  <- c("Uncertain significance")

call_counts <- function(df, label) {
  tibble(
    Cohort                              = label,
    `Total oncogenicity-classified rows` = nrow(df),
    `Benign / Likely benign`             = sum(df$somatic_calling %in% benign_terms),
    `Uncertain significance`             = sum(df$somatic_calling %in% uncertain_terms),
    `Oncogenic / Likely oncogenic`       = sum(df$somatic_calling %in% c("Oncogenic", "Likely oncogenic")),
    `Other / unrecognised term`          = nrow(df) -
      sum(df$somatic_calling %in% c(benign_terms, uncertain_terms, "Oncogenic", "Likely oncogenic"))
  )
}

calls_summary <- bind_rows(
  call_counts(oncogenicity_rows, "All samples"),
  oncogenicity_rows %>% group_split(cohort) %>% map_dfr(~ call_counts(.x, unique(.x$cohort)))
)

message(
  "\nOncogenicity call counts:\n",
  paste(capture.output(print(calls_summary, width = Inf)), collapse = "\n")
)

benign_uncertain_detail <- oncogenicity_rows %>%
  filter(somatic_calling %in% c(benign_terms, uncertain_terms)) %>%
  select(sample, cohort, gene, clinvar_variation_id, somatic_calling, cancer_type_condition) %>%
  arrange(cohort, sample, gene)

## ============================================================
## 5. Write Excel workbook
## ============================================================
write_xlsx(
  list(
    "ID Coverage"            = id_coverage,
    "Oncogenicity Call Counts" = calls_summary,
    "Benign+Uncertain Detail" = benign_uncertain_detail
  ),
  path = out_xlsx
)

message("\nSaved summary workbook to ", out_xlsx)