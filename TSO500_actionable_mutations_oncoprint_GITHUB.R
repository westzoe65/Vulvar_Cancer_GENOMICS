###############################################
# TSO500 ACTIONABILITY (OncoKB + ESCAT) ONCOPRINT AND SUMMARY
###############################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(janitor)
  library(writexl)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## =========================
## 0) Palettes
## =========================
na_col <- "#ECEFF1"

ALT_COL <- "#1E88E5"

ESCAT_COL <- "#FB8C00"

db_pal <- c("OncoKB" = ALT_COL, "ESCAT" = ESCAT_COL)

base_pal <- c(
  "#E53935","#8E24AA","#1E88E5","#00ACC1",
  "#43A047","#FB8C00","#6D4C41","#757575",
  "#3949AB","#039BE5","#00897B","#7CB342",
  "#F4511E","#D81B60","#546E7A"
)

mut_pal <- c(
  "Missense"               = "#1E88E5",
  "Inframe"                = "#00ACC1",
  "Fusions"                = "#8E24AA",
  "Frameshift"             = "#D81B60",
  "Nonsense"               = "#F48FB1",
  "Duplication"            = "#26A69A",
  "Deletion"                = "#5E35B1",
  "Intron"                 = "#90CAF9",
  "Synonymous"             = "#B0BEC5",
  "Wild-type"              = "#90CAF9",
  "Amplification"          = "#1E88E5",
  "Loss"                   = "#B71C1C",
  "Oncogenic Mutation"     = "#039BE5",
  "Tumour Mutation Burden" = "#FF6F00",
  "Splice"                 = "#43A047",
  "Truncating"             = "#FB8C00",
  "SNV"                    = "#D81B60"
)

broad_pal <- c(
  "Thoracic"         = "#185FA5",
  "Breast"           = "#D4537E",
  "Gastrointestinal" = "#BA7517",
  "Hepatobiliary"    = "#854F0B",
  "Pancreatic"       = "#3B6D11",
  "Gynaecological"   = "#993556",
  "Genitourinary"    = "#0F6E56",
  "Head & Neck"      = "#534AB7",
  "Skin"             = "#D85A30",
  "CNS"              = "#7F77DD",
  "Haematological"   = "#A32D2D",
  "Histiocytic"      = "#1D9E75",
  "Sarcoma"          = "#6D4C41",
  "Endocrine"        = "#EF9F27",
  "Solid Tumour"     = "#43A047",
  "Other / Rare"     = "#546E7A",
  "Unknown"          = "#B0BEC5"
)

tc_strip_pal <- c(
  "Chemotherapy"     = "#EF9A9A",
  "Immunotherapy"    = "#A5D6A7",
  "Targeted Therapy" = "#90CAF9",
  "Combined"         = "#CE93D8"
)

broad_sort_order <- c(
  "Solid Tumour"=1,"Head & Neck"=2,"Gynaecological"=3,
  "Thoracic"=4,"Breast"=5,"Gastrointestinal"=6,
  "Hepatobiliary"=7,"Pancreatic"=8,"Genitourinary"=9,
  "Skin"=10,"CNS"=11,"Haematological"=12,
  "Histiocytic"=13,"Sarcoma"=14,"Endocrine"=15,
  "Other / Rare"=16,"Unknown"=17
)

risk_pal <- c(
  "Early-stage (node-negative)" = "#E3F2FD",
  "Late-stage (node-positive)"  = "#C2185B"
)

ROW_FONT_SIZE    <- 11
ANN_FONT_SIZE    <- 10
SAMPLE_FONT_SIZE <- 9
TITLE_FONT_SIZE  <- 13

## =========================
## 1) Input / output
## =========================
infile  <- "<ONCOKB_ESCAT_FILE>"
out_dir <- "<OUTPUT_DIR>"
ts      <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_pdf_denovo <- file.path(out_dir, paste0("OncoKB_Oncoprint_DeNovo_",     ts, ".pdf"))
out_png_denovo <- file.path(out_dir, paste0("OncoKB_Oncoprint_DeNovo_",     ts, ".png"))
out_pdf_recur  <- file.path(out_dir, paste0("OncoKB_Oncoprint_Recurrence_", ts, ".pdf"))
out_png_recur  <- file.path(out_dir, paste0("OncoKB_Oncoprint_Recurrence_", ts, ".png"))
out_xlsx    <- file.path(out_dir, paste0("Actionable_Gene_Summary_", ts, ".xlsx"))
out_bar_pdf <- file.path(out_dir, paste0("SolidTumour_ActionableGenes_", ts, ".pdf"))
out_bar_png <- file.path(out_dir, paste0("SolidTumour_ActionableGenes_", ts, ".png"))

staging_file   <- "<STAGE_FILE>"
meta_id_col    <- "randomisation_no"
meta_stage_col <- "stage"

file_char_stage_recur_sort <- "<RECURRENT_PATIENT_CHARACTERISTICS_FILE>"
id_candidates_sort <- c("Randomisation No.", "Randomization No.", "Sample", "ID", "Patient_no.", "Patient No.")

normalize_pq_id <- function(x) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  ifelse(str_detect(x, "^\\d+$"), paste0("pq", as.integer(x)), x)
}

read_event_annotation <- function(file, id_col_candidates, value_col_candidates) {
  d <- read_xlsx(file, .name_repair="minimal") %>% clean_names()
  id_clean  <- janitor::make_clean_names(id_col_candidates)
  val_clean <- janitor::make_clean_names(value_col_candidates)
  id_col  <- intersect(id_clean,  names(d))[1]
  val_col <- intersect(val_clean, names(d))[1]
  if (is.na(id_col) || is.na(val_col)) {
    warning("Could not find id/value columns in: ", file,
            " - columns present: ", paste(names(d), collapse=", "))
    return(tibble(sample=character(), value=character()))
  }
  d %>%
    transmute(sample = normalize_pq_id(.data[[id_col]]),
              value  = str_squish(as.character(.data[[val_col]]))) %>%
    mutate(value = na_if(value, ""))
}

align_to_cols <- function(df, cols) {
  v <- setNames(df$value, df$sample)
  unname(v[normalize_pq_id(cols)])
}

stage_to_rank <- function(x) {
  x <- toupper(str_squish(as.character(x)))
  x[x %in% "IA"]   <- "1A"
  x[x %in% "IB"]   <- "1B"
  x[x %in% "IIIA"] <- "3A"
  x[x %in% "IIIB"] <- "3B"
  x[x %in% "IIIC"] <- "3C"
  match(x, c("1A","1B","3A","3B","3C"))
}

recur_pq_ids <- normalize_pq_id(c(7, 14, 21, 33, 46, 47))

## =========================
## 2) Read data
## =========================
find_col <- function(pat, cols) {
  h <- cols[str_detect(cols, regex(pat, ignore_case=TRUE))]
  if (length(h)==0) NA_character_ else h[1]
}

read_oncokb_df <- function(path) {
  d <- read_xlsx(path, sheet=1) %>% clean_names()
  col_map <- list(
    level_classification = find_col("level_class|classification", names(d)),
    database              = find_col("^database$",                 names(d)),
    level                = find_col("^level$",                    names(d)),
    therapy              = find_col("^therapy$|^drug$",           names(d)),
    therapy_class        = find_col("^therapy_class$|^drug_class",names(d)),
    therapy_subclass     = find_col("^therapy_subclass$|^drug_sub",names(d)),
    cancer_type          = find_col("cancer.*type|type.*cancer",  names(d)),
    cancer_class         = find_col("cancer.*class",              names(d)),
    broad_category       = find_col("broad",                      names(d)),
    gene                 = find_col("^gene$",                     names(d)),
    mutation             = find_col("^mutation$",                 names(d))
  )
  for (logical_name in names(col_map)) {
    raw_name <- col_map[[logical_name]]
    if (!is.na(raw_name) && raw_name %in% names(d) && raw_name != logical_name)
      names(d)[names(d) == raw_name] <- logical_name
  }
  if (!"broad_category" %in% names(d)) d$broad_category <- NA_character_
  if (!"cancer_class"   %in% names(d)) d$cancer_class   <- NA_character_
  if (!"database"       %in% names(d)) d$database       <- NA_character_
  d
}

df <- read_oncokb_df(infile)
cat("OncoKB+ESCAT clean names:\n"); print(names(df)); cat("\n")
cat("Database values found:\n"); print(table(df$database, useNA="ifany")); cat("\n")

df_recur <- df

classify_database <- function(x) {
  x <- toupper(trimws(as.character(x)))
  case_when(
    str_detect(x, "ONCOKB") ~ "OncoKB",
    str_detect(x, "ESCAT")  ~ "ESCAT",
    TRUE ~ NA_character_
  )
}

## =========================
## 3) Sample columns per cohort
## =========================
all_pq_cols <- names(df)[str_detect(names(df), "^pq\\d+$")]
if (length(all_pq_cols)==0) stop("No PQ columns found in OncoKB workbook (`infile`).")

cat("All pq columns in workbook:\n"); print(all_pq_cols)
cat("Recurrence IDs sought:\n");      print(recur_pq_ids)

recurrent_cols <- intersect(all_pq_cols, recur_pq_ids)
denovo_cols    <- setdiff(all_pq_cols, recur_pq_ids)
if (length(denovo_cols)==0) stop("No de novo PQ columns after removing recurrence samples.")

if (length(recurrent_cols)==0)
  warning("None of the recurrence IDs (", paste(recur_pq_ids, collapse=", "),
          ") matched a workbook column - check the pq column names printed above. ",
          "The recurrence oncoprint will be skipped.")

cat("De novo pq columns:   ", length(denovo_cols),    "  ->  ", paste(denovo_cols, collapse=", "),    "\n")
cat("Recurrence pq columns:", length(recurrent_cols), "  ->  ", paste(recurrent_cols, collapse=", "), "\n\n")

stage_denovo_sort <- read_event_annotation(staging_file,               id_candidates_sort, c("Stage"))
stage_recur_sort  <- read_event_annotation(file_char_stage_recur_sort, id_candidates_sort, c("Stage"))

denovo_rank <- stage_to_rank(align_to_cols(stage_denovo_sort, denovo_cols))
denovo_cols <- denovo_cols[order(denovo_rank, as.integer(str_remove(denovo_cols,"^pq")))]

if (length(recurrent_cols) > 0) {
  recurrent_rank <- stage_to_rank(align_to_cols(stage_recur_sort, recurrent_cols))
  recurrent_cols <- recurrent_cols[order(recurrent_rank, as.integer(str_remove(recurrent_cols,"^pq")))]
}

sample_cols      <- denovo_cols
n_total_patients <- length(sample_cols)
sample_names     <- paste0("Patient ", str_remove(sample_cols, "^pq"))

###############################################################################
## PART A - ONCOPRINT
###############################################################################

## =========================
## 4) Clean and filter (oncoprint)
## =========================
make_df2 <- function(df_src, cols) {
  d <- df_src %>%
    mutate(across(all_of(cols), ~suppressWarnings(as.numeric(.x)))) %>%
    mutate(
      gene_clean    = toupper(trimws(as.character(gene))),
      therapy_clean = trimws(as.character(therapy)),
      tc_clean      = trimws(as.character(therapy_class)),
      cancer_type_c = trimws(as.character(cancer_type)),
      broad_cat     = trimws(as.character(broad_category)),
      mut_class     = trimws(as.character(mutation)),
      level_str     = trimws(as.character(level)),
      lc_str        = trimws(as.character(level_classification)),
      level_num     = suppressWarnings(as.numeric(trimws(as.character(level)))),
      db_clean      = classify_database(database)
    ) %>%
    mutate(
      broad_cat     = if_else(is.na(broad_cat)|broad_cat %in% c("","NA"), "Unknown", broad_cat),
      broad_cat     = if_else(broad_cat=="Pan-Cancer", "Solid Tumour", broad_cat),
      therapy_clean = if_else(is.na(therapy_clean)|therapy_clean=="","Unknown Therapy",therapy_clean),
      tc_clean      = if_else(is.na(tc_clean)|tc_clean=="","Targeted Therapy",tc_clean),
      level_num     = if_else(is.na(level_num), 99, level_num),
      is_tmb    = str_detect(gene_clean, regex("tumou?r.mutation.burden|\\bTMB\\b", ignore_case=TRUE)),
      is_wt     = str_detect(mut_class,  regex("wild.?type|\\bWT\\b", ignore_case=TRUE)),
      mut_class = case_when(
        is_tmb ~ "Tumour Mutation Burden",
        is_wt  ~ "Wild-type",
        is.na(mut_class) ~ "SNV",
        TRUE   ~ mut_class
      ),
      gene_clean = case_when(
        is_tmb ~ "TMB",
        is.na(gene_clean)|gene_clean=="" ~ "UNKNOWN",
        TRUE ~ gene_clean
      )
    ) %>%
    filter(!str_detect(lc_str, regex("resist", ignore_case=TRUE))) %>%
    filter(gene_clean != "TMB") %>%
    filter(!is_wt) %>%
    filter(!is.na(db_clean)) %>%
    filter(
      (db_clean=="OncoKB" & level_num %in% c(1,2,3)) |
        (db_clean=="ESCAT")
    )
  
  he <- rowSums(d[, cols, drop=FALSE]==2, na.rm=TRUE) > 0
  d[he, , drop=FALSE]
}

df2 <- make_df2(df, denovo_cols)
if (nrow(df2)==0) stop("No eligible de novo rows.")
cat("De novo broad_cat after filter:\n")
print(sort(table(df2$broad_cat), decreasing=TRUE)); cat("\n")
cat("De novo rows by Database after filter:\n")
print(table(df2$db_clean)); cat("\n")

df2_recur <- if (is.null(df_recur) || length(recurrent_cols)==0) NULL else make_df2(df_recur, recurrent_cols)
if (!is.null(df2_recur)) {
  cat("Recurrence broad_cat after filter:\n")
  print(sort(table(df2_recur$broad_cat), decreasing=TRUE)); cat("\n")
  cat("Recurrence rows by Database after filter:\n")
  print(table(df2_recur$db_clean)); cat("\n")
}

## =========================
## 5) Tidy up mutation classes (both cohorts)
## =========================
tidy_mut_class <- function(d) {
  if (is.null(d)) return(d)
  d %>% mutate(mut_class = if_else(mut_class %in% names(mut_pal), mut_class, "SNV"))
}
df2       <- tidy_mut_class(df2)
df2_recur <- tidy_mut_class(df2_recur)

## =========================
## 6) One row per gene, per cohort
## =========================
build_cohort_elig_mat <- function(df2_src, cols, col_labels) {
  if (is.null(df2_src) || length(cols) == 0 || nrow(df2_src) == 0)
    return(matrix("", nrow=0, ncol=0))
  
  gf <- df2_src %>%
    select(gene_clean, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to="sample_col", values_to="v") %>%
    filter(!is.na(v), v==2) %>%
    distinct(gene_clean, sample_col) %>%
    count(gene_clean, name="n_alt") %>%
    arrange(desc(n_alt), gene_clean)
  
  grows <- gf$gene_clean
  if (length(grows)==0)
    return(matrix("", nrow=0, ncol=length(cols), dimnames=list(NULL, col_labels)))
  
  m <- matrix("", nrow=length(grows), ncol=length(cols),
              dimnames=list(grows, col_labels))
  for (gi in seq_along(grows)) {
    g      <- grows[gi]
    g_rows <- which(df2_src$gene_clean == g)
    for (j in seq_along(cols)) {
      types <- character(0)
      for (i in g_rows) {
        v <- df2_src[[cols[j]]][i]
        if (!is.na(v) && v==2) types <- c(types, df2_src$db_clean[i])
      }
      types <- unique(types)
      if (length(types)>0) m[gi, j] <- paste(types, collapse=";")
    }
  }
  m
}

denovo_labels <- paste0("Patient ", str_remove(denovo_cols,    "^pq"))
recur_labels  <- paste0("Patient ", str_remove(recurrent_cols, "^pq"))

elig_mat_denovo <- build_cohort_elig_mat(df2,       denovo_cols,    denovo_labels)
elig_mat_recur  <- build_cohort_elig_mat(df2_recur, recurrent_cols, recur_labels)

cat("De novo oncoprint genes:   ", nrow(elig_mat_denovo), "\n")
cat("Recurrence oncoprint genes:", nrow(elig_mat_recur),  "\n")

get_type <- function(x) {
  x <- as.character(x)
  if (is.na(x) || x=="") character(0) else strsplit(x, ";")[[1]]
}

## =========================
## 7) alter_fun: OncoKB and ESCAT colours
## =========================
alter_fun <- list(
  background = function(x, y, w, h)
    grid.rect(x, y, w*0.92, h*0.55, gp=gpar(fill=na_col, col=NA)),
  OncoKB = function(x, y, w, h)
    grid.rect(x, y, w*0.90, h*0.55, gp=gpar(fill=ALT_COL, col=NA)),
  ESCAT = function(x, y, w, h)
    grid.rect(x, y, w*0.45, h*0.55, gp=gpar(fill=ESCAT_COL, col=NA))
)
col <- db_pal

## =========================
## 8) Clinical and patient characteristic annotations
## =========================
to_numeric_safe <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "N/A", "Not recorded")] <- NA_character_
  suppressWarnings(as.numeric(str_extract(x, "[-+]?[0-9]*\\.?[0-9]+")))
}

recode_tumour_size <- function(x) case_when(
  x=="0" ~ "<2 cm", x=="1" ~ "2-4 cm", x=="2" ~ ">4 cm", TRUE ~ "Not recorded"
)
recode_nodes <- function(x) case_when(
  x=="0" ~ "0", x=="1" ~ "1", x=="2" ~ "\u22652", TRUE ~ "Not recorded"
)
recode_lvsi <- function(x) case_when(
  x %in% c("0","No","no","NO","Absent","absent")  ~ "Absent",
  x %in% c("1","Yes","yes","YES","Present","present") ~ "Present", TRUE ~ "Not recorded"
)
recode_doi <- function(x) case_when(
  x=="0" ~ "<2 mm", x=="1" ~ "\u22652 mm", TRUE ~ "Not recorded"
)
recode_hpv <- function(x) { x <- tolower(x); case_when(
  x %in% c("positive","pos","1","yes") ~ "Positive",
  x %in% c("negative","neg","0","no")  ~ "Negative", TRUE ~ "Not recorded"
) }
recode_p53 <- function(x) { x <- tolower(x); case_when(
  x %in% c("mutant","mut","1","mutated","abnormal") ~ "Abnormal",
  x %in% c("wild-type","wildtype","wt","0")          ~ "Wild-type", TRUE ~ "Not recorded"
) }

file_msi_denovo  <- "<MSI_FILE>"
file_tmb_denovo  <- "<TMB_FILE>"
file_char_denovo <- "<PATIENT_CHARACTERISTICS_FILE>"
file_age_denovo  <- "<AGE_FILE>"

file_msi_recur         <- "<RECURRENT_MSI_FILE>"
file_tmb_recur         <- "<RECURRENT_TMB_FILE>"
file_char_stage_recur  <- "<RECURRENT_PATIENT_CHARACTERISTICS_FILE>"
file_age_recur         <- "<RECURRENT_AGE_FILE>"

id_candidates <- c("Randomisation No.", "Randomization No.", "Sample", "ID", "Patient_no.", "Patient No.")

msi_denovo <- read_event_annotation(file_msi_denovo, id_candidates, c("MSI","MSI score","MSI_Score","MSI_Event"))
msi_recur  <- read_event_annotation(file_msi_recur,  id_candidates, c("MSI","MSI score","MSI_Score","MSI_Event"))

tmb_denovo <- read_event_annotation(file_tmb_denovo, id_candidates, c("TMB","Total TMB","Total_TMB","TMB_Event"))
tmb_recur  <- read_event_annotation(file_tmb_recur,  id_candidates, c("TMB","Total TMB","Total_TMB","TMB_Event"))

stage_denovo <- read_event_annotation(staging_file,          id_candidates, c("Stage"))
stage_recur  <- read_event_annotation(file_char_stage_recur, id_candidates, c("Stage"))

recode_stage_vec <- function(x) {
  s <- toupper(str_squish(as.character(x)))
  s[s %in% "IA"]   <- "1A"
  s[s %in% "IB"]   <- "1B"
  s[s %in% "IIIA"] <- "3A"
  s[s %in% "IIIB"] <- "3B"
  s[s %in% "IIIC"] <- "3C"
  s[!(s %in% c("1A","1B","3A","3B","3C"))] <- NA_character_
  factor(s, levels=c("1A","1B","3A","3B","3C"))
}

read_char_file <- function(file) {
  d <- read_xlsx(file, .name_repair="minimal") %>% clean_names()
  pt_col <- intersect(names(d), c("patient_no","patient_number","patient","id",
                                  "randomisation_no","randomization_no","sample"))[1]
  if (is.na(pt_col)) { warning("No patient ID column found in: ", file); return(NULL) }
  pull_col <- function(cands) {
    c2 <- intersect(cands, names(d))[1]
    if (is.na(c2)) rep(NA_character_, nrow(d)) else as.character(d[[c2]])
  }
  tibble(
    sample              = normalize_pq_id(d[[pt_col]]),
    `Tumour size`       = recode_tumour_size(pull_col(c("tumour_size","tumor_size","size"))),
    `Positive nodes`    = recode_nodes(pull_col(c("positive_nodes","node_count","nodes"))),
    LVSI                = recode_lvsi(pull_col(c("lvsi","lvi"))),
    `Depth of invasion` = recode_doi(pull_col(c("depth_of_invasion","doi","depth"))),
    HPV                 = recode_hpv(pull_col(c("hpv","hpv_status"))),
    p53                 = recode_p53(pull_col(c("p53","tp53")))
  ) %>% distinct(sample, .keep_all=TRUE)
}

char_denovo <- read_char_file(file_char_denovo)
char_recur  <- read_char_file(file_char_stage_recur)

align_char_cohort <- function(char_df, cols, col_name, levels_order) {
  v   <- setNames(char_df[[col_name]], char_df$sample)
  out <- unname(v[normalize_pq_id(cols)])
  out[is.na(out)] <- "Not recorded"
  factor(out, levels=levels_order)
}

read_age_file <- function(file) {
  d <- read_xlsx(file, .name_repair="minimal") %>% clean_names()
  id_col  <- intersect(names(d), c("patient_no","patient_number","patient","id",
                                   "randomisation_no","randomization_no","sample"))[1]
  val_col <- intersect(names(d), c("age","age_at_diagnosis","age_event",
                                   "age_at_diagnosis_2","age_diagnosis"))[1]
  if (is.na(id_col) || is.na(val_col)) { warning("No age columns found in: ", file); return(NULL) }
  tibble(sample = normalize_pq_id(d[[id_col]]),
         age    = suppressWarnings(as.numeric(as.character(d[[val_col]])))) %>%
    distinct(sample, .keep_all=TRUE)
}

age_denovo <- read_age_file(file_age_denovo)
age_recur  <- read_age_file(file_age_recur)

make_cohort_clinical <- function(cols, msi_src, tmb_src, stage_src, char_src, age_src) {
  msi_v <- to_numeric_safe(align_to_cols(msi_src, cols)); msi_v <- ifelse(is.na(msi_v)|!is.finite(msi_v), 0, msi_v)
  tmb_v <- to_numeric_safe(align_to_cols(tmb_src, cols)); tmb_v <- ifelse(is.na(tmb_v)|!is.finite(tmb_v), 0, tmb_v)
  ac <- function(name, lv) align_char_cohort(char_src, cols, name, lv)
  list(
    msi        = msi_v,
    tmb        = tmb_v,
    stage      = recode_stage_vec(align_to_cols(stage_src, cols)),
    tumoursize = ac("Tumour size",       c("<2 cm","2-4 cm",">4 cm","Not recorded")),
    doi        = ac("Depth of invasion", c("<2 mm","\u22652 mm","Not recorded")),
    hpv        = ac("HPV",               c("Negative","Positive","Not recorded")),
    p53        = ac("p53",               c("Wild-type","Abnormal","Not recorded")),
    lvsi       = ac("LVSI",              c("Absent","Present","Not recorded")),
    nodes      = ac("Positive nodes",    c("0","1","\u22652","Not recorded")),
    age        = age_src$age[match(normalize_pq_id(cols), age_src$sample)]
  )
}

clin_denovo <- make_cohort_clinical(denovo_cols,    msi_denovo, tmb_denovo, stage_denovo, char_denovo, age_denovo)
clin_recur  <- if (length(recurrent_cols)==0) NULL else
  make_cohort_clinical(recurrent_cols, msi_recur,  tmb_recur,  stage_recur,  char_recur,  age_recur)

msi_max <- 35; msi_thresh <- 30
tmb_max <- 45; tmb_thresh <- 10

stage_cols <- c("1A"="#E8F5E9","1B"="#E3F2FD","3A"="#F8BBD0","3B"="#F06292","3C"="#C2185B")
col_tumour_size <- c("<2 cm"="#BBDEFB","2-4 cm"="#42A5F5",">4 cm"="#1565C0","Not recorded"=na_col)
col_doi         <- c("<2 mm"="#EDE7F6","\u22652 mm"="#6A1B9A","Not recorded"=na_col)
col_lvsi        <- c("Absent"="#FFEBEE","Present"="#B71C1C","Not recorded"=na_col)
col_nodes       <- c("0"="#E8F5E9","1"="#66BB6A","\u22652"="#1B5E20","Not recorded"=na_col)
col_hpv         <- c("Negative"="#E0F7FA","Positive"="#00695C","Not recorded"=na_col)
col_p53         <- c("Wild-type"="#FFF3E0","Abnormal"="#E65100","Not recorded"=na_col)

make_bottom_anno <- function(clin) {
  age_v    <- clin$age
  age_plot <- ifelse(is.na(age_v)|!is.finite(age_v), 0, age_v)
  age_max  <- if (all(is.na(age_v))) 100 else ceiling(max(age_v, na.rm=TRUE)/10)*10
  if (!is.finite(age_max) || age_max<=0) age_max <- 100
  age_axis_at <- seq(0, age_max, by=20)
  age_thresh  <- if (all(is.na(age_v))) NA_real_ else mean(age_v, na.rm=TRUE)
  
  ha <- HeatmapAnnotation(
    MSI = anno_barplot(
      clin$msi, ylim=c(0,msi_max), gp=gpar(fill="#26A69A", col=NA),
      border=TRUE, border_gp=gpar(col="black", lwd=0.8), height=unit(2.2,"cm"),
      axis_param=list(at=c(0,10,20,30), labels=c("0","10","20","30"), side="left", gp=gpar(fontsize=ANN_FONT_SIZE))
    ),
    TMB = anno_barplot(
      clin$tmb, ylim=c(0,tmb_max), gp=gpar(fill="#D81B60", col=NA),
      border=TRUE, border_gp=gpar(col="black", lwd=0.8), height=unit(2.2,"cm"),
      axis_param=list(at=c(0,10,20,30,40), labels=c("0","10","20","30","40"), side="left", gp=gpar(fontsize=ANN_FONT_SIZE))
    ),
    Stage               = clin$stage,
    `Tumour size`       = clin$tumoursize,
    `Depth of invasion` = clin$doi,
    HPV                 = clin$hpv,
    p53                 = clin$p53,
    LVSI                = clin$lvsi,
    `Positive nodes`    = clin$nodes,
    Age = anno_points(
      age_plot, ylim=c(0,age_max), pch=16, size=unit(2,"mm"),
      gp=gpar(col="black", fill="black"), border=TRUE, border_gp=gpar(col="black", lwd=0.8),
      height=unit(2.2,"cm"),
      axis_param=list(at=age_axis_at, labels=as.character(age_axis_at), side="left", gp=gpar(fontsize=ANN_FONT_SIZE))
    ),
    spacer = anno_empty(border=FALSE, height=unit(6,"mm")),
    col = list(
      Stage=stage_cols, `Tumour size`=col_tumour_size, `Depth of invasion`=col_doi,
      HPV=col_hpv, p53=col_p53, LVSI=col_lvsi, `Positive nodes`=col_nodes
    ),
    which="column",
    gap=unit(c(rep(2,10), 0), "mm"),
    gp=gpar(col="white", lwd=0.8),
    border=FALSE,
    annotation_name_side="left",
    annotation_name_gp=gpar(fontsize=ANN_FONT_SIZE, fontface="plain"),
    show_annotation_name=c(
      MSI=TRUE, TMB=TRUE, Stage=TRUE, `Tumour size`=TRUE, `Depth of invasion`=TRUE,
      HPV=TRUE, p53=TRUE, LVSI=TRUE, `Positive nodes`=TRUE, Age=TRUE,
      spacer=FALSE
    ),
    simple_anno_size=unit(0.32,"cm"),
    annotation_height=unit.c(
      unit(2.2,"cm"), unit(2.2,"cm"),
      unit(0.32,"cm"), unit(0.32,"cm"), unit(0.32,"cm"),
      unit(0.32,"cm"), unit(0.32,"cm"), unit(0.32,"cm"), unit(0.32,"cm"),
      unit(2.2,"cm"), unit(6,"mm")
    )
  )
  list(anno=ha, age_thresh=age_thresh, age_max=age_max)
}

## =========================
## 9-11) Build and draw one oncoprint per cohort
## =========================
ROW_HEIGHT <- unit(5, "mm")

MASTER_SAMPLES  <- 39
MASTER_WIDTH_IN <- 22
LEGEND_WIDTH_IN <- 4.0
COL_WIDTH_IN    <- (MASTER_WIDTH_IN - LEGEND_WIDTH_IN) / MASTER_SAMPLES

ROWBAR_WIDTH_IN <- 1.6

## =========================
## Grouped bar helpers
## =========================
compute_actionable_counts <- function(elig_mat) {
  if (nrow(elig_mat)==0 || ncol(elig_mat)==0) {
    return(matrix(0, nrow=ncol(elig_mat), ncol=2,
                  dimnames=list(colnames(elig_mat), c("OncoKB","ESCAT"))))
  }
  onco_counts  <- apply(elig_mat, 2, function(col) sum(str_detect(col, "OncoKB")))
  escat_counts <- apply(elig_mat, 2, function(col) sum(str_detect(col, "ESCAT")))
  cbind(OncoKB = onco_counts, ESCAT = escat_counts)
}

compute_row_actionable_counts <- function(elig_mat) {
  if (nrow(elig_mat)==0 || ncol(elig_mat)==0) {
    return(matrix(0, nrow=0, ncol=2, dimnames=list(NULL, c("OncoKB","ESCAT"))))
  }
  onco_counts  <- apply(elig_mat, 1, function(row) sum(str_detect(row, "OncoKB")))
  escat_counts <- apply(elig_mat, 1, function(row) sum(str_detect(row, "ESCAT")))
  cbind(OncoKB = onco_counts, ESCAT = escat_counts)
}

build_top_bar <- function(elig_mat) {
  counts <- compute_actionable_counts(elig_mat)
  HeatmapAnnotation(
    `Actionable genes` = anno_barplot(
      counts,
      beside     = TRUE,
      bar_width  = 0.8,
      gp         = gpar(fill = c(ALT_COL, ESCAT_COL), col = NA),
      border     = TRUE,
      height     = unit(2, "cm"),
      axis_param = list(side = "left", labels_rot = 0)
    ),
    spacer = anno_empty(border=FALSE, height=unit(5,"mm")),
    show_legend=FALSE,
    show_annotation_name=c(`Actionable genes`=TRUE, spacer=FALSE),
    annotation_name_side="left",
    annotation_name_rot=90,
    annotation_name_gp=gpar(fontsize=ANN_FONT_SIZE, fontface="plain")
  )
}

mut_legend <- Legend(
  title="Database", at=c("OncoKB","ESCAT"), labels=c("OncoKB","ESCAT"),
  legend_gp=gpar(fill=c(ALT_COL, ESCAT_COL)), ncol=1
)

draw_cohort_oncoprint <- function(elig_mat, clin, cohort_title, out_pdf, out_png) {
  if (is.null(clin) || nrow(elig_mat)==0 || ncol(elig_mat)==0) {
    warning("No actionable genes for cohort '", cohort_title, "' - skipping its oncoprint.")
    return(invisible(NULL))
  }
  
  ba <- make_bottom_anno(clin)
  
  op <- oncoPrint(
    elig_mat,
    alter_fun               = alter_fun,
    alter_fun_is_vectorized = FALSE,
    get_type                = get_type,
    col                     = col,
    row_gap                 = unit(0.8,"mm"),
    column_title            = NULL,
    top_annotation          = build_top_bar(elig_mat),
    bottom_annotation       = ba$anno,
    right_annotation        = rowAnnotation(
      `Patients` = anno_barplot(
        compute_row_actionable_counts(elig_mat),
        which      = "row",
        beside     = TRUE,
        bar_width  = 0.8,
        gp         = gpar(fill = c(ALT_COL, ESCAT_COL), col = NA),
        border     = TRUE,
        width      = unit(3, "cm"),
        axis_param = list(side = "top", labels_rot = 0)
      ),
      annotation_name_side = "top",
      annotation_name_gp   = gpar(fontsize=ANN_FONT_SIZE, fontface="plain")
    ),
    row_names_gp            = gpar(fontsize=ROW_FONT_SIZE),
    row_names_side          = "right",
    column_names_gp         = gpar(fontsize=SAMPLE_FONT_SIZE),
    show_pct                = FALSE,
    show_row_names          = TRUE,
    remove_empty_rows       = TRUE,
    remove_empty_columns    = FALSE,
    show_column_names       = TRUE,
    height                  = ROW_HEIGHT * nrow(elig_mat)
  )
  
  pdf_w <- COL_WIDTH_IN * ncol(elig_mat) + LEGEND_WIDTH_IN + ROWBAR_WIDTH_IN
  pdf_h <- max(14, nrow(elig_mat)*0.20 + 8)
  
  draw_one <- function() {
    draw(
      op,
      heatmap_legend_side    = "right",
      annotation_legend_side = "right",
      annotation_legend_list = list(mut_legend),
      padding                = unit(c(2, 8, 2, 2), "mm")
    )
    decorate_annotation("Actionable genes", {
      grid.lines(x=unit(c(0,1),"npc"), y=unit(c(0,0),"native"),
                 gp=gpar(col="black", lwd=1.5))
    })
    decorate_annotation("MSI", {
      grid.lines(x=unit(c(0,1),"npc"), y=unit(msi_thresh/msi_max,"npc"),
                 gp=gpar(col="black", lty=2, lwd=1.2))
    })
    decorate_annotation("TMB", {
      grid.lines(x=unit(c(0,1),"npc"), y=unit(tmb_thresh/tmb_max,"npc"),
                 gp=gpar(col="black", lty=2, lwd=1.2))
    })
    if (is.finite(ba$age_thresh)) {
      decorate_annotation("Age", {
        grid.lines(x=unit(c(0,1),"npc"), y=unit(ba$age_thresh/ba$age_max,"npc"),
                   gp=gpar(col="black", lty=2, lwd=1.2))
      })
    }
  }
  
  while (dev.cur() > 1) dev.off()
  pdf(out_pdf, width=pdf_w, height=pdf_h)
  tryCatch(draw_one(), error=function(e) cat("Draw error (", cohort_title, "):", conditionMessage(e), "\n"))
  dev.off()
  cat("Oncoprint PDF (", cohort_title, "):\n ", out_pdf, "\n")
  
  while (dev.cur() > 1) dev.off()
  png(out_png, width=pdf_w*150, height=pdf_h*150, res=150)
  tryCatch(draw_one(), error=function(e) cat("Draw error (", cohort_title, "):", conditionMessage(e), "\n"))
  dev.off()
  cat("Oncoprint PNG (", cohort_title, "):\n ", out_png, "\n")
  
  invisible(op)
}

out_pdf_denovo <- file.path(out_dir, paste0("OncoKB_Oncoprint_DeNovo_",     ts, ".pdf"))
out_png_denovo <- file.path(out_dir, paste0("OncoKB_Oncoprint_DeNovo_",     ts, ".png"))
out_pdf_recur  <- file.path(out_dir, paste0("OncoKB_Oncoprint_Recurrence_", ts, ".pdf"))
out_png_recur  <- file.path(out_dir, paste0("OncoKB_Oncoprint_Recurrence_", ts, ".png"))

draw_cohort_oncoprint(elig_mat_denovo, clin_denovo, "De novo",    out_pdf_denovo, out_png_denovo)
draw_cohort_oncoprint(elig_mat_recur,  clin_recur,  "Recurrence", out_pdf_recur,  out_png_recur)

###############################################################################
## PART B - ACTIONABLE GENE SUMMARY (Excel) + SOLID TUMOUR BAR CHART
###############################################################################

## =========================
## 12) Clean & filter (actionable gene summary)
## =========================
df_gene <- df %>%
  mutate(across(all_of(sample_cols), ~suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    gene_clean = toupper(trimws(as.character(gene))),
    broad_cat  = trimws(as.character(broad_category)),
    lc_str     = trimws(as.character(level_classification)),
    level_str  = trimws(as.character(level)),
    level_num  = suppressWarnings(as.numeric(trimws(as.character(level)))),
    mut_str    = trimws(as.character(mutation)),
    therapy_clean = trimws(as.character(therapy)),
    db_clean   = classify_database(database)
  ) %>%
  mutate(
    broad_cat  = if_else(is.na(broad_cat)|broad_cat %in% c("","NA"), "Unknown", broad_cat),
    broad_cat  = if_else(broad_cat=="Pan-Cancer", "Solid Tumour", broad_cat),
    is_tmb     = str_detect(gene_clean, regex("tumou?r.mutation.burden|\\bTMB\\b", ignore_case=TRUE)),
    is_wt_req  = str_detect(mut_str, regex("wild.?type|\\bWT\\b", ignore_case=TRUE)),
    therapy_clean = if_else(is.na(therapy_clean)|therapy_clean=="", "Unknown Therapy", therapy_clean),
    gene_clean = case_when(
      is_tmb ~ "TMB",
      is.na(gene_clean)|gene_clean=="" ~ "UNKNOWN",
      TRUE ~ gene_clean
    )
  ) %>%
  filter(!str_detect(lc_str, regex("resist", ignore_case=TRUE))) %>%
  filter(gene_clean != "TMB") %>%
  filter(!is_wt_req) %>%
  filter(!is.na(db_clean)) %>%
  filter(
    (db_clean=="OncoKB" & !is.na(level_num) & level_num %in% c(1,2,3)) |
      (db_clean=="ESCAT")
  )

has_elig_gene <- rowSums(df_gene[, sample_cols, drop=FALSE]==2, na.rm=TRUE) > 0
df_gene <- df_gene[has_elig_gene, ]
if (nrow(df_gene)==0) stop("No eligible/actionable rows after filtering (actionable gene summary).")

## =========================
## 13) Collapse to unique gene x patient actionable calls, per database
## =========================
gene_patient_hits <- df_gene %>%
  select(gene_clean, broad_cat, db_clean, all_of(sample_cols)) %>%
  pivot_longer(cols=all_of(sample_cols), names_to="patient", values_to="value") %>%
  filter(value == 2) %>%
  distinct(gene_clean, broad_cat, db_clean, patient)

## =========================
## 13b) Sample-level detail: one row per sample per targetable drug
## =========================
sample_drug_hits <- df_gene %>%
  select(gene_clean, broad_cat, therapy_clean, level_str, db_clean, all_of(sample_cols)) %>%
  pivot_longer(cols=all_of(sample_cols), names_to="patient", values_to="value") %>%
  filter(value == 2) %>%
  distinct(patient, gene_clean, therapy_clean, broad_cat, level_str, db_clean) %>%
  arrange(patient, gene_clean, therapy_clean) %>%
  transmute(
    Patient           = patient,
    Gene              = gene_clean,
    `Targetable Drug` = therapy_clean,
    `Broad Category`  = broad_cat,
    Database          = db_clean,
    `Level`           = level_str
  )

## =========================
## 14) Summarise: n actionable patients + % frequency, per gene, PER DATABASE
## =========================
gene_broad_lookup <- gene_patient_hits %>%
  count(gene_clean, broad_cat, sort=TRUE) %>%
  group_by(gene_clean) %>%
  slice_max(n, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  select(gene_clean, broad_cat)

gene_summary <- gene_patient_hits %>%
  distinct(gene_clean, db_clean, patient) %>%
  count(gene_clean, db_clean, name="n_actionable") %>%
  left_join(gene_broad_lookup, by="gene_clean") %>%
  mutate(
    total_patients = n_total_patients,
    pct_frequency  = round(100 * n_actionable / n_total_patients, 1),
    broad_sort     = broad_sort_order[broad_cat],
    broad_sort     = if_else(is.na(broad_sort), 99, broad_sort)
  ) %>%
  arrange(broad_sort, desc(pct_frequency), gene_clean) %>%
  transmute(
    Gene              = gene_clean,
    Database          = db_clean,
    `Broad Category`  = broad_cat,
    `n Actionable`    = n_actionable,
    `Total Patients`  = total_patients,
    `% Frequency`     = pct_frequency
  )

## =========================
## 15) Summary counts (overall + per broad category), PER DATABASE
## =========================
overall_summary <- tibble(
  `Metric` = c(
    "Total actionable genes - OncoKB (Level 1-3, WT-required excluded)",
    "Total actionable genes - ESCAT (all levels, WT-required excluded)",
    "Total patients screened"
  ),
  `Value`  = c(
    n_distinct(gene_summary$Gene[gene_summary$Database=="OncoKB"]),
    n_distinct(gene_summary$Gene[gene_summary$Database=="ESCAT"]),
    n_total_patients
  )
)

by_category_summary <- gene_summary %>%
  group_by(`Broad Category`, Database) %>%
  summarise(`n Actionable Genes` = n_distinct(Gene), .groups="drop") %>%
  arrange(desc(`n Actionable Genes`))

## =========================
## 17) Solid Tumour-only sheet, PER DATABASE
## =========================
solid_tumour_summary <- gene_patient_hits %>%
  filter(broad_cat == "Solid Tumour") %>%
  distinct(gene_clean, db_clean, patient) %>%
  count(gene_clean, db_clean, name="n_actionable") %>%
  mutate(
    total_patients = n_total_patients,
    pct_frequency  = round(100 * n_actionable / n_total_patients, 1)
  ) %>%
  arrange(desc(pct_frequency), gene_clean) %>%
  transmute(
    Gene             = gene_clean,
    Database         = db_clean,
    `n Actionable`   = n_actionable,
    `Total Patients` = total_patients,
    `% Frequency`    = pct_frequency
  )

## =========================
## 18) Write workbook
## =========================
write_xlsx(
  list(
    "Summary"                       = overall_summary,
    "By Broad Category"             = by_category_summary,
    "All Actionable Genes"          = gene_summary,
    "Solid Tumour Genes"            = solid_tumour_summary,
    "Sample-Level Actionable Drugs" = sample_drug_hits
  ),
  path = out_xlsx
)
cat("Excel summary written to:\n", out_xlsx, "\n")

## =========================
## 19) Bar chart: solid tumour actionable genes and TMB-high, by stage
## =========================================================================

if (!file.exists(staging_file)) {
  stop("Staging file not found at:\n", staging_file)
}

meta <- read_xlsx(staging_file, sheet=1) %>% clean_names()

if (!meta_id_col %in% names(meta) || !meta_stage_col %in% names(meta)) {
  stop(
    "`meta_id_col` ('", meta_id_col, "') or `meta_stage_col` ('", meta_stage_col,
    "') not found in staging file.\n",
    "Columns available: ", paste(names(meta), collapse=", ")
  )
}

meta_clean <- meta %>%
  transmute(
    id_raw     = as.character(.data[[meta_id_col]]),
    stage_raw  = trimws(as.character(.data[[meta_stage_col]]))
  ) %>%
  mutate(
    id_num  = suppressWarnings(as.integer(str_extract(id_raw, "\\d+"))),
    patient = if_else(is.na(id_num), NA_character_, paste0("pq", id_num))
  ) %>%
  filter(!is.na(patient)) %>%
  mutate(
    stage_group = str_extract(toupper(str_remove_all(stage_raw, "\\s")), "^(1A|1B|3A|3B|3C|3)"),
    risk_class = case_when(
      stage_group %in% c("1A","1B")          ~ "Early-stage (node-negative)",
      stage_group %in% c("3","3A","3B","3C") ~ "Late-stage (node-positive)",
      TRUE ~ NA_character_
    )
  )

unmapped <- meta_clean %>% filter(is.na(risk_class)) %>% distinct(stage_raw)
if (nrow(unmapped) > 0) {
  cat("NOTE: the following stage/node values did not map to a risk class ",
      "and were excluded from the bar chart:\n")
  print(unmapped)
}

meta_clean <- meta_clean %>% filter(!is.na(risk_class))

meta_clean <- meta_clean %>% filter(patient %in% sample_cols)
n_by_risk  <- meta_clean %>% count(risk_class, name="n_patients")

solid_hits <- gene_patient_hits %>%
  filter(broad_cat == "Solid Tumour") %>%
  inner_join(meta_clean, by="patient") %>%
  distinct(gene_clean, patient, risk_class, db_clean)

if (nrow(solid_hits) == 0) stop("No Solid Tumour actionable hits with a mapped risk class - check the metadata join.")

df_tmb <- df %>%
  mutate(across(all_of(sample_cols), ~suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    gene_str  = toupper(trimws(as.character(gene))),
    broad_cat = trimws(as.character(broad_category)),
    lc_str    = trimws(as.character(level_classification)),
    level_num = suppressWarnings(as.numeric(trimws(as.character(level)))),
    is_tmb    = str_detect(gene_str, regex("tumou?r.mutation.burden|\\bTMB\\b", ignore_case=TRUE)),
    db_clean  = classify_database(database)
  ) %>%
  mutate(
    broad_cat = if_else(is.na(broad_cat)|broad_cat %in% c("","NA"), "Unknown", broad_cat),
    broad_cat = if_else(broad_cat=="Pan-Cancer", "Solid Tumour", broad_cat)
  ) %>%
  filter(is_tmb) %>%
  filter(!str_detect(lc_str, regex("resist", ignore_case=TRUE))) %>%
  filter(!is.na(db_clean)) %>%
  filter(
    (db_clean=="OncoKB" & !is.na(level_num) & level_num %in% c(1,2,3)) |
      (db_clean=="ESCAT")
  ) %>%
  filter(broad_cat == "Solid Tumour")

tmb_patient_hits <- df_tmb %>%
  select(db_clean, all_of(sample_cols)) %>%
  pivot_longer(cols=all_of(sample_cols), names_to="patient", values_to="value") %>%
  filter(value == 2) %>%
  distinct(patient, db_clean)

tmb_hits <- tmb_patient_hits %>%
  inner_join(meta_clean, by="patient") %>%
  distinct(patient, risk_class, db_clean)

tmb_bar_row <- tmb_hits %>%
  count(risk_class, db_clean, name="n_actionable") %>%
  right_join(
    expand.grid(risk_class = names(risk_pal), db_clean = c("OncoKB","ESCAT"),
                stringsAsFactors = FALSE),
    by = c("risk_class","db_clean")
  ) %>%
  mutate(n_actionable = if_else(is.na(n_actionable), 0L, n_actionable)) %>%
  left_join(n_by_risk, by="risk_class") %>%
  mutate(
    pct        = round(100 * n_actionable / n_patients, 0),
    gene_clean = "TMB-high"
  ) %>%
  select(gene_clean, risk_class, db_clean, n_actionable, n_patients, pct)

bar_data <- solid_hits %>%
  count(gene_clean, risk_class, db_clean, name="n_actionable") %>%
  left_join(n_by_risk, by="risk_class") %>%
  mutate(pct = round(100 * n_actionable / n_patients, 0)) %>%
  complete(gene_clean, risk_class = names(risk_pal), db_clean = c("OncoKB","ESCAT"),
           fill=list(n_actionable=0, pct=0)) %>%
  bind_rows(tmb_bar_row)

gene_order <- bar_data %>%
  group_by(gene_clean) %>%
  summarise(tot = sum(pct), .groups="drop") %>%
  arrange(desc(tot)) %>%
  pull(gene_clean)

bar_data <- bar_data %>%
  mutate(
    gene_clean = factor(gene_clean, levels=gene_order),
    risk_class = factor(risk_class, levels=names(risk_pal)),
    db_clean   = factor(db_clean,   levels=c("OncoKB","ESCAT"))
  )

p_bar <- ggplot(bar_data, aes(x=gene_clean, y=pct, fill=db_clean)) +
  geom_col(position=position_dodge(width=0.8), width=0.75, colour=NA) +
  geom_text(
    aes(label=paste0(pct,"%"), colour=db_clean, group=db_clean),
    position=position_dodge(width=0.8), vjust=-0.4, size=2.6,
    fontface="bold", show.legend=FALSE
  ) +
  scale_fill_manual(values=db_pal, name="Database") +
  scale_colour_manual(values=db_pal, guide="none") +
  scale_y_continuous(limits=c(0, max(bar_data$pct, na.rm=TRUE)*1.15), expand=c(0,0)) +
  facet_wrap(~risk_class, ncol=1, scales="free_x") +
  labs(x=NULL, y="Percentage (%)") +
  theme_classic(base_size=10) +
  theme(
    axis.text.x     = element_text(angle=45, hjust=1, face="italic"),
    legend.position = "right",
    legend.title    = element_text(face="bold"),
    strip.text      = element_text(face="bold")
  )

ggsave(out_bar_pdf, p_bar, width=max(7, n_distinct(bar_data$gene_clean)*0.55+2), height=9)
ggsave(out_bar_png, p_bar, width=max(7, n_distinct(bar_data$gene_clean)*0.55+2), height=9, dpi=300)

cat("Bar chart written to:\n", out_bar_pdf, "\n", out_bar_png, "\n")

cat("\nDone!\n")