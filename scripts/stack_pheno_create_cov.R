#!/usr/bin/env Rscript
# stack_pheno_create_cov.R
#
# Two supported input modes:
#   1) GCTA per-trait files (legacy): --pheno_dir=/path/to/dir
#      Each file has 3 columns: ID, FID, trait_value
#
#   2) Wide processed CSV (recommended): --processed_file=/path/to/processed_data_ready.csv
#      Must include an ID column (default: rfid) and regressedlr_* trait columns.
#
# In both modes, trait/cohort definitions are taken from:
#   --pheno_cohort_file (pheno_cohort_project_dict.csv)

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}
has_flag <- function(flag) any(args == flag)

print_usage <- function() {
  cat(
"Usage:
  # Mode A (recommended): stack from processed_data_ready.csv
  Rscript stack_pheno_create_cov.R \
    --processed_file=/path/to/processed_data_ready.csv \
    --pheno_cohort_file=/path/to/pheno_cohort_project_dict.csv \
    --output_pheno_file=/path/to/combined_phenotype.csv \
    --output_cov_file=/path/to/indicator_covariates.csv

  # Mode B (legacy): stack from per-trait GCTA phenotype files
  Rscript stack_pheno_create_cov.R \
    --pheno_dir=/path/to/pheno_dir \
    --pheno_cohort_file=/path/to/pheno_cohort_project_dict.csv \
    --output_pheno_file=/path/to/combined_phenotype.csv \
    --output_cov_file=/path/to/indicator_covariates.csv

Arguments:
  --processed_file
      Wide CSV with columns: rfid (or --id_col) and regressedlr_* trait columns.
      If provided, this mode is used.

  --id_col
      ID column name in processed file (default: rfid)

  --pheno_dir
      Directory with per-trait GCTA phenotype files (ID, FID, trait_value).
      Used only if --processed_file is not provided.

  --pheno_cohort_file
      CSV with columns: pheno_file, cohort_num, project_name, N

  --output_pheno_file
      Output stacked phenotype CSV: ID, trait_value

  --output_cov_file
      Output incidence-matrix covariate CSV: ID, cohort1, cohort2, ...

Options:
  --help
      Print this message and exit.
"
  )
}

if (has_flag("--help") || length(args) == 0) {
  print_usage()
  quit(status = 0)
}

processed_file    <- get_arg("--processed_file", "")
id_col            <- get_arg("--id_col", "rfid")
pheno_dir         <- get_arg("--pheno_dir", "")
pheno_cohort_file <- get_arg("--pheno_cohort_file")
output_pheno_file <- get_arg("--output_pheno_file")
output_cov_file   <- get_arg("--output_cov_file")

missing <- c(
  if (is.na(pheno_cohort_file) || pheno_cohort_file == "") "--pheno_cohort_file" else NULL,
  if (is.na(output_pheno_file) || output_pheno_file == "") "--output_pheno_file" else NULL,
  if (is.na(output_cov_file)   || output_cov_file   == "") "--output_cov_file" else NULL
)

# Must provide either processed_file or pheno_dir
if (is.na(processed_file)) processed_file <- ""
if (is.na(pheno_dir)) pheno_dir <- ""

if (processed_file == "" && pheno_dir == "") {
  missing <- c(missing, "--processed_file (or --pheno_dir)")
}

if (length(missing) > 0) {
  cat("ERROR: Missing required arguments:", paste(missing, collapse = ", "), "\n\n")
  print_usage()
  quit(status = 1)
}

if (!file.exists(pheno_cohort_file)) stop("pheno_cohort_file does not exist: ", pheno_cohort_file)

dir.create(dirname(output_pheno_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(output_cov_file),   showWarnings = FALSE, recursive = TRUE)

# Read cohort info (defines cohort order + which traits to stack)
pheno_cohort_info <- fread(pheno_cohort_file)
if (!all(c("pheno_file", "cohort_num") %in% names(pheno_cohort_info))) {
  stop("pheno_cohort_file must contain columns: pheno_file, cohort_num")
}

# Initialize outputs
combined_pheno <- data.table(ID = character(), trait_value = numeric())
combined_cov   <- data.table(ID = character())  # start with IDs only; add cohort cols iteratively

# ---------------------------------------------------------------------
# Mode A: processed file (wide)
# ---------------------------------------------------------------------
if (processed_file != "") {
  if (!file.exists(processed_file)) stop("processed_file does not exist: ", processed_file)

  message("[stack] Mode=processed (wide CSV): ", processed_file)
  dt <- fread(processed_file)

  if (!(id_col %in% names(dt))) {
    stop("processed_file is missing id_col '", id_col, "'. Available columns include: ",
         paste(head(names(dt), 20), collapse = ", "))
  }

  # Ensure ID is character
  dt[, (id_col) := as.character(get(id_col))]

  for (i in 1:nrow(pheno_cohort_info)) {
    cohort <- pheno_cohort_info$cohort_num[i]

    # mapping has pheno_file like: regressedlr_<trait>.txt
    trait_col <- sub("\\.txt$", "", pheno_cohort_info$pheno_file[i])

    if (!(trait_col %in% names(dt))) {
      stop("processed_file missing trait column: ", trait_col)
    }

    ph <- dt[!is.na(get(trait_col)), .(ID = get(id_col), trait_value = as.numeric(get(trait_col)))]

    # Stack
    combined_pheno <- rbind(combined_pheno, ph, use.names = TRUE)

    # Incidence covariate for this cohort
    cohort_cov <- unique(ph[, .(ID)])
    cohort_cov[, (cohort) := 1L]
    combined_cov <- merge(combined_cov, cohort_cov, by = "ID", all = TRUE)
  }
}

# ---------------------------------------------------------------------
# Mode B: per-trait GCTA phenotype files
# ---------------------------------------------------------------------
if (processed_file == "" && pheno_dir != "") {
  if (!dir.exists(pheno_dir)) stop("pheno_dir does not exist: ", pheno_dir)
  message("[stack] Mode=gcta (per-trait files): ", pheno_dir)

  for (i in 1:nrow(pheno_cohort_info)) {
    cohort <- pheno_cohort_info$cohort_num[i]
    pheno_file <- file.path(pheno_dir, pheno_cohort_info$pheno_file[i])
    if (!file.exists(pheno_file)) stop("Missing phenotype file: ", pheno_file)

    # Most of our files are headerless; we force column names.
    pheno_data <- fread(pheno_file, col.names = c("ID", "FID", "trait_value"))
    pheno_data[, ID := as.character(ID)]
    pheno_data <- pheno_data[!is.na(trait_value)]

    combined_pheno <- rbind(combined_pheno, pheno_data[, .(ID, trait_value)], use.names = TRUE)

    cohort_cov <- unique(pheno_data[, .(ID)])
    cohort_cov[, (cohort) := 1L]
    combined_cov <- merge(combined_cov, cohort_cov, by = "ID", all = TRUE)
  }
}

# Fill NA -> 0 for cohort indicator columns
for (cc in setdiff(names(combined_cov), "ID")) {
  set(combined_cov, which(is.na(combined_cov[[cc]])), cc, 0L)
}

# Order columns: ID then cohort1..cohortK (numeric order)
cohort_cols <- grep("^cohort", names(combined_cov), value = TRUE)
cohort_cols <- cohort_cols[order(as.integer(gsub("^cohort", "", cohort_cols)))]
setcolorder(combined_cov, c("ID", cohort_cols))

# Write files
fwrite(combined_pheno, output_pheno_file, sep = ",", col.names = TRUE)
fwrite(combined_cov,   output_cov_file,   sep = ",", col.names = TRUE)

cat("Wrote stacked phenotype file: ", output_pheno_file, "\n", sep = "")
cat("Wrote incidence covariate file: ", output_cov_file, "\n", sep = "")
cat("Rows in combined_pheno: ", nrow(combined_pheno), "\n", sep = "")
cat("Unique IDs in covariates: ", nrow(combined_cov), "\n", sep = "")

rm(list = ls())
