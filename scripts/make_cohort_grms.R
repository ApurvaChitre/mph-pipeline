#!/usr/bin/env Rscript

# make_cohort_grms.R
# ------------------------------------------------------------
# Create MPH custom GRMs for multiple cohorts/traits by masking
# a full cross-cohort GRM (gg) into:
#   - N single-cohort GRMs  (cohort1 ... cohortN)
#   - N*(N-1)/2 cross-cohort GRMs (cohort1_cohort2, ...)
#   - (N-1) diagonal residual matrices (cohort1_diagonal ... cohort(N-1)_diagonal)
#
# This script is intended to be run inside a MPH run directory that
# contains:
#   pheno/indicator_covariates.csv
#   grm/gg.grm.*
#
# It can be run from anywhere as long as you pass --workdir.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

as_int <- function(x, default = NA_integer_) {
  if (is.na(x) || x == "") return(default)
  suppressWarnings(as.integer(x))
}

workdir      <- get_arg("--workdir", ".")
cov_file     <- get_arg("--cov_file", file.path(workdir, "pheno", "indicator_covariates.csv"))
gg_prefix    <- get_arg("--gg_prefix", file.path(workdir, "grm", "gg"))
outdir       <- get_arg("--outdir", file.path(workdir, "grm"))
mph_functs   <- get_arg("--mph_functs", Sys.getenv("MPH_FUNCTS_R", "/tscc/projects/ps-palmer/software/local/src/mph/mph_functs.R"))
n_cohorts_in <- as_int(get_arg("--n_cohorts", NA_character_), default = NA_integer_)

# Validate paths
workdir <- normalizePath(workdir, mustWork = TRUE)

if (!file.exists(cov_file)) stop("Missing covariate file: ", cov_file)
if (!file.exists(paste0(gg_prefix, ".grm.bin"))) stop("Missing gg GRM: ", paste0(gg_prefix, ".grm.bin"))
if (!file.exists(mph_functs)) stop("Missing mph_functs.R: ", mph_functs)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Read inputs
cohort_data <- read.csv(cov_file, stringsAsFactors = FALSE)

# Infer number of cohorts if not provided: count columns named cohort*
cohort_cols <- grep("^cohort\\d+$", names(cohort_data), value = TRUE)
if (length(cohort_cols) == 0) {
  stop("No cohort columns found in cov_file. Expected columns like cohort1, cohort2, ...")
}

n_cohorts <- n_cohorts_in
if (is.na(n_cohorts)) n_cohorts <- length(cohort_cols)

cohorts <- paste0("cohort", 1:n_cohorts)

missing_cols <- setdiff(cohorts, names(cohort_data))
if (length(missing_cols) > 0) {
  stop(
    "Covariate file is missing expected cohort columns: ",
    paste(missing_cols, collapse = ", "),
    "\nFound cohort-like columns: ",
    paste(cohort_cols, collapse = ", ")
  )
}

# Load MPH helper I/O
source(mph_functs)

# Read the GRM and force symmetry
matrix_data <- read_grm(gg_prefix)
matrix_data <- matrix_data + t(matrix_data)
diag(matrix_data) <- diag(matrix_data) / 2

# Generate cohort combinations
cohort_combinations <- combn(cohorts, 2, simplify = FALSE)

# Mask helpers
process_matrix <- function(matrix, cohort_data, cohort1, cohort2 = NULL) {
  result <- matrix

  if (is.null(cohort2)) {
    cohort1_ids <- cohort_data$ID[cohort_data[[cohort1]] == 1]
    mask <- outer(rownames(result) %in% cohort1_ids, colnames(result) %in% cohort1_ids)
    result[!mask] <- 0
  } else {
    cohort1_ids <- cohort_data$ID[cohort_data[[cohort1]] == 1]
    cohort2_ids <- cohort_data$ID[cohort_data[[cohort2]] == 1]
    mask1 <- outer(rownames(result) %in% cohort1_ids, colnames(result) %in% cohort2_ids)
    mask2 <- outer(rownames(result) %in% cohort2_ids, colnames(result) %in% cohort1_ids)
    result[!(mask1 | mask2)] <- 0
  }

  return(result)
}

create_diagonal_matrix <- function(matrix, cohort_data, cohort) {
  result <- matrix * 0
  cohort_ids <- cohort_data$ID[cohort_data[[cohort]] == 1]
  diag_indices <- which(rownames(result) %in% cohort_ids)
  result[cbind(diag_indices, diag_indices)] <- 1
  return(result)
}

# Case 1: Single cohort matrices
for (cohort in cohorts) {
  result <- process_matrix(matrix_data, cohort_data, cohort)
  write_grm(file.path(outdir, cohort), colnames(result), result)
}

# Case 2: Cohort pair matrices
for (combo in cohort_combinations) {
  combo_name <- paste(combo, collapse = "_")
  result <- process_matrix(matrix_data, cohort_data, combo[1], combo[2])
  write_grm(file.path(outdir, combo_name), colnames(result), result)
}

# Case 3: Diagonal matrices for each cohort (all but last; MPH includes final one automatically)
if (length(cohorts) >= 2) {
  for (cohort in cohorts[-length(cohorts)]) {
    diagonal_name <- paste0(cohort, "_diagonal")
    result <- create_diagonal_matrix(matrix_data, cohort_data, cohort)
    write_grm(file.path(outdir, diagonal_name), colnames(result), result)
  }
}

cat("Wrote custom GRMs to: ", normalizePath(outdir), "\n", sep = "")

