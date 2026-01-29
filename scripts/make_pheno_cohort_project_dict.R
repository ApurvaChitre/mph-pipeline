#!/usr/bin/env Rscript
# make_pheno_cohort_project_dict.R
# Sequential cohort numbering.
# include/exclude traits can be:
#   --include_traits=trait1,trait2,...
#   --exclude_traits=traitA,traitB,...
# OR a file path prefixed with @:
#   --include_traits=@/path/to/traits.txt
#   --exclude_traits=@/path/to/traits.txt

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

normalize_trait <- function(x) {
  x <- str_trim(x)
  x <- x[x != ""]
  str_remove(x, "^regressedlr_")
}

read_trait_spec <- function(spec) {
  # Returns character vector of traits from:
  # - empty/NA -> character(0)
  # - "@/path/file.txt" -> one per line
  # - "a,b,c" -> comma-separated
  if (is.na(spec) || spec == "") return(character(0))
  
  # File mode if starts with "@"
  if (startsWith(spec, "@")) {
    path <- sub("^@", "", spec)
    if (!file.exists(path)) stop("Trait file not found: ", path)
    x <- readLines(path, warn = FALSE)
    x <- str_trim(x)
    x <- x[x != ""]
    return(x)
  }
  
  # Comma-separated mode
  x <- unlist(str_split(spec, ","))
  x <- str_trim(x)
  x <- x[x != ""]
  x
}

infer_project_name <- function(trait) {
  toks <- unlist(str_split(trait, "_"))
  if (length(toks) == 0) return(NA_character_)
  
  if (toks[1] == "u01") {
    return(paste(toks[1:min(3, length(toks))], collapse = "_"))
  }
  if (toks[1] == "p50") {
    if (length(toks) >= 4 && str_detect(toks[4], "^20\\d\\d$")) {
      return(paste(toks[1:4], collapse = "_"))
    }
    if (length(toks) >= 4 && toks[4] == "open") {
      return(paste(toks[1:4], collapse = "_"))
    }
    return(paste(toks[1:min(3, length(toks))], collapse = "_"))
  }
  
  paste(toks[1:min(3, length(toks))], collapse = "_")
}

processed_file <- get_arg("--processed")
out_dir        <- get_arg("--outdir", ".")
out_name       <- get_arg("--outname", "pheno_cohort_project_dict.csv")
include_spec   <- get_arg("--include_traits", "")
exclude_spec   <- get_arg("--exclude_traits", "")

if (is.na(processed_file) || processed_file == "") {
  stop("Missing --processed=/path/to/processed_data_ready.csv")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, out_name)

pheno <- read_csv(processed_file, col_types = cols(.default = col_guess(), rfid = col_character()))
if (!("rfid" %in% names(pheno))) stop("Processed file must contain column 'rfid': ", processed_file)

reg_cols <- names(pheno)[str_detect(names(pheno), "^regressedlr_")]
if (length(reg_cols) == 0) stop("No regressedlr_ columns found in processed file: ", processed_file)

all_traits <- str_remove(reg_cols, "^regressedlr_")

include_traits <- normalize_trait(read_trait_spec(include_spec))
exclude_traits <- normalize_trait(read_trait_spec(exclude_spec))

# Default: all traits in file order
traits_keep <- all_traits

# If include list given: use that order (then filter to existing)
if (length(include_traits) > 0) {
  traits_keep <- include_traits
}

# Apply exclusion
if (length(exclude_traits) > 0) {
  traits_keep <- setdiff(traits_keep, exclude_traits)
}

# Keep only traits that exist in this file
traits_keep <- traits_keep[traits_keep %in% all_traits]
if (length(traits_keep) == 0) stop("After filtering, no traits remain.")

# If not using include list, keep file column order
if (length(include_traits) == 0) {
  traits_keep <- all_traits[all_traits %in% traits_keep]
}

reg_cols_keep <- paste0("regressedlr_", traits_keep)

# Non-NA counts per regressedlr column
N_vec <- vapply(reg_cols_keep, function(cc) sum(!is.na(pheno[[cc]])), integer(1))

dict <- tibble(
  pheno_file   = paste0(reg_cols_keep, ".txt"),
  cohort_num   = paste0("cohort", seq_along(traits_keep)),
  project_name = vapply(traits_keep, infer_project_name, character(1)),
  N            = as.integer(N_vec)
)

write_csv(dict, out_file)
cat("Wrote:", out_file, "\n")
cat("Traits included:", nrow(dict), "\n")
cat("Total N (sum across traits):", sum(dict$N, na.rm = TRUE), "\n")

# Examples:
# Rscript make_pheno_cohort_project_dict.R --processed=... --outdir=... \
#   --include_traits=u01_tom_jhou_session1_locomotor1,p50_hao_chen_open_field_first_15
#
# Rscript make_pheno_cohort_project_dict.R --processed=... --outdir=... \
#   --exclude_traits=p50_shelly_flagel_2014_hab1_total_distance,u01_peter_kalivas_oft_distance_1_first_15
#
# File mode:
#   --include_traits=@/path/to/traits_keep.txt
#   --exclude_traits=@/path/to/traits_skip.txt
