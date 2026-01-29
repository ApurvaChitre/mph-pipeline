#!/usr/bin/env Rscript

# Compare GenomicSEM-ready MPH outputs (.RData) between two directories.
#
# This is meant as a regression test: if the matrices (and their dimnames)
# differ between an "old" run and a "new" run, GenomicSEM may behave
# differently or error out.

suppressPackageStartupMessages({
  library(optparse)
})

opt_list <- list(
  make_option(c("--old_dir"), type = "character",
              help = "Directory containing the reference/old GenomicSEM outputs", metavar = "DIR"),
  make_option(c("--new_dir"), type = "character",
              help = "Directory containing the new GenomicSEM outputs to compare", metavar = "DIR"),
  make_option(c("--files"), type = "character",
              default = "MPH_genomicSEM.RData,S_matrix_MPH.RData,V_matrix_MPH.RData,I_matrix_MPH.RData,N_matrix_MPH.RData",
              help = "Comma-separated list of .RData files to compare [default: %default]",
              metavar = "STR"),
  make_option(c("--out_tsv"), type = "character", default = "",
              help = "Optional: write a TSV summary to this path", metavar = "FILE"),
  make_option(c("--tolerance"), type = "double", default = 1e-10,
              help = "Numeric tolerance for max absolute difference [default: %default]",
              metavar = "FLOAT"),
  make_option(c("--quiet"), action = "store_true", default = FALSE,
              help = "Reduce console output")
)

opt <- parse_args(OptionParser(option_list = opt_list))

say <- function(...) {
  if (!isTRUE(opt$quiet)) cat(..., "\n")
}

die <- function(...) {
  stop(paste0(...), call. = FALSE)
}

if (is.null(opt$old_dir) || opt$old_dir == "") die("Missing --old_dir")
if (is.null(opt$new_dir) || opt$new_dir == "") die("Missing --new_dir")
if (!dir.exists(opt$old_dir)) die("old_dir does not exist: ", opt$old_dir)
if (!dir.exists(opt$new_dir)) die("new_dir does not exist: ", opt$new_dir)

files <- trimws(strsplit(opt$files, ",", fixed = TRUE)[[1]])
files <- files[nzchar(files)]
if (length(files) == 0) die("No files specified via --files")

load_rdata <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  as.list(e)
}

is_matrix_like <- function(x) {
  is.matrix(x) || (!is.null(dim(x)) && length(dim(x)) == 2)
}

safe_dim <- function(x) {
  d <- dim(x)
  if (is.null(d)) return(NA_character_)
  paste0(d, collapse = "x")
}

safe_dn_state <- function(x) {
  dn <- dimnames(x)
  if (is.null(dn)) return("NULL")
  rn <- dn[[1]]
  cn <- dn[[2]]
  rn_state <- if (is.null(rn)) "row=NULL" else paste0("row=", length(rn))
  cn_state <- if (is.null(cn)) "col=NULL" else paste0("col=", length(cn))
  paste(rn_state, cn_state, sep = ",")
}

max_abs_diff <- function(a, b) {
  # assumes dims already checked
  if (!is.numeric(a) || !is.numeric(b)) return(NA_real_)
  # handle Matrix class or other matrix-likes
  am <- tryCatch(as.matrix(a), error = function(e) NULL)
  bm <- tryCatch(as.matrix(b), error = function(e) NULL)
  if (is.null(am) || is.null(bm)) return(NA_real_)
  suppressWarnings(max(abs(am - bm), na.rm = TRUE))
}

results <- list()

add_row <- function(file, object, status, details,
                    class_old = NA_character_, class_new = NA_character_,
                    dim_old = NA_character_, dim_new = NA_character_,
                    dimnames_old = NA_character_, dimnames_new = NA_character_,
                    max_abs = NA_real_) {
  results[[length(results) + 1]] <<- data.frame(
    file = file,
    object = object,
    status = status,
    details = details,
    class_old = class_old,
    class_new = class_new,
    dim_old = dim_old,
    dim_new = dim_new,
    dimnames_old = dimnames_old,
    dimnames_new = dimnames_new,
    max_abs_diff = max_abs,
    stringsAsFactors = FALSE
  )
}

compare_any <- function(a, b, file, obj_path) {
  # recurse on lists
  if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b)) {
    na <- names(a)
    nb <- names(b)
    if (!setequal(na, nb)) {
      add_row(file, obj_path, "DIFF",
              paste0("List element names differ. old={", paste(sort(na), collapse = ","), "} new={", paste(sort(nb), collapse = ","), "}"),
              class_old = paste(class(a), collapse = ","),
              class_new = paste(class(b), collapse = ","))
    }
    common <- intersect(na, nb)
    for (nm in common) {
      compare_any(a[[nm]], b[[nm]], file, paste0(obj_path, "$", nm))
    }
    return(invisible(NULL))
  }

  class_old <- paste(class(a), collapse = ",")
  class_new <- paste(class(b), collapse = ",")
  dim_old <- safe_dim(a)
  dim_new <- safe_dim(b)
  dn_old <- safe_dn_state(a)
  dn_new <- safe_dn_state(b)

  # full attribute-aware equality (captures dimnames differences)
  ae <- tryCatch(all.equal(a, b, check.attributes = TRUE), error = function(e) e$message)
  ae_ok <- isTRUE(ae)

  max_abs <- NA_real_
  if (is_matrix_like(a) && is_matrix_like(b) && identical(dim(a), dim(b))) {
    max_abs <- max_abs_diff(a, b)
  }

  # determine status
  status <- if (ae_ok) "OK" else "DIFF"
  details <- if (ae_ok) "" else {
    if (is.character(ae)) paste(ae, collapse = " ; ") else "Not equal"
  }

  # If numeric difference is tiny but attributes differ, make that obvious
  if (!ae_ok && is.finite(max_abs) && max_abs <= opt$tolerance) {
    details <- paste0(details, " | numeric within tol (", opt$tolerance, ")")
  }

  add_row(file, obj_path, status, details,
          class_old = class_old, class_new = class_new,
          dim_old = dim_old, dim_new = dim_new,
          dimnames_old = dn_old, dimnames_new = dn_new,
          max_abs = max_abs)
}

for (f in files) {
  old_path <- file.path(opt$old_dir, f)
  new_path <- file.path(opt$new_dir, f)

  if (!file.exists(old_path)) {
    add_row(f, "<file>", "MISSING", paste0("Missing in old_dir: ", old_path))
    next
  }
  if (!file.exists(new_path)) {
    add_row(f, "<file>", "MISSING", paste0("Missing in new_dir: ", new_path))
    next
  }

  say("[compare] ", f)

  old_objs <- load_rdata(old_path)
  new_objs <- load_rdata(new_path)

  old_names <- names(old_objs)
  new_names <- names(new_objs)

  if (!setequal(old_names, new_names)) {
    add_row(f, "<objects>", "DIFF",
            paste0("Top-level object names differ. old={", paste(sort(old_names), collapse = ","), "} new={", paste(sort(new_names), collapse = ","), "}"))
  }

  common <- intersect(old_names, new_names)
  for (nm in common) {
    compare_any(old_objs[[nm]], new_objs[[nm]], f, nm)
  }
}

res <- do.call(rbind, results)

# Print a compact summary
n_ok <- sum(res$status == "OK")
n_diff <- sum(res$status == "DIFF")
n_missing <- sum(res$status == "MISSING")

say("\n=== SUMMARY ===")
say("OK: ", n_ok, " | DIFF: ", n_diff, " | MISSING: ", n_missing)

if (n_diff > 0 || n_missing > 0) {
  say("\n--- Differences (first 50) ---")
  print(utils::head(res[res$status != "OK", c("file", "object", "status", "details", "dim_old", "dim_new", "dimnames_old", "dimnames_new", "max_abs_diff")], 50), row.names = FALSE)
} else {
  say("No differences detected.")
}

if (!is.null(opt$out_tsv) && opt$out_tsv != "") {
  dir.create(dirname(opt$out_tsv), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(res, file = opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  say("\nWrote TSV summary: ", opt$out_tsv)
}

if (n_diff > 0 || n_missing > 0) {
  quit(status = 1)
}
