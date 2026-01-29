#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(stringr)
})

opt_list <- list(
  make_option(c("--mph_vc_file"), type = "character",
              help = "MPH variance-component output (.mq.vc.csv)", metavar = "FILE"),
  make_option(c("--mapping_file"), type = "character",
              help = "Mapping file (pheno_cohort_project_dict.csv)", metavar = "FILE"),
  make_option(c("--gcta_h2_file"), type = "character", default = "",
              help = "GCTA heritability.tsv (optional)", metavar = "FILE"),
  make_option(c("--gcta_rg_file"), type = "character", default = "",
              help = "GCTA genetic_correlation_melted_table.csv (optional)", metavar = "FILE"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory", metavar = "DIR"),
  make_option(c("--prefix"), type = "character", default = "mph_vs_gcta",
              help = "Prefix for output files [default: %default]", metavar = "STR"),
  make_option(c("--mph_metric"), type = "character", default = "var",
              help = "Trait-level MPH metric for plot: 'var' (genetic variance) or 'h2' (heritability via var_g/(var_g+resid)). NOTE: h2 SE not computed [default: %default]",
              metavar = "var|h2"),
  make_option(c("--quiet"), action = "store_true", default = FALSE,
              help = "Reduce console output")
)

opt <- parse_args(OptionParser(option_list = opt_list))

msg <- function(...) {
  if (!isTRUE(opt$quiet)) cat(..., "\n")
}

stop_if_missing <- function(path, label) {
  if (is.null(path) || path == "" || !file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
}

stop_if_missing(opt$mph_vc_file, "--mph_vc_file")
stop_if_missing(opt$mapping_file, "--mapping_file")
if (is.null(opt$outdir) || opt$outdir == "") stop("Missing --outdir")

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

msg("[compare] Reading mapping file: ", opt$mapping_file)
map <- fread(opt$mapping_file)
if (!"pheno_file" %in% names(map)) {
  stop("mapping_file must have a column named 'pheno_file' (got: ", paste(names(map), collapse = ", "), ")")
}

map[, trait_full := gsub("\\.txt$", "", basename(pheno_file))]
map[, trait_short := gsub("^regressedlr_", "", trait_full)]
N <- nrow(map)
if (N < 2) stop("Need at least 2 traits; found N=", N)

trait_full <- map$trait_full
trait_short <- map$trait_short

# Preserve mapping order explicitly (merges may reorder rows)
trait_full_order <- trait_full

msg("[compare] N traits inferred from mapping: ", N)


msg("[compare] Reading MPH VC file: ", opt$mph_vc_file)
# check.names=TRUE makes duplicate column names unique (e.g., appending .1)
vc <- fread(opt$mph_vc_file, check.names = TRUE)
if (!"vc_name" %in% names(vc)) stop("mph_vc_file must contain a 'vc_name' column")
if (!"var" %in% names(vc)) stop("mph_vc_file must contain a 'var' column")
if (!"seV" %in% names(vc)) stop("mph_vc_file must contain a 'seV' column")

vc[, vc_base := basename(vc_name)]

num_genetic <- N * (N + 1) / 2
if (nrow(vc) < num_genetic) {
  stop("mph_vc_file has ", nrow(vc), " rows but expected at least ", num_genetic,
       " (N*(N+1)/2) for N=", N, ".")
}

g <- vc[1:num_genetic]


# ------------------------------------------------------------
# Helper: locate the covariance-matrix column for a component
# In MPH output, covariance matrix columns often appear twice;
# the second copy becomes unique with a .1 suffix.
# ------------------------------------------------------------
find_cov_col <- function(vcb) {
  # Prefer the second copy (.1) when present
  pat1 <- paste0("(^|\\.)", vcb, "\\.1$")
  cand <- grep(pat1, names(vc), value = TRUE)
  if (length(cand) == 1) return(cand)

  # Fall back: exact end-match without .1
  pat0 <- paste0("(^|\\.)", vcb, "$")
  cand <- grep(pat0, names(vc), value = TRUE)
  if (length(cand) == 1) return(cand)

  stop("Could not find a covariance column for component '", vcb, "'.")
}

cov_cols <- vapply(g$vc_base, find_cov_col, character(1))

# Covariance matrix among genetic components (rows/cols in same order as g$vc_base)
V <- as.matrix(g[, ..cov_cols])
rownames(V) <- g$vc_base
colnames(V) <- g$vc_base


# ------------------------------------------------------------
# Build MPH trait-level table
# ------------------------------------------------------------
mph_trait <- data.table(
  trait_full = trait_full,
  trait_short = trait_short,
  cohort = paste0("cohort", seq_len(N))
)

for (i in seq_len(N)) {
  comp <- paste0("cohort", i)
  row <- g[vc_base == comp]
  if (nrow(row) != 1) stop("Could not uniquely find MPH component: ", comp)
  mph_trait[i, mph_var := row$var]
  mph_trait[i, mph_var_se := row$seV]
}

# Optional MPH heritability (assumes residuals parameterized as err + cohortX_diagonal)
if (opt$mph_metric == "h2") {
  # need err and diagonal components
  err_row <- vc[vc_base == "err"]
  if (nrow(err_row) != 1) stop("Could not uniquely find residual component 'err' in mph_vc_file")
  err_var <- err_row$var

  mph_trait[, resid_var := err_var]
  for (i in seq_len(N - 1)) {
    dname <- paste0("cohort", i, "_diagonal")
    drow <- vc[vc_base == dname]
    if (nrow(drow) != 1) stop("Could not uniquely find residual diagonal component: ", dname)
    mph_trait[i, resid_var := err_var + drow$var]
  }

  mph_trait[, mph_h2 := mph_var / (mph_var + resid_var)]
  # NOTE: we do not compute SE for mph_h2 (needs full covariance for (var, diag, err))
}


# ------------------------------------------------------------
# Optional: read GCTA heritability file
# ------------------------------------------------------------
gcta_h2 <- NULL
if (!is.null(opt$gcta_h2_file) && opt$gcta_h2_file != "") {
  stop_if_missing(opt$gcta_h2_file, "--gcta_h2_file")
  msg("[compare] Reading GCTA heritability: ", opt$gcta_h2_file)

  h2raw <- fread(opt$gcta_h2_file, sep = "\t", header = TRUE, check.names = FALSE)
  # first column is trait name (often blank header)
  setnames(h2raw, 1, "trait")

  h2col <- grep("V\\(G\\)/Vp", names(h2raw), value = TRUE)
  if (length(h2col) == 0) {
    # sometimes check.names may slightly alter; try broader
    h2col <- grep("Vp", names(h2raw), value = TRUE)
  }
  secol <- grep("heritability_SE", names(h2raw), value = TRUE)

  if (length(h2col) == 0) stop("Could not find column 'V(G)/Vp' in gcta_h2_file")
  if (length(secol) == 0) stop("Could not find column 'heritability_SE' in gcta_h2_file")

  gcta_h2 <- h2raw[, .(
    trait = trait,
    gcta_h2 = as.numeric(get(h2col[1])),
    gcta_h2_se = as.numeric(get(secol[1]))
  )]

  mph_trait <- merge(mph_trait, gcta_h2, by.x = "trait_full", by.y = "trait", all.x = TRUE)
  # Re-order to match the original mapping order
  mph_trait <- mph_trait[match(trait_full_order, mph_trait$trait_full)]
}


# ------------------------------------------------------------
# Compute MPH genetic correlations + SE
# ------------------------------------------------------------
correlation_and_se <- function(var1, var2, cov12, cov3x3) {
  if (is.na(var1) || is.na(var2) || is.na(cov12)) return(list(r = NA_real_, se = NA_real_))
  if (var1 <= 0 || var2 <= 0) return(list(r = NA_real_, se = NA_real_))

  r <- cov12 / sqrt(var1 * var2)

  # Delta method using 3x3 sampling covariance of (var1, var2, cov12)
  V11 <- cov3x3[1, 1]
  V22 <- cov3x3[2, 2]
  V33 <- cov3x3[3, 3]
  V12 <- cov3x3[1, 2]
  V13 <- cov3x3[1, 3]
  V23 <- cov3x3[2, 3]

  if (is.na(V11) || is.na(V22) || is.na(V33)) return(list(r = r, se = NA_real_))
  if (cov12 == 0) return(list(r = r, se = NA_real_))

  var_r <- r^2 * (
    (V11 / (4 * var1^2)) +
      (V22 / (4 * var2^2)) +
      (V33 / (cov12^2)) +
      (V12 / (2 * var1 * var2)) -
      (V13 / (var1 * cov12)) -
      (V23 / (var2 * cov12))
  )

  if (!is.finite(var_r) || var_r < 0) {
    return(list(r = r, se = NA_real_))
  }

  list(r = r, se = sqrt(var_r))
}


pairs <- t(combn(seq_len(N), 2))
mph_rg <- data.table(
  i = pairs[, 1],
  j = pairs[, 2]
)
mph_rg[, trait1 := trait_full[i]]
mph_rg[, trait2 := trait_full[j]]
mph_rg[, trait1_short := trait_short[i]]
mph_rg[, trait2_short := trait_short[j]]

mph_rg[, comp1 := paste0("cohort", i)]
mph_rg[, comp2 := paste0("cohort", j)]
mph_rg[, comp12 := paste0("cohort", i, "_cohort", j)]

mph_rg[, var1 := g[match(comp1, vc_base), var]]
mph_rg[, var2 := g[match(comp2, vc_base), var]]
mph_rg[, cov12 := g[match(comp12, vc_base), var]]

mph_rg[, `:=`(rg_mph = NA_real_, rg_mph_se = NA_real_)]
for (k in seq_len(nrow(mph_rg))) {
  c1 <- mph_rg$comp1[k]
  c2 <- mph_rg$comp2[k]
  c12 <- mph_rg$comp12[k]
  cov3 <- V[c(c1, c2, c12), c(c1, c2, c12), drop = FALSE]
  res <- correlation_and_se(mph_rg$var1[k], mph_rg$var2[k], mph_rg$cov12[k], cov3)
  mph_rg$rg_mph[k] <- res$r
  mph_rg$rg_mph_se[k] <- res$se
}

# key for joining
mph_rg[, pair_key := paste(pmin(trait1, trait2), pmax(trait1, trait2), sep = "__")]


# ------------------------------------------------------------
# Optional: read GCTA genetic-correlation file
# ------------------------------------------------------------
gcta_rg <- NULL
if (!is.null(opt$gcta_rg_file) && opt$gcta_rg_file != "") {
  stop_if_missing(opt$gcta_rg_file, "--gcta_rg_file")
  msg("[compare] Reading GCTA genetic correlations: ", opt$gcta_rg_file)

  rgraw <- fread(opt$gcta_rg_file, check.names = FALSE)

  # Sometimes there is a leading unnamed index column
  if (!all(c("trait1", "trait2") %in% names(rgraw))) {
    # try dropping first column
    rgraw <- rgraw[, -1, with = FALSE]
  }
  if (!all(c("trait1", "trait2") %in% names(rgraw))) {
    stop("gcta_rg_file must contain columns trait1 and trait2")
  }

  if (!"genetic_correlation" %in% names(rgraw)) stop("gcta_rg_file missing 'genetic_correlation' column")
  if (!"rG_SE" %in% names(rgraw)) stop("gcta_rg_file missing 'rG_SE' column")

  gcta_rg <- rgraw[, .(
    trait1 = trait1,
    trait2 = trait2,
    rg_gcta = as.numeric(genetic_correlation),
    rg_gcta_se = as.numeric(rG_SE)
  )]
  gcta_rg[, pair_key := paste(pmin(trait1, trait2), pmax(trait1, trait2), sep = "__")]

  mph_rg <- merge(mph_rg, gcta_rg, by = "pair_key", all.x = TRUE)
}


# ------------------------------------------------------------
# Write tables
# ------------------------------------------------------------
trait_out <- file.path(opt$outdir, paste0(opt$prefix, "_traits.tsv"))
rg_out <- file.path(opt$outdir, paste0(opt$prefix, "_rg_pairs.tsv"))

fwrite(mph_trait, trait_out, sep = "\t", na = "NA")
fwrite(mph_rg, rg_out, sep = "\t", na = "NA")

msg("[compare] Wrote: ", trait_out)
msg("[compare] Wrote: ", rg_out)


# ------------------------------------------------------------
# Plot helper (base R) for two-series errorbar plot
# ------------------------------------------------------------
plot_two_series <- function(xlabels, y1, se1, y2, se2,
                            legend1 = "GCTA", legend2 = "MPH",
                            main = "", ylab = "", outfile = "") {
  if (outfile == "") stop("outfile required")

  n <- length(xlabels)
  x <- seq_len(n)
  offset <- 0.15

  # -----------------------------
  # Plot sizing: auto-scale output size so long/vertical x-axis labels
  # don't crush the plotting region.
  #
  # Defaults were tuned for TSCC-style long trait names.
  # You can always re-run with a different prefix/outdir if you want
  # to keep multiple versions.
  # -----------------------------
  max_len <- max(nchar(xlabels), na.rm = TRUE)

  # Width: scale with number of labels (pair plots get wide fast)
  width_px <- max(2400, 220 * n)
  width_px <- min(width_px, 12000)  # cap to avoid absurdly huge files

  # Height: scale with label length (rotated labels need bottom space)
  height_px <- max(900, 600 + 18 * max_len)
  height_px <- min(height_px, 3000)

  # Bottom margin in "lines"; scale a bit with label length.
  bottom_mar <- min(30, max(10, 8 + ceiling(max_len / 3)))

  # Axis text size: shrink slightly when there are lots of labels
  cex_ax <- if (n > 35) 0.55 else if (n > 25) 0.60 else if (n > 15) 0.70 else 0.80

  ymax <- max(c(y1 + se1, y2 + se2), na.rm = TRUE)
  ymin <- min(c(y1 - se1, y2 - se2), na.rm = TRUE)

  png(outfile, width = width_px, height = height_px, res = 200)
  par(mar = c(bottom_mar, 6, 4, 2) + 0.1)

  plot(x, y1, type = "n",
       ylim = c(ymin, ymax),
       xaxt = "n", xlab = "", ylab = ylab,
       main = main)

  axis(1, at = x, labels = xlabels, las = 2, cex.axis = cex_ax)
  grid(nx = NA, ny = NULL)

  # Series 1
  points(x - offset, y1, pch = 16)
  arrows(x - offset, y1 - se1, x - offset, y1 + se1,
         angle = 90, code = 3, length = 0.03)

  # Series 2
  points(x + offset, y2, pch = 15)
  arrows(x + offset, y2 - se2, x + offset, y2 + se2,
         angle = 90, code = 3, length = 0.03)

  legend("topleft", legend = c(legend1, legend2),
         pch = c(16, 15), bty = "n")

  dev.off()
}


# Trait-level plot (requires GCTA h2)
if (!is.null(gcta_h2)) {
  xlabels <- mph_trait$trait_short
  y1 <- mph_trait$gcta_h2
  se1 <- mph_trait$gcta_h2_se

  if (opt$mph_metric == "h2") {
    y2 <- mph_trait$mph_h2
    # no SE for h2 currently
    se2 <- rep(NA_real_, length(y2))
    main <- "Comparison of SNP h\u00b2 (GCTA) and h\u00b2 (MPH) by Trait"
    ylab <- "Heritability (h\u00b2)"
    legend2 <- "MPH h\u00b2"
  } else {
    y2 <- mph_trait$mph_var
    se2 <- mph_trait$mph_var_se
    main <- "Comparison of SNP h\u00b2 (GCTA) and Var (MPH) by Trait"
    ylab <- "Heritability (h\u00b2) and Variance"
    legend2 <- "Var (MPH)"
  }

  out_png <- file.path(opt$outdir, paste0(opt$prefix, "_h2_var_compare.png"))
  plot_two_series(xlabels, y1, se1, y2, se2,
                  legend1 = "SNP h\u00b2 (GCTA)", legend2 = legend2,
                  main = main, ylab = ylab, outfile = out_png)
  msg("[compare] Wrote plot: ", out_png)
}


# rG plot (requires GCTA rg)
if (!is.null(gcta_rg)) {
  # Build x-axis labels
  mph_rg[, pair_label := paste0(trait1_short, "_", trait2_short)]
  xlabels <- mph_rg$pair_label
  y1 <- mph_rg$rg_gcta
  se1 <- mph_rg$rg_gcta_se
  y2 <- mph_rg$rg_mph
  se2 <- mph_rg$rg_mph_se

  out_png <- file.path(opt$outdir, paste0(opt$prefix, "_rg_compare.png"))
  plot_two_series(xlabels, y1, se1, y2, se2,
                  legend1 = "GCTA", legend2 = "MPH",
                  main = "Comparison of genetic correlation estimates: MPH vs GCTA",
                  ylab = "Genetic correlation estimates",
                  outfile = out_png)
  msg("[compare] Wrote plot: ", out_png)
}


msg("[compare] Done.")
