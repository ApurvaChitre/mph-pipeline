#!/usr/bin/env bash
# -------------------------------------------------------------------
# MPH pipeline configuration (TSCC / Slurm)
#
# Copy this file, edit paths, then run:
#   bash bin/submit_mph_pipeline.sh --config config/my_run.sh
# -------------------------------------------------------------------

# ----- Required: where to write outputs (will create subdirs inside) -----
# Example:
#   MPH_DIR="/tscc/projects/ps-palmer/$USER/paper_2026/mph"
MPH_DIR="/tscc/projects/ps-palmer/YOUR_USER/paper_2026/mph"

# -------------------------------------------------------------------
# Phenotype inputs
# -------------------------------------------------------------------
# Choose ONE phenotype source:
#   - processed : stack regressedlr_* columns directly from PROCESSED_FILE (recommended)
#   - gcta      : use pre-formatted per-trait GCTA phenotype files in PHENO_DIR
PHENO_SOURCE="processed"

# Wide processed file used for:
#   1) building the mapping file (pheno_cohort_project_dict.csv)
#   2) stacking phenotypes (if PHENO_SOURCE="processed")
PROCESSED_FILE="/tscc/projects/ps-palmer/YOUR_USER/paper_2026/individual_projects/processed_data_ready.csv"

# Only used if PHENO_SOURCE="gcta"
# Directory containing per-trait phenotype files named like:
#   regressedlr_<trait>.txt
# Each file should have 3 columns: ID, FID, trait_value (headerless is OK).
PHENO_DIR="/tscc/projects/ps-palmer/YOUR_USER/paper_2026/individual_projects/data/pheno"

# ----- Trait selection (optional) -----
# 7-trait example: exclude these two traits
EXCLUDE_TRAITS="regressedlr_meinhardt_open_field_test_track_length_total_cm,regressedlr_r01_giordano_open_field_15_min_distance_m"
EXCLUDE_TRAITS_FILE=""   # alternatively: "config/exclude_traits.example.txt" (one trait per line; with or without regressedlr_ prefix)

# Similarly for include traits (if set, ONLY these are used, in this order)
INCLUDE_TRAITS=""
INCLUDE_TRAITS_FILE=""   # one trait per line; with or without regressedlr_ prefix

# -------------------------------------------------------------------
# Genotype inputs
# -------------------------------------------------------------------
# Option A (default): provide the full genotype set; pipeline will subset it to phenotype IDs
GENO_BFILE="/tscc/projects/ps-palmer/YOUR_USER/paper_2026/individual_projects/genotypes/genotypes"

# Option B: if you ALREADY have a subset genotype prefix (bed/bim/fam),
# set it here and the pipeline will use it (symlink/copy into MPH_DIR/genotypes/subset_genotypes.*)
SUBSET_GENO_BFILE=""

# Autosomes only by default
CHR_RANGE="1-20"

# Use --keep-fam (1-column file of IDs) OR --keep (2-column FID IID file)
PLINK_KEEP_FLAG="--keep-fam"
PLINK_BIN="plink"

# -------------------------------------------------------------------
# Software paths
# -------------------------------------------------------------------
MPH_BIN="/tscc/projects/ps-palmer/software/local/src/mph/mph"
MPH_FUNCTS_R="/tscc/projects/ps-palmer/software/local/src/mph/mph_functs.R"

# Optional: conda env that has R + required packages
# (data.table, readr, dplyr, stringr, optparse, Matrix)
CONDA_ENV="genomic_sem"

# Optional: path to conda.sh for non-interactive shells (recommended on TSCC).
# If conda activate fails inside Slurm, set this explicitly:
#   CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
CONDA_SH=""

# -------------------------------------------------------------------
# Slurm defaults (edit to match your allocation)
# -------------------------------------------------------------------
SLURM_ACCOUNT="csd795"
SLURM_PARTITION="condo"
SLURM_QOS="condo"

# Resources (override if needed)
PREP_TIME="1:00:00"
PREP_CPUS="4"
PREP_MEM="16G"

GRM_TIME="3:00:00"
GRM_CPUS="8"
GRM_MEM="64G"

COHORT_GRMS_TIME="1:00:00"
COHORT_GRMS_CPUS="8"
COHORT_GRMS_MEM="32G"

REML_TIME="3:00:00"
REML_CPUS="8"
REML_MEM="32G"

GSEM_TIME="0:30:00"
GSEM_CPUS="2"
GSEM_MEM="8G"

# -------------------------------------------------------------------
# Optional QC: compare MPH estimates to GCTA results (if you have them)
# -------------------------------------------------------------------
# If enabled, this runs *after* genomicsem_inputs and writes outputs to:
#   ${MPH_DIR}/qc_gcta/
RUN_GCTA_COMPARE="0"   # set to 1 to enable

# Paths to precomputed GCTA outputs (formats shown in docs)
GCTA_H2_FILE=""        # e.g. /path/to/heritability.tsv
GCTA_RG_FILE=""        # e.g. /path/to/genetic_correlation_melted_table.csv

# MPH trait-level metric to compare to GCTA h2:
#   - var : use MPH genetic variance estimates (best when phenotypes are standardized)
#   - h2  : compute var_g / (var_g + var_e) using MPH diagonal components
MPH_GCTA_METRIC="var"

# QC resources
GCTA_COMPARE_TIME="0:30:00"
GCTA_COMPARE_CPUS="1"
GCTA_COMPARE_MEM="4G"
