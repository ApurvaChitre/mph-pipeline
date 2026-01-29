# Scripts

This directory contains the **command-line** scripts used by the Slurm workflow:

- `make_pheno_cohort_project_dict.R`
  - Builds `pheno/pheno_cohort_project_dict.csv` from your `processed_data_ready.csv`
  - Supports `--include_traits` / `--exclude_traits` (comma list or `@file.txt`)

- `stack_pheno_create_cov.R`
  - Stacks phenotypes into a single `combined_phenotype.csv`
  - Creates the cohort indicator matrix `indicator_covariates.csv`
  - Supports **two input modes**:
    - `--processed_file=processed_data_ready.csv` (recommended; no per-trait files needed)
    - `--pheno_dir=/path/to/gcta_pheno_dir` (legacy per-trait GCTA files)

- `make_cohort_grms.R`
  - Masks the full GRM (`grm/gg`) into all required per-cohort and cross-cohort GRMs

- `generate_matrices_I_matrix_no_overlap.R`
  - Converts MPH REML output (`reml/all.mq.vc.csv`) into `S`, `V`, `I`, `N`, `m`
    in a format usable by GenomicSEM.

- `compare_mph_vs_gcta.R`
  - Optional QC helper: compares MPH outputs against GCTA outputs (if you have them)
  - Writes tables + simple base-R plots

- `compare_genomicsem_outputs.R`
  - Regression test helper: compares two directories of GenomicSEM-ready `.RData`
    outputs to confirm object structure + **dimnames** match exactly.

These scripts are called automatically by `bin/submit_mph_pipeline.sh`.
