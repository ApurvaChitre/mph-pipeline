# Running the MPH pipeline (TSCC / Slurm)

This repo contains both:
- **Background/explanation** (main `README.md`)
- A **coworker-friendly, automated pipeline** that runs on **TSCC (Slurm)** via `sbatch` job dependencies.

If you are **not** on TSCC / Slurm, you can still run the workflow by executing the scripts in order (see “Manual run (no Slurm)” below).

---

## Environment / installation

See `docs/INSTALL.md` for options (recommended: conda env). A quick sanity check:

```bash
bash bin/check_deps.sh config/my_run.sh
```

### Optional but recommended: run the Slurm debug job first

This is a quick preflight that runs on a compute node and checks:
- your config can be sourced
- output directories can be created
- conda env activation works in non-interactive Slurm jobs (if `CONDA_ENV` is set)
- required commands and R packages are available

```bash
bash bin/submit_debug_env.sh --config config/my_run.sh
```

Then open the log in `${MPH_DIR}/logs/` (it will say PASS/FAIL at the bottom).

### Conda note (important on TSCC)

Slurm jobs run in **non-interactive shells**, so `conda activate` often requires sourcing `conda.sh`.

If your jobs can’t activate the env, set this in your config:

```bash
# on the login node:
conda info --base
# then:
CONDA_SH="/path/from/that/command/etc/profile.d/conda.sh"
```


## Starting point (what you need)

Minimum inputs:

1) **Phenotypes (wide CSV)**
- `processed_data_ready.csv` containing:
  - `rfid` column (sample IDs)
  - one column per trait starting with `regressedlr_...`

2) **Genotypes (PLINK)**
- A PLINK prefix: `GENO_BFILE=/path/to/genotypes` (expects `.bed/.bim/.fam`)
- Autosomes only (default in config is `CHR_RANGE="1-20"`)

Outputs you care about:

- `reml/all.mq.vc.csv` (MPH variance components + SEs + covariances of estimates)
- `gsem/MPH_genomicSEM.RData` (list with `S`, `V`, `I`, `N`, `m` for GenomicSEM)

---

## Quickstart (one command on TSCC)

1) Clone the repo:
```bash
git clone <your_repo_url>
cd <repo>
```

2) Copy/edit config:
```bash
cp config/config.example.sh config/my_run.sh
nano config/my_run.sh
```

3) Submit the full pipeline (Slurm job chain):
```bash
bash bin/submit_mph_pipeline.sh --config config/my_run.sh
```

Monitor:
```bash
squeue -u $USER
```

Logs:
- `${MPH_DIR}/logs/<jobname>_<jobid>.out`
- `${MPH_DIR}/logs/<jobname>_<jobid>.err`

---

## Optional QC: compare MPH vs GCTA (h² and rG)

If you have **precomputed GCTA results** for the same trait set, you can have the pipeline run a small QC step at the end to compare:

- Trait-level estimates (GCTA h² vs MPH genetic variance or MPH h²)
- Pairwise genetic correlations (rG)

Enable it in your config:

```bash
RUN_GCTA_COMPARE="1"
GCTA_H2_FILE="/path/to/heritability.tsv"
GCTA_RG_FILE="/path/to/genetic_correlation_melted_table.csv"

# If phenotypes are standardized (variance ~1), MPH genetic variance ≈ h².
# Otherwise you can try MPH_GCTA_METRIC="h2".
MPH_GCTA_METRIC="var"
```

Outputs:

- Tables: `${MPH_DIR}/qc_gcta/mph_vs_gcta_trait_table.tsv`, `${MPH_DIR}/qc_gcta/mph_vs_gcta_rg_table.tsv`
- Plots: `${MPH_DIR}/qc_gcta/mph_vs_gcta_h2_or_var.png`, `${MPH_DIR}/qc_gcta/mph_vs_gcta_rg.png`

This QC stage is *optional* so coworkers can still run MPH → GenomicSEM inputs even if they don’t have GCTA outputs.

## 7-trait example (recommended)

This repo’s examples are set up for the **7-trait run** by excluding 2 traits:

- `regressedlr_meinhardt_open_field_test_track_length_total_cm`
- `regressedlr_r01_giordano_open_field_15_min_distance_m`

You can do this either in the config directly:
```bash
EXCLUDE_TRAITS="regressedlr_meinhardt_open_field_test_track_length_total_cm,regressedlr_r01_giordano_open_field_15_min_distance_m"
```

…or by using the example file:
```bash
EXCLUDE_TRAITS_FILE="config/exclude_traits.example.txt"
```

If you find it easier to explicitly list the traits you want to run (instead of excluding a few), you can use **include mode**. When you use an include list, the order you provide defines cohort numbering (cohort1, cohort2, …):

```bash
# comma-separated
INCLUDE_TRAITS="u01_tom_jhou_session1_locomotor1,p50_hao_chen_open_field_first_15"

# or file mode (one trait per line)
INCLUDE_TRAITS_FILE="/path/to/traits_keep.txt"
```

The pipeline automatically infers **N traits** from the resulting mapping file and generates the GenomicSEM matrices accordingly.

---

## Important note about sample overlap (MPH requirement)

MPH requires that the GRM sample set and phenotype/covariate sample set **match**.

- If the genotype set contains **extra samples** (not in phenotype IDs), MPH will error.
- If some phenotype IDs are **missing in genotypes**, MPH will also error.

**This pipeline enforces matching** by creating a keep-list from the phenotype IDs and subsetting genotypes to exactly those IDs (unless you explicitly provide an already-subset genotype prefix).

---

## Phenotype input options

### Option A (default): use `processed_data_ready.csv` only
This is the lowest-dependency approach: you do **not** need per-trait GCTA phenotype files.

Config:
- `PHENO_SOURCE="processed"`
- `PROCESSED_FILE=/path/to/processed_data_ready.csv`

The pipeline will:
1) create `pheno/pheno_cohort_project_dict.csv`
2) stack the `regressedlr_*` columns into `pheno/combined_phenotype.csv`
3) create `pheno/indicator_covariates.csv`

### Option B (legacy): per-trait GCTA phenotype files
If you already have per-trait GCTA phenotype files (3 columns: `ID FID trait_value`) you can use them.

Config:
- `PHENO_SOURCE="gcta"`
- `PHENO_DIR=/path/to/pheno_dir` containing files named like `regressedlr_<trait>.txt`

---

## Trait selection for the phenotype-prep step (mapping file)

The mapping file `pheno/pheno_cohort_project_dict.csv` is created by `scripts/make_pheno_cohort_project_dict.R`.
It scans `processed_data_ready.csv` for columns that start with `regressedlr_` and assigns sequential cohort names (`cohort1`, `cohort2`, …) to the selected traits.

You can control trait selection in two ways:

- **Include list**: keep only the listed traits (**and the order you provide becomes the cohort order**).
- **Exclude list**: start from all `regressedlr_*` traits and drop a subset.

Both `--include_traits` and `--exclude_traits` accept:
- comma-separated traits, **or**
- “file mode” with a leading `@` (one trait per line).

Trait names can be provided **with or without** the `regressedlr_` prefix.

```bash
# Examples:
# Rscript scripts/make_pheno_cohort_project_dict.R --processed=... --outdir=... \
#   --include_traits=u01_tom_jhou_session1_locomotor1,p50_hao_chen_open_field_first_15
#
# Rscript scripts/make_pheno_cohort_project_dict.R --processed=... --outdir=... \
#   --exclude_traits=p50_shelly_flagel_2014_hab1_total_distance,u01_peter_kalivas_oft_distance_1_first_15
#
# File mode:
#   --include_traits=@/path/to/traits_keep.txt
#   --exclude_traits=@/path/to/traits_skip.txt
```

In the **automated pipeline**, you set these via config variables:
- `INCLUDE_TRAITS` / `INCLUDE_TRAITS_FILE`
- `EXCLUDE_TRAITS` / `EXCLUDE_TRAITS_FILE`

If you provide both include and exclude lists, the pipeline applies the include list first, then removes any excluded traits.

---

## Genotype subsetting options

### Default behavior (recommended)
The pipeline will create:
- `${MPH_DIR}/genotypes/subset_ids.txt` (from the phenotype IDs)
- `${MPH_DIR}/genotypes/subset_genotypes.*` (PLINK subset; autosomes only)

### If you *already* have a subset genotype set
If you already created the PLINK subset elsewhere, set:
```bash
SUBSET_GENO_BFILE="/path/to/your_subset_prefix"
```

The pipeline will symlink/copy it into:
- `${MPH_DIR}/genotypes/subset_genotypes.*`

So downstream steps still work unchanged.

---

## Stage-by-stage run (Slurm or manual)
You can run each stage one-by-one.

**On TSCC/Slurm**, always pass your config file:

```bash
CONFIG_FILE="$(realpath config/my_run.sh)"
sbatch --export=ALL,CONFIG_FILE="$CONFIG_FILE" slurm/prep_inputs.sbatch
```

If you are **not** on Slurm, run the commands shown below directly (and ignore the `sbatch ...` lines).


If you’re not on TSCC/Slurm, you can run the stages manually in order.

(Assuming you are in the repo root, and you’ve set the same variables as the config file.)

1) Mapping:
```bash
Rscript scripts/make_pheno_cohort_project_dict.R \
  --processed "$PROCESSED_FILE" \
  --outdir "$MPH_DIR/pheno" \
  --outname pheno_cohort_project_dict.csv \
  --exclude_traits="$EXCLUDE_TRAITS"   # or --include_traits=... (see above)
```

2) Stack phenotypes + covariates:
```bash
Rscript scripts/stack_pheno_create_cov.R --processed_file "$PROCESSED_FILE" --pheno_cohort_file "$MPH_DIR/pheno/pheno_cohort_project_dict.csv" --output_pheno_file "$MPH_DIR/pheno/combined_phenotype.csv" --output_cov_file "$MPH_DIR/pheno/indicator_covariates.csv"
```

3) Subset genotypes:
```bash
cut -d',' -f1 "$MPH_DIR/pheno/indicator_covariates.csv" | tail -n +2 > "$MPH_DIR/genotypes/subset_ids.txt"

plink --bfile "$GENO_BFILE" --keep-fam "$MPH_DIR/genotypes/subset_ids.txt" --chr 1-20 --make-bed --out "$MPH_DIR/genotypes/subset_genotypes"
```

4) Make GRM:
```bash
sbatch --export=ALL,CONFIG_FILE="$CONFIG_FILE" slurm/make_grm.sbatch   # or run the command inside it
```

5) Create custom GRMs:
```bash
Rscript scripts/make_cohort_grms.R --workdir "$MPH_DIR" --mph_functs "$MPH_FUNCTS_R"
```

6) MPH REML:
```bash
sbatch --export=ALL,CONFIG_FILE="$CONFIG_FILE" slurm/reml.sbatch       # or run the command inside it
```

7) GenomicSEM matrices:
```bash
sbatch --export=ALL,CONFIG_FILE="$CONFIG_FILE" slurm/genomicsem_inputs.sbatch
```

---

## Troubleshooting

- If the prep stage fails with a “sample mismatch” error, it means your genotype set does not contain exactly the same sample IDs as the phenotype set. Fix the mismatch (or regenerate the subset) and re-run.
- If you want to re-run cleanly, delete the output directory:
```bash
rm -rf "$MPH_DIR"
```
and submit again.

