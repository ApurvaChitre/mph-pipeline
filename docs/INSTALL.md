# Installing dependencies

This pipeline is designed to run on **TSCC (Slurm)**, where MPH is already installed at:

- `/tscc/projects/ps-palmer/software/local/src/mph/mph`
- `/tscc/projects/ps-palmer/software/local/src/mph/mph_functs.R`

You **do not** need to install MPH yourself on TSCC (just set `MPH_BIN` / `MPH_FUNCTS_R` in your config).

## Minimum requirements

You need these available in your batch job environment:

- `Rscript`
- R packages: `data.table`, `optparse`, `dplyr`, `readr`, `stringr`, `Matrix`
- `plink` (or set `PLINK_BIN=/path/to/plink`)
- standard unix tools: `awk`, `cut`, `wc`, `perl`

## Recommended: use a conda/mamba env (recommended)

If you want a consistent environment across coworkers, use the included example:

```bash
# If you have mamba, use it (usually much faster):
mamba env create -f env/environment.example.yml -n genomic_sem

# Otherwise conda works too:
# conda env create -f env/environment.example.yml -n genomic_sem

conda activate genomic_sem
```

Then in your run config set:

```bash
CONDA_ENV="genomic_sem"
```

### Conda activation inside Slurm jobs

TSCC Slurm jobs run in **non-interactive** shells. To make `conda activate` work reliably, the pipeline sources `conda.sh`.

- You should **not** execute `conda.sh` directly (it may not be executable).
- You should **source** it, e.g.:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
```

If your conda install is not in `$HOME/miniconda3`, set this in your run config:

```bash
# run on login node:
conda info --base
# then set in config:
CONDA_SH="/path/from/that/command/etc/profile.d/conda.sh"
```

### Exporting YOUR working environment to YAML

If you already have a working environment (e.g. `genomic_sem`) and want to share it:

```bash
conda activate genomic_sem

# Option A: smaller, reproducible (recommended)
conda env export --from-history | sed '/^prefix:/d' > env/environment.yml

# Option B: fully pinned (bigger)
conda env export --no-builds | sed '/^prefix:/d' > env/environment.full.yml
```

Commit `env/environment.yml` (or `environment.full.yml`) if you want coworkers to recreate it exactly.

### Troubleshooting env creation on TSCC

If you see errors like:

- `SafetyError: ... appears to be corrupted`
- `ClobberError: ... incompatible packages due to a shared path`

Common fixes:

1) **Use strict channel priority** (helps avoid clobbering between channels):

```bash
conda config --set channel_priority strict
```

2) **Avoid mixing `defaults` with `conda-forge`** in the YAML (this repo’s example uses only `conda-forge` + `bioconda`).

3) **Clear the package cache** if a package download is corrupted:

```bash
conda clean --all
```

Then retry the create command.

## If you prefer NOT to use conda

That's fine. Just ensure `Rscript` and the required packages are available in batch jobs (e.g. via TSCC modules),
and set in your run config:

```bash
CONDA_ENV=""
```

## Sanity-check your environment

Run:

```bash
bash bin/check_deps.sh config/my_run.sh
```

It will check for required commands and required R packages.
