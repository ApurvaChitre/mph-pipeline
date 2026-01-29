#!/usr/bin/env bash
#SBATCH --job-name=mph_cohort_grms
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --account=csd795
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --export=ALL

set -eo pipefail

: "${CONFIG_FILE:?CONFIG_FILE not set (pass via sbatch --export=ALL,CONFIG_FILE=/path/to/config.sh)}"
# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${MPH_DIR:?MPH_DIR must be set in config}"
: "${MPH_FUNCTS_R:?MPH_FUNCTS_R must be set in config}"

# Repo directory (see note in slurm/prep_inputs.sbatch)
REPO_DIR="${PIPELINE_REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
if [[ ! -d "${REPO_DIR}" ]]; then
  echo "ERROR: REPO_DIR does not exist: ${REPO_DIR}" >&2
  echo "Hint: submit via bin/submit_mph_pipeline.sh (it exports PIPELINE_REPO_DIR)." >&2
  exit 2
fi
cd "${REPO_DIR}"

# Optional conda env activation (for R packages)
# NOTE: avoid sourcing ~/.bashrc in Slurm jobs (can break in non-interactive shells).
if [[ -n "${CONDA_ENV:-}" ]]; then
  CONDA_SH_PATH="${CONDA_SH:-}"

  if [[ -z "${CONDA_SH_PATH}" ]]; then
    # Common defaults
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
  fi

  if [[ ! -f "${CONDA_SH_PATH}" ]]; then
    echo "ERROR: CONDA_ENV is set (${CONDA_ENV}) but conda.sh could not be found." >&2
    echo "Set CONDA_SH in your config (run 'conda info --base' on the login node)." >&2
    exit 2
  fi

  # shellcheck disable=SC1090
  source "${CONDA_SH_PATH}"
  conda activate "${CONDA_ENV}"

  echo "[env] Activated conda env: ${CONDA_ENV}"
  echo "[env] Rscript: $(command -v Rscript || echo NA)"
fi

command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found. Set CONDA_ENV in config or load an R module." >&2; exit 1; }


GG_PREFIX="${MPH_DIR}/grm/gg"
COV_FILE="${MPH_DIR}/pheno/indicator_covariates.csv"

if [[ ! -f "${GG_PREFIX}.grm.bin" ]]; then
  echo "ERROR: Missing gg GRM: ${GG_PREFIX}.grm.bin (did make_grm step finish?)" >&2
  exit 1
fi
if [[ ! -f "${COV_FILE}" ]]; then
  echo "ERROR: Missing covariate file: ${COV_FILE} (did prep step finish?)" >&2
  exit 1
fi


echo "[cohort_grms] Creating custom GRMs in ${MPH_DIR}/grm"
Rscript "${REPO_DIR}/scripts/make_cohort_grms.R" \
  --workdir="${MPH_DIR}" \
  --mph_functs="${MPH_FUNCTS_R}"

echo "[cohort_grms] Done."

