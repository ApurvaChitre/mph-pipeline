#!/usr/bin/env bash
#SBATCH --job-name=mph_genomicsem_inputs
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --account=csd795
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --export=ALL

set -eo pipefail

: "${CONFIG_FILE:?CONFIG_FILE not set (pass via sbatch --export=ALL,CONFIG_FILE=/path/to/config.sh)}"
# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${MPH_DIR:?MPH_DIR must be set in config}"

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


VC_FILE="${MPH_DIR}/reml/all.mq.vc.csv"
MAPPING_FILE="${MPH_DIR}/pheno/pheno_cohort_project_dict.csv"
BIM_FILE="${MPH_DIR}/genotypes/subset_genotypes.bim"
OUT_DIR="${MPH_DIR}/gsem"

if [[ -f "${OUT_DIR}/MPH_genomicSEM.RData" ]]; then
  echo "[gsem] Found existing output: ${OUT_DIR}/MPH_genomicSEM.RData. Skipping."
  exit 0
fi

if [[ ! -f "${VC_FILE}" ]]; then
  echo "ERROR: Missing vc file: ${VC_FILE} (did REML finish?)" >&2
  exit 1
fi
if [[ ! -f "${MAPPING_FILE}" ]]; then
  echo "ERROR: Missing mapping file: ${MAPPING_FILE}" >&2
  exit 1
fi
if [[ ! -f "${BIM_FILE}" ]]; then
  echo "ERROR: Missing BIM file: ${BIM_FILE}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}" "${MPH_DIR}/logs"

# infer N traits from mapping file (header + N rows)
N=$(($(wc -l < "${MAPPING_FILE}") - 1))
if [[ "${N}" -lt 1 ]]; then
  echo "ERROR: Could not infer number of traits from ${MAPPING_FILE}" >&2
  exit 1
fi

# build "cohort1,cohort2,..."
traits=""
for i in $(seq 1 "${N}"); do
  if [[ -z "$traits" ]]; then
    traits="cohort${i}"
  else
    traits="${traits},cohort${i}"
  fi
done


echo "[gsem] Generating GenomicSEM matrices (N=${N}) ..."
Rscript "${REPO_DIR}/scripts/generate_matrices_I_matrix_no_overlap.R" \
  --num_traits "${N}" \
  --traits "${traits}" \
  --vc_file "${VC_FILE}" \
  --mapping_file "${MAPPING_FILE}" \
  --bim_file "${BIM_FILE}" \
  --output_dir "${OUT_DIR}"

echo "[gsem] Done."

