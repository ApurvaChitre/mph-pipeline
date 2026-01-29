#!/usr/bin/env bash

#SBATCH --job-name=mph_gcta_compare

# NOTE: Resources are set at submission time by bin/submit_mph_pipeline.sh

set -eo pipefail

if [[ -z "${CONFIG_FILE:-}" ]]; then
  echo "ERROR: CONFIG_FILE not set (pass via sbatch --export=ALL,CONFIG_FILE=/path/to/config.sh)" >&2
  exit 2
fi

echo "[mph] gcta_compare starting"
echo "HOST=$(hostname)"
echo "CONFIG_FILE=${CONFIG_FILE}"

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

# Repo directory (see note in slurm/prep_inputs.sbatch)
REPO_DIR="${PIPELINE_REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
if [[ ! -d "${REPO_DIR}" ]]; then
  echo "ERROR: REPO_DIR does not exist: ${REPO_DIR}" >&2
  echo "Hint: submit via bin/submit_mph_pipeline.sh (it exports PIPELINE_REPO_DIR)." >&2
  exit 2
fi
cd "${REPO_DIR}"

# Activate conda env if requested
if [[ -n "${CONDA_ENV:-}" ]]; then
  CONDA_SH_CANDIDATES=(
    "${CONDA_SH:-}"
    "$HOME/miniconda3/etc/profile.d/conda.sh"
    "$HOME/anaconda3/etc/profile.d/conda.sh"
    "/opt/conda/etc/profile.d/conda.sh"
  )
  for c in "${CONDA_SH_CANDIDATES[@]}"; do
    if [[ -n "${c}" && -f "${c}" ]]; then
      # shellcheck source=/dev/null
      source "${c}"
      conda activate "${CONDA_ENV}" || {
        echo "ERROR: Failed to activate conda env: ${CONDA_ENV}" >&2
        exit 2
      }
      break
    fi
  done
fi

VC_FILE="${MPH_DIR}/reml/all.mq.vc.csv"
MAP_FILE="${MPH_DIR}/pheno/pheno_cohort_project_dict.csv"
OUTDIR="${MPH_DIR}/qc_gcta"

mkdir -p "${OUTDIR}"

if [[ -z "${GCTA_H2_FILE:-}" && -z "${GCTA_RG_FILE:-}" ]]; then
  echo "[mph] RUN_GCTA_COMPARE is enabled but GCTA_H2_FILE/GCTA_RG_FILE are empty. Nothing to do." >&2
  exit 2
fi

echo "[mph] MPH vc file: ${VC_FILE}"
echo "[mph] Mapping file: ${MAP_FILE}"
echo "[mph] GCTA h2 file: ${GCTA_H2_FILE:-<unset>}"
echo "[mph] GCTA rg file: ${GCTA_RG_FILE:-<unset>}"
echo "[mph] OUTDIR: ${OUTDIR}"

Rscript "${REPO_DIR}/scripts/compare_mph_vs_gcta.R" \
  --mph_vc_file "${VC_FILE}" \
  --mapping_file "${MAP_FILE}" \
  --gcta_h2_file "${GCTA_H2_FILE:-}" \
  --gcta_rg_file "${GCTA_RG_FILE:-}" \
  --outdir "${OUTDIR}" \
  --prefix "mph_vs_gcta" \
  --mph_metric "${MPH_GCTA_METRIC:-var}"

echo "[mph] gcta_compare done"
