#!/usr/bin/env bash
#SBATCH --job-name=mph_reml
#SBATCH --time=3:00:00
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
: "${MPH_BIN:?MPH_BIN must be set in config}"

GRM_DIR="${MPH_DIR}/grm"
PHENO_FILE="${MPH_DIR}/pheno/combined_phenotype.csv"
COV_FILE="${MPH_DIR}/pheno/indicator_covariates.csv"
MAPPING_FILE="${MPH_DIR}/pheno/pheno_cohort_project_dict.csv"

OUT_DIR="${MPH_DIR}/reml"
OUT_PREFIX="${OUT_DIR}/all"
GRM_LIST="${GRM_DIR}/chr.grms.txt"

if [[ -f "${OUT_PREFIX}.mq.vc.csv" ]]; then
  echo "[reml] Found existing REML output: ${OUT_PREFIX}.mq.vc.csv. Skipping."
  exit 0
fi

mkdir -p "${OUT_DIR}" "${MPH_DIR}/logs"

if [[ ! -f "${MAPPING_FILE}" ]]; then
  echo "ERROR: Missing mapping file: ${MAPPING_FILE}" >&2
  exit 1
fi

# Infer number of cohorts/traits from mapping file
N=$(($(wc -l < "${MAPPING_FILE}") - 1))
if [[ "${N}" -lt 1 ]]; then
  echo "ERROR: Could not infer number of traits from ${MAPPING_FILE}" >&2
  exit 1
fi

# -------------------------------------------------------------------
# Build a deterministic GRM list (do NOT rely on file timestamps)
# Order is important: genetic components first, then residual diagonals.
#
# Genetic components:
#   - cohort1..cohortN
#   - cohort1_cohort2..cohort(N-1)_cohortN
# Residual components:
#   - cohort1_diagonal..cohort(N-1)_diagonal  (MPH adds final diagonal automatically)
# -------------------------------------------------------------------
{
  # Single-cohort genetic GRMs
  for i in $(seq 1 "${N}"); do
    echo "${GRM_DIR}/cohort${i}"
  done

  # Cross-cohort genetic GRMs
  if [[ "${N}" -ge 2 ]]; then
    for i in $(seq 1 $((N-1))); do
      for j in $(seq $((i+1)) "${N}"); do
        echo "${GRM_DIR}/cohort${i}_cohort${j}"
      done
    done
  fi

  # Residual diagonal GRMs (all but last)
  if [[ "${N}" -ge 2 ]]; then
    for i in $(seq 1 $((N-1))); do
      echo "${GRM_DIR}/cohort${i}_diagonal"
    done
  fi
} > "${GRM_LIST}"

# Verify all referenced GRMs exist
missing=0
while read -r prefix; do
  if [[ ! -f "${prefix}.grm.bin" ]]; then
    echo "ERROR: Missing GRM component: ${prefix}.grm.bin" >&2
    missing=1
  fi
done < "${GRM_LIST}"

if [[ "${missing}" -ne 0 ]]; then
  echo "ERROR: One or more GRMs are missing. Did make_cohort_grms finish successfully?" >&2
  exit 1
fi

echo "[reml] GRM list written to: ${GRM_LIST}"
echo "[reml] Running MPH --reml ..."
"${MPH_BIN}" --reml   --grm_list "${GRM_LIST}"   --phenotype "${PHENO_FILE}"   --trait trait_value   --covariate_file "${COV_FILE}"   --covariate_names all   --output "${OUT_PREFIX}"

echo "[reml] Done."
