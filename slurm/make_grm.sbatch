#!/usr/bin/env bash
#SBATCH --job-name=mph_make_grm
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
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

SUBSET_PREFIX="${MPH_DIR}/genotypes/subset_genotypes"
BIM_FILE="${SUBSET_PREFIX}.bim"
SNP_INFO="${MPH_DIR}/genotypes/snp_info.csv"
OUT_PREFIX="${MPH_DIR}/grm/gg"

if [[ -f "${OUT_PREFIX}.grm.bin" ]]; then
  echo "[make_grm] Found existing GRM: ${OUT_PREFIX}.grm.bin. Skipping."
  exit 0
fi

if [[ ! -f "$BIM_FILE" ]]; then
  echo "ERROR: Missing BIM file: $BIM_FILE (did prep step run?)" >&2
  exit 1
fi

mkdir -p "${MPH_DIR}/grm" "${MPH_DIR}/logs" "${MPH_DIR}/genotypes"

echo "[make_grm] Writing SNP info to: $SNP_INFO"
perl -e 'print "SNP\n"; while(<>){@c=split /\s+/; print "$c[1]\n"; }' < "$BIM_FILE" > "$SNP_INFO"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "[make_grm] Running MPH --make_grm ..."
"${MPH_BIN}" \
  --make_grm \
  --binary_genotype "${SUBSET_PREFIX}" \
  --snp_info "${SNP_INFO}" \
  --num_threads "${THREADS}" \
  --out "${OUT_PREFIX}"

echo "[make_grm] Done."

