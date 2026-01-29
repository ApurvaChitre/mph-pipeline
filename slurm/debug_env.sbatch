#!/usr/bin/env bash
#SBATCH --job-name=mph_debug_env
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=csd795
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --export=ALL

# "Soft strict" mode: fail on errors, but don't die on unset variables (conda init can be noisy).
set -eo pipefail

echo "=== MPH DEBUG ENV ==="
echo "host=$(hostname)"
echo "date=$(date)"
echo "pwd=$(pwd)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "CONFIG_FILE=${CONFIG_FILE:-<unset>}"

if [[ -z "${CONFIG_FILE:-}" ]]; then
  echo "ERROR: CONFIG_FILE not set. Submit with: sbatch --export=ALL,CONFIG_FILE=/path/to/config.sh ..." >&2
  exit 2
fi
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: CONFIG_FILE not found: $CONFIG_FILE" >&2
  exit 2
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"

echo ""
echo "== Config summary =="
echo "MPH_DIR=${MPH_DIR:-<unset>}"
echo "PHENO_SOURCE=${PHENO_SOURCE:-<unset>}"
echo "PROCESSED_FILE=${PROCESSED_FILE:-<unset>}"
echo "PHENO_DIR=${PHENO_DIR:-<unset>}"
echo "GENO_BFILE=${GENO_BFILE:-<unset>}"
echo "SUBSET_GENO_BFILE=${SUBSET_GENO_BFILE:-<unset>}"
echo "PLINK_BIN=${PLINK_BIN:-plink}"
echo "MPH_BIN=${MPH_BIN:-<unset>}"
echo "MPH_FUNCTS_R=${MPH_FUNCTS_R:-<unset>}"
echo "CONDA_ENV=${CONDA_ENV:-<none>}"
echo "CONDA_SH=${CONDA_SH:-<auto>}"
echo ""

errs=0
check_file () {
  local f="$1"
  local label="$2"
  if [[ -z "$f" ]]; then
    echo "MISSING: $label (empty)"
    errs=1
    return 0
  fi
  if [[ ! -e "$f" ]]; then
    echo "MISSING: $label -> $f"
    errs=1
    return 0
  fi
  echo "OK: $label -> $f"
}

check_cmd () {
  local c="$1"
  if ! command -v "$c" >/dev/null 2>&1; then
    echo "MISSING: command '$c' on PATH"
    errs=1
    return 0
  fi
  echo "OK: $c -> $(command -v "$c")"
}

# MPH_DIR + subdirs
if [[ -z "${MPH_DIR:-}" ]]; then
  echo "ERROR: MPH_DIR is not set in config." >&2
  exit 2
fi
mkdir -p "$MPH_DIR"/{logs,pheno,genotypes,tmp} || { echo "ERROR: cannot create subdirs in MPH_DIR=$MPH_DIR" >&2; exit 3; }
echo "OK: created/verified MPH_DIR subdirs:"
ls -ld "$MPH_DIR" "$MPH_DIR/logs" "$MPH_DIR/pheno" "$MPH_DIR/genotypes" "$MPH_DIR/tmp"

echo ""
echo "== Input file checks =="

if [[ "${PHENO_SOURCE:-processed}" == "processed" ]]; then
  check_file "${PROCESSED_FILE:-}" "PROCESSED_FILE"
elif [[ "${PHENO_SOURCE}" == "gcta" ]]; then
  if [[ -d "${PHENO_DIR:-}" ]]; then
    echo "OK: PHENO_DIR -> ${PHENO_DIR}"
  else
    echo "MISSING: PHENO_DIR -> ${PHENO_DIR:-<unset>}"
    errs=1
  fi
fi

# genotype prefix checks
if [[ -n "${SUBSET_GENO_BFILE:-}" ]]; then
  check_file "${SUBSET_GENO_BFILE}.bed" "SUBSET_GENO_BFILE.bed"
  check_file "${SUBSET_GENO_BFILE}.bim" "SUBSET_GENO_BFILE.bim"
  check_file "${SUBSET_GENO_BFILE}.fam" "SUBSET_GENO_BFILE.fam"
else
  check_file "${GENO_BFILE:-}.bed" "GENO_BFILE.bed"
  check_file "${GENO_BFILE:-}.bim" "GENO_BFILE.bim"
  check_file "${GENO_BFILE:-}.fam" "GENO_BFILE.fam"
fi

echo ""
echo "== Tooling checks (before conda) =="
check_cmd awk
check_cmd cut
check_cmd wc
check_cmd perl
check_cmd "${PLINK_BIN:-plink}" || true
check_cmd Rscript || true
check_cmd python || true

echo ""
echo "== Conda activation (if requested) =="

if [[ -n "${CONDA_ENV:-}" ]]; then
  CONDA_SH_PATH="${CONDA_SH:-}"
  if [[ -z "${CONDA_SH_PATH}" ]]; then
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
  fi

  if [[ ! -f "${CONDA_SH_PATH}" ]]; then
    echo "ERROR: CONDA_ENV is set (${CONDA_ENV}) but conda.sh could not be found." >&2
    echo "Set CONDA_SH in config (run 'conda info --base' on login node)." >&2
    errs=1
  else
    # shellcheck disable=SC1090
    source "${CONDA_SH_PATH}"
    if conda activate "${CONDA_ENV}"; then
      echo "OK: activated conda env: ${CONDA_ENV}"
      echo "PATH=$PATH"
      check_cmd Rscript
      check_cmd python || true
      check_cmd "${PLINK_BIN:-plink}" || true

      echo ""
      echo "== R package sanity check =="
      if command -v Rscript >/dev/null 2>&1; then
        Rscript - <<'RS'
pkgs <- c("data.table","optparse","dplyr","readr","stringr","Matrix")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) {
  cat("Missing R packages:\n")
  cat(paste0("  - ", missing), sep="\n")
  quit(status=1)
} else {
  cat("All required R packages are installed.\n")
}
RS
        if [[ $? -ne 0 ]]; then
          errs=1
        fi
      fi
    else
      echo "ERROR: conda activate failed for env: ${CONDA_ENV}" >&2
      errs=1
    fi
  fi
fi

echo ""
echo "== MPH binary checks =="
if [[ -n "${MPH_BIN:-}" ]]; then
  if [[ -x "${MPH_BIN}" ]]; then
    echo "OK: MPH_BIN executable -> ${MPH_BIN}"
  else
    echo "MISSING/NOT EXECUTABLE: MPH_BIN -> ${MPH_BIN}"
    errs=1
  fi
else
  echo "MISSING: MPH_BIN is not set in config"
  errs=1
fi

if [[ -n "${MPH_FUNCTS_R:-}" ]]; then
  if [[ -f "${MPH_FUNCTS_R}" ]]; then
    echo "OK: MPH_FUNCTS_R -> ${MPH_FUNCTS_R}"
  else
    echo "MISSING: MPH_FUNCTS_R -> ${MPH_FUNCTS_R}"
    errs=1
  fi
fi

echo ""
if [[ $errs -ne 0 ]]; then
  echo "DEBUG RESULT: FAIL (see messages above)" >&2
  exit 1
fi

echo "DEBUG RESULT: PASS"
echo "=== DONE ==="
