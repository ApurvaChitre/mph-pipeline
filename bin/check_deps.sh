#!/usr/bin/env bash
set -eo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash bin/check_deps.sh /path/to/config.sh" >&2
  exit 1
fi

CONFIG_FILE="$1"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: config file not found: $CONFIG_FILE" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"

# If CONDA_ENV is set, try to activate it for consistent R/plink availability.
if [[ -n "${CONDA_ENV:-}" ]]; then
  CONDA_SH_PATH="${CONDA_SH:-}"
  if [[ -z "${CONDA_SH_PATH}" ]]; then
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
  fi

  if [[ -f "${CONDA_SH_PATH}" ]]; then
    set +u
    # shellcheck disable=SC1090
    source "${CONDA_SH_PATH}"
    conda activate "${CONDA_ENV}" || true
  else
    echo "WARN: CONDA_ENV is set (${CONDA_ENV}) but CONDA_SH not found. Skipping activation." >&2
  fi
fi

echo "== Checking required paths/vars =="

: "${MPH_BIN:?MPH_BIN must be set in config}"
: "${MPH_FUNCTS_R:?MPH_FUNCTS_R must be set in config}"
: "${PLINK_BIN:=plink}"

if [[ ! -x "$MPH_BIN" ]]; then
  echo "ERROR: MPH_BIN not executable: $MPH_BIN" >&2
  exit 1
fi
if [[ ! -f "$MPH_FUNCTS_R" ]]; then
  echo "ERROR: MPH_FUNCTS_R not found: $MPH_FUNCTS_R" >&2
  exit 1
fi

echo "MPH_BIN:      $MPH_BIN"
echo "MPH_FUNCTS_R: $MPH_FUNCTS_R"
echo "PLINK_BIN:    $PLINK_BIN"
echo "CONDA_ENV:    ${CONDA_ENV:-<none>}"

echo ""
echo "== Checking commands on PATH =="

need_cmd() {
  local c="$1"
  if ! command -v "$c" >/dev/null 2>&1; then
    echo "ERROR: missing command on PATH: $c" >&2
    return 1
  fi
  echo "OK: $c -> $(command -v "$c")"
}

# If PLINK_BIN is an absolute path, check executable; otherwise check PATH
if [[ "$PLINK_BIN" == /* ]]; then
  [[ -x "$PLINK_BIN" ]] || { echo "ERROR: PLINK_BIN not executable: $PLINK_BIN" >&2; exit 1; }
  echo "OK: PLINK_BIN -> $PLINK_BIN"
else
  need_cmd "$PLINK_BIN"
fi

need_cmd Rscript
need_cmd awk
need_cmd cut
need_cmd wc
need_cmd perl

echo ""
echo "== Checking R packages =="

Rscript - <<'RS'
pkgs <- c("data.table","optparse","dplyr","readr","stringr","Matrix")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) {
  cat("Missing R packages:\n")
  cat(paste0("  - ", missing), sep="\n")
  cat("\n\nInstall them in your conda env (recommended) or your R library.\n")
  quit(status=1)
} else {
  cat("All required R packages are installed.\n")
}
RS

echo ""
echo "All checks passed."
