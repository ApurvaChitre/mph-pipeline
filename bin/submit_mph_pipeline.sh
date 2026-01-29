#!/usr/bin/env bash
set -eo pipefail

usage() {
  cat <<'EOF'
Submit the full MPH workflow as a Slurm job chain.

Usage:
  bash bin/submit_mph_pipeline.sh --config /path/to/config.sh

This will submit (with dependencies):
  1) prep inputs
  2) make GRM (gg)
  3) make cohort GRMs
  4) MPH REML
  5) generate GenomicSEM inputs

EOF
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-c)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIG_FILE" ]]; then
  echo "ERROR: --config is required" >&2
  exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: config file not found: $CONFIG_FILE" >&2
  exit 1
fi

# Make config path absolute for compute nodes
CONFIG_FILE="$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$CONFIG_FILE")"

# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${MPH_DIR:?MPH_DIR must be set in config}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${MPH_DIR}/logs"
mkdir -p "$LOG_DIR"

# Common sbatch options from config (optional)
SBATCH_COMMON=()
if [[ -n "${SLURM_ACCOUNT:-}" ]]; then SBATCH_COMMON+=(--account="$SLURM_ACCOUNT"); fi
if [[ -n "${SLURM_PARTITION:-}" ]]; then SBATCH_COMMON+=(--partition="$SLURM_PARTITION"); fi
if [[ -n "${SLURM_QOS:-}" ]]; then SBATCH_COMMON+=(--qos="$SLURM_QOS"); fi

# Always write logs to MPH_DIR/logs
SBATCH_COMMON+=(--output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err")

submit() {
  local script="$1"
  local dep="${2:-}"
  shift 2 || true

  # IMPORTANT:
  # Slurm executes a *spooled copy* of each .sbatch script (e.g. under /cm/local/apps/slurm/var/spool/...)
  # so those scripts cannot reliably infer the repo root from "$0" / "${BASH_SOURCE[0]}".
  # We therefore:
  #   1) force the working directory to the repo root via --chdir
  #   2) pass PIPELINE_REPO_DIR as an exported variable
  local export_vars="ALL,CONFIG_FILE=${CONFIG_FILE},PIPELINE_REPO_DIR=${REPO_DIR}"
  local args=(sbatch --parsable --chdir="$REPO_DIR" --export="$export_vars")
  args+=("${SBATCH_COMMON[@]}")

  if [[ -n "$dep" ]]; then
    # dep may contain a federation suffix (e.g. 12345;cluster) or whitespace; sanitize
    dep="${dep%%;*}"
    dep="${dep%% *}"
    args+=(--dependency="afterok:${dep}")
  fi

  # Any extra sbatch args
  if [[ $# -gt 0 ]]; then
    args+=("$@")
  fi

  args+=("$script")

  # IMPORTANT: print debug to STDERR so command substitution captures ONLY the jobid
  echo "+ ${args[*]}" >&2

  local out
  out="$("${args[@]}")"
  # Slurm may return "jobid;cluster"; keep only the numeric jobid
  out="${out%%;*}"
  out="${out%% *}"
  echo "$out"
}

# Stage-specific resource overrides (optional)
PREP_ARGS=()
if [[ -n "${PREP_TIME:-}" ]]; then PREP_ARGS+=(--time="$PREP_TIME"); fi
if [[ -n "${PREP_CPUS:-}" ]]; then PREP_ARGS+=(--cpus-per-task="$PREP_CPUS"); fi
if [[ -n "${PREP_MEM:-}" ]]; then PREP_ARGS+=(--mem="$PREP_MEM"); fi

GRM_ARGS=()
if [[ -n "${GRM_TIME:-}" ]]; then GRM_ARGS+=(--time="$GRM_TIME"); fi
if [[ -n "${GRM_CPUS:-}" ]]; then GRM_ARGS+=(--cpus-per-task="$GRM_CPUS"); fi
if [[ -n "${GRM_MEM:-}" ]]; then GRM_ARGS+=(--mem="$GRM_MEM"); fi

COHORT_ARGS=()
if [[ -n "${COHORT_GRMS_TIME:-}" ]]; then COHORT_ARGS+=(--time="$COHORT_GRMS_TIME"); fi
if [[ -n "${COHORT_GRMS_CPUS:-}" ]]; then COHORT_ARGS+=(--cpus-per-task="$COHORT_GRMS_CPUS"); fi
if [[ -n "${COHORT_GRMS_MEM:-}" ]]; then COHORT_ARGS+=(--mem="$COHORT_GRMS_MEM"); fi

REML_ARGS=()
if [[ -n "${REML_TIME:-}" ]]; then REML_ARGS+=(--time="$REML_TIME"); fi
if [[ -n "${REML_CPUS:-}" ]]; then REML_ARGS+=(--cpus-per-task="$REML_CPUS"); fi
if [[ -n "${REML_MEM:-}" ]]; then REML_ARGS+=(--mem="$REML_MEM"); fi

GSEM_ARGS=()
if [[ -n "${GSEM_TIME:-}" ]]; then GSEM_ARGS+=(--time="$GSEM_TIME"); fi
if [[ -n "${GSEM_CPUS:-}" ]]; then GSEM_ARGS+=(--cpus-per-task="$GSEM_CPUS"); fi
if [[ -n "${GSEM_MEM:-}" ]]; then GSEM_ARGS+=(--mem="$GSEM_MEM"); fi

# Optional: GCTA comparison stage (QC)
GCTA_ARGS=()
if [[ -n "${GCTA_COMPARE_TIME:-}" ]]; then GCTA_ARGS+=(--time="$GCTA_COMPARE_TIME"); fi
if [[ -n "${GCTA_COMPARE_CPUS:-}" ]]; then GCTA_ARGS+=(--cpus-per-task="$GCTA_COMPARE_CPUS"); fi
if [[ -n "${GCTA_COMPARE_MEM:-}" ]]; then GCTA_ARGS+=(--mem="$GCTA_COMPARE_MEM"); fi

echo "Submitting MPH pipeline with config: $CONFIG_FILE"
echo "Repo: $REPO_DIR"
echo "Outputs: $MPH_DIR"
echo

jid_prep=$(submit "${REPO_DIR}/slurm/prep_inputs.sbatch" "" "${PREP_ARGS[@]}")
jid_grm=$(submit "${REPO_DIR}/slurm/make_grm.sbatch" "${jid_prep}" "${GRM_ARGS[@]}")
jid_cohort=$(submit "${REPO_DIR}/slurm/make_cohort_grms.sbatch" "${jid_grm}" "${COHORT_ARGS[@]}")
jid_reml=$(submit "${REPO_DIR}/slurm/reml.sbatch" "${jid_cohort}" "${REML_ARGS[@]}")
jid_gsem=$(submit "${REPO_DIR}/slurm/genomicsem_inputs.sbatch" "${jid_reml}" "${GSEM_ARGS[@]}")

jid_gcta_compare=""
if [[ "${RUN_GCTA_COMPARE:-0}" == "1" ]]; then
  # Only submit if at least one of the GCTA inputs is provided.
  if [[ -n "${GCTA_H2_FILE:-}" || -n "${GCTA_RG_FILE:-}" ]]; then
    jid_gcta_compare=$(submit "${REPO_DIR}/slurm/gcta_compare.sbatch" "${jid_gsem}" "${GCTA_ARGS[@]}")
  else
    echo "NOTE: RUN_GCTA_COMPARE=1 but GCTA_H2_FILE and GCTA_RG_FILE are empty; skipping QC stage." >&2
  fi
fi

cat <<EOF

Submitted job chain:
  prep_inputs:          ${jid_prep}
  make_grm:             ${jid_grm}
  make_cohort_grms:     ${jid_cohort}
  reml:                 ${jid_reml}
  genomicsem_inputs:    ${jid_gsem}
  gcta_compare (QC):     ${jid_gcta_compare:-<not submitted>}

Logs:
  ${LOG_DIR}/<jobname>_<jobid>.out
  ${LOG_DIR}/<jobname>_<jobid>.err

Monitor:
  squeue -u ${USER}

Key outputs (after completion):
  ${MPH_DIR}/reml/all.mq.vc.csv
  ${MPH_DIR}/gsem/MPH_genomicSEM.RData
  ${MPH_DIR}/qc_gcta/ (if gcta_compare is enabled)

EOF

