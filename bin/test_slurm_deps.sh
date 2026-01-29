#!/usr/bin/env bash
set -eo pipefail

usage() {
  cat <<'EOF'
Quick Slurm dependency sanity test.

Usage:
  bash bin/test_slurm_deps.sh --config /path/to/config.sh

Submits two tiny jobs:
  job1: prints hostname/date, sleeps 5s
  job2: runs after job1 succeeds (afterok)

Writes logs to: $MPH_DIR/logs/
EOF
}

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-c)
      CONFIG_FILE="$2"; shift 2;;
    --help|-h)
      usage; exit 0;;
    *)
      echo "Unknown argument: $1" >&2
      usage; exit 1;;
  esac
done

if [[ -z "${CONFIG_FILE}" ]]; then
  echo "ERROR: --config is required" >&2
  exit 1
fi
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: config file not found: ${CONFIG_FILE}" >&2
  exit 1
fi

CONFIG_FILE="$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${CONFIG_FILE}")"
# shellcheck disable=SC1090
source "${CONFIG_FILE}"
: "${MPH_DIR:?MPH_DIR must be set in config}"

LOG_DIR="${MPH_DIR}/logs"
mkdir -p "${LOG_DIR}"

SBATCH_COMMON=()
if [[ -n "${SLURM_ACCOUNT:-}" ]]; then SBATCH_COMMON+=(--account="${SLURM_ACCOUNT}"); fi
if [[ -n "${SLURM_PARTITION:-}" ]]; then SBATCH_COMMON+=(--partition="${SLURM_PARTITION}"); fi
if [[ -n "${SLURM_QOS:-}" ]]; then SBATCH_COMMON+=(--qos="${SLURM_QOS}"); fi

submit_wrap() {
  local name="$1"; shift
  local dep="${1:-}"; shift || true
  local cmd="$1"; shift || true

  local args=(sbatch --parsable --job-name="${name}" --output="${LOG_DIR}/${name}_%j.out" --error="${LOG_DIR}/${name}_%j.err")
  args+=("${SBATCH_COMMON[@]}")
  if [[ -n "${dep}" ]]; then
    dep="${dep%%;*}"; dep="${dep%% *}"
    args+=(--dependency="afterok:${dep}")
  fi
  args+=(--wrap "${cmd}")

  echo "+ ${args[*]}" >&2
  local out
  out="$("${args[@]}")"
  out="${out%%;*}"; out="${out%% *}"
  echo "${out}"
}

echo "Submitting dependency sanity test..."
job1=$(submit_wrap "mph_deptest_1" "" "echo JOB1_START; hostname; date; sleep 5; echo JOB1_DONE")
job2=$(submit_wrap "mph_deptest_2" "${job1}" "echo JOB2_START; hostname; date; echo JOB2_DONE")

echo
echo "Submitted:"
echo "  job1: ${job1}"
echo "  job2: ${job2} (afterok:${job1})"
echo
echo "Monitor:"
echo "  squeue -u ${USER}"
echo "  sacct -j ${job1},${job2} --format=JobID,JobName%20,State,ExitCode,Elapsed"
echo
echo "Logs:"
echo "  ${LOG_DIR}/mph_deptest_1_<jobid>.out"
echo "  ${LOG_DIR}/mph_deptest_2_<jobid>.out"
