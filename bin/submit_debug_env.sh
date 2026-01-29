#!/usr/bin/env bash
set -eo pipefail

usage() {
  echo "Usage: bash bin/submit_debug_env.sh --config /path/to/config.sh" >&2
}

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIG_FILE" ]]; then
  echo "ERROR: --config is required" >&2
  usage
  exit 1
fi

CONFIG_FILE="$(realpath "$CONFIG_FILE")"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: config not found: $CONFIG_FILE" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${MPH_DIR:?MPH_DIR must be set in config}"

mkdir -p "$MPH_DIR/logs"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="${REPO_DIR}/slurm/debug_env.sbatch"

# Slurm defaults (can still be overridden by user via sbatch args if needed)
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-}"
QOS="${SLURM_QOS:-}"

# Ensure the debug job runs from the repo root (and knows where the repo is)
export_vars="ALL,CONFIG_FILE=${CONFIG_FILE},PIPELINE_REPO_DIR=${REPO_DIR}"
args=(sbatch --parsable --chdir="$REPO_DIR" --export="$export_vars")
[[ -n "$ACCOUNT"   ]] && args+=(--account="$ACCOUNT")
[[ -n "$PARTITION" ]] && args+=(--partition="$PARTITION")
[[ -n "$QOS"       ]] && args+=(--qos="$QOS")

args+=(--output="$MPH_DIR/logs/%x_%j.out" --error="$MPH_DIR/logs/%x_%j.err")

echo "+ ${args[*]} $SCRIPT" >&2
jid="$("${args[@]}" "$SCRIPT")"
jid="${jid%%;*}"
jid="${jid%% *}"
echo "Submitted debug job: $jid"
echo "Logs: $MPH_DIR/logs/mph_debug_env_${jid}.out / .err"
