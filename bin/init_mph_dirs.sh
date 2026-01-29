#!/usr/bin/env bash
set -eo pipefail

# Create the standard MPH run directory structure.
#
# Usage:
#   bash init_mph_dirs.sh /path/to/mph_run
#
# If no argument is provided, creates ./mph

MPH_DIR="${1:-mph}"

subdirs=(grm genotypes reml pheno gsem logs tmp)

mkdir -p "$MPH_DIR"
for sd in "${subdirs[@]}"; do
  mkdir -p "$MPH_DIR/$sd"
done

echo "Created MPH directory structure under: $MPH_DIR"
printf "  - %s\n" "${subdirs[@]}"
