#!/bin/bash

# Main directory name
main_dir="mph"

# Subdirectories
subdirs=("grm" "genotypes" "reml" "pheno" "scripts")

# Create the main directory
mkdir -p "$main_dir"

# Create subdirectories within the main directory
for subdir in "${subdirs[@]}"; do
  mkdir -p "$main_dir/$subdir"
done

echo "Directory structure created under $main_dir:"
