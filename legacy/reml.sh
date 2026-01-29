#!/bin/bash

#SBATCH --job-name=reml               # Job name
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Number of tasks (processes)
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=32G                    # Total memory for the job
#SBATCH --account=csd795             # Account name
#SBATCH --partition=condo            # Partition name
#SBATCH --qos=condo                  # Quality of Service
#SBATCH --export=ALL                 # Export all environment variables



# Define directories and file paths
GRM_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/grm"
PHENO_FILE="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/pheno/combined_phenotype.csv"
COVARIATE_FILE="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/pheno/indicator_covariates.csv"
OUTPUT_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/reml"
GRM_LIST="${GRM_DIR}/chr.grms.txt"
OUTPUT_FILE="${OUTPUT_DIR}/all"


# Generate the GRM list
ls -tr ${GRM_DIR}/*.grm.bin | grep -v "gg.grm.bin" | sed 's/\.grm\.bin//' > ${GRM_LIST}

# Load necessary modules 
conda activate genomic_sem

# Run the mph_ccc_kalivas command
/tscc/projects/ps-palmer/software/local/src/mph/mph --reml \
  --grm_list "${GRM_LIST}" \
  --phenotype "${PHENO_FILE}" \
  --trait trait_value \
  --covariate_file "${COVARIATE_FILE}" \
  --covariate_names all \
  --output "${OUTPUT_FILE}"
