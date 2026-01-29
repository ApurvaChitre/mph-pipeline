#!/bin/bash

#SBATCH --job-name=grm_creation         # Job name
#SBATCH --time=3:00:00                  # Time limit hrs:min:sec
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem-per-cpu=8G                # Memory per CPU core
#SBATCH --account=csd795                # Account name
#SBATCH --partition=condo               # Partition name
#SBATCH --qos=condo                     # Quality of Service
#SBATCH --export=ALL                    # Export all environment variables

# Define directories
BASE_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/"
GENOTYPE_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/genotypes"
GRM_DIR="$BASE_DIR/grm"
SNP_INFO_FILE="$BASE_DIR/genotypes/snp_info.csv"
GENOTYPE_FILE="$GENOTYPE_DIR/subset_genotypes.bim"
GENOTYPE_BASENAME="$(basename "$GENOTYPE_FILE" .bim)"
OUTPUT_FILE="$GRM_DIR/gg"

# Change to the directory where the job was submitted
cd $SLURM_SUBMIT_DIR

# Ensure directories exist
mkdir -p "$GRM_DIR"

# Step 1: Generate SNP info
perl -e '
    print "SNP\n"; 
    while(<>) {
        @c = split /\s+/; 
        print "$c[1]\n";
    }
' < "$GENOTYPE_FILE" > "$SNP_INFO_FILE"

# Step 2: Create GRM
/tscc/projects/ps-palmer/software/local/src/mph/mph \
    --make_grm \
    --binary_genotype "$GENOTYPE_DIR/$GENOTYPE_BASENAME" \
    --snp_info "$SNP_INFO_FILE" \
    --num_threads 14 \
    --out "$OUTPUT_FILE"
