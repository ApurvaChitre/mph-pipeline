#!/bin/bash

#SBATCH --job-name=cohort_grms       # Job name
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Number of tasks (processes)
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=32G                    # Total memory for the job
#SBATCH --account=csd795             # Account name
#SBATCH --partition=condo            # Partition name
#SBATCH --qos=condo                  # Quality of Service
#SBATCH --export=ALL                 # Export all environment variables


conda activate genomic_sem

# Run the R script
Rscript /tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph/scripts/make_cohort_grms.R --n_cohorts=8
