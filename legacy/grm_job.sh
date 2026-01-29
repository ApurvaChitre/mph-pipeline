#!/bin/bash

#SBATCH --job-name=grm_creation  # Job name
#SBATCH --time=3:00:00                # Time limit hrs:min:sec (increased walltime)
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust as necessary)
#SBATCH --mem-per-cpu=8G               # Memory per CPU core
#SBATCH --account=csd795               # Account name
#SBATCH --partition=condo              # Partition name
#SBATCH --qos=condo                    # Quality of Service
#SBATCH --export=ALL                   # Export all environment variables
# Change to the directory where the job was submitted
cd $SLURM_SUBMIT_DIR

/tscc/projects/ps-palmer/software/local/src/mph/mph --make_grm --binary_genotype /tscc/projects/ps-palmer/apurva/locomotor/mph/genotypes/subset_genotypes --snp_info snp_info.csv --num_threads 14 --out gg