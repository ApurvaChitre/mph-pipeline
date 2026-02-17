# Genetic Covariance and Variance-Covariance Matrix Estimation for Heterogeneous Stock Rats

This repository provides resources and scripts to estimate genetic covariance and variance-covariance matrices required for the **Genomic SEM R package**. The estimation is performed using **[MPH (MINQUE for Partitioning Heritability)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11093526/)**, customized specifically for research involving **Heterogeneous Stock (HS) rats**.

**Pipeline automation (TSCC/Slurm):** For an end-to-end runner (with `sbatch` dependencies and a single config file), see `docs/RUNNING_PIPELINE.md`.

> **Note:** The automated runner in this repo is designed for the **TSCC HPC cluster (Slurm)**. If you are running on a different system (no Slurm), you can still use the underlying scripts and run the stages manually (see `docs/RUNNING_PIPELINE.md`).

### Scope of the Code
This repository is designed for cases where **traits are measured in different cohorts**, rather than in a single cohort with multiple traits. In this context, cohort does not refer to a phenotyping batch. Specifically:
- **Cohorts represent distinct Heterogeneous Stock rats (HS) projects arising from collaborations with the Palmer lab, with each project focusing on a unique set of traits.**
- For example, **Traits A, B, and C are measured in cohorts 1, 2, and 3, respectively.**

The goal is to estimate the genetic covariance matrix between the three traits, accounting for the cohort structure and ensuring accurate variance-covariance modeling.

### Goal of the Analysis

The objective is to compute two key matrices for use in the **Genomic SEM R package**:

1. **Matrix S (Genetic Covariance Matrix):**
   - A symmetric $n \times n$ matrix for $n$ traits.
   - **Diagonal**: Genetic variance for each trait.
   - **Off-diagonal**: Genetic covariances between traits.

2. **Matrix V (Variance-Covariance Matrix of Parameter Estimates):**
   - A symmetric $\frac{n(n+1)}{2} \times \frac{n(n+1)}{2}$ matrix.
   - **Diagonal**: Squared standard errors of genetic variances.
   - **Off-diagonal**: Covariances of standard errors for genetic variances and covariances.

This method overcomes limitations of **LDSC (Linkage Disequilibrium Score Regression)**, which provides out-of-bound estimates for HS rats.

### Correspondence Highlights with the author of MPH, Jicai Jiang

The function for the problem has not been hard-coded in MPH, but MPH does have a flexible (though requiring more data manipulation) approach to doing so. You can stack the phenotypes across cohorts as a single trait. By customizing genomic relationship matrices, the estimation problem for multiple cohorts is reduced to a single-variate multi-component model.

Suppose **G** is a genomic relationship matrix across three cohorts (1-3):

$$
\mathbf{G} =
\begin{pmatrix}
G_{11} & G_{12} & G_{13} \\
G_{21} & G_{22} & G_{23} \\
G_{31} & G_{32} & G_{33}
\end{pmatrix}
$$

Suppose $\mathbf{y}_1, \mathbf{y}_2, \mathbf{y}_3$ are phenotypes for cohorts 1-3. 

$$
\begin{pmatrix} 
\mathbf{y}_1 \\ 
\mathbf{y}_2 \\ 
\mathbf{y}_3 
\end{pmatrix} 
$$

The **variance-covariance matrix** for phenotypes measured across three cohorts is as follows:

![Variance-Covariance Matrix](formula.jpg)

- **Matrix on the Left**: Represents genetic variance ($\sigma^2_g$) and genetic covariance between traits across cohorts.
  - **Diagonal elements**: Genetic variance for each cohort.
  - **Off-diagonal elements**: Genetic covariance between cohorts.

- **Matrix on the Right**: Represents residual variance ($\sigma^2_e$).
  - **Diagonal elements**: Residual variance for each cohort.
  - **Off-diagonal elements**: Zero, indicating no residual covariance between cohorts.

This can be expanded as follows:
![Expanded Variance-Covariance Matrix](expanded_formula.jpg)

### Overview of steps for the analysis (3 Cohorts)

The whole analysis (for 3 cohorts) includes a few major steps:

1. **Compute G Across Cohorts**  
   Compute the genomic relationship matrix (G) across all cohorts.

2. **Create Custom Matrices**  
   - Generate the six genomic matrices and the 3 residual matrices (shown above).  
   - Only 3 residual matrices are required because MPH automatically includes a diagonal matrix as the sixth one.  
   - Use the provided R scripts for manipulating GRMs [GRM Scripts](https://jiang18.github.io/mph/scripts/#grm-inputoutput)
   - This step should be quick if the total number of rats across cohorts is around 10k.

3. **List Custom Relationship Matrices**  
   - Include all 9 custom relationship matrices in the GRM list file for MPH.

4. **Run MPH**  
   - Execute MPH, which will output the matrices **S** and **V** in the `.mq.vc.csv` file.

## NEW INSTRUCTIONS TO RUN THE PIPELINE

**📝 Note:** Ensure that the **RFIDs/Samples** overlap between the **genotype** and **phenotype** data. MPH will error out even if the genotype file contains additional samples not present in the phenotype data, even if there is overlap. This approach is intended for autosomes only.

**Quickstart (TSCC/Slurm):**
```bash
cp config/config.example.sh config/my_run.sh
nano config/my_run.sh
bash bin/submit_mph_pipeline.sh --config config/my_run.sh
```

For full details (inputs, options like using an already-subset genotype set, and manual non-Slurm execution), see `docs/RUNNING_PIPELINE.md`.

### File Formats

For details about the required file formats, please refer to:

- The **[MPH Manual](https://jiang18.github.io/mph/options/)**, which provides comprehensive documentation on all input and output file structures.
- The Palmer Lab TSCC directory at `/tscc/projects/ps-palmer/s3/data/genomic_sem`, which contains example files and datasets for genomic SEM analyses.

These resources will provide all necessary guidance on preparing and using the required files for analysis.

### Acknowledgments

Special thanks to **[Jicai Jiang](https://cals.ncsu.edu/animal-science/people/jicai-jiang/)** for his invaluable guidance on using the **MPH framework** and providing the code constructing custom matrices.

For more details, visit the **[MPH manual](https://jiang18.github.io/mph/)**.
