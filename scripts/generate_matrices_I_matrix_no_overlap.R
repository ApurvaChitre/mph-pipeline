#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(optparse)
library(Matrix)


# Set up command-line argument parsing
option_list <- list(
  make_option(c("-n", "--num_traits"), type = "integer", default = 6, 
              help = "Number of traits/cohorts in the analysis [default: %default]"),
  make_option(c("-t", "--traits"), type = "character", default = "cohort1,cohort2,cohort3,cohort4,cohort5,cohort6",
              help = "Comma-separated list of trait names [default: %default]"),
  make_option(c("-v", "--vc_file"), type = "character", default = "all.mq.vc.csv",
              help = "Variance-covariance file path [default: %default]"),
  make_option(c("-m", "--mapping_file"), type = "character", default = "pheno_cohort_project_dict.csv",
              help = "Mapping file path [default: %default]"),
  make_option(c("-b", "--bim_file"), type = "character", default = NULL,
              help = "Path to the .bim file for SNP count [required]"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Output directory for matrices [default: current directory]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))



# Assign arguments to variables
N <- opt$num_traits
traits <- unlist(strsplit(opt$traits, ","))
vc_file <- opt$vc_file
mapping_file <- opt$mapping_file
bim_file <- opt$bim_file
output_dir <- opt$output_dir
# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE)

# Read input files
df_vc <- read.csv(vc_file, header = TRUE, stringsAsFactors = FALSE)
df_mapping <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)



# Calculate m value (number of SNPs) from .bim file
m_value <- as.numeric(system(paste("wc -l", bim_file, "| awk '{print $1}'"), intern = TRUE))
cat("m value (number of SNPs):", m_value, "\n")


# Process the Variance-Covariance (VC) Data
df_vc$vc_name <- basename(df_vc$vc_name)

# Step 1: Subset the dataframe to keep only n(n+1)/2 rows
num_rows_to_keep <- (N * (N + 1)) / 2  # e.g., For N=6, this is 21 rows
df_subset <- df_vc[1:num_rows_to_keep, ]

# Step 2: Extract initial columns
df_initial <- df_subset %>%
  select(vc_name, var, seV)

# Step 3: Extract columns after skipping
columns_to_skip <- 10 + num_rows_to_keep + (N - 1)
df_after_skip <- df_subset[, (columns_to_skip + 1):(columns_to_skip + num_rows_to_keep)]
colnames(df_after_skip) <- df_vc$vc_name[1:num_rows_to_keep]

# Combine initial and additional columns
df_combined <- cbind(df_initial, df_after_skip)

# Generate S_matrix
# 1. Initialize the Genetic Covariance Matrix (Matrix S)
S_matrix <- matrix(0, nrow = N, ncol = N)

# Fill the diagonal with variances (rows 1 to 6)
diag(S_matrix) <- df_combined$var[1:N]

# Fill the off-diagonal with covariances (rows 7 to 21)
index <- N + 1
for (i in 1:(N-1)) {
  for (j in (i+1):N) {
    S_matrix[i, j] <- df_combined$var[index]
    S_matrix[j, i] <- df_combined$var[index]
    index <- index + 1
  }
}

# Save S_matrix
save(S_matrix, file = file.path(output_dir, "S_matrix_MPH.RData"))

# Generate V_matrix
num_traits <- length(traits)
num_elements <- num_traits * (num_traits + 1) / 2  # Number of elements in the matrix

# Initialize the matrix with zeros
vc_matrix <- matrix(0, nrow = num_elements, ncol = num_elements)
# Create labels for the rows and columns based on variances (Vii) and covariances (Cij)
labels <- c()
for (i in 1:num_traits) {
  for (j in i:num_traits) {
    if (i == j) {
      labels <- c(labels, paste0("V", i))  # Diagonal (variance) elements
    } else {
      labels <- c(labels, paste0("C", i, j))  # Off-diagonal (covariance) elements
    }
  }
}

# Assign labels to the rows and columns of the matrix
rownames(vc_matrix) <- labels
colnames(vc_matrix) <- labels

# Fill the V_matrix
# Function to extract cohort names from row/column names
extract_cohort <- function(name) {
  if (grepl("^V", name)) {
    # Variance case (Vii -> cohort[i])
    cohort_idx <- as.numeric(gsub("V", "", name))  # Extract numeric part after "V"
    return(traits[cohort_idx])
  } else if (grepl("^C", name)) {
    # Covariance case (Cij -> cohort[i]_cohort[j])
    cohort_indices <- as.numeric(strsplit(sub("C", "", name), "")[[1]])  # Split after "C"
    return(paste0(traits[cohort_indices[1]], "_", traits[cohort_indices[2]]))
  }
  return(NA)
}

# Loop to print the row and column names and construct cohort names
for (i in 1:num_elements) {
  for (j in i:num_elements) {
    row_name <- rownames(vc_matrix)[i]
    col_name <- colnames(vc_matrix)[j]
    
    # Extract the corresponding cohort names using the extract_cohort function
    row_cohort <- extract_cohort(row_name)
    col_cohort <- extract_cohort(col_name)
    
    # Printing the constructed row and column cohort names
    print(paste("Row:", row_name, "->", row_cohort, "| Column:", col_name, "->", col_cohort))
    
    # Fill diagonal (variance) elements with squared standard errors
    if (i == j) {
      vc_matrix[i, j] <- df_combined$seV[which(df_combined$vc_name == row_cohort)]^2
    }
    # Fill off-diagonal elements
    else {
      vc_matrix[i, j] <- df_combined[which(df_combined$vc_name == row_cohort), col_cohort]
      vc_matrix[j, i] <- df_combined[which(df_combined$vc_name == row_cohort), col_cohort]
    }
  }
}



# ===========================================
# (AFTER your existing S and V code)
# ===========================================


# --- Generate I_matrix (identity: no sample overlap, λGC = 1) -------------
I_matrix           <- diag(1,  nrow = N,  ncol = N)   # 1’s on the diagonal, 0 elsewhere
dimnames(I_matrix) <- list(traits, traits)            # keep trait names for clarity
save(I_matrix, file = file.path(output_dir, "I_matrix_MPH.RData"))
cat("Identity I_matrix generated and saved.\n")
# --------------------------------------------------------------------------




# Read the mapping file
df_mapping <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)

# Extract sample sizes and traits
sample_sizes <- df_mapping$N
traits <- gsub("regressedlr_", "", df_mapping$pheno_file)
traits <- gsub(".txt", "", traits)

# Number of traits
num_traits <- length(traits)

# Initialize the N matrix
N_matrix <- matrix(0, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))

# Fill the diagonal with sample sizes
diag(N_matrix) <- sample_sizes

# Fill the off-diagonal with geometric means of sample sizes
for (i in 1:(num_traits - 1)) {
  for (j in (i + 1):num_traits) {
    geometric_mean <- sqrt(sample_sizes[i] * sample_sizes[j])
    N_matrix[i, j] <- geometric_mean
    N_matrix[j, i] <- geometric_mean
  }
}

# Save the N matrix
save(N_matrix, file = file.path(output_dir, "N_matrix_MPH.RData"))

cat("N_matrix generated and saved.\n")
# Print the N matrix for verification
print(N_matrix)

# Function to check if the matrix has any unfilled (NA) elements
check_empty_elements <- function(matrix) {
  empty_elements <- which(is.na(matrix), arr.ind = TRUE)
  if (length(empty_elements) == 0) {
    print("No empty elements found in the matrix.")
  } else {
    print(paste("Empty elements found at positions:", toString(empty_elements)))
  }
}

# Function to check for symmetry
check_symmetry <- function(matrix) {
  for (i in 1:nrow(matrix)) {
    for (j in i:ncol(matrix)) {
      if (!all.equal(matrix[i, j], matrix[j, i])) {
        print(paste("Matrix is not symmetric at positions:", i, j))
        return(FALSE)
      }
    }
  }
  print("Matrix is symmetric.")
  return(TRUE)
}

# Call the functions to check the matrix
check_empty_elements(vc_matrix)  # Check for empty elements
check_symmetry(vc_matrix)  # Check for symmetry

# Save V_matrix
save(vc_matrix, file = file.path(output_dir, "V_matrix_MPH.RData"))




rownames(I_matrix) <- NULL
colnames(I_matrix) <- NULL
dimnames(I_matrix) <- NULL


# Extract the lower triangle of N_matrix, including the diagonal
N_vector <- N_matrix[lower.tri(N_matrix, diag = TRUE)]

# Create a 1-row matrix from the extracted vector
N_single_row <- matrix(N_vector, nrow = 1)

# Remove dimension names to match LDSCoutput$N
dimnames(N_single_row) <- NULL

dimnames(vc_matrix) <- NULL

# Assign NULL to rownames and trait names to colnames
rownames(S_matrix) <- NULL
colnames(S_matrix) <- traits

mph_out <- list(
  V = vc_matrix,        # Variance-Covariance Matrix
  S = S_matrix,         # Genetic Covariance Matrix
  I = I_matrix,         # Cross-Trait Intercept Matrix
  N = N_single_row,     # Single-row vector of sample sizes
  m = m_value           # Number of SNPs from the .bim file
)


save(mph_out, file = file.path(output_dir, "MPH_genomicSEM.RData"))

# Print completion message
cat("Matrices generated and saved to:", output_dir, "\n")
