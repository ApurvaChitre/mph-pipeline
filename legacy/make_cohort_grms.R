# Read the cohort data
cohort_data <- read.csv("indicator_covariates.csv", stringsAsFactors = FALSE)

# Read the matrix
#####CHANGE TO MPH DIRECTORY TO READ IN THE AUXILLIARY SCRIPTS OR DOWNLOAD FROM HERE: https://jiang18.github.io/mph/scripts/
source("/home/wolftech/jjiang26/projects/ans590/demo/mph_functs.R")
matrix_data <- read_grm("gg")
matrix_data = matrix_data + t(matrix_data)
diag(matrix_data) = diag(matrix_data)/2

# Function to process matrix based on cohorts
process_matrix <- function(matrix, cohort_data, cohort1, cohort2 = NULL) {
  result <- matrix
  
  if (is.null(cohort2)) {
    # Case 1: Keep elements corresponding to cohort1 individuals
    cohort1_ids <- cohort_data$ID[cohort_data[[cohort1]] == 1]
    mask <- outer(rownames(result) %in% cohort1_ids, colnames(result) %in% cohort1_ids)
    result[!mask] <- 0
  } else {
    # Case 2: Keep elements for cohort1-cohort2 interactions
    cohort1_ids <- cohort_data$ID[cohort_data[[cohort1]] == 1]
    cohort2_ids <- cohort_data$ID[cohort_data[[cohort2]] == 1]
    mask1 <- outer(rownames(result) %in% cohort1_ids, colnames(result) %in% cohort2_ids)
    mask2 <- outer(rownames(result) %in% cohort2_ids, colnames(result) %in% cohort1_ids)
    result[!(mask1 | mask2)] <- 0
  }
  
  return(result)
}

# Generate all combinations of cohorts
cohorts <- paste0("cohort", 1:6)
cohort_combinations <- combn(cohorts, 2, simplify = FALSE)

# Case 1: Single cohort matrices
for (cohort in cohorts) {
  result <- process_matrix(matrix_data, cohort_data, cohort)
  write_grm(paste0("grm/", cohort), colnames(result), result)
  
}

# Case 2: Cohort pair matrices
for (combo in cohort_combinations) {
  combo_name <- paste(combo, collapse = "_")
  result <- process_matrix(matrix_data, cohort_data, combo[1], combo[2])
  write_grm(paste0("grm/", combo_name), colnames(result), result)

}

# Function to create diagonal matrix for a cohort
create_diagonal_matrix <- function(matrix, cohort_data, cohort) {
  result <- matrix * 0  # Create a zero matrix of the same size
  cohort_ids <- cohort_data$ID[cohort_data[[cohort]] == 1]
  diag_indices <- which(rownames(result) %in% cohort_ids)
  result[cbind(diag_indices, diag_indices)] <- 1
  return(result)
}

# Case 3: Diagonal matrices for each cohort
for (cohort in cohorts[1:5]) {  # We only need 5 matrices as per the request
  diagonal_name <- paste0(cohort, "_diagonal")
  result <- create_diagonal_matrix(matrix_data, cohort_data, cohort)
  write_grm(paste0("grm/", diagonal_name), colnames(result), result)

}
