library(data.table)

# Define paths
pheno_dir <- "/tscc/projects/ps-palmer/apurva/locomotor/individual_projects_seq/data/pheno"
pheno_cohort_file <- "/tscc/projects/ps-palmer/apurva/locomotor/mph/individual_projects_seq/pheno_cohort_project_dict.csv"
output_pheno_file <- "/tscc/projects/ps-palmer/apurva/locomotor/mph/individual_projects_seq/combined_phenotype.csv"
output_cov_file <- "/tscc/projects/ps-palmer/apurva/locomotor/mph/individual_projects_seq/indicator_covariates.csv"

# Read the pheno cohort info
pheno_cohort_info <- fread(pheno_cohort_file)

# Initialize data tables for phenotype and covariate files
combined_pheno <- data.table(ID = character(), trait_value = numeric())
combined_cov <- data.table(ID = character())

# Initialize columns for combined_cov
for (i in 1:nrow(pheno_cohort_info)) {
  cohort <- pheno_cohort_info$cohort_num[i]
  combined_cov[, (cohort) := 0]
}

# Process each phenotype file
for (i in 1:nrow(pheno_cohort_info)) {
  # Read the phenotype file
  pheno_file <- file.path(pheno_dir, pheno_cohort_info$pheno_file[i])
  pheno_data <- fread(pheno_file, col.names = c("ID", "FID", "trait_value"))
  
  # Remove rows with NA in the trait_value column
  pheno_data <- pheno_data[!is.na(trait_value)]
  
  # Add cohort-specific indicator
  cohort <- pheno_cohort_info$cohort_num[i]
  pheno_data[, (cohort) := 1]
  
  # Merge with combined phenotype data
  combined_pheno <- rbind(combined_pheno, pheno_data[, .(ID, trait_value)])
  
  # Prepare cohort-specific covariate data
  cohort_cov <- pheno_data[, .(ID)]
  cohort_cov[, (cohort) := 1]
  
  # Merge with combined covariate data
  combined_cov <- merge(combined_cov, cohort_cov, by = "ID", all = TRUE)
}

# Replace NAs in covariate file with 0
combined_cov[is.na(combined_cov)] <- 0

# Merge cohort columns with .x and .y suffixes
cohort_cols <- grep("^cohort", names(combined_cov), value = TRUE)
combined_cov <- combined_cov[, lapply(.SD, max, na.rm = TRUE), by = ID, .SDcols = cohort_cols]




# Identify columns that end with ".y"
y_cols <- grep("\\.y$", names(combined_cov), value = TRUE)

# Subset to keep only these columns, including the ID column if necessary
combined_cov <- combined_cov[, c("ID", y_cols), with = FALSE]

# Print the resulting data.table to check the columns
print(combined_cov)


# Remove the ".y" suffix from column names
setnames(combined_cov, old = names(combined_cov), new = gsub("\\.y$", "", names(combined_cov)))

# Print the resulting data.table to check the column names
print(combined_cov)



# Write the combined phenotype and covariate files
fwrite(combined_pheno, output_pheno_file, sep = ",", col.names = TRUE)
fwrite(combined_cov, output_cov_file, sep = ",", col.names = TRUE)
