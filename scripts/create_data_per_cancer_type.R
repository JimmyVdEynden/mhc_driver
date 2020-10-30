# Script used to generate "data/data_per_cancer_type.rds"

# # Imports

library(tidyverse)

# # Data

# ## Mutation status and affinity

# Load affinity and mutation data
HLA_matrix <- readRDS("data/HLA_matrix.rds")
mut_matrix <- readRDS("data/mut_matrix.rds")

# Exclude the 3 columns (mutations) that have an NA value for all samples
idx_NA_col <- which(colSums(is.na(HLA_matrix)) == nrow(HLA_matrix))
mut_matrix <- mut_matrix[,-idx_NA_col]
HLA_matrix <- HLA_matrix[,-idx_NA_col]

# Remove rows (samples) that contain NA value
# 1540 patients have NA values for all mutations and are excluded
idx_NA <- which(rowSums(is.na(HLA_matrix)) > 0)
mut_matrix <- mut_matrix[-idx_NA,]
HLA_matrix <- HLA_matrix[-idx_NA,]

# ## Cancer types

# Generate a list that contains the corresponding cancer type for each row in the HLA / mutation matrices
cancer_types_all = readRDS("downloads/TCGA/clin/TCGA_cancer_id.rds")
cancer_types = cancer_types_all[
    # Extract the patient specific part from the TCGA barcode
    mut_matrix %>% rownames %>% str_sub(1, 12)
  ]

# Check whether "cancer_types" has the right dimensions
stopifnot(length(cancer_types) == nrow(mut_matrix))
stopifnot(length(cancer_types) == nrow(HLA_matrix))

# # Generate RDS

# Group mutation status and affinity data per cancer type
data_per_cancer_type = cancer_types %>%
  # Ignore the existing row names. We use the indices to match cancer types with mut_matrix and HLA_matrices 
  unname %>%
  # Convert to tibble. Add current row ids as a column.
  as_tibble_col("cancer_type") %>% rowid_to_column %>%
  # Remove samples where no cancer type is found
  filter(!is.na(cancer_type)) %>%
  # Get list with the original row numbers per cancer type
  group_by(cancer_type) %>% summarise(list(rowid), .groups = 'drop') %>% deframe %>%
  # Select rows for both matrices
  map(~list(mutations = mut_matrix[.,], affinities = HLA_matrix[.,]))

# Add pan cancer as a "cancer type"
data_per_cancer_type$PAN = list(mutations = mut_matrix, affinities = HLA_matrix)

saveRDS(data_per_cancer_type, "data/data_per_cancer_type.rds")
