#####################################################################
# Extract Read Counts by Condition
#####################################################################
# This script extracts and saves read counts grouped by experimental
# conditions from a DESeq2 analysis

# Load required libraries
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)

# Create output directory if it doesn't exist
dir.create("./data", showWarnings = FALSE)

# Load the saved DESeq2 object
dds <- readRDS("./data/dds_object.rds")

# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Get sample information
sample_info <- colData(dds)

# Get DMEM and NB sample names
dmem_samples <- rownames(sample_info)[sample_info$condition == "DMEM"]
nb_samples <- rownames(sample_info)[sample_info$condition == "NB"]

# Extract counts for each condition
dmem_counts <- normalized_counts[, dmem_samples, drop=FALSE]
nb_counts <- normalized_counts[, nb_samples, drop=FALSE]

# Calculate mean counts per condition
dmem_mean_counts <- rowMeans(dmem_counts)
nb_mean_counts <- rowMeans(nb_counts)

# Get gene symbols for ENSEMBL IDs
gene_symbols <- mapIds(org.Hs.eg.db,
                      keys = rownames(normalized_counts),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Create data frames for each condition with all samples and mean
dmem_results <- as.data.frame(dmem_counts) %>%
  mutate(
    ensembl_id = rownames(dmem_counts),
    gene_symbol = gene_symbols,
    mean_count = dmem_mean_counts
  )

nb_results <- as.data.frame(nb_counts) %>%
  mutate(
    ensembl_id = rownames(nb_counts),
    gene_symbol = gene_symbols,
    mean_count = nb_mean_counts
  )

# Save results to CSV files
write.csv(dmem_results, "./data/DMEM_condition_counts.csv", row.names = FALSE)
write.csv(nb_results, "./data/NB_condition_counts.csv", row.names = FALSE)

# Create a combined summary with mean counts for both conditions
combined_summary <- data.frame(
  ensembl_id = rownames(normalized_counts),
  gene_symbol = gene_symbols,
  DMEM_mean = dmem_mean_counts,
  NB_mean = nb_mean_counts
)

# Save combined summary
write.csv(combined_summary, "./data/condition_means_summary.csv", row.names = FALSE)

# Print summary information
cat("Files created:\n")
cat("1. DMEM_condition_counts.csv - Individual and mean counts for DMEM condition\n")
cat("2. NB_condition_counts.csv - Individual and mean counts for NB condition\n")
cat("3. condition_means_summary.csv - Summary of mean counts for both conditions\n")

cat("\nNumber of genes analyzed:", nrow(normalized_counts), "\n")
cat("Number of DMEM samples:", length(dmem_samples), "\n")
cat("Number of NB samples:", length(nb_samples), "\n") 