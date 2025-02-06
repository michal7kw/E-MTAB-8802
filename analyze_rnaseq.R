#####################################################################
# RNA-Seq Differential Expression Analysis Pipeline
#####################################################################
# This script performs differential expression analysis on RNA-seq data
# using DESeq2. It processes HTSeq count files, performs statistical
# analysis, and generates various visualizations and result files.

# Create output directory if it doesn't exist
dir.create("./data", showWarnings = FALSE)
dir.create("./data/plots", showWarnings = FALSE)

# Install and load required Bioconductor packages
# BiocManager is the package manager for Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install core analysis packages:
# - DESeq2: Main package for differential expression analysis
# - pheatmap: For creating heatmaps
# - ggplot2: For creating plots
# - EnhancedVolcano: For creating enhanced volcano plots
# - org.Hs.eg.db: Human genome annotation database
BiocManager::install(c("DESeq2", "pheatmap", "ggplot2", "EnhancedVolcano", "org.Hs.eg.db"))

# Load all required libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(dplyr)

# Set working directory to data location
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/E-MTAB-8802")  # Changed from "./data" to "." since we'll use "./data" for output

# Get list of all HTSeq count files in the directory
count_files <- list.files("./data", pattern="*_HTSeq_formatted.count$", full.names=TRUE)

# Remove failed sample
count_files <- count_files[!grepl("3373-EN-12", count_files)]

#####################################################################
# Data Import and Preprocessing
#####################################################################

# Check if any count files were found
if (length(count_files) == 0) {
    stop("No HTSeq count files found in the current directory. 
          Files should end with '_HTSeq_formatted.count'")
}

# Function to read and process HTSeq count files
# Parameters:
#   file: Path to HTSeq count file
# Returns:
#   A data frame with gene counts and cleaned sample names
read_htseq <- function(file) {
    if (!file.exists(file)) {
        stop(sprintf("File not found: %s", file))
    }
    
    counts <- try(read.table(file, stringsAsFactors = FALSE))
    if (inherits(counts, "try-error")) {
        stop(sprintf("Error reading file: %s", file))
    }
    
    if (ncol(counts) != 2) {
        stop(sprintf("Invalid file format in %s: Expected 2 columns", file))
    }
    
    # Remove version numbers from ENSEMBL IDs (e.g., ENSG00000001.4 -> ENSG00000001)
    rownames(counts) <- gsub("\\.[0-9]+$", "", counts$V1)
    # Keep only count column
    counts <- counts[,2,drop=FALSE]
    # Clean up sample names by removing file extension
    colnames(counts) <- gsub("_HTSeq_formatted.count$", "", basename(file))
    return(counts)
}

# Read all count files with error handling
count_list <- list()
for (file in count_files) {
    tryCatch({
        count_list[[file]] <- read_htseq(file)
    }, error = function(e) {
        warning(sprintf("Error processing file %s: %s", file, e$message))
    })
}

if (length(count_list) == 0) {
    stop("No valid count data could be loaded")
}

# Combine into count matrix
countData <- do.call(cbind, count_list)

# Verify data
if (nrow(countData) == 0 || ncol(countData) == 0) {
    stop("Empty count matrix created")
}

# Print dimensions for verification
cat("Count matrix dimensions:", dim(countData), "\n")

# Extract sample order
actual_sample_order <- colnames(countData)
print("Sample order:")
print(actual_sample_order)

#####################################################################
# Experimental Design Setup
#####################################################################

# Create metadata for samples
# Assumes samples ending with -1, -2, or -3 are condition B, others are condition A
sample_conditions <- data.frame(
    sample = actual_sample_order,
    condition = ifelse(grepl("-(1|2|3)$", actual_sample_order), 
                      "NB",  # Changed from conditionB to NB
                      "DMEM")  # Changed from conditionA to DMEM
)

# Create sample information data frame for DESeq2
sampleInfo <- data.frame(
    condition = factor(sample_conditions$condition),
    row.names = sample_conditions$sample
)

# Print sample information for verification
print("Sample Information:")
print(sampleInfo)

#####################################################################
# DESeq2 Analysis
#####################################################################

# Create DESeq2 object with validation
tryCatch({
    dds <- DESeqDataSetFromMatrix(
        countData = countData,
        colData = sampleInfo,
        design = ~ condition
    )
}, error = function(e) {
    stop(sprintf("Error creating DESeq2 object: %s", e$message))
})

# Filter out genes with low counts (less than 10 reads total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform differential expression analysis
dds <- DESeq(dds)

# Save the DESeq2 object for further analysis
saveRDS(dds, "./data/dds_object.rds")

# Extract results and sort by adjusted p-value
res <- results(dds, contrast=c("condition", "DMEM", "NB"))
res <- res[order(res$padj),]

# Map ENSEMBL IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                      keys = rownames(res),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
res$gene_symbol <- gene_symbols

# Save normalized counts for plotting
normalized_counts <- counts(dds, normalized=TRUE)
# Clean up column names to match original sample names
colnames(normalized_counts) <- gsub("\\.", "-", colnames(normalized_counts))
write.csv(normalized_counts, "./data/normalized_counts.csv")

# Save raw counts
raw_counts <- counts(dds, normalized=FALSE)
# Clean up column names to match original sample names
colnames(raw_counts) <- gsub("\\.", "-", colnames(raw_counts))
write.csv(raw_counts, "./data/raw_counts.csv")

# Save complete results table with more detailed information
results_df <- as.data.frame(res) %>%
  mutate(
    ensembl_id = rownames(.),
    baseMean = round(baseMean, 2),
    log2FoldChange = round(log2FoldChange, 3),
    padj = round(padj, 6)
  )
write.csv(results_df, "./data/differential_expression_results.csv")

# Save significant genes (adjusted p-value < 0.05 and |log2FC| > 1)
sig_genes <- subset(results_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig_genes, "./data/significant_genes.csv")

# Save sample information
write.csv(sampleInfo, "./data/sample_info.csv", row.names=TRUE)

#####################################################################
# Generate Plots
#####################################################################

# Generate MA plot
tryCatch({
    png("./data/plots/MA_plot.png")
    plotMA(res, ylim=c(-2,2))
    dev.off()
}, error = function(e) {
    warning(sprintf("Error generating MA plot: %s", e$message))
})

# Create enhanced volcano plot
tryCatch({
    png("./data/plots/volcano_plot.png", width=1200, height=1000)
    EnhancedVolcano(res,
        lab = res$gene_symbol,
        x = 'log2FoldChange',
        y = 'padj',
        title = 'Differential Expression',
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 3.0,
        labSize = 6.0)
    dev.off()
}, error = function(e) {
    warning(sprintf("Error generating volcano plot: %s", e$message))
})

# Generate heatmap
tryCatch({
    topGenes <- head(order(res$padj), 50)
    mat <- counts(dds, normalized=TRUE)[topGenes,]
    rownames(mat) <- gene_symbols[rownames(mat)]
    png("./data/plots/heatmap.png", width=1200, height=1000)
    pheatmap(mat,
             scale="row",
             show_rownames=TRUE,
             main="Top 50 Differentially Expressed Genes")
    dev.off()
}, error = function(e) {
    warning(sprintf("Error generating heatmap: %s", e$message))
})

# Print summary statistics of the differential expression results
summary(res) 