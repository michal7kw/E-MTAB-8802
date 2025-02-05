#####################################################################
# RNA-Seq Differential Expression Analysis Pipeline
#####################################################################
# This script performs differential expression analysis on RNA-seq data
# using DESeq2. It processes HTSeq count files, performs statistical
# analysis, and generates various visualizations and result files.

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

# Set working directory to data location
setwd("./data")

# Get list of all HTSeq count files in the directory
count_files <- list.files(pattern="*_HTSeq_formatted.count$")

#####################################################################
# Data Import and Preprocessing
#####################################################################

# Function to read and process HTSeq count files
# Parameters:
#   file: Path to HTSeq count file
# Returns:
#   A data frame with gene counts and cleaned sample names
read_htseq <- function(file) {
    counts <- read.table(file, stringsAsFactors = FALSE)
    # Remove version numbers from ENSEMBL IDs (e.g., ENSG00000001.4 -> ENSG00000001)
    rownames(counts) <- gsub("\\.[0-9]+$", "", counts$V1)
    # Keep only count column
    counts <- counts[,2,drop=FALSE]
    # Clean up sample names by removing file extension
    colnames(counts) <- gsub("_HTSeq_formatted.count$", "", basename(file))
    return(counts)
}

# Read all count files and combine them into a single matrix
count_list <- lapply(count_files, read_htseq)
countData <- do.call(cbind, count_list)

# Extract and print sample order for verification
actual_sample_order <- gsub("_HTSeq_formatted.count$", "", colnames(countData))
print("Actual sample order in countData:")
print(actual_sample_order)

#####################################################################
# Experimental Design Setup
#####################################################################

# Create metadata for samples
# Assumes samples ending with -1, -2, or -3 are condition B, others are condition A
sample_conditions <- data.frame(
    sample = actual_sample_order,
    condition = ifelse(grepl("-(1|2|3)$", actual_sample_order), 
                      "conditionB", 
                      "conditionA")
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

# Create DESeq2 object with count data and experimental design
dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = sampleInfo,
    design = ~ condition  # Formula specifying experimental design
)

# Filter out genes with low counts (less than 10 reads total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform differential expression analysis
dds <- DESeq(dds)

# Extract results and sort by adjusted p-value
res <- results(dds)
res <- res[order(res$padj),]

# Map ENSEMBL IDs to gene symbols for better interpretability
gene_symbols <- mapIds(org.Hs.eg.db,
                      keys = rownames(res),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
res$gene_symbol <- gene_symbols

#####################################################################
# Results Output and Visualization
#####################################################################

# Save complete results table
write.csv(as.data.frame(res), file="differential_expression_results.csv")

# Generate MA plot
# Shows log2 fold changes vs mean normalized counts
png("MA_plot.png")
plotMA(res, ylim=c(-2,2))
dev.off()

# Create enhanced volcano plot
# Shows statistical significance vs fold change
png("volcano_plot.png", width=1200, height=1000)
EnhancedVolcano(res,
    lab = res$gene_symbol,
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Differential Expression',
    pCutoff = 0.05,    # Significance threshold
    FCcutoff = 1,      # Fold change threshold
    pointSize = 3.0,
    labSize = 6.0)
dev.off()

# Generate heatmap of top 50 differentially expressed genes
topGenes <- head(order(res$padj), 50)
mat <- counts(dds, normalized=TRUE)[topGenes,]
rownames(mat) <- gene_symbols[rownames(mat)]
png("heatmap.png", width=1200, height=1000)
pheatmap(mat,
         scale="row",           # Scale by row to show relative changes
         show_rownames=TRUE,
         main="Top 50 Differentially Expressed Genes")
dev.off()

# Extract significantly differentially expressed genes
# Criteria: adjusted p-value < 0.05 and absolute log2 fold change > 1
sigGenes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sigGenes), file="significant_genes.csv")

# Print summary statistics of the differential expression results
summary(res) 