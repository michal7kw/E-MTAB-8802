# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(org.Hs.eg.db)  # Add this for gene symbol mapping
library(tibble)  # Add this for rownames_to_column function

#####################################################################
# Read normalized counts from both sources
#####################################################################

# Read normalized counts from DESeq2 analysis
normalized_counts <- read.csv("./data/normalized_counts.csv", row.names=1)

# Read table4 data
table4 <- read.csv("./data/table4.csv", stringsAsFactors=FALSE)

#####################################################################
# List of genes to compare
#####################################################################

genes_of_interest <- c(
    "ENSG00000085563",
    "ENSG00000125257",
    "ENSG00000106688",
    "ENSG00000117394",
    "ENSG00000131389",
    "ENSG00000103257",
    "ENSG00000138821"
)

#####################################################################
# Map Ensembl IDs to gene symbols
#####################################################################

# Map ENSEMBL IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                      keys = genes_of_interest,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Create a mapping dataframe
gene_mapping <- data.frame(
    ensembl_id = genes_of_interest,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE
)

# Print the mapping for verification
cat("\nGene ID to Symbol mapping:\n")
print(gene_mapping)

#####################################################################
# Process normalized counts data
#####################################################################

# Extract and reshape normalized counts for our genes
normalized_data <- normalized_counts[genes_of_interest,] %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_mapping, by = "ensembl_id") %>%
    pivot_longer(
        -c(ensembl_id, gene_symbol),
        names_to = "sample",
        values_to = "deseq2_count"
    ) %>%
    mutate(
        # Remove leading X (if any) and replace dots with hyphens
        sample = gsub("^X", "", sample),
        sample = gsub("\\.", "-", sample),
        # Assign condition: if sample ends in -1, -2, or -3 then NB, else DMEM
        condition = ifelse(grepl("-(1|2|3)$", sample), "NB", "DMEM")
    )

#####################################################################
# Process table4 data
#####################################################################

# Extract normalized counts from table4
normalized_cols <- grep("normalized", names(table4), value = TRUE)

table4_data <- table4 %>%
    dplyr::filter(ensemblid %in% genes_of_interest) %>%
    dplyr::select(ensemblid, hgnc_symbol, all_of(normalized_cols)) %>%
    pivot_longer(
        cols = all_of(normalized_cols),
        names_to = "sample",
        values_to = "table4_count"
    ) %>%
    mutate(
        # Clean up sample names:
        sample = gsub(" \\(normalized\\)$", "", sample),
        sample = gsub("\\.\\.\\.", "-", sample),
        sample = gsub("\\.", "-", sample),
        sample = gsub("DMEM_F12-", "", sample),
        sample = gsub("TF-NB-", "", sample),
        # Remove any trailing "-normalized" text that remains (e.g. "--normalized-")
        sample = sub("[-]+normalized[-]+$", "", sample),
        # Assign condition: if sample ends in -1, -2, or -3 then NB, else DMEM
        condition = ifelse(grepl("-(1|2|3)$", sample), "NB", "DMEM")
    ) %>%
    dplyr::rename(ensembl_id = ensemblid)

# Print debug information
cat("\nDESeq2 normalized samples and conditions:\n")
print(normalized_data %>% 
    dplyr::select(sample, condition) %>% 
    distinct())

cat("\nTable4 samples and conditions:\n")
print(table4_data %>% 
    dplyr::select(sample, condition) %>% 
    distinct())

#####################################################################
# Compare the data
#####################################################################

comparison_data <- normalized_data %>%
    group_by(ensembl_id, gene_symbol, condition) %>%
    summarise(
        deseq2_mean = mean(deseq2_count),
        deseq2_n = n(),
        .groups = "drop"
    ) %>%
    left_join(
        table4_data %>%
            group_by(ensembl_id, condition) %>%
            summarise(
                table4_mean = mean(table4_count),
                table4_n = n(),
                .groups = "drop"
            ),
        by = c("ensembl_id", "condition")
    )

# Print summary of the comparison data
cat("\nSummary of comparison data:\n")
print(comparison_data %>% 
    dplyr::select(gene_symbol, condition, deseq2_n, table4_n))

# Create comparison plot
comparison_plot <- ggplot(comparison_data, aes(x = deseq2_mean, y = table4_mean)) +
    geom_point(aes(color = condition), size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    facet_wrap(~paste0(gene_symbol, " (", ensembl_id, ")"), scales = "free") +
    scale_color_manual(values = c("DMEM" = "blue", "NB" = "red")) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "top"
    ) +
    labs(
        x = "DESeq2 normalized counts",
        y = "Table 4 normalized counts",
        title = "Comparison of normalized counts between DESeq2 and Table 4",
        color = "Condition"
    )

# Save the plot
ggsave("./data/plots/count_comparison_plot.png", 
       comparison_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Print correlation statistics
cat("\nCorrelation statistics between DESeq2 and Table 4 counts:\n")
correlations <- comparison_data %>%
    group_by(gene_symbol, ensembl_id) %>%
    summarise(
        correlation = cor(deseq2_mean, table4_mean),
        .groups = "drop"
    )
print(correlations) 