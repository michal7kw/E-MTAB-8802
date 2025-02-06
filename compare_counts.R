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

# Read raw counts from DESeq2 analysis
raw_counts <- read.csv("./data/raw_counts.csv", row.names=1)

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
# Process raw counts data from both sources
#####################################################################

# Extract and reshape raw counts from DESeq2 analysis
raw_deseq2_data <- raw_counts[genes_of_interest,] %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_mapping, by = "ensembl_id") %>%
    pivot_longer(
        -c(ensembl_id, gene_symbol),
        names_to = "sample",
        values_to = "deseq2_count"
    ) %>%
    mutate(
        condition = ifelse(grepl("^NB", sample), "NB", "DMEM")
    )

# Process raw counts from table4
table4_raw_data <- table4 %>%
    dplyr::filter(ensemblid %in% genes_of_interest) %>%
    dplyr::select(ensemblid, hgnc_symbol, 
                  DMEM_10_raw, DMEM_11_raw, 
                  NB_1_raw, NB_2_raw, NB_3_raw) %>%
    pivot_longer(
        cols = ends_with("_raw"),
        names_to = "sample",
        values_to = "table4_count"
    ) %>%
    mutate(
        sample = sub("_raw$", "", sample),
        condition = ifelse(grepl("^NB", sample), "NB", "DMEM")
    ) %>%
    dplyr::rename(ensembl_id = ensemblid)

#####################################################################
# Process normalized counts data from both sources
#####################################################################

# Extract and reshape normalized counts from DESeq2 analysis
normalized_deseq2_data <- normalized_counts[genes_of_interest,] %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_mapping, by = "ensembl_id") %>%
    pivot_longer(
        -c(ensembl_id, gene_symbol),
        names_to = "sample",
        values_to = "deseq2_count"
    ) %>%
    mutate(
        condition = ifelse(grepl("^NB", sample), "NB", "DMEM")
    )

# Process normalized counts from table4
table4_normalized_data <- table4 %>%
    dplyr::filter(ensemblid %in% genes_of_interest) %>%
    dplyr::select(ensemblid, hgnc_symbol, 
                  DMEM_10_normalized, DMEM_11_normalized,
                  NB_1_normalized, NB_2_normalized, NB_3_normalized) %>%
    pivot_longer(
        cols = ends_with("_normalized"),
        names_to = "sample",
        values_to = "table4_count"
    ) %>%
    mutate(
        sample = sub("_normalized$", "", sample),
        condition = ifelse(grepl("^NB", sample), "NB", "DMEM")
    ) %>%
    dplyr::rename(ensembl_id = ensemblid)

#####################################################################
# Create comparison data for both raw and normalized counts
#####################################################################

# Function to prepare comparison data
prepare_comparison_data <- function(deseq2_data, table4_data) {
    deseq2_data %>%
        group_by(ensembl_id, gene_symbol, condition) %>%
        summarise(
            deseq2_mean = mean(deseq2_count),
            deseq2_sd = sd(deseq2_count),
            deseq2_n = n(),
            .groups = "drop"
        ) %>%
        left_join(
            table4_data %>%
                group_by(ensembl_id, condition) %>%
                summarise(
                    table4_mean = mean(table4_count),
                    table4_sd = sd(table4_count),
                    table4_n = n(),
                    .groups = "drop"
                ),
            by = c("ensembl_id", "condition")
        )
}

raw_comparison <- prepare_comparison_data(raw_deseq2_data, table4_raw_data)
normalized_comparison <- prepare_comparison_data(normalized_deseq2_data, table4_normalized_data)

#####################################################################
# Create plots
#####################################################################

# Function to create comparison plot
create_comparison_plot <- function(comparison_data, title_prefix) {
    plot_data <- comparison_data %>%
        pivot_longer(
            cols = c(deseq2_mean, table4_mean),
            names_to = "source",
            values_to = "count"
        ) %>%
        mutate(source = ifelse(source == "deseq2_mean", "HTSeq", "Table 4"))
    
    ggplot(plot_data, 
           aes(x = condition, y = count, fill = source)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_wrap(~paste0(gene_symbol, " (", ensembl_id, ")"), 
                  scales = "free_y") +
        scale_fill_manual(values = c("HTSeq" = "#1f77b4", "Table 4" = "#ff7f0e")) +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            legend.position = "top"
        ) +
        labs(
            x = "Condition",
            y = "Counts",
            title = paste(title_prefix, "counts between HTSeq and Table 4"),
            fill = "Data Source"
        )
}

# Create and save raw counts comparison plot
raw_plot <- create_comparison_plot(raw_comparison, "Comparison of raw")
ggsave("./data/plots/raw_count_comparison_plot.png", 
       raw_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Create and save normalized counts comparison plot
normalized_plot <- create_comparison_plot(normalized_comparison, "Comparison of normalized")
ggsave("./data/plots/normalized_count_comparison_plot.png", 
       normalized_plot, 
       width = 12, 
       height = 8, 
       dpi = 300) 