# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)  # Add this for textGrob
library(Hmisc) # Add this for stat_summary

#####################################################################
# Load Results from CSV files
#####################################################################

# Load the complete differential expression results
diff_expr <- read.csv("./data/differential_expression_results.csv", row.names=1)

# Load normalized counts and sample information
norm_counts <- read.csv("./data/normalized_counts.csv", row.names=1)
sample_info <- read.csv("./data/sample_info.csv", row.names=1)

# Clean up sample names in normalized counts
colnames(norm_counts) <- gsub("^X", "", colnames(norm_counts))
colnames(norm_counts) <- gsub("\\.", "-", colnames(norm_counts))

# Remove failed sample if present
if("3373-EN-12" %in% colnames(norm_counts)) {
    norm_counts <- norm_counts[, !colnames(norm_counts) %in% "3373-EN-12"]
}
if("3373-EN-12" %in% rownames(sample_info)) {
    sample_info <- sample_info[!rownames(sample_info) %in% "3373-EN-12", ]
}

#####################################################################
# Create PCA Plot (Panel a)
#####################################################################

# Perform PCA on normalized counts
pca_result <- prcomp(t(log2(norm_counts + 1)), scale=TRUE)
pca_data <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2],
    condition = sample_info$condition
)

# Calculate variance
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x=PC1, y=PC2, color=condition)) +
    geom_point(size=4) +
    xlab(sprintf("PC1: %.1f%% variance", var_explained[1])) +
    ylab(sprintf("PC2: %.1f%% variance", var_explained[2])) +
    theme_bw() +
    scale_color_manual(values=c("DMEM"="blue", "NB"="red")) +
    theme(legend.position="top") +
    ggtitle("(a)")

#####################################################################
# Create Volcano Plot (Panel b)
#####################################################################

# Create volcano plot from differential expression results
volcano_plot <- ggplot(diff_expr, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=ifelse(padj < 0.05,
                               ifelse(abs(log2FoldChange) > 5,
                                      "FC > 5",
                                      "FC > 2"),
                               "Non-significant")),
               alpha=0.6, size=1) +
    scale_color_manual(values=c("red", "blue", "grey"),
                      name="Significance",
                      labels=c("FC > 5", "FC > 2", "Non-significant")) +
    theme_bw() +
    xlab("log2(fold change)") +
    ylab("-log10(FDR-adjusted p-value)") +
    theme(legend.position="top") +
    ggtitle("(b)")

#####################################################################
# Create Transporter Expression Bar Plots (Panel c)
#####################################################################

# List of transporters to analyze
transporters <- c("ABCB1", "ABCC4", "SLC1A1", "SLC2A1", 
                 "SLC6A6", "SLC7A5", "SLC39A8")

# Get normalized counts for transporters
transporter_data <- diff_expr %>%
    mutate(ensembl_id = rownames(diff_expr)) %>%
    filter(gene_symbol %in% transporters) %>%
    left_join(
        as.data.frame(norm_counts) %>%
            mutate(ensembl_id = rownames(norm_counts)) %>%
            pivot_longer(-ensembl_id,
                        names_to = "sample",
                        values_to = "count"),
        by = "ensembl_id"
    ) %>%
    left_join(
        data.frame(
            sample = rownames(sample_info),
            condition = sample_info$condition,
            stringsAsFactors = FALSE
        ),
        by = "sample"
    )

# Create bar plot for each transporter with mean ± SD
transporter_plot <- ggplot(transporter_data, 
                          aes(x=condition, y=count, fill=condition)) +
    # Bar for mean value
    stat_summary(fun = mean, geom = "bar", 
                position = position_dodge(width = 0.9),
                na.rm = TRUE) +
    # Error bars for mean ± standard deviation
    stat_summary(fun.data = mean_sdl, 
                geom = "errorbar", 
                position = position_dodge(width = 0.9), 
                width = 0.2,
                na.rm = TRUE) +
    # Individual data points
    geom_point(position = position_dodge(width = 0.9), 
               size = 2, 
               alpha = 0.6,
               na.rm = TRUE) +
    facet_wrap(~gene_symbol, scales="free_y", ncol=4) +
    scale_fill_manual(values=c("DMEM"="#E0E0E0", "NB"="#808080")) + # Lighter grey colors
    theme_bw() +
    theme(
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, size = 12),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = "top"
    ) +
    ylab("Normalized counts") +
    xlab("") +
    ggtitle("(c)")

# Add debugging print for final data
print("Structure of transporter_data:")
print(str(transporter_data))
print("Unique conditions in transporter_data:")
print(unique(transporter_data$condition))

#####################################################################
# Combine and Save Plots
#####################################################################

# Save individual plots
ggsave("./data/plots/pca_plot.png", pca_plot, width=6, height=5, dpi=300)
ggsave("./data/plots/volcano_plot.png", volcano_plot, width=6, height=5, dpi=300)
ggsave("./data/plots/transporter_plot.png", transporter_plot, width=12, height=8, dpi=300)

# Combine all plots into one figure
combined_plot <- grid.arrange(
    pca_plot, volcano_plot, transporter_plot,
    layout_matrix = rbind(c(1,2), c(3,3)),
    heights = c(1, 1.5),
    top = textGrob(
        "Basal media alters the transcriptome and transporter activity in CC3-derived brain microvascular endothelial cell-like cells",
        gp = gpar(fontsize = 12),
        just = "center"
    )
)

ggsave("./data/plots/combined_plot.png", combined_plot, width=12, height=10, dpi=300)

# Print summary statistics for transporters
transporter_stats <- transporter_data %>%
    group_by(gene_symbol, condition) %>%
    summarise(
        Mean = mean(count),
        SD = sd(count),
        N = n(),
        .groups = "drop"
    )
print(transporter_stats) 