# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Read table4 data using consistent options with compare_counts.R
table4 <- read.csv("./data/table4.csv", stringsAsFactors = FALSE)
cat("\nDebug: Original table4 data (first few rows):\n")
print(head(table4))

# Clean up ensemblid column by removing version numbers (e.g., '.16') if applicable
if ("ensemblid" %in% colnames(table4)) {
  table4 <- table4 %>%
    mutate(ensemblid = sub("\\..*", "", ensemblid))
  cat("\nDebug: Cleaned table4 data (first few rows):\n")
  print(head(table4))
}

# List of transporters we want to analyze
transporters <- c("ABCB1", "ABCC4", "SLC1A1", "SLC2A1", 
                 "SLC6A6", "SLC7A5", "SLC39A8")

# Filter for our transporters of interest
transporter_data <- table4 %>%
  filter(hgnc_symbol %in% transporters)

# Extract and reshape normalized count columns
normalized_counts <- transporter_data %>%
  select(hgnc_symbol, 
         contains("normalized")) %>%
  pivot_longer(
    cols = contains("normalized"),
    names_to = "sample",
    values_to = "count"
  ) %>%
  # Clean up sample names and add condition
  mutate(
    condition = case_when(
      grepl("DMEM", sample) ~ "DMEM",
      grepl("TF-NB", sample) ~ "NB"
    )
  ) %>%
  # Clean and convert condition to factor
  mutate(
    condition = trimws(condition),
    condition = toupper(condition),
    condition = factor(condition, levels = c("DMEM", "NB"))
  )

# Add debug prints
print("Debug: Checking condition values and factor levels:")
print(table(normalized_counts$condition))
print(levels(normalized_counts$condition))

# Create the plot
transporter_plot <- ggplot(normalized_counts, 
                          aes(x = condition, y = count, fill = condition)) +
  # Bar for mean value
  stat_summary(fun = mean, 
              geom = "bar", 
              position = position_dodge(width = 0.9)) +
  # Error bars for mean ± standard deviation
  stat_summary(fun.data = mean_sdl, 
              geom = "errorbar", 
              position = position_dodge(width = 0.9), 
              width = 0.2) +
  # Individual data points
  geom_point(position = position_dodge(width = 0.9), 
             size = 2, 
             color = "black") +
  # Facet by gene with free y-scales
  facet_wrap(~hgnc_symbol, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("DMEM" = "#404040", "NB" = "#CCCCCC")) +
  # Customize theme to match analysis plot
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0)
  ) +
  # Labels
  labs(
    y = "Normalized counts",
    x = "",
    title = "(c)"
  )

# Add debugging prints
print("Debug: Color mapping that will be used:")
print(c("DMEM" = "#404040", "NB" = "#CCCCCC"))

print("Debug: Checking if plot object contains correct color mapping:")
print(transporter_plot$scales$scales[[1]]$palette(2))

# Save the plot
ggsave("./data/plots/transporter_plot_from_table4.png", 
       transporter_plot, 
       width = 10, 
       height = 8, 
       dpi = 300)

# After saving, add:
print("Debug: Saving plot to:")
print(normalizePath("./data/plots/transporter_plot_from_table4.png"))
print("Debug: File saved. Checking if file exists:")
print(file.exists("./data/plots/transporter_plot_from_table4.png"))
print("Debug: File modification time:")
print(file.mtime("./data/plots/transporter_plot_from_table4.png"))