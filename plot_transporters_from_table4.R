# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Read the table4 data
table4 <- read.csv("./data/table4.csv", row.names=1)

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
  )

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
  # Use grey scale for DMEM (filled) and NB (empty)
  scale_fill_manual(values = c("DMEM" = "grey40", "NB" = "white")) +
  # Customize theme
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
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

# Print summary statistics
summary_stats <- normalized_counts %>%
  group_by(hgnc_symbol, condition) %>%
  summarise(
    mean = mean(count),
    sd = sd(count),
    n = n(),
    .groups = "drop"
  )
print(summary_stats)

# Save the plot
ggsave("./data/plots/transporter_plot_from_table4.png", 
       transporter_plot, 
       width = 10, 
       height = 8, 
       dpi = 300) 