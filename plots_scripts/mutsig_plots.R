#!/usr/bin/env Rscript

# Script for plotting clinical data: bar chart for mutational signatures and primary tumor
# Usage: Rscript mutsig_plots.R
# Cristina Pastor, October 2025

library(tidyverse)

# CSV file
df <- read_csv("/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos.csv")

# Output directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/"
if(!dir.exists(output_dir)) dir.create(output_dir)

# Detect mutational signatures and tumor columns
tumor_cols <- grep("^type_", names(df), value = TRUE)
sig_cols <- grep("^mutational_signatures_", names(df), value = TRUE)

# Pivot longer the information
df_long <- df %>%
  pivot_longer(
    cols = matches("^(type_|mutational_signatures_)"),
    names_to = c(".value", "tumor_num"),
    names_pattern = "(type_|mutational_signatures_)(\\d+)"
  ) %>%
  filter(!is.na(type_) & type_ != "") %>%
  rename(tumor = type_, signatures = mutational_signatures_) %>%
  separate_rows(signatures, sep = ",") %>% # separate mutational signatures by ','
  mutate(signatures = trimws(signatures)) %>%   # remove white spaces
  filter(!is.na(signatures) & signatures != "")

# Count tumors per signature
df_plot <- df_long %>%
  group_by(signatures, tumor) %>%
  summarise(num_tumores = n(), .groups = "drop")

# Personalized colors per type cancer
color_map <- c(
  "malignant colorectal neoplasm" = "#377EB8",
  "malignant neoplasm of appendix" = "#7030A0",
  "malignant tumor of breast" = "#dd60a2ff", 
  "malignant tumor of kidney" = "#FFFF33", 
  "malignant tumor of lung" = "#006400",
  "malignant tumor of meninges" = "#A65628",
  "malignant tumor of ovary" = "#984EA3", 
  "malignant tumor of pancreas" = "#800000",
  "malignant tumor of prostate" = "#00FFFF", 
  "malignant tumor of skin" = "#2E8B57", 
  "malignant tumor of stomach" = "#8DA0CB", 
  "malignant tumor of testis" = "#E41A1C",
  "malignant tumor of thyroid gland" = "#FF7F00",
  "sarcoma" = "#666666", 
  "malignant neoplasm of liver" = "#008080",
  "malignant neoplasm of urinary bladder" = "#A6D854",
  "malignant tumor of ampulla of vater" = "#FFFAC8",
  "malignant tumor of mediastinum" = "#B2ABD2",
  "pheochromocytoma" = "#1d4a6bff"
)

# ggplot
gg_bar <- ggplot(df_plot, aes(x = signatures, y = num_tumores, fill = tumor)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = color_map,
    labels = scales::wrap_format(25)
  ) +
  labs(title = "Distribution of mutational signatures across tumor types", 
        x = "Mutational signatures", 
        y = "Number of cases", 
        fill = "Tumor type") +
  theme_minimal(base_family = "Times") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
    axis.text.x = element_text(angle = 68, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.box.margin = margin(t = 20),
    plot.margin = margin(40, 40, 40, 40),
    
  )


# Save PDF
pdf_path <- file.path(output_dir, "bar_mut_sig.pdf")
ggsave(pdf_path, gg_bar, width = 20, height = 15, useDingbats = FALSE)


cat("Bar chart saved:", pdf_path, "\n")

cat("\n##### PIPELINE COMPLETED #####\n")