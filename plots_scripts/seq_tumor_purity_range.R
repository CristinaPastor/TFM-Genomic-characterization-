#!/usr/bin/env Rscript

# Script for plotting clinical data
# Usage: Rscript seq_tumor_purity_range.R
# Cristina Pastor, November 2025


library(tidyverse)
library(ggplot2)

df <- read.csv("/home/vant/Documentos/TFM/impact_data/cleaned_casos_with_seq.csv")


df <- df[(df$sequenced_tumor == 'yes' & df$sequenced_germline == 'yes' & !is.na(df$survival)), ]

# Detect columns
tumor_cols <- grep("^type_", colnames(df), value = TRUE)
tumor_purity_cols <- grep("^tumor_purity_\\d+$", colnames(df), value = TRUE)
case_id_column <- grep("^case_id$", colnames(df), value = TRUE)
tmb_cols <- grep("^TMB_\\d+_maftools$", colnames(df), value = TRUE)

# Get indices from column names
indices <- sub(".*_(\\d+)_maftools$", "\\1", tmb_cols)

# Join type_i and tumor_purity_i columns into a long format (where TMB_column is not empty)
df_long <- lapply(indices, function(i) {
  df %>%
    select(
      case_id = all_of(case_id_column),
      tumor_type = any_of(paste0("type_", i)),
      tumor_purity = any_of(paste0("tumor_purity_", i)),
      TMB = any_of(paste0("TMB_", i, "_maftools"))
    ) %>%
    filter(!is.na(TMB) & !is.na(tumor_type)) %>%
    # make ranges of tumor purity (by 10) and if there is no tumor purity value, categorize in "Unknown"
    mutate(
      tumor_purity_clean = str_squish(tumor_purity),
      tumor_purity_clean = str_replace_all(tumor_purity_clean, "%", ""),
      tumor_purity_clean = str_replace_all(tumor_purity_clean, ",", "."),
      tumor_purity_clean = as.numeric(tumor_purity_clean),
      tumor_purity_range = case_when(
        is.na(tumor_purity_clean) ~ "Unknown",
        tumor_purity_clean < 10 ~ "<10%",
        tumor_purity_clean >= 10 & tumor_purity_clean < 20 ~ "10-19%",
        tumor_purity_clean >= 20 & tumor_purity_clean < 30 ~ "20-29%",
        tumor_purity_clean >= 30 & tumor_purity_clean < 40 ~ "30-39%",
        tumor_purity_clean >= 40 & tumor_purity_clean < 50 ~ "40-49%",
        tumor_purity_clean >= 50 & tumor_purity_clean < 60 ~ "50-59%",
        tumor_purity_clean >= 60 & tumor_purity_clean < 70 ~ "60-69%",
        tumor_purity_clean >= 70 & tumor_purity_clean < 80 ~ "70-79%",
        tumor_purity_clean >= 80 & tumor_purity_clean < 90 ~ "80-89%",
        tumor_purity_clean >= 90 ~ "90-100%"
      )
    )
}) %>%
  bind_rows()

# Calculate mean and median tumor purity (excluding NAs)
purity_stats <- df_long %>%
  summarise(
    mean_tumor_purity = mean(tumor_purity_clean, na.rm = TRUE),
    median_tumor_purity = median(tumor_purity_clean, na.rm = TRUE),
    n_with_purity = sum(!is.na(tumor_purity_clean)),
    n_total = n()
  )

cat("\nTumor purity statistics:\n")
cat("Mean tumor purity:", round(purity_stats$mean_tumor_purity, 2), "\n")
cat("Median tumor purity:", round(purity_stats$median_tumor_purity, 2), "\n")
cat("Samples with purity:", purity_stats$n_with_purity, "/", purity_stats$n_total, "\n")


# Count how many tumors per tumor purity range
purity_counts <- df_long %>%
  group_by(tumor_purity_range) %>%
  summarise(n = n(), .groups = 'drop')

# Count how many tumors per tumor purity range and tumor type
purity_type_counts <- df_long %>%
  group_by(tumor_purity_range, tumor_type) %>%
  summarise(n = n(), .groups = 'drop')

# Color per tumor type
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

# Plotting
gg_purity <- ggplot(purity_type_counts, aes(x = tumor_purity_range, y = n, fill = tumor_type)) +
    geom_col(position = "stack") +
    scale_fill_manual(
        values = color_map,
        labels = scales::wrap_format(25)
    ) + 
    # n inside bars
    geom_text(
        aes(label = ifelse(n > 1, n, "")), 
        position = position_stack(vjust = 0.5), 
        size = 6
    ) +
    # n on top of bars
    geom_text(
        data = purity_counts,
        aes(x = tumor_purity_range, y = n, label = n),
        vjust = -0.5,
        size = 8,
        inherit.aes = FALSE,
        color = "grey30"
    ) +
    labs(
        title = "Tumor purity ranges per tumor type",
        x = "Tumor purity range",
        y = "Number of tumors",
        fill = "Tumor type"
    ) +
    theme_bw() +
    theme(
        text = element_text(family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
        axis.text.x = element_text(size = 24, margin = margin(t = 10)),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 26, margin = margin(t = 10)),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 24),
        legend.position = "bottom",
        legend.margin = margin(30, 0, 0, 0),
        plot.margin = margin(40, 40, 40, 40)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)))

pdf_purity <- file.path("/home/vant/Documentos/TFM/impact_data/tumor_purity_ranges_per_type.pdf")
ggsave(pdf_purity, plot = gg_purity, width = 21, height = 22)
cat("\n##### PIPELINE COMPLETED #####\n")