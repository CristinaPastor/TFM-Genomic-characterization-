#!/usr/bin/env Rscript

# Script for plotting clinical data
# Usage: Rscript tumor_type_seq_plot.R
# Cristina Pastor, October 2025

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# Read CSVs
data <- read_csv("/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_seq.csv")

# Output directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/tumor_type_seq_plots/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory created:", output_dir, "\n")
} else {
  cat("Directory exists:", output_dir, "\n")
}


data <- data[(data$sequenced_tumor == 'yes' & data$sequenced_germline == 'yes' & !is.na(data$survival)), ]


# Detect columns
tumor_cols <- grep("^type_", colnames(data), value = TRUE)
tmb_cols <- grep("^TMB_\\d+_maftools$", colnames(data), value = TRUE)
age_columns <- grep("^diagnosis_age_\\d+$", colnames(data), value = TRUE)
sex_column <- grep("^sex$", colnames(data), value = TRUE)
case_id_column <- grep("^case_id$", colnames(data), value = TRUE)

if (length(tmb_cols) == 0) stop("No valid TMB maftools columns found.")

# Extract numeric indexes
indices <- sub(".*_(\\d+)_maftools$", "\\1", tmb_cols)

# Join type_i and TMB_i
long_data_tumor_seq <- lapply(indices, function(i) {
  data %>%
      select(
        case_id = any_of(case_id_column),
        sex = any_of(sex_column),
        tumor_type = any_of(paste0("type_", i)),
        TMB = any_of(paste0("TMB_", i, "_maftools"))
      ) %>%
      mutate(
        TMB = suppressWarnings(as.numeric(gsub(",", ".", gsub("[^0-9.,]", "", TMB))))
      ) %>%
      filter(!is.na(TMB))
  }) %>%
  bind_rows() 
  
# Count tumors per type
tumor_sex_counts <- long_data_tumor_seq %>%
  count(tumor_type, sex, name = "n")

tumor_order <- tumor_sex_counts %>%
  group_by(tumor_type) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  pull(tumor_type)

tumor_sex_counts <- tumor_sex_counts %>%
  mutate(tumor_type = factor(tumor_type, levels = tumor_order))


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


color_sex <- c(
    "female" = "#9155c2",
    "male" = "#55c26f"
)


# bar distribution of primary tumors
gg_bar <- ggplot(tumor_sex_counts, aes(x = tumor_type, y = n, fill = sex)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(
    values = color_sex
  ) +
  labs(
    title = "Distribution of sequenced tumor types by sex", 
    x = "Tumor type",
    y = "Number of tumors",
    fill = "Sex"
  ) +
  geom_text(
    aes(label = n),
    position = position_dodge2(width = 0.8, preserve = "single"),
    vjust = -0.3,
    size = 8
  ) +
  theme_minimal(base_family = "Times") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
    axis.text.x = element_text(angle = 68, hjust = 1, size = 23),
    axis.text.y = element_text(size = 23),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.box.margin = margin(t = 20),
    plot.margin = margin(40, 40, 40, 40)
  )

# Save PDF
pdf_path <- file.path(output_dir, "bar_tumor_type_seq.pdf")
ggsave(pdf_path, gg_bar, width = 16, height = 13, useDingbats = FALSE)
cat("PDF saved:", pdf_path, "\n")





# Bar chart with sequenced tumors group by age interval
long_data_age <- lapply(indices, function(i) {
  data %>%
    select(
      tumor_type = any_of(paste0("type_", i)),
      age = any_of(paste0("diagnosis_age_", i)),
      TMB = any_of(paste0("TMB_", i, "_maftools"))
    ) %>%
    # mutate(age = suppressWarnings(as.numeric(age))) %>%
    filter(!is.na(tumor_type), !is.na(age), !is.na(TMB))
}) %>%
  bind_rows()


age_min <- floor(min(long_data_age$age, na.rm = TRUE) / 5) * 5
age_max <- ceiling(max(long_data_age$age, na.rm = TRUE) / 5) * 5

long_data_age <- long_data_age %>%
  mutate(
    age_group = cut(
      age,
      breaks = seq(age_min, age_max, by = 10), # age interval 10
      right = FALSE
    )
  )


df_bar <- long_data_age %>%
  count(age_group, tumor_type, name = "n") %>%
  mutate(age_group = droplevels(age_group))

df_bar_totals <- df_bar %>%
  group_by(age_group) %>%
  summarise(ntotal = sum(n), .groups = "drop")

gg_bar <- ggplot(df_bar, aes(x = age_group, y = n, fill = tumor_type)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = color_map,
    labels = scales::wrap_format(25)
  ) +
  # n inside bars
  geom_text(
    aes(label = ifelse(n > 1, n, "")), # use an emtpy label for n <= 1 (otherwise numbers are not in their place)
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  # n on top of bars
  geom_text(
    data = df_bar_totals,
    aes(x = age_group, y = ntotal, label = paste0("n = ", ntotal)),
    vjust = -0.5,
    size = 8,
    inherit.aes = FALSE,
    color = "grey30"
  ) +
  labs(
    title = "Sequenced tumor types by age group",
    x = "Age at diagnosis",
    y = "Number of tumors",
    fill = "Tumor type"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24),
    legend.position = "bottom",
    plot.margin = margin(40, 40, 40, 40)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)))



pdf_bar <- file.path(output_dir, "bar_tumor_type_by_age_group_10.pdf")
ggsave(pdf_bar, gg_bar, width = 21, height = 16, useDingbats = FALSE)
cat("PDF saved:", pdf_bar, "\n")


cat("\n##### PIPELINE COMPLETED #####\n")