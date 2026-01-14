#!/usr/bin/env Rscript

# Script for plotting clinical data
# Usage: Rscript seq_tumor_subtypes.R
# Cristina Pastor, November 2025

library(tidyverse)
library(ggplot2)

df <- read.csv("/home/vant/Documentos/TFM/impact_data/cleaned_casos_with_seq.csv")


df <- df[(df$sequenced_tumor == 'yes' & df$sequenced_germline == 'yes' & !is.na(df$survival)), ]

# Detect columns
tumor_cols <- grep("^type_", colnames(df), value = TRUE)
subtype_cols <- grep("^subtype_", colnames(df), value = TRUE)
case_id_column <- grep("^case_id$", colnames(df), value = TRUE)
tmb_cols <- grep("^TMB_\\d+_maftools$", colnames(df), value = TRUE)

# Get indices from column names
indices <- sub(".*_(\\d+)_maftools$", "\\1", tmb_cols)

rename_as_unknown <- c(
  "malignant tumor of breast",
  "malignant colorectal neoplasm",
  "malignant tumor of ovary",
  "malignant tumor of prostate",
  "malignant tumor of lung",
  "malignant neoplasm of liver",
  "malignant tumor of ampulla of vater",
  "pheochromocytoma"
)

# Join type_i and subtype_i columns into a long format (where TMB_column is not empty)
long_data <- lapply(indices, function(i) {
  df %>%
    select(
      case_id = all_of(case_id_column),
      tumor_type = any_of(paste0("type_", i)),
      tumor_subtype = any_of(paste0("subtype_", i)),
      TMB = any_of(paste0("TMB_", i, "_maftools"))
    ) %>%
    mutate(
      tumor_subtype = str_trim(tumor_subtype),
      tumor_subtype = na_if(tumor_subtype, ""),
      tumor_subtype = replace_na(tumor_subtype, "unknown"),
      tumor_subtype = ifelse(tumor_subtype %in% rename_as_unknown, "unknown", tumor_subtype)
    ) %>%
    filter(!is.na(TMB) & !is.na(tumor_type))
}) %>%
  bind_rows()

cat("\nTumor subtype counts:\n")
count(long_data, tumor_subtype, sort = TRUE)

# Unique tumor subtypes per tumor type
unique_long_data <- long_data %>%
  mutate(
    tumor_subtype = str_trim(tumor_subtype),
    tumor_subtype = na_if(tumor_subtype, "")
  ) %>%
  group_by(tumor_type) %>%
  summarise(
    subtypes = paste(sort(unique(replace_na(tumor_subtype, "unknown"))),
                    collapse = ", "),
    .groups = "drop"
  )

table_path <- "/home/vant/Documentos/TFM/impact_data/unique_tumor_subtypes_per_type.csv"
write.csv(unique_long_data, table_path, row.names = FALSE)

# Counts of tumor subtypes per tumor type
subtype_counts <- long_data %>%
  group_by(tumor_type, tumor_subtype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(tumor_type, desc(count))

# counts per tumor type
type_counts <- long_data %>%
  group_by(tumor_type) %>%
  summarise(count = n(), .groups = 'drop')

# subtype colors
subtype_colors <- c(

  # Breast (#dd60a2ff)
  "ductal carcinoma in situ of breast"           = "#B8447F",
  "infiltrating carcinoma with ductal and lobular features" = "#DD60A2",
  "infiltrating duct carcinoma of breast"        = "#F08BBE",
  "infiltrating lobular carcinoma of breast"     = "#A23363",
  "neuroendocrine carcinoma of breast"           = "#8C1D4F",
  "papillary carcinoma of breast"                = "#F5C4DE",

  # Colorectal (#377EB8)
  "familial multiple polyposis syndrome"        = "#1F5FA5",
  "malignant tumor of colon"                    = "#4F8EDC",
  "primary adenocarcinoma of colon"             = "#9BBCE5",
  
  # Ovary (#984EA3)
  "clear cell adenocarcinoma of ovary"           = "#6A2C91",
  "endometrioid carcinoma ovary"                 = "#8A4FA8",
  "primary high grade serous adenocarcinoma of ovary" = "#B07AC2",
  "primary low grade serous adenocarcinoma of ovary"  = "#C9A3D6",
  "primary mucinous adenocarcinoma of ovary"     = "#5A1F73",
  "teratoma of ovary"                            = "#D6B8E3",
  
  # Thyroid (#FF7F00)
  "papillary thyroid carcinoma"                  = "#CC6600",
  "papillary thyroid microcarcinoma"             = "#FF9F40",

  # Kidney (#FFFF33)
  "chromophobe renal cell carcinoma"             = "#CCCC00",
  "clear cell carcinoma of kidney"               = "#FFFF33",
  "papillary renal cell carcinoma"               = "#999900",

  # Skin (#2E8B57)
  "basal cell carcinoma of skin"                 = "#1F6B42",
  "malignant melanoma of skin"                   = "#4FA97A",

  # Pancreas (#800000)
  "pancreatic ductal adenocarcinoma"             = "#660000",
  "primary malignant neuroendocrine neoplasm of pancreas" = "#A33333",

  # Prostate (#00FFFF)
  "adenocarcinoma of prostate"                   = "#00BFBF",

  # Lung (#006400)
  "carcinoid tumor of lung"                      = "#004B00",
  "spindle cell sarcoma"                         = "#4F8F4F",

  # Sarcoma (#666666)
  "sarcoma"                                      = "#666666",

  # Appendix (#7030A0)
  "neuroendocrine tumor of appendix"            = "#5A2385",
  "primary mucinous adenocarcinoma of appendix" = "#9A6BC4",

  # Liver (#008080)
  "hepatocellular carcinoma"                    = "#006666",

  # Stomach (#8DA0CB)
  "adenocarcinoma of stomach"                    = "#6C7FA8",
  "diffuse adenocarcinoma of stomach"            = "#AAB8DB",

  # Bladder (#A6D854)
  "papillary carcinoma of bladder"              = "#A6D854",

  # Ampulla of Vater (#FFFAC8)

  # CNS (#A65628)
  "primary malignant meningioma"                 = "#A65628",

  # Testis (#E41A1C)
  "seminoma"                                     = "#E41A1C",

  # Adrenal (#1d4a6bff)

  "unknown"                                      = "#B0B0B0"
)

# Order by tumor type counts and order legend by subtypes
subtype_counts <- subtype_counts %>%
  mutate(
    tumor_subtype = factor(tumor_subtype, levels = names(subtype_colors))
  ) %>%
  arrange(tumor_subtype) %>%
  left_join(type_counts, by = "tumor_type", suffix = c("", "_total")) %>%
  mutate(
    tumor_type = fct_reorder(tumor_type, count_total, .desc = TRUE)
  ) %>%
  select(-count_total)


# Plotting
gg_bar <- ggplot(subtype_counts, aes(x = tumor_type, y = count, fill = tumor_subtype)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    labels = scales::wrap_format(23),
    values = subtype_colors
  ) +
  # n inside bars
  geom_text(
    aes(label = ifelse(count>1, count, "")), 
    position = position_stack(vjust = 0.5), 
    size = 6
  ) +
  # n on top of bars
  geom_text(
    data = type_counts,
    aes(x = tumor_type, y = count, label = count),
    vjust = -0.5,
    size = 8,
    inherit.aes = FALSE,
    color = "grey30"
  ) +
  labs(
    title = "Tumor subtypes per tumor type",
    x = "Tumor type",
    y = "Number of tumors",
    fill = "Tumor subtype"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
    axis.text.x = element_text(angle = 68, hjust = 1, size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24),
    legend.position = "bottom",
    plot.margin = margin(40, 40, 40, 40)
) + 
scale_y_continuous(expand = expansion(mult = c(0, 0.08)))

pdf_bar <- file.path("/home/vant/Documentos/TFM/impact_data/tumor_subtypes_per_type.pdf")
ggsave(pdf_bar, plot = gg_bar, width = 22, height = 25)

cat("\n##### PIPELINE COMPLETED #####\n")