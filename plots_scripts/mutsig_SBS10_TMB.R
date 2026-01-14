#!/usr/bin/env Rscript

# Script to test correlation between mutational signature SBS10 presence and higher TMB values
# Cristina Pastor, December 2025

library(tidyverse)
library(ggplot2)
library(ggpubr)

df <- read.csv("/home/vant/Documentos/TFM/impact_data/cleaned_casos_with_seq.csv")


df <- df[(df$sequenced_tumor == 'yes' & df$sequenced_germline == 'yes' & !is.na(df$survival)), ]

# Detect columns
case_id_column <- grep("^case_id$", colnames(df), value = TRUE)
tmb_cols <- grep("^TMB_\\d+_maftools$", colnames(df), value = TRUE)
mutsig_cols <- grep("^mutational_signatures_\\d+$", colnames(df), value = TRUE)

# Get indices from column names
indices <- sub(".*_(\\d+)_maftools$", "\\1", tmb_cols)

# Join mutational_signatures_i and TMB_i_maftools columns into a long format (where TMB_column is not empty)
df_long <- lapply(indices, function(i) {
  df %>%
    select(
      case_id = all_of(case_id_column),
      mutational_signatures = any_of(paste0("mutational_signatures_", i)),
      TMB = any_of(paste0("TMB_", i, "_maftools"))
    ) %>%
    filter(!is.na(TMB) & !is.na(mutational_signatures))
}) %>%
  bind_rows()

# Identify samples with SBS10 signature
df_long <- df_long %>%
  mutate(
    has_SBS10 = str_detect(mutational_signatures, "SBS10a|SBS10b")
  ) %>%
  mutate(
    SBS10_status = case_when(
      has_SBS10 ~ "SBS10 Present",
      TRUE ~ "SBS10 Absent"
    )
  )

# Plot TMB distribution by SBS10 status
ggplot(df_long, aes(x = SBS10_status, y = TMB, fill = SBS10_status)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(
    title = "TMB distribution by SBS10 signature status",
    x = "SBS10 signature status",
    y = "Tumor Mutational Burden (log10 scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 19, hjust = 0.5, margin = margin(b = 35)),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 17, margin = margin(t = 10)),
    axis.title.y = element_text(size = 17, margin = margin(r = 10)),
    legend.title = element_text(size = 17, margin = margin(b = 10)),
    legend.text = element_text(size = 15),
    legend.position = "none"
    )

# Save plot
output_path <- "/home/vant/Documentos/TFM/impact_data/"
ggsave(paste0(output_path, "TMB_by_SBS10_status.png"), width = 8, height = 7)

# Statistical test
library(car)
library(rstatix)
library(ggpubr)
library(nortest)

# Check normality (anderson-darling test and shapiro-wilk test)
normality_test <- df_long %>%
  group_by(SBS10_status) %>%
  summarise(
    anderson_darling = ad.test(TMB)$p.value,
    shapiro_wilk = shapiro.test(TMB)$p.value
  )

# print normality test results
normality_test %>%
  rowwise() %>%
  mutate(
    AD_normal = ifelse(anderson_darling >= 0.05, "Normal", "Not normal"),
    SW_normal = ifelse(shapiro_wilk >= 0.05, "Normal", "Not normal")
  ) %>%
  select(SBS10_status, AD_normal, SW_normal) %>%
  print()

# check homogeneity of variances (levene's test)
levene_test <- leveneTest(TMB ~ SBS10_status, data = df_long)

# print levene's test result
levene_p <- levene_test$`Pr(>F)`[1]

if (levene_p >= 0.05) {
  cat("\nVariances are homogeneous (Levene p =", round(levene_p, 4), ")\n")
} else {
  cat("\nVariances are NOT homogeneous (Levene p =", round(levene_p, 4), ")\n")
}

# since normality is not met, but variances are homogeneous, use wilcoxon test
wilcox_test <- wilcox.test(TMB ~ SBS10_status, data = df_long)

wilcox_p <- wilcox_test$p.value

cat("\nWilcoxon test p-value:", wilcox_p, "\n")


if (wilcox_p < 0.05) {
  cat("\nSignificant difference in TMB between SBS10 status groups (Wilcoxon p =", round(wilcox_p, 4), ")\n")
} else {
  cat("\nNo significant difference in TMB between SBS10 status groups (Wilcoxon p =", round(wilcox_p, 4), ")\n")
}

# save all statistical results to a text file
output_stats_file <- file.path(output_path, "TMB_SBS10_statistical_results.txt")
sink(output_stats_file)
cat("Normality Test Results:\n")
print(normality_test %>%
  rowwise() %>%
  mutate(
    AD_normal = ifelse(anderson_darling >= 0.05, "Normal", "Not normal"),
    SW_normal = ifelse(shapiro_wilk >= 0.05, "Normal", "Not normal")
  ) %>%
  select(SBS10_status, AD_normal, SW_normal))

cat("\nLevene's Test for Homogeneity of Variances:\n")
cat("Levene p-value:", round(levene_p, 4), "\n")
if (levene_p >= 0.05) {
  cat("Variances are homogeneous.\n")
} else {
  cat("Variances are NOT homogeneous.\n")
}

cat("\nWilcoxon Test Results:\n")
cat("Wilcoxon p-value:", wilcox_p, "\n")
if (wilcox_p < 0.05) {
  cat("Significant difference in TMB between SBS10 status groups.\n")
} else {
  cat("No significant difference in TMB between SBS10 status groups.\n")
}
sink()

cat("\n##### PIPELINE COMPLETED #####\n")