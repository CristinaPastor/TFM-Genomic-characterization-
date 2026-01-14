#!/usr/bin/env Rscript

# Script for seeing if samples with higher diagnosis age have a higher proportion of clock-like signatures (SBS1 and SBS5)
# Cristina Pastor, December 2025


library(tidyverse)


df <- read.csv("/home/vant/Documentos/TFM/impact_data/cleaned_casos_with_seq.csv")


df <- df[(df$sequenced_tumor == 'yes' & df$sequenced_germline == 'yes' & !is.na(df$survival)), ]

# Detect columns
case_id_column <- grep("^case_id$", colnames(df), value = TRUE)
tmb_cols <- grep("^TMB_\\d+_maftools$", colnames(df), value = TRUE)
mutsig_cols <- grep("^mutational_signatures_\\d+$", colnames(df), value = TRUE)
diagnosis_age_cols <- grep("^diagnosis_age_\\d+$", colnames(df), value = TRUE)

# Get indices from column names
indices <- sub(".*_(\\d+)_maftools$", "\\1", tmb_cols)

# # Reshape data to long format
df_long <- df %>%
  pivot_longer(
    cols = matches("(mutational_signatures|diagnosis_age)_\\d+"),
    names_to = c(".value", "tumor_id"),
    names_pattern = "(.*)_(\\d+)"
  ) %>%
  filter(!is.na(diagnosis_age) & !is.na(mutational_signatures))

df_clock <- df_long %>%
  mutate(
    has_SBS1 = str_detect(mutational_signatures, "SBS1"), 
    has_SBS5 = str_detect(mutational_signatures, "SBS5"),
    age_group = cut(diagnosis_age,
                    breaks = c(0,30,40,50,60,70,80,90),
                    right = FALSE)
  )

# calculate proportions of SBS1 and SBS5 per age group
prop_counts <- df_clock %>%
  group_by(age_group) %>%
  summarise(
    SBS1 = mean(has_SBS1),
    SBS5 = mean(has_SBS5)
  ) %>%
  pivot_longer(cols = c(SBS1, SBS5),
               names_to = "clock_like",
               values_to = "proportion")
# Plotting
gg_plot <- ggplot(prop_counts, aes(x = age_group, y = proportion, fill = clock_like)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of clock-like mutational signatures by diagnosis age",
    x = "Diagnosis age group",
    y = "Proportion of tumors",
    fill = "Clock-like signature"
  ) +
  theme(
    axis.text.x = element_text(size = 18, margin = margin(b = 20)),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.title.x = element_text(size = 20, margin = margin(t = 20)),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  theme_minimal()

output_path <- "/home/vant/Documentos/TFM/impact_data/"
ggsave(paste0(output_path, "clock_like_sigs_by_age_group.png"), gg_plot, width = 8, height = 5)


cat("\n##### PIPELINE COMPLETED #####\n")