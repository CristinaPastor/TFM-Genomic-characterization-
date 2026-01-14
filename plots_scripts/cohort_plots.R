#!/usr/bin/env Rscript

# Cristina Pastor, October 2025

library(readr)
library(tidyverse)
library(ggplot2)

# CSV file
df <- read_csv(
    "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_seq.csv"
    )

# Output directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/cohort_plots/"
if(!dir.exists(output_dir)) dir.create(output_dir)


df <- df[!is.na(df$survival), ]
df <- df[df$sequenced_germline == "yes", ]


# Sex pie
df_sex <- df %>%
  group_by(sex) %>%
  summarise(num_sex = n(), .groups = "drop") 

color_sex <- c(
    "female" = "#9155c2",
    "male" = "#55c26f"
)

pie_sex <- ggplot(df_sex, aes(x = "", y = num_sex, fill = sex)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = paste0(num_sex, " (", round((num_sex / sum(num_sex)) * 100, 1), "%)")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_fill_manual(values = color_sex) +
  coord_polar(theta = "y") +
  labs(title = "Sex distribution", fill = "Sex") +
  theme_void() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    plot.margin = margin(40, 40, 40, 40),
    panel.spacing = unit(0, "mm")
  )

# Save PDF
pdf_sex <- file.path(output_dir, "pie_sex.pdf")
ggsave(pdf_sex, pie_sex, width = 7, height = 6, useDingbats = FALSE)
cat("PDF saved:", pdf_sex, "\n")


# Affected pie
df_affected <- df %>%
  group_by(affected) %>%
  summarise(num_affected = n(), .groups = "drop") 

color_affected <- c(
    "yes" = "#a5c255",
    "no" = "#c25578"
)

pie_affected <- ggplot(df_affected, aes(x = "", y = num_affected, fill = affected)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = paste0(num_affected, " (", round((num_affected / sum(num_affected)) * 100, 1), "%)")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_fill_manual(values = color_affected) +
  coord_polar(theta = "y") +
  labs(title = "Affected status distribution", fill = "Affected status") +
  theme_void() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    plot.margin = margin(40, 40, 40, 40),
    panel.spacing = unit(0, "mm")
  )

# Save PDF
pdf_affected <- file.path(output_dir, "pie_affected.pdf")
ggsave(pdf_affected, pie_affected, width = 7, height = 6, useDingbats = FALSE)
cat("PDF saved:", pdf_affected, "\n")


# Pie for distribution of families according to the number of affected individuals
df_family <- df %>%
  group_by(family_id) %>%
  summarise(
    family_case = if_else(
      any(family_case == "yes", na.rm = TRUE), 
      "Families with >1 affected", 
      "Families with 1 affected"
    ),
    .groups = "drop"
  ) %>%
  count(family_case, name = "num_family_case")

color_family_case <- c(
    "Families with >1 affected" = "#bdb543ff",
    "Families with 1 affected" = "#e682d0ff"
)

pie_family <- ggplot(df_family, aes(x = "", y = num_family_case, fill = family_case)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = paste0(num_family_case, " (", round((num_family_case / sum(num_family_case)) * 100, 1), "%)")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_fill_manual(
    values = color_family_case,
    labels = scales::wrap_format(25)
  ) +
  coord_polar(theta = "y") +
  labs(
    title = stringr::str_wrap(
      "Family distribution by number of affected individuals", 
      width = 30
    ),
    fill = "Family category"
  ) +
  theme_void() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.2, size = 20, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    legend.key.size = unit(0.9, "cm"),
    plot.margin = margin(40, 40, 40, 40),
    panel.spacing = unit(0, "mm")
  )

# Save PDF
pdf_family <- file.path(output_dir, "pie_family.pdf")
ggsave(pdf_family, pie_family, width = 7, height = 6, useDingbats = FALSE)
cat("PDF saved:", pdf_family, "\n")


# survival status distribution pie
df_survival <- df %>%
  group_by(survival) %>%
  summarise(num_survival = n(), .groups = "drop") 

color_survival <- c(
    "alive" = "#55c2a7",
    "dead" = "#c255a7"
)

pie_survival <- ggplot(df_survival, aes(x = "", y = num_survival, fill = survival)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = paste0(num_survival, " (", round((num_survival / sum(num_survival)) * 100, 1), "%)")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_fill_manual(values = color_survival) +
  coord_polar(theta = "y") +
  labs(title = "Survival distribution", fill = "Survival status") +
  theme_void() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    plot.margin = margin(40, 40, 40, 40),
    panel.spacing = unit(0, "mm")
  )

# Save PDF
pdf_survival <- file.path(output_dir, "pie_survival.pdf")
ggsave(pdf_survival, pie_survival, width = 7, height = 6, useDingbats = FALSE)
cat("PDF saved:", pdf_survival, "\n")

cat("\n##### PIPELINE COMPLETED #####\n")