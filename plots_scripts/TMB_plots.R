#!/usr/bin/env Rscript

# Script for plotting TMB
# Usage: Rscript TMB_plots.R
# Cristina Pastor, October 2025

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Read data
data <- read_csv("/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_seq.csv")

# Filter data for sequenced tumor and germline with survival information
data <- data %>%
  filter(sequenced_tumor == "yes" & sequenced_germline == 'yes' & !is.na(survival))

# Define and create output directory
if (!dir.exists("/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/TMB_plots")) {
  dir.create("/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/TMB_plots", recursive = TRUE)
}

output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/TMB_plots/"

# Reshape TMB data to long format
tmb_long <- data %>%
  select(
    starts_with("TMB_"),
    sequenced_tumor
  ) %>%
  mutate(across(matches("^TMB_"), ~ as.numeric(str_replace(as.character(.x), ",", ".")))) %>%
  pivot_longer(
    cols = matches("^TMB_[0-9]+(_mutscape|_maftools)?$"),
    names_to = "tmb_source",
    values_to = "TMB"
  ) %>%
  filter(!is.na(TMB)) %>%
  mutate(
    method = case_when(
      str_detect(tmb_source, "maftools") ~ "Maftools",
      str_detect(tmb_source, "mutscape") ~ "MutScape",
      TRUE ~ "Base"
    )
  )

# Filter for visualization: TMB values > 150 are not shown
tmb_plot <- tmb_long %>%
  filter(TMB <= 150)

# Colors
tmb_colors <- c(
  "Base" = "#1B9E77",
  "Maftools" = "#7570B3",
  "MutScape" = "#D95F02"
)

# Boxplot
boxplot_TMB_method <- ggplot(tmb_plot, aes(x = method, y = TMB, fill = method)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = tmb_colors) +
  labs(
    title = "Distribution of TMB by method",
    subtitle = "Kruskal-Wallis p = 0.055 (ns)",
    x = "TMB method",
    y = "TMB (mutations/Mb)",
    caption = "Outliers with TMB > 150 were excluded from visualization"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none",
    plot.margin = margin(40, 40, 40, 40),
    plot.caption = element_text(size = 16, hjust = 0, margin = margin(t = 10))
  )

# Save PDF
output_file <- file.path(output_dir, "TMB_boxplot_distribution_by_method.pdf")

ggsave(filename = output_file,
  plot = boxplot_TMB_method,
  width = 8,
  height = 6,
  useDingbats = FALSE
)

cat("PDF saved in:", output_dir, "TMB_boxplot_distribution_by_method.pdf\n")




# Barplot for TMB distribution per tumor type and method
# Transform TMB data to long format
tmb_list <- list()

# Max number of tumors per individual
max_tumors <- 7

for(i in 1:max_tumors){
  type_col <- paste0("type_", i)
  tmb_base <- paste0("TMB_", i)
  tmb_mutscape <- paste0("TMB_", i, "_mutscape")
  tmb_maftools <- paste0("TMB_", i, "_maftools")
  
  if(all(c(type_col, tmb_base, tmb_mutscape, tmb_maftools) %in% colnames(data))){
    tmp <- data %>%
      select(all_of(c(type_col, tmb_base, tmb_mutscape, tmb_maftools))) %>%
      filter(!is.na(.data[[type_col]])) %>%
      mutate(
        Base = as.numeric(str_replace(.data[[tmb_base]], ",", ".")),
        MutScape = as.numeric(str_replace(.data[[tmb_mutscape]], ",", ".")),
        Maftools = as.numeric(str_replace(.data[[tmb_maftools]], ",", "."))
      ) %>%
      select(type = all_of(type_col), Base, MutScape, Maftools) %>%
      pivot_longer(cols = c(Base, MutScape, Maftools), names_to = "method", values_to = "TMB") %>%
      filter(!is.na(TMB) & TMB <= 150)
    
    tmb_list[[i]] <- tmp
  }
}

# Combine all tumors
tmb_long <- bind_rows(tmb_list)


# Boxplot
boxplot_TMB <- ggplot(
  tmb_long, 
  aes(x = reorder(type, TMB, FUN = median, na.rm = TRUE), y = TMB, fill = method)
  ) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.9, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = tmb_colors) +
  labs(
    title = "Distribution of TMB by tumor type and method",
    x = "Tumor type",
    y = "TMB (mutations/Mb)",
    caption = "Values TMB > 150 not shown"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 68, hjust = 1, size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    plot.caption = element_text(size = 16, hjust = 0),
    plot.margin = margin(40, 40, 40, 40)
  )

# Save PDF
output_file <- file.path(output_dir, "TMB_boxplot_by_tumor_method.pdf")

ggsave(output_file, boxplot_TMB, width = 10, height = 12, units = "in", useDingbats = FALSE)

cat("PDF saved in:", output_file, "\n")
cat("##### PIPELINE FINISHED #####\n")

