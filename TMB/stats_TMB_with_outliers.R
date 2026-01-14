#!/usr/bin/env Rscript

# Script for statistics analyses of TMB
# Usage: Rscript stats_TMB.R
# Author: Cristina Pastor, October 2025

# Load libraries

library(dplyr)
library(tidyr)
library(stringr)
library(car)       # for leveneTest
library(rstatix)   
library(readr)
library(ggplot2)
library(dplyr)
library(purrr)
library(nortest)      # lillie.test, ad.test, cvm.test
library(ggpubr)       # qqplot
library(FSA)

# Read CSV
tmb <- read_csv("/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_TMBman.csv")

# Output directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/TMB/statistics_output_with_outliers/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Output directory ready:", output_dir, "\n")

# Data preparation
# Detect TMB columns (including unknown) and type columns
tmb_cols <- grep("TMB", names(tmb), value = TRUE)
type_cols <- grep("^type_", names(tmb), value = TRUE)

# Extract TMB IDs to match with type_ columns
tmb_ids <- unique(str_extract(tmb_cols, "\\d+"))

# Convert TMB columns to numeric
tmb <- tmb %>%
  mutate(across(
    all_of(tmb_cols),
    ~ as.numeric(str_replace_all(., ",", ".")) 
  ))
  

# Create a list to store results
long_list <- list()

for (id in tmb_ids) {
  # TMB columns and type column for this id
  type_col <- paste0("type_", id)
  cols <- c(
    paste0("TMB_", id),
    paste0("TMB_", id, "_maftools"),
    paste0("TMB_", id, "_mutscape")
  )

  # Filter the columns that exist in the dataframe
  cols <- cols[cols %in% names(tmb)]
  if (!(type_col %in% names(tmb))) next

  # Convert to long format for that set
  df_long <- tmb %>%
    select(case_id, all_of(type_col), all_of(cols)) %>%
    pivot_longer(
      cols = all_of(cols),
      names_to = "tmp_col",
      values_to = "TMB"
    ) %>%
    mutate(
      method = case_when(
        str_detect(tmp_col, "maftools") ~ "maftools",
        str_detect(tmp_col, "mutscape") ~ "mutscape",
        TRUE ~ "base"
      ),
      type_tumor = .data[[type_col]]
    ) %>%
    select(case_id, type_tumor, method, TMB)
  
  long_list[[id]] <- df_long
}

# Combine all long dataframes
tmb_long <- bind_rows(long_list)

# Statistics between TMB method
print("\n\nSTATISTICS OF TMB METHOD\n\n")

# 1. Normality test (Shapiro-Wilk, Kolmogorov–Smirnov, K-S corregido: Lilliefors y Anderson-Darling)
# filter methods with at least 3 non-NA values
multi_normality_tests <- function(df, var = "TMB", group_vars, min_n = 3) {
  

  safe_test <- function(x, test) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) < min_n || sd(x) == 0) return(NA_real_)
    
    tryCatch({
      switch(test,
             "shapiro" = shapiro.test(x)$p.value,
             "ks"      = ks.test(x, "pnorm", mean(x), sd(x))$p.value,
             "lillie"  = lillie.test(x)$p.value,
             "ad"      = ad.test(x)$p.value,
             NA_real_)
    }, error = function(e) NA_real_)
  }
  
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n        = sum(!is.na(.data[[var]])),
      shapiro_p = safe_test(.data[[var]], "shapiro"),
      ks_p      = safe_test(.data[[var]], "ks"),
      lillie_p  = safe_test(.data[[var]], "lillie"),
      ad_p      = safe_test(.data[[var]], "ad"),
      .groups = "drop"
    ) %>%
    mutate(
      decision_AD = case_when(
        n < min_n ~ "No calculable (n < 3)",
        ad_p > 0.05 ~ "H0 is not rejected --> Normal", # Anderson-Darling is used for final decision
        TRUE ~ "H0 is rejected --> No Normal"
      )
    )
}

# Normality per method
normality_method <- multi_normality_tests(tmb_long, var = 'TMB', group_vars = c('method'))

cat("\nNormality test:\n")
cat("Anderson-Darling test p-value:", normality_method$ad_p, "\n")

# 1.1 Plot histogram
# Graph distribution and normal curve per method
plot_normality <- function(df, var = "TMB", group = "method", output_dir = ".", filename_prefix = "TMB_normality_plot") {

  # Filter NAs
  df <- df %>% filter(!is.na(.data[[var]]), !is.na(.data[[group]]))
  
  # Asegurate that group is factor
  df[[group]] <- as.factor(df[[group]])

  # Calculate mean and sd by group
  df_stats <- df %>%
    group_by(.data[[group]]) %>%
    summarise(mean = mean(.data[[var]], na.rm = TRUE),
              sd = sd(.data[[var]], na.rm = TRUE),
              min_val = min(.data[[var]], na.rm = TRUE),
              max_val = max(.data[[var]], na.rm = TRUE),
              .groups = "drop")

  # Create a dataframe with the normal curve for each group
  df_norm <- df_stats %>%
    group_by(.data[[group]]) %>%
    do({
      g <- .
      x_vals <- seq(g$min_val, g$max_val, length.out = 200)
      y_vals <- dnorm(x_vals, mean = g$mean, sd = g$sd)
      data.frame(x = x_vals, y = y_vals) %>% # data frame with x and y
        mutate(!!group := g[[group]])  # add the dynamic group column with mutate
    }) %>% ungroup()

  # faceted histogram with density and normal curve
  p <- ggplot(df, aes(x = .data[[var]])) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "gray80", color = "black") +
    geom_density(color = "red", linewidth = 1) +
    geom_line(data = df_norm, aes(x = x, y = y), color = "blue", linewidth = 1, linetype = "dashed") +
    facet_wrap(as.formula(paste("~", group)), scales = "free") +
    labs(
      title = paste("Distribution of", var, "by", group, "\nRed line: real density | Blue line: theoretical normal"),
      x = var,
      y = "Density"
    ) +
    theme_minimal(base_size = 22)

  # save plot
  filename <- file.path(output_dir, paste0(filename_prefix, "_by_", group, ".pdf"))
  ggsave(filename, plot = p, width = 12, height = 6)
  message("Gráfico facetado guardado en: ", filename)

  return(p)
}


plot_normality(tmb_long, var = "TMB", group = "method", output_dir = output_dir)

# 1.2 Plot QQ-plot by method
plot_qq_method <- function(df, var = "TMB", group = "method", output_dir = ".", filename_prefix = "QQ_TMB_by_method") {
  df <- df %>% filter(!is.na(.data[[var]]), !is.na(.data[[group]]))
  df[[group]] <- as.factor(df[[group]])
  
  p <- ggplot(df, aes(sample = .data[[var]])) +
    stat_qq() +
    stat_qq_line(color = "blue", linetype = "dashed") +
    facet_wrap(as.formula(paste("~", group)), scales = "free") +
    labs(title = "TMB QQ-plot per method", x = "Theoretical quantiles", y = "Observed quantiles") +
    theme_minimal(base_size = 13)
  
  filename <- file.path(output_dir, paste0(filename_prefix, ".png"))
  ggsave(filename, plot = p, width = 12, height = 6)
  message("QQ-plot per method saved in: ", filename)
  return(p)
}

plot_qq_method(tmb_long, var = "TMB", group = "method", output_dir = output_dir)

# 2. Homogeneity of variances
levene_res_method <- leveneTest(TMB ~ method, data = tmb_long) # returns a dataframe with several columns
levene_p_method <- levene_res_method[1, "Pr(>F)"] # the p-value is the first row of the column Pr(>F)
bartlett_p_method <- bartlett.test(TMB ~ method, data = tmb_long)$p.value

levene_text_method <- ifelse(levene_p_method > 0.05,
                      'H0 is not rejected: there is homoscedasticity of variances between TMB methods',
                      'H0 is rejected: there is no homoscedasticity of variances between TMB methods')

bartlett_text_method <- ifelse(bartlett_p_method > 0.05,
                      'H0 is not rejected: there is homoscedasticity of variances between TMB methods',
                      'H0 is rejected: there is no homoscedasticity of variances between TMB methods')

cat("\nTest of homogeneity of variances:\n")
cat("Levene test p-value:", levene_p_method, "\n")
cat(levene_text_method, '\n')
cat("Bartlett test p-value:", bartlett_p_method, "\n")
cat(bartlett_text_method, '\n')

# 3. One-way ANOVA (if assumptions are met) or Kruskal-Wallis (if assumptions are not met)
if (all(normality_method$ad_p > 0.05) & levene_p_method > 0.05) {
  cat("\nData meets assumptions (normality and homocedasticity) -> One-way ANOVA\n")
  res_method <- aov(TMB ~ method, data = tmb_long)
  summary_res <- summary(res_method)
  test_used <- "ANOVA"
} else {
  cat("\nData does not meet assumptions (normality and/or homocedasticity) -> Kruskal–Wallis\n")
  res_method <- kruskal.test(TMB ~ method, data = tmb_long)
  summary_res <- res_method
  test_used <- "Kruskal–Wallis"
}

cat("\nTest result of the global test (", test_used, "):\n", sep = "")
print(summary_res)

# 4. Save results
results_file_method_txt <- file.path(output_dir, "stats_TMB_by_method.txt")
results_file_method_rds <- file.path(output_dir, "stats_TMB_by_method.RDS")

sink(results_file_method_txt)
cat("### STATISTICAL RESULTS COMPARING TMB METHODS ###\n\n")

cat("## Normality test (Shapiro-Wilk, KS, Lilliefors, Anderson-Darling)\n")
print(normality_method)
cat("\nNote: the final decision is based on Anderson–Darling (column 'decision_AD').\n\n")

cat("\nTest of homogeneity of variances:\n")
cat("Levene p-value:", levene_p_method, "-->", levene_text_method)
cat("\nBartlett p-value:", bartlett_p_method, "-->", bartlett_text_method, "\n")
cat("\nNote: the final decision is based on the Levene test\n\n")

cat("## Global test of comparison between methods\n")
cat("Type of test:", test_used, "\n\n")
print(summary_res)

# p-value interpretation
p_global <- if(test_used == "ANOVA") {
  summary_res[[1]][["Pr(>F)"]][1]  # ANOVA: first row, Pr(>F)
} else if(test_used == "Kruskal–Wallis") {
  summary_res$p.value               # Kruskal-Wallis: p.value
} else {
  NA
}

if (!is.na(p_global)) {
  if (p_global < 0.05) {
    cat("\nInterpretation: p-value =", round(p_global,4), 
        "-> There are significant differences between methods.\n")
  } else {
    cat("\nInterpretation: p-value =", round(p_global,4), 
        "-> There are no significant differences between methods.\n")
  }
}

sink()


saveRDS(list(
  normality = normality_method,
  levene_p = levene_p_method,
  bartlett_p = bartlett_p_method,
  levene_text = levene_text_method,
  bartlett_text = bartlett_text_method,
  test_used = test_used,
  summary = summary_res
), results_file_method_rds)

cat("\nResults saved in:\n",
    results_file_method_txt, "\n",
    results_file_method_rds, "\n")


# Statistics between TMB method and tumor type --> filtered, with at least 5 cases per group
print("\n\nSTATISTICS OF TMB METHOD AND TUMOR TYPE FILTERING (min 5 cases per group)\n\n")

tmb_long_filtered <- tmb_long %>%
  filter(!is.na(TMB)) %>%
  group_by(type_tumor, method) %>%
  filter(n() >= 5) %>%
  ungroup()

# 1. Normality test (Shapiro-Wilk, Kolmogorov–Smirnov, K-S corregido: Lilliefors y Anderson-Darling)
# 1.1 Primary tumor type-method
normality_tumor_filtered <- multi_normality_tests(df = tmb_long_filtered, var = 'TMB', group_vars = c('type_tumor', 'method'))


# 1.2 Plot histogram
# By tumor type and method pairs: one facet per tumor and different colors per method
plot_normality_tumor_color <- function(df, var = "TMB", output_dir = ".", filename_prefix = "TMB_normality_plot_by_tumor_filtered") {
  
  df <- df %>% 
    filter(!is.na(TMB), !is.na(type_tumor), !is.na(method))
  
  df_stats <- df %>%
    group_by(type_tumor, method) %>%
    summarise(
      mean = mean(.data[[var]], na.rm = TRUE),
      sd   = sd(.data[[var]], na.rm = TRUE),
      min_val = min(.data[[var]], na.rm = TRUE),
      max_val = max(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # dataframe with normal curves per tumor-method
  df_norm <- df_stats %>%
    group_by(type_tumor, method) %>%
    do({
      g <- .
      x_vals <- seq(g$min_val, g$max_val, length.out = 200)
      y_vals <- dnorm(x_vals, mean = g$mean, sd = g$sd)
      data.frame(x = x_vals, y = y_vals)
    })
  
  # Define linetypes for methods
  linetypes <- c("base" = "dashed", "maftools" = "solid", "mutscape" = "dotted")
  
  p <- ggplot(df, aes(x = .data[[var]], fill = method)) +
    geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 25, alpha = 0.3, color = "black") +
    geom_density(aes(y = after_stat(density)), color = "black", linewidth = 1.2) +  # observed density
    geom_line(
      data = df_norm,
      aes(x = x, y = y, linetype = method),
      color = "blue",
      linewidth = 0.9
    ) + # theoretical normal curves with different linetypes
    scale_linetype_manual(values = linetypes) +
    facet_wrap(~ type_tumor, scales = "free") +
    labs(
      title = "TMB distribution per type of tumor\nDensity observed (black) vs theoretical normal (blue lines)",
      x = var, y = "Density",
      fill = "Method", linetype = "Method"
    ) +
    theme_minimal(base_size = 22) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
  
  filename <- file.path(output_dir, paste0(filename_prefix, ".pdf"))
  ggsave(filename, plot = p, width = 14, height = 10)
  message("Plot saved in: ", filename)
  
  return(p)
}

plot_normality_tumor_color(tmb_long_filtered, output_dir = output_dir)

# 1.3 QQ-plot By tumor type and method pairs: one facet per tumor and different colors per method
plot_qq_tumor_method <- function(df, var = "TMB", output_dir = ".", filename_prefix = "QQ_TMB_by_tumor_method_filtered") {
  df <- df %>% filter(!is.na(TMB), !is.na(type_tumor), !is.na(method))
  
  p <- ggplot(df, aes(sample = TMB, color = method)) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") +
    facet_wrap(~ type_tumor, scales = "free") +
    labs(title = "TMB QQ-plot per tumor and method",
         x = "Theoretical quantiles", y = "Observed quantiles",
         color = "Method") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  
  filename <- file.path(output_dir, paste0(filename_prefix, ".png"))
  ggsave(filename, plot = p, width = 14, height = 10)
  message("QQ-plot per tumor and method saved in: ", filename)
  return(p)
}

plot_qq_tumor_method(tmb_long_filtered, output_dir = output_dir)

# 2. Homogeneity of variances (by method already done)
levene_res_tumor_filt <- leveneTest(TMB ~ type_tumor, data = tmb_long_filtered)
levene_p_tumor_filt <- levene_res_tumor_filt[1, "Pr(>F)"] 
bartlett_p_tumor_filt <- bartlett.test(TMB ~ type_tumor, data = tmb_long_filtered)$p.value

levene_text_tumor_filt <- ifelse(levene_p_tumor_filt > 0.05,
                      'H0 is not rejected: there is homoscedasticity of variances',
                      'H0 is rejected: there is no homoscedasticity of variances')

bartlett_text_tumor_filt <- ifelse(bartlett_p_tumor_filt > 0.05,
                      'H0 is not rejected: there is homoscedasticity of variances',
                      'H0 is rejected: there is no homoscedasticity of variances')

cat("\nHomogeneity test of variances:\n")
cat("Levene test p-value:", levene_p_tumor_filt, "\n")
cat(levene_text_tumor_filt, '\n')
cat("Bartlett test p-value:", bartlett_p_tumor_filt, "\n")
cat(bartlett_text_tumor_filt, '\n')

# 3. multiple ANOVA if the assumptions are met
if (all(normality_tumor_filtered$ad_p > 0.05) & levene_p_tumor_filt > 0.05) {
  cat("\nData follows normality and homoscedasticity -> Multifactorial ANOVA\n")
  anova_mult_res_filt <- aov(TMB ~ type_tumor * method, data = tmb_long_filtered)
  summary_anova_mult_filt <- summary(anova_mult_res_filt)
  test_used <- "multifactorial ANOVA"
} else {
  cat("\nData does not meet normality and/or homoscedasticity assumptions")
  summary_anova_mult_filt <- NA
  test_used <- "Not applicable"
}

cat("\nGlobal result of the test (", test_used, "):\n", sep = "")
print(summary_anova_mult_filt)

# 4. Save results
results_file_tumor_filt_txt <- file.path(output_dir, "stats_TMB_by_tumor_method_filtered.txt")
results_file_tumor_filt_rds <- file.path(output_dir, "stats_TMB_by_tumor_method_filtered.RDS")

sink(results_file_tumor_filt_txt)
cat("### STATISTIC RESULTS OF TMB BY TUMOR AND METHOD (only with >=5 cases per tumor-method) ###\n\n")

cat("## Normality Test (Shapiro-Wilk, KS, Lilliefors, Anderson-Darling)\n")
print(normality_tumor_filtered)
cat("\nNote: the final decision is based on Anderson–Darling (column 'decision_AD').\n\n")

cat("## Homogeneity test of variances\n")
cat("Levene p-value:", levene_p_tumor_filt, "-->", levene_text_tumor_filt, "\n")
cat("Bartlett p-value:", bartlett_p_tumor_filt, "-->", bartlett_text_tumor_filt, "\n")
cat("\nNote: the final decision is based on the Levene test\n\n")

cat("## Global test of comparison between tumors and methods\n")
cat("Type of test:", test_used, "\n\n")
if (!is.na(summary_anova_mult_filt)) print(summary_anova_mult_filt)
# p-value interpretation
if (!is.na(summary_anova_mult_filt) & test_used == "multifactorial ANOVA") {
  # Extract p-value from the summary table of the multifactorial aov
  # summary(aov) returns a list of data.frames; the first row is usually 'type_tumor'
  p_global <- summary_anova_mult_filt[[1]][["Pr(>F)"]][1] 
  
  if (!is.na(p_global)) {
    if (p_global < 0.05) {
      cat("\nInterpretation: p-value =", round(p_global,4), 
          "-> There are significant differences between the groups.\n")
    } else {
      cat("\nInterpretation: p-value =", round(p_global,4), 
          "-> There are no significant differences between the groups.\n")
    }
  }
}

sink()

saveRDS(list(
  normality = normality_tumor_filtered,
  levene_p = levene_p_tumor_filt,
  bartlett_p = bartlett_p_tumor_filt,
  levene_text = levene_text_tumor_filt,
  bartlett_text = bartlett_text_tumor_filt,
  test_used = test_used,
  summary = summary_anova_mult_filt
), results_file_tumor_filt_rds)

cat("\nResultados guardados en:\n",
    results_file_tumor_filt_txt, "\n",
    results_file_tumor_filt_rds, "\n")




# Kruskal-Wallis between pairs with >= 5 cases
# Multifactorial ANOVA cannot be performed, so I will check for differences between tumors within each method using Kruskal-Wallis
cat("\nStatistical analysis by Kruskal-Wallis between tumor-method pairs with >= 5 cases\n")

kruskal_results <- tmb_long_filtered %>%
  group_by(method) %>%
  summarise(
    tumors_included = paste(unique(type_tumor), collapse = ", "),
    statistic = kruskal.test(TMB ~ type_tumor, data = pick(everything()))$statistic,
    df = kruskal.test(TMB ~ type_tumor, data = pick(everything()))$parameter,
    p_kw = kruskal.test(TMB ~ type_tumor, data = pick(everything()))$p.value,
    decision = ifelse(
      kruskal.test(TMB ~ type_tumor, data = pick(everything()))$p.value < 0.05,
      "H0 is rejected: there are significant differences between tumors",
      "H0 is not rejected: there are no significant differences between tumors"
    )
  )


write.table(kruskal_results, file = file.path(output_dir, "KruskalWallis_pairs_filtered.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# POST-HOC
# As three methods are significant, we do post hoc to see which tumors within each method are different
# Filter the methods where there are significant differences
significant_methods <- kruskal_results %>%
  filter(p_kw < 0.05) %>%
  pull(method)

# List to save the Dunn results
dunn_list <- list()

for (m in significant_methods) {
  cat("\n\n### Dunn's test for method:", m, "###\n")
  
  df_method <- tmb_long_filtered %>% filter(method == m)
  
  # Dunn test with Bonferroni correction
  dunn_res <- dunnTest(TMB ~ type_tumor, data = df_method, method = "bonferroni")
  
  # Save the results in the list
  dunn_list[[m]] <- dunn_res$res
  
  # Show on screen
  print(dunn_res$res)
}

# Save combined results in a single data frame
dunn_combined <- bind_rows(dunn_list, .id = "method")

# p-value interpretation: if p.adj < 0.05 --> significant differences
dunn_combined <- dunn_combined %>%
  mutate(
    decision = ifelse(P.adj < 0.05,
                      "Significant difference between tumors",
                      "Not significant")
  )

# Save results
dunn_file <- file.path(output_dir, "Dunn_posthoc_pairs_filtered.txt")
write.table(dunn_combined, file = dunn_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResults of Dunn's test saved in:\n", dunn_file, "\n")

# Create dataframe summary

# Results of Dunn are in a list by method, e.g.:
# dunn_results$base, dunn_results$maftools, dunn_results$mutscape
# with columns: Comparison, Z, P.unadj, P.adj

# First generate statistical summary by method and tumor
cat("\nGenerating summary dataframe\n")
tmb_summary <- tmb_long_filtered %>%
  group_by(method, type_tumor) %>%
  summarise(
    n = n(),
    mean_TMB = mean(TMB, na.rm = TRUE),
    median_TMB = median(TMB, na.rm = TRUE),
    sd_TMB = sd(TMB, na.rm = TRUE),
    min_TMB = min(TMB, na.rm = TRUE),
    max_TMB = max(TMB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(method) %>%
  arrange(method, desc(median_TMB)) %>%
  mutate(rank_median = row_number()) %>%
  ungroup()

# add Dunn significance information
add_dunn_significance <- function(summary_df, dunn_df) {
  summary_df <- summary_df %>%
    rowwise() %>%
    mutate(
      significant_vs = paste(
        dunn_df$Comparison[dunn_df$P.adj < 0.05 & grepl(type_tumor, dunn_df$Comparison)],
        collapse = "; "
      )
    ) %>%
    ungroup()
  return(summary_df)
}

tmb_summary_extended <- tmb_summary %>%
  group_by(method) %>%
  group_modify(~ add_dunn_significance(.x, dunn_combined %>% filter(method == .y$method[1]))) %>%
  ungroup()

# Save summary dataframe
cat("\nSaving summary dataframe in", output_dir, "\n")
write.csv(tmb_summary_extended, file.path(output_dir, "summary_pairs_filtered.csv"), row.names = FALSE)

tmb_summary_extended

cat("\n\n##### PIPELINE COMPLETED #####\n\n")