#!/usr/bin/env Rscript

# Script to detect mutually exclusive, co-occuring and altered gene-sets (or top) with maftools (somaticInteractions function)
# somaticInteractions: performs Pair-wise Fisher’s Exact test to detect mutually exclusive or co-occuring events
# It is corrected by mutation rate per kb (long genes have more probabilities to mutate)
# Usage: Rscript somatic_interactions.R
# Cristina Pastor, October 2025

library(maftools)
library(tidyverse)
library(data.table)

file_all_snv <- read.maf("/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_snvs_maf_tumors.maf")
bed_file <- fread("/storage/scratch01/groups/bu/dcarrero_common/databases/gencode.v44.annotation_reformat.bed")

# Gene length (necessary for correcting by gene length)
# Not necessary to remove the header, fread should detect it and remove it

# Calculate gene length
gene_length <- bed_file[ , .(Length = sum(Stop - Start)), by = Gene_Symbol]

gene.summary <- getGeneSummary(file_all_snv)

# Merge maf and gene length by gene name
gene.summary <- merge(gene.summary, gene_length,
                by.x = "Hugo_Symbol",
                by.y = "Gene_Symbol",
                all.x = TRUE) # keep all genes in MAF file (if there is no length it will be NA)

# Calculate mutation rate per kb
gene.summary$mut.rate.kb <- gene.summary$total / (gene.summary$Length / 1000)

# Somatic interactions
# Directories
som_int_output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/somatic_interactions/output_length_corrected/"
if (!dir.exists(som_int_output_dir)) dir.create(som_int_output_dir, recursive = TRUE)

plot_output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/somatic_interactions_length_corrected/"
if(!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)

# Select top 50, 40, 25 and 20 genes (higher mutation rate per kb)
top50_genes <- head(gene.summary[order(-gene.summary$mut.rate.kb), ]$Hugo_Symbol, 50)
top40_genes <- head(gene.summary[order(-gene.summary$mut.rate.kb), ]$Hugo_Symbol, 40)
top25_genes <- head(gene.summary[order(-gene.summary$mut.rate.kb), ]$Hugo_Symbol, 25)
top20_genes <- head(gene.summary[order(-gene.summary$mut.rate.kb), ]$Hugo_Symbol, 20)

# Analyze somatic interactions
cat("\n-----Analyzing somatic interactions (full)...\n")
som_int_res <- somaticInteractions(
    maf = file_all_snv,
    genes = top50_genes,
    pvalue = c(0.05, 0.01),
    returnAll = TRUE,
    plotPadj = TRUE
)
cat("\n-----Analysis complete.\n")

# Save complete table
write.table(som_int_res,
            file = file.path(som_int_output_dir, "somatic_interactions_stats_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n-----Full statistics table (top 50 genes) saved in", som_int_output_dir, "\n")

# Graph top20
som_int_file_20 <- file.path(plot_output_dir, "somatic_interactions_top20.pdf")
pdf(som_int_file_20, width = 9, height = 7)

somaticInteractions(
    maf = file_all_snv,
    genes = top20_genes,
    pvalue = c(0.05, 0.01),
    returnAll = FALSE,
    plotPadj = TRUE,
    fontSize = 1,
    leftMar = 8,
    topMar = 8
)
dev.off()
cat("\n-----Top 20 plot saved.\n")

# Graph top 40
som_int_file_40 <- file.path(plot_output_dir, "somatic_interactions_top40.pdf")
pdf(som_int_file_40, width = 12, height = 8)

somaticInteractions(
    maf = file_all_snv,
    genes = top40_genes,
    pvalue = c(0.05, 0.01),
    returnAll = FALSE,
    plotPadj = TRUE,
    fontSize = 1.5,
    leftMar = 8,
    topMar = 8
)
dev.off()
cat("\n-----Top 40 plot saved.\n")



# Save only significant pairs (pAdj < 0.05)
cat("\n-----Saving significant pairs (pAdj<0.05) table...\n")
sig_pairs <- som_int_res[som_int_res$pAdj < 0.05, ] %>%
  arrange(pAdj) # show more significant pairs at the top

write.table(sig_pairs, 
            file = file.path(som_int_output_dir, "somatic_interactions_significant_pairs_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n-----Significant pairs (pAdj<0.05) saved in directory", som_int_output_dir, "\n")

# create and save mutual exclusive and co-occuring tables
cat("\n-----Creating and saving mutual exclusive and co-occuring tables...\n")
mutexcl <- som_int_res[som_int_res$Event == "Mutually_Exclusive", ]
cooccur <- som_int_res[som_int_res$Event == "Co_Occurence", ]

write.table(mutexcl, file = file.path(som_int_output_dir, "mutually_exclusive_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(cooccur, file = file.path(som_int_output_dir, "co_occuring_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n-----Mutual exclusive and co-occuring tables saved in directory", som_int_output_dir, "\n")


cat("\n### PIPELINE FINISHED ###\n")