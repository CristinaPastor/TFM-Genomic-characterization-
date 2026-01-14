#!/usr/bin/env Rscript

# Script for plotting VAF
# Usage: Rscript TMB_plots.R
# Cristina Pastor, December 2025

library(maftools) 


# file
file_all_snvs <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_snvs_maf_tumors.maf"
# read file
snvs_maf <- read.maf(maf = file_all_snvs, verbose = FALSE) 


# Define output directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/" 

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE) 


genes <- c(
            # driver genes
            "BRCA1", "BRCA2", "ATM", "KRAS", "BRAF", "NRAS", "PIK3CA", "PALB2", "CHEK2", 
            "ERBB2", "APC", "EGFR", "ALK", "IDH1", "KIT", "PDGFRA", "TP53", 
            "PTEN", "FLT3", "NPM1", "DNMT3A", "POLE", "VHL", "RET", "IDH2", "FLCN", "NF1", "HRAS",
            # CPGs
            "ARID1A", "RBM10", "ASXL1", "TP53", "ZFHX3", "NOTCH1", "PTPN14", "ARID3A", "TGFBR2",
            "KMT2C", "ARID1B", "CHD2", "HLA-A", "KMT2D", "KMT2A", "MN1", "LRP1B", "ZFP36L2",
            "LRP5", "PIK3CA", "IRS2", "HNF1A", "APC", "ATR", "MSH6"
            )

# filter unique genes
genes <- unique(genes)

# create VAF values: t_alt_count / t_depth
snvs_maf@data$VAF <- with(
  snvs_maf@data,
  t_alt_count / t_depth
)

VAF_somatic <- file.path(output_dir, "somatic_VAF.pdf") 


pdf(VAF_somatic, width = 6, height = 8)

# Graph plotVaf
plotVaf(
    maf = snvs_maf,
    vafCol = "VAF",
    genes = genes,
    keepGeneOrder = TRUE,
    orderByMedian = TRUE,
    showN = TRUE,
    flip = TRUE,
    height = 8,
    width = 6
)

dev.off()

cat("\n##### PIPELINE COMPLETED #####\n")