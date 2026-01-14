#!/usr/bin/env Rscript

# Script for merging all the MAF files from tumors with SNVs in one MAF file
# Usage: Rscript merge_all_snvs_maf_tumors.R
# Cristina Pastor, November 2025

library(maftools)

# path to families directories
families_dir <- "/storage/scratch01/groups/bu/impact_vuscan/vcf2maf/output/" # for maf tumor files


# Iterate each family in order to find all MAF files 
cat("\n-----Starting to iterate each family...\n")
maf_files <- list.files(
  path = families_dir,
  pattern = "\\.vep\\.maf$",
  recursive = TRUE,
  full.names = TRUE
)

# read all maf files
cat("\n-----Starting to read all tumor maf files...\n")

# maf_list is a list of objects MAF, one per tumor
maf_list <- lapply(maf_files, function(f) {
    print(paste("Reading:", f))
    read.maf(maf = f)
}) 


# Merge all maf files in one
cat("\n-----Merging all MAF files\n")
cohort_maf <- merge_mafs(maf_list)

# Save cohort_maf
# Define directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# Create path for MAF combined file
output_maf_path <- file.path(output_dir, "all_snvs_maf_tumors.maf")


# Save MAFs
write.table(
  cohort_maf@data,
  file = output_maf_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n-----All tumor MAF files merged and saved in:", output_maf_path, "\n")

cat("\n\n##### PIPELINE COMPLETED #####\n\n")
