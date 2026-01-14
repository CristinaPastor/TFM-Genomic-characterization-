# Script to calculate TMB with maftools from a MAF file (already filtered)
# Usage: Rscript calculate_TMB.R <input_filtered_maf_file> <output_TMB.txt> <captured_size_bp>
# Cristina Pastor, September 2025

# define arguments
args <- commandArgs(trailingOnly = TRUE)
maf_file <- args[1]
output_file <- args[2]
captured_size_bp <- as.numeric(args[3])

# Check inputs
if(!file.exists(maf_file)) stop("MAF file not found")
if(is.na(captured_size_bp) || captured_size_bp <= 0) stop("Invalid captured size")

library(maftools)

# change working directory to output directory to avoid Rplots.pdf being created elsewhere
setwd(dirname(output_file))

# read MAF and calculate TMB
maf <- read.maf(maf = maf_file)

captured_size_mb <- captured_size_bp / 1e6
tmb_results <- tmb(maf = maf, captureSize = captured_size_mb)

tmb_value <- tmb_results$total_perMB[1]
write(tmb_value, file = output_file)

# save full TMB table
full_table_file <- sub("\\.txt$", "_full.tsv", output_file)
write.table(tmb_results, full_table_file, sep = "\t", quote = FALSE, row.names = FALSE)

# save summary PDF
pdf_file <- sub("\\.txt$", "_summary.pdf", output_file)
pdf(pdf_file)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

cat("[INFO] TMB, full table, and summary PDF saved in", dirname(output_file), "\n")
