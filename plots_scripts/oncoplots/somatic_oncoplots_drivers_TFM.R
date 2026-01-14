#!/usr/bin/env Rscript

# Script for plotting oncoplots with selected driver genes
# Cristina Pastor, November 2025


library(maftools) 
library(data.table) 

# Files 
file_all_snvs <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_snvs_maf_tumors.maf"
file_all_snvs_VAF80 <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_snvs_80_maf_tumors.maf"
file_all_snvs_filtered <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_snvs_filtered_maf_tumors.maf" 
file_all_cnvs <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_cnvs_tumors.csv"
file_clinical <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/clinical_annotations.tsv" 



# read files
cat("\n--- Reading files...\n") 
snvs_maf <- read.maf(maf = file_all_snvs, verbose = FALSE) 
snvs_VAF80_maf <- read.maf(maf = file_all_snvs_VAF80, verbose = FALSE)
snvs_filtered_maf <- read.maf(maf = file_all_snvs_filtered, verbose = FALSE)
cnvs_file <- read.delim(file_all_cnvs, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
clinical_annotations <- fread(file_clinical, sep = "\t", header = TRUE, stringsAsFactors = FALSE)



# Subset only needed columns 
clinical_annotations_subset <- clinical_annotations[, .(Tumor_Sample_Barcode, tumor_primary, tumor_purity_range, unique_case)] 

# Change format of CNVs
cnvs_file$CN <- toupper(cnvs_file$CN)

cnvs_file$CN[cnvs_file$CN == "AMP"] <- "Amp"
cnvs_file$CN[cnvs_file$CN == "DUP"] <- "Amp"
cnvs_file$CN[cnvs_file$CN == "DEL,DUP"] <- "Amp"

cnvs_file$CN[cnvs_file$CN == "DEL"] <- "Del"


colnames(cnvs_file)[colnames(cnvs_file) == "Gene"] <- "Hugo_Symbol"
colnames(cnvs_file)[colnames(cnvs_file) == "Sample_name"] <- "Tumor_Sample_Barcode"
colnames(cnvs_file)[colnames(cnvs_file) == "CN"] <- "Variant_Classification"


# persolanized colors
mutation_colors <- c(
  "Missense_Mutation" = "#E69F00",
  "In_Frame_Del" = "#D55E00",
  "Frame_Shift_Del" = "#009E73",
  "Splice_Site" = "#CC79A7",
  "Frame_Shift_Ins" = "#F0E442",
  "Nonsense_Mutation" = "#56B4E9",
  "In_Frame_Ins" = "#0072B2",
  "Multi_Hit" = "#999999",
  "Amp" = "#246930", 
  "Del" = "#6e0b1b"
)

annotation_colors <- list(
  # Tumors
  "tumor_primary" = c(
  "malignant colorectal neoplasm" = "#377EB8",
  "malignant neoplasm of appendix" = "#7030A0",
  "malignant tumor of breast" = "#dd60a2ff", 
  "malignant tumor of kidney" = "#FFFF33", 
  "malignant tumor of lung" = "#006400",
  "malignant tumor of meninges" = "#A65628",
  "malignant tumor of ovary" = "#984EA3", 
  "malignant tumor of pancreas" = "#800000",
  "malignant tumor of prostate" = "#00FFFF", 
  "malignant tumor of skin" = "#2E8B57", 
  "malignant tumor of stomach" = "#8DA0CB", 
  "malignant tumor of testis" = "#E41A1C",
  "malignant tumor of thyroid gland" = "#FF7F00",
  "sarcoma" = "#666666", 
  "malignant neoplasm of liver" = "#008080",
  "malignant neoplasm of urinary bladder" = "#A6D854",
  "malignant tumor of ampulla of vater" = "#FFFAC8",
  "malignant tumor of mediastinum" = "#B2ABD2",
  "pheochromocytoma" = "#1d4a6bff"
  ),

  "unique_case" = c(
    "yes" = "#5FA15F",
    "no" = "#C75948"
  ),

  "tumor_purity_range" = c(
    "0-10" = "#FFE6CC",
    "11-20" = "#FEC999",
    "21-30" = "#FDB866",
    "31-40" = "#FCAB33",
    "41-50" = "#F89E1A",
    "51-60" = "#F49000",
    "61-70" = "#EB8400",
    "71-80" = "#E17800",
    "81-90" = "#D76B00",
    "91-100" = "#CC5500",
    "unknown" = "#B3B3B3"
  )
)


# 1. Oncoplot with only snvs and driver genes
###### 1.1 NO FILTERED SNVs
# Define output directory
output_dir <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/oncoplots_TFM/" 

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE) 


driver_genes <- c("BRCA1", "BRCA2", "ATM", "KRAS", "BRAF", "NRAS", "PIK3CA", "PALB2", "CHEK2", 
                "ERBB2", "APC", "EGFR", "ALK", "IDH1", "KIT", "PDGFRA", "TP53", "TMPRSS1-ERG", 
                "PTEN", "FLT3", "NPM1", "DNMT3A", "POLE", "VHL", "RET", "IDH2", "FLCN", "NF1", "HRAS")

# Make oncoplot 
cat("\n-----Making oncoplot with SNVs and driver genes...\n") 

# --Tumor grouped
# Open device 
pdf(file.path(output_dir, "somatic_oncoplot_drivers_tumor.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = snvs_maf,
  genes = driver_genes,
  clinicalFeatures = c("tumor_primary", "unique_case", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()



cat("\n-----Oncoplot with SNVs and driver genes saved in ", output_dir, "\n")


# --Unique case grouped
# Open device 
pdf(file.path(output_dir, "somatic_oncoplot_drivers_unique.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = snvs_maf,
  genes = driver_genes,
  clinicalFeatures = c("unique_case", "tumor_primary", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs and driver genes saved in ", output_dir, "\n")

###### 1.2 FILTERED SNVs WITH VAF > 80
# Define output directory
output_dir_VAF80 <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/oncoplots_VAF80_TFM/" 

if(!dir.exists(output_dir_VAF80)) dir.create(output_dir_VAF80, recursive = TRUE) 

# Make oncoplot 
cat("\n-----Making oncoplot with SNVs (VAF > 80) and driver genes ...\n") 

# --Tumor grouped
# Open device 
pdf(file.path(output_dir_VAF80, "somatic_oncoplot_VAF80_drivers_tumor.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = snvs_VAF80_maf,
  genes = driver_genes,
  clinicalFeatures = c("tumor_primary", "unique_case", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs (VAF > 80) and driver genes saved in ", output_dir_VAF80, "\n")


# --Unique case grouped
# Open device 
pdf(file.path(output_dir_VAF80, "somatic_oncoplot_VAF80_drivers_unique.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = snvs_VAF80_maf,
  genes = driver_genes,
  clinicalFeatures = c("unique_case", "tumor_primary", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs (VAF > 80) and driver genes saved in ", output_dir_VAF80, "\n")


###### 1.3 FILTERED SNVs WITH DD-score
# Define output directory
output_dir_filtered <- "/storage/scratch01/groups/bu/impact_vuscan/characterization_plots/oncoplots_filtered_TFM/" 

if(!dir.exists(output_dir_filtered)) dir.create(output_dir_filtered, recursive = TRUE) 

# Make oncoplot 
cat("\n-----Making oncoplot with SNVs (dd-score) and driver genes...\n") 

# --Tumor grouped
# Open device 
pdf(file.path(output_dir_filtered, "somatic_oncoplot_filtered_drivers_tumor.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = snvs_filtered_maf,
  genes = driver_genes,
  clinicalFeatures = c("tumor_primary", "unique_case", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs (dd-score) and driver genes saved in ", output_dir_filtered, "\n")

# --Unique case grouped
# Open device 
pdf(file.path(output_dir_filtered, "somatic_oncoplot_filtered_drivers_unique.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = snvs_filtered_maf,
  genes = driver_genes,
  clinicalFeatures = c("unique_case", "tumor_primary", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs (dd-score) and driver genes saved in ", output_dir_filtered, "\n")


# 2. Oncoplot with SNVs + CNVs and driver genes

###### 2.1 NO FILTERED SNVs
# IDs that have CNVs
cnv_samples <- unique(cnvs_file$Tumor_Sample_Barcode)
cat("\nNumber of cases with CNVs:", length(cnv_samples), "\n")


combined_snvs_cnvs <- read.maf(
  maf = "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/snvs_with_cnvs.maf",
  cnTable = cnvs_file,
  verbose = FALSE
)

# Make oncoplot 
cat("\n-----Making oncoplot with SNVs + CNVs and driver genes...\n") 

# --Tumor grouped
# Open device 
pdf(file.path(output_dir, "somatic_oncoplot_CNVs_SNVs_drivers_tumor.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = combined_snvs_cnvs,
  genes = driver_genes,
  clinicalFeatures = c("tumor_primary", "unique_case", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs + CNVs and driver genes saved in ", output_dir, "\n")

# --Unique case grouped
# Open device 
pdf(file.path(output_dir, "somatic_oncoplot_CNVs_SNVs_drivers_unique.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = combined_snvs_cnvs,
  genes = driver_genes,
  clinicalFeatures = c("unique_case", "tumor_primary", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs + CNVs and driver genes saved in ", output_dir, "\n")



###### 2.2 FILTERED SNVs WITH VAF > 80

combined_snvs_cnvs_80 <- read.maf(
  maf = "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/snvs_80_with_cnvs.maf",
  cnTable = cnvs_file,
  verbose = FALSE
)


# Make oncoplot 
cat("\n-----Making oncoplot with SNVs + CNVs (VAF > 80) and driver genes...\n") 

# --Tumor grouped
# Open device 
pdf(file.path(output_dir_VAF80, "somatic_oncoplot_CNVs_SNVs_VAF80_drivers_tumor.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = combined_snvs_cnvs_80,
  genes = driver_genes,
  clinicalFeatures = c("tumor_primary", "unique_case", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs + CNVs (VAF > 80) and driver genes saved in ", output_dir_VAF80, "\n")


# --Unique case grouped
# Open device 
pdf(file.path(output_dir_VAF80, "somatic_oncoplot_CNVs_SNVs_VAF80_drivers_unique.pdf"), width = 22, height = 14)

par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = combined_snvs_cnvs_80,
  genes = driver_genes,
  clinicalFeatures = c("unique_case", "tumor_primary", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
)


dev.off()

cat("\n-----Oncoplot with SNVs + CNVs (VAF > 80) and driver genes saved in ", output_dir_VAF80, "\n")


# 3. Oncoplot with CNVs
# Create dummy MAF (to represent the CNVs a MAF file is needed)
dummy_maf_df <- data.frame(
  Hugo_Symbol = "DUMMY",
  Chromosome = "1",
  Start_Position = 1,
  End_Position = 1,
  Strand = "+",
  Variant_Classification = "Missense_Mutation",
  Variant_Type = "SNP",
  Reference_Allele = "A",
  Tumor_Seq_Allele2 = "T",
  Tumor_Sample_Barcode = cnv_samples[1], # first sample with dummy mutation
  stringsAsFactors = FALSE
)

# Read MAF
dummy_maf_with_cnv <- read.maf(
  maf = dummy_maf_df,
  cnTable = cnvs_file,
  verbose = FALSE
)

# Remove DUMMY mutation
dummy_maf_with_cnv@data <- dummy_maf_with_cnv@data[
  dummy_maf_with_cnv@data$Hugo_Symbol != "DUMMY", 
]

cat("\nCases in the CNVs oncoplot:", length(unique(dummy_maf_with_cnv@data$Tumor_Sample_Barcode)), "\n")

# Make oncoplot 
cat("\n-----Making oncoplot with CNVs and driver genes...\n") 

# --Tumor grouped
pdf(file.path(output_dir, "somatic_oncoplot_CNVs_drivers_tumor.pdf"), width = 22, height = 14)
par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = dummy_maf_with_cnv,
  genes = driver_genes,
  clinicalFeatures = c("tumor_primary", "unique_case", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = FALSE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
  )

dev.off()
cat("\n-----Oncoplot with CNVs and driver genes saved in", output_dir, "\n")

# --Unique case grouped
pdf(file.path(output_dir, "somatic_oncoplot_CNVs_drivers_unique.pdf"), width = 22, height = 14)
par(
    mar = c(6, 10, 12, 6),
    oma = c(4, 0, 2, 0),
    cex = 1,
    cex.axis = 1.4,
    cex.lab  = 1.4
)

oncoplot(
  maf = dummy_maf_with_cnv,
  genes = driver_genes,
  clinicalFeatures = c("unique_case", "tumor_primary", "tumor_purity_range"),
  annotationDat = clinical_annotations_subset,
  sortByAnnotation = TRUE,
  draw_titv = FALSE,
  gene_mar = 12, # margin in the genes
  fontSize = 1, # gene names size
  legendFontSize = 1.4,
  titleFontSize = 1.8,
  annotationFontSize = 1.4,
  annotationColor = annotation_colors, # clinic info
  colors = mutation_colors # mutations
  )

dev.off()
cat("\n-----Oncoplot with CNVs and driver genes saved in", output_dir, "\n")

cat("\n##### PIPELINE COMPLETED #####\n")