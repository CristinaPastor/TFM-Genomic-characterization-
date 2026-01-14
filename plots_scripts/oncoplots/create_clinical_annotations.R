#!/usr/bin/env Rscript

# Script for creating clinical annotations df for oncoplots
# Usage: Rscript create_clinical_annotations.R
# Cristina Pastor, November 2025

library(data.table)
library(maftools)
library(tidyverse)

# Create clinical data: bands of clinical annotations
# load csv clinical data
clinical <- read.csv("/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_suspected.csv")
setDT(clinical)

# Change column name family_case to unique_case
setnames(clinical, "family_case", "unique_case")
# Change values 'no' to 'yes' and viceversa
clinical <- clinical %>%
  mutate(
    unique_case = case_when(
      unique_case == "no" ~ "yes",
      unique_case == "yes" ~ "no",
      TRUE ~ unique_case
    )
  )

# Automatic cases: individuals with just one tumor (one suspected_sporadic_* column filled)
suspected_cols <- grep("^suspected_sporadic_[0-9]+$", names(clinical), value = TRUE)


# Detect tumor_purity_* columns
purity_cols <- grep("^tumor_purity_[0-9]+$", names(clinical), value = TRUE)

# Clean: remove % and convert to numeric (tumor_purity format is for example 80%)
clinical[, (purity_cols) := lapply(.SD, function(x) {
  x <- gsub("%", "", x)      # remove %
  x <- trimws(x)
  x[x == ""] <- NA           # empty to NA
  as.numeric(x)
}), .SDcols = purity_cols]


auto_cases <- clinical[
  , {
    filled <- suspected_cols[sapply(.SD, function(x) !is.na(x) & x != "")]

    # If there is only one suspected_col filled
    if (length(filled) == 1) {
      col <- filled

      # Asign tumor_suffix by default
      suffix <- "02"
      # Exceptions
      if (case_id == "1591-1591-4impact") suffix <- "06"
      if (case_id == "1871-1591-4impact") suffix <- "04"
      

      list(
        case_id = case_id,        
        unique_case = unique_case,
        suspected_col = col,
        tumor_suffix = suffix,
        suspicion_sporadic = get(col),
        tumor_purity = get(sub("suspected_sporadic_", "tumor_purity_", col)) / 100, # values will be 0.80
        tumor_purity_percent = get(sub("suspected_sporadic_", "tumor_purity_", col)), # the values will be 80
        tumor_primary = get(sub("suspected_sporadic_", "type_", col))
      )
    }
  },
  .SDcols = suspected_cols,
  by = case_id
]



# Create auxiliar table for manual cases:
manual_map <- data.table(
  case_id = c(
    "1837-1837-4impact", "1837-1837-4impact",
    "1884-1884-4impact", "1884-1884-4impact",
    "2560-2560-4impact", "2560-2560-4impact",
    "2748-2748-4impact", "2748-2748-4impact",
    "3700-3700-4impact", "3700-3700-4impact"
  ),
  suspected_col = c(
    "suspected_sporadic_2", "suspected_sporadic_1",
    "suspected_sporadic_1", "suspected_sporadic_2",
    "suspected_sporadic_1", "suspected_sporadic_2",
    "suspected_sporadic_1", "suspected_sporadic_2",
    "suspected_sporadic_1", "suspected_sporadic_2"
  ),
  tumor_suffix = c(
    "02", "04",
    "04", "02",
    "04", "02",
    "04", "02",
    "04", "02"
  ),
  unique_case = c(
    "no", "no",
    "no", "no",
    "no", "no",
    "no", "no",
    "no", "no"
  ),
  tumor_purity = c(
    0.9, 0.3,
    0.95, 0.75,
    0.3, 0.3,
    0.8, 0.7,
    "unknown", "unknown"
  ),
  tumor_purity_percent = c(
    90, 30,
    95, 75,
    30, 30,
    80, 70,
    "unknown", "unknown"
  )
)

# Extract suspicion and type using indicated column
manual_map[, suspicion_sporadic := mapply(function(case, col)
  clinical[case_id == case, get(col)], case_id, suspected_col)]

manual_map[, tumor_primary := mapply(function(case, col) {
  if (grepl("unknown", col)) {
    return("unknown")
  } else {
    return(clinical[case_id == case, get(sub("suspected_sporadic", "type", col))])
  }
}, case_id, suspected_col)]



# Check case_id is character in both tables
auto_cases[, case_id := as.character(case_id)]
manual_map[, case_id := as.character(case_id)]


head(auto_cases[, .(case_id, tumor_suffix, unique_case, tumor_purity, tumor_purity_percent)])
head(manual_map[, .(case_id, tumor_suffix, unique_case, tumor_purity, tumor_purity_percent)])


# Combine automatic and manual
combined_clin <- rbindlist(list(auto_cases, manual_map), fill = TRUE)


combined_clin[, purity_num := suppressWarnings(as.numeric(tumor_purity_percent))]

combined_clin[!is.na(purity_num),
              tumor_purity_range :=
                cut(purity_num,
                    breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                    include.lowest = TRUE,
                    right = TRUE,
                    labels=c("0-10","11-20","21-30","31-40","41-50",
                             "51-60","61-70","71-80","81-90","91-100"))]

combined_clin[is.na(purity_num),
              tumor_purity_range := "unknown"]


# If multiple columns related to case_id exist (e.g., i.case_id or case_id.1), coalesce them
cid_cols <- grep("^case_id(\\.|$|_)?", names(combined_clin), value = TRUE)
if (length(cid_cols) > 1) {
  # Create a definitive case_id by taking the first non-NA and non-empty value per row
  combined_clin[, case_id := {
    vals <- .SD
    # combine across columns: take the first non-empty value in order
    apply(vals, 1, function(r) {
      r <- as.character(r)
      r <- r[!is.na(r) & r != ""]
      if (length(r) == 0) return(NA_character_) else return(r[1])
    })
  }, .SDcols = cid_cols]
  # remove old columns except for the newly constructed 'case_id'
  cols_to_remove <- setdiff(cid_cols, "case_id")
  if (length(cols_to_remove) > 0) combined_clin[, (cols_to_remove) := NULL]
}


# If the same case_id appears in both tables, prioritize the manually curated version
# (this does not remove rows with different tumor_suffix)
combined_clin <- combined_clin[
  order(fifelse(case_id %in% manual_map$case_id, 2, 1)),  # prioridad a manuales
  .SD[.N],
  by = .(case_id, tumor_suffix)
]

# Remove potential columns created by rbind (e.g., i.case_id)
if ("i.case_id" %in% names(combined_clin)) combined_clin[, i.case_id := NULL]
# also remove columns like case_id.1 if they exist
caseid_like <- grep("^case_id[\\.0-9]+$", names(combined_clin), value = TRUE)
if (length(caseid_like) > 0) combined_clin[, (caseid_like) := NULL]

# Create Tumor_Sample_Barcode
combined_clin[, Tumor_Sample_Barcode := paste0(case_id, "-", tumor_suffix)]


combined_clin <- combined_clin[combined_clin$case_id != "2132-1521-4impact", ]

# Save final clinical table
output_dir_clin_tsv <- "/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/"

clinical_outfile <- file.path(output_dir_clin_tsv, "clinical_annotations.tsv")
fwrite(combined_clin, clinical_outfile, sep = "\t")

# Save in RDS format
saveRDS(combined_clin, file.path(output_dir_clin_tsv, "clinical_annotations.rds"))

cat("\n-----Clinical annotation table and RDS saved in:", clinical_outfile, "\n")


cat("\n\n##### PIPELINE COMPLETED #####\n\n")
