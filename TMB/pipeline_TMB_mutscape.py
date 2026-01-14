#!/usr/bin/env python3

"""
Pipeline to calculate TMB using MutScape and add it to the clinical CSV.
All cases requiring manual review are logged in manual_cases_to_check_mutscape.txt.
Usage: python3 pipeline_TMB_mutscape.py > pipeline_TMB_mutscape_log.txt 2> pipeline_TMB_mutscape_error.txt
Author: Cristina Pastor, October 2025
"""

import os
import subprocess
import pandas as pd
import shutil
import re

# Define paths
families_dir = "/storage/scratch01/groups/bu/impact_vuscan/vcf2maf/output/"  # Folder with subfolders per family_id
mutscape_script = "/storage/scratch01/users/cpastor/software/MutScape/mutscape/mafAnalysis.py"
tmb_bp_file = "/storage/scratch01/users/cpastor/0350_test/captured_regions/HyperExomeV2_capture_targets_size.txt"
output_base = "/storage/scratch01/groups/bu/impact_vuscan/TMB/mutscape/output"

clinical_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_TMB.csv"
manual_txt = "/storage/scratch01/groups/bu/impact_vuscan/TMB/manual_cases_to_check_mutscape.txt"


print("Starting TMB pipeline with MutScape...")

# Read TMB value in bp
with open(tmb_bp_file) as f:
    tmb_bp = f.read().strip()
print(f"TMB value in base pairs loaded: {tmb_bp}")

# Load clinical CSV
print(f"Loading clinical CSV: {clinical_csv}")
df_clinical = pd.read_csv(clinical_csv)

# Backup original CSV (the one with only maftools TMBs) before modifying
backup_csv = clinical_csv.replace(".csv", "_maftools.csv")
shutil.copy(clinical_csv, backup_csv)
print(f"Backup of clinical CSV created: {backup_csv}")

# Create output folder if it doesn't exist
os.makedirs(output_base, exist_ok=True)
print(f"Output base folder: {output_base}")

# Lists to store manual review cases and automatic updates
manual_cases = []
auto_count = 0
manual_count = 0
total_cases = 0
multi_tumor_manuals = 0

# Iterate over each family folder
for fam_id in os.listdir(families_dir):
    fam_path = os.path.join(families_dir, fam_id)
    if not os.path.isdir(fam_path):
        continue

    print(f"\nProcessing family {fam_id}...")
    fam_output = os.path.join(output_base, fam_id)
    os.makedirs(fam_output, exist_ok=True)

    # DataFrame to accumulate TSVs per family
    df_family_tsv = pd.DataFrame()

    # Find all MAF files for the family
    maf_files = [os.path.join(fam_path, f) for f in os.listdir(fam_path) if f.endswith(".vep.maf")]
    if not maf_files:
        print(f"  No MAF files found for family {fam_id}, skipping.")
        continue

    # Run MutScape for each MAF
    for maf in maf_files:
        sample_base = os.path.splitext(os.path.basename(maf))[0]  # eg: 1837-1837-4impact-02
        maf_output_dir = os.path.join(fam_output, sample_base)
        os.makedirs(maf_output_dir, exist_ok=True)

        print(f"  Running MutScape on {maf} in folder {maf_output_dir}...")
        cmd = [
            "python3", mutscape_script,
            "-f", maf,
            "-tmb", tmb_bp,
            "-o", maf_output_dir,
            "-p", maf_output_dir
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"    Successfully processed {maf}")
        except subprocess.CalledProcessError as e:
            print(f"    ERROR: MutScape failed for {maf}: {e}")
            continue

        # Read TSV after MutScape
        tsv_file = os.path.join(maf_output_dir, "TMB_analysis.tsv")
        if not os.path.exists(tsv_file):
            print(f"    WARNING: {tsv_file} not generated")
            continue

        df_tsv = pd.read_csv(tsv_file, sep="\t")

        # Clean 'sample' column
        df_tsv["sample"] = df_tsv["sample"].apply(lambda x: os.path.splitext(os.path.basename(str(x)))[0])
        for ext in [".annotated.somatic.renamed", ".vep", ".somatic", ".maf", ".vcf"]:
            df_tsv["sample"] = df_tsv["sample"].str.replace(ext, "", regex=False)

        df_tsv["case_id"] = df_tsv["sample"].apply(lambda x: "-".join(x.split("-")[:-1]))

        # accumulate in the family DataFrame
        df_family_tsv = pd.concat([df_family_tsv, df_tsv], ignore_index=True)

        print("\n[DEBUG] df_family_tsv for now:")
        print(df_family_tsv[["case_id", "sample", "TMB_nonsynonmous"]])


    # Insert TMB into clinical CSV per individual
    counts_per_individual = df_family_tsv.groupby("case_id")["sample"].nunique().to_dict()

    for case_id_check, maf_count in counts_per_individual.items():
        mask_ind = df_family_tsv["case_id"] == case_id_check
        samples_ind = df_family_tsv.loc[mask_ind, "sample"].tolist()

        # Case 0: Multiple tumors sequenced: manual
        if maf_count > 1:
            for sample in samples_ind:
                manual_cases.append(f"{fam_id}\t{case_id_check}\t{sample}\tMULTIPLE_TUMORS_SEQUENCED")
                manual_count += 1
            continue

        # Only one tumor sequenced
        sample = samples_ind[0]
        tmb_value = df_family_tsv.loc[df_family_tsv["sample"] == sample, "TMB_nonsynonmous"].values[0]
        mask_clinical = df_clinical["case_id"] == case_id_check

        if not mask_clinical.any():
            manual_cases.append(f"{fam_id}\t{case_id_check}\t{sample}\tNO_CASE_IN_CSV")
            manual_count += 1
            continue

        # Extract the row in clinical CSV
        row = df_clinical.loc[mask_clinical].iloc[0]

        # Identify filled type and TMB columns
        filled_tmb_cols = [
            c for c in df_clinical.columns
            if re.match(r"^TMB_\d+$", c)  # Only TMB_1, TMB_2, etc. (without suffixes)
            and pd.notna(row[c])
            and str(row[c]).strip() not in ["", "NA", "NaN"]
        ]


        filled_type_cols = [
            c for c in df_clinical.columns
            if c.startswith("type_")
            and pd.notna(row[c]) and str(row[c]).strip() not in ["", "NA", "NaN", "."]
        ]

        print(f"DEBUG {fam_id} {case_id_check}: types={filled_type_cols}, tmbs={filled_tmb_cols}, maf_count={maf_count}")

        # CASE 1: Already a filled TMB_n: put after
        if len(filled_tmb_cols) == 1:
            n = filled_tmb_cols[0].split("_")[1]
            base_col = f"TMB_{n}"
            new_col = f"{base_col}_mutscape"

            if new_col not in df_clinical.columns:
                insert_pos = df_clinical.columns.get_loc(base_col) + 1
                df_clinical.insert(insert_pos, new_col, pd.NA)

            df_clinical.loc[mask_clinical, new_col] = tmb_value
            auto_count += 1
            total_cases += 1
            continue

        # Case 2: No TMB_n but only one type_n filled: put after
        if not filled_tmb_cols and len(filled_type_cols) == 1:
            n = filled_type_cols[0].split("_")[1]
            base_col = f"TMB_{n}"
            new_col = f"{base_col}_mutscape"

            if new_col not in df_clinical.columns:
                if base_col in df_clinical.columns:
                    insert_pos = df_clinical.columns.get_loc(base_col) + 1
                else:
                    insert_pos = len(df_clinical.columns)
                df_clinical.insert(insert_pos, new_col, pd.NA)

            df_clinical.loc[mask_clinical, new_col] = tmb_value
            auto_count += 1
            total_cases += 1
            continue

        # Case 3: Multiple type_n and no TMB_n filled: manual
        if not filled_tmb_cols and len(filled_type_cols) > 1:
            manual_cases.append(f"{fam_id}\t{case_id_check}\t{sample}\tMULTIPLE_TYPES_NO_TMB")
            manual_count += 1
            total_cases += 1
            continue

        # Other cases: manual
        manual_cases.append(f"{fam_id}\t{case_id_check}\t{sample}\tAMBIGUOUS_CASE")
        manual_count += 1
        total_cases += 1



    # Flag multi-tumor cases within the family
    counts_per_individual = df_family_tsv.groupby("case_id")["sample"].nunique().to_dict()

    for case_id_check, count in counts_per_individual.items():
        if count > 1:
            print(f"  Family {fam_id}: case_id {case_id_check} has multiple tumors ({count}) → manual")
            samples = df_family_tsv.loc[df_family_tsv["case_id"] == case_id_check, "sample"].tolist()

            for sample in samples:
                entry = f"{fam_id}\t{case_id_check}\t{sample}\tMULTIPLE_TUMORS_SEQUENCED"
                if entry not in manual_cases:  #avoid duplicates
                    manual_cases.append(entry)
                    manual_count += 1
            multi_tumor_manuals += 1


# Save updated clinical CSV
df_clinical.to_csv(clinical_csv, index=False)
print(f"\nClinical CSV updated with MutScape TMB values: {clinical_csv}")

# Save manual cases
with open(manual_txt, "w") as f:
    f.write("family_id\tcase_id\tsample\treason\n")
    f.write("\n".join(manual_cases))
print(f"Manual review cases saved to: {manual_txt}")

# Summary
print("\n--- PIPELINE SUMMARY ---")
print(f"Total cases processed: {total_cases}")
print(f"Cases added automatically: {auto_count}")
unique_manuals = set(manual_cases)  # evita duplicados globales
print(f"Cases requiring manual review: {len(unique_manuals)} (including multi-tumor cases)")
print(f"   Manual detected during processing: {len(unique_manuals) - multi_tumor_manuals}")
print(f"   Manual due to multi-tumor: {multi_tumor_manuals}")
print("-------------------------")
print("\nPipeline completed")
