#!/usr/bin/env python3
"""
Script to add tumor_id in cleaned_casos.csv
Usage: python add_tumor_id.py
Cristina Pastor, October 2025
"""

import os
import pandas as pd
import re

# Rutas
families_dir = "/storage/scratch01/groups/bu/impact_vuscan/vcf2maf/output/"
clinical_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_TMBman.csv"
out_clinical_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_tumor_id.csv"
manual_txt = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/manual_tumor_id.txt"

# Load clinical csv
print(f"Loading clinical CSV: {clinical_csv}")
df_clinical = pd.read_csv(clinical_csv)

# Create dictionary case_id: row in the CSV
clinical_rows = {}
for idx, row in df_clinical.iterrows():
    cid = row["case_id"]  # case_id of each row
    clinical_rows[cid] = idx  # idx is the row number (position in the CSV)


manual_cases = []

# Iterate over each family
for fam_id in os.listdir(families_dir):
    fam_path = os.path.join(families_dir, fam_id)
    if not os.path.isdir(fam_path):
        continue

    # Only .vep.maf files
    maf_files = [os.path.join(fam_path, f) for f in os.listdir(fam_path) if f.endswith(".vep.maf")]
    if not maf_files:
        print(f"[INFO] No MAF files found for family {fam_id}...")
        continue

    maf_info = []

    # Extract case_id and tumor_code from each maf
    for maf in maf_files:
        maf_filename = os.path.basename(maf)
        m = re.match(r'^(\d+-\d+-.*)-(\d+)\.vep\.maf$', maf_filename)

        if not m:
            print(f"[INFO] Filename did not match regex: {maf_filename}")
            continue

        case_id = m.group(1)
        tumor_code = m.group(2)

        maf_info.append({
            "case_id": case_id,
            "tumor_code": tumor_code
        })

    if not maf_info:
        continue

    # Aggregate tumor_codes by case_id
    case_id_tumor_code = {}
    for m in maf_info:
        case_id = m["case_id"]
        tumor_code = m["tumor_code"]

        if case_id not in case_id_tumor_code:
            case_id_tumor_code[case_id] = []

        case_id_tumor_code[case_id].append(tumor_code)

    print(f"[INFO] Family {fam_id} → {case_id_tumor_code}")

    # detect relevant columns
    type_cols = {}        # n --> type_n
    tmb_cols = {}         # n --> TMB_n
    tmb_unknown_cols = {} # n --> TMB_unknown_n

    for col in df_clinical.columns:
        m_type = re.match(r"type_(\d+)$", col)
        if m_type:
            type_cols[int(m_type.group(1))] = col

        m_tmb = re.match(r"TMB_(\d+)$", col)
        if m_tmb:
            tmb_cols[int(m_tmb.group(1))] = col

        m_tmb_unknown = re.match(r"TMB_unknown_(\d+)$", col)
        if m_tmb_unknown:
            tmb_unknown_cols[int(m_tmb_unknown.group(1))] = col

    # order columns by n
    type_cols = dict(sorted(type_cols.items()))
    tmb_cols = dict(sorted(tmb_cols.items()))
    tmb_unknown_cols = dict(sorted(tmb_unknown_cols.items()))

    # inicialize manual cases for this family

    ### Case 1: individuals with multiple tumors: manual
    for case_id, tumors in case_id_tumor_code.items():
        if len(tumors) > 1:
            manual_cases.append(f"{case_id}\tMULTIPLE_TUMORS")

    ### function for filled cells
    def is_filled(v):
        if pd.isna(v):
            return False
        if isinstance(v, str) and v.strip() == "":
            return False
        return True

    ### Case 2: individuals with a single tumor: assign tumor_id_n
    for case_id, tumors in case_id_tumor_code.items():

        # Skip manual cases
        if case_id in [m.split("\t")[0] for m in manual_cases]:
            continue

        # Only individuals with 1 tumor
        if len(tumors) != 1:
            continue

        tumor_code = tumors[0]
        tumor_id_value = f"{case_id}-{tumor_code}"

        # Obtain row in CSV
        row_idx = clinical_rows.get(case_id)
        if row_idx is None:
            manual_cases.append(f"{case_id}\tNO_CLINICAL_ROW")
            continue

        # case 2A: only a type_n filled: use that n
        filled_type = [n for n, col in type_cols.items() if is_filled(df_clinical.at[row_idx, col])]
        if len(filled_type) == 1:
            n = filled_type[0]

            tumor_col = f"tumor_id_{n}"
            if tumor_col not in df_clinical.columns:
                insert_pos = df_clinical.columns.get_loc(type_cols[n])
                df_clinical.insert(insert_pos, tumor_col, "")

            df_clinical.at[row_idx, tumor_col] = tumor_id_value
            continue

        # case 2B: multiple type_n filled: use TMB_n
        if len(filled_type) > 1:
            filled_tmb = [n for n, col in tmb_cols.items() if is_filled(df_clinical.at[row_idx, col])]

            if len(filled_tmb) == 1:
                n = filled_tmb[0]

                tumor_col = f"tumor_id_{n}"
                if tumor_col not in df_clinical.columns:
                    insert_pos = df_clinical.columns.get_loc(type_cols[n])
                    df_clinical.insert(insert_pos, tumor_col, "")

                df_clinical.at[row_idx, tumor_col] = tumor_id_value
                continue

            if len(filled_tmb) > 1:
                manual_cases.append(f"{case_id}\tMULTIPLE_TMB")
                continue

        # Case 2C: TMB_unknown_n filled
        filled_tmb_unknown = [n for n, col in tmb_unknown_cols.items() if is_filled(df_clinical.at[row_idx, col])]
        if len(filled_tmb_unknown) == 1:
            n = filled_tmb_unknown[0]

            tumor_col = f"tumor_unknown_id_{n}"
            if tumor_col not in df_clinical.columns:
                insert_pos = df_clinical.columns.get_loc(tmb_unknown_cols[n])
                df_clinical.insert(insert_pos, tumor_col, "")

            df_clinical.at[row_idx, tumor_col] = tumor_id_value
            continue

        if len(filled_tmb_unknown) > 1:
            manual_cases.append(f"{case_id}\tMULTIPLE_TMB_UNKNOWN")
            continue

        # case 2D: ambiguous case: manual
        manual_cases.append(f"{case_id}\tAMBIGUOUS_SINGLE")

# Save updated clinical CSV
print(f"Saving {os.path.basename(out_clinical_csv)} in {out_clinical_csv}...")
df_clinical.to_csv(out_clinical_csv, index=False)

# save manual cases to check
print(f"Saving manual cases to check in {manual_txt}...")
with open(manual_txt, "w") as f:
    f.write("\n".join(manual_cases))

print("##### PIPELINE FINISHED #####")

