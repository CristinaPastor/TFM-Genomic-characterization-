#!/usr/bin/env python3
"""
Script to manually add TMB values to cleaned_casos.csv
Automatically creates the TMB_unknown columns in the correct order.
Usage: python manual_add_TMB.py
Cristina Pastor, October 2025
"""

import pandas as pd
import os

# File paths
CSV_FILE = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_TMB.csv"
OUTPUT_FILE = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_TMBman.csv"
MUTSCAPE_BASE = "/storage/scratch01/groups/bu/impact_vuscan/TMB/mutscape/output/"
MAFTOOLS_BASE = "/storage/scratch01/groups/bu/impact_vuscan/TMB/maftools/output/"

# Load CSV
df = pd.read_csv(CSV_FILE)

# Functions to get TMB values from external files
def get_from_mutscape(sample_path):
    """
    Read TMB_nonsynonmous from a Mutscape TMB_analysis.tsv file.

    sample_path format: "path/to/TMB_analysis.tsv"
    Example: "1837/1837-1837-4impact-02.vep/TMB_analysis.tsv"
    """
    try:
        # Extract the tumor name from the path (last folder before TMB_analysis.tsv)
        tumor_folder = os.path.basename(os.path.dirname(sample_path))
        sample_name = os.path.basename(os.path.dirname(sample_path)).replace(".vep", "").strip()
        
        if not os.path.exists(sample_path):
            print(f"[WARN] Not found: {sample_path}")
            return None
        
        # read TSV
        tsv = pd.read_csv(sample_path, sep="\t", comment="#")
        tsv.columns = [c.strip() for c in tsv.columns]
        if "sample" not in tsv.columns:
            print(f"[WARN] Columns found: {tsv.columns} in {sample_path}")
        tsv["sample"] = tsv["sample"].astype(str).str.strip()
        
        # Find the row matching the tumor_folder in the 'sample' column
        row = tsv.loc[tsv["sample"].str.startswith(sample_name), "TMB_nonsynonmous"]
        if row.empty:
            print(f"[WARN] No TMB_nonsynonmous found for {tumor_folder} in {sample_path}")
            return None
        
        return float(row.iloc[0])
    
    except Exception as e:
        print(f"[WARN] Error reading Mutscape for {sample_path}: {e}")
        return None


def get_from_maftools(sample_path):
    """
    Read TMB value from a Maftools _TMB_maftools.txt file.

    sample_path format: "path/to/sample_TMB_maftools.txt"
    Example: "1837/1837-1837-4impact-02_TMB_maftools.txt"
    """
    try:
        # sample_path example: "1837/1837-1837-4impact-02"
        # Build full path to the TMB_maftools.txt file
        base_dir = "/storage/scratch01/groups/bu/impact_vuscan/TMB/maftools/output"
        fam_id = os.path.basename(os.path.dirname(sample_path))
        sample_name = os.path.basename(sample_path)

        maftools_path = os.path.join(base_dir, fam_id, f"{sample_name}_TMB_maftools.txt")

        if not os.path.exists(maftools_path):
            print(f"[WARN] Not found: {maftools_path}")
            return None

        with open(maftools_path) as f:
            value = f.read().strip()
            if not value:
                print(f"[WARN] Empty TMB file: {maftools_path}")
                return None
            return float(value)

    except Exception as e:
        print(f"[WARN] Error reading maftools for {sample_path}: {e}")
        return None

# Create base columns TMB_unknown_1 and TMB_unknown_2 if not present
for base_col in ["TMB_unknown_1", "TMB_unknown_2"]:
    if base_col not in df.columns:
        df[base_col] = pd.NA
        print(f"[INFO] Column created: {base_col} (at the end of the CSV)")

# Function to insert new column right after another
def insert_column_after(df, new_col, after_col):
    """Insert a new empty column immediately after another."""
    if new_col in df.columns:
        return df
    cols = list(df.columns)
    if after_col in cols:
        idx = cols.index(after_col) + 1
        cols.insert(idx, new_col)
        df[new_col] = pd.NA
        df = df[cols]
        print(f"[INFO] Column created: {new_col} (after {after_col})")
    else:
        df[new_col] = pd.NA
        print(f"[INFO] Column created: {new_col} (at the end, {after_col} not found)")
    return df

def insert_column_at(df, new_col, idx):
    """Insert a new empty column at a specific index position."""
    if new_col in df.columns:
        return df
    cols = list(df.columns)
    cols.insert(idx, new_col)
    df[new_col] = pd.NA
    df = df[cols]
    print(f"[INFO] Column created: {new_col} (at position {idx})")
    return df


def ensure_column_position(df, base_col, tool=None):
    """
    Create TMB_unknown_?, TMB_unknown_?_mutscape and TMB_unknown_?_maftools columns
    at the end of the CSV in this order:
        base --> mutscape --> maftools
    It does NOT touch TMB_1, TMB_2, etc.
    """
    # Only apply to TMB_unknown columns
    if not base_col.startswith("TMB_unknown"):
        return df

    cols = list(df.columns)

    # Ensure base column exists
    if base_col not in cols:
        df = insert_column_at(df, base_col, len(cols))
        cols = list(df.columns)

    # No tool specified --> stop here
    if tool is None:
        return df

    # Determine new column name and position
    col_name = f"{base_col}_{tool}"
    if col_name in cols:
        return df  # already exists

    base_idx = cols.index(base_col)
    insert_idx = len(cols)  # default: end of file

    if tool == "mutscape":
        insert_idx = base_idx + 1
    elif tool == "maftools":
        mutscape_col = f"{base_col}_mutscape"
        insert_idx = cols.index(mutscape_col) + 1 if mutscape_col in cols else base_idx + 1

    df = insert_column_at(df, col_name, insert_idx)
    return df



# External data sources
sources = {
    "mutscape": get_from_mutscape,
    "maftools": get_from_maftools,
}

# Cases to add manually
# Format: ("case_id", "TMB_column", "source", "FAM_ID/sample_id_in_source")
manual_tmbs = [
    # Mutscape
    ("1837-1837-4impact", "TMB_2_mutscape", "mutscape", "9.355873144734987"), #-02 liver
    ("1837-1837-4impact", "TMB_1_mutscape", "mutscape", "4.307406942873039"), #-04 breast
    ("1884-1884-4impact", "TMB_1_mutscape", "mutscape", "4.955833794488335"), #-04 sarcoma
    ("1884-1884-4impact", "TMB_2_mutscape", "mutscape", "8.753762496806498"), #-02 mama 
    ("2079-2077-4impact", "TMB_2_mutscape", "mutscape", "10.374829625844738"), # mama
    ("2560-2560-4impact", "TMB_1_mutscape", "mutscape", "0.9494821755795407"), #-04 stomach
    ("2560-2560-4impact", "TMB_2_mutscape", "mutscape", "17.438050688082786"), #-02 liver
    ("2748-2748-4impact", "TMB_1_mutscape", "mutscape", "18.2254261507585"), #-04 renal último
    ("2748-2748-4impact", "TMB_2_mutscape", "mutscape", "9.957983792663477"), #-02 pheochromocytoma
    ("3018-3018-4impact", "TMB_1_mutscape", "mutscape", "3.682138193101146"), # sarcoma
    ("3700-3700-4impact", "TMB_1_mutscape", "mutscape", "1.55159282350803"), #-04 bladder
    ("3700-3700-4impact", "TMB_2_mutscape", "mutscape", "1.5747509253514334"), #-02 renal

    # Maftools 
    ("1837-1837-4impact", "TMB_2_maftools", "maftools", "1837/1837-1837-4impact-02"), #-02 liver 9.355873 
    ("1837-1837-4impact", "TMB_1_maftools", "maftools", "1837/1837-1837-4impact-04"), #-04 breast 4.307407
    ("1884-1884-4impact", "TMB_1_maftools", "maftools", "1884/1884-1884-4impact-04"), #-04 sarcoma 4.955834
    ("1884-1884-4impact", "TMB_2_maftools", "maftools", "1884/1884-1884-4impact-02"), #-02 mama 8.753762
    ("2079-2077-4impact", "TMB_2_maftools", "maftools", "2077/2079-2077-4impact-02"), # 10.37483
    ("2560-2560-4impact", "TMB_1_maftools", "maftools", "2560/2560-2560-4impact-04"), #-04 stomach 0.9494822
    ("2560-2560-4impact", "TMB_2_maftools", "maftools", "2560/2560-2560-4impact-02"), #-02 liver 17.43805
    ("2748-2748-4impact", "TMB_1_maftools", "maftools", "2748/2748-2748-4impact-04"), #-04 renal 18.22543
    ("2748-2748-4impact", "TMB_2_maftools", "maftools", "2748/2748-2748-4impact-02"), #-02 pheochromocytoma 9.957984
    ("3018-3018-4impact", "TMB_1_maftools", "maftools", "3018/3018-3018-4impact-02"), # 3.682138
    ("3700-3700-4impact", "TMB_1_maftools", "maftools", "3700/3700-3700-4impact-04"), #-04 1.551593
    ("3700-3700-4impact", "TMB_2_maftools", "maftools", "3700/3700-3700-4impact-02") #-02 1.574751
]

# Process each case
for case_id, column, source, sample_path in manual_tmbs:
    if source not in sources:
        print(f"[WARN] Unknown source: {source}")
        continue

    getter = sources[source]
    # Build full path if it's from mutscape or maftools
    if source == "mutscape":
        full_path = os.path.join(MUTSCAPE_BASE, sample_path)
    elif source == "maftools":
        full_path = os.path.join(MAFTOOLS_BASE, sample_path)
    else:
        full_path = sample_path

    # Check if sample_path is a number
    try:
        value = float(sample_path)
        print(f"[INFO] Using manual numeric TMB value {value} for {case_id} ({source})")
    except ValueError:
        value = getter(full_path)

    if value is None:
        print(f"[WARN] No value found for {sample_path} in {source}.")
        continue

    parts = column.split("_")
    if parts[1] == "unknown":
        base_col = "_".join(parts[:3])  # TMB_unknown_1
    else:
        base_col = "_".join(parts[:2])  # TMB_1, TMB_2, etc.

    df = ensure_column_position(df, base_col, source)

    mask = df["case_id"] == case_id
    if not mask.any():
        print(f"[WARN] Case {case_id} not found in CSV.")
        continue

    df.loc[mask, column] = value

# Reorder columns
maftools_cols = sorted([c for c in df.columns if c.startswith("TMB_unknown_") and "maftools" in c])
mutscape_cols = sorted([c for c in df.columns if c.startswith("TMB_unknown_") and "mutscape" in c])
tmb_order = ["TMB_unknown_1", "TMB_unknown_2"] + maftools_cols + mutscape_cols

# Apply reordering
unknown_cols = [c for c in tmb_order if c in df.columns]
unknown_df = df[unknown_cols]
main_df = df.drop(columns=unknown_cols, errors="ignore")
df = pd.concat([main_df, unknown_df], axis=1)

# Save updated CSV file
df.to_csv(OUTPUT_FILE, index=False)
print(f"Updated file saved to: {OUTPUT_FILE}")

print("\n##### PIPELINE COMPLETED #####\n")