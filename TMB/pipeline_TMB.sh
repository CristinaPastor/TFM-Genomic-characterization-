#!/bin/bash

# Pipeline for TMB calculation using maftools, filtering non-synonymous variants, and adding it to cleaned_casos.csv grouped by individual.
# All cases that require manual review (according to the criteria) are logged in manual_cases_to_check_maftools.txt.
# Usage: bash pipeline_TMB.sh
# Cristina Pastor, October 2025

set -euo pipefail

# define directories and files
ROOT_DIR="/storage/scratch01/groups/bu/impact_vuscan/vcf2maf/output/"
OUT_DIR_MAFTOOLS="/storage/scratch01/groups/bu/impact_vuscan/TMB/maftools/output/"
CLINIC_CSV="/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos.csv"
CAPTURED_SIZE_FILE="/storage/scratch01/users/cpastor/0350_test/captured_regions/HyperExomeV2_capture_targets_size.txt"

TMB_DIR="/storage/scratch01/groups/bu/impact_vuscan/TMB/"
LOG_FILE="$TMB_DIR/pipeline_TMB_maftools.log"
ERR_FILE="$TMB_DIR/pipeline_TMB_maftools.err"

# stdout to log, stderr to error (also to terminal)
exec > >(tee -a "$LOG_FILE")
exec 2> >(tee -a "$ERR_FILE" >&2)

mkdir -p "$OUT_DIR_MAFTOOLS"

# Scripts
FILTER_SCRIPT="/storage/scratch01/groups/bu/impact_vuscan/TMB/scripts/filter_maf_nonsyn.py"
TMB_SCRIPT="/storage/scratch01/groups/bu/impact_vuscan/TMB/scripts/calculate_TMB_maftools.R"

# Prepare working CSV: create a copy to modify
WORKING_CSV="${CLINIC_CSV%/*}/cleaned_casos_with_TMB.csv"
if [[ ! -f "$WORKING_CSV" ]]; then
    cp "$CLINIC_CSV" "$WORKING_CSV"
    echo "[INFO] Created working cleaned_casos copy: $WORKING_CSV"
fi

# File to log manual cases
MANUAL_LOG="$TMB_DIR/manual_cases_to_check_maftools.txt"
echo "" > "$MANUAL_LOG"

AUTO_ADDED=0
MANUAL_COUNT=0

declare -A CASE_FILES_MAFTOOLS=()

# Iterate each family
for fam in "$ROOT_DIR"/*; do
    [[ -d "$fam" ]] || continue
    CASE_FILES_MAFTOOLS=()
    FAMID=$(basename "$fam")
    echo "== Processing family: $FAMID =="

    # Search MAF files for the family
    mapfile -t MAF_LIST < <(find "$fam" -type f -name "*-${FAMID}-4impact-0?.vep.maf")

    if [[ ${#MAF_LIST[@]} -eq 0 ]]; then
        echo "  [WARNING] No MAF found for family $FAMID"
        continue
    fi
    echo "  -> Found ${#MAF_LIST[@]} MAF(s)"

    mkdir -p "$OUT_DIR_MAFTOOLS/$FAMID"

    # Process each MAF file found
    for MAF in "${MAF_LIST[@]}"; do
        echo "    -> Processing MAF: $MAF"
        BASENAME=$(basename "$MAF" .vep.maf)

        FAM_OUT="$OUT_DIR_MAFTOOLS/$FAMID"
        mkdir -p "$FAM_OUT"

        # Filter non-synonymous variants
        FILTERED_MAF="$FAM_OUT/${BASENAME}_filtered.maf"
        echo "      -> Filtering MAF for maftools"
        python3 "$FILTER_SCRIPT" "$MAF" "$FILTERED_MAF"

        # Calculate TMB using maftools
        if [[ ! -s "$CAPTURED_SIZE_FILE" ]]; then
            echo "[ERROR] Captured size file missing or empty: $CAPTURED_SIZE_FILE"
            exit 1
        fi
        CAPTURED_SIZE=$(cat "$CAPTURED_SIZE_FILE")

        TMB_OUT="$FAM_OUT/${BASENAME}_TMB_maftools.txt"
        echo "      -> Calculating TMB (maftools)"
        Rscript "$TMB_SCRIPT" "$FILTERED_MAF" "$TMB_OUT" "$CAPTURED_SIZE"
        echo "         TMB saved in: $TMB_OUT"

        # Save mapping case_id_base -> files
        case_id_full=$(basename "$TMB_OUT" | cut -d"_" -f1)   # eg: 0350-0350-4impact-02
        case_id_base=$(echo "$case_id_full" | awk -F'-' '{print $1"-"$2"-"$3}') # eg: 0350-0350-4impact
        CASE_FILES_MAFTOOLS["$case_id_base"]+="$TMB_OUT "
    done

    # Add TMB to CSV grouped by individual
    ALL_CIDS=($(printf "%s\n" "${!CASE_FILES_MAFTOOLS[@]}" | sort -u))

    for cid in "${ALL_CIDS[@]}"; do
        FILES_MAFTOOLS=(${CASE_FILES_MAFTOOLS[$cid]:-})

        # Only one MAF:
        if [[ ${#FILES_MAFTOOLS[@]} -eq 1 ]]; then
            echo "  -> Adding TMB (maftools) for $cid"
            TMB_FILE="${FILES_MAFTOOLS[0]}"

            read AUTO_ADD_PY MANUAL_PY CASE_STATUS < <(python3 - <<PYCODE
import pandas as pd
import re, sys

csv_path = "$WORKING_CSV"
df = pd.read_csv(csv_path)

# Read TMB value
tmb_file = "$TMB_FILE"
with open(tmb_file) as f:
    tmb_val = float(f.read().strip())

sample = "$cid"
base_id = re.sub(r"-\d+$", "", sample)

mask = df["case_id"] == base_id
auto_added = 0
manual_count = 0
status = "auto"

if mask.sum() == 0:
    manual_count += 1
    status = "manual:no_caseid"
else:
    row = df.loc[mask].iloc[0]
    type2_val = row["type_2"] if "type_2" in row else None

    # Case 1: single clinical tumor (empty type_2)
    if pd.isna(type2_val) or str(type2_val).strip() == "":
        col_name = "TMB_1_maftools"
        if col_name not in df.columns:
            # insert after TMB_1 if exists, else at end
            insert_pos = list(df.columns).index("TMB_1") + 1 if "TMB_1" in df.columns else len(df.columns)
            df.insert(insert_pos, col_name, pd.NA)
        df.loc[mask, col_name] = tmb_val
        auto_added += 1

    # Case 2: multiple tumors but only one MAF
    else:
        tmb_cols = sorted([c for c in df.columns if re.match(r"TMB_\d+$", c)])
        filled_cols = [c for c in tmb_cols if not pd.isna(row[c])]

        if len(filled_cols) == 1:
            # There is only one filled TMB_X: assign TMB to that same number
            base_col = filled_cols[0]
            col_name = f"{base_col}_maftools"
            if col_name not in df.columns:
                insert_pos = list(df.columns).index(base_col) + 1
                df.insert(insert_pos, col_name, pd.NA)
            df.loc[mask, col_name] = tmb_val
            auto_added += 1
        else:
            # If there is no filled TMB_X or more than one filled TMB_X: manual
            manual_count += 1
            status = "manual:multi_tumor"

df.to_csv(csv_path, index=False)
print(f"{auto_added} {manual_count} {status}")

PYCODE
)

            AUTO_ADDED=$((AUTO_ADDED + AUTO_ADD_PY))
            MANUAL_COUNT=$((MANUAL_COUNT + MANUAL_PY))

            if [[ "$CASE_STATUS" == manual:* ]]; then
                echo "$cid [maftools] - $CASE_STATUS" >> "$MANUAL_LOG"
            fi

        # Multiple MAFs: manual
        elif [[ ${#FILES_MAFTOOLS[@]} -gt 1 ]]; then
            echo "  [INFO] Multiple MAFs for $cid -> manual update"
            for f in "${FILES_MAFTOOLS[@]}"; do
                echo "$cid - $(basename "$f") [maftools] - manual:multiple_mafs" >> "$MANUAL_LOG"
            done
            MANUAL_COUNT=$((MANUAL_COUNT + ${#FILES_MAFTOOLS[@]}))
        fi
    done

    echo "------- Family $FAMID finished -------"
done

echo "------------------- Pipeline completed --------------------"
echo "[INFO] Cases requiring manual TMB insertion saved in $MANUAL_LOG"
echo "[INFO] Total cases auto-added: $AUTO_ADDED"
echo "[INFO] Total cases requiring manual check: $MANUAL_COUNT"
