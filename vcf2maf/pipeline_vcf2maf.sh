#!/bin/bash

# Pipeline to convert all tumor VCFs to MAFs
# Usage: bash pipeline_vcf2maf.sh
# Cristina Pastor, September 2025

set -euo pipefail

# Define directories and files
ROOT_DIR="/storage/scratch01/groups/bu/impact_vuscan/familias"
OUT_DIR="/storage/scratch01/groups/bu/impact_vuscan/vcf2maf/output2/"

VCF2MAF="/storage/scratch01/users/cpastor/software/vcf2maf-1.6.22/vcf2maf.pl"
REF_FASTA="/storage/scratch01/users/cpastor/genomes/hg38/hg38.fa"
VEP_PATH="/home/cpastor/miniforge3/envs/vep/bin"

mkdir -p "$OUT_DIR"

# Iterate each family
for fam in "$ROOT_DIR"/*; do
    [[ -d "$fam" ]] || continue
    FAM_ID=$(basename "$fam") # remove all the path
    echo "--- Processing family: $FAM_ID ..."

    PRIOR_DIR="$fam/output/prioritization"
    # If PRIOR_DIR does not exist warn me 
    if [[ ! -d "$PRIOR_DIR" ]]; then
        echo "--[WARNING] No prioritization dir found for family $FAM_ID"
        continue
    fi

    # Search for all the cases in prioritization directory 
    # exclude files containing filtered and v2
    mapfile -t CASES < <(find "$PRIOR_DIR" -maxdepth 1 -type f \
        -name "*-${FAM_ID}-4impact-0?.annotated.somatic*.vcf" ! -name "*filtered.vcf" ! -name "*v2*" \ 
        -printf "%f\n" | sed -E "s/-${FAM_ID}-4impact.*//" | sort -u)

    # if no VCF files are found warn me
    if [[ ${#CASES[@]} -eq 0 ]]; then
        echo "--[WARNING] No cases found for family $FAM_ID"
        continue
    fi

    VCF_ARRAY=()

    for CASE in "${CASES[@]}"; do
        echo "--- Processing case: $CASE ..."

        # Search for all the case's VCFs (there can be several tumors)
        # exclude files containing filtered and v2
        mapfile -t NORMAL_VCFS < <(find "$PRIOR_DIR" -maxdepth 1 -type f \
            -name "${CASE}-${FAM_ID}-4impact-0?.annotated.somatic*.vcf" ! -name "*filtered.vcf" ! -name "*v2*" | sort)

        # If there are no VCFs skip the case
        if [[ ${#NORMAL_VCFS[@]} -eq 0 ]]; then
            echo "--[WARNING] No VCFs found for case $CASE in family $FAM_ID"
            continue
        fi

        # Add all the VCFs found in VCF_ARRAY array
        for vcf in "${NORMAL_VCFS[@]}"; do
            echo "------Found VCF: $vcf"
            VCF_ARRAY+=("$vcf")
        done

    done

    # Process each VCF file in VCF_ARRAY
    for VCF in "${VCF_ARRAY[@]}"; do
        echo "---Processing VCF: $VCF"

        # Create output directory for each family with its id
        FAM_OUT="$OUT_DIR/$FAM_ID"
        mkdir -p "$FAM_OUT"

        # Define basename
        BASENAME=$(basename "$VCF")

        if [[ "$BASENAME" =~ ^([0-9]{4}-[0-9]{4}-4impact-0[0-9]) ]]; then # basename has to have 4 numbers-4 numbers-4impact-two numbers (0350-0350-4impact-02)
            SAMPLE_BASE="${BASH_REMATCH[1]}"
        else
            echo "-----[WARNING] Cannot parse sample base from $BASENAME"
            continue
        fi

        # Define output MAF file name and tumor and normal ids
        MAF_OUT="$FAM_OUT/${SAMPLE_BASE}.vep.maf"
        TUMOR_ID="$SAMPLE_BASE"
        NORMAL_ID="${SAMPLE_BASE%-0?}-01" # remove '-0X' and add '-01'

        echo "-----[CHECK] Tumor ID: $TUMOR_ID"
        echo "-----[CHECK] Normal ID: $NORMAL_ID"
        echo "-----[CHECK] MAF file name: $MAF_OUT"

        # Run vcf2maf to convert VCF to MAF
        echo "-----Running vcf2maf..."
        perl "$VCF2MAF" \
            --input-vcf "$VCF" \
            --output-maf "$MAF_OUT" \
            --ref-fasta "$REF_FASTA" \
            --vep-path "$VEP_PATH" \
            --species "homo_sapiens" \
            --ncbi-build "GRCh38" \
            --tumor-id "$TUMOR_ID" \
            --normal-id "$NORMAL_ID" \
            --tmp-dir "$FAM_OUT" \
            --vep-forks=16

        echo "-----MAF generated: $MAF_OUT"
    done
done

echo "##### PIPELINE COMPLETED #####"
