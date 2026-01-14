#!/usr/bin/env python

# Annotate and merge CNVs (from tumors) using gencode file
# Cristina Pastor, November 2025

import sys
import os
import argparse
import gzip
import re
import glob

parser = argparse.ArgumentParser(description="Annotate CNV file.")
parser.add_argument("--gencode", dest="gencode", required=False, default="/storage/scratch01/groups/bu/dcarrero_common/databases/gencode.v44.annotation_reformat.bed", help="Input gencode file with gene names and coordinates.")
parser.add_argument("--outFile", dest="outFile", required=False, default="/storage/scratch01/groups/bu/impact_vuscan/cohort_characterization/oncoplots/all_cnvs_tumors.csv", help="Output annotated CNV bed file with all samples.")


base_dir = "/storage/scratch01/groups/bu/impact_vuscan/familias/"

args = parser.parse_args()


# Load gencode bed
genes = []
with open(args.gencode) as bed:
    for line in bed:
        if line.strip() and not line.startswith("#"):
            chrom, start, end, strand, gene = line.strip().split()[:5]
            genes.append((chrom.replace("chr", ""), int(start), int(end), gene))

# Look for files
# Remove germinal CNVs files (01)
cnv_basic_files = [
    f for f in glob.glob(os.path.join(base_dir, "????", "*.CNVCalls.Tumor_vs_Normal.tsv"))
    if not os.path.basename(f).split("-")[-1].startswith("01")
]

cnv_ann_files = [
    f for f in glob.glob(os.path.join(base_dir, "????", "*.CNVCalls*ann.tag.tsv"))
    if not os.path.basename(f).split("-")[-1].startswith("01")
]

cnv_ann_files += [
    f for f in glob.glob(os.path.join(base_dir, "????", "*.CNVs.PON.annotated.tsv"))
    if not os.path.basename(f).split("-")[-1].startswith("01")
]

cnv_ann_files += [
    f for f in glob.glob(os.path.join(base_dir, "????", "*.CNVs.annotated_parsed.tsv"))
    if not os.path.basename(f).split("-")[-1].startswith("01")
]

# Files that will be processed:
processed_files = set(cnv_basic_files + cnv_ann_files)

with open(args.outFile, "w") as out:
    out.write("Gene\tSample_name\tCN\n")

    # Process files .CNVCalls.Tumor_vs_Normal.tsv (Navarra y Galicia)
    for cnv_file in cnv_basic_files:
        # Skip .ann.tag.tsv or .annotated_parsed.tsv files
        if cnv_file.endswith(".ann.tag.tsv") or cnv_file.endswith(".annotated_parsed.tsv"):
            continue

        sample_name = os.path.basename(cnv_file).split(".CNVCalls.Tumor_vs_Normal.tsv")[0]
        print(f"[BÁSICO] Procesando {sample_name}...")

        with open(cnv_file) as cnv:
            header = cnv.readline()
            for line in cnv:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                chrom, start, end, cn_type = cols[1], int(cols[2]), int(cols[3]), cols[5]
                cn_type = cn_type.replace("<", "").replace(">", "")
                cn_type = cn_type.replace("<", "").replace(">", "")  # quita los símbolos

                # Divide with coma if there are various types
                types = cn_type.split(",")  # eg: ["DEL","DUP"]

                for t in types:
                    t_clean = t.strip()
                    if t_clean == "DEL":
                        final_cn = "Del"
                    elif t_clean == "DUP":
                        final_cn = "Amp"
                    else:
                        print(t_clean) # skip LOH
                        continue

                for gchrom, gstart, gend, gene in genes:
                    if gchrom == chrom and not (end < gstart or start > gend):
                        out.write(f"{gene}\t{sample_name}\t{cn_type}\n")

    # Process files .CNVCalls.Tumor_vs_Normal.ann.tag.tsv (CNAG)
    for ann_file in cnv_ann_files:
        sample_name = os.path.basename(ann_file).split(".")[0]
        print(f"[ANNOTADED] Processing {sample_name}...")

        with open(ann_file) as ann:
            header = ann.readline().strip().split("\t")
            try:
                chrom_idx = header.index("SV_chrom")
                start_idx = header.index("SV_start")
                end_idx = header.index("SV_end")
                type_idx = header.index("SV_type")
            except ValueError:
                print(f"{ann_file} does not have expected columns, it is omitted")
                continue

            for line in ann:
                cols = line.strip().split("\t")
                if len(cols) <= max(chrom_idx, start_idx, end_idx, type_idx):
                    continue
                chrom, start, end, sv_type = cols[chrom_idx], int(cols[start_idx]), int(cols[end_idx]), cols[type_idx]
                if sv_type == "DEL":
                    cn_type = "Del"
                elif sv_type == "DUP":
                    cn_type = "Amp"
                else:
                    print(sv_type)

                for gchrom, gstart, gend, gene in genes:
                    if gchrom == chrom and not (end < gstart or start > gend):
                        out.write(f"{gene}\t{sample_name}\t{cn_type}\n")

# Check not processed CNV files
all_cnv_files = glob.glob(os.path.join(base_dir, "????", "*CNV*tsv"))
unprocessed = [f for f in all_cnv_files if f not in processed_files]

if unprocessed:
    print("\nFiles with 'CNV' in the name that were not processed:")
    for f in unprocessed:
        print("   ", f)
else:
    print("\nAll files with 'CNV' in the name were processed")