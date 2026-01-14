"""
Clean MAF file (select non-synonymous) to calculate Tumor Mutation Burden (TMB).

Cristina Pastor, September 2025
"""

import sys
import pandas as pd 

# Check input arguments
if len(sys.argv) != 3:
    print("Usage: python filter_maf_nonsyn.py <input_maf_file.maf> <output_filtered_maf>")
    sys.exit(1)

# The input file (.maf obtained after ejecuting vcf2maf) is the second argument
input_file = sys.argv[1]
# Define the output_file
output_file = sys.argv[2]

# Keep header lines
header_lines = []
with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            header_lines.append(line)
        else:
            break

# Read the .maf file
if input_file.endswith(".maf"):
    maf = pd.read_csv(input_file, sep = "\t", comment = "#")
else:
    print("Unsupported file format. Use .maf")
    sys.exit(1)

# Filter non-synonymous mutations in 'Variant_Classification' column:
non_synonymous_maf = maf[maf['Variant_Classification'].isin([
    "Splice_Site",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonstop_Mutation",
    "In_Frame_Ins",
    "In_Frame_Del",
    "Missense_Mutation",
    "Translation_Start_Site"
])]

# Save filtered maf with header lines and filtered non-synonymous mutations
with open(output_file, "w") as f:
    f.writelines(header_lines) 
    non_synonymous_maf.to_csv(f, sep="\t", index=False)


print(f"Filtered maf saved to: {output_file}")