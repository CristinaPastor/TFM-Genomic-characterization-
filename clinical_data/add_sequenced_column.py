import os
import pandas as pd

# Paths 
tabla_maestra_path = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_ped.csv"
base_path = "/storage/scratch01/groups/bu/impact_vuscan/familias/"

# Load table
df = pd.read_csv(tabla_maestra_path)

# Add new columns
germline_flags = []
tumor_flags = []

for _, row in df.iterrows():
    case_id = row["case_id"]
    family_id = str(row["family_id"]).zfill(4)  # confirm family_id is zero-padded to 4 digits
    family_path = os.path.join(base_path, family_id)

    # File paths
    germline_file = os.path.join(family_path, f"{case_id}-01.g.vcf.gz")
    tumor_file_02 = os.path.join(family_path, f"{case_id}-02.vcf.gz")
    tumor_file_04 = os.path.join(family_path, f"{case_id}-04.vcf.gz")
    tumor_file_06 = os.path.join(family_path, f"{case_id}-06.vcf.gz")

    # Flags
    germline_flags.append("yes" if os.path.exists(germline_file) else "no")
    tumor_flags.append("yes" if os.path.exists(tumor_file_02) or os.path.exists(tumor_file_04) or os.path.exists(tumor_file_06) else "no")

# Insert after "affected"
df.insert(df.columns.get_loc("affected") + 1, "sequenced_germline", germline_flags)
df.insert(df.columns.get_loc("affected") + 2, "sequenced_tumor", tumor_flags)

# Save result
output_path = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_seq.csv"
df.to_csv(output_path, index=False)

print(f"File saved to: {output_path}")
print("##### PIPELINE COMPLETED #####")