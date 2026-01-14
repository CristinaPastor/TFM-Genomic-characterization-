import csv
import os

input_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_TMBman.csv"
output_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_ped.csv"

rows = []

with open(input_csv, newline="") as infile:
    reader = csv.DictReader(infile)
    fieldnames = reader.fieldnames
    
    # Insert new column right after "affected"
    if "affected" not in fieldnames:
        raise ValueError("Column 'affected' not found in CSV.")
    idx = fieldnames.index("affected") + 1
    new_fieldnames = fieldnames[:idx] + ["affected_ped"] + fieldnames[idx:]
    
    for row in reader:
        case_id = row["case_id"]
        family_id = row["family_id"]
        ped_file = f"/storage/scratch01/groups/bu/impact_vuscan/familias/{family_id}/{family_id}.ped"

        affected_status = "NA"
        if os.path.exists(ped_file):
            with open(ped_file) as f:
                lines = f.readlines()
            
            # If there is a header "aff", we ignore it
            if lines and lines[0].strip().split()[-1].lower() == "aff":
                data_lines = lines[1:]
            else:
                data_lines = lines
            
            for line in data_lines:
                parts = line.strip().split()
                if not parts:
                    continue
                if parts[0] == f"{case_id}-01":
                    last_col = parts[-1]
                    if last_col == "1":
                        affected_status = "no"
                    elif last_col == "2":
                        affected_status = "yes"
                    break
        
        row["affected_ped"] = affected_status
        rows.append(row)

# save new CSV with columns in correct order
with open(output_csv, "w", newline="") as outfile:
    writer = csv.DictWriter(outfile, fieldnames=new_fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow({col: row.get(col, "") for col in new_fieldnames})

print(f"File saved to: {output_csv}")
