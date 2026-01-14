"""
Clean the .csv file from IMPACT platform

Cristina Pastor, July 2025
"""
import sys
import pandas as pd 
import os
import csv
from collections import defaultdict

# Check input arguments
if len(sys.argv) < 2:
    print("Usage: python clean_impact_csv.py <input_impact_file.csv>")
    sys.exit(1)

# The input file (.csv converted from .xls downloaded from IMPACT platform) is the second argument
input_file = sys.argv[1]
# Define the output_file
output_file = "cleaned_" + os.path.basename(input_file)

# Read the .csv file
if input_file.endswith(".csv"):
    df = pd.read_csv(input_file)
else:
    print("Unsupported file format. Use .csv")
    sys.exit(1)


# 1. SELECT COLUMNS 
# Define columns of interest
columns_of_interest = [
    'IdCaso',
    'PerfilIndividuo',
    'Sexo',
    'CasoFamiliar',
    'Survival',
    'PruebasComplementarias', # Where HBOC is (PC)
    'DatosGeneticosAnalisis', # In some cases there is commercial panel info (DGA)
    'TratamientosEIntervencionesTabla', # Where treatment is (TEIT)
    'DiagnosisAndTumorInformation' # (DATI)
]

# Check if all columns exist
missing = []
for col in columns_of_interest:
    if col not in df.columns:
        missing.append(col)

if missing:
    print(f"These columns were not found in the file: {missing}")
    sys.exit(1)

# Select columns of interest and save
df_filtered = df[columns_of_interest]


# 2. CLEANING COLUMN functions
def clean_PC(cell):
    """
    Function for cleaning PC column (PruebasComplementarias)
    """
    # If the cell is empty or NaN, return an empty dictionary
    if pd.isna(cell):
        return {"hboc": None, "other_tests": None}

    # Ensure the cell is a string
    if not isinstance(cell, str):
        cell = str(cell)

    # Split the input string into tests by '@#'
    tests = cell.split('@#')

    hboc_info = None
    other_tests = []

    for test in tests:
        try: 
            test_name = test.split("Prueba: ")[1].split("@")[0].lower().strip()
        except IndexError:
            continue

        if test_name == "hboc":
        # If there are results, extract them
            hboc_info = test.split("Hallazgos:")[1].lower().strip()
        # If other tests were made
        else:
            other_tests.append(test.lower().replace("@", "-"))

    # Create 'info' dictionary
    info = {
        "hboc": hboc_info,
        "other_tests": " | ".join(other_tests) if other_tests else None
    }

    return info

def clean_DGA(cell):
    """
    Function for cleaning DGA column (DatosGeneticosAnalisis)
    """
    # If the cell is empty or NaN, return an empty dictionary
    if pd.isna(cell):
        return {}
    
    # Ensure the cell is a string
    if not isinstance(cell, str):
        cell = str(cell)

    # Split the input string into parts by '@#'
    parts = cell.split('@#')

    # Replace '@' with '-' in each part, then join them with ' | '
    clean_parts = " | ".join(part.replace("@", "-").replace(";", ",").lower().strip() for part in parts)

    # If the list is empty, returns an empty string
    # For tumor data when the tumor is unknown and for extra information
    info = {"extra_info": clean_parts if parts else ""}

    return info


def clean_TEIT(cell):
    """
    Function for cleaning TEIT column (TratamientosEIntervencionesTabla)
    """
    # If the cell is empty or NaN, return an empty dictionary
    if pd.isna(cell):
        return {}

    # Ensure the cell is a string
    if not isinstance(cell, str):
        cell = str(cell)

    # Split the input string into treatments by '@#'
    treatments = cell.split('@#')

    # Replace '@' with '-' in each treatment, then join them with ' | '
    clean_treatments = " | ".join(treatment.replace("@", "-").lower().strip() for treatment in treatments)

    # If the list is empty, returns an empty string
    info_treatments = {"treatments": clean_treatments if treatments else ""}

    return info_treatments


def clean_DATI(cell):
    """
    Function for cleaning DATI column (DiagnosisAndTumorInformation)
    """
    # If the cell is empty or NaN, return an empty dictionary
    if pd.isna(cell):
        return {}

    # Ensure the cell is a string
    if not isinstance(cell, str):
        cell = str(cell)

    # Split the input string into diagnoses by '@#'
    diagnoses = cell.split('@#')
    data_DATI = {}

    # Loops through 'diagnoses' with a counter (idx) starting at 1
    for idx, diag in enumerate(diagnoses, 1):
        # Split all 'diagnoses' into 'parts' by '@'
        parts = diag.strip().split('@')
        # Variables are initialized as empty strings
        tumor_info = {
            "type": "",
            "subtype": "",
            "diagnosis_age": "",
            "molecular_classification": "",
            "tumor_purity": "",
            "mutational_signatures": "",
            "TMB": ""
        }
        #type = subtype = diagnosis_age = molecular_classification = ""
        #tumor_purity = mutational_signatures = TMB = ""

        for part in parts:
            part = part.strip()
            # Remove the label 'Type:' and keep only the tumor 'type'
            if part.startswith("Type:"):
                tumor_info["type"] = part.replace("Type:", "").strip()
                #type = part.replace("Type:", "").strip()
            # Remove the label 'SUBTYPE:' and keep only the tumor 'subtype'
            elif part.startswith("SUBTYPE:"):
                tumor_info["subtype"] = part.replace("SUBTYPE:", "").strip()
                #subtype = part.replace("SUBTYPE:", "").strip()
            # Remove the label 'DIAGNOSIS AGE:' and keep only the number
            elif part.startswith("DIAGNOSIS AGE:"):
                tumor_info["diagnosis_age"] = part.replace("DIAGNOSIS AGE:", "").strip()
                #diagnosis_age = part.replace("DIAGNOSIS AGE:", "").strip()
            # Remove the label 'MOLECULAR CLASSIFICATION:'
            elif part.startswith("MOLECULAR CLASSIFICATION:"):
                tumor_info["molecular_classification"] = part.replace("MOLECULAR CLASSIFICATION:", "").strip()
                #molecular_classification = part.replace("MOLECULAR CLASSIFICATION:", "").strip()
            # Remove the label 'TUMOR PURITY:'
            elif part.startswith("TUMOR PURITY:"):
                tumor_info["tumor_purity"] = part.replace("TUMOR PURITY:", "").strip()
                #tumor_purity = part.replace("TUMOR PURITY:", ""). strip()
            # Remove the label 'MUTATIONAL SIGNATURES:'
            elif part.startswith("MUTATIONAL SIGNATURES:"):
                tumor_info["mutational_signatures"] = part.replace("MUTATIONAL SIGNATURES:", "").strip()
                #mutational_signatures = part.replace("MUTATIONAL SIGNATURES:", "").strip()
            # Remove the label 'TMB' and keep only the number
            elif part.startswith("TMB"):
                tumor_info["TMB"] = part.replace("TMB:", "").strip()
                #TMB = part.replace("TMB:", "").strip()

        if tumor_info["type"]:
            for k, v in tumor_info.items():
                if v:  # solo columnas no vacías
                    data_DATI[f"{k}_{idx}"] = v

        # Only append if there's a tumor type
        #if tumor_info["type"]:
         #   data_DATI.append(tumor_info)
    
        # If a 'type' was found (not empty)
        """ if type:
            data_DATI[f"type_{idx}"] = type
            if subtype:
                data_DATI[f"subtype_{idx}"] = subtype
            if diagnosis_age:
                data_DATI[f"diagnosis_age_{idx}"] = diagnosis_age
            if molecular_classification:
                data_DATI[f"molecular_classification_{idx}"] = molecular_classification
            if tumor_purity:
                data_DATI[f"tumor_purity_{idx}"] = tumor_purity
            if mutational_signatures:
                data_DATI[f"mutational_signatures_{idx}"] = mutational_signatures
            if TMB:
                data_DATI[f"TMB_{idx}"] = TMB """

    return data_DATI


# List of columns DATI and DGA
columns_to_clean = ["PruebasComplementarias","DatosGeneticosAnalisis", 
                    "TratamientosEIntervencionesTabla", "DiagnosisAndTumorInformation"] 

# Loop through the 'columns_to_clean'
for col in columns_to_clean:
    if col in df_filtered.columns:
        # Apply the clean_PC function and expand the dictionary into separate columns
        if col == "PruebasComplementarias":
            cleaned = df_filtered[col].apply(clean_PC)
            expanded_cols = pd.DataFrame(cleaned.tolist())
        # Apply the clean_DGA function and expand the dictionary into separate columns
        elif col == 'DatosGeneticosAnalisis':
            cleaned = df_filtered[col].apply(clean_DGA)
            expanded_cols = pd.DataFrame(cleaned.tolist())
        # Apply the clean_TEIT function and expand the dictionary into separate columns
        elif col == 'TratamientosEIntervencionesTabla':
            cleaned = df_filtered[col].apply(clean_TEIT)
            expanded_cols = pd.DataFrame(cleaned.tolist())
        # Apply the clean_DATI function and expand the dictionary into separate columns
        elif col == "DiagnosisAndTumorInformation":
            cleaned = df_filtered[col].apply(clean_DATI)
            expanded_cols = pd.DataFrame(cleaned.tolist())

            # Order columns per sufix
            suffixes = sorted({
                int(c.split("_")[-1])
                for c in expanded_cols.columns if "_" in c and c.split("_")[-1].isdigit()
            })
            fields_order = [
                "type",
                "subtype",
                "diagnosis_age",
                "molecular_classification",
                "tumor_purity",
                "mutational_signatures",
                "TMB"
            ]
            ordered_cols = [
                f"{field}_{suf}"
                for suf in suffixes
                for field in fields_order
                if f"{field}_{suf}" in expanded_cols.columns
            ]
            expanded_cols = expanded_cols.reindex(columns=ordered_cols)

        # Replaces the original column (col) in 'df_filtered' with the new expanded columns from the 'cleaned' data
        df_filtered = pd.concat([df_filtered.drop(columns=col), expanded_cols], axis=1)


# 3. CLEANING
# Rename columns (Spanish to English) 
column_translation = {
    "IdCaso": "case_id",
    "PerfilIndividuo": "affected",
    "Sexo": "sex",
    "CasoFamiliar": "family_case",
    "Survival": "survival",
    "hboc": "HBOC"
}

df_filtered.rename(columns=column_translation, inplace=True)


# Function for cleaning 'survival' column
def clean_survival(cell):
    if pd.isna(cell):
        return ""
    
    # Remove "Status:" and "@#". Clean spaces and lower the letters
    return cell.replace("Status:", "").replace("@#", "").strip().lower()

# Apply clean_survival function
if "survival" in df_filtered.columns:
    df_filtered["survival"] = df_filtered["survival"].apply(clean_survival)

# Remove survival data from 0445-0445, 0541-0541, 3412-3412, 3629-3629, 3672-3672, 3734-3734, 3769-3769, 4076-4076 because they are NOT FROM THE COHORT
ids_to_remove = ["0445-0445-4impact", "0541-0541-4impact", "3412-3412-4impact", "3629-3629-4impact", "3672-3672-4impact", "3734-3734-4impact", "3769-3769-4impact", "4076-4076-4impact"]
df_filtered.loc[df_filtered["case_id"].isin(ids_to_remove), "survival"] = ""

# 4. LOWERCASE all the cell values (string) except from 'mutational_signatures' and 'molecular_classification' columns
cols_to_exclude = [col for col in df_filtered.columns 
                   if (("mutational_signatures" in col or "molecular_classification" in col) and "tissue" not in col)]

# Selects only columns with type object (usually strings)
for col in df_filtered.select_dtypes(include='object').columns:
    if col not in cols_to_exclude:
        df_filtered[col] = df_filtered[col].str.lower().str.strip()

# 5. TRANSLATE column 'sex'
if 'sex' in df_filtered.columns:
    df_filtered['sex'] = df_filtered['sex'].replace('hombre', 'male')
    df_filtered['sex'] = df_filtered['sex'].replace('mujer', 'female')

# 6. REPLACE 'probando' and 'familiar afecto' for 'yes' and 'familiar' for 'no' in column 'affected', previously named 'perfil_individuo'
if 'affected' in df_filtered.columns:
    df_filtered['affected'] = df_filtered['affected'].replace('probando', 'yes')
    df_filtered['affected'] = df_filtered['affected'].replace('familiar afecto', 'yes')
    df_filtered['affected'] = df_filtered['affected'].replace('familiar', 'no')

# 7. TRANSLATE column 'family_case', previously named 'caso_familiar'
if 'family_case' in df_filtered.columns:
    df_filtered['family_case'] = df_filtered['family_case'].replace('si', 'yes')

# 8. CREATE 'family_id' column and INSERT it after 'case_id' column
df_filtered.insert(loc=1, column='family_id', value=df_filtered['case_id'].str.split('-').str[1].str.zfill(4))

# 9. SAVE cleaned csv file
df_filtered.to_csv(output_file, index=False, quoting=csv.QUOTE_NONNUMERIC)
print(f"Cleaned file saved to: {output_file}")