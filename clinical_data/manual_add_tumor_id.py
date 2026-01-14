#!/usr/bin/env python3
"""
Script to manually add tumor_ids to cleaned_casos_with_tumor_id.csv in those cases with more than one tumor or ambiguous (no TMB base calculated)
Usage: python manual_add_TMB.py
Cristina Pastor, October 2025
"""

import pandas as pd
import os

clinical_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_tumor_id.csv"
out_clinical_csv = "/storage/scratch01/groups/bu/impact_vuscan/clinical_data/cleaned_casos_with_tumor_id_man.csv"

# Open clinical csv
print(f"Loading clinical CSV: {clinical_csv}")
df_clinical = pd.read_csv(clinical_csv)

print("Adding manual tumor_id...")
df_clinical.loc[df_clinical["case_id"] == "1884-1884-4impact", "tumor_id_1"] = "1884-1884-4impact-04"
df_clinical.loc[df_clinical["case_id"] == "1884-1884-4impact", "tumor_id_2"] = "1884-1884-4impact-02"
df_clinical.loc[df_clinical["case_id"] == "2560-2560-4impact", "tumor_id_1"] = "2560-2560-4impact-04"
df_clinical.loc[df_clinical["case_id"] == "2560-2560-4impact", "tumor_id_2"] = "2560-2560-4impact-02"
df_clinical.loc[df_clinical["case_id"] == "3018-3018-4impact", "tumor_id_1"] = "3018-3018-4impact-02"
df_clinical.loc[df_clinical["case_id"] == "3700-3700-4impact", "tumor_id_1"] = "3700-3700-4impact-04"
df_clinical.loc[df_clinical["case_id"] == "3700-3700-4impact", "tumor_id_2"] = "3700-3700-4impact-02"
df_clinical.loc[df_clinical["case_id"] == "1837-1837-4impact", "tumor_id_2"] = "1837-1837-4impact-02"
df_clinical.loc[df_clinical["case_id"] == "1837-1837-4impact", "tumor_id_1"] = "1837-1837-4impact-04"
df_clinical.loc[df_clinical["case_id"] == "2748-2748-4impact", "tumor_id_1"] = "2748-2748-4impact-04"
df_clinical.loc[df_clinical["case_id"] == "2748-2748-4impact", "tumor_id_2"] = "2748-2748-4impact-02"
df_clinical.loc[df_clinical["case_id"] == "2079-2077-4impact", "tumor_id_2"] = "2079-2077-4impact-02"

print(f"Saving updated CSV in {out_clinical_csv}...")
df_clinical.to_csv(out_clinical_csv, index=False)


print("##### PIPELINE COMPLETED #####")