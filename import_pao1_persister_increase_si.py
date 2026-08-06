import pandas as pd
import sqlite3

csv_file = "PAO1_Persister_Cell_TnSeq_10fold_Increase_Survival_Index.csv"
db_file = "pseudomonas_new.db"
table_name = "PAO1_Persister_TnSeq_10fold_Increase_SI"

df = pd.read_csv(csv_file, skiprows=1)

df.columns = [
    "PAO1_Locus_Tag",
    "Gene",
    "Protein_ID",
    "Read_Count",
    "Rep1_Input_Sites",
    "Rep1_Cipro_Sites",
    "Rep2_Input_Sites",
    "Rep2_Cipro_Sites",
    "Rep3_Input_Sites",
    "Rep3_Cipro_Sites",
    "Rep1_Survival_Index",
    "Rep2_Survival_Index",
    "Rep3_Survival_Index",
    "Mean_Survival_Index"
]

conn = sqlite3.connect(db_file)
df.to_sql(table_name, conn, if_exists="replace", index=False)
conn.close()

print("Imported successfully:", table_name)
print("Rows:", len(df))
print("Columns:", list(df.columns))