import pandas as pd
import sqlite3

csv_file = "PAO1_Tobramycin_TnSeq.csv"
db_file = "pseudomonas_new.db"
table_name = "PAO1_Tobramycin_TnSeq"

df = pd.read_csv(csv_file, skiprows=4)

df.columns = [
    "PAO1_Locus_Tag",
    "Gene",
    "Function",
    "Gene_Length_bp",
    "Pre_Growth_Hits",
    "Pre_Growth_Reads",
    "Growth_No_Tobramycin_Hits",
    "Growth_No_Tobramycin_Reads",
    "Growth_With_Tobramycin_Hits",
    "Growth_With_Tobramycin_Reads",
    "Selection_Ratio",
    "Previous_Study_Identified",
    "Previous_Study_MIC",
    "This_Study_Mutants_Available",
    "This_Study_Mutants_Tested",
    "This_Study_MIC"
]

conn = sqlite3.connect(db_file)
df.to_sql(table_name, conn, if_exists="replace", index=False)
conn.close()

print("Imported successfully:", table_name)
print("Rows:", len(df))
print("Columns:", list(df.columns))