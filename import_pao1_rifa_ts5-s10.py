import pandas as pd
import sqlite3
import re

DB = "pseudomonas_new.db"

files = {
    "S5_PolymyxinB_1h": "PAO1_PolymycinB_Rifampicin_Table_S5.csv",
    "S6_PolymyxinB_4h": "PAO1_PolymycinB_Rifampicin_Table_S6.csv",
    "S7_PolymyxinB_24h": "PAO1_PolymycinB_Rifampicin_Table_S7.csv",
    "S8_PolymyxinB_Rifampicin_1h": "PAO1_PolymycinB_Rifampicin_Table_S8.csv",
    "S9_PolymyxinB_Rifampicin_4h": "PAO1_PolymycinB_Rifampicin_Table_S9.csv",
    "S10_PolymyxinB_Rifampicin_24h": "PAO1_PolymycinB_Rifampicin_Table_S10.csv",
}

def clean_col(col):
    col = str(col).strip()
    col = re.sub(r"\s+", "_", col)
    col = col.replace("-", "_")
    col = col.replace("/", "_")
    col = col.replace("(", "")
    col = col.replace(")", "")
    col = col.replace(".", "")
    return col

conn = sqlite3.connect(DB)

for table_name, file in files.items():
    print(f"\nImporting {file}")

    df = pd.read_csv(file)
    df.columns = [clean_col(c) for c in df.columns]

    sqlite_table = "PAO1_PolymycinB_Rifampicin_" + table_name

    df.to_sql(sqlite_table, conn, if_exists="replace", index=False)

    print("Saved table:", sqlite_table)
    print("Rows:", len(df))
    print("Columns:", df.columns.tolist())

conn.close()

print("\nImport complete.")