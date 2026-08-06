import pandas as pd
import sqlite3

file_name = "KP13_PolymyxinB_DEG_Table.csv"
db_name = "pseudomonas_new.db"
table_name = "KP13_PolymyxinB_DEG_Table"

df = pd.read_csv("KP13_PolymyxinB_DEG_Table.csv", skiprows=1, usecols=range(8))

# Clean column names
df.columns = (
    df.columns
    .str.strip()
    .str.replace(" ", "_")
    .str.replace("-", "_")
    .str.replace("/", "_")
)
print(df.head())
print("\nRows:", len(df))
print("\nColumns:")
print(df.columns.tolist())

conn = sqlite3.connect("pseudomanas_new.db")

df.to_sql(
    "KP13_PolymyxinB_DEG_Table",
    conn, if_exists="replace", index=False
)
conn.close()

print("KP13_PolymyxinB_DEG_Table imported successfully")