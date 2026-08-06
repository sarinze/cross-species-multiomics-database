import pandas as pd
import sqlite3

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql(
    "SELECT * FROM KPHS_K56_COLISTIN_COMBINATION_DEG",
    conn
)

conn.close()

df.to_csv(
    "KPHS_K56_COLISTIN_COMBINATION_DEG_CLEAN.csv",
    index=False
)

print("Clean CSV converted successfully")
print("Rows:", len(df))
print("Columns:", len(df.columns))
print(df.columns.tolist())