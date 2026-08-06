import sqlite3
import pandas as pd

DB = "pseudomonas_new.db"

conn = sqlite3.connect(DB)

tables = pd.read_sql(
    "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;",
    conn
)

print("========== TABLES ==========")
print(tables.to_string(index=False))

for table in tables["name"]:
    if "polymy" in table.lower() or "rif" in table.lower():
        print(f"\n========== {table} ==========")
        df = pd.read_sql(f"SELECT * FROM {table} LIMIT 5", conn)
        print(df.head())
        print("Columns:")
        print(df.columns.tolist())

conn.close()