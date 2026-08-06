import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

table = "KLEBSIELLA_PSEUDOMONAS_MASTER_LINKER"

df = pd.read_sql_query(f"SELECT * FROM {table} LIMIT 10", conn)

print("Columns:")
print(df.columns.tolist())

print("\nFirst 10 rows:")
print(df)

print("\nTotal rows:")
count = pd.read_sql_query(f"SELECT COUNT(*) AS row_count FROM {table}", conn)
print(count)

conn.close()