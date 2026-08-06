import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql_query("""
SELECT *
FROM KPPR1_Genome_Annotation
LIMIT 10
""", conn)

print("Columns:")
print(df.columns.tolist())

print("\nFirst 10 rows:")
print(df)

conn.close()