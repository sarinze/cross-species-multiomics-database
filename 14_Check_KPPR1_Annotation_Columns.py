import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql_query("""
SELECT *
FROM KPPR1_Genome_Annotation
LIMIT 20
""", conn)

conn.close()

print("Columns:")
print(list(df.columns))

print("\nFirst 20 rows:")
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
print(df)