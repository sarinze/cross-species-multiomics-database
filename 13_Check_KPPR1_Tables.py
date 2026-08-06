import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

tables = pd.read_sql_query("""
SELECT name
FROM sqlite_master
WHERE type='table'
ORDER BY name
""", conn)

print(tables.to_string(index=False))

conn.close()