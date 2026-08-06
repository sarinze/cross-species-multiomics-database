import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

# Total rows in SQLite
total = pd.read_sql_query("""
SELECT COUNT(*) AS Total_Rows
FROM Klebsiella_TnSeq
""", conn)

# Distinct Gene_ID values
distinct = pd.read_sql_query("""
SELECT COUNT(DISTINCT Gene_ID) AS Unique_Genes
FROM Klebsiella_TnSeq
""", conn)

# Check for duplicate Gene_IDs
duplicates = pd.read_sql_query("""
SELECT
    Gene_ID,
    COUNT(*) AS Number_of_Rows
FROM Klebsiella_TnSeq
GROUP BY Gene_ID
HAVING COUNT(*) > 1
ORDER BY Number_of_Rows DESC
LIMIT 20
""", conn)

print("Total rows:")
print(total)

print("\nUnique Gene_IDs:")
print(distinct)

print("\nGenes appearing more than once:")
print(duplicates)

conn.close()