import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql_query("""
SELECT
    Protein_Accession,
    COUNT(*) AS Count
FROM KPPR1_Genome_Annotation
WHERE Protein_Accession IS NOT NULL
GROUP BY Protein_Accession
HAVING COUNT(*) > 1
ORDER BY Count DESC
""", conn)

print("Duplicated KPPR1 protein accessions:")
print(df)

conn.close()