import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql_query("""
SELECT
    Protein_Accession,
    Protein_Name,
    Full_FASTA_Header
FROM kppr1_protein_fasta
LIMIT 30
""", conn)

conn.close()

pd.set_option("display.max_colwidth", 300)
print(df)