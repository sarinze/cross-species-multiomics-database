import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

mgh = pd.read_sql_query("""
SELECT
    Protein_Accession,
    Protein_Name
FROM Kleb_MGH_78578_Genome_Annotation
WHERE Protein_Accession IS NOT NULL
""", conn)

print("Rows with protein accessions:", len(mgh))
print(mgh.head())

conn.close()

print("\nThis table has accessions, but not sequences.")
print("So we need to recover the WP_ FASTA file used to build MGH78578_CLEAN_DB, or re-download/export the RefSeq protein FASTA.")
