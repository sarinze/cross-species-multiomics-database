import pandas as pd
import sqlite3
import re
from pathlib import Path

INPUT_FILE = "KPPR1_Protein_faa_RefSeq.csv"
DB_FILE = "pseudomonas_new.db"
TABLE_NAME = "kppr1_refseq_protein_fasta"

if not Path(INPUT_FILE).exists():
    raise FileNotFoundError(f"Could not find file: {INPUT_FILE}")

raw = pd.read_csv(INPUT_FILE)

lines = raw.iloc[:, 0].dropna().astype(str).str.strip().tolist()

rows = []
current_header = None
current_sequence = []

def make_row(header, sequence_lines):
    sequence = "".join(sequence_lines)

    parts = header.split(" ", 1)
    protein_accession = parts[0].replace(">", "").strip()
    description = parts[1].strip() if len(parts) > 1 else ""

    organism_match = re.search(r"\[(.*?)\]", description)
    organism = organism_match.group(1) if organism_match else ""

    protein_name = re.sub(r"\s*\[.*?\]\s*$", "", description).strip()

    return {
        "Protein_Accession": protein_accession,
        "Protein_Name": protein_name,
        "Organism": organism,
        "Protein_Length": len(sequence),
        "Protein_Sequence": sequence,
        "Full_FASTA_Header": header
    }

for line in lines:
    if line.startswith(">"):
        if current_header is not None:
            rows.append(make_row(current_header, current_sequence))
        current_header = line
        current_sequence = []
    else:
        current_sequence.append(line)

if current_header is not None:
    rows.append(make_row(current_header, current_sequence))

df = pd.DataFrame(rows)

print(df.head())
print()
print("Rows:", len(df))
print("Columns:", list(df.columns))

df.to_csv("KPPR1_RefSeq_Protein_FASTA_Clean.csv", index=False)

conn = sqlite3.connect(DB_FILE)
df.to_sql(TABLE_NAME, conn, if_exists="replace", index=False)
conn.close()

print()
print("Saved CSV: KPPR1_RefSeq_Protein_FASTA_Clean.csv")
print(f"Imported SQLite table: {TABLE_NAME}")