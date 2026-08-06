import pandas as pd

INPUT_FILE = "KPPR1_RefSeq_Protein_FASTA_Clean.csv"
OUTPUT_FILE = "KPPR1_proteins.faa"

df = pd.read_csv(INPUT_FILE)

print("Columns:")
print(df.columns.tolist())
print("\nFirst rows:")
print(df.head())

# Adjust these names if your printed columns are different
id_col = "Protein_Accession"
seq_col = "Protein_Sequence"

df = df.dropna(subset=[id_col, seq_col])

with open(OUTPUT_FILE, "w") as f:
    for _, row in df.iterrows():
        protein_id = str(row[id_col]).strip()
        seq = str(row[seq_col]).strip().replace(" ", "").replace("\n", "")
        f.write(f">{protein_id}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

print("\nFASTA created successfully")
print("Output file:", OUTPUT_FILE)
print("Sequences:", len(df))