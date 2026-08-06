import pandas as pd

cols = [
    "query", "subject", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

forward = pd.read_csv("KPPR1_vs_MGH78578.tsv", sep="\t", names=cols)
reverse = pd.read_csv("MGH78578_vs_KPPR1.tsv", sep="\t", names=cols)

print("Forward BLAST: KPPR1 -> MGH78578")
print(forward.head(20)[["query", "subject", "pident", "bitscore"]])

print("\nReverse BLAST: MGH78578 -> KPPR1")
print(reverse.head(20)[["query", "subject", "pident", "bitscore"]])

print("\nForward query examples:")
print(forward["query"].drop_duplicates().head(10).tolist())

print("\nForward subject examples:")
print(forward["subject"].drop_duplicates().head(10).tolist())

print("\nReverse query examples:")
print(reverse["query"].drop_duplicates().head(10).tolist())

print("\nReverse subject examples:")
print(reverse["subject"].drop_duplicates().head(10).tolist())