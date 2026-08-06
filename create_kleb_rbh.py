import pandas as pd
import sqlite3
cols = [
    "query",
    "subject",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
]
mgh_vs_nctc = pd.read_csv("MGH78578_vs_NCTC5055_REFSEQ.tsv", sep="\t", names=cols)
nctc_vs_mgh = pd.read_csv("NCTC5055_vs_MGH78578_REFSEQ.tsv", sep="\t", names=cols)

mgh_best = (
    mgh_vs_nctc.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)
nctc_best = (
    nctc_vs_mgh.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)
rbh = mgh_best.merge(
    nctc_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("_MGH78578_to_NCTC5055", "_NCTC5055_to_MGH78578")
)
orthologs = rbh[[
    "query_MGH78578_to_NCTC5055",
    "subject_MGH78578_to_NCTC5055",
    "pident_MGH78578_to_NCTC5055",
    "length_MGH78578_to_NCTC5055",
    "evalue_MGH78578_to_NCTC5055",
    "bitscore_MGH78578_to_NCTC5055"
]]
orthologs.columns = [
    "MGH78578_SEQUENCE",
    "NCTC5055_SEQUENCE",
    "PERCENT_IDENTITY",
    "ALIGNMENT_LENGTH",
    "EVALUE",
    "BITSCORE"
]
orthologs.to_csv("MGH78578_NCTC5055_PROT_RBH_Orthologs.csv", index=False)

print("Number of reciprocal best hits:", len(orthologs))
print(orthologs.head())
print("Saved file: MGH78578_NCTC5055_PROT_RBH_Orthologs.csv")
conn = sqlite3.connect("pseudomonas_new.db")
orthologs.to_sql("MGH78578_NCTC5055_PROT_RBH_LINKER", conn, if_exists="replace", index=False)
conn.close()
print("Imported into SQLite as table: MGH78578_NCTC5055_PROT_RBH_LINKER")


mgh_annot = pd.read_csv("Kleb_MGH_78578_Genome_Annotation.csv")
nctc_annot = pd.read_csv("NCTC5055_Genome_Annotation.csv")

final_linker = orthologs.merge(
    mgh_annot,
    left_on="MGH78578_SEQUENCE",
    right_on="Protein_Accession",
    how="left"
)
final_linker = final_linker.merge(
    nctc_annot,
    left_on="NCTC5055_SEQUENCE",
    right_on="Protein_Accession",
    how="left",
    suffixes=("_MGH78578", "_NCTC5055")
    )

final_linker.to_csv("MGH78578_NCTC5055_Final_Locus_Linker.csv", index=False)

print("Saved final linker table:")
print("MGH78578_NCTC5055_Final_Locus_Linker.csv")

conn = sqlite3.connect("pseudomonas_new.db")

final_linker.to_sql("MGH78578_NCTC5055_FINAL_LINKER", conn, if_exists="replace", index=False)
conn.close()
print("Imported into SQLite as:")
print("MGH78578_NCTC5055_FINAL_LINKER")

print("MGH columns:")
print(mgh_annot.columns.tolist())

print("\nNCTC columns:")
print(nctc_annot.columns.tolist())

print("\nMGH Protein_Accession preview:")
print(mgh_annot["Protein_Accession"].head(10))

print("\nNCTC Protein_Accession preview:")
print(nctc_annot["Protein_Accession"].head(10))


def make_rbh(forward_file, reverse_file, strain_a, strain_b, output_csv, sqlite_table):
    print(f"\nRunning RBH for {strain_a} vs {strain_b}")
    forward = pd.read_csv(forward_file, sep="\t", names=cols)
    reverse = pd.read_csv(reverse_file, sep="\t", names=cols)
    forward_best = (
        forward.sort_values("bitscore", ascending=False)
        .drop_duplicates("query")
    )
    reverse_best = (
        reverse.sort_values("bitscore", ascending=False)
        .drop_duplicates("query")
    )
    rbh = forward_best.merge(
        reverse_best,
        left_on=["query", "subject"],
        right_on=["subject", "query"],
        suffixes=(f"_{strain_a}_to_{strain_b}", f"_{strain_b}_to_{strain_a}")
    )
    print("Forward best hits:", len(forward_best))
    print("Reverse best hits:", len(reverse_best))
    print("Reciprocal best hits after merge:", len(rbh))
    orthologs = rbh[[
        f"query_{strain_a}_to_{strain_b}",
        f"subject_{strain_a}_to_{strain_b}",
        f"pident_{strain_a}_to_{strain_b}",
        f"length_{strain_a}_to_{strain_b}",
        f"evalue_{strain_a}_to_{strain_b}",
        f"bitscore_{strain_a}_to_{strain_b}"
    ]]
    orthologs.columns = [
        f"{strain_a}_SEQUENCE",
        f"{strain_b}_SEQUENCE",
        "PERCENT_IDENTITY",
        "ALIGNMENT_LENGTH",
        "EVALUE",
        "BITSCORE"
    ]
    orthologs.to_csv(output_csv, index=False)
    conn = sqlite3.connect("pseudomonas_new.db")
    orthologs.to_sql(sqlite_table, conn, if_exists="replace", index=False)
    conn.close()
    print("Number of reciprocal best hits:", len(orthologs))
    print("Saved file:", output_csv)
    print("Imported into SQLite as table:", sqlite_table)
    print(orthologs.head())

make_rbh(
    forward_file="Kleb_HS11286_vs_NCTC5055_REFSEQ.tsv",
    reverse_file="NCTC5055_vs_Kleb_HS11286_REFSEQ.tsv",
    strain_a="Kleb_HS11286",
    strain_b="NCTC5055",
    output_csv="Kleb_HS11286_NCTC5055_PROT_RBH_Orthologs.csv",
    sqlite_table="Kleb_HS11286_NCTC5055_PROT_RBH_LINKER",
)
make_rbh(
    forward_file="Kleb_HS11286_vs_MGH78578_REFSEQ.tsv",
    reverse_file="MGH78578_vs_Kleb_HS11286_REFSEQ.tsv",
    strain_a="Kleb_HS11286",
    strain_b="MGH78578",
    output_csv="Kleb_HS11286_MGH78578_PROT_RBH_Orthologs.csv",
    sqlite_table="Kleb_HS11286_MGH78578_PROT_RBH_LINKER",
)