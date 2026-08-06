import pandas as pd
import sqlite3

FILES = [
    ("PA14", "MGH78578", "PA14_MGH78578_ANNOTATED_LINKER.csv"),
    ("PA14", "NCTC5055", "PA14_NCTC5055_ANNOTATED_LINKER.csv"),
    ("PA14", "HS11286", "PA14_HS11286_ANNOTATED_LINKER.csv"),
    ("PAO1", "MGH78578", "PAO1_MGH78578_ANNOTATED_LINKER.csv"),
    ("PAO1", "NCTC5055", "PAO1_NCTC5055_ANNOTATED_LINKER.csv"),
    ("PAO1", "HS11286", "PAO1_HS11286_ANNOTATED_LINKER.csv"),
]

def get_col(df, possible_names):
    for name in possible_names:
        if name in df.columns:
            return df[name]
    return pd.Series([None] * len(df))

all_tables = []

for pseudo_strain, kleb_strain, filename in FILES:
    print(f"Loading {filename}")
    df = pd.read_csv(filename)

    out = pd.DataFrame()
    out["Pseudomonas_Strain"] = [pseudo_strain] * len(df)
    out["Klebsiella_Strain"] = [kleb_strain] * len(df)

    if pseudo_strain == "PA14":
        out["Pseudomonas_Locus_Tag"] = get_col(df, ["PA14_New_Locus_Tag"])
        out["Pseudomonas_Original_Locus_Tag"] = get_col(df, ["PA14_Original_Locus_Tag"])
        out["Pseudomonas_Gene_Name"] = get_col(df, ["PA14_Gene_Name"])
        out["Pseudomonas_Symbol"] = get_col(df, ["PA14_Symbol"])
        out["Pseudomonas_Gene_Type"] = get_col(df, ["PA14_Gene_Type"])
        out["Pseudomonas_Protein_Name"] = get_col(df, ["PA14_Protein_Name"])
        out["Pseudomonas_Protein_Accession"] = get_col(df, ["PA14_Protein_Accession"])
        out["Pseudomonas_Protein_Length"] = get_col(df, ["PA14_Protein_Length"])
    else:
        out["Pseudomonas_Locus_Tag"] = get_col(df, ["PAO1_Locus_Tag"])
        out["Pseudomonas_Original_Locus_Tag"] = [None] * len(df)
        out["Pseudomonas_Gene_Name"] = get_col(df, ["PAO1_Gene_Name"])
        out["Pseudomonas_Symbol"] = get_col(df, ["PAO1_Symbol"])
        out["Pseudomonas_Gene_Type"] = get_col(df, ["PAO1_Gene_Type"])
        out["Pseudomonas_Protein_Name"] = get_col(df, ["PAO1_Protein_Name"])
        out["Pseudomonas_Protein_Accession"] = get_col(df, ["PAO1_Protein_Accession"])
        out["Pseudomonas_Protein_Length"] = get_col(df, ["PAO1_Protein_Length"])

    out["Klebsiella_Locus_Tag"] = get_col(df, [f"{kleb_strain}_Locus_Tag"])
    out["Klebsiella_Gene_Name"] = get_col(df, [f"{kleb_strain}_Gene_Name"])
    out["Klebsiella_Symbol"] = get_col(df, [f"{kleb_strain}_Symbol"])
    out["Klebsiella_Gene_Type"] = get_col(df, [f"{kleb_strain}_Gene_Type"])
    out["Klebsiella_Protein_Name"] = get_col(df, [f"{kleb_strain}_Protein_Name"])
    out["Klebsiella_Protein_Accession"] = get_col(df, [f"{kleb_strain}_Protein_Accession"])
    out["Klebsiella_Protein_Length"] = get_col(df, [f"{kleb_strain}_Protein_Length"])

    out["Percent_Identity"] = get_col(df, [
        f"pident_{pseudo_strain}_to_{kleb_strain}",
        f"pident{pseudo_strain}_to_{kleb_strain}",
    ])

    out["Alignment_Length"] = get_col(df, [
        f"length_{pseudo_strain}_to_{kleb_strain}",
        f"length{pseudo_strain}_to_{kleb_strain}",
    ])

    out["Evalue"] = get_col(df, [
        f"evalue_{pseudo_strain}_to_{kleb_strain}",
        f"evalue{pseudo_strain}_to_{kleb_strain}",
    ])

    out["Bitscore"] = get_col(df, [
        f"bitscore_{pseudo_strain}_to_{kleb_strain}",
        f"bitscore{pseudo_strain}_to_{kleb_strain}",
    ])

    out["Source_File"] = [filename] * len(df)

    all_tables.append(out)

master = pd.concat(all_tables, ignore_index=True)

master.to_csv("KLEBSIELLA_PSEUDOMONAS_MASTER_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
master.to_sql(
    "KLEBSIELLA_PSEUDOMONAS_MASTER_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("\nMASTER LINKER CREATED")
print("Rows:", len(master))

print("\nPseudomonas strains:")
print(master["Pseudomonas_Strain"].unique())

print("\nKlebsiella strains:")
print(master["Klebsiella_Strain"].unique())

print("\nRows by comparison:")
print(master.groupby(["Pseudomonas_Strain", "Klebsiella_Strain"]).size())

print("\nSaved CSV: KLEBSIELLA_PSEUDOMONAS_MASTER_LINKER.csv")
print("Imported SQLite Table: KLEBSIELLA_PSEUDOMONAS_MASTER_LINKER")