import sqlite3
import pandas as pd

DB = "pseudomonas_new.db"

INPUT_TABLE = "Master_Multiomics_Matrix_With_KPPR1_TnSeq"
OUTPUT_TABLE = "Master_Multiomics_Matrix_With_Polymyxin_Rifampicin"
OUTPUT_CSV = "Master_Multiomics_Matrix_With_Polymyxin_Rifampicin.csv"

PAO1_ID = "Pseudomonas_Locus_Tag"

datasets = {
    "PAO1_PolymycinB_Rifampicin_S5_PolymyxinB_1h": "PB_1h",
    "PAO1_PolymycinB_Rifampicin_S6_PolymyxinB_4h": "PB_4h",
    "PAO1_PolymycinB_Rifampicin_S7_PolymyxinB_24h": "PB_24h",
    "PAO1_PolymycinB_Rifampicin_S8_PolymyxinB_Rifampicin_1h": "PB_RIF_1h",
    "PAO1_PolymycinB_Rifampicin_S9_PolymyxinB_Rifampicin_4h": "PB_RIF_4h",
    "PAO1_PolymycinB_Rifampicin_S10_PolymyxinB_Rifampicin_24h": "PB_RIF_24h",
}

conn = sqlite3.connect(DB)

master = pd.read_sql(f"SELECT * FROM {INPUT_TABLE}", conn)

print("========== INPUT MASTER ==========")
print("Rows:", len(master))
print("Columns:", len(master.columns))

master[PAO1_ID] = master[PAO1_ID].astype(str).str.strip()

final = master.copy()

for table, prefix in datasets.items():
    print(f"\n========== ADDING {prefix} ==========")

    df = pd.read_sql(f"SELECT * FROM {table}", conn)

    # keep only real columns, removing Unnamed garbage columns
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

    df["Locus_Tag"] = df["Locus_Tag"].astype(str).str.strip()

    keep = df[["Locus_Tag", "Log2FC", "adj_P_Value"]].copy()

    keep = keep.rename(
        columns={
            "Log2FC": f"PAO1_{prefix}_Log2FC",
            "adj_P_Value": f"PAO1_{prefix}_Adjusted_P_Value",
        }
    )

    keep = keep.drop_duplicates(subset=["Locus_Tag"])

    before_cols = len(final.columns)

    final = final.merge(
        keep,
        left_on=PAO1_ID,
        right_on="Locus_Tag",
        how="left",
    )

    final = final.drop(columns=["Locus_Tag"])

    print("Rows after merge:", len(final))
    print("Columns added:", len(final.columns) - before_cols)
    print("Non-null Log2FC:", final[f"PAO1_{prefix}_Log2FC"].notna().sum())
    print("Non-null adjusted p-value:", final[f"PAO1_{prefix}_Adjusted_P_Value"].notna().sum())

print("\n========== FINAL QC ==========")
print("Final rows:", len(final))
print("Final columns:", len(final.columns))

new_cols = [c for c in final.columns if c.startswith("PAO1_PB")]
print("\nNew polymyxin/rifampicin columns:")
print(new_cols)

print("\nNon-null counts:")
print(final[new_cols].notna().sum())

final.to_csv(OUTPUT_CSV, index=False)
final.to_sql(OUTPUT_TABLE, conn, if_exists="replace", index=False)

conn.close()

print("\n========== COMPLETE ==========")
print("Output CSV:", OUTPUT_CSV)
print("SQLite table:", OUTPUT_TABLE)