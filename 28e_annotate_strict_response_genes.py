# ======================================================================
# STEP 28E: PA14-CENTRED ANNOTATION OF STRICT CONSERVED/DIVERGENT GENES
# Uses Final_Cluster_Gene_Catalogue.csv as the explicit annotation source
# ======================================================================

from pathlib import Path
import pandas as pd
import numpy as np


print("\n" + "=" * 78)
print("STEP 28E: PA14-CENTRED STRICT RESPONSE GENE ANNOTATION")
print("=" * 78 + "\n")


# ----------------------------------------------------------------------
# 1. FILE PATHS
# ----------------------------------------------------------------------

PROJECT_DIRECTORY = Path.cwd()

STRICT_DIRECTORY = (
    PROJECT_DIRECTORY /
    "Strict_Conserved_Response_Analysis_Outputs"
)

STRICT_CONSERVED_FILE = (
    STRICT_DIRECTORY /
    "Strict_Conserved_Response_Genes.csv"
)

STRICT_DIVERGENT_FILE = (
    STRICT_DIRECTORY /
    "Strict_Divergent_Response_Genes.csv"
)

GENE_CATALOGUE_FILE = (
    PROJECT_DIRECTORY /
    "Final_Cluster_Gene_Catalogue.csv"
)

OUTPUT_DIRECTORY = (
    PROJECT_DIRECTORY /
    "Strict_PA14_Annotated_Response_Outputs"
)

OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# 2. OUTPUT FILES
# ----------------------------------------------------------------------

CONSERVE_OUTPUT = (
    OUTPUT_DIRECTORY /
    "Strict_Conserved_Response_Genes_PA14_Annotated.csv"
)

DIVERGENT_OUTPUT = (
    OUTPUT_DIRECTORY /
    "Strict_Divergent_Response_Genes_PA14_Annotated.csv"
)

COMBINED_OUTPUT = (
    OUTPUT_DIRECTORY /
    "Strict_Combined_Response_Genes_PA14_Annotated.csv"
)

UNMATCHED_OUTPUT = (
    OUTPUT_DIRECTORY /
    "Strict_Response_Genes_Unmatched_to_Catalogue.csv"
)

EXCEL_OUTPUT = (
    OUTPUT_DIRECTORY /
    "Strict_PA14_Annotated_Response_Genes.xlsx"
)


# ----------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ----------------------------------------------------------------------

def clean_identifier(series):
    """Standardise locus tags before matching."""
    return (
        series.astype("string")
        .str.strip()
        .replace(
            {
                "": pd.NA,
                "nan": pd.NA,
                "None": pd.NA,
                "<NA>": pd.NA,
            }
        )
    )


def first_non_missing(series):
    """Return the first available value in a group."""
    values = series.dropna()

    if len(values) == 0:
        return pd.NA

    return values.iloc[0]


def combine_gene_label(row):
    """
    Prefer the existing Gene_Label from the catalogue.
    If absent, construct:
    PA14 locus tag | gene name | Klebsiella locus tag
    """

    existing = row.get("Gene_Label", pd.NA)

    if pd.notna(existing):
        existing = str(existing).strip()

        if existing and existing.lower() not in {"nan", "none"}:
            return existing

    values = []

    for column in [
        "Pseudomonas_Original_Locus_Tag",
        "Pseudomonas_Gene_Name",
        "Klebsiella_Locus_Tag",
    ]:
        value = row.get(column, pd.NA)

        if pd.notna(value):
            value = str(value).strip()

            if (
                value
                and value.lower() not in {"nan", "none"}
                and value not in values
            ):
                values.append(value)

    return " | ".join(values) if values else "Unlabelled orthologue"


def colistin_status(row):
    """Create a readable strict colistin-response description."""

    response_class = str(
        row.get("Cross_Species_Response_Class", "")
    ).strip()

    pseudomonas = pd.to_numeric(
        pd.Series([
            row.get(
                "Pseudomonas_Colistin_Consensus_Log2FC",
                np.nan
            )
        ]),
        errors="coerce"
    ).iloc[0]

    klebsiella = pd.to_numeric(
        pd.Series([
            row.get(
                "Klebsiella_Colistin_Log2FC",
                np.nan
            )
        ]),
        errors="coerce"
    ).iloc[0]

    if response_class == "Conserved":
        if pd.notna(pseudomonas) and pseudomonas > 0:
            return "Conserved induced"

        if pd.notna(pseudomonas) and pseudomonas < 0:
            return "Conserved repressed"

        return "Conserved"

    if response_class == "Divergent":
        if pd.notna(pseudomonas) and pd.notna(klebsiella):
            if pseudomonas > 0 and klebsiella < 0:
                return (
                    "Divergent: induced in PA14 / "
                    "repressed in K56"
                )

            if pseudomonas < 0 and klebsiella > 0:
                return (
                    "Divergent: repressed in PA14 / "
                    "induced in K56"
                )

        return "Divergent"

    return response_class or "Not classified"


# ----------------------------------------------------------------------
# 4. CHECK INPUT FILES
# ----------------------------------------------------------------------

required_files = [
    STRICT_CONSERVED_FILE,
    STRICT_DIVERGENT_FILE,
    GENE_CATALOGUE_FILE,
]

missing_files = [
    str(path)
    for path in required_files
    if not path.exists()
]

if missing_files:
    raise FileNotFoundError(
        "\nThe following required files were not found:\n  "
        + "\n  ".join(missing_files)
    )


# ----------------------------------------------------------------------
# 5. LOAD STRICT RESULTS
# ----------------------------------------------------------------------

strict_conserved = pd.read_csv(
    STRICT_CONSERVED_FILE,
    low_memory=False
)

strict_divergent = pd.read_csv(
    STRICT_DIVERGENT_FILE,
    low_memory=False
)

strict_conserved["Strict_Analysis_Group"] = "Conserved"
strict_divergent["Strict_Analysis_Group"] = "Divergent"

strict_combined = pd.concat(
    [strict_conserved, strict_divergent],
    ignore_index=True,
    sort=False
)

print(
    f"Loaded strict conserved genes: "
    f"{len(strict_conserved)}"
)

print(
    f"Loaded strict divergent genes: "
    f"{len(strict_divergent)}"
)


# ----------------------------------------------------------------------
# 6. LOAD THE FINAL GENE CATALOGUE
# ----------------------------------------------------------------------

catalogue = pd.read_csv(
    GENE_CATALOGUE_FILE,
    low_memory=False
)

required_catalogue_columns = [
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "Gene_Label",
    "Final_Cluster",
]

missing_catalogue_columns = [
    column
    for column in required_catalogue_columns
    if column not in catalogue.columns
]

if missing_catalogue_columns:
    raise KeyError(
        "\nThe catalogue is missing these required columns:\n  "
        + "\n  ".join(missing_catalogue_columns)
    )

if "Klebsiella_Locus_Tag" not in strict_combined.columns:
    raise KeyError(
        "\nThe strict output files do not contain "
        "'Klebsiella_Locus_Tag'."
    )

print(
    f"Loaded Final_Cluster_Gene_Catalogue.csv: "
    f"{len(catalogue)} rows"
)


# ----------------------------------------------------------------------
# 7. CLEAN MATCHING IDENTIFIERS
# ----------------------------------------------------------------------

strict_combined["Klebsiella_Locus_Tag"] = clean_identifier(
    strict_combined["Klebsiella_Locus_Tag"]
)

catalogue["Klebsiella_Locus_Tag"] = clean_identifier(
    catalogue["Klebsiella_Locus_Tag"]
)

catalogue["Pseudomonas_Original_Locus_Tag"] = clean_identifier(
    catalogue["Pseudomonas_Original_Locus_Tag"]
)


# ----------------------------------------------------------------------
# 8. CHECK PA14 ANNOTATION AVAILABILITY
# ----------------------------------------------------------------------

catalogue_rows_with_pa14 = int(
    catalogue["Pseudomonas_Original_Locus_Tag"]
    .notna()
    .sum()
)

catalogue_rows_without_pa14 = int(
    catalogue["Pseudomonas_Original_Locus_Tag"]
    .isna()
    .sum()
)

print(
    f"Catalogue rows with Pseudomonas locus tags: "
    f"{catalogue_rows_with_pa14}"
)

print(
    f"Catalogue rows without Pseudomonas locus tags: "
    f"{catalogue_rows_without_pa14}"
)


# ----------------------------------------------------------------------
# 9. COLLAPSE CATALOGUE TO ONE ROW PER KLEBSIELLA LOCUS TAG
# ----------------------------------------------------------------------

annotation_columns = [
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Gene_Name",
    "Gene_Label",
    "Final_Cluster",
]

# Include useful measurements from the catalogue only when present.
optional_catalogue_columns = [
    "Colistin_RNAseq_Fold_Change",
    "Colistin_Riboseq_Fold_Change",
    "Tobramycin_RNAseq_Fold_Change",
    "Tobramycin_Riboseq_Fold_Change",
    "PAO1_Persister_All_Mean_SI",
    "K56_vs_Colistin_Log2FC",
    "K56_vs_Combination_Log2FC",
    "Colistin_vs_Combination_Log2FC",
    "KPPR1_TnSeq_Log2FC_Output_Input",
    "Measured_Omics_Count",
]

annotation_columns += [
    column
    for column in optional_catalogue_columns
    if column in catalogue.columns
]

catalogue_subset = catalogue[
    ["Klebsiella_Locus_Tag"] + annotation_columns
].copy()

catalogue_subset = catalogue_subset[
    catalogue_subset["Klebsiella_Locus_Tag"].notna()
].copy()

aggregation = {
    column: first_non_missing
    for column in annotation_columns
}

catalogue_lookup = (
    catalogue_subset
    .groupby(
        "Klebsiella_Locus_Tag",
        as_index=False,
        dropna=False
    )
    .agg(aggregation)
)

print(
    f"Unique Klebsiella catalogue identifiers: "
    f"{len(catalogue_lookup)}"
)


# ----------------------------------------------------------------------
# 10. MERGE STRICT RESULTS WITH PA14 ANNOTATIONS
# ----------------------------------------------------------------------

# Rename catalogue measurements that already exist in the strict files.
rename_map = {}

for column in catalogue_lookup.columns:
    if (
        column != "Klebsiella_Locus_Tag"
        and column in strict_combined.columns
    ):
        rename_map[column] = f"{column}_Catalogue"

catalogue_lookup = catalogue_lookup.rename(
    columns=rename_map
)

annotated = strict_combined.merge(
    catalogue_lookup,
    how="left",
    on="Klebsiella_Locus_Tag",
    validate="many_to_one"
)


# ----------------------------------------------------------------------
# 11. RECONCILE CLUSTER COLUMN
# ----------------------------------------------------------------------

if "Final_Cluster_Catalogue" in annotated.columns:
    if "Final_Cluster" in annotated.columns:
        annotated["Final_Cluster"] = (
            annotated["Final_Cluster"]
            .combine_first(
                annotated["Final_Cluster_Catalogue"]
            )
        )
    else:
        annotated["Final_Cluster"] = (
            annotated["Final_Cluster_Catalogue"]
        )

    annotated = annotated.drop(
        columns=["Final_Cluster_Catalogue"]
    )


# ----------------------------------------------------------------------
# 12. CREATE PA14-CENTRED REPRESENTATIVE LABELS
# ----------------------------------------------------------------------

annotated["Representative_Gene_Label"] = annotated.apply(
    combine_gene_label,
    axis=1
)

annotated["Colistin_Conservation_Status"] = annotated.apply(
    colistin_status,
    axis=1
)


# ----------------------------------------------------------------------
# 13. REPORT MATCHING SUCCESS
# ----------------------------------------------------------------------

matched_to_catalogue = int(
    annotated["Gene_Label"].notna().sum()
)

matched_to_pa14 = int(
    annotated["Pseudomonas_Original_Locus_Tag"]
    .notna()
    .sum()
)

matched_to_gene_name = int(
    annotated["Pseudomonas_Gene_Name"]
    .notna()
    .sum()
)

print("\n" + "-" * 78)
print("ANNOTATION MATCHING RESULTS")
print("-" * 78)

print(
    f"Strict genes matched to catalogue: "
    f"{matched_to_catalogue} of {len(annotated)}"
)

print(
    f"Strict genes with Pseudomonas locus tags: "
    f"{matched_to_pa14} of {len(annotated)}"
)

print(
    f"Strict genes with Pseudomonas gene names: "
    f"{matched_to_gene_name} of {len(annotated)}"
)


# ----------------------------------------------------------------------
# 14. IDENTIFY UNMATCHED ORTHOLOGUES
# ----------------------------------------------------------------------

unmatched = annotated[
    annotated["Pseudomonas_Original_Locus_Tag"].isna()
].copy()

unmatched.to_csv(
    UNMATCHED_OUTPUT,
    index=False
)


# ----------------------------------------------------------------------
# 15. ORDER PA14-CENTRED OUTPUT COLUMNS
# ----------------------------------------------------------------------

preferred_columns = [
    "Representative_Gene_Label",
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "Gene_Label",
    "Final_Cluster",
    "Colistin_Conservation_Status",
    "Strict_Analysis_Group",
    "Colistin_RNAseq_Fold_Change",
    "Colistin_Riboseq_Fold_Change",
    "PA14_Colistin_RNAseq_Log2FC",
    "PA14_Colistin_Riboseq_Log2FC",
    "Pseudomonas_Colistin_Consensus_Log2FC",
    "K56_vs_Colistin_Log2FC",
    "Klebsiella_Colistin_Log2FC",
    "Cross_Species_Response_Class",
    "Response_Direction",
    "Cross_Species_Response_Strength",
    "Measured_Omics_Count",
]

preferred_columns = [
    column
    for column in preferred_columns
    if column in annotated.columns
]

remaining_columns = [
    column
    for column in annotated.columns
    if column not in preferred_columns
]

annotated = annotated[
    preferred_columns + remaining_columns
].copy()


# ----------------------------------------------------------------------
# 16. SPLIT CONSERVED AND DIVERGENT OUTPUTS
# ----------------------------------------------------------------------

annotated_conserved = annotated[
    annotated["Strict_Analysis_Group"] == "Conserved"
].copy()

annotated_divergent = annotated[
    annotated["Strict_Analysis_Group"] == "Divergent"
].copy()


# ----------------------------------------------------------------------
# 17. SAVE RESULTS
# ----------------------------------------------------------------------

annotated_conserved.to_csv(
    CONSERVE_OUTPUT,
    index=False
)

annotated_divergent.to_csv(
    DIVERGENT_OUTPUT,
    index=False
)

annotated.to_csv(
    COMBINED_OUTPUT,
    index=False
)

with pd.ExcelWriter(
    EXCEL_OUTPUT,
    engine="openpyxl"
) as writer:

    annotated_conserved.to_excel(
        writer,
        sheet_name="Strict_Conserved",
        index=False
    )

    annotated_divergent.to_excel(
        writer,
        sheet_name="Strict_Divergent",
        index=False
    )

    annotated.to_excel(
        writer,
        sheet_name="Combined",
        index=False
    )

    unmatched.to_excel(
        writer,
        sheet_name="Unmatched_PA14",
        index=False
    )


# ----------------------------------------------------------------------
# 18. FINAL SUMMARY
# ----------------------------------------------------------------------

print("\nSaved output files:")

print(f"  {CONSERVE_OUTPUT}")
print(f"  {DIVERGENT_OUTPUT}")
print(f"  {COMBINED_OUTPUT}")
print(f"  {UNMATCHED_OUTPUT}")
print(f"  {EXCEL_OUTPUT}")

if matched_to_pa14 < len(annotated):
    print(
        "\nNOTE: Some strict genes did not recover a Pseudomonas "
        "locus tag. This means those Klebsiella entries either lack "
        "a mapped Pseudomonas orthologue in the catalogue or use a "
        "different Klebsiella identifier format."
    )

print("\n" + "=" * 78)
print("STEP 28E COMPLETE")
print("=" * 78 + "\n")