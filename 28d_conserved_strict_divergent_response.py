import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


print("\n" + "=" * 78)
print("STEP 28D: ANNOTATE STRICT CONSERVED AND DIVERGENT RESPONSE GENES")
print("=" * 78 + "\n")


# ----------------------------------------------------------------
# 1. FILE SETTINGS
# ----------------------------------------------------------------

# Run this script from the main project folder.
PROJECT_DIRECTORY = Path.cwd()

# Step 28A strict output folder.
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

# Separate Step 28D output folder.
OUTPUT_DIRECTORY = (
    PROJECT_DIRECTORY /
    "Strict_Annotated_Response_Analysis_Outputs"
)

OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

ANNOTATED_CONSERVED_FILE = (
    OUTPUT_DIRECTORY /
    "Strict_Conserved_Response_Genes_Annotated.csv"
)

ANNOTATED_DIVERGENT_FILE = (
    OUTPUT_DIRECTORY /
    "Strict_Divergent_Response_Genes_Annotated.csv"
)

ANNOTATED_COMBINED_FILE = (
    OUTPUT_DIRECTORY /
    "Strict_Conserved_and_Divergent_Genes_Annotated.csv"
)

ANNOTATION_AUDIT_FILE = (
    OUTPUT_DIRECTORY /
    "Strict_Annotation_Source_Audit.csv"
)

UNMATCHED_FILE = (
    OUTPUT_DIRECTORY /
    "Strict_Genes_Without_Annotation_Matches.csv"
)

EXCEL_FILE = (
    OUTPUT_DIRECTORY /
    "Strict_Annotated_Response_Analysis.xlsx"
)


# ----------------------------------------------------------------
# 2. EXPECTED STRICT COUNTS
# ----------------------------------------------------------------

EXPECTED_CONSERVED_COUNT = 10
EXPECTED_DIVERGENT_COUNT = 8


# ----------------------------------------------------------------
# 3. CANDIDATE ANNOTATION FILES
# ----------------------------------------------------------------

# The script searches these first. Files that do not exist are skipped.
PREFERRED_ANNOTATION_FILES = [
    "Final_Cluster_Gene_Catalogue.csv",
    "Final_Biological_Cluster_Characterization_DETAILED.csv",
    "Final_Biological_Cluster_Characterization.csv",
    "Gene_Level_Cluster_Characterization.csv",
    "Cluster_Gene_Level_Characterization.csv",
    "Final_Clustered_Matrix.csv",
]

# It will then inspect other CSV files in the project folder if needed.
EXCLUDED_DIRECTORY_NAMES = {
    "Strict_Conserved_Response_Analysis_Outputs",
    "Strict_Annotated_Response_Analysis_Outputs",
}


# ----------------------------------------------------------------
# 4. COLUMN NAME CANDIDATES
# ----------------------------------------------------------------

COLUMN_CANDIDATES = {
    "Representative_Gene_Label": [
        "Representative_Gene_Label",
        "Representative_Gene",
        "Representative_Label",
        "Gene_Label",
    ],
    "PA14_Locus_Tag": [
        "PA14_Locus_Tag",
        "PA14_locus_tag",
        "Pseudomonas_Locus_Tag",
        "Pseudomonas_Gene_ID",
        "PA14_Gene_ID",
    ],
    "Gene_Name": [
        "Gene_Name",
        "gene_name",
        "Gene",
        "gene",
        "Symbol",
        "symbol",
        "PA14_Gene_Name",
    ],
    "Klebsiella_Locus_Tag": [
        "Klebsiella_Locus_Tag",
        "MGH78578_Locus_Tag",
        "K56_Locus_Tag",
        "KPHS_Locus_Tag",
        "Klebsiella_Gene_ID",
    ],
    "Product": [
        "Gene_Product",
        "gene_product",
        "Product",
        "product",
        "Description",
        "description",
        "Protein_Name",
        "protein_name",
        "Annotation",
        "annotation",
    ],
    "Final_Cluster": [
        "Final_Cluster",
        "Cluster",
        "Cluster_ID",
        "Hierarchical_Cluster",
    ],
    "Biological_Fingerprint": [
        "Biological_Fingerprint",
        "Biological_Cluster_Fingerprint",
        "Cluster_Biological_Fingerprint",
        "Final_Biological_Fingerprint",
    ],
    "TnSeq_Conservation_Status": [
        "TnSeq_Conservation_Status",
        "Tnseq_Conservation_Status",
        "TnSeq_Status",
        "Fitness_Conservation_Status",
    ],
    "Outlier_Status": [
        "Outlier_Status",
        "Cluster_Outlier_Status",
        "Is_Outlier",
    ],
}


# ----------------------------------------------------------------
# 5. HELPER FUNCTIONS
# ----------------------------------------------------------------

def normalise_column_name(name):
    """Normalise a column name for tolerant comparisons."""
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def find_column(dataframe, possible_names):
    """Find the first exact or normalised match for a column."""
    for name in possible_names:
        if name in dataframe.columns:
            return name

    normalised_lookup = {
        normalise_column_name(column): column
        for column in dataframe.columns
    }

    for name in possible_names:
        key = normalise_column_name(name)
        if key in normalised_lookup:
            return normalised_lookup[key]

    return None


def standardise_columns(dataframe):
    """
    Add standard annotation columns without deleting original columns.
    """
    result = dataframe.copy()

    for standard_name, candidates in COLUMN_CANDIDATES.items():
        source_column = find_column(result, candidates)

        if source_column is not None:
            if standard_name not in result.columns:
                result[standard_name] = result[source_column]

    return result


def clean_identifier(series):
    """Clean identifiers while preserving missing values."""
    cleaned = series.astype("string").str.strip()
    cleaned = cleaned.replace(
        {
            "": pd.NA,
            "nan": pd.NA,
            "None": pd.NA,
            "<NA>": pd.NA,
        }
    )
    return cleaned


def first_non_missing(series):
    """Return the first non-missing value in a grouped series."""
    non_missing = series.dropna()

    if len(non_missing) == 0:
        return pd.NA

    return non_missing.iloc[0]


def read_csv_safely(path):
    """Read a CSV with tolerant encoding handling."""
    try:
        return pd.read_csv(path, low_memory=False)
    except UnicodeDecodeError:
        return pd.read_csv(
            path,
            encoding="latin-1",
            low_memory=False
        )
    except Exception as error:
        warnings.warn(f"Could not read {path}: {error}")
        return None


def discover_csv_files():
    """Find likely annotation CSVs in the project directory."""
    discovered = []
    seen = set()

    for filename in PREFERRED_ANNOTATION_FILES:
        candidate = PROJECT_DIRECTORY / filename

        if candidate.exists():
            resolved = candidate.resolve()
            if resolved not in seen:
                discovered.append(candidate)
                seen.add(resolved)

    for candidate in PROJECT_DIRECTORY.rglob("*.csv"):
        if any(
            directory_name in candidate.parts
            for directory_name in EXCLUDED_DIRECTORY_NAMES
        ):
            continue

        resolved = candidate.resolve()

        if resolved in seen:
            continue

        discovered.append(candidate)
        seen.add(resolved)

    return discovered


def score_annotation_table(dataframe):
    """
    Score a table based on how useful it is for gene-level annotation.
    """
    standardised = standardise_columns(dataframe)

    score = 0

    weighted_columns = {
        "Klebsiella_Locus_Tag": 10,
        "PA14_Locus_Tag": 8,
        "Representative_Gene_Label": 8,
        "Gene_Name": 6,
        "Product": 4,
        "Final_Cluster": 3,
        "Biological_Fingerprint": 5,
        "TnSeq_Conservation_Status": 4,
        "Outlier_Status": 3,
    }

    for column, weight in weighted_columns.items():
        if column in standardised.columns:
            score += weight

    return score


def collapse_annotation_table(dataframe, key_column):
    """
    Collapse duplicate annotations to one row per identifier.
    """
    working = dataframe.copy()

    if key_column not in working.columns:
        return None

    working[key_column] = clean_identifier(working[key_column])
    working = working[working[key_column].notna()].copy()

    if len(working) == 0:
        return None

    useful_columns = [
        column for column in [
            key_column,
            "Representative_Gene_Label",
            "PA14_Locus_Tag",
            "Gene_Name",
            "Klebsiella_Locus_Tag",
            "Product",
            "Final_Cluster",
            "Biological_Fingerprint",
            "TnSeq_Conservation_Status",
            "Outlier_Status",
        ]
        if column in working.columns
    ]

    working = working[useful_columns]

    aggregation = {
        column: first_non_missing
        for column in useful_columns
        if column != key_column
    }

    collapsed = (
        working
        .groupby(key_column, as_index=False, dropna=False)
        .agg(aggregation)
    )

    return collapsed


def fill_from_lookup(
    target,
    lookup,
    target_key,
    lookup_key,
    columns_to_fill,
    source_label
):
    """
    Fill missing annotation columns in target using a lookup table.
    """
    if (
        lookup is None
        or target_key not in target.columns
        or lookup_key not in lookup.columns
    ):
        return target, 0

    left = target.copy()
    right = lookup.copy()

    left[target_key] = clean_identifier(left[target_key])
    right[lookup_key] = clean_identifier(right[lookup_key])

    right_columns = [
        lookup_key
    ] + [
        column for column in columns_to_fill
        if column in right.columns
    ]

    right = right[right_columns].drop_duplicates(
        subset=[lookup_key],
        keep="first"
    )

    rename_map = {
        column: f"{column}__lookup"
        for column in right_columns
        if column != lookup_key
    }

    right = right.rename(columns=rename_map)

    merged = left.merge(
        right,
        how="left",
        left_on=target_key,
        right_on=lookup_key,
        suffixes=("", "__key")
    )

    matched_mask = pd.Series(False, index=merged.index)

    for column in columns_to_fill:
        lookup_column = f"{column}__lookup"

        if lookup_column not in merged.columns:
            continue

        if column not in merged.columns:
            merged[column] = pd.NA

        before_missing = merged[column].isna()
        merged[column] = merged[column].combine_first(
            merged[lookup_column]
        )

        newly_filled = before_missing & merged[column].notna()
        matched_mask = matched_mask | newly_filled

        merged = merged.drop(columns=[lookup_column])

    extra_key_columns = [
        column for column in merged.columns
        if column.endswith("__key")
    ]

    if lookup_key != target_key and lookup_key in merged.columns:
        merged = merged.drop(columns=[lookup_key])

    if extra_key_columns:
        merged = merged.drop(columns=extra_key_columns)

    if "Annotation_Source" not in merged.columns:
        merged["Annotation_Source"] = pd.NA

    source_needed = matched_mask & merged["Annotation_Source"].isna()
    merged.loc[source_needed, "Annotation_Source"] = source_label

    return merged, int(matched_mask.sum())


def make_representative_label(row):
    """
    Build a readable PA14 | gene name | Klebsiella label.
    """
    components = []

    for column in [
        "PA14_Locus_Tag",
        "Gene_Name",
        "Klebsiella_Locus_Tag",
    ]:
        value = row.get(column, pd.NA)

        if pd.notna(value):
            value = str(value).strip()

            if value and value.lower() not in {"nan", "none"}:
                if value not in components:
                    components.append(value)

    if components:
        return " | ".join(components)

    existing = row.get("Representative_Gene_Label", pd.NA)

    if pd.notna(existing):
        return str(existing)

    return "Unlabelled gene"


def classify_colistin_status(row):
    """
    Convert the strict Step 28A class into a readable status.
    """
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


def infer_outlier_status(dataframe):
    """
    Use supplied outlier labels where available. Otherwise infer an
    outlier cluster as one containing three or fewer genes in the
    complete 6,970-gene classification table.
    """
    result = dataframe.copy()

    if "Outlier_Status" not in result.columns:
        result["Outlier_Status"] = pd.NA

    normalised_existing = (
        result["Outlier_Status"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    valid_existing = normalised_existing.isin(
        {"OUTLIER", "NOT OUTLIER"}
    )

    result.loc[
        valid_existing,
        "Outlier_Status"
    ] = normalised_existing[valid_existing]

    # Prefer the complete strict classification file for cluster sizes.
    full_classification_file = (
        STRICT_DIRECTORY /
        "Strict_Cross_Species_Response_Classification.csv"
    )

    cluster_sizes = None

    if full_classification_file.exists():
        full_df = read_csv_safely(full_classification_file)

        if (
            full_df is not None
            and "Final_Cluster" in full_df.columns
        ):
            cluster_sizes = (
                full_df["Final_Cluster"]
                .value_counts(dropna=False)
                .to_dict()
            )

    if cluster_sizes is None and "Final_Cluster" in result.columns:
        cluster_sizes = (
            result["Final_Cluster"]
            .value_counts(dropna=False)
            .to_dict()
        )

    if cluster_sizes is not None:
        missing_status = result["Outlier_Status"].isna()

        result.loc[
            missing_status,
            "Outlier_Status"
        ] = result.loc[
            missing_status,
            "Final_Cluster"
        ].map(
            lambda cluster: (
                "OUTLIER"
                if cluster_sizes.get(cluster, np.inf) <= 3
                else "NOT OUTLIER"
            )
        )

    result["Outlier_Status"] = (
        result["Outlier_Status"]
        .fillna("NOT AVAILABLE")
    )

    return result


# ----------------------------------------------------------------
# 6. LOAD STRICT STEP 28A RESULTS
# ----------------------------------------------------------------

missing_strict_files = [
    str(path)
    for path in [
        STRICT_CONSERVED_FILE,
        STRICT_DIVERGENT_FILE,
    ]
    if not path.exists()
]

if missing_strict_files:
    raise FileNotFoundError(
        "\nThe following Step 28A strict files were not found:\n  "
        + "\n  ".join(missing_strict_files)
        + "\n\nRun Step 28A first, or place its output folder in "
          "the project directory."
    )

strict_conserved = read_csv_safely(
    STRICT_CONSERVED_FILE
)

strict_divergent = read_csv_safely(
    STRICT_DIVERGENT_FILE
)

if strict_conserved is None or strict_divergent is None:
    raise RuntimeError(
        "The strict Step 28A files could not be read."
    )

strict_conserved["Strict_Analysis_Group"] = "Conserved"
strict_divergent["Strict_Analysis_Group"] = "Divergent"

strict_combined = pd.concat(
    [strict_conserved, strict_divergent],
    ignore_index=True,
    sort=False
)

strict_combined = standardise_columns(strict_combined)

if "Klebsiella_Locus_Tag" not in strict_combined.columns:
    raise KeyError(
        "\nThe strict files do not contain a recognised "
        "Klebsiella locus-tag column.\n"
        f"Available columns:\n{list(strict_combined.columns)}"
    )

strict_combined["Klebsiella_Locus_Tag"] = clean_identifier(
    strict_combined["Klebsiella_Locus_Tag"]
)

print(
    f"Loaded strict conserved genes: {len(strict_conserved)}"
)
print(
    f"Loaded strict divergent genes: {len(strict_divergent)}"
)


# ----------------------------------------------------------------
# 7. VERIFY STRICT COUNTS
# ----------------------------------------------------------------

if len(strict_conserved) != EXPECTED_CONSERVED_COUNT:
    warnings.warn(
        f"Expected {EXPECTED_CONSERVED_COUNT} strict conserved genes, "
        f"but found {len(strict_conserved)}."
    )

if len(strict_divergent) != EXPECTED_DIVERGENT_COUNT:
    warnings.warn(
        f"Expected {EXPECTED_DIVERGENT_COUNT} strict divergent genes, "
        f"but found {len(strict_divergent)}."
    )


# ----------------------------------------------------------------
# 8. DISCOVER AND SCORE ANNOTATION TABLES
# ----------------------------------------------------------------

candidate_files = discover_csv_files()
annotation_tables = []
audit_rows = []

for path in candidate_files:
    table = read_csv_safely(path)

    if table is None:
        continue

    standardised = standardise_columns(table)
    score = score_annotation_table(table)

    useful_columns = [
        column for column in [
            "Representative_Gene_Label",
            "PA14_Locus_Tag",
            "Gene_Name",
            "Klebsiella_Locus_Tag",
            "Product",
            "Final_Cluster",
            "Biological_Fingerprint",
            "TnSeq_Conservation_Status",
            "Outlier_Status",
        ]
        if column in standardised.columns
    ]

    audit_rows.append(
        {
            "Source_File": str(
                path.relative_to(PROJECT_DIRECTORY)
            ),
            "Rows": len(standardised),
            "Columns": len(standardised.columns),
            "Annotation_Score": score,
            "Useful_Columns": " | ".join(useful_columns),
        }
    )

    if score > 0:
        annotation_tables.append(
            {
                "path": path,
                "data": standardised,
                "score": score,
            }
        )

annotation_tables = sorted(
    annotation_tables,
    key=lambda item: item["score"],
    reverse=True
)

audit_df = pd.DataFrame(audit_rows)

if len(audit_df) > 0:
    audit_df = audit_df.sort_values(
        "Annotation_Score",
        ascending=False
    )

audit_df.to_csv(
    ANNOTATION_AUDIT_FILE,
    index=False
)

print(
    f"\nCandidate annotation tables inspected: "
    f"{len(candidate_files)}"
)
print(
    f"Tables containing useful annotation columns: "
    f"{len(annotation_tables)}"
)


# ----------------------------------------------------------------
# 9. MERGE GENE-LEVEL ANNOTATIONS
# ----------------------------------------------------------------

annotated = strict_combined.copy()

annotation_columns_to_fill = [
    "Representative_Gene_Label",
    "PA14_Locus_Tag",
    "Gene_Name",
    "Klebsiella_Locus_Tag",
    "Product",
    "Biological_Fingerprint",
    "TnSeq_Conservation_Status",
    "Outlier_Status",
]

merge_audit = []

for item in annotation_tables:
    path = item["path"]
    table = item["data"]

    # First choice: exact Klebsiella locus-tag match.
    lookup = collapse_annotation_table(
        table,
        "Klebsiella_Locus_Tag"
    )

    annotated, matched = fill_from_lookup(
        target=annotated,
        lookup=lookup,
        target_key="Klebsiella_Locus_Tag",
        lookup_key="Klebsiella_Locus_Tag",
        columns_to_fill=annotation_columns_to_fill,
        source_label=str(path.relative_to(PROJECT_DIRECTORY)),
    )

    merge_audit.append(
        {
            "Source_File": str(
                path.relative_to(PROJECT_DIRECTORY)
            ),
            "Merge_Key": "Klebsiella_Locus_Tag",
            "Rows_With_New_Annotations": matched,
        }
    )


# ----------------------------------------------------------------
# 10. MERGE CLUSTER-LEVEL BIOLOGICAL FINGERPRINTS
# ----------------------------------------------------------------

cluster_columns_to_fill = [
    "Biological_Fingerprint",
    "Outlier_Status",
]

for item in annotation_tables:
    path = item["path"]
    table = item["data"]

    if "Final_Cluster" not in table.columns:
        continue

    cluster_lookup = collapse_annotation_table(
        table,
        "Final_Cluster"
    )

    annotated, matched = fill_from_lookup(
        target=annotated,
        lookup=cluster_lookup,
        target_key="Final_Cluster",
        lookup_key="Final_Cluster",
        columns_to_fill=cluster_columns_to_fill,
        source_label=str(path.relative_to(PROJECT_DIRECTORY)),
    )

    merge_audit.append(
        {
            "Source_File": str(
                path.relative_to(PROJECT_DIRECTORY)
            ),
            "Merge_Key": "Final_Cluster",
            "Rows_With_New_Annotations": matched,
        }
    )


# ----------------------------------------------------------------
# 11. BUILD PUBLICATION-READY COLUMNS
# ----------------------------------------------------------------

annotated["Representative_Gene_Label"] = annotated.apply(
    make_representative_label,
    axis=1
)

annotated["Colistin_Conservation_Status"] = annotated.apply(
    classify_colistin_status,
    axis=1
)

if "TnSeq_Conservation_Status" not in annotated.columns:
    annotated["TnSeq_Conservation_Status"] = (
        "Not available in source tables"
    )
else:
    annotated["TnSeq_Conservation_Status"] = (
        annotated["TnSeq_Conservation_Status"]
        .fillna("Not available in source tables")
    )

if "Biological_Fingerprint" not in annotated.columns:
    annotated["Biological_Fingerprint"] = (
        "Not available in source tables"
    )
else:
    annotated["Biological_Fingerprint"] = (
        annotated["Biological_Fingerprint"]
        .fillna("Not available in source tables")
    )

annotated = infer_outlier_status(annotated)

if "Product" not in annotated.columns:
    annotated["Product"] = pd.NA


# ----------------------------------------------------------------
# 12. SELECT AND ORDER OUTPUT COLUMNS
# ----------------------------------------------------------------

publication_columns = [
    "Representative_Gene_Label",
    "PA14_Locus_Tag",
    "Gene_Name",
    "Klebsiella_Locus_Tag",
    "Product",
    "Final_Cluster",
    "Biological_Fingerprint",
    "Colistin_Conservation_Status",
    "TnSeq_Conservation_Status",
    "Outlier_Status",
    "Strict_Analysis_Group",
    "Colistin_RNAseq_Fold_Change",
    "Colistin_Riboseq_Fold_Change",
    "K56_vs_Colistin_Log2FC",
    "PA14_Colistin_RNAseq_Log2FC",
    "PA14_Colistin_Riboseq_Log2FC",
    "Pseudomonas_Colistin_Consensus_Log2FC",
    "Klebsiella_Colistin_Log2FC",
    "Cross_Species_Response_Class",
    "Response_Direction",
    "Cross_Species_Response_Strength",
    "Annotation_Source",
]

publication_columns = [
    column
    for column in publication_columns
    if column in annotated.columns
]

annotated = annotated[publication_columns].copy()

annotated_conserved = annotated[
    annotated["Strict_Analysis_Group"] == "Conserved"
].copy()

annotated_divergent = annotated[
    annotated["Strict_Analysis_Group"] == "Divergent"
].copy()


# ----------------------------------------------------------------
# 13. IDENTIFY UNMATCHED GENES
# ----------------------------------------------------------------

unmatched_mask = (
    annotated["PA14_Locus_Tag"].isna()
    if "PA14_Locus_Tag" in annotated.columns
    else pd.Series(True, index=annotated.index)
)

if "Gene_Name" in annotated.columns:
    unmatched_mask = (
        unmatched_mask
        &
        annotated["Gene_Name"].isna()
    )

unmatched = annotated[unmatched_mask].copy()

unmatched.to_csv(
    UNMATCHED_FILE,
    index=False
)


# ----------------------------------------------------------------
# 14. SAVE ANNOTATED OUTPUTS
# ----------------------------------------------------------------

annotated_conserved.to_csv(
    ANNOTATED_CONSERVED_FILE,
    index=False
)

annotated_divergent.to_csv(
    ANNOTATED_DIVERGENT_FILE,
    index=False
)

annotated.to_csv(
    ANNOTATED_COMBINED_FILE,
    index=False
)

merge_audit_df = pd.DataFrame(merge_audit)

with pd.ExcelWriter(
    EXCEL_FILE,
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
        sheet_name="Combined_Annotated",
        index=False
    )

    unmatched.to_excel(
        writer,
        sheet_name="Unmatched_Genes",
        index=False
    )

    audit_df.to_excel(
        writer,
        sheet_name="Source_Audit",
        index=False
    )

    merge_audit_df.to_excel(
        writer,
        sheet_name="Merge_Audit",
        index=False
    )


# ----------------------------------------------------------------
# 15. FINAL REPORT
# ----------------------------------------------------------------

pa14_annotated_count = (
    int(annotated["PA14_Locus_Tag"].notna().sum())
    if "PA14_Locus_Tag" in annotated.columns
    else 0
)

gene_name_count = (
    int(annotated["Gene_Name"].notna().sum())
    if "Gene_Name" in annotated.columns
    else 0
)

fingerprint_count = int(
    (
        annotated["Biological_Fingerprint"]
        != "Not available in source tables"
    ).sum()
)

tnseq_count = int(
    (
        annotated["TnSeq_Conservation_Status"]
        != "Not available in source tables"
    ).sum()
)

print("\n" + "-" * 78)
print("ANNOTATION RESULTS")
print("-" * 78)

print(f"Strict conserved genes exported: {len(annotated_conserved)}")
print(f"Strict divergent genes exported: {len(annotated_divergent)}")
print(f"Genes with PA14 locus tags: {pa14_annotated_count}")
print(f"Genes with gene names: {gene_name_count}")
print(f"Genes with biological fingerprints: {fingerprint_count}")
print(f"Genes with Tn-seq conservation labels: {tnseq_count}")
print(f"Genes lacking PA14/gene-name matches: {len(unmatched)}")

print("\nSaved annotated files:")
print(f"  {ANNOTATED_CONSERVED_FILE}")
print(f"  {ANNOTATED_DIVERGENT_FILE}")
print(f"  {ANNOTATED_COMBINED_FILE}")
print(f"  {UNMATCHED_FILE}")
print(f"  {ANNOTATION_AUDIT_FILE}")
print(f"  {EXCEL_FILE}")

if pa14_annotated_count == 0:
    print(
        "\nIMPORTANT: No PA14 locus tags were recovered. "
        "Place Final_Cluster_Gene_Catalogue.csv or another "
        "gene-level ortholog annotation table in the project folder "
        "and rerun this script."
    )

print("\n" + "=" * 78)
print("STEP 28D COMPLETE")
print("=" * 78 + "\n")