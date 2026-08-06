# ================================================================
# STEP 28A: STRICT CONSERVED AND DIVERGENT RESPONSE ANALYSIS
# Figure 3.7 outputs
# ================================================================

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


print("\n" + "=" * 72)
print("STEP 28A: STRICT CONSERVED AND DIVERGENT RESPONSE ANALYSIS")
print("=" * 72 + "\n")


# ----------------------------------------------------------------
# 1. FILE SETTINGS
# ----------------------------------------------------------------

INPUT_FILE = "Final_Clustered_Matrix.csv"

OUTPUT_DIRECTORY = "Strict_Conserved_Response_Analysis_Outputs"

CONSERVE_FILE = os.path.join(
    OUTPUT_DIRECTORY,
    "Strict_Conserved_Response_Genes.csv"
)

DIVERGENT_FILE = os.path.join(
    OUTPUT_DIRECTORY,
    "Strict_Divergent_Response_Genes.csv"
)

ALL_RESULTS_FILE = os.path.join(
    OUTPUT_DIRECTORY,
    "Strict_Cross_Species_Response_Classification.csv"
)

SUMMARY_FILE = os.path.join(
    OUTPUT_DIRECTORY,
    "Strict_Conserved_Divergent_Response_Summary.csv"
)

HEATMAP_DATA_FILE = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure_3_7_Heatmap_Data.csv"
)

SCATTER_FIGURE = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure_3_7A_Cross_Species_Response_Scatterplot.png"
)

HEATMAP_FIGURE = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure_3_7B_Conserved_Divergent_Heatmap.png"
)

COMBINED_FIGURE = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure_3_7_Conserved_and_Divergent_Response_Genes.png"
)


# ----------------------------------------------------------------
# 2. ANALYSIS SETTINGS
# ----------------------------------------------------------------

# Strict response threshold:
# Both organisms must show an absolute response of at least 1 log2FC.
LOG2FC_THRESHOLD = 1.0

# Minimum number of Pseudomonas colistin measurements required.
# Setting this to 2 means that both RNA-seq and Ribo-seq must be available.
MIN_PSEUDOMONAS_MEASUREMENTS = 2

# Expected values from the original analysis.
EXPECTED_CONSERVED_COUNT = 7
EXPECTED_DIVERGENT_COUNT = 15

# Figure resolution.
FIGURE_DPI = 600


# ----------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ----------------------------------------------------------------

def find_column(dataframe, possible_names, required=True):
    """
    Find the first matching column from a list of possible names.
    This accommodates slight naming differences between script versions.
    """

    for name in possible_names:
        if name in dataframe.columns:
            return name

    if required:
        raise KeyError(
            "\nNone of the expected columns were found:\n"
            f"{possible_names}\n\n"
            f"Available columns are:\n{list(dataframe.columns)}"
        )

    return None


def safely_convert_numeric(series):
    """
    Convert a pandas Series to numeric values.
    Non-numeric values are converted to NaN.
    """

    return pd.to_numeric(series, errors="coerce")


def fold_change_to_log2(series):
    """
    Convert an ordinary fold-change column to signed log2 fold change.

    Positive fold changes remain positive.
    Negative values are interpreted as signed fold changes when present.

    Values equal to zero become missing because log2(0) is undefined.
    """

    numeric_series = safely_convert_numeric(series)

    output = pd.Series(
        np.nan,
        index=numeric_series.index,
        dtype=float
    )

    positive_mask = numeric_series > 0
    negative_mask = numeric_series < 0

    output.loc[positive_mask] = np.log2(
        numeric_series.loc[positive_mask]
    )

    # This section retains the sign when a dataset already contains
    # signed fold-change values rather than ratios.
    output.loc[negative_mask] = -np.log2(
        np.abs(numeric_series.loc[negative_mask]) + 1
    )

    return output


def prepare_response_column(dataframe, column_name):
    """
    Return a log2-scale response column.

    Columns containing 'log2' are treated as already log2-transformed.
    Other response columns are assessed before conversion.
    """

    values = safely_convert_numeric(dataframe[column_name])

    lower_name = column_name.lower()

    if "log2" in lower_name:
        print(f"Using as log2FC: {column_name}")
        return values

    # The PA14 columns in the current matrix contain signed values,
    # including values below zero. They therefore behave like
    # log-scale response values despite some names containing
    # 'Fold_Change' or 'FC'.
    if values.min(skipna=True) < 0:
        print(
            f"Using as signed response/log2-scale values: {column_name}"
        )
        return values

    print(f"Converting ordinary fold change to log2FC: {column_name}")
    return fold_change_to_log2(values)


def choose_gene_label_columns(dataframe):
    """
    Identify useful gene identifier and annotation columns.
    """

    identifier_candidates = [
        "PA14_Locus_Tag",
        "PA14_locus_tag",
        "PA14_Gene",
        "PA14_Gene_Name",
        "Gene",
        "gene",
        "Gene_Name",
        "gene_name",
        "Locus_Tag",
        "locus_tag",
        "PAO1_Locus_Tag",
        "Klebsiella_Locus_Tag",
        "MGH78578_Locus_Tag"
    ]

    annotation_candidates = [
        "Gene_Product",
        "gene_product",
        "Product",
        "product",
        "Description",
        "description",
        "Protein_Name",
        "protein_name",
        "Annotation",
        "annotation"
    ]

    identifier_column = find_column(
        dataframe,
        identifier_candidates,
        required=False
    )

    annotation_column = find_column(
        dataframe,
        annotation_candidates,
        required=False
    )

    if identifier_column is None:
        dataframe["Analysis_Gene_ID"] = [
            f"Gene_{i + 1}" for i in range(len(dataframe))
        ]
        identifier_column = "Analysis_Gene_ID"

    return identifier_column, annotation_column


def classify_response(row):
    """
    Classify each gene using strict cross-species criteria.

    Conserved:
        Pseudomonas and Klebsiella both respond above threshold
        and have the same response direction.

    Divergent:
        Pseudomonas and Klebsiella both respond above threshold
        but have opposite response directions.

    Below threshold:
        Comparable measurements exist, but one or both organisms
        do not reach the strict response threshold.

    Not comparable:
        Required response measurements are missing.
    """

    pseudomonas = row["Pseudomonas_Colistin_Consensus_Log2FC"]
    klebsiella = row["Klebsiella_Colistin_Log2FC"]

    if pd.isna(pseudomonas) or pd.isna(klebsiella):
        return "Not comparable"

    p_responsive = abs(pseudomonas) >= LOG2FC_THRESHOLD
    k_responsive = abs(klebsiella) >= LOG2FC_THRESHOLD

    if not (p_responsive and k_responsive):
        return "Below strict threshold"

    if np.sign(pseudomonas) == np.sign(klebsiella):
        return "Conserved"

    return "Divergent"


# ----------------------------------------------------------------
# 4. LOAD INPUT DATA
# ----------------------------------------------------------------

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(
        f"\nInput file not found: {INPUT_FILE}\n"
        "Place the script in the same folder as "
        "Final_Clustered_Matrix.csv."
    )

os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")


# ----------------------------------------------------------------
# 5. DETECT REQUIRED COLUMNS
# ----------------------------------------------------------------

pseudomonas_rnaseq_column = find_column(
    df,
    [
        "Colistin_RNAseq_Fold_Change",
        "Colistin_RNAseq_FC",
        "Colistin_RNAseq_Log2FC",
        "PA14_Colistin_RNAseq_Fold_Change",
        "PA14_Colistin_RNAseq_Log2FC"
    ]
)

pseudomonas_riboseq_column = find_column(
    df,
    [
        "Colistin_Riboseq_Fold_Change",
        "Colistin_Riboseq_FC",
        "Colistin_Riboseq_Log2FC",
        "PA14_Colistin_Riboseq_Fold_Change",
        "PA14_Colistin_Riboseq_Log2FC"
    ]
)

klebsiella_colistin_column = find_column(
    df,
    [
        "K56_vs_Colistin_Log2FC",
        "K56_Colistin_Log2FC",
        "Klebsiella_Colistin_Log2FC",
        "Colistin_K56_Log2FC"
    ]
)

gene_id_column, annotation_column = choose_gene_label_columns(df)

print("\nDetected columns:")
print(f"  Gene identifier: {gene_id_column}")
print(f"  Annotation: {annotation_column}")
print(f"  PA14 colistin RNA-seq: {pseudomonas_rnaseq_column}")
print(f"  PA14 colistin Ribo-seq: {pseudomonas_riboseq_column}")
print(f"  K56 colistin response: {klebsiella_colistin_column}")


# ----------------------------------------------------------------
# 6. PREPARE RESPONSE VALUES
# ----------------------------------------------------------------

df["PA14_Colistin_RNAseq_Log2FC"] = prepare_response_column(
    df,
    pseudomonas_rnaseq_column
)

df["PA14_Colistin_Riboseq_Log2FC"] = prepare_response_column(
    df,
    pseudomonas_riboseq_column
)

df["Klebsiella_Colistin_Log2FC"] = prepare_response_column(
    df,
    klebsiella_colistin_column
)

pseudomonas_columns = [
    "PA14_Colistin_RNAseq_Log2FC",
    "PA14_Colistin_Riboseq_Log2FC"
]

df["Pseudomonas_Measurement_Count"] = (
    df[pseudomonas_columns]
    .notna()
    .sum(axis=1)
)

# Require both PA14 RNA-seq and Ribo-seq measurements.
valid_pseudomonas_mask = (
    df["Pseudomonas_Measurement_Count"]
    >= MIN_PSEUDOMONAS_MEASUREMENTS
)

df["Pseudomonas_Colistin_Consensus_Log2FC"] = np.nan

# Median is used as the consensus because it is less affected by
# unusually large values than the arithmetic mean.
df.loc[
    valid_pseudomonas_mask,
    "Pseudomonas_Colistin_Consensus_Log2FC"
] = (
    df.loc[valid_pseudomonas_mask, pseudomonas_columns]
    .median(axis=1)
)


# ----------------------------------------------------------------
# 7. REQUIRE CONSISTENCY BETWEEN RNA-SEQ AND RIBO-SEQ
# ----------------------------------------------------------------

df["PA14_RNA_Ribo_Same_Direction"] = (
    np.sign(df["PA14_Colistin_RNAseq_Log2FC"])
    ==
    np.sign(df["PA14_Colistin_Riboseq_Log2FC"])
)

# When one PA14 layer increases and the other decreases, the gene is
# not considered a clear Pseudomonas consensus response.
inconsistent_mask = (
    valid_pseudomonas_mask
    &
    ~df["PA14_RNA_Ribo_Same_Direction"]
)

df.loc[
    inconsistent_mask,
    "Pseudomonas_Colistin_Consensus_Log2FC"
] = np.nan

print(
    "\nGenes with inconsistent PA14 RNA-seq/Ribo-seq directions:",
    int(inconsistent_mask.sum())
)


# ----------------------------------------------------------------
# 8. CLASSIFY CONSERVED AND DIVERGENT RESPONSES
# ----------------------------------------------------------------

df["Cross_Species_Response_Class"] = df.apply(
    classify_response,
    axis=1
)

df["Response_Direction"] = np.where(
    df["Pseudomonas_Colistin_Consensus_Log2FC"] > 0,
    "Increased",
    np.where(
        df["Pseudomonas_Colistin_Consensus_Log2FC"] < 0,
        "Decreased",
        "Neutral"
    )
)

df["Cross_Species_Response_Strength"] = (
    np.abs(df["Pseudomonas_Colistin_Consensus_Log2FC"])
    +
    np.abs(df["Klebsiella_Colistin_Log2FC"])
)

conserved_df = df[
    df["Cross_Species_Response_Class"] == "Conserved"
].copy()

divergent_df = df[
    df["Cross_Species_Response_Class"] == "Divergent"
].copy()

conserved_df = conserved_df.sort_values(
    "Cross_Species_Response_Strength",
    ascending=False
)

divergent_df = divergent_df.sort_values(
    "Cross_Species_Response_Strength",
    ascending=False
)

conserved_count = len(conserved_df)
divergent_count = len(divergent_df)

print("\nStrict classification results:")
print(f"  Conserved genes: {conserved_count}")
print(f"  Divergent genes: {divergent_count}")


# ----------------------------------------------------------------
# 9. SELECT OUTPUT COLUMNS
# ----------------------------------------------------------------

output_columns = []

for column in [
    gene_id_column,
    annotation_column,
    "Final_Cluster",
    pseudomonas_rnaseq_column,
    pseudomonas_riboseq_column,
    klebsiella_colistin_column,
    "PA14_Colistin_RNAseq_Log2FC",
    "PA14_Colistin_Riboseq_Log2FC",
    "Pseudomonas_Colistin_Consensus_Log2FC",
    "Klebsiella_Colistin_Log2FC",
    "PA14_RNA_Ribo_Same_Direction",
    "Cross_Species_Response_Class",
    "Response_Direction",
    "Cross_Species_Response_Strength"
]:
    if column is not None and column in df.columns:
        if column not in output_columns:
            output_columns.append(column)


# ----------------------------------------------------------------
# 10. SAVE TABLES
# ----------------------------------------------------------------

conserved_df[output_columns].to_csv(
    CONSERVE_FILE,
    index=False
)

divergent_df[output_columns].to_csv(
    DIVERGENT_FILE,
    index=False
)

df[output_columns].to_csv(
    ALL_RESULTS_FILE,
    index=False
)

summary_df = pd.DataFrame(
    {
        "Response_Category": [
            "Conserved",
            "Divergent",
            "Below strict threshold",
            "Not comparable"
        ],
        "Number_of_Genes": [
            conserved_count,
            divergent_count,
            int(
                (
                    df["Cross_Species_Response_Class"]
                    == "Below strict threshold"
                ).sum()
            ),
            int(
                (
                    df["Cross_Species_Response_Class"]
                    == "Not comparable"
                ).sum()
            )
        ]
    }
)

summary_df.to_csv(
    SUMMARY_FILE,
    index=False
)

print("\nSaved output tables:")
print(f"  {CONSERVE_FILE}")
print(f"  {DIVERGENT_FILE}")
print(f"  {ALL_RESULTS_FILE}")
print(f"  {SUMMARY_FILE}")


# ----------------------------------------------------------------
# 11. PREPARE FIGURE DATA
# ----------------------------------------------------------------

figure_gene_df = pd.concat(
    [conserved_df, divergent_df],
    axis=0,
    ignore_index=True
)

figure_gene_df["Figure_Group"] = (
    figure_gene_df["Cross_Species_Response_Class"]
)

# Build readable labels.
figure_gene_df["Figure_Gene_Label"] = (
    figure_gene_df[gene_id_column]
    .fillna("Unknown gene")
    .astype(str)
)

# Ensure duplicate labels remain distinguishable.
duplicated_labels = figure_gene_df[
    "Figure_Gene_Label"
].duplicated(keep=False)

figure_gene_df.loc[
    duplicated_labels,
    "Figure_Gene_Label"
] = (
    figure_gene_df.loc[
        duplicated_labels,
        "Figure_Gene_Label"
    ]
    +
    "_"
    +
    figure_gene_df.loc[
        duplicated_labels
    ].index.astype(str)
)

heatmap_columns = [
    "PA14_Colistin_RNAseq_Log2FC",
    "PA14_Colistin_Riboseq_Log2FC",
    "Klebsiella_Colistin_Log2FC"
]

heatmap_labels = [
    "PA14 RNA-seq",
    "PA14 Ribo-seq",
    "K56 RNA-seq"
]

heatmap_output_columns = [
    gene_id_column,
    "Cross_Species_Response_Class"
] + heatmap_columns

figure_gene_df[heatmap_output_columns].to_csv(
    HEATMAP_DATA_FILE,
    index=False
)


# ----------------------------------------------------------------
# 12. PANEL A: CROSS-SPECIES SCATTERPLOT
# ----------------------------------------------------------------

comparable_df = df[
    df["Pseudomonas_Colistin_Consensus_Log2FC"].notna()
    &
    df["Klebsiella_Colistin_Log2FC"].notna()
].copy()

fig, ax = plt.subplots(figsize=(8.2, 7.2))

below_threshold = comparable_df[
    comparable_df["Cross_Species_Response_Class"]
    == "Below strict threshold"
]

ax.scatter(
    below_threshold["Pseudomonas_Colistin_Consensus_Log2FC"],
    below_threshold["Klebsiella_Colistin_Log2FC"],
    s=16,
    alpha=0.28,
    label="Other comparable genes"
)

ax.scatter(
    conserved_df["Pseudomonas_Colistin_Consensus_Log2FC"],
    conserved_df["Klebsiella_Colistin_Log2FC"],
    s=55,
    marker="o",
    edgecolors="black",
    linewidths=0.5,
    label=f"Conserved (n = {conserved_count})"
)

ax.scatter(
    divergent_df["Pseudomonas_Colistin_Consensus_Log2FC"],
    divergent_df["Klebsiella_Colistin_Log2FC"],
    s=65,
    marker="^",
    edgecolors="black",
    linewidths=0.5,
    label=f"Divergent (n = {divergent_count})"
)

ax.axhline(
    0,
    linewidth=0.8,
    linestyle="--"
)

ax.axvline(
    0,
    linewidth=0.8,
    linestyle="--"
)

ax.axhline(
    LOG2FC_THRESHOLD,
    linewidth=0.6,
    linestyle=":"
)

ax.axhline(
    -LOG2FC_THRESHOLD,
    linewidth=0.6,
    linestyle=":"
)

ax.axvline(
    LOG2FC_THRESHOLD,
    linewidth=0.6,
    linestyle=":"
)

ax.axvline(
    -LOG2FC_THRESHOLD,
    linewidth=0.6,
    linestyle=":"
)

ax.set_xlabel(
    "Pseudomonas colistin consensus response (log2FC)",
    fontsize=11
)

ax.set_ylabel(
    "Klebsiella K56 colistin response (log2FC)",
    fontsize=11
)

ax.set_title(
    "A  Cross-species colistin-response comparison",
    loc="left",
    fontsize=12,
    fontweight="bold"
)

ax.legend(
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

fig.savefig(
    SCATTER_FIGURE,
    dpi=FIGURE_DPI,
    bbox_inches="tight"
)

plt.close(fig)

print(f"Saved scatterplot: {SCATTER_FIGURE}")


# ----------------------------------------------------------------
# 13. PANEL B: RESPONSE HEATMAP
# ----------------------------------------------------------------

if len(figure_gene_df) > 0:

    heatmap_matrix = (
        figure_gene_df
        .set_index("Figure_Gene_Label")[heatmap_columns]
    )

    # Keep conserved genes together and divergent genes together.
    row_labels = heatmap_matrix.index.tolist()

    heatmap_height = max(
        5.5,
        len(heatmap_matrix) * 0.32
    )

    fig, ax = plt.subplots(
        figsize=(7.4, heatmap_height)
    )

    image = ax.imshow(
        heatmap_matrix.values,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-max(
            2,
            np.nanmax(np.abs(heatmap_matrix.values))
        ),
        vmax=max(
            2,
            np.nanmax(np.abs(heatmap_matrix.values))
        )
    )

    ax.set_xticks(
        np.arange(len(heatmap_labels))
    )

    ax.set_xticklabels(
        heatmap_labels,
        rotation=30,
        ha="right"
    )

    ax.set_yticks(
        np.arange(len(row_labels))
    )

    ax.set_yticklabels(
        row_labels,
        fontsize=8
    )

    ax.set_title(
        "B  Conserved and divergent gene-response profiles",
        loc="left",
        fontsize=12,
        fontweight="bold"
    )

    # Horizontal separator between conserved and divergent genes.
    if conserved_count > 0 and divergent_count > 0:
        ax.axhline(
            conserved_count - 0.5,
            linewidth=1.5
        )

        ax.text(
            -0.65,
            (conserved_count - 1) / 2,
            "Conserved",
            rotation=90,
            va="center",
            ha="center",
            fontweight="bold"
        )

        ax.text(
            -0.65,
            conserved_count + (divergent_count - 1) / 2,
            "Divergent",
            rotation=90,
            va="center",
            ha="center",
            fontweight="bold"
        )

    colour_bar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.035,
        pad=0.03
    )

    colour_bar.set_label(
        "Response (log2FC)"
    )

    fig.tight_layout()

    fig.savefig(
        HEATMAP_FIGURE,
        dpi=FIGURE_DPI,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(f"Saved heatmap: {HEATMAP_FIGURE}")

else:
    warnings.warn(
        "No conserved or divergent genes were identified, "
        "so the heatmap was not generated."
    )


# ----------------------------------------------------------------
# 14. COMBINED FIGURE 3.7
# ----------------------------------------------------------------

if len(figure_gene_df) > 0:

    number_of_rows = len(figure_gene_df)

    combined_height = max(
        7.2,
        number_of_rows * 0.30
    )

    fig = plt.figure(
        figsize=(15.5, combined_height)
    )

    grid = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[1.15, 0.85],
        wspace=0.48
    )

    # -------------------------
    # Combined panel A
    # -------------------------

    ax1 = fig.add_subplot(grid[0, 0])

    ax1.scatter(
        below_threshold[
            "Pseudomonas_Colistin_Consensus_Log2FC"
        ],
        below_threshold[
            "Klebsiella_Colistin_Log2FC"
        ],
        s=15,
        alpha=0.25,
        label="Other comparable genes"
    )

    ax1.scatter(
        conserved_df[
            "Pseudomonas_Colistin_Consensus_Log2FC"
        ],
        conserved_df[
            "Klebsiella_Colistin_Log2FC"
        ],
        s=55,
        marker="o",
        edgecolors="black",
        linewidths=0.5,
        label=f"Conserved (n = {conserved_count})"
    )

    ax1.scatter(
        divergent_df[
            "Pseudomonas_Colistin_Consensus_Log2FC"
        ],
        divergent_df[
            "Klebsiella_Colistin_Log2FC"
        ],
        s=65,
        marker="^",
        edgecolors="black",
        linewidths=0.5,
        label=f"Divergent (n = {divergent_count})"
    )

    ax1.axhline(
        0,
        linewidth=0.8,
        linestyle="--"
    )

    ax1.axvline(
        0,
        linewidth=0.8,
        linestyle="--"
    )

    for threshold in [
        LOG2FC_THRESHOLD,
        -LOG2FC_THRESHOLD
    ]:
        ax1.axhline(
            threshold,
            linewidth=0.6,
            linestyle=":"
        )

        ax1.axvline(
            threshold,
            linewidth=0.6,
            linestyle=":"
        )

    ax1.set_xlabel(
        "Pseudomonas colistin consensus response (log2FC)"
    )

    ax1.set_ylabel(
        "Klebsiella K56 colistin response (log2FC)"
    )

    ax1.set_title(
        "A  Cross-species response comparison",
        loc="left",
        fontsize=12,
        fontweight="bold"
    )

    ax1.legend(
        frameon=False,
        loc="best"
    )

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # -------------------------
    # Combined panel B
    # -------------------------

    ax2 = fig.add_subplot(grid[0, 1])

    heatmap_matrix = (
        figure_gene_df
        .set_index("Figure_Gene_Label")[heatmap_columns]
    )

    maximum_value = max(
        2,
        np.nanmax(np.abs(heatmap_matrix.values))
    )

    image = ax2.imshow(
        heatmap_matrix.values,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-maximum_value,
        vmax=maximum_value
    )

    ax2.set_xticks(
        np.arange(len(heatmap_labels))
    )

    ax2.set_xticklabels(
        heatmap_labels,
        rotation=30,
        ha="right"
    )

    ax2.set_yticks(
        np.arange(len(heatmap_matrix.index))
    )

    ax2.set_yticklabels(
        heatmap_matrix.index,
        fontsize=8
    )

    ax2.set_title(
        "B  Gene-level response profiles",
        loc="left",
        fontsize=12,
        fontweight="bold"
    )

    if conserved_count > 0 and divergent_count > 0:
        ax2.axhline(
            conserved_count - 0.5,
            linewidth=1.5
        )

    colour_bar = fig.colorbar(
        image,
        ax=ax2,
        fraction=0.045,
        pad=0.04
    )

    colour_bar.set_label(
        "Response (log2FC)"
    )

    fig.savefig(
        COMBINED_FIGURE,
        dpi=FIGURE_DPI,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(f"Saved combined Figure 3.7: {COMBINED_FIGURE}")


# ----------------------------------------------------------------
# 15. VERIFY ORIGINAL COUNTS
# ----------------------------------------------------------------

print("\n" + "-" * 72)
print("LEGACY COUNT VERIFICATION")
print("-" * 72)

if (
    conserved_count == EXPECTED_CONSERVED_COUNT
    and
    divergent_count == EXPECTED_DIVERGENT_COUNT
):
    print(
        "SUCCESS: The strict analysis reproduced the original result:"
    )
    print(
        f"  {conserved_count} conserved genes"
    )
    print(
        f"  {divergent_count} divergent genes"
    )

else:
    print(
        "WARNING: The strict analysis did not reproduce the expected "
        "7 conserved and 15 divergent genes."
    )

    print(
        f"\nObserved: {conserved_count} conserved and "
        f"{divergent_count} divergent."
    )

    print(
        "\nThis means at least one of the following differs from the "
        "first script:"
    )

    print(
        "  1. The response threshold"
    )

    print(
        "  2. The columns compared"
    )

    print(
        "  3. Whether RNA-seq and Ribo-seq agreement was required"
    )

    print(
        "  4. Whether missing values were excluded"
    )

    print(
        "  5. Whether an ordinary fold-change column was converted "
        "to log2FC"
    )

    print(
        "\nDo not alter the threshold simply to force the expected counts. "
        "Compare these settings with the original script or original CSV "
        "outputs first."
    )


print("\n" + "=" * 72)
print("STEP 28A COMPLETE")
print("=" * 72 + "\n")