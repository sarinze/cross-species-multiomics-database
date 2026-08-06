import pandas as pd
import matplotlib.pyplot as plt
import os

print("\n========== STEP 27: SCATTERPLOT ANALYSIS ==========\n")

INPUT_FILE = "Final_Clustered_Matrix.csv"
BIOLOGY_FILE = "Final_Biological_Cluster_Characterization.csv"

OUTPUT_DIR = "Scatterplot_Analysis_Outputs"
SUMMARY_OUTPUT = "Scatterplot_Analysis_Summary.csv"
OUTLIER_OUTPUT = "Scatterplot_Outlier_Genes.csv"

CLUSTER_COL = "Final_Cluster"

os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(INPUT_FILE)
biology = pd.read_csv(BIOLOGY_FILE)

print(f"Loaded matrix: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# Add biological fingerprint to each gene
df = df.merge(
    biology[[CLUSTER_COL, "Biological_Fingerprint", "Outlier_Status"]],
    on=CLUSTER_COL,
    how="left"
)

# ---------------------------------------------------
# Scatterplot comparisons
# ---------------------------------------------------

comparisons = [
    {
        "name": "Pseudomonas_Colistin_RNAseq_vs_Tobramycin_RNAseq",
        "x": "Colistin_RNAseq_Fold_Change",
        "y": "Tobramycin_RNAseq_Fold_Change",
        "question": "Do Pseudomonas genes responding to colistin also respond to tobramycin?"
    },
    {
        "name": "Pseudomonas_Colistin_RNAseq_vs_Colistin_Riboseq",
        "x": "Colistin_RNAseq_Fold_Change",
        "y": "Colistin_Riboseq_Fold_Change",
        "question": "Does colistin transcriptional response match ribosome-level response?"
    },
    {
        "name": "Pseudomonas_Tobramycin_RNAseq_vs_Tobramycin_Riboseq",
        "x": "Tobramycin_RNAseq_Fold_Change",
        "y": "Tobramycin_Riboseq_Fold_Change",
        "question": "Does tobramycin transcriptional response match ribosome-level response?"
    },
    {
        "name": "Pseudomonas_Colistin_RNAseq_vs_Persister_TnSeq",
        "x": "Colistin_RNAseq_Fold_Change",
        "y": "PAO1_Persister_All_Mean_SI",
        "question": "Do colistin-responsive genes also influence Pseudomonas persister fitness?"
    },
    {
        "name": "Pseudomonas_Tobramycin_RNAseq_vs_Persister_TnSeq",
        "x": "Tobramycin_RNAseq_Fold_Change",
        "y": "PAO1_Persister_All_Mean_SI",
        "question": "Do tobramycin-responsive genes also influence Pseudomonas persister fitness?"
    },
    {
        "name": "Pseudomonas_Colistin_vs_Klebsiella_K56_Colistin",
        "x": "Colistin_RNAseq_Fold_Change",
        "y": "K56_vs_Colistin_Log2FC",
        "question": "Do orthologous genes show similar colistin responses in Pseudomonas and Klebsiella?"
    },
    {
        "name": "Klebsiella_K56_Colistin_vs_K56_Combination",
        "x": "K56_vs_Colistin_Log2FC",
        "y": "K56_vs_Combination_Log2FC",
        "question": "How does Klebsiella K56 colistin response compare with combination response?"
    },
    {
        "name": "Klebsiella_K56_Response_vs_KPPR1_TnSeq",
        "x": "K56_vs_Combination_Log2FC",
        "y": "KPPR1_TnSeq_Log2FC_Output_Input",
        "question": "Does Klebsiella antibiotic response relate to KPPR1 mutant fitness?"
    },
    {
        "name": "Pseudomonas_Persister_TnSeq_vs_KPPR1_TnSeq",
        "x": "PAO1_Persister_All_Mean_SI",
        "y": "KPPR1_TnSeq_Log2FC_Output_Input",
        "question": "Do Pseudomonas persister fitness and Klebsiella KPPR1 fitness show similar trends?"
    }
]

# ---------------------------------------------------
# Helper functions
# ---------------------------------------------------

def get_gene_label(row):
    for col in ["Gene_Label", "Pseudomonas_Gene_Name", "Klebsiella_Gene_Name", "Pseudomonas_Original_Locus_Tag", "Klebsiella_Locus_Tag"]:
        if col in row and pd.notna(row[col]) and str(row[col]).strip() != "":
            return str(row[col])
    return "Unknown"


def make_scatterplot(data, x_col, y_col, name, question):
    plot_df = data[[x_col, y_col, CLUSTER_COL, "Gene_Label"]].copy()
    plot_df = plot_df.dropna(subset=[x_col, y_col])

    if plot_df.empty:
        print(f"Skipping {name}: no valid data.")
        return None, None

    # Outlier score based on distance from origin
    plot_df["Scatter_Outlier_Score"] = (plot_df[x_col].abs() + plot_df[y_col].abs())

    top_outliers = (
        plot_df
        .sort_values("Scatter_Outlier_Score", ascending=False)
        .head(15)
        .copy()
    )

    # Quadrant counts
    both_high = ((plot_df[x_col] > 1) & (plot_df[y_col] > 1)).sum()
    x_high_y_low = ((plot_df[x_col] > 1) & (plot_df[y_col] < -1)).sum()
    x_low_y_high = ((plot_df[x_col] < -1) & (plot_df[y_col] > 1)).sum()
    both_low = ((plot_df[x_col] < -1) & (plot_df[y_col] < -1)).sum()
    neutral = len(plot_df) - both_high - x_high_y_low - x_low_y_high - both_low

    plt.figure(figsize=(10, 8))

    scatter = plt.scatter(
        plot_df[x_col],
        plot_df[y_col],
        c=plot_df[CLUSTER_COL],
        cmap="tab20",
        s=15,
        alpha=0.65
    )

    plt.axhline(0, linewidth=0.8)
    plt.axvline(0, linewidth=0.8)

    plt.axhline(1, linestyle="--", linewidth=0.7)
    plt.axhline(-1, linestyle="--", linewidth=0.7)
    plt.axvline(1, linestyle="--", linewidth=0.7)
    plt.axvline(-1, linestyle="--", linewidth=0.7)

    # Label top 10 outliers
    for _, row in top_outliers.head(10).iterrows():
        label = get_gene_label(row)
        plt.text(
            row[x_col],
            row[y_col],
            label,
            fontsize=7,
            alpha=0.8
        )

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(name.replace("_", " "))

    cbar = plt.colorbar(scatter)
    cbar.set_label("Final Cluster")

    plt.tight_layout()

    output_plot = os.path.join(OUTPUT_DIR, f"{name}.png")
    plt.savefig(output_plot, dpi=300)
    plt.close()

    summary_row = {
        "Comparison": name,
        "Question": question,
        "X_Column": x_col,
        "Y_Column": y_col,
        "N_Genes_Plotted": len(plot_df),
        "Both_High_Count": both_high,
        "X_High_Y_Low_Count": x_high_y_low,
        "X_Low_Y_High_Count": x_low_y_high,
        "Both_Low_Count": both_low,
        "Neutral_or_Mixed_Count": neutral,
        "Plot_File": output_plot
    }

    outlier_rows = []

    for _, row in top_outliers.iterrows():
        outlier_rows.append({
            "Comparison": name,
            "Gene_Label": get_gene_label(row),
            "Final_Cluster": row[CLUSTER_COL],
            "X_Column": x_col,
            "X_Value": row[x_col],
            "Y_Column": y_col,
            "Y_Value": row[y_col],
            "Scatter_Outlier_Score": row["Scatter_Outlier_Score"]
        })

    return summary_row, outlier_rows


# ---------------------------------------------------
# Run all scatterplots
# ---------------------------------------------------

summary_rows = []
all_outlier_rows = []

for comp in comparisons:
    x = comp["x"]
    y = comp["y"]

    if x not in df.columns or y not in df.columns:
        print(f"Skipping {comp['name']}: missing column.")
        continue

    print(f"\nCreating scatterplot: {comp['name']}")

    summary_row, outlier_rows = make_scatterplot(
        df,
        x,
        y,
        comp["name"],
        comp["question"]
    )

    if summary_row is not None:
        summary_rows.append(summary_row)
        all_outlier_rows.extend(outlier_rows)

summary_df = pd.DataFrame(summary_rows)
outlier_df = pd.DataFrame(all_outlier_rows)

summary_path = os.path.join(OUTPUT_DIR, SUMMARY_OUTPUT)
outlier_path = os.path.join(OUTPUT_DIR, OUTLIER_OUTPUT)

summary_df.to_csv(summary_path, index=False)
outlier_df.to_csv(outlier_path, index=False)

print(f"\nSaved summary: {summary_path}")
print(f"Saved outlier genes: {outlier_path}")

print("\nScatterplot summary:")
print(summary_df.to_string(index=False))

print("\n========== STEP 27 COMPLETE ==========")