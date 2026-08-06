import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

print("\n========== STEP 26: PCA AND CLUSTER VISUALIZATION ==========\n")

INPUT_FILE = "Final_Clustered_Matrix.csv"
BIOLOGY_FILE = "Final_Biological_Cluster_Characterization.csv"

PCA_COORDS_OUTPUT = "PCA_Cluster_Coordinates.csv"
PCA_LOADINGS_OUTPUT = "PCA_Feature_Loadings.csv"
PCA_PLOT_OUTPUT = "PCA_Cluster_Visualization.png"
PCA_CENTROID_OUTPUT = "PCA_Cluster_Centroids.csv"
PCA_CENTROID_PLOT_OUTPUT = "PCA_Cluster_Centroid_Visualization.png"

CLUSTER_COL = "Final_Cluster"

df = pd.read_csv(INPUT_FILE)
biology = pd.read_csv(BIOLOGY_FILE)

print(f"Loaded matrix: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# -----------------------------
# Define omics columns only
# -----------------------------

omics_cols = [
    "Colistin_RNAseq_Fold_Change",
    "Colistin_Riboseq_Fold_Change",
    "Tobramycin_RNAseq_Fold_Change",
    "Tobramycin_Riboseq_Fold_Change",
    "PAO1_Persister_All_Rep1_SI",
    "PAO1_Persister_All_Rep2_SI",
    "PAO1_Persister_All_Rep3_SI",
    "PAO1_Persister_All_Mean_SI",
    "K56_vs_Colistin_Log2FC",
    "K56_vs_Combination_Log2FC",
    "Colistin_vs_Combination_Log2FC",
    "KPPR1_TnSeq_Fold_Change",
    "KPPR1_TnSeq_Log2FC_Output_Input"
]

missing = [c for c in omics_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing omics columns: {missing}")

X = df[omics_cols].copy()

# Fill missing values with 0 because missing means not measured/no signal after matrix preparation
X = X.fillna(0)

# -----------------------------
# Scale data
# -----------------------------

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# -----------------------------
# PCA
# -----------------------------

pca = PCA(n_components=2)
pca_result = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame({
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
    CLUSTER_COL: df[CLUSTER_COL],
    "Gene_Label": df["Gene_Label"],
    "Pseudomonas_Original_Locus_Tag": df["Pseudomonas_Original_Locus_Tag"],
    "Pseudomonas_Gene_Name": df["Pseudomonas_Gene_Name"],
    "Klebsiella_Locus_Tag": df["Klebsiella_Locus_Tag"],
    "Klebsiella_Gene_Name": df["Klebsiella_Gene_Name"],
})

# Add biological fingerprint
pca_df = pca_df.merge(
    biology[[CLUSTER_COL, "Biological_Fingerprint", "Dominant_Driver", "Outlier_Status"]],
    on=CLUSTER_COL,
    how="left"
)

pca_df.to_csv(PCA_COORDS_OUTPUT, index=False)

print(f"\nSaved: {PCA_COORDS_OUTPUT}")

print("\nExplained variance:")
print(f"PC1: {pca.explained_variance_ratio_[0] * 100:.2f}%")
print(f"PC2: {pca.explained_variance_ratio_[1] * 100:.2f}%")
print(f"Total PC1 + PC2: {pca.explained_variance_ratio_[:2].sum() * 100:.2f}%")

# -----------------------------
# PCA feature loadings
# -----------------------------

loadings = pd.DataFrame(
    pca.components_.T,
    columns=["PC1_Loading", "PC2_Loading"],
    index=omics_cols
).reset_index().rename(columns={"index": "Feature"})

loadings["PC1_Abs_Loading"] = loadings["PC1_Loading"].abs()
loadings["PC2_Abs_Loading"] = loadings["PC2_Loading"].abs()

loadings.to_csv(PCA_LOADINGS_OUTPUT, index=False)

print(f"Saved: {PCA_LOADINGS_OUTPUT}")

print("\nTop features contributing to PC1:")
print(
    loadings.sort_values("PC1_Abs_Loading", ascending=False)
    [["Feature", "PC1_Loading"]]
    .head(10)
    .to_string(index=False)
)

print("\nTop features contributing to PC2:")
print(
    loadings.sort_values("PC2_Abs_Loading", ascending=False)
    [["Feature", "PC2_Loading"]]
    .head(10)
    .to_string(index=False)
)

# -----------------------------
# PCA scatter plot: all genes
# -----------------------------

plt.figure(figsize=(11, 8))

scatter = plt.scatter(
    pca_df["PC1"],
    pca_df["PC2"],
    c=pca_df[CLUSTER_COL],
    cmap="tab20",
    s=14,
    alpha=0.65
)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
plt.title("PCA Visualization of Multi-Omics Clusters")

cbar = plt.colorbar(scatter)
cbar.set_label("Final Cluster")

plt.tight_layout()
plt.savefig(PCA_PLOT_OUTPUT, dpi=300)
plt.close()

print(f"Saved: {PCA_PLOT_OUTPUT}")

# -----------------------------
# Cluster centroid plot
# -----------------------------

centroids = (
    pca_df
    .groupby(CLUSTER_COL)
    .agg(
        PC1=("PC1", "mean"),
        PC2=("PC2", "mean"),
        Cluster_Size=("Gene_Label", "count")
    )
    .reset_index()
)

centroids = centroids.merge(
    biology[[CLUSTER_COL, "Biological_Fingerprint", "Dominant_Driver", "Outlier_Status"]],
    on=CLUSTER_COL,
    how="left"
)

centroids.to_csv(PCA_CENTROID_OUTPUT, index=False)

plt.figure(figsize=(12, 9))

plt.scatter(
    centroids["PC1"],
    centroids["PC2"],
    s=centroids["Cluster_Size"] * 0.8,
    alpha=0.65
)

for _, row in centroids.iterrows():
    plt.text(
        row["PC1"],
        row["PC2"],
        str(int(row[CLUSTER_COL])),
        fontsize=9,
        ha="center",
        va="center"
    )

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
plt.title("PCA Cluster Centroids with Cluster Sizes")

plt.tight_layout()
plt.savefig(PCA_CENTROID_PLOT_OUTPUT, dpi=300)
plt.close()

print(f"Saved: {PCA_CENTROID_OUTPUT}")
print(f"Saved: {PCA_CENTROID_PLOT_OUTPUT}")

print("\nCluster centroid summary:")
print(
    centroids[
        [
            CLUSTER_COL,
            "Cluster_Size",
            "PC1",
            "PC2",
            "Dominant_Driver",
            "Outlier_Status",
            "Biological_Fingerprint"
        ]
    ]
    .sort_values("Cluster_Size", ascending=False)
    .to_string(index=False)
)

print("\n========== STEP 26 COMPLETE ==========")