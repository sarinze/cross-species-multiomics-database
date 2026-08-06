import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram

print("\n========== STEP 25: CLUSTER HEATMAPS ==========\n")

INPUT_FILE = "Final_Clustered_Matrix.csv"

# Output file names
OUTPUT_MEAN_HEATMAP = "Mean_MultiOmics_Cluster_Profile_Heatmap.png"
OUTPUT_ZSCORE_HEATMAP = "Mean_MultiOmics_Zscore_Cluster_Heatmap.png"
OUTPUT_ROW_ZSCORE_HEATMAP = "Gene_Level_MultiOmics_Row_Zscore_Heatmap.png"
OUTPUT_DENDROGRAM = "Hierarchical_Cluster_Dendrogram.png"

df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

if "Final_Cluster" not in df.columns:
    raise ValueError("Final_Cluster column not found.")

metadata_keywords = [
    "gene", "locus", "tag", "name", "description",
    "product", "protein", "accession", "cluster"
]

metadata_cols = [
    col for col in df.columns
    if any(keyword in col.lower() for keyword in metadata_keywords)
]

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

omics_cols = [
    col for col in numeric_cols
    if col not in metadata_cols
    and col not in ["Final_Cluster", "Measured_Omics_Count"]
]

print(f"\nMetadata columns detected: {len(metadata_cols)}")
print(f"Omics columns used for heatmaps: {len(omics_cols)}")
print(omics_cols)

# -----------------------------
# Cluster mean profile matrix
# -----------------------------
cluster_profile = df.groupby("Final_Cluster")[omics_cols].mean()

cluster_sizes = df["Final_Cluster"].value_counts().sort_index()

# Order clusters by size, largest first
cluster_order = cluster_sizes.sort_values(ascending=False).index.tolist()
cluster_profile = cluster_profile.loc[cluster_order]

cluster_profile.to_csv("Cluster_Profile_Matrix_Mean.csv")
print("\nSaved: Cluster_Profile_Matrix_Mean.csv")

# -----------------------------
# Heatmap helper
# -----------------------------
def make_heatmap(data, output_file, title, colorbar_label):
    plt.figure(figsize=(14, 10))

    im = plt.imshow(
        data,
        aspect="auto",
        cmap="coolwarm",
        vmin=-3 if "Zscore" in output_file or "Zscore" in title else None,
        vmax=3 if "Zscore" in output_file or "Zscore" in title else None
    )

    plt.colorbar(im, label=colorbar_label)

    plt.xticks(
        ticks=np.arange(len(data.columns)),
        labels=data.columns,
        rotation=90,
        fontsize=14
    )

    plt.yticks(
        ticks=np.arange(len(data.index)),
        labels=data.index,
        fontsize=14
    )

    plt.xlabel("Multi-Omics Datasets")
    plt.ylabel("Hierarchical Gene Clusters")
    plt.title(title, fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved heatmap: {output_file}")

# -----------------------------
# Heatmap 1: mean raw values
# -----------------------------
make_heatmap(
    cluster_profile,
    OUTPUT_MEAN_HEATMAP,
    "Mean Multi-Omics Profiles of Hierarchical Gene Clusters",
    "Mean Omics Value"
)

# -----------------------------
# Heatmap 2: z-scored cluster means
# -----------------------------
cluster_profile_z = pd.DataFrame(
    StandardScaler().fit_transform(cluster_profile),
    index=cluster_profile.index,
    columns=cluster_profile.columns
)

cluster_profile_z.to_csv("Cluster_Profile_Matrix_Zscore.csv")
print("Saved: Cluster_Profile_Matrix_Zscore.csv")

make_heatmap(
    cluster_profile_z,
    OUTPUT_ZSCORE_HEATMAP,
    "Mean Multi-Omics Z-Score Profiles of Hierarchical Gene Clusters",
    "Mean Cluster Z-score"
)

# -----------------------------
# Heatmap 3: gene-level row z-score
# -----------------------------
gene_matrix = df[omics_cols].replace([np.inf, -np.inf], np.nan).fillna(0)

row_mean = gene_matrix.mean(axis=1)
row_std = gene_matrix.std(axis=1).replace(0, np.nan)

gene_row_z = gene_matrix.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0)
gene_row_z["Final_Cluster"] = df["Final_Cluster"]

gene_row_z = gene_row_z.sort_values("Final_Cluster")

plt.figure(figsize=(14, 10))

im = plt.imshow(
    gene_row_z[omics_cols],
    aspect="auto",
    cmap="coolwarm",
    vmin=-3,
    vmax=3
)

plt.colorbar(im, label="Row Z-score")

plt.xticks(
    ticks=np.arange(len(omics_cols)),
    labels=omics_cols,
    rotation=90,
    fontsize=14
)

plt.yticks([])
plt.xlabel("Multi-Omics Datasets")
plt.ylabel("Genes Ordered by Hierarchical Cluster")
plt.title(
    "Gene-Level Multi-Omics Heatmap (Genes Ordered by Hierarchical Cluster)",
    fontsize=14
)

plt.tight_layout()
plt.savefig(OUTPUT_ROW_ZSCORE_HEATMAP, dpi=300)
plt.close()

print(f"Saved heatmap: {OUTPUT_ROW_ZSCORE_HEATMAP}")

# -----------------------------
# Cluster dendrogram
# -----------------------------
Z = linkage(cluster_profile_z, method="ward")

plt.figure(figsize=(10, 8))

dendrogram(
    Z,
    labels=cluster_profile_z.index.astype(str).tolist(),
    leaf_rotation=0
)

plt.xlabel("Hierarchical Gene Clusters")
plt.ylabel("Distance")
plt.title("Hierarchical Relationships Among Gene Clusters", fontsize=14)

plt.tight_layout()
plt.savefig(OUTPUT_DENDROGRAM, dpi=300)
plt.close()

print(f"Saved dendrogram: {OUTPUT_DENDROGRAM}")

# -----------------------------
# Save long-format data
# -----------------------------
long_df = cluster_profile.reset_index().melt(
    id_vars="Final_Cluster",
    var_name="Omics_Feature",
    value_name="Mean_Value"
)

long_df.to_csv("Cluster_Heatmap_Data_Long.csv", index=False)
print("Saved: Cluster_Heatmap_Data_Long.csv")

print("\nCluster sizes:")
print(cluster_sizes.sort_values(ascending=False).to_string())

print("\n========== STEP 25 COMPLETE ==========\n")