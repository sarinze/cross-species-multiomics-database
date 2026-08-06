import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

print("\n========== STEP 26b: PCA ANNOTATED BIOLOGICAL CLUSTER PLOT ==========\n")

INPUT_FILE = "Final_Clustered_Matrix.csv"
BIOLOGY_FILE = "Final_Biological_Cluster_Characterization.csv"

OUTPUT_CENTROIDS = "PCA_Annotated_Biological_Centroids.csv"
OUTPUT_PLOT = "PCA_Annotated_Biological_Cluster_Plot.png"

CLUSTER_COL = "Final_Cluster"

df = pd.read_csv(INPUT_FILE)
biology = pd.read_csv(BIOLOGY_FILE)

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

X = df[omics_cols].fillna(0)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame({
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
    CLUSTER_COL: df[CLUSTER_COL],
    "Gene_Label": df["Gene_Label"]
})

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
    biology[
        [
            CLUSTER_COL,
            "Biological_Fingerprint",
            "Dominant_Driver",
            "Secondary_Driver",
            "Outlier_Status"
        ]
    ],
    on=CLUSTER_COL,
    how="left"
)

def simplify_fingerprint(fp):
    """
    Makes a short label for the PCA plot.
    """
    parts = str(fp).split("|")
    parts = [p.strip() for p in parts]

    important = []

    for p in parts:
        if p.startswith("HIGH "):
            important.append(p)
        elif p.startswith("MOD HIGH "):
            important.append(p)
        elif p.startswith("LOW "):
            important.append(p)
        elif p.startswith("MOD LOW "):
            important.append(p)
        elif "OUTLIER" in p:
            important.append("OUTLIER")

    if len(important) == 0:
        return "NEUTRAL PROFILE"

    return "\n".join(important[:3])


centroids["Short_Biological_Label"] = centroids["Biological_Fingerprint"].apply(simplify_fingerprint)

centroids["Plot_Label"] = centroids.apply(
    lambda row: f"C{int(row[CLUSTER_COL])}\n{row['Short_Biological_Label']}",
    axis=1
)

centroids.to_csv(OUTPUT_CENTROIDS, index=False)

print(f"Saved: {OUTPUT_CENTROIDS}")

print("\nExplained variance:")
print(f"PC1: {pca.explained_variance_ratio_[0] * 100:.2f}%")
print(f"PC2: {pca.explained_variance_ratio_[1] * 100:.2f}%")
print(f"Total PC1 + PC2: {pca.explained_variance_ratio_[:2].sum() * 100:.2f}%")

plt.figure(figsize=(15, 11))

plt.scatter(
    centroids["PC1"],
    centroids["PC2"],
    s=centroids["Cluster_Size"] * 0.9,
    alpha=0.6
)

for _, row in centroids.iterrows():
    plt.text(
        row["PC1"],
        row["PC2"],
        row["Plot_Label"],
        fontsize=8,
        ha="center",
        va="center"
    )

plt.axhline(0, linewidth=0.7)
plt.axvline(0, linewidth=0.7)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
plt.title("PCA Cluster Centroids Annotated by Biological Fingerprint")

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300)
plt.close()

print(f"Saved: {OUTPUT_PLOT}")

print("\nAnnotated centroid summary:")
print(
    centroids[
        [
            CLUSTER_COL,
            "Cluster_Size",
            "PC1",
            "PC2",
            "Short_Biological_Label",
            "Outlier_Status"
        ]
    ].to_string(index=False)
)

print("\n========== STEP 26b COMPLETE ==========")