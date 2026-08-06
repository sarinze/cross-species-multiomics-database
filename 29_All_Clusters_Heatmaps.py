# ============================================================
# STEP 25B: ALL CLUSTERS HEATMAPS
# Creates:
# 1. Standard heatmap of all cluster profiles
# 2. Clustered heatmap / clustermap with dendrogram
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# INPUT FILE
# ------------------------------------------------------------

INPUT_FILE = "Cluster_Profile_Matrix_Mean.csv"

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

df = pd.read_csv(INPUT_FILE)

print("========== ALL CLUSTERS HEATMAP ==========")
print("Loaded:", INPUT_FILE)
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])

# ------------------------------------------------------------
# SET CLUSTER COLUMN AS INDEX
# ------------------------------------------------------------

df = df.set_index("Final_Cluster")

# Make sure cluster order is numerical
df = df.sort_index()

# ------------------------------------------------------------
# 1. STANDARD HEATMAP
# ------------------------------------------------------------

plt.figure(figsize=(18, 10))

sns.heatmap(
    df,
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    linecolor="lightgrey",
    annot=False,
    cbar_kws={"label": "Mean Omics Value"}
)

plt.title("Average Multi-Omics Profile Across All Clusters", fontsize=18, weight="bold")
plt.xlabel("Omics Features", fontsize=12)
plt.ylabel("Final Cluster", fontsize=12)

plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)

plt.tight_layout()

plt.savefig("All_Clusters_Standard_Heatmap.png", dpi=300)
print("Saved: All_Clusters_Standard_Heatmap.png")

plt.show()

# ------------------------------------------------------------
# 2. CLUSTERED HEATMAP / CLUSTERMAP
# ------------------------------------------------------------

g = sns.clustermap(
    df,
    cmap="RdBu_r",
    center=0,
    figsize=(18, 10),
    linewidths=0.3,
    linecolor="lightgrey",
    cbar_kws={"label": "Mean Omics Value"},
    dendrogram_ratio=0.15
)

g.fig.suptitle(
    "Hierarchical Clustered Heatmap of Multi-Omics Cluster Profiles",
    fontsize=18,
    weight="bold",
    y=1.03
)

plt.savefig("All_Clusters_Clustered_Heatmap.png", dpi=300, bbox_inches="tight")
print("Saved: All_Clusters_Clustered_Heatmap.png")

plt.show()

print("========== COMPLETE ==========")