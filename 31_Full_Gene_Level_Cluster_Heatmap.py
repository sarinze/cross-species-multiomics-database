# ============================================================
# STEP 31: FULL GENE-LEVEL CLUSTER HEATMAP
# Shows ALL genes together, grouped by Final_Cluster
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# INPUT FILE
# ------------------------------------------------------------

INPUT_FILE = "Final_Clustered_Matrix.csv"

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

df = pd.read_csv(INPUT_FILE)

print("========== FULL GENE-LEVEL CLUSTER HEATMAP ==========")
print("Loaded:", INPUT_FILE)
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])

# ------------------------------------------------------------
# CHECK FINAL_CLUSTER EXISTS
# ------------------------------------------------------------

if "Final_Cluster" not in df.columns:
    raise ValueError("ERROR: Final_Cluster column not found in the file.")

# ------------------------------------------------------------
# SORT GENES BY CLUSTER
# ------------------------------------------------------------

df = df.sort_values("Final_Cluster")

# ------------------------------------------------------------
# IDENTIFY METADATA COLUMNS TO EXCLUDE
# ------------------------------------------------------------

metadata_cols = [
    "Gene_Label",
    "Final_Cluster",
    "Cluster_Label",
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name"
]

metadata_cols = [col for col in metadata_cols if col in df.columns]

# ------------------------------------------------------------
# KEEP ONLY NUMERIC OMICS COLUMNS
# ------------------------------------------------------------

omics_cols = [
    col for col in df.columns
    if col not in metadata_cols
    and pd.api.types.is_numeric_dtype(df[col])
]

# Remove non-response numeric summary columns
omics_cols = [
    col for col in omics_cols
    if col not in ["Measured_Omics_Count"]
]

print("Omics columns used:", len(omics_cols))
print(omics_cols)

# ------------------------------------------------------------
# CREATE HEATMAP MATRIX
# ------------------------------------------------------------

heatmap_data = df[omics_cols]

# Replace missing values with 0 for plotting
heatmap_data = heatmap_data.fillna(0)

# ------------------------------------------------------------
# CREATE FULL HEATMAP
# ------------------------------------------------------------

plt.figure(figsize=(18, 30))

sns.heatmap(
    heatmap_data,
    cmap="RdBu_r",
    center=0,
    yticklabels=False,
    xticklabels=True,
    cbar_kws={"label": "Omics Value"}
)

plt.title(
    "Full Gene-Level Multi-Omics Heatmap Grouped by Cluster",
    fontsize=18,
    weight="bold"
)

plt.xlabel("Omics Features")
plt.ylabel("Genes grouped by Final Cluster")

plt.xticks(rotation=45, ha="right")

plt.tight_layout()

OUTPUT = "Full_Gene_Level_All_Clusters_Heatmap.png"
plt.savefig(OUTPUT, dpi=300)

print("Saved:", OUTPUT)
plt.show()

print("========== COMPLETE ==========")