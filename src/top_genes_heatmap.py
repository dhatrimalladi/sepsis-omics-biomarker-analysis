import pandas as pd
import matplotlib.pyplot as plt

expression_file = "data/expression_sepsis_d1_vs_healthy.csv"
top_genes_file = "results/top_50_differentially_expressed_genes.csv"
output_file = "results/figures/top_50_genes_heatmap.png"

# Load expression data
expression = pd.read_csv(expression_file)

# Remove GEO end marker
expression = expression[
    expression["ID_REF"] != "!series_matrix_table_end"
]

# Load top 50 differential-expression results
top = pd.read_csv(top_genes_file)

# Keep probe IDs that are present in the expression matrix
probe_ids = top["probe_id"].tolist()
probe_ids = [
    probe for probe in probe_ids
    if probe in expression["ID_REF"].values
]

# Extract expression values
heatmap_data = expression[
    expression["ID_REF"].isin(probe_ids)
].set_index("ID_REF")

# Keep the same order as the ranked top-50 table
heatmap_data = heatmap_data.loc[probe_ids]

# Z-score each gene across samples
heatmap_z = heatmap_data.sub(
    heatmap_data.mean(axis=1), axis=0
).div(
    heatmap_data.std(axis=1), axis=0
)

plt.figure(figsize=(14, 10))

plt.imshow(
    heatmap_z,
    aspect="auto",
    interpolation="nearest"
)

plt.colorbar(label="Expression Z-score")
plt.xlabel("Samples")
plt.ylabel("Top 50 probes")
plt.title("Top 50 Differentially Expressed Probes")

plt.xticks(
    range(len(heatmap_z.columns)),
    heatmap_z.columns,
    rotation=90,
    fontsize=6
)

plt.yticks(
    range(len(heatmap_z.index)),
    heatmap_z.index,
    fontsize=6
)

plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()

print("Heatmap saved to:", output_file)
print("Probes plotted:", len(probe_ids))
