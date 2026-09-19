import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

expression_file = "data/expression_sepsis_d1_vs_healthy.csv"
metadata_file = "data/sample_metadata.csv"

# Load expression data
df = pd.read_csv(expression_file)

# Set probe IDs as row names
expression = df.set_index("ID_REF")

# Transpose: samples become rows, probes become columns
X = expression.T

# Load metadata
metadata = pd.read_csv(metadata_file)
metadata = metadata.set_index("sample_id")

# Keep metadata in the same order as expression samples
groups = metadata.loc[X.index, "group"]

# Standardize probe measurements
X_scaled = StandardScaler().fit_transform(X)

# PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame(
    principal_components,
    index=X.index,
    columns=["PC1", "PC2"]
)

pca_df["group"] = groups

# Plot
plt.figure(figsize=(10, 7))

for group in pca_df["group"].unique():
    subset = pca_df[pca_df["group"] == group]
    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        s=70
    )

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
plt.title("PCA of Sepsis D1 and Healthy Samples")
plt.legend()
plt.tight_layout()

output_file = "results/figures/pca_sepsis_vs_healthy.png"
plt.savefig(output_file, dpi=300)
plt.close()

print(f"PC1 variance explained: {pca.explained_variance_ratio_[0] * 100:.2f}%")
print(f"PC2 variance explained: {pca.explained_variance_ratio_[1] * 100:.2f}%")
print(f"PCA figure saved to: {output_file}")
