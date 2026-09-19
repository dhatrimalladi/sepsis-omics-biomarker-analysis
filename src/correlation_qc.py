import pandas as pd
import matplotlib.pyplot as plt

expression_file = "data/expression_sepsis_d1_vs_healthy.csv"

# Load expression data
df = pd.read_csv(expression_file)

# Remove probe IDs
expression = df.set_index("ID_REF")

# Samples as rows
samples = expression.T

# Calculate Pearson correlation between samples
correlation = samples.T.corr()

# Print correlation of the candidate sample with all other samples
candidate = "GSM4075794"

print(f"\nCorrelations for {candidate}:")
print(correlation[candidate].sort_values().to_string())

# Average correlation with all other samples
mean_correlation = correlation[candidate].drop(candidate).mean()

print(f"\nMean correlation of {candidate} with other samples: {mean_correlation:.4f}")

# Create heatmap
plt.figure(figsize=(12, 10))

plt.imshow(
    correlation,
    aspect="auto",
    interpolation="nearest"
)

plt.colorbar(label="Pearson correlation")

plt.xticks(
    range(len(correlation.columns)),
    correlation.columns,
    rotation=90,
    fontsize=6
)

plt.yticks(
    range(len(correlation.index)),
    correlation.index,
    fontsize=6
)

plt.title("Sample-to-Sample Expression Correlation")
plt.tight_layout()

output_file = "results/figures/sample_correlation_heatmap.png"
plt.savefig(output_file, dpi=300)
plt.close()

print(f"\nCorrelation heatmap saved to: {output_file}")
