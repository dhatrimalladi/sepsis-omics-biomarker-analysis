import pandas as pd
import matplotlib.pyplot as plt

# Load expression matrix
file_path = "data/expression_sepsis_d1_vs_healthy.csv"
df = pd.read_csv(file_path)

# Remove probe ID column
expression = df.drop(columns=["ID_REF"])

# Create boxplot
plt.figure(figsize=(16, 6))
expression.boxplot(rot=90)

plt.title("Expression Distribution Across Samples")
plt.xlabel("Samples")
plt.ylabel("Expression Value")

plt.tight_layout()

# Save figure
output_file = "results/figures/expression_distribution.png"
plt.savefig(output_file, dpi=300)
plt.close()

print(f"QC figure saved to: {output_file}")
