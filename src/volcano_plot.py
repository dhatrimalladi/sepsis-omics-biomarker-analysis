import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input_file = "results/differential_expression_results.csv"
output_file = "results/figures/volcano_plot.png"

df = pd.read_csv(input_file)

# Calculate -log10(FDR)
df["neg_log10_fdr"] = -np.log10(df["fdr"].clip(lower=1e-300))

# Define significance
df["significant"] = (
    (df["fdr"] < 0.05) &
    (df["mean_difference"].abs() >= 1)
)

plt.figure(figsize=(10, 7))

plt.scatter(
    df["mean_difference"],
    df["neg_log10_fdr"],
    s=8,
    alpha=0.5
)

sig = df[df["significant"]]

plt.scatter(
    sig["mean_difference"],
    sig["neg_log10_fdr"],
    s=10,
    alpha=0.7
)

plt.axvline(1, linestyle="--")
plt.axvline(-1, linestyle="--")
plt.axhline(-np.log10(0.05), linestyle="--")

plt.xlabel("Mean expression difference (Sepsis D1 - Healthy)")
plt.ylabel("-log10(FDR)")
plt.title("Differential Expression: Sepsis D1 vs Healthy")

plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()

print("Significant probes:", len(sig))
print("Volcano plot saved to:", output_file)
