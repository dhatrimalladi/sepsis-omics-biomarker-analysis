import pandas as pd
from scipy.stats import ttest_ind

expression_file = "data/expression_sepsis_d1_vs_healthy.csv"
metadata_file = "data/sample_metadata.csv"
annotation_file = "data/probe_gene_annotation.csv"

# Load expression data
expr = pd.read_csv(expression_file)

# Remove the GEO end marker
expr = expr[expr["ID_REF"] != "!series_matrix_table_end"].copy()

# Load metadata
meta = pd.read_csv(metadata_file)
meta = meta[meta["group"].isin(["Sepsis_D1", "Healthy"])]

sepsis_samples = meta.loc[meta["group"] == "Sepsis_D1", "sample_id"].tolist()
healthy_samples = meta.loc[meta["group"] == "Healthy", "sample_id"].tolist()

# Keep only samples present in the expression matrix
sepsis_samples = [x for x in sepsis_samples if x in expr.columns]
healthy_samples = [x for x in healthy_samples if x in expr.columns]

results = []

for _, row in expr.iterrows():
    probe = row["ID_REF"]

    sepsis_values = pd.to_numeric(row[sepsis_samples], errors="coerce")
    healthy_values = pd.to_numeric(row[healthy_samples], errors="coerce")

    sepsis_mean = sepsis_values.mean()
    healthy_mean = healthy_values.mean()

    # Difference on the log-expression scale
    mean_difference = sepsis_mean - healthy_mean

    _, p_value = ttest_ind(
        sepsis_values,
        healthy_values,
        equal_var=False,
        nan_policy="omit"
    )

    results.append({
        "probe_id": probe,
        "sepsis_mean": sepsis_mean,
        "healthy_mean": healthy_mean,
        "mean_difference": mean_difference,
        "p_value": p_value
    })

results = pd.DataFrame(results)

# Multiple-testing correction using Benjamini-Hochberg FDR
results = results.sort_values("p_value").reset_index(drop=True)
n = len(results)
results["rank"] = range(1, n + 1)
results["fdr"] = (results["p_value"] * n / results["rank"]).clip(upper=1)

# Add gene annotation
annotation = pd.read_csv(annotation_file)

results = results.merge(
    annotation[["probe_id", "gene_symbol", "entrez_gene_id"]],
    on="probe_id",
    how="left"
)

# Put important columns first
results = results[
    [
        "probe_id",
        "gene_symbol",
        "entrez_gene_id",
        "sepsis_mean",
        "healthy_mean",
        "mean_difference",
        "p_value",
        "fdr"
    ]
]

output_file = "results/differential_expression_results.csv"
results.to_csv(output_file, index=False)

print("Sepsis samples:", len(sepsis_samples))
print("Healthy samples:", len(healthy_samples))
print("Probes analyzed:", len(results))
print("Results saved to:", output_file)

print("\nTop 10 results:")
print(results.head(10).to_string(index=False))
