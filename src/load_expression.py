import gzip
import pandas as pd

expression_file = "data/GSE137340_series_matrix.txt.gz"
metadata_file = "data/sample_metadata.csv"

# Read sample metadata
metadata = pd.read_csv(metadata_file)

# Keep only Sepsis_D1 and Healthy samples
selected = metadata[metadata["group"].isin(["Sepsis_D1", "Healthy"])]

selected_samples = selected["sample_id"].tolist()

print(f"Selected samples: {len(selected_samples)}")
print(selected["group"].value_counts())

# Find the expression matrix header and data
with gzip.open(expression_file, "rt") as f:
    lines = f.readlines()

header_index = None

for i, line in enumerate(lines):
    if line.startswith('"ID_REF"'):
        header_index = i
        break

if header_index is None:
    raise ValueError("Expression matrix header not found")

# Read expression matrix
expression = pd.read_csv(
    expression_file,
    sep="\t",
    compression="gzip",
    skiprows=header_index
)

# Keep ID_REF and selected samples
columns_to_keep = ["ID_REF"] + selected_samples
expression = expression[columns_to_keep]

# Save selected expression matrix
output_file = "data/expression_sepsis_d1_vs_healthy.csv"
expression.to_csv(output_file, index=False)

print(f"Expression matrix saved to: {output_file}")
print(f"Matrix shape: {expression.shape}")
