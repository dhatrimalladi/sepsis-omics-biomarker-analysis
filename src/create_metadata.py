import gzip
import re
import pandas as pd

file_path = "data/GSE137340_series_matrix.txt.gz"

with gzip.open(file_path, "rt") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = re.findall(r'"([^"]+)"', line)
        elif line.startswith("!Sample_geo_accession"):
            gsm_ids = re.findall(r'"([^"]+)"', line)
            break

metadata = pd.DataFrame({
    "sample_id": gsm_ids,
    "sample_title": titles
})

def assign_group(title):
    if title.endswith("D1"):
        return "Sepsis_D1"
    elif title == "ctrl_sm" or title.startswith("C"):
        return "Healthy"
    else:
        return "Exclude"

metadata["group"] = metadata["sample_title"].apply(assign_group)

metadata.to_csv("data/sample_metadata.csv", index=False)

print(metadata.to_string(index=False))
print("\nGroup counts:")
print(metadata["group"].value_counts())
