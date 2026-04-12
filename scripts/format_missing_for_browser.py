import pandas as pd
import os

input_file = "/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt"
output_file = "/mnt/temp2/kronckbm/gvf_output/missing_baseline_browser.csv"

pmids = []
with open(input_file, "r") as f:
    for line in f:
        pmid = line.strip()
        if pmid:
            pmids.append(pmid)

df = pd.DataFrame({"PMID": pmids, "Status": ["missing"] * len(pmids)})
df.to_csv(output_file, index=False)
print(f"Created {output_file} with {len(df)} rows")
