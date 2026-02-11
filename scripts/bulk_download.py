import os

import requests

# Replace these with actual values or arguments
pmids_file = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/scripts/sample_pmids.txt"
output_dir = "output_directory"

# Ensure output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

with open(pmids_file, "r") as f:
    pmids = f.read().splitlines()

for pmid in pmids:
    url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/bioC/?pmid={pmid}"
    response = requests.get(url)
    if response.status_code == 200:
        # Save the XML data to the output directory
        with open(os.path.join(output_dir, f"{pmid}.xml"), "wb") as out_file:
            out_file.write(response.content)
    else:
        print(f"Failed to retrieve {pmid}: {response.status_code}")
