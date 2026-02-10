#!/usr/bin/env python3
"""Analyze gold standard Excel file for recall metrics."""
import pandas as pd
import os
import glob
import re
import json

# Read gold standard Excel file
xls_path = "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
df = pd.read_excel(xls_path)

print("=== GOLD STANDARD ANALYSIS ===\n")
print(f"Total rows (carriers): {len(df)}")
print(f"\nColumns: {list(df.columns)}")

# Identify key columns
pmid_col = None
variant_col = None
lqt_col = None

for col in df.columns:
    col_lower = str(col).lower()
    if 'pmid' in col_lower or 'pubmed' in col_lower:
        pmid_col = col
    if 'variant' in col_lower or 'mutation' in col_lower:
        if variant_col is None:
            variant_col = col
    if 'lqt' in col_lower or 'affect' in col_lower or 'phenotype' in col_lower:
        lqt_col = col

print(f"\nIdentified columns:")
print(f"  PMID column: {pmid_col}")
print(f"  Variant column: {variant_col}")
print(f"  LQT column: {lqt_col}")

# Show sample data
print(f"\n=== Sample Data (first 5 rows) ===")
print(df.head())

# Get unique PMIDs
if pmid_col:
    # PMIDs might be in multiple columns - let's check for any column containing PMIDs
    pmid_values = df[pmid_col].dropna().unique()
    print(f"\n=== PMID Statistics ===")
    print(f"Unique PMIDs in '{pmid_col}': {len(pmid_values)}")
    print(f"Sample PMIDs: {list(pmid_values[:10])}")
else:
    # Search for PMIDs in all columns
    print("\nSearching for PMID-like columns...")
    for col in df.columns:
        sample_vals = df[col].dropna().head()
        print(f"  {col}: {list(sample_vals)[:3]}")

# Count unique variants
if variant_col:
    variants = df[variant_col].dropna().unique()
    print(f"\n=== Variant Statistics ===")
    print(f"Unique variants: {len(variants)}")
    print(f"Sample variants: {list(variants[:10])}")

# Count LQT-affected
if lqt_col:
    print(f"\n=== LQT Status ===")
    print(df[lqt_col].value_counts())

# Get all downloaded PMIDs
download_dir = "/mnt/temp2/kronckbm/gvf_output/verified_downloads_20260208"
downloaded_pmids = set()
for pdf in glob.glob(os.path.join(download_dir, "*.pdf")):
    basename = os.path.basename(pdf).replace('.pdf', '')
    # Extract PMID (might have suffixes like _supp)
    pmid = re.split(r'[_-]', basename)[0]
    if pmid.isdigit():
        downloaded_pmids.add(int(pmid))

print(f"\n=== Downloaded Papers ===")
print(f"Total PDFs: {len(glob.glob(os.path.join(download_dir, '*.pdf')))}")
print(f"Unique PMIDs: {len(downloaded_pmids)}")

# Save analysis
analysis = {
    "total_rows": len(df),
    "columns": list(df.columns),
    "pmid_col": pmid_col,
    "variant_col": variant_col,
    "lqt_col": lqt_col,
    "downloaded_pmids": len(downloaded_pmids)
}

with open("gold_standard_analysis.json", "w") as f:
    json.dump(analysis, f, indent=2, default=str)
print("\nSaved analysis to gold_standard_analysis.json")
