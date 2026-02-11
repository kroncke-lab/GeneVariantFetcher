#!/usr/bin/env python3
"""Extract all unique variants from gold standard Excel."""

import pandas as pd
import json
import re

# Read the gold standard Excel
xls_path = "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
df = pd.read_excel(xls_path)

print("=== GOLD STANDARD COLUMN NAMES ===")
print(df.columns.tolist())
print()

print("=== FIRST 5 ROWS ===")
print(df.head(5).to_string())
print()

# Find the variant column (usually the main identifier)
# Look for Mutation or Variant columns
for col in df.columns:
    print(f"{col}: {df[col].iloc[0] if pd.notna(df[col].iloc[0]) else 'NaN'}")
print()

# Extract all unique variant names from likely columns
variant_columns = []
for col in df.columns:
    col_lower = col.lower()
    if any(
        term in col_lower for term in ["variant", "mut", "change", "notation", "allele"]
    ):
        variant_columns.append(col)

print(f"=== POTENTIAL VARIANT COLUMNS ===")
print(variant_columns)

# Get all unique values from these columns
all_variants = set()
for col in variant_columns:
    unique_vals = df[col].dropna().unique()
    for v in unique_vals:
        if isinstance(v, str) and len(v) > 0:
            all_variants.add(v.strip())

print(f"\n=== ALL UNIQUE VARIANTS (first 50) ===")
sorted_variants = sorted(all_variants)
for v in sorted_variants[:50]:
    print(f"  {v}")
print(f"\nTotal unique: {len(all_variants)}")

# Also get the "Mutation" column specifically if it exists
if "Mutation" in df.columns:
    mutations = df["Mutation"].dropna().unique()
    print(f"\n=== UNIQUE FROM 'Mutation' COLUMN ===")
    for m in sorted(mutations)[:30]:
        print(f"  {m}")
    print(f"Total: {len(mutations)}")

# Save all unique variants to JSON
output = {
    "columns_found": variant_columns,
    "total_unique_variants": len(all_variants),
    "variants": sorted(list(all_variants)),
}

with open("utils/gold_standard_variants.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\n=== SAVED TO utils/gold_standard_variants.json ===")
