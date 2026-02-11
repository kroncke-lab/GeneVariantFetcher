#!/usr/bin/env python3
"""Find PMIDs that need to be downloaded."""

import sys

sys.path.insert(0, "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher")

import sqlite3
import pandas as pd
import re


def normalize_pmid(p):
    if pd.isna(p):
        return None
    s = str(p).strip()
    s = re.sub(r"\.0$", "", s)
    s = re.sub(r"\D+", "", s)
    return s if s else None


# Load Excel PMIDs
excel_path = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
df = pd.read_excel(excel_path)
excel_pmids = set(normalize_pmid(p) for p in df["PMID"].unique() if normalize_pmid(p))

# Load SQLite PMIDs
conn = sqlite3.connect("/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db")
cursor = conn.cursor()
cursor.execute("SELECT DISTINCT pmid FROM papers")
sqlite_pmids = set(row[0] for row in cursor.fetchall() if row[0])
conn.close()

# Find missing PMIDs
missing = sorted(excel_pmids - sqlite_pmids, key=lambda x: int(x) if x.isdigit() else 0)

print(f"Excel PMIDs: {len(excel_pmids)}")
print(f"SQLite PMIDs: {len(sqlite_pmids)}")
print(f"Missing PMIDs: {len(missing)}")

# Save to file
out_path = "/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt"
with open(out_path, "w") as f:
    for pmid in missing:
        f.write(f"{pmid}\n")
print(f"\nSaved missing PMIDs to: {out_path}")

# Count variant entries per missing PMID (to prioritize)
print("\nTop 20 missing PMIDs by Excel variant count:")
missing_with_counts = []
for pmid in missing:
    count = len(df[df["PMID"].apply(normalize_pmid) == pmid])
    missing_with_counts.append((pmid, count))

missing_with_counts.sort(key=lambda x: -x[1])
total_missing_entries = sum(c for _, c in missing_with_counts)
print(f"Total variant entries in missing PMIDs: {total_missing_entries}")
for pmid, count in missing_with_counts[:20]:
    print(f"  PMID {pmid}: {count} variants")
