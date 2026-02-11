#!/usr/bin/env python3
"""Analyze PMID overlap between Excel and SQLite."""

import sys

sys.path.insert(0, "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher")

import sqlite3
import pandas as pd

# Load Excel PMIDs
excel_path = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
df = pd.read_excel(excel_path)
excel_pmids = set(
    str(p).strip() for p in df["PMID"].dropna().unique() if str(p).strip()
)

# Load SQLite PMIDs
conn = sqlite3.connect("/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db")
cursor = conn.cursor()
cursor.execute("SELECT DISTINCT pmid FROM papers")
sqlite_pmids = set(row[0] for row in cursor.fetchall())

# Load variants-only PMIDs (papers with actual variant data)
cursor.execute("SELECT DISTINCT pmid FROM variant_papers")
variants_pmids = set(row[0] for row in cursor.fetchall())

conn.close()

print(f"Excel unique PMIDs: {len(excel_pmids)}")
print(f"SQLite papers (total): {len(sqlite_pmids)}")
print(f"SQLite papers (with variants): {len(variants_pmids)}")

overlap = excel_pmids & variants_pmids
print(f"\nPMID overlap (SQLite has variants for these Excel PMIDs): {len(overlap)}")
print(f"Coverage: {len(overlap) / len(excel_pmids) * 100:.1f}% of Excel PMIDs")

missing_from_sqlite = excel_pmids - sqlite_pmids
print(f"\nExcel PMIDs NOT in SQLite at all: {len(missing_from_sqlite)}")

# PMIDs that are in SQLite but have no variants
in_sqlite_no_variants = sqlite_pmids - variants_pmids
print(f"SQLite papers with NO variants extracted: {len(in_sqlite_no_variants)}")

# Show some missing PMIDs
print("\nSample of Excel PMIDs NOT in SQLite:")
for pmid in sorted(missing_from_sqlite)[:20]:
    print(f"  {pmid}")

# Show variant count per matching PMID
conn = sqlite3.connect("/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db")
cursor = conn.cursor()

print("\n--- PMID-level variant comparison ---")
for pmid in sorted(overlap)[:10]:
    cursor.execute(
        """
        SELECT COUNT(*) FROM variant_papers WHERE pmid = ?
    """,
        (pmid,),
    )
    sqlite_count = cursor.fetchone()[0]

    excel_count = len(df[df["PMID"].astype(str).str.strip() == pmid])
    print(f"PMID {pmid}: Excel={excel_count}, SQLite={sqlite_count}")

conn.close()
