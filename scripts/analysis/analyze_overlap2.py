#!/usr/bin/env python3
"""Analyze PMID overlap between Excel and SQLite (fixed PMID parsing)."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher')

import sqlite3
import pandas as pd
import re

def normalize_pmid(p):
    """Normalize PMID - strip decimals and whitespace."""
    if pd.isna(p):
        return None
    s = str(p).strip()
    # Remove .0 if present
    s = re.sub(r'\.0$', '', s)
    # Remove any non-digit
    s = re.sub(r'\D+', '', s)
    return s if s else None

# Load Excel PMIDs
excel_path = '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls'
df = pd.read_excel(excel_path)
excel_pmids = set(normalize_pmid(p) for p in df['PMID'].unique() if normalize_pmid(p))

# Load SQLite PMIDs
conn = sqlite3.connect('/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db')
cursor = conn.cursor()

cursor.execute("SELECT DISTINCT pmid FROM papers")
sqlite_pmids = set(row[0] for row in cursor.fetchall() if row[0])

cursor.execute("SELECT DISTINCT pmid FROM variant_papers")
variants_pmids = set(row[0] for row in cursor.fetchall() if row[0])

print(f"Excel unique PMIDs: {len(excel_pmids)}")
print(f"SQLite papers (total): {len(sqlite_pmids)}")
print(f"SQLite papers (with variants): {len(variants_pmids)}")

overlap = excel_pmids & variants_pmids
print(f"\nPMID overlap: {len(overlap)}")
print(f"Coverage: {len(overlap)/len(excel_pmids)*100:.1f}% of Excel PMIDs")

missing_from_sqlite = excel_pmids - sqlite_pmids
no_extraction = excel_pmids - variants_pmids
print(f"Excel PMIDs NOT in SQLite at all: {len(missing_from_sqlite)}")
print(f"Excel PMIDs without variant extraction: {len(no_extraction)}")

# For overlapping PMIDs, count variants
total_excel_in_overlap = 0
total_sqlite_in_overlap = 0

print("\n--- Sample PMID-level variant comparison ---")
for pmid in sorted(overlap)[:15]:
    cursor.execute("SELECT COUNT(*) FROM variant_papers WHERE pmid = ?", (pmid,))
    sqlite_count = cursor.fetchone()[0]
    
    excel_count = len(df[df['PMID'].apply(normalize_pmid) == pmid])
    total_excel_in_overlap += excel_count
    total_sqlite_in_overlap += sqlite_count
    print(f"PMID {pmid}: Excel={excel_count}, SQLite={sqlite_count}")

print(f"\n--- Summary for overlapping PMIDs ---")
print(f"Total Excel entries in overlap: {total_excel_in_overlap}")
print(f"Total SQLite variants in overlap: {total_sqlite_in_overlap}")

# Full count
full_excel = 0
full_sqlite = 0
for pmid in overlap:
    cursor.execute("SELECT COUNT(*) FROM variant_papers WHERE pmid = ?", (pmid,))
    full_sqlite += cursor.fetchone()[0]
    full_excel += len(df[df['PMID'].apply(normalize_pmid) == pmid])

print(f"\nFull overlap analysis:")
print(f"Total Excel entries for overlapping PMIDs: {full_excel}")
print(f"Total SQLite variants for overlapping PMIDs: {full_sqlite}")

conn.close()
