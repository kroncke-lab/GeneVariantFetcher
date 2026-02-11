#!/usr/bin/env python3
"""Comprehensive GVF Recall Analysis."""

import pandas as pd
import os
import glob
import re
import json
from collections import defaultdict

# Read gold standard Excel file
xls_path = "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
df = pd.read_excel(xls_path)

print("=" * 70)
print("GVF RECALL ANALYSIS - COMPREHENSIVE REPORT")
print("=" * 70)

# Get all downloaded PMIDs and check for supplements
download_dir = "/mnt/temp2/kronckbm/gvf_output/verified_downloads_20260208"
downloaded_files = glob.glob(os.path.join(download_dir, "*.pdf"))

pmid_to_files = defaultdict(list)
for pdf in downloaded_files:
    basename = os.path.basename(pdf)
    # Check if it's a supplement
    is_supp = "supp" in basename.lower() or "_s" in basename.lower()
    # Extract PMID
    pmid_match = re.match(r"^(\d+)", basename)
    if pmid_match:
        pmid = int(pmid_match.group(1))
        pmid_to_files[pmid].append({"file": basename, "is_supplement": is_supp})

# Analyze downloads
full_text_only = 0
full_text_plus_supp = 0
for pmid, files in pmid_to_files.items():
    has_main = any(not f["is_supplement"] for f in files)
    has_supp = any(f["is_supplement"] for f in files)
    if has_main and has_supp:
        full_text_plus_supp += 1
    elif has_main:
        full_text_only += 1

# Get gold standard PMIDs
gold_pmids = set()
df_with_pmid = df.dropna(subset=["PMID"])
for pmid in df_with_pmid["PMID"].unique():
    try:
        gold_pmids.add(int(float(pmid)))
    except (ValueError, TypeError):
        pass

downloaded_pmids = set(pmid_to_files.keys())
missing_pmids = gold_pmids - downloaded_pmids
extra_pmids = downloaded_pmids - gold_pmids

print(f"\n{'=' * 70}")
print("PHASE 1: PAPER RECALL")
print(f"{'=' * 70}")
print(f"Gold Standard PMIDs:     {len(gold_pmids)}")
print(f"Downloaded PMIDs:        {len(downloaded_pmids)}")
print(f"  - Full text only:      {full_text_only}")
print(f"  - Full text + supps:   {full_text_plus_supp}")
print(f"Missing PMIDs:           {len(missing_pmids)}")
print(f"Extra PMIDs (not in GS): {len(extra_pmids)}")
print(
    f"\nPaper Recall Rate:       {len(downloaded_pmids & gold_pmids)}/{len(gold_pmids)} = {100 * len(downloaded_pmids & gold_pmids) / len(gold_pmids):.1f}%"
)

# Variant analysis - which variants are in downloaded papers
print(f"\n{'=' * 70}")
print("PHASE 2: VARIANT RECALL (Paper-Level)")
print(f"{'=' * 70}")

# Group variants by PMID
df_with_pmid["PMID_int"] = df_with_pmid["PMID"].apply(
    lambda x: int(float(x)) if pd.notna(x) else None
)

variants_in_downloaded = set()
variants_in_missing = set()
total_variants = set(df["Variant"].dropna().unique())

for _, row in df_with_pmid.iterrows():
    pmid = row["PMID_int"]
    variant = row["Variant"]
    if pd.isna(variant):
        continue
    if pmid in downloaded_pmids:
        variants_in_downloaded.add(variant)
    else:
        variants_in_missing.add(variant)

# Some variants appear in both downloaded and missing papers
variants_only_in_missing = variants_in_missing - variants_in_downloaded

print(f"Total unique variants (gold standard):   {len(total_variants)}")
print(f"Variants in downloaded papers:           {len(variants_in_downloaded)}")
print(f"Variants only in missing papers:         {len(variants_only_in_missing)}")
print(
    f"\nVariant Coverage (by paper):             {len(variants_in_downloaded)}/{len(total_variants)} = {100 * len(variants_in_downloaded) / len(total_variants):.1f}%"
)

# Carrier analysis
print(f"\n{'=' * 70}")
print("PHASE 3: CARRIER RECALL")
print(f"{'=' * 70}")

total_carriers = len(df)
carriers_in_downloaded = len(
    df_with_pmid[df_with_pmid["PMID_int"].isin(downloaded_pmids)]
)
carriers_in_missing = total_carriers - carriers_in_downloaded

# Count carriers with no PMID
carriers_no_pmid = len(df) - len(df_with_pmid)

print(f"Total carriers (gold standard):          {total_carriers}")
print(f"Carriers in downloaded papers:           {carriers_in_downloaded}")
print(f"Carriers in missing papers:              {carriers_in_missing}")
print(f"Carriers with no PMID:                   {carriers_no_pmid}")
print(
    f"\nCarrier Coverage:                        {carriers_in_downloaded}/{total_carriers} = {100 * carriers_in_downloaded / total_carriers:.1f}%"
)

# LQT-Affected analysis
print(f"\n{'=' * 70}")
print("PHASE 4: LQT-AFFECTED RECALL")
print(f"{'=' * 70}")

# LQT >= 1 means affected
df["LQT_affected"] = df["LQT"].apply(lambda x: x >= 1 if pd.notna(x) else False)
total_lqt_affected = df["LQT_affected"].sum()

df_with_pmid["LQT_affected"] = df_with_pmid["LQT"].apply(
    lambda x: x >= 1 if pd.notna(x) else False
)
lqt_in_downloaded = df_with_pmid[df_with_pmid["PMID_int"].isin(downloaded_pmids)][
    "LQT_affected"
].sum()

print(f"Total LQT-affected (gold standard):      {int(total_lqt_affected)}")
print(f"LQT-affected in downloaded papers:       {int(lqt_in_downloaded)}")
print(
    f"\nLQT-Affected Coverage:                   {int(lqt_in_downloaded)}/{int(total_lqt_affected)} = {100 * lqt_in_downloaded / total_lqt_affected:.1f}%"
)

# Summary table
print(f"\n{'=' * 70}")
print("SUMMARY TABLE")
print(f"{'=' * 70}")
print(f"| {'Metric':<30} | {'Current':<15} | {'Coverage':<10} |")
print(f"|{'-' * 32}|{'-' * 17}|{'-' * 12}|")
print(
    f"| {'Paper Recall':<30} | {len(downloaded_pmids & gold_pmids)}/{len(gold_pmids):<13} | {100 * len(downloaded_pmids & gold_pmids) / len(gold_pmids):.1f}%{' ' * 5} |"
)
print(f"| {'  - Full text only':<30} | {full_text_only:<15} | {' ' * 10} |")
print(
    f"| {'  - Full text + supplements':<30} | {full_text_plus_supp:<15} | {' ' * 10} |"
)
print(
    f"| {'Variant Coverage (by paper)':<30} | {len(variants_in_downloaded)}/{len(total_variants):<11} | {100 * len(variants_in_downloaded) / len(total_variants):.1f}%{' ' * 5} |"
)
print(
    f"| {'Carrier Coverage':<30} | {carriers_in_downloaded}/{total_carriers:<11} | {100 * carriers_in_downloaded / total_carriers:.1f}%{' ' * 5} |"
)
print(
    f"| {'LQT-Affected Coverage':<30} | {int(lqt_in_downloaded)}/{int(total_lqt_affected):<11} | {100 * lqt_in_downloaded / total_lqt_affected:.1f}%{' ' * 5} |"
)

# List missing PMIDs
print(f"\n{'=' * 70}")
print(f"MISSING PMIDs ({len(missing_pmids)} papers)")
print(f"{'=' * 70}")
missing_sorted = sorted(missing_pmids)
for i, pmid in enumerate(missing_sorted):
    if i < 20:
        # Get count of carriers for this PMID
        carrier_count = len(df_with_pmid[df_with_pmid["PMID_int"] == pmid])
        print(f"  {pmid} - {carrier_count} carriers")
    elif i == 20:
        print(f"  ... and {len(missing_sorted) - 20} more")
        break

# Save detailed results
results = {
    "paper_recall": {
        "gold_standard": len(gold_pmids),
        "downloaded": len(downloaded_pmids),
        "overlap": len(downloaded_pmids & gold_pmids),
        "missing": len(missing_pmids),
        "full_text_only": full_text_only,
        "full_text_plus_supp": full_text_plus_supp,
    },
    "variant_recall": {
        "total_variants": len(total_variants),
        "variants_in_downloaded": len(variants_in_downloaded),
        "variants_only_in_missing": len(variants_only_in_missing),
    },
    "carrier_recall": {
        "total_carriers": total_carriers,
        "carriers_in_downloaded": carriers_in_downloaded,
        "carriers_no_pmid": carriers_no_pmid,
    },
    "lqt_recall": {
        "total_affected": int(total_lqt_affected),
        "affected_in_downloaded": int(lqt_in_downloaded),
    },
    "missing_pmids": sorted(missing_pmids),
}

with open("recall_analysis_results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nDetailed results saved to recall_analysis_results.json")
