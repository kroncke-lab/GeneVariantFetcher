#!/usr/bin/env python3
"""Check if specific expected variants are in the extracted set."""

import sys
from pathlib import Path
import json

repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

from utils.variant_normalizer import normalize_variant

# Load recall data
with open("/mnt/temp2/kronckbm/gvf_output/per_variant_recall_20260210_v3.json") as f:
    recall_data = json.load(f)

matched = set(recall_data["variants"]["matched"])
missed = set(recall_data["variants"]["missed"])
novel = set(recall_data["variants"]["novel"])
all_extracted = matched | novel

# Variants we expect might be in novel but should match gold
expected_extractions = [
    # From unmatched analysis - "normalization matches"
    "p.Gly262Alafs*98",
    "p.Arg176Trp",
    "p.Val644Phe",
    "p.Glu58Lys",
    "p.Gly903Arg",
    "P72Q",
    "T613M",
    "A561V",
    "W1001X",
    "S660L",
    "R863X",
]

print("=== CHECKING EXPECTED VARIANTS ===")
print(f"Extracted variants total: {len(all_extracted)}")
print(f"  Matched: {len(matched)}")
print(f"  Novel: {len(novel)}")
print()

for v in expected_extractions:
    v_upper = v.upper()
    normalized = normalize_variant(v, "KCNH2")
    in_matched = v in matched or v_upper in [x.upper() for x in matched]
    in_novel = v in novel or v_upper in [x.upper() for x in novel]
    in_missed = normalized in missed or normalized.upper() in [
        x.upper() for x in missed
    ]

    status = "MATCHED" if in_matched else ("IN NOVEL" if in_novel else "NOT EXTRACTED")
    gold_status = (
        "IN GOLD(missed)"
        if in_missed
        else ("IN GOLD(matched)" if normalized in matched else "NOT IN GOLD")
    )

    print(f"  {v:<25} -> {normalized:<15} {status:<15} {gold_status}")

# Check if the extracted form of G262fsX is in novel
print()
print("=== SEARCHING FOR SPECIFIC FORMS ===")
search_terms = [
    "262",
    "G262",
    "GLY262",
    "Gly262",
    "176",
    "R176",
    "Arg176",
    "ARG176",
    "660",
    "S660",
]

for term in search_terms:
    found = [v for v in novel if term.upper() in v.upper()][:5]
    if found:
        print(f"  '{term}' found in novel: {found}")
    else:
        print(f"  '{term}' NOT found in novel")

# Check what forms of the missed variants exist in all_extracted
print()
print("=== CROSS-CHECKING MISSED VS EXTRACTED ===")
sample_missed = ["G262fsX", "W1001X", "S660L", "P72Q", "R863X"]

for m in sample_missed:
    m_upper = m.upper()
    # Look for any extracted variant that normalizes to this
    matching_extracted = []
    for ext in all_extracted:
        ext_norm = normalize_variant(ext, "KCNH2")
        if ext_norm.upper() == m_upper:
            matching_extracted.append(ext)

    if matching_extracted:
        print(f"  {m}: FOUND -> {matching_extracted}")
    else:
        print(f"  {m}: NOT FOUND in extracted")
