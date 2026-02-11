#!/usr/bin/env python3
"""
Comprehensive variant matching using aliases.

This script:
1. Loads all extracted variants (from novel + matched)
2. Normalizes each using the alias dictionary
3. Matches against gold standard
4. Reports improved recall
"""

import json
import re
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Amino acid mappings
AA_3_TO_1 = {
    "Ala": "A",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y",
    "Ter": "*",
    "Stop": "*",
    "Xaa": "X",
}


def normalize_to_canonical(variant, aliases):
    """Normalize a variant to canonical form using alias dictionary."""
    if not variant:
        return None

    v = str(variant).strip()
    original = v

    # Try direct lookup first
    v_upper = v.upper()
    if v_upper in aliases:
        return aliases[v_upper]

    # Clean up common issues
    # Remove trailing asterisks (DMS artifacts): A561V* -> A561V
    if v.endswith("*") and len(v) >= 2 and v[-2].isalpha():
        v = v[:-1]

    # Remove trailing annotations: R176W(het) -> R176W
    v = re.sub(r"\([^)]*\)$", "", v)
    v = re.sub(r"/[Ww][Tt]$", "", v)
    v = v.strip()

    v_upper = v.upper()
    if v_upper in aliases:
        return aliases[v_upper]

    # Try without p. prefix
    if v_upper.startswith("P."):
        stripped = v_upper[2:]
        if stripped in aliases:
            return aliases[stripped]

    # Try three-letter to single-letter conversion
    # Pattern: p.Ala561Val or Ala561Val
    m = re.match(r"^(?:P\.)?([A-Z][A-Z][A-Z])(\d+)([A-Z][A-Z][A-Z])$", v_upper)
    if m:
        ref_3 = m.group(1).capitalize()
        alt_3 = m.group(3).capitalize()
        ref_1 = AA_3_TO_1.get(ref_3, "")
        alt_1 = AA_3_TO_1.get(alt_3, "")
        if ref_1 and alt_1:
            single = f"{ref_1}{m.group(2)}{alt_1}"
            if single.upper() in aliases:
                return aliases[single.upper()]
            return single

    # Handle frameshift: p.Ala193fs*46, Gly262Alafs*98
    m = re.match(
        r"^(?:P\.)?([A-Z][A-Z][A-Z])(\d+)(?:[A-Z][A-Z][A-Z])?FS[\*X]?\d*$", v_upper
    )
    if m:
        ref_3 = m.group(1).capitalize()
        ref_1 = AA_3_TO_1.get(ref_3, "")
        if ref_1:
            canonical = f"{ref_1}{m.group(2)}fsX"
            if canonical.upper() in aliases:
                return aliases[canonical.upper()]
            return canonical.upper()

    # Handle single-letter frameshift: A193fsX10, A193fs*46
    m = re.match(r"^(?:P\.)?([A-Z])(\d+)FS[\*X]?\d*$", v_upper)
    if m:
        canonical = f"{m.group(1)}{m.group(2)}FSX"
        if canonical in aliases:
            return aliases[canonical]
        return f"{m.group(1)}{m.group(2)}fsX"

    # Handle stop codons: p.Arg864Ter, R864*, R864X
    m = re.match(r"^(?:P\.)?([A-Z][A-Z][A-Z])(\d+)(TER|\*|X|STOP)$", v_upper)
    if m:
        ref_3 = m.group(1).capitalize()
        ref_1 = AA_3_TO_1.get(ref_3, "")
        if ref_1:
            canonical = f"{ref_1}{m.group(2)}X"
            if canonical.upper() in aliases:
                return aliases[canonical.upper()]
            return canonical

    m = re.match(r"^(?:P\.)?([A-Z])(\d+)(TER|\*|X|STOP)$", v_upper)
    if m:
        canonical = f"{m.group(1)}{m.group(2)}X"
        if canonical in aliases:
            return aliases[canonical]
        return canonical

    # Single-letter format: A561V
    m = re.match(r"^(?:P\.)?([A-Z])(\d+)([A-Z])$", v_upper)
    if m:
        canonical = f"{m.group(1)}{m.group(2)}{m.group(3)}"
        if canonical in aliases:
            return aliases[canonical]
        return canonical

    # Return uppercase original if no normalization possible
    return v_upper


def main():
    # Load alias dictionary
    with open("utils/kcnh2_variant_aliases.json") as f:
        alias_data = json.load(f)
        aliases = alias_data["aliases"]

    # Load recall data
    with open(
        "/mnt/temp2/kronckbm/gvf_output/per_variant_recall_20260210_v3.json"
    ) as f:
        recall_data = json.load(f)

    # Load the gold standard Excel to get ALL extracted variants
    xls_path = "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
    df = pd.read_excel(xls_path)

    # Get gold standard variants (EXCLUDE?=0)
    include_mask = df["EXCLUDE?"] == 0
    gold_variants = set(df[include_mask]["Variant"].dropna().unique())
    gold_all = set(df["Variant"].dropna().unique())

    print("=" * 70)
    print("COMPREHENSIVE VARIANT MATCHING WITH ALIASES")
    print("=" * 70)
    print(f"Gold standard (include only): {len(gold_variants)}")
    print(f"Gold standard (all): {len(gold_all)}")
    print(f"Aliases in dictionary: {len(aliases)}")

    # Build canonical gold standard lookup
    gold_canonical = {}  # canonical -> original
    for v in gold_variants:
        canonical = normalize_to_canonical(v, aliases)
        gold_canonical[canonical] = v

    print(f"Unique canonical gold standard: {len(gold_canonical)}")

    # Get all GVF extracted variants
    matched = set(recall_data["variants"]["matched"])
    novel = set(recall_data["variants"]["novel"])
    all_extracted = matched | novel

    print(f"\nGVF extracted variants: {len(all_extracted)}")
    print(f"  Already matched: {len(matched)}")
    print(f"  Novel (unmatched): {len(novel)}")

    # Now, for each extracted variant, normalize and try to match
    new_matches = []
    old_matches = []
    still_unmatched = []

    for ext in all_extracted:
        canonical = normalize_to_canonical(ext, aliases)

        if canonical in gold_canonical:
            if ext in matched:
                old_matches.append((ext, gold_canonical[canonical], "original"))
            else:
                new_matches.append((ext, gold_canonical[canonical], canonical))
        else:
            # Try case-insensitive match
            found = False
            for gc, gv in gold_canonical.items():
                if gc.upper() == canonical.upper():
                    if ext in matched:
                        old_matches.append((ext, gv, "case-insensitive"))
                    else:
                        new_matches.append((ext, gv, canonical))
                    found = True
                    break
            if not found:
                still_unmatched.append((ext, canonical))

    print(f"\n" + "=" * 70)
    print("MATCHING RESULTS")
    print("=" * 70)
    print(f"Original matches (already working): {len(old_matches)}")
    print(f"NEW matches via normalization: {len(new_matches)}")
    print(f"Still unmatched: {len(still_unmatched)}")

    # Show new matches
    print(f"\n--- NEW MATCHES ({len(new_matches)}) ---")
    for ext, gold, canonical in sorted(new_matches, key=lambda x: x[1]):
        print(f"  {ext:25} -> {gold:15} (canonical: {canonical})")

    # Calculate new recall
    matched_gold = set()
    for ext, gold, _ in old_matches + new_matches:
        matched_gold.add(gold)

    original_recall = len(matched) / len(gold_variants) * 100
    new_recall = len(matched_gold) / len(gold_variants) * 100

    print(f"\n" + "=" * 70)
    print("RECALL SUMMARY")
    print("=" * 70)
    print(
        f"Original unique variant recall: {len(matched)}/{len(gold_variants)} = {original_recall:.1f}%"
    )
    print(
        f"New unique variant recall: {len(matched_gold)}/{len(gold_variants)} = {new_recall:.1f}%"
    )
    print(
        f"Improvement: +{len(new_matches)} variants (+{(new_recall - original_recall):.1f}%)"
    )

    # Analyze what's still missed
    print(f"\n" + "=" * 70)
    print("ANALYSIS OF STILL MISSED VARIANTS")
    print("=" * 70)

    missed_gold = gold_variants - matched_gold
    print(f"Gold standard variants still not found: {len(missed_gold)}")

    # Categorize
    cdna = [v for v in missed_gold if v.startswith("c.")]
    ivs = [v for v in missed_gold if v.upper().startswith("IVS")]
    exon = [v for v in missed_gold if "EXON" in v.upper() or "EX" in v.upper()]
    chrom = [
        v for v in missed_gold if "DEL" in v.upper() and ("Q" in v.upper() or "(" in v)
    ]
    protein = [v for v in missed_gold if v not in cdna + ivs + exon + chrom]

    print(f"  cDNA notation: {len(cdna)}")
    if cdna:
        print(f"    Examples: {cdna[:5]}")
    print(f"  IVS/splice: {len(ivs)}")
    if ivs:
        print(f"    Examples: {ivs[:5]}")
    print(f"  Exon deletions/dups: {len(exon)}")
    if exon:
        print(f"    Examples: {exon[:5]}")
    print(f"  Chromosomal: {len(chrom)}")
    if chrom:
        print(f"    Examples: {chrom[:5]}")
    print(f"  Protein variants: {len(protein)}")
    if protein:
        print(f"    Examples: {sorted(protein)[:20]}")

    # Save results for use in normalizer
    output = {
        "original_matched": len(matched),
        "new_matches": len(new_matches),
        "total_matched": len(matched_gold),
        "gold_standard_count": len(gold_variants),
        "original_recall": original_recall,
        "new_recall": new_recall,
        "new_match_details": [
            {"extracted": e, "gold": g, "canonical": c} for e, g, c in new_matches
        ],
        "still_missed": sorted(list(missed_gold)),
    }

    with open("utils/matching_improvement_analysis.json", "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n=== Saved analysis to utils/matching_improvement_analysis.json ===")


if __name__ == "__main__":
    main()
