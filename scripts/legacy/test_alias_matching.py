#!/usr/bin/env python3
"""
Test alias matching to see how many missed variants can be recovered.
"""

import json
import re
from pathlib import Path

# Load the alias dictionary
with open("utils/kcnh2_variant_aliases.json") as f:
    data = json.load(f)
    aliases = data["aliases"]

# Load the recall data
with open("/mnt/temp2/kronckbm/gvf_output/per_variant_recall_20260210_v3.json") as f:
    recall_data = json.load(f)

# Get the gold standard variants (from the recall data)
gold_matched = set(recall_data["variants"]["matched"])
gold_missed = set(recall_data["variants"]["missed"])
gold_all = gold_matched | gold_missed

# Get the GVF extracted variants (novel = extracted but not matched)
gvf_novel = recall_data["variants"]["novel"]

print("=" * 60)
print("ALIAS MATCHING TEST")
print("=" * 60)
print(f"Gold standard variants: {len(gold_all)}")
print(f"  Currently matched: {len(gold_matched)}")
print(f"  Currently missed: {len(gold_missed)}")
print(f"GVF novel variants: {len(gvf_novel)}")
print(f"Aliases in dictionary: {len(aliases)}")

# Now, for each GVF novel variant, see if it maps to a missed gold standard variant
new_matches = []
still_missed = []

# First, build a set of canonical forms for missed variants
missed_canonicals = set()
for v in gold_missed:
    v_upper = v.upper()
    if v_upper in aliases:
        missed_canonicals.add(aliases[v_upper])
    else:
        missed_canonicals.add(v_upper)

print(f"\nMissed variants canonical forms: {len(missed_canonicals)}")
print(f"  Sample: {list(missed_canonicals)[:10]}")

# Now check each novel GVF variant
for novel in gvf_novel:
    novel_upper = novel.upper()

    # Check if this novel variant maps to a canonical form
    if novel_upper in aliases:
        canonical = aliases[novel_upper]
        # Check if this canonical is in the missed set
        if canonical in missed_canonicals:
            new_matches.append((novel, canonical))
        elif canonical in gold_all or canonical.upper() in [
            g.upper() for g in gold_all
        ]:
            # Already matched
            pass
        else:
            # Novel - not in gold at all
            pass
    else:
        # Try alternate normalization
        # Strip p. prefix and try again
        if novel_upper.startswith("P."):
            stripped = novel_upper[2:]
            if stripped in aliases:
                canonical = aliases[stripped]
                if canonical in missed_canonicals:
                    new_matches.append((novel, canonical))

print("\n" + "=" * 60)
print("POTENTIAL NEW MATCHES VIA ALIASING")
print("=" * 60)
print(f"Found {len(new_matches)} potential new matches")

for novel, canonical in new_matches[:30]:
    print(f"  {novel} -> {canonical}")

if len(new_matches) > 30:
    print(f"  ... and {len(new_matches) - 30} more")

# Calculate new potential recall
new_matched_count = len(gold_matched) + len(new_matches)
new_recall = new_matched_count / len(gold_all) * 100

print("\n" + "=" * 60)
print("RECALL IMPROVEMENT")
print("=" * 60)
print(
    f"Original matched: {len(gold_matched)} / {len(gold_all)} = {len(gold_matched) / len(gold_all) * 100:.1f}%"
)
print(f"New potential matches: {len(new_matches)}")
print(f"New total matched: {new_matched_count} / {len(gold_all)} = {new_recall:.1f}%")
print(
    f"Improvement: +{len(new_matches)} variants (+{len(new_matches) / len(gold_all) * 100:.1f}%)"
)

# Also check: how many missed variants have aliases that weren't found?
print("\n" + "=" * 60)
print("REMAINING MISSED ANALYSIS")
print("=" * 60)

matched_canonicals = set(c for _, c in new_matches)
remaining_missed = [
    v for v in gold_missed if v.upper() not in [c.upper() for c in matched_canonicals]
]

print(f"Still missed after aliasing: {len(remaining_missed)}")

# Categorize remaining missed
cdna_missed = [v for v in remaining_missed if v.startswith("c.")]
ivs_missed = [v for v in remaining_missed if v.upper().startswith("IVS")]
exon_missed = [
    v
    for v in remaining_missed
    if "EXON" in v.upper() or "EX" in v.upper() or "DEL" in v.upper()
]
other_missed = [
    v for v in remaining_missed if v not in cdna_missed + ivs_missed + exon_missed
]

print(f"  cDNA notation: {len(cdna_missed)}")
print(f"  IVS/splice: {len(ivs_missed)}")
print(f"  Exon/large del: {len(exon_missed)}")
print(f"  Other protein: {len(other_missed)}")

print("\nOther missed (sample):")
for v in sorted(other_missed)[:20]:
    print(f"    {v}")
