#!/usr/bin/env python3
"""
Recalculate variant recall using the updated normalizer with comprehensive aliases.

This script:
1. Loads gold standard variants from Excel
2. Loads all raw extracted variants from JSON files in the output directory
3. Normalizes both gold standard and extracted variants using the updated normalizer
4. Calculates new recall metrics
"""

import sys
from pathlib import Path
import json
import pandas as pd
from collections import defaultdict

# Add the repo root to path
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

from utils.variant_normalizer import normalize_variant

def load_gold_standard():
    """Load gold standard variants from Excel."""
    xls_path = repo_root / "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
    df = pd.read_excel(xls_path)
    
    # Filter to included variants
    include_mask = df['EXCLUDE?'] == 0
    gold_df = df[include_mask]
    
    # Get unique variants and normalize
    gold_variants = {}  # normalized -> original
    for v in gold_df['Variant'].dropna().unique():
        v_str = str(v).strip()
        normalized = normalize_variant(v_str, 'KCNH2')
        gold_variants[normalized] = v_str
    
    return gold_variants, gold_df

def load_extracted_variants():
    """Load all extracted variants from the GVF output."""
    output_dir = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2/20260128_210249")
    
    all_extracted = set()
    
    # Load per-file extractions
    extraction_files = list(output_dir.glob("**/extractions_*.json"))
    for f in extraction_files:
        try:
            with open(f) as fp:
                data = json.load(fp)
                if 'variants' in data:
                    for v in data['variants']:
                        if isinstance(v, dict):
                            # Try different keys for variant name
                            for key in ['protein_change', 'protein_notation', 'variant', 'name']:
                                if key in v and v[key]:
                                    all_extracted.add(str(v[key]).strip())
                        elif isinstance(v, str):
                            all_extracted.add(v.strip())
        except Exception as e:
            print(f"Warning: Could not load {f}: {e}")
    
    # Also load aggregated results if available
    agg_files = list(output_dir.glob("**/aggregated_*.json"))
    for f in agg_files:
        try:
            with open(f) as fp:
                data = json.load(fp)
                for v in data.get('variants', []):
                    if isinstance(v, dict):
                        for key in ['protein_change', 'protein_notation', 'variant', 'name']:
                            if key in v and v[key]:
                                all_extracted.add(str(v[key]).strip())
                    elif isinstance(v, str):
                        all_extracted.add(v.strip())
        except Exception as e:
            print(f"Warning: Could not load {f}: {e}")
    
    return all_extracted

def load_recall_json():
    """Load the existing recall JSON with already-matched variants."""
    recall_path = Path("/mnt/temp2/kronckbm/gvf_output/per_variant_recall_20260210_v3.json")
    with open(recall_path) as f:
        return json.load(f)

def main():
    print("=" * 70)
    print("RECALCULATING RECALL WITH UPDATED NORMALIZER")
    print("=" * 70)
    
    # Load gold standard
    gold_variants, gold_df = load_gold_standard()
    print(f"Gold standard unique variants: {len(gold_variants)}")
    
    # Load existing recall data
    recall_data = load_recall_json()
    matched_original = set(recall_data['variants']['matched'])
    novel_original = set(recall_data['variants']['novel'])
    
    # All extracted = matched + novel
    all_extracted = matched_original | novel_original
    print(f"Total extracted variants: {len(all_extracted)}")
    print(f"  Previously matched: {len(matched_original)}")
    print(f"  Previously novel: {len(novel_original)}")
    
    # Normalize all extracted variants
    extracted_normalized = {}  # normalized -> [originals]
    for v in all_extracted:
        norm = normalize_variant(v, 'KCNH2')
        if norm not in extracted_normalized:
            extracted_normalized[norm] = []
        extracted_normalized[norm].append(v)
    
    print(f"Unique normalized extracted: {len(extracted_normalized)}")
    
    # Find matches
    new_matches = []
    still_missed = []
    
    for gold_norm, gold_orig in gold_variants.items():
        if gold_norm in extracted_normalized:
            new_matches.append({
                'gold': gold_orig,
                'gold_normalized': gold_norm,
                'extracted_forms': extracted_normalized[gold_norm]
            })
        else:
            # Try case-insensitive
            found = False
            for ext_norm, ext_origs in extracted_normalized.items():
                if ext_norm.upper() == gold_norm.upper():
                    new_matches.append({
                        'gold': gold_orig,
                        'gold_normalized': gold_norm,
                        'extracted_forms': ext_origs
                    })
                    found = True
                    break
            if not found:
                still_missed.append(gold_orig)
    
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Original recall: {len(matched_original)}/{len(gold_variants)} = {len(matched_original)/len(gold_variants)*100:.1f}%")
    print(f"New matches with updated normalizer: {len(new_matches)}")
    print(f"Still missed: {len(still_missed)}")
    print(f"New recall: {len(new_matches)}/{len(gold_variants)} = {len(new_matches)/len(gold_variants)*100:.1f}%")
    
    # Show NEW matches (not in original matched set)
    new_recovered = []
    for m in new_matches:
        if m['gold'] not in matched_original:
            new_recovered.append(m)
    
    print()
    print(f"=== NEWLY RECOVERED VARIANTS ({len(new_recovered)}) ===")
    for m in sorted(new_recovered, key=lambda x: x['gold'])[:30]:
        print(f"  {m['gold']:<20} <- {m['extracted_forms'][0]}")
    
    if len(new_recovered) > 30:
        print(f"  ... and {len(new_recovered) - 30} more")
    
    # Categorize still missed
    print()
    print(f"=== STILL MISSED ({len(still_missed)}) ===")
    cdna = [v for v in still_missed if v.startswith('c.')]
    ivs = [v for v in still_missed if v.upper().startswith('IVS')]
    exon = [v for v in still_missed if 'EXON' in v.upper() or 'EX' in v.upper()]
    chrom = [v for v in still_missed if 'DEL' in v.upper() and ('Q' in v.upper() or '(' in v)]
    protein = [v for v in still_missed if v not in cdna + ivs + exon + chrom]
    
    print(f"  cDNA notation: {len(cdna)}")
    print(f"  IVS/splice: {len(ivs)}")
    print(f"  Exon deletions/dups: {len(exon)}")
    print(f"  Chromosomal: {len(chrom)}")
    print(f"  Protein variants not extracted: {len(protein)}")
    
    if protein:
        print(f"\n  Sample of missed protein variants:")
        for v in sorted(protein)[:30]:
            print(f"    {v}")
    
    # Save results
    output = {
        "date": "2026-02-10",
        "description": "Recall with updated normalizer and comprehensive aliases",
        "gold_standard_count": len(gold_variants),
        "original_matched": len(matched_original),
        "original_recall": len(matched_original) / len(gold_variants) * 100,
        "new_matched": len(new_matches),
        "new_recall": len(new_matches) / len(gold_variants) * 100,
        "improvement": len(new_matches) - len(matched_original),
        "newly_recovered": [m['gold'] for m in new_recovered],
        "still_missed": still_missed,
        "missed_by_category": {
            "cdna": cdna,
            "ivs": ivs,
            "exon": exon,
            "chromosomal": chrom,
            "protein": protein
        }
    }
    
    output_path = repo_root / "utils/recall_improvement_results.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"\n=== Saved results to {output_path} ===")
    

if __name__ == '__main__':
    main()
