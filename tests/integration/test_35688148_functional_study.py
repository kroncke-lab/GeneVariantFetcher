#!/usr/bin/env python3
"""Test extraction for PMID 35688148 (massively parallel functional assay paper)

This tests the fix for functional study detection - the LLM should now:
1. Recognize this as a functional study
2. NOT extract assay replicates (n=cells) as patient carriers
3. Set penetrance_data counts to null for most variants
4. Only extract actual clinical data from LQTS/Unaffected columns
"""

import sys
import os
import json
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Only load .env if needed - shell env vars take precedence
if not os.getenv("NCBI_EMAIL"):
    from dotenv import load_dotenv

    load_dotenv()

from pipeline.extraction import ExpertExtractor
from utils.models import Paper


def main():
    pmid = "35688148"
    gene = "KCNH2"
    output_dir = Path("output/KCNH2/20260114_204346")

    # Read the DATA_ZONES.md file (this is what the extractor will use)
    data_zones_path = output_dir / "pmc_fulltext" / f"{pmid}_DATA_ZONES.md"
    full_text_path = output_dir / "pmc_fulltext" / f"{pmid}_FULL_CONTEXT.md"

    if data_zones_path.exists():
        content = data_zones_path.read_text()
        print(f"Using DATA_ZONES.md: {len(content):,} chars")
    elif full_text_path.exists():
        content = full_text_path.read_text()
        print(f"Using FULL_CONTEXT.md: {len(content):,} chars")
    else:
        print("No content found!")
        sys.exit(1)

    # Create Paper object
    paper = Paper(
        pmid=pmid,
        title="A massively parallel assay accurately discriminates between functionally normal and abnormal variants in a hotspot domain of KCNH2.",
        abstract="",
        full_text=content,
        gene_symbol=gene,
    )

    print(f"\nTesting extraction for PMID {pmid}")
    print(f"Gene: {gene}")
    print("-" * 50)

    # Run extraction with updated prompts
    extractor = ExpertExtractor(fulltext_dir=str(output_dir / "pmc_fulltext"))
    result = extractor.extract(paper)

    if not result.success:
        print(f"Extraction failed: {result.error}")
        sys.exit(1)

    data = result.extracted_data
    print(f"\nExtraction successful!")
    print(f"Model used: {result.model_used}")

    # Check extraction metadata
    meta = data.get("extraction_metadata", {})
    print(f"\n=== EXTRACTION METADATA ===")
    print(f"Study type: {meta.get('study_type', 'NOT SET')}")
    print(f"Total variants: {meta.get('total_variants_found', 0)}")
    print(f"Confidence: {meta.get('extraction_confidence', 'N/A')}")
    print(f"Notes: {meta.get('notes', 'N/A')}")

    # Check variants
    variants = data.get("variants", [])
    print(f"\n=== VARIANTS ({len(variants)} total) ===")

    # Show first 10 variants
    for i, v in enumerate(variants[:10]):
        print(f"\n  Variant {i+1}: {v.get('protein_notation', 'N/A')}")
        pdata = v.get("penetrance_data", {})
        carriers = pdata.get("total_carriers_observed")
        affected = pdata.get("affected_count")
        print(f"    Carriers: {carriers if carriers is not None else 'NULL'}")
        print(f"    Affected: {affected if affected is not None else 'NULL'}")
        print(f"    Source: {v.get('source_location', 'N/A')}")

    # Count how many have null penetrance data (expected for functional study)
    null_penetrance = sum(
        1
        for v in variants
        if v.get("penetrance_data", {}).get("total_carriers_observed") is None
    )

    print(f"\n=== SUMMARY ===")
    print(f"Total variants: {len(variants)}")
    print(f"Variants with NULL penetrance: {null_penetrance}")
    print(f"Variants with numeric penetrance: {len(variants) - null_penetrance}")

    # Check if the fix worked
    study_type = meta.get("study_type", "")
    if "functional" in study_type.lower():
        print("\n✓ SUCCESS: Study correctly identified as FUNCTIONAL")
    else:
        print(f"\n⚠ WARNING: Study type is '{study_type}' - expected 'functional'")

    # For a functional study, most variants should have null penetrance
    if null_penetrance > len(variants) * 0.5:
        print(
            f"✓ SUCCESS: {null_penetrance}/{len(variants)} variants have NULL penetrance (expected for functional study)"
        )
    else:
        print(
            f"⚠ WARNING: Only {null_penetrance}/{len(variants)} variants have NULL penetrance"
        )
        print(
            "   This suggests assay replicates may still be extracted as patient counts"
        )

    # Save result
    output_json = Path(f"tests/test_35688148_extraction_fixed.json")
    with open(output_json, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\nSaved to: {output_json}")


if __name__ == "__main__":
    main()
