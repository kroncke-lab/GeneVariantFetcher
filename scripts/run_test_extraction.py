#!/usr/bin/env python3
"""
Test extraction run for validating Phase 3 features.
Run: 20260206_143433
"""

import sys
from pathlib import Path

from harvesting.orchestrator import PMCHarvester
from utils.run_manifest import RunManifest
from utils.variant_normalizer import VariantNormalizer, normalize_variant_list

OUTPUT_DIR = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2_test_20260206")

# Test PMIDs (known KCNH2 papers with variants)
TEST_PMIDS = ["11320260", "15466642", "17905336", "19841300", "20031618"]


def main():
    print("=" * 60)
    print("GVF Phase 3 Test Run - 20260206_143433")
    print("=" * 60)

    # Create run manifest
    manifest = RunManifest(OUTPUT_DIR, "KCNH2")
    manifest.set_config(test_run=True, pmid_count=len(TEST_PMIDS))
    manifest.save()
    print(f"\n✓ Created manifest: {manifest.run_id[:8]}...")

    # Initialize harvester
    print("\nInitializing harvester...")
    harvester = PMCHarvester(str(OUTPUT_DIR), gene_symbol="KCNH2")
    print("✓ Harvester initialized with circuit breakers")
    print(f"  Elsevier circuit: {harvester.elsevier_circuit.state}")
    print(f"  Wiley circuit: {harvester.wiley_circuit.state}")
    print(f"  Springer circuit: {harvester.springer_circuit.state}")

    print(f"\nTest PMIDs: {TEST_PMIDS}")

    # Process each PMID
    print("\n" + "-" * 60)
    print("Processing PMIDs...")
    print("-" * 60)

    results = []
    for pmid in TEST_PMIDS:
        print(f"\n[{pmid}] Processing...")
        try:
            # Just test the download/conversion (not full extraction)
            result = harvester.process_pmid(pmid)
            results.append(
                {"pmid": pmid, "status": "success", "result": str(result)[:100]}
            )
            print(f"[{pmid}] ✓ Success")
        except Exception as e:
            results.append({"pmid": pmid, "status": "error", "error": str(e)[:100]})
            print(f"[{pmid}] ✗ Error: {str(e)[:50]}")

    # Test variant normalization
    print("\n" + "-" * 60)
    print("Testing variant normalization...")
    print("-" * 60)

    test_variants = [
        {"protein_notation": "A561V"},
        {"protein_notation": "p.G628S"},
        {"protein_notation": "Arg534Cys"},
        {"protein_notation": "T613M"},
        {"protein_notation": "A561V"},  # duplicate
    ]

    normalized, stats = normalize_variant_list(test_variants, "KCNH2")
    print("\nNormalization stats:")
    print(f"  Input: {stats['total_input']}")
    print(f"  Output: {stats['total_output']}")
    print(f"  Normalized: {stats['normalized']}")
    print(f"  Duplicates removed: {stats['duplicates_removed']}")

    print("\nNormalized variants:")
    for v in normalized:
        print(f"  {v.get('canonical_key')}")

    # Finalize manifest
    manifest.update_statistics(
        pmids_processed=len(TEST_PMIDS),
        successful=len([r for r in results if r["status"] == "success"]),
        failed=len([r for r in results if r["status"] == "error"]),
    )
    manifest.finalize()
    print("\n✓ Manifest finalized")

    # Summary
    print("\n" + "=" * 60)
    print("TEST RUN COMPLETE")
    print("=" * 60)
    print(f"Output: {OUTPUT_DIR}")
    print(f"Manifest: {OUTPUT_DIR / 'run_manifest.json'}")
    print(f"Status files: {OUTPUT_DIR / 'pmid_status/'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
