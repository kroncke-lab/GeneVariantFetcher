#!/usr/bin/env python3
"""
Phase 1: Batch Download Missing PMIDs
Scripts to download the TOP 10 high-impact PMIDs for GVF pipeline completion.
"""

import os
import sys
from pathlib import Path

# Ensure this script runs from any directory
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from harvesting.orchestrator import PMCHarvester


def main():
    # Setup paths
    gvf_repo = Path("/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/")
    output_base = Path("/mnt/temp2/kronckbm/gvf_output/")
    output_dir = output_base / "download_batch_20260207"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Top 10 high-impact PMIDs for testing
    with open("/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt", "r") as f:
        high_impact_pmids = [line.strip() for line in f]

    print("=" * 80)
    print("GVF PHASE 1: DOWNLOADING HIGH-IMPACT PMIDs")
    print("=" * 80)
    print(f"Total PMIDs to process: {len(high_impact_pmids)}")
    print(f"Output directory: {output_dir.absolute()}")
    print("Target PMIDs:", " ".join(high_impact_pmids))
    print()

    # Initialize harvester
    harvester = PMCHarvester(output_dir=str(output_dir), gene_symbol="batch_download")

    # Download all PMIDs with harvest orchestrator
    try:
        harvester.harvest(
            pmids=high_impact_pmids,
            delay=3.0,  # Increased delay for robustness
            run_scout=False,  # Do post-processing separately if needed
        )

        print("\n" + "=" * 80)
        print("DOWNLOAD BATCH COMPLETE!")
        print("=" * 80)
        print("Results saved to:")
        print(f"  - Harvest output: {output_dir}")
        print(f"  - Success log: {output_dir}/successful_downloads.csv")
        print(f"  - Paywall log: {output_dir}/paywalled_missing.csv")
        print(f"  - Download log: {output_base}/download_log_20260206.md")

    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        return False

    return True


if __name__ == "__main__":
    # Ensure we're in the GVF directory
    # os.chdir("/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/")
    success = main()
    sys.exit(0 if success else 1)
