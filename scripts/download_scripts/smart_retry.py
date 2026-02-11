#!/usr/bin/env python3
"""
Smart retry: Find PMIDs that passed filter but weren't downloaded, and retry them.
Uses the actual _FULL_CONTEXT.md files to determine what was downloaded.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from harvesting import PMCHarvester


def main():
    base_dir = Path("/home/kronckbm/gvf_output/KCNH2/20260128_210249")
    pmc_dir = base_dir / "pmc_fulltext"

    # Get all PMIDs from the PMID list (passed relevance filter)
    pmid_file = base_dir / "KCNH2_pmids.txt"
    all_pmids = set()
    if pmid_file.exists():
        with open(pmid_file) as f:
            all_pmids = set(line.strip() for line in f if line.strip().isdigit())
    print(f"Total PMIDs in list: {len(all_pmids)}")

    # Get PMIDs that were actually downloaded (have _FULL_CONTEXT.md)
    downloaded = set()
    for md_file in pmc_dir.glob("*_FULL_CONTEXT.md"):
        pmid = md_file.stem.replace("_FULL_CONTEXT", "")
        downloaded.add(pmid)
    print(f"Already downloaded: {len(downloaded)}")

    # Find missing PMIDs
    missing = all_pmids - downloaded
    print(f"Missing PMIDs to retry: {len(missing)}")

    if not missing:
        print("Nothing to retry!")
        return

    # Get PMIDs that passed filter (from filtered_out.csv, invert it)
    filtered_out_file = base_dir / "pmid_status" / "filtered_out.csv"
    filtered_out = set()
    if filtered_out_file.exists():
        with open(filtered_out_file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split(",")
                if parts and parts[0].isdigit():
                    filtered_out.add(parts[0])

    # Only retry PMIDs that weren't filtered out
    retry_pmids = missing - filtered_out
    print(f"PMIDs to retry (passed filter but not downloaded): {len(retry_pmids)}")

    if not retry_pmids:
        print("All missing PMIDs were filtered out - nothing to retry")
        return

    # Create harvester and retry
    harvester = PMCHarvester(output_dir=str(pmc_dir), gene_symbol="KCNH2")

    successful = 0
    failed = 0
    retry_list = sorted(retry_pmids)

    for idx, pmid in enumerate(retry_list, 1):
        print(f"\n[{idx}/{len(retry_list)}] Retrying PMID {pmid}...")

        try:
            success, result, _ = harvester.download_pmid(pmid)
            if success:
                successful += 1
                print("  ✅ Success!")
            else:
                failed += 1
                print(f"  ❌ Failed: {result}")
        except Exception as e:
            failed += 1
            print(f"  ❌ Error: {e}")

        time.sleep(2)

    print(f"\n{'=' * 60}")
    print("Retry complete!")
    print(f"  ✅ Successful: {successful}")
    print(f"  ❌ Failed: {failed}")
    print(f"  New total downloaded: {len(downloaded) + successful}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
