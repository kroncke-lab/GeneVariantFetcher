#!/usr/bin/env python3
"""
Retry downloading papers that failed in the initial harvest.
Focuses on papers most likely to succeed on retry.
"""

import csv
import sys
import time
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from harvesting import PMCHarvester


def get_failed_pmids(paywalled_csv: Path, success_csv: Path) -> list:
    """Get unique PMIDs that failed but might succeed on retry."""

    # Get already successful PMIDs
    successful = set()
    if success_csv.exists():
        with open(success_csv) as f:
            reader = csv.reader(f)
            next(reader, None)  # skip header
            for row in reader:
                if row:
                    successful.add(row[0])

    # Get failed PMIDs with their reasons
    failed_reasons = {}
    if paywalled_csv.exists():
        with open(paywalled_csv) as f:
            reader = csv.reader(f)
            next(reader, None)  # skip header
            for row in reader:
                if row and row[0].isdigit():
                    pmid = row[0]
                    reason = row[1] if len(row) > 1 else ""
                    if pmid not in failed_reasons:
                        failed_reasons[pmid] = []
                    failed_reasons[pmid].append(reason)

    # Filter to PMIDs worth retrying
    # Skip: truly paywalled, junk content, validation failures
    skip_patterns = [
        "No PMCID found, not free full text",
        "Junk content detected",
        "Content from non-article domain",
        "Content too short",
        "Missing paper structure",
    ]

    retry_pmids = []
    for pmid, reasons in failed_reasons.items():
        if pmid in successful:
            continue

        # Check if all reasons are skip-worthy
        all_skip = all(
            any(skip in reason for skip in skip_patterns) for reason in reasons
        )

        if not all_skip:
            retry_pmids.append(pmid)

    return retry_pmids


def main():
    base_dir = Path("/home/kronckbm/gvf_output/KCNH2/20260128_210249")
    pmc_dir = base_dir / "pmc_fulltext"

    paywalled_csv = pmc_dir / "paywalled_missing.csv"
    success_csv = pmc_dir / "successful_downloads.csv"

    # Get PMIDs to retry
    retry_pmids = get_failed_pmids(paywalled_csv, success_csv)
    print(f"Found {len(retry_pmids)} PMIDs worth retrying")

    if not retry_pmids:
        print("No PMIDs to retry")
        return

    # Create harvester pointing to same output directory
    harvester = PMCHarvester(output_dir=str(pmc_dir), gene_symbol="KCNH2")

    # Retry each PMID
    successful = 0
    failed = 0

    for idx, pmid in enumerate(retry_pmids, 1):
        print(f"\n[{idx}/{len(retry_pmids)}] Retrying PMID {pmid}...")

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

        # Rate limiting
        time.sleep(2)

    print(f"\n{'=' * 60}")
    print("Retry complete!")
    print(f"  ✅ Successful: {successful}")
    print(f"  ❌ Failed: {failed}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
