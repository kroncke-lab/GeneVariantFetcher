#!/usr/bin/env python3
"""
Track failed downloads and extraction failures from GVF runs.
Outputs a summary of what didn't work and why.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path


def analyze_run(run_dir: Path):
    """Analyze a GVF run directory for failures."""
    run_dir = Path(run_dir)

    # Load workflow summary if exists
    summary_files = list(run_dir.glob("*_workflow_summary.json"))
    summary = {}
    if summary_files:
        with open(summary_files[0]) as f:
            summary = json.load(f)

    # Count discovered PMIDs
    pmid_file = list(run_dir.glob("*_pmids.txt"))
    discovered = set()
    if pmid_file:
        discovered = set(pmid_file[0].read_text().strip().split("\n"))

    # Check papers directory for downloaded content
    papers_dir = run_dir / "papers"
    downloaded = set()
    download_status = defaultdict(list)

    if papers_dir.exists():
        for paper_dir in papers_dir.iterdir():
            if paper_dir.is_dir():
                pmid = paper_dir.name
                full_context = paper_dir / "FULL_CONTEXT.md"
                data_zones = paper_dir / "DATA_ZONES.md"

                if full_context.exists() and full_context.stat().st_size > 100:
                    downloaded.add(pmid)
                    download_status["has_fulltext"].append(pmid)
                elif data_zones.exists() and data_zones.stat().st_size > 100:
                    downloaded.add(pmid)
                    download_status["has_data_zones_only"].append(pmid)
                else:
                    download_status["empty_or_missing"].append(pmid)

    # Check extractions
    extractions_dir = run_dir / "extractions"
    extracted = set()
    extraction_failures = []
    extraction_successes = []

    if extractions_dir.exists():
        for ext_file in extractions_dir.glob("*.json"):
            try:
                with open(ext_file) as f:
                    data = json.load(f)
                pmid = (
                    ext_file.stem.split("_PMID_")[1]
                    if "_PMID_" in ext_file.stem
                    else ext_file.stem
                )
                extracted.add(pmid)

                # Check if extraction found variants
                variants = data.get("variants", [])
                if variants:
                    extraction_successes.append((pmid, len(variants)))
                else:
                    extraction_failures.append((pmid, "no_variants_found"))
            except Exception as e:
                extraction_failures.append((ext_file.stem, f"parse_error: {e}"))

    # Calculate gaps
    not_downloaded = discovered - downloaded - extracted
    downloaded_not_extracted = downloaded - extracted

    # Print report
    print("=" * 80)
    print(f"GVF RUN ANALYSIS: {run_dir.name}")
    print("=" * 80)

    print("\n📊 SUMMARY:")
    print(f"  Discovered PMIDs:     {len(discovered)}")
    print(f"  Downloaded (any):     {len(downloaded)}")
    print(f"  Extracted:            {len(extracted)}")
    print(f"  With variants:        {len(extraction_successes)}")

    print("\n❌ GAPS:")
    print(f"  Not downloaded:       {len(not_downloaded)}")
    print(f"  Downloaded, not extracted: {len(downloaded_not_extracted)}")
    print(f"  Extracted, no variants:    {len(extraction_failures)}")

    print("\n📁 DOWNLOAD STATUS:")
    for status, pmids in download_status.items():
        print(f"  {status}: {len(pmids)}")

    # Save detailed failure lists
    failures_dir = run_dir / "failure_tracking"
    failures_dir.mkdir(exist_ok=True)

    # Not downloaded
    if not_downloaded:
        with open(failures_dir / "not_downloaded.txt", "w") as f:
            f.write("\n".join(sorted(not_downloaded)))
        print(
            f"\n📝 Saved {len(not_downloaded)} not-downloaded PMIDs to failure_tracking/not_downloaded.txt"
        )

    # Downloaded but not extracted
    if downloaded_not_extracted:
        with open(failures_dir / "downloaded_not_extracted.txt", "w") as f:
            f.write("\n".join(sorted(downloaded_not_extracted)))
        print(
            f"📝 Saved {len(downloaded_not_extracted)} downloaded-not-extracted PMIDs to failure_tracking/downloaded_not_extracted.txt"
        )

    # Extraction failures
    if extraction_failures:
        with open(failures_dir / "extraction_failures.json", "w") as f:
            json.dump(extraction_failures, f, indent=2)
        print(
            f"📝 Saved {len(extraction_failures)} extraction failures to failure_tracking/extraction_failures.json"
        )

    # Top extractions by variant count
    if extraction_successes:
        top_10 = sorted(extraction_successes, key=lambda x: x[1], reverse=True)[:10]
        print("\n🏆 TOP 10 EXTRACTIONS BY VARIANT COUNT:")
        for pmid, count in top_10:
            print(f"  PMID {pmid}: {count} variants")

    return {
        "discovered": len(discovered),
        "downloaded": len(downloaded),
        "extracted": len(extracted),
        "with_variants": len(extraction_successes),
        "not_downloaded": len(not_downloaded),
        "downloaded_not_extracted": len(downloaded_not_extracted),
        "extraction_failures": len(extraction_failures),
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python track_failures.py <run_directory>")
        print(
            "Example: python track_failures.py /mnt/temp2/kronckbm/gvf_output/KCNH2/20260202_173749/"
        )
        sys.exit(1)

    run_dir = Path(sys.argv[1])
    if not run_dir.exists():
        print(f"Error: {run_dir} does not exist")
        sys.exit(1)

    analyze_run(run_dir)
