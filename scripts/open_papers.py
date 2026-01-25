#!/usr/bin/env python3
"""
Opens paper URLs in your default browser for manual download.
Run this, download PDFs manually, then run organize_downloads.py.

Usage:
    python scripts/open_papers.py paywalled_missing.csv --max 10
"""

import csv
import sys
import time
import webbrowser
from pathlib import Path


def get_unique_pmids(csv_path: Path) -> list:
    """Extract unique PMIDs and their best URL from CSV."""
    pmids = {}
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pmid = row.get("PMID", "").strip()
            url = row.get("URL", "").strip()
            if pmid and url:
                # Skip non-article URLs
                if any(
                    skip in url
                    for skip in [
                        "clinicaltrials.gov",
                        "coriell.org",
                        "lens.org",
                        "medlineplus.gov",
                    ]
                ):
                    continue
                # Prefer DOI URLs
                if pmid not in pmids or "doi.org" in url:
                    pmids[pmid] = url
    return list(pmids.items())


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Open paper URLs in browser")
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--max", type=int, default=10, help="Max papers to open")
    parser.add_argument("--start", type=int, default=0, help="Start at paper N")
    parser.add_argument(
        "--delay", type=float, default=2.0, help="Seconds between opens"
    )
    parser.add_argument("--pmids", type=str, help="Specific PMIDs (comma-separated)")
    args = parser.parse_args()

    papers = get_unique_pmids(args.csv_file)

    if args.pmids:
        specific = set(args.pmids.split(","))
        papers = [(p, u) for p, u in papers if p in specific]

    papers = papers[args.start : args.start + args.max]

    print(f"\nOpening {len(papers)} papers in your browser...")
    print("Download each PDF manually (click download button or Cmd+S)")
    print("Files go to ~/Downloads, then run organize_downloads.py\n")

    for i, (pmid, url) in enumerate(papers):
        print(f"[{i+1}/{len(papers)}] PMID {pmid}: {url[:60]}...")
        webbrowser.open(url)
        if i < len(papers) - 1:
            time.sleep(args.delay)

    print(f"\nDone! Download PDFs from the {len(papers)} browser tabs.")
    print(f"PMIDs opened: {', '.join(p for p, u in papers)}")


if __name__ == "__main__":
    main()
