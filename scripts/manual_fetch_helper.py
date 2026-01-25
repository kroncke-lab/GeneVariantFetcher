#!/usr/bin/env python3
"""
Manual Fetch Helper - Opens papers in your real browser for manual download.

This avoids CAPTCHA issues by using your actual browser (with cookies, history, etc.)
instead of Playwright automation which gets detected.

Usage:
    python scripts/manual_fetch_helper.py paywalled_missing.csv --output-dir pmc_fulltext

Workflow:
    1. Opens each paper URL in your default browser
    2. You download the PDF manually (Cmd+S or click download)
    3. Press Enter to continue to next paper, 's' to skip, 'q' to quit
    4. After you're done, run --organize to move files from Downloads to target dir
"""

import csv
import os
import re
import shutil
import subprocess
import sys
import webbrowser
from datetime import datetime
from pathlib import Path
from typing import Optional

# Common download locations
DOWNLOAD_DIRS = [
    Path.home() / "Downloads",
    Path.home() / "Desktop",
]


def get_unique_pmids(csv_path: Path) -> list:
    """Extract unique PMIDs and their first URL from CSV."""
    pmids = {}
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pmid = row.get("PMID", "").strip()
            url = row.get("URL", "").strip()
            if pmid and pmid not in pmids and url:
                # Skip non-article URLs
                if "clinicaltrials.gov" in url or "coriell.org" in url:
                    continue
                pmids[pmid] = url
    return list(pmids.items())


def find_recent_downloads(pmid: str, since: datetime, extensions: list = None) -> list:
    """Find files downloaded since a given time that might be for this PMID."""
    if extensions is None:
        extensions = [".pdf", ".html", ".htm"]

    found = []
    for dl_dir in DOWNLOAD_DIRS:
        if not dl_dir.exists():
            continue
        for f in dl_dir.iterdir():
            if not f.is_file():
                continue
            # Check modification time
            mtime = datetime.fromtimestamp(f.stat().st_mtime)
            if mtime < since:
                continue
            # Check extension
            if f.suffix.lower() not in extensions:
                continue
            found.append(f)

    return sorted(found, key=lambda f: f.stat().st_mtime, reverse=True)


def organize_downloads(output_dir: Path, download_log: Path):
    """Move downloaded files from Downloads to output directory."""
    if not download_log.exists():
        print(f"No download log found at {download_log}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    with open(download_log, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            pmid, filepath = parts
            src = Path(filepath)
            if not src.exists():
                print(f"  File not found: {src}")
                continue

            # Determine destination name
            ext = src.suffix.lower()
            if ext == ".pdf":
                dst = output_dir / f"{pmid}_supplement.pdf"
            else:
                dst = output_dir / f"PMID_{pmid}_FULL_CONTEXT{ext}"

            shutil.move(str(src), str(dst))
            print(f"  Moved: {src.name} -> {dst.name}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Manual fetch helper for paywalled papers"
    )
    parser.add_argument("csv_file", type=Path, help="CSV file with PMID,URL columns")
    parser.add_argument("--output-dir", type=Path, help="Directory to save files")
    parser.add_argument(
        "--organize", action="store_true", help="Organize previously downloaded files"
    )
    parser.add_argument("--start-at", type=int, default=0, help="Start at paper N")
    parser.add_argument("--pmids", type=str, help="Comma-separated specific PMIDs")
    args = parser.parse_args()

    if not args.csv_file.exists():
        print(f"CSV file not found: {args.csv_file}")
        sys.exit(1)

    # Set up paths
    output_dir = args.output_dir or args.csv_file.parent
    download_log = output_dir / "manual_download_log.txt"

    if args.organize:
        organize_downloads(output_dir, download_log)
        return

    # Get papers to process
    papers = get_unique_pmids(args.csv_file)

    if args.pmids:
        specific = set(args.pmids.split(","))
        papers = [(p, u) for p, u in papers if p in specific]

    papers = papers[args.start_at :]

    print(f"\n{'='*60}")
    print(f"Manual Fetch Helper")
    print(f"{'='*60}")
    print(f"Papers to process: {len(papers)}")
    print(f"Output directory: {output_dir}")
    print(f"\nInstructions:")
    print(f"  1. Each paper will open in your browser")
    print(f"  2. Download the PDF (Cmd+S or click download button)")
    print(f"  3. Press Enter when done, 's' to skip, 'q' to quit")
    print(f"{'='*60}\n")

    downloaded = []

    for i, (pmid, url) in enumerate(papers):
        print(f"\n[{i+1}/{len(papers)}] PMID {pmid}")
        print(f"  URL: {url[:80]}...")

        # Record start time
        start_time = datetime.now()

        # Open in default browser
        webbrowser.open(url)

        # Wait for user
        response = (
            input("  Press Enter when downloaded, 's' to skip, 'q' to quit: ")
            .strip()
            .lower()
        )

        if response == "q":
            print("\nQuitting...")
            break
        elif response == "s":
            print("  Skipped")
            continue

        # Look for recently downloaded files
        recent = find_recent_downloads(pmid, start_time)
        if recent:
            print(f"  Found {len(recent)} recent download(s):")
            for j, f in enumerate(recent[:3]):
                print(f"    {j+1}. {f.name}")

            choice = input(
                "  Enter number to use (or Enter for first, 'n' for none): "
            ).strip()

            if choice == "n":
                continue
            elif choice.isdigit() and 1 <= int(choice) <= len(recent):
                selected = recent[int(choice) - 1]
            else:
                selected = recent[0]

            # Log the download
            with open(download_log, "a") as f:
                f.write(f"{pmid}\t{selected}\n")
            downloaded.append((pmid, selected))
            print(f"  Logged: {selected.name}")
        else:
            print("  No recent downloads found. Continuing...")

    print(f"\n{'='*60}")
    print(f"Session complete!")
    print(f"  Downloaded: {len(downloaded)} papers")
    print(f"  Log file: {download_log}")
    print(f"\nTo organize downloads into {output_dir}, run:")
    print(f"  python {__file__} {args.csv_file} --output-dir {output_dir} --organize")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
