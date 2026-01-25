#!/usr/bin/env python3
"""
Organizes downloaded PDFs from ~/Downloads into the target directory.
Matches files by PMID patterns in filenames.

Usage:
    python scripts/organize_downloads.py TARGET_DIR --pmids 12345,67890
"""

import re
import shutil
import sys
from datetime import datetime, timedelta
from pathlib import Path


def find_downloads(pmids: list, since_minutes: int = 60) -> dict:
    """Find recently downloaded files that might match PMIDs."""
    downloads = Path.home() / "Downloads"
    cutoff = datetime.now() - timedelta(minutes=since_minutes)

    found = {}

    for f in downloads.iterdir():
        if not f.is_file():
            continue
        if f.suffix.lower() not in [".pdf", ".html", ".htm"]:
            continue

        # Check if recent
        mtime = datetime.fromtimestamp(f.stat().st_mtime)
        if mtime < cutoff:
            continue

        # Try to match PMID in filename
        fname = f.stem.lower()
        for pmid in pmids:
            if pmid in fname:
                if pmid not in found:
                    found[pmid] = []
                found[pmid].append(f)
                break

    return found


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Organize downloaded papers")
    parser.add_argument("target_dir", type=Path)
    parser.add_argument("--pmids", type=str, required=True, help="PMIDs to look for")
    parser.add_argument(
        "--since", type=int, default=60, help="Look at files from last N minutes"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Show what would be done"
    )
    args = parser.parse_args()

    pmids = [p.strip() for p in args.pmids.split(",")]
    args.target_dir.mkdir(parents=True, exist_ok=True)

    print(f"Looking for {len(pmids)} PMIDs in ~/Downloads (last {args.since} min)...")

    found = find_downloads(pmids, args.since)

    if not found:
        print("\nNo matching files found in Downloads.")
        print("Make sure you downloaded PDFs and the filenames contain the PMID.")

        # Show recent PDFs
        downloads = Path.home() / "Downloads"
        recent = sorted(
            [f for f in downloads.iterdir() if f.suffix.lower() == ".pdf"],
            key=lambda f: f.stat().st_mtime,
            reverse=True,
        )[:10]

        if recent:
            print("\nRecent PDFs in Downloads:")
            for f in recent:
                print(f"  {f.name}")
        return

    print(f"\nFound files for {len(found)} PMIDs:\n")

    moved = 0
    for pmid, files in found.items():
        for f in files:
            dest_name = f"PMID_{pmid}_supplement{f.suffix}"
            dest = args.target_dir / dest_name

            print(f"  {f.name} -> {dest_name}")

            if not args.dry_run:
                shutil.move(str(f), str(dest))
                moved += 1

    if args.dry_run:
        print(f"\nDry run - would move {len(found)} files")
    else:
        print(f"\nMoved {moved} files to {args.target_dir}")


if __name__ == "__main__":
    main()
