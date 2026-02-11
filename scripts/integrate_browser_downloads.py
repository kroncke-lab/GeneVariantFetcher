#!/usr/bin/env python3
"""
Integrate browser_downloads into the main pipeline.

This script:
1. Organizes browser_downloads by PMID
2. Converts main PDFs to markdown (FULL_CONTEXT.md)
3. Moves supplements to {PMID}_supplements/
4. Runs Data Scout to create DATA_ZONES.md
5. Prepares for extraction

Usage:
    python scripts/integrate_browser_downloads.py [--dry-run]
"""

import argparse
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path

# Add project root
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.settings import get_settings
from harvesting.format_converters import FormatConverter
from pipeline.data_scout import GeneticDataScout


def parse_browser_filename(filename):
    """Parse {PMID}_{rest} or {PMID}_Supp_{N}_{rest} format."""
    parts = filename.split("_")
    if not parts[0].isdigit():
        return None, None, None

    pmid = parts[0]

    if len(parts) >= 3 and parts[1] == "Supp":
        # Supplement file: {PMID}_Supp_{N}_{filename}
        return pmid, "supplement", filename
    else:
        # Main file: {PMID}_{filename}
        return pmid, "main", filename


def main():
    parser = argparse.ArgumentParser(description="Integrate browser downloads")
    parser.add_argument("--dry-run", action="store_true", help="Show what would happen")
    parser.add_argument(
        "--source", default="gvf_output/browser_downloads", help="Source directory"
    )
    parser.add_argument(
        "--target",
        default="gvf_output/KCNH2/20260128_210249/pmc_fulltext",
        help="Target pmc_fulltext directory",
    )
    parser.add_argument("--gene", default="KCNH2", help="Gene symbol")
    args = parser.parse_args()

    source = PROJECT_ROOT / args.source
    target = PROJECT_ROOT / args.target

    if not source.exists():
        print(f"❌ Source not found: {source}")
        return 1

    print("=" * 70)
    print("BROWSER DOWNLOADS INTEGRATION")
    print("=" * 70)
    print(f"Source: {source}")
    print(f"Target: {target}")
    print(f"Gene: {args.gene}")
    print(f"Mode: {'DRY RUN' if args.dry_run else 'LIVE'}")
    print("=" * 70)

    # Parse all files
    files_by_pmid = defaultdict(lambda: {"main": [], "supplement": []})

    for f in os.listdir(source):
        filepath = source / f
        if filepath.is_dir():
            continue
        pmid, ftype, fname = parse_browser_filename(f)
        if pmid:
            files_by_pmid[pmid][ftype].append(filepath)

    print(f"\nFound {len(files_by_pmid)} unique PMIDs")

    # Check which already have FULL_CONTEXT
    existing = set()
    for f in os.listdir(target):
        if f.endswith("_FULL_CONTEXT.md"):
            existing.add(f.replace("_FULL_CONTEXT.md", ""))

    new_pmids = [p for p in files_by_pmid if p not in existing]
    update_pmids = [p for p in files_by_pmid if p in existing]

    print(f"  - NEW (no FULL_CONTEXT): {len(new_pmids)}")
    print(f"  - UPDATE (add supplements): {len(update_pmids)}")

    if args.dry_run:
        print("\n### DRY RUN - Actions that would be taken ###\n")

        for pmid in new_pmids:
            files = files_by_pmid[pmid]
            print(f"PMID {pmid}:")
            for f in files["main"]:
                print(f"  MAIN: {f.name} → {pmid}_FULL_CONTEXT.md")
            for f in files["supplement"]:
                print(f"  SUPP: {f.name} → {pmid}_supplements/")

        for pmid in update_pmids:
            files = files_by_pmid[pmid]
            if files["supplement"]:
                print(f"PMID {pmid} (update):")
                for f in files["supplement"]:
                    print(f"  ADD SUPP: {f.name} → {pmid}_supplements/")

        print("\n✅ Dry run complete. Run without --dry-run to execute.")
        return 0

    # Initialize converter
    converter = FormatConverter()

    # Process NEW PMIDs
    processed = 0
    for pmid in new_pmids:
        files = files_by_pmid[pmid]
        print(f"\nProcessing {pmid}...")

        # Create supplements directory
        supp_dir = target / f"{pmid}_supplements"
        supp_dir.mkdir(exist_ok=True)

        # Convert main file(s) to markdown
        full_text = ""
        for main_file in files["main"]:
            print(f"  Converting: {main_file.name}")
            try:
                ext = main_file.suffix.lower()
                if ext == ".pdf":
                    full_text += converter.pdf_to_markdown(main_file)
                elif ext == ".docx":
                    full_text += converter.docx_to_markdown(main_file)
                elif ext == ".doc":
                    full_text += converter.doc_to_markdown(main_file)
                elif ext in [".html", ".htm"]:
                    full_text += main_file.read_text(errors="ignore")
                else:
                    full_text += f"[File: {main_file.name}]\n"
                full_text += "\n\n"
            except Exception as e:
                print(f"    ⚠ Error converting: {e}")

        # Process supplements
        supp_text = ""
        for idx, supp_file in enumerate(files["supplement"], 1):
            print(f"  Supplement: {supp_file.name}")
            # Copy to supplements dir
            dest = supp_dir / supp_file.name
            shutil.copy2(supp_file, dest)

            # Also convert to text
            try:
                ext = supp_file.suffix.lower()
                supp_text += f"\n\n# SUPPLEMENTAL FILE {idx}: {supp_file.name}\n\n"
                if ext == ".pdf":
                    supp_text += converter.pdf_to_markdown(supp_file)
                elif ext in [".xlsx", ".xls"]:
                    supp_text += converter.excel_to_markdown(supp_file)
                elif ext == ".docx":
                    supp_text += converter.docx_to_markdown(supp_file)
                elif ext == ".doc":
                    supp_text += converter.doc_to_markdown(supp_file)
                else:
                    supp_text += f"[File available: {supp_file.name}]\n"
            except Exception as e:
                print(f"    ⚠ Error converting supplement: {e}")

        # Write FULL_CONTEXT.md
        combined = f"# MAIN TEXT\n\n{full_text}\n\n{supp_text}"
        output_file = target / f"{pmid}_FULL_CONTEXT.md"
        output_file.write_text(combined, encoding="utf-8")
        print(f"  ✓ Wrote: {output_file.name} ({len(combined):,} chars)")
        processed += 1

    # Process UPDATE PMIDs (add supplements to existing)
    for pmid in update_pmids:
        files = files_by_pmid[pmid]
        if not files["supplement"]:
            continue

        print(f"\nUpdating {pmid} (adding supplements)...")
        supp_dir = target / f"{pmid}_supplements"
        supp_dir.mkdir(exist_ok=True)

        for supp_file in files["supplement"]:
            dest = supp_dir / supp_file.name
            if not dest.exists():
                shutil.copy2(supp_file, dest)
                print(f"  ✓ Added: {supp_file.name}")

    print(f"\n{'=' * 70}")
    print("INTEGRATION COMPLETE")
    print(f"  New FULL_CONTEXT.md files: {processed}")
    print(f"  Updated with supplements: {len(update_pmids)}")
    print("\nNext step: Run Data Scout on new files")
    print(f"  python -m cli.scout {target} --gene {args.gene}")
    print(f"{'=' * 70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
