#!/usr/bin/env python3
"""Test harvesting for PMID 19716085 (KCNH2)"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from harvesting import PMCHarvester
from pathlib import Path

def main():
    pmid = "19716085"
    output_dir = Path("test_kcnh2_19716085")

    print(f"Testing PMID: {pmid}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)

    # Create harvester
    harvester = PMCHarvester(output_dir=str(output_dir))

    # Process the PMID
    success, result = harvester.process_pmid(pmid)

    if success:
        print(f"\n✓ SUCCESS: Created {result}")

        # Read and show content summary
        with open(result, 'r') as f:
            content = f.read()

        print(f"\nContent length: {len(content):,} characters")
        print(f"Line count: {len(content.splitlines()):,} lines")

        # Check for supplements
        if "## Supplementary" in content or "supplement" in content.lower():
            print("✓ Contains supplementary material references")

        # Show first part of content
        print("\n" + "=" * 50)
        print("FIRST 3000 CHARACTERS:")
        print("=" * 50)
        print(content[:3000])

        # Show what sections exist
        print("\n" + "=" * 50)
        print("SECTION HEADERS FOUND:")
        print("=" * 50)
        for line in content.splitlines():
            if line.startswith('#'):
                print(line)

        # Check for supplement files
        supp_dir = output_dir / f"{pmid}_supplements"
        if supp_dir.exists():
            supp_files = list(supp_dir.iterdir())
            print(f"\n✓ Supplement directory exists with {len(supp_files)} files:")
            for f in supp_files:
                print(f"  - {f.name}")
        else:
            print(f"\n⚠ No supplement directory at {supp_dir}")

    else:
        print(f"\n✗ FAILED: {result}")

        # Check logs
        for log_file in ["paywalled_missing.csv", "successful_downloads.csv"]:
            log_path = output_dir / log_file
            if log_path.exists():
                print(f"\nLog file {log_file}:")
                with open(log_path) as f:
                    print(f.read())

if __name__ == "__main__":
    main()
