#!/usr/bin/env python3
"""
Post-download processing: convert PDFs to markdown and run data scout.

Run after browser_download.py to prepare files for extraction.

Usage:
    python tests/run_post_download.py
"""

import subprocess
import sys
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
VENV_PYTHON = PROJECT_ROOT / "venv" / "bin" / "python"
CSV_FILE = (
    PROJECT_ROOT
    / "tests"
    / "test_run_output"
    / "KCNH2"
    / "downloads"
    / "test_pmids_download.csv"
)
TARGET_DIR = (
    PROJECT_ROOT / "tests" / "test_run_output" / "KCNH2" / "downloads" / "pmc_fulltext"
)


def main():
    print("=" * 70)
    print("POST-DOWNLOAD PROCESSING")
    print("=" * 70)
    print(f"Target: {TARGET_DIR}")
    print("=" * 70)
    print()

    # Build command for fetch_manager in process-only mode
    cmd = [
        str(VENV_PYTHON),
        "-m",
        "cli.fetch_manager",
        str(CSV_FILE),
        "--target-dir",
        str(TARGET_DIR),
        "--process-only",  # Skip download loop, just process
        "--convert",  # Convert PDFs to markdown
        "--run-scout",  # Run data scout
        "--gene",
        "KCNH2",
    ]

    print(f"Command: {' '.join(cmd)}")
    print()

    result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))

    if result.returncode == 0:
        print("\n" + "=" * 70)
        print("POST-PROCESSING COMPLETE")
        print("=" * 70)
        print(f"\nFiles in: {TARGET_DIR}")
        print("  - *_FULL_CONTEXT.md: Combined paper text")
        print("  - *_DATA_ZONES.md: Condensed high-value sections")
        print("\nNext: Run extraction pipeline")

    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
