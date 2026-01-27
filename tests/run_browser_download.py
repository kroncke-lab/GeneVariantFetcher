#!/usr/bin/env python3
"""
Automated browser download of test PMIDs.

Run with Claude Code cowork mode for hands-free paper downloading.
Uses Vanderbilt VPN for institutional access.

Usage:
    python tests/run_browser_download.py [--max N] [--use-claude]
"""

import sys
import subprocess
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
    import argparse

    parser = argparse.ArgumentParser(description="Run browser-based paper downloads")
    parser.add_argument("--max", type=int, help="Maximum papers to process")
    parser.add_argument(
        "--use-claude",
        action="store_true",
        help="Use Claude API to find download links",
    )
    parser.add_argument(
        "--headless", action="store_true", help="Run without visible browser"
    )
    parser.add_argument(
        "--retry-paywall", action="store_true", help="Retry only paywalled papers"
    )
    parser.add_argument(
        "--retry-captcha", action="store_true", help="Retry CAPTCHA-blocked papers"
    )
    parser.add_argument("--pmids", type=str, help="Comma-separated specific PMIDs")
    args = parser.parse_args()

    # Ensure target directory exists
    TARGET_DIR.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [
        str(VENV_PYTHON),
        "-m",
        "cli.browser_fetch",
        str(CSV_FILE),
        "--target-dir",
        str(TARGET_DIR),
        "--use-profile",  # Use real Chrome profile for VPN/cookies
        "--slow-mo",
        "300",  # Slow down to avoid detection
    ]

    if args.max:
        cmd.extend(["--max", str(args.max)])
    if args.use_claude:
        cmd.append("--use-claude")
    if args.headless:
        cmd.append("--headless")
    if args.retry_paywall:
        cmd.append("--retry-paywall")
    if args.retry_captcha:
        cmd.append("--retry-captcha")
    if args.pmids:
        cmd.extend(["--pmids", args.pmids])

    print("=" * 70)
    print("BROWSER PAPER DOWNLOAD")
    print("=" * 70)
    print(f"CSV: {CSV_FILE}")
    print(f"Target: {TARGET_DIR}")
    print(f"Command: {' '.join(cmd)}")
    print("=" * 70)
    print()

    # Run browser fetch
    result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))

    if result.returncode == 0:
        print("\n" + "=" * 70)
        print("DOWNLOAD PHASE COMPLETE")
        print("=" * 70)
        print(f"\nFiles saved to: {TARGET_DIR}")
        print("\nNext steps:")
        print("  1. Check paywalled_missing.csv for papers that need manual download")
        print("  2. Run extraction: python tests/run_test_pipeline.py")

    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
