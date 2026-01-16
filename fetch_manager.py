#!/usr/bin/env python3
"""
Fetch Manager: Semi-Manual Literature Retrieval Workflow

This script assists with manual downloading of papers that couldn't be automatically
fetched. It automates the file management and organization steps while you manually
download from publisher websites.

Workflow:
1. Opens each paper URL in your browser (DOI preferred, falls back to PubMed)
2. Waits for you to manually download PDF/supplements
3. Detects new files in your Downloads folder
4. Renames and moves them to the target directory with proper naming
5. Tracks progress in the CSV file

Usage:
    python fetch_manager.py <csv_file> [options]

    # Use a paywalled_missing.csv from a workflow run
    python fetch_manager.py output/GENE/20250101_120000/pmc_fulltext/paywalled_missing.csv

    # Use a custom CSV with PMID and DOI columns
    python fetch_manager.py my_papers.csv --doi-column DOI

    # Specify a custom target directory
    python fetch_manager.py papers.csv --target-dir ~/MyData/Papers
"""

import os
import sys
import time
import shutil
import webbrowser
import logging
from pathlib import Path
from datetime import datetime
from typing import Optional, Set, List, Tuple

import pandas as pd

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


# Default configuration
DEFAULT_DOWNLOADS_DIR = Path.home() / "Downloads"
DEFAULT_TARGET_DIR = Path.home() / "OneDrive - VUMC" / "Kroncke_Lab" / "Manual_Retrieval"

# File extensions to track
PDF_EXTENSIONS = {'.pdf'}
SPREADSHEET_EXTENSIONS = {'.xlsx', '.xls', '.csv', '.tsv'}
DOC_EXTENSIONS = {'.docx', '.doc'}
TEXT_EXTENSIONS = {'.txt'}
ARCHIVE_EXTENSIONS = {'.zip', '.rar', '.gz'}
ALL_TRACKED_EXTENSIONS = (
    PDF_EXTENSIONS
    | SPREADSHEET_EXTENSIONS
    | DOC_EXTENSIONS
    | TEXT_EXTENSIONS
    | ARCHIVE_EXTENSIONS
)

# Files to ignore
IGNORED_PATTERNS = {
    '.DS_Store',
    '.crdownload',  # Chrome partial download
    '.part',        # Firefox partial download
    '.tmp',
    '.download',
    'desktop.ini',
    'Thumbs.db',
}


def get_downloads_snapshot(downloads_dir: Path) -> Set[Path]:
    """
    Get a snapshot of current files in the Downloads directory.

    Args:
        downloads_dir: Path to Downloads folder

    Returns:
        Set of file paths currently in Downloads
    """
    files = set()
    if not downloads_dir.exists():
        logger.warning(f"Downloads directory does not exist: {downloads_dir}")
        return files

    for item in downloads_dir.iterdir():
        if item.is_file():
            # Skip hidden files and ignored patterns
            if item.name.startswith('.'):
                continue
            if any(pattern in item.name for pattern in IGNORED_PATTERNS):
                continue
            files.add(item)

    return files


def identify_new_files(before: Set[Path], after: Set[Path]) -> List[Path]:
    """
    Identify files that appeared between before and after snapshots.

    Args:
        before: Snapshot before download
        after: Snapshot after download

    Returns:
        List of new file paths
    """
    new_files = after - before

    # Filter out partial downloads and temporary files
    valid_files = []
    for f in new_files:
        if f.suffix.lower() in ALL_TRACKED_EXTENSIONS:
            valid_files.append(f)
        else:
            logger.debug(f"Ignoring file with untracked extension: {f.name}")

    return sorted(valid_files, key=lambda x: x.stat().st_mtime)


def generate_new_filename(pmid: str, original_file: Path, file_index: int = 0) -> str:
    """
    Generate a standardized filename based on PMID and file type.

    Args:
        pmid: PubMed ID
        original_file: Original file path
        file_index: Index for multiple files of same type

    Returns:
        New filename string
    """
    suffix = original_file.suffix.lower()

    if suffix in PDF_EXTENSIONS:
        base_name = f"{pmid}_Main_Text"
        if file_index > 0:
            base_name = f"{pmid}_Supplement_{file_index}"
        return f"{base_name}.pdf"

    elif suffix in SPREADSHEET_EXTENSIONS:
        base_name = f"{pmid}_Supp_Data"
        if file_index > 0:
            base_name = f"{pmid}_Supp_Data_{file_index}"
        return f"{base_name}{suffix}"

    elif suffix in DOC_EXTENSIONS:
        base_name = f"{pmid}_Supplement_Doc"
        if file_index > 0:
            base_name = f"{pmid}_Supplement_Doc_{file_index}"
        return f"{base_name}{suffix}"

    elif suffix in TEXT_EXTENSIONS:
        base_name = f"{pmid}_Supplement_Text"
        if file_index > 0:
            base_name = f"{pmid}_Supplement_Text_{file_index}"
        return f"{base_name}{suffix}"

    elif suffix in ARCHIVE_EXTENSIONS:
        base_name = f"{pmid}_Supplement_Archive"
        if file_index > 0:
            base_name = f"{pmid}_Supplement_Archive_{file_index}"
        return f"{base_name}{suffix}"

    else:
        # Fallback for other file types
        return f"{pmid}_file_{file_index}{suffix}"


def move_and_rename_files(
    new_files: List[Path],
    pmid: str,
    target_dir: Path
) -> List[Tuple[Path, Path]]:
    """
    Rename and move new files to target directory.

    Args:
        new_files: List of new files to process
        pmid: PubMed ID for naming
        target_dir: Destination directory

    Returns:
        List of (source, destination) tuples for successfully moved files
    """
    moved_files = []
    pdf_count = 0
    spreadsheet_count = 0
    doc_count = 0
    text_count = 0
    archive_count = 0

    # Ensure target directory exists
    target_dir.mkdir(parents=True, exist_ok=True)

    for source_file in new_files:
        suffix = source_file.suffix.lower()

        # Determine file index for multiple files of same type
        if suffix in PDF_EXTENSIONS:
            file_index = pdf_count
            pdf_count += 1
        elif suffix in SPREADSHEET_EXTENSIONS:
            file_index = spreadsheet_count
            spreadsheet_count += 1
        elif suffix in DOC_EXTENSIONS:
            file_index = doc_count
            doc_count += 1
        elif suffix in TEXT_EXTENSIONS:
            file_index = text_count
            text_count += 1
        elif suffix in ARCHIVE_EXTENSIONS:
            file_index = archive_count
            archive_count += 1
        else:
            file_index = 0

        new_name = generate_new_filename(pmid, source_file, file_index)
        dest_file = target_dir / new_name

        # Handle existing files
        if dest_file.exists():
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            stem = dest_file.stem
            suffix = dest_file.suffix
            dest_file = target_dir / f"{stem}_{timestamp}{suffix}"

        try:
            shutil.move(str(source_file), str(dest_file))
            logger.info(f"  Moved: {source_file.name} -> {dest_file.name}")
            moved_files.append((source_file, dest_file))
        except Exception as e:
            logger.error(f"  Failed to move {source_file.name}: {e}")

    return moved_files


def construct_url(
    row: pd.Series,
    doi_column: Optional[str] = None,
    pmid_column: str = "PMID",
) -> str:
    """
    Construct URL for a paper from DOI or PMID.

    Args:
        row: DataFrame row with paper info
        doi_column: Name of DOI column (if present)

    Returns:
        URL string
    """
    # Try DOI first (preferred)
    if doi_column and doi_column in row.index:
        doi = row[doi_column]
        if pd.notna(doi) and str(doi).strip():
            doi = str(doi).strip()
            # Handle various DOI formats
            if doi.startswith('http'):
                return doi
            elif doi.startswith('doi:'):
                doi = doi[4:].strip()
            elif doi.startswith('10.'):
                pass  # Already clean DOI
            return f"https://doi.org/{doi}"

    # Check for URL column (from paywalled_missing.csv)
    if 'URL' in row.index and pd.notna(row['URL']):
        return str(row['URL']).strip()

    # Fall back to PubMed
    pmid_value = row.get(pmid_column, row.get("PMID", ""))
    pmid = str(pmid_value).strip()
    if not pmid or pmid.lower() == "nan":
        logger.warning("Missing PMID value in column '%s'", pmid_column)
        return ""
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def wait_for_user_action() -> str:
    """
    Wait for user input to continue or skip.

    Returns:
        User input string (empty for Enter, 'n' to skip, 'q' to quit)
    """
    print("\n  Press ENTER when done downloading, 'n' to skip, 'q' to quit: ", end='', flush=True)
    try:
        response = input().strip().lower()
        return response
    except EOFError:
        return 'q'
    except KeyboardInterrupt:
        print()
        return 'q'


def run_fetch_manager(
    csv_file: Path,
    downloads_dir: Path = DEFAULT_DOWNLOADS_DIR,
    target_dir: Path = DEFAULT_TARGET_DIR,
    doi_column: Optional[str] = None,
    pmid_column: str = 'PMID',
    status_column: str = 'Status',
    start_index: int = 0,
):
    """
    Main workflow loop for semi-manual paper fetching.

    Args:
        csv_file: Path to CSV file with paper list
        downloads_dir: Path to Downloads folder to monitor
        target_dir: Path to move renamed files to
        doi_column: Name of DOI column (optional)
        pmid_column: Name of PMID column
        status_column: Name of status column for tracking progress
        start_index: Row index to start from (for resuming)
    """
    # Load the CSV
    logger.info(f"Loading CSV: {csv_file}")
    df = pd.read_csv(csv_file)

    # Validate required columns
    if pmid_column not in df.columns:
        logger.error(f"Required column '{pmid_column}' not found in CSV")
        logger.info(f"Available columns: {list(df.columns)}")
        sys.exit(1)

    # Add status column if not present
    if status_column not in df.columns:
        df[status_column] = ''

    # Validate directories
    if not downloads_dir.exists():
        logger.error(f"Downloads directory does not exist: {downloads_dir}")
        sys.exit(1)

    target_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Target directory: {target_dir}")

    # Count pending papers
    df[status_column] = df[status_column].fillna("")
    pending_mask = ~df[status_column].isin(["done", "skipped"])
    total_pending = pending_mask.sum()
    logger.info(f"Found {total_pending} papers to process (out of {len(df)} total)")

    # Process each paper
    processed = 0
    for idx, row in df.iterrows():
        # Skip already done or skipped
        if str(row[status_column]).strip().lower() in {"done", "skipped"}:
            continue

        # Skip if before start index
        if idx < start_index:
            continue

        pmid = str(row.get(pmid_column, "")).strip()
        if not pmid or pmid.lower() == "nan":
            logger.warning(
                "Row %s missing PMID in column '%s'; marking as skipped.",
                idx,
                pmid_column,
            )
            df.at[idx, status_column] = 'skipped'
            df.to_csv(csv_file, index=False)
            continue

        url = construct_url(row, doi_column, pmid_column=pmid_column)

        processed += 1
        remaining = total_pending - processed + 1

        while True:
            print("\n" + "=" * 60)
            print(f"Paper {processed}/{total_pending} (remaining: {remaining})")
            print(f"  PMID: {pmid}")
            print(f"  URL:  {url}")
            if 'Reason' in row.index and pd.notna(row['Reason']):
                print(f"  Note: {row['Reason']}")
            print("=" * 60)

            # Take 'before' snapshot
            before_snapshot = get_downloads_snapshot(downloads_dir)

            # Open URL in browser
            if url:
                logger.info("Opening in browser...")
                webbrowser.open(url)
            else:
                logger.warning("No URL available for PMID %s; skipping.", pmid)
                df.at[idx, status_column] = 'skipped'
                df.to_csv(csv_file, index=False)
                break

            # Wait for user
            response = wait_for_user_action()

            if response == 'q':
                logger.info("Quitting...")
                df.to_csv(csv_file, index=False)
                logger.info(f"Progress saved to {csv_file}")
                return

            if response == 'n':
                logger.info("Skipping this paper...")
                df.at[idx, status_column] = 'skipped'
                df.to_csv(csv_file, index=False)
                break

            # Take 'after' snapshot
            # Small delay to ensure files are fully written
            time.sleep(0.5)
            after_snapshot = get_downloads_snapshot(downloads_dir)

            # Identify new files
            new_files = identify_new_files(before_snapshot, after_snapshot)

            if not new_files:
                logger.warning("No new files detected in Downloads folder")
                print("  No files found. Enter 'r' to retry, 's' to skip, or ENTER to mark done anyway: ", end='', flush=True)
                retry_response = input().strip().lower()

                if retry_response == 'r':
                    continue
                if retry_response == 's':
                    df.at[idx, status_column] = 'skipped'
                    df.to_csv(csv_file, index=False)
                    break
                # Otherwise mark as done (user confirmed no downloads needed)
            else:
                logger.info(f"Found {len(new_files)} new file(s):")
                for f in new_files:
                    logger.info(f"  - {f.name}")

                # Move and rename files
                moved = move_and_rename_files(new_files, pmid, target_dir)

                if moved:
                    logger.info(f"Successfully processed {len(moved)} file(s)")

            # Mark as done and save
            df.at[idx, status_column] = 'done'
            df.to_csv(csv_file, index=False)
            logger.info("Progress saved")
            break

    print("\n" + "=" * 60)
    logger.info("All papers processed!")
    logger.info(f"Files saved to: {target_dir}")
    print("=" * 60)


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Semi-manual literature retrieval workflow assistant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a paywalled_missing.csv from workflow output
  python fetch_manager.py output/GENE/run/pmc_fulltext/paywalled_missing.csv

  # Use a custom CSV with DOI column
  python fetch_manager.py my_papers.csv --doi-column DOI

  # Specify custom Downloads and target directories
  python fetch_manager.py papers.csv --downloads ~/Downloads --target ~/Papers

  # Resume from a specific row (0-indexed)
  python fetch_manager.py papers.csv --start 10
        """
    )

    parser.add_argument(
        "csv_file",
        type=Path,
        help="Path to CSV file containing PMID (and optionally DOI) columns"
    )
    parser.add_argument(
        "--downloads", "-d",
        type=Path,
        default=DEFAULT_DOWNLOADS_DIR,
        help=f"Path to Downloads folder to monitor (default: {DEFAULT_DOWNLOADS_DIR})"
    )
    parser.add_argument(
        "--target-dir", "-t",
        type=Path,
        default=DEFAULT_TARGET_DIR,
        help=f"Target directory for renamed files (default: {DEFAULT_TARGET_DIR})"
    )
    parser.add_argument(
        "--doi-column",
        type=str,
        default=None,
        help="Name of DOI column in CSV (if present, DOI URLs are preferred)"
    )
    parser.add_argument(
        "--pmid-column",
        type=str,
        default="PMID",
        help="Name of PMID column in CSV (default: PMID)"
    )
    parser.add_argument(
        "--status-column",
        type=str,
        default="Status",
        help="Name of status column for progress tracking (default: Status)"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Row index to start from (for resuming interrupted sessions)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate CSV file exists
    if not args.csv_file.exists():
        logger.error(f"CSV file not found: {args.csv_file}")
        sys.exit(1)

    # Run the workflow
    run_fetch_manager(
        csv_file=args.csv_file,
        downloads_dir=args.downloads,
        target_dir=args.target_dir,
        doi_column=args.doi_column,
        pmid_column=args.pmid_column,
        status_column=args.status_column,
        start_index=args.start,
    )


if __name__ == "__main__":
    main()
