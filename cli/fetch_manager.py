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
6. (Optional) Converts downloaded files to markdown format for extraction
7. (Optional) Runs Data Scout to create DATA_ZONES.md files

Usage:
    python fetch_manager.py <csv_file> [options]

    # Use a paywalled_missing.csv from a workflow run
    python fetch_manager.py output/GENE/20250101_120000/pmc_fulltext/paywalled_missing.csv

    # Use a custom CSV with PMID and DOI columns
    python fetch_manager.py my_papers.csv --doi-column DOI

    # Specify a custom target directory
    python fetch_manager.py papers.csv --target-dir ~/MyData/Papers

    # Convert downloaded files to markdown and run scout (for re-extraction)
    python fetch_manager.py paywalled_missing.csv --convert --run-scout --gene SCN5A
"""

import logging
import shutil
import sys
import time
import webbrowser
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Set, Tuple

import pandas as pd

# Import converters and scout (optional - for --convert and --run-scout flags)
try:
    from harvesting.format_converters import FormatConverter

    CONVERTER_AVAILABLE = True
except ImportError:
    CONVERTER_AVAILABLE = False
    FormatConverter = None

try:
    from config.settings import get_settings
    from pipeline.data_scout import GeneticDataScout

    SCOUT_AVAILABLE = True
except ImportError:
    SCOUT_AVAILABLE = False
    GeneticDataScout = None
    get_settings = None

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# Default configuration
DEFAULT_DOWNLOADS_DIR = Path.home() / "Downloads"
# Default target directory - use current working directory to avoid hardcoded paths
# Users can override with --target-dir flag
DEFAULT_TARGET_DIR = Path.cwd() / "manual_retrieval"

# File extensions to track
PDF_EXTENSIONS = {".pdf"}
SPREADSHEET_EXTENSIONS = {".xlsx", ".xls", ".csv", ".tsv"}
DOC_EXTENSIONS = {".docx", ".doc"}
TEXT_EXTENSIONS = {".txt"}
ARCHIVE_EXTENSIONS = {".zip", ".rar", ".gz"}
ALL_TRACKED_EXTENSIONS = (
    PDF_EXTENSIONS
    | SPREADSHEET_EXTENSIONS
    | DOC_EXTENSIONS
    | TEXT_EXTENSIONS
    | ARCHIVE_EXTENSIONS
)

# Files to ignore
IGNORED_PATTERNS = {
    ".DS_Store",
    ".crdownload",  # Chrome partial download
    ".part",  # Firefox partial download
    ".tmp",
    ".download",
    "desktop.ini",
    "Thumbs.db",
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
            if item.name.startswith("."):
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
    new_files: List[Path], pmid: str, target_dir: Path
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


def convert_to_markdown(pmid: str, target_dir: Path) -> Optional[str]:
    """
    Convert downloaded files for a PMID to unified markdown format.

    Looks for files matching the PMID in target_dir and creates
    a {PMID}_FULL_CONTEXT.md file combining all content.

    Args:
        pmid: PubMed ID
        target_dir: Directory containing downloaded files

    Returns:
        Path to created markdown file, or None if conversion failed
    """
    if not CONVERTER_AVAILABLE:
        logger.error(
            "FormatConverter not available. Install harvesting module dependencies."
        )
        return None

    import shutil

    converter = FormatConverter()
    markdown_parts = []

    # Find all files for this PMID
    pmid_files = list(target_dir.glob(f"{pmid}_*"))
    if not pmid_files:
        logger.warning(f"No files found for PMID {pmid} in {target_dir}")
        return None

    # Sort files: main text first, then supplements
    main_files = [
        f for f in pmid_files if "Main_Text" in f.name or "main" in f.name.lower()
    ]
    supp_files = [
        f
        for f in pmid_files
        if f not in main_files and f.suffix.lower() not in [".md", ".json"]
    ]

    # Create supplements folder if there are supplement files
    supplements_dir = None
    if supp_files:
        supplements_dir = target_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)

    # Process main text file(s)
    for idx, file_path in enumerate(main_files):
        ext = file_path.suffix.lower()
        if ext == ".pdf":
            content = converter.pdf_to_markdown(file_path)
            markdown_parts.append(f"# MAIN TEXT\n\n{content}")
        elif ext == ".docx":
            content = converter.docx_to_markdown(file_path)
            markdown_parts.append(f"# MAIN TEXT\n\n{content}")
        elif ext == ".doc":
            content = converter.doc_to_markdown(file_path)
            markdown_parts.append(f"# MAIN TEXT\n\n{content}")
        elif ext in [".txt", ".md"]:
            content = file_path.read_text(encoding="utf-8", errors="ignore")
            markdown_parts.append(f"# MAIN TEXT\n\n{content}")

    # Process supplemental files
    for idx, file_path in enumerate(supp_files, 1):
        ext = file_path.suffix.lower()
        content = ""

        if ext == ".pdf":
            content = converter.pdf_to_markdown(file_path)
        elif ext in [".xlsx", ".xls"]:
            content = converter.excel_to_markdown(file_path)
        elif ext == ".docx":
            content = converter.docx_to_markdown(file_path)
        elif ext == ".doc":
            content = converter.doc_to_markdown(file_path)
        elif ext in [".txt", ".csv"]:
            try:
                content = file_path.read_text(encoding="utf-8", errors="ignore")
            except Exception as e:
                content = f"[Error reading file: {e}]"

        if content:
            markdown_parts.append(
                f"\n\n# SUPPLEMENTAL FILE {idx}: {file_path.name}\n\n{content}"
            )

        # Move supplement file to supplements folder
        if supplements_dir and file_path.exists():
            dest = supplements_dir / file_path.name
            try:
                shutil.move(str(file_path), str(dest))
                logger.debug(f"Moved {file_path.name} to {pmid}_supplements/")
            except Exception as e:
                logger.warning(f"Could not move {file_path.name}: {e}")

    if not markdown_parts:
        logger.warning(f"No convertible content found for PMID {pmid}")
        return None

    # Write unified markdown
    unified_content = "\n".join(markdown_parts)
    output_file = target_dir / f"{pmid}_FULL_CONTEXT.md"
    output_file.write_text(unified_content, encoding="utf-8")

    logger.info(f"Created {output_file.name} ({len(unified_content)} characters)")
    if supplements_dir:
        logger.info(f"Organized {len(supp_files)} supplements into {pmid}_supplements/")
    return str(output_file)


def run_data_scout(pmid: str, target_dir: Path, gene_symbol: str) -> Optional[str]:
    """
    Run the Genetic Data Scout on a FULL_CONTEXT.md file.

    Args:
        pmid: PubMed ID
        target_dir: Directory containing the FULL_CONTEXT.md file
        gene_symbol: Gene symbol for relevance scoring

    Returns:
        Path to created DATA_ZONES.md file, or None if scout failed
    """
    if not SCOUT_AVAILABLE:
        logger.error(
            "GeneticDataScout not available. Install pipeline module dependencies."
        )
        return None

    full_context_file = target_dir / f"{pmid}_FULL_CONTEXT.md"
    if not full_context_file.exists():
        logger.warning(f"No FULL_CONTEXT.md found for PMID {pmid}")
        return None

    try:
        settings = get_settings()
        scout = GeneticDataScout(
            gene_symbol=gene_symbol,
            min_relevance_score=settings.scout_min_relevance if settings else 0.3,
            max_zones=settings.scout_max_zones if settings else 30,
        )

        unified_content = full_context_file.read_text(encoding="utf-8")
        report = scout.scan(unified_content, pmid=pmid)

        # Write zone metadata JSON
        zones_json_path = target_dir / f"{pmid}_DATA_ZONES.json"
        zones_json_path.write_text(scout.to_json(report))

        # Write condensed markdown
        zones_md_path = target_dir / f"{pmid}_DATA_ZONES.md"
        zones_md_path.write_text(scout.format_markdown(report, unified_content))

        logger.info(
            f"Scout: {report.zones_kept}/{report.total_zones_found} zones kept "
            f"({report.compression_ratio:.0%} compression)"
        )
        return str(zones_md_path)

    except Exception as e:
        logger.error(f"Data scout failed for PMID {pmid}: {e}")
        return None


def process_completed_downloads(
    csv_file: Path,
    target_dir: Path,
    gene_symbol: Optional[str] = None,
    convert: bool = False,
    run_scout: bool = False,
    pmid_column: str = "PMID",
    status_column: str = "Status",
) -> Tuple[int, int]:
    """
    Process completed downloads: convert to markdown and/or run scout.

    Args:
        csv_file: Path to CSV file with paper list
        target_dir: Directory containing downloaded files
        gene_symbol: Gene symbol for scout (required if run_scout=True)
        convert: Whether to convert files to markdown
        run_scout: Whether to run data scout
        pmid_column: Name of PMID column
        status_column: Name of status column

    Returns:
        Tuple of (converted_count, scouted_count)
    """
    df = pd.read_csv(csv_file)

    if status_column not in df.columns:
        logger.error(f"Status column '{status_column}' not found in CSV")
        return 0, 0

    # Find completed downloads
    done_mask = df[status_column].str.lower() == "done"
    done_pmids = df.loc[done_mask, pmid_column].astype(str).tolist()

    if not done_pmids:
        logger.info("No completed downloads to process")
        return 0, 0

    logger.info(f"Processing {len(done_pmids)} completed downloads...")

    converted_count = 0
    scouted_count = 0

    for pmid in done_pmids:
        pmid = pmid.strip()
        if not pmid or pmid.lower() == "nan":
            continue

        # Convert to markdown if requested
        if convert:
            result = convert_to_markdown(pmid, target_dir)
            if result:
                converted_count += 1

        # Run scout if requested
        if run_scout:
            if not gene_symbol:
                logger.warning("Gene symbol required for scout. Use --gene flag.")
                break
            result = run_data_scout(pmid, target_dir, gene_symbol)
            if result:
                scouted_count += 1

    return converted_count, scouted_count


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
            if doi.startswith("http"):
                return doi
            elif doi.startswith("doi:"):
                doi = doi[4:].strip()
            elif doi.startswith("10."):
                pass  # Already clean DOI
            else:
                # DOI doesn't start with expected prefix, try anyway
                pass
            return f"https://doi.org/{doi}"

    # Check for URL column (from paywalled_missing.csv)
    if "URL" in row.index and pd.notna(row["URL"]):
        return str(row["URL"]).strip()

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
    print(
        "\n  Press ENTER when done downloading, 'n' to skip, 'q' to quit: ",
        end="",
        flush=True,
    )
    try:
        response = input().strip().lower()
        return response
    except EOFError:
        return "q"
    except KeyboardInterrupt:
        print()
        return "q"


def run_fetch_manager(
    csv_file: Path,
    downloads_dir: Path = DEFAULT_DOWNLOADS_DIR,
    target_dir: Optional[Path] = None,
    doi_column: Optional[str] = None,
    pmid_column: str = "PMID",
    status_column: str = "Status",
    start_index: int = 0,
    convert: bool = False,
    run_scout: bool = False,
    gene_symbol: Optional[str] = None,
    process_only: bool = False,
):
    """
    Main workflow loop for semi-manual paper fetching.

    Args:
        csv_file: Path to CSV file with paper list
        downloads_dir: Path to Downloads folder to monitor
        target_dir: Path to move renamed files to (auto-detected from CSV if None)
        doi_column: Name of DOI column (optional)
        pmid_column: Name of PMID column
        status_column: Name of status column for tracking progress
        start_index: Row index to start from (for resuming)
        convert: Whether to convert downloaded files to markdown
        run_scout: Whether to run data scout on converted files
        gene_symbol: Gene symbol for scout (required if run_scout=True)
        process_only: Skip download loop, just process existing completed downloads
    """
    # Auto-detect target directory from CSV path if not specified
    # If CSV is in pmc_fulltext/, use that directory
    if target_dir is None:
        csv_parent = csv_file.parent
        if csv_parent.name == "pmc_fulltext" or "pmc_fulltext" in str(csv_parent):
            target_dir = csv_parent
            logger.info(f"Auto-detected target directory: {target_dir}")
        else:
            target_dir = DEFAULT_TARGET_DIR

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
        df[status_column] = ""

    # Validate directories
    if not downloads_dir.exists():
        logger.error(f"Downloads directory does not exist: {downloads_dir}")
        sys.exit(1)

    target_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Target directory: {target_dir}")

    # Validate scout requirements
    if run_scout and not gene_symbol:
        logger.error("--gene is required when using --run-scout")
        sys.exit(1)

    if run_scout and not SCOUT_AVAILABLE:
        logger.error("GeneticDataScout not available. Install pipeline dependencies.")
        sys.exit(1)

    if convert and not CONVERTER_AVAILABLE:
        logger.error("FormatConverter not available. Install harvesting dependencies.")
        sys.exit(1)

    # Count pending papers
    df[status_column] = df[status_column].fillna("")
    pending_mask = ~df[status_column].isin(["done", "skipped"])
    total_pending = pending_mask.sum()
    logger.info(f"Found {total_pending} papers to process (out of {len(df)} total)")

    # If process_only mode, skip download loop and just convert/scout
    if process_only:
        logger.info("Process-only mode: skipping download loop")
        if convert or run_scout:
            converted, scouted = process_completed_downloads(
                csv_file=csv_file,
                target_dir=target_dir,
                gene_symbol=gene_symbol,
                convert=convert,
                run_scout=run_scout,
                pmid_column=pmid_column,
                status_column=status_column,
            )
            print(f"\nProcessing complete: {converted} converted, {scouted} scouted")
        return

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
            df.at[idx, status_column] = "skipped"
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
            if "Reason" in row.index and pd.notna(row["Reason"]):
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
                df.at[idx, status_column] = "skipped"
                df.to_csv(csv_file, index=False)
                break

            # Wait for user
            response = wait_for_user_action()

            if response == "q":
                logger.info("Quitting...")
                df.to_csv(csv_file, index=False)
                logger.info(f"Progress saved to {csv_file}")
                return

            if response == "n":
                logger.info("Skipping this paper...")
                df.at[idx, status_column] = "skipped"
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
                print(
                    "  No files found. Enter 'r' to retry, 's' to skip, or ENTER to mark done anyway: ",
                    end="",
                    flush=True,
                )
                retry_response = input().strip().lower()

                if retry_response == "r":
                    continue
                if retry_response == "s":
                    df.at[idx, status_column] = "skipped"
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
            df.at[idx, status_column] = "done"
            df.to_csv(csv_file, index=False)
            logger.info("Progress saved")
            break

    print("\n" + "=" * 60)
    logger.info("All papers processed!")
    logger.info(f"Files saved to: {target_dir}")
    print("=" * 60)

    # Post-processing: convert and/or scout
    if convert or run_scout:
        print("\n" + "=" * 60)
        logger.info("Running post-download processing...")
        converted, scouted = process_completed_downloads(
            csv_file=csv_file,
            target_dir=target_dir,
            gene_symbol=gene_symbol,
            convert=convert,
            run_scout=run_scout,
            pmid_column=pmid_column,
            status_column=status_column,
        )
        print(f"Post-processing complete: {converted} converted, {scouted} scouted")
        print("=" * 60)

        if convert:
            print("\nNext steps:")
            print("  - Re-run extraction using GUI 'folder job' feature")
            print("  - Or run: python automated_workflow.py ... (with existing folder)")


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

  # Convert downloaded files to markdown and run scout (for re-extraction)
  python fetch_manager.py paywalled_missing.csv --convert --run-scout --gene SCN5A

  # Only process already-completed downloads (no download loop)
  python fetch_manager.py paywalled_missing.csv --process-only --convert --run-scout --gene TTR
        """,
    )

    parser.add_argument(
        "csv_file",
        type=Path,
        help="Path to CSV file containing PMID (and optionally DOI) columns",
    )
    parser.add_argument(
        "--downloads",
        "-d",
        type=Path,
        default=DEFAULT_DOWNLOADS_DIR,
        help=f"Path to Downloads folder to monitor (default: {DEFAULT_DOWNLOADS_DIR})",
    )
    parser.add_argument(
        "--target-dir",
        "-t",
        type=Path,
        default=None,
        help="Target directory for renamed files (auto-detected from CSV path if in pmc_fulltext/)",
    )
    parser.add_argument(
        "--doi-column",
        type=str,
        default=None,
        help="Name of DOI column in CSV (if present, DOI URLs are preferred)",
    )
    parser.add_argument(
        "--pmid-column",
        type=str,
        default="PMID",
        help="Name of PMID column in CSV (default: PMID)",
    )
    parser.add_argument(
        "--status-column",
        type=str,
        default="Status",
        help="Name of status column for progress tracking (default: Status)",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Row index to start from (for resuming interrupted sessions)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "--convert",
        "-c",
        action="store_true",
        help="Convert downloaded files to FULL_CONTEXT.md format for extraction",
    )
    parser.add_argument(
        "--run-scout",
        action="store_true",
        help="Run Data Scout to create DATA_ZONES.md files (requires --gene)",
    )
    parser.add_argument(
        "--gene",
        "-g",
        type=str,
        default=None,
        help="Gene symbol for Data Scout relevance scoring (required with --run-scout)",
    )
    parser.add_argument(
        "--process-only",
        action="store_true",
        help="Skip download loop, only process already-completed downloads",
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
        convert=args.convert,
        run_scout=args.run_scout,
        gene_symbol=args.gene,
        process_only=args.process_only,
    )


if __name__ == "__main__":
    main()
