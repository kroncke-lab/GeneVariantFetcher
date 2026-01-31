#!/usr/bin/env python3
"""
Standalone CLI for the Data Scout stage.

Processes full-context markdown files through GeneticDataScout to identify
high-value data zones (tables and text sections) containing patient-level
or variant-specific information.

Usage:
    # From a directory of downloaded papers
    python -m cli.scout --input downloads/ --output scout_results/ --gene KCNQ1

    # From a download manifest
    python -m cli.scout --input downloads/manifest.json --output scout_results/ --gene BRCA1

    # With custom manifest output path
    python -m cli.scout --input downloads/ --output results/ --gene SCN5A \
        --manifest-out results/scout_manifest.json
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from pipeline.data_scout import GeneticDataScout
from utils.manifest import Manifest, ManifestEntry, Status, Stage

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def find_input_files(input_path: Path, manifest: Optional[Manifest] = None) -> list[Path]:
    """
    Find input files to process.

    If manifest provided, uses only PMIDs with SUCCESS status.
    Otherwise, globs for *_FULL_CONTEXT.md files.

    Args:
        input_path: Directory containing downloaded papers
        manifest: Optional download manifest to filter by

    Returns:
        List of paths to *_FULL_CONTEXT.md files
    """
    if manifest:
        # Get successful PMIDs from manifest
        successful = manifest.get_successful()
        files = []
        for entry in successful:
            # Look for the FULL_CONTEXT file
            pattern = f"{entry.pmid}_FULL_CONTEXT.md"
            matches = list(input_path.glob(pattern))
            if matches:
                files.append(matches[0])
            else:
                logger.warning(f"PMID {entry.pmid}: No FULL_CONTEXT.md file found")
        return files
    else:
        # Glob for all FULL_CONTEXT files
        return list(input_path.glob("*_FULL_CONTEXT.md"))


def extract_pmid_from_filename(filename: str) -> Optional[str]:
    """Extract PMID from filename like '12345678_FULL_CONTEXT.md'."""
    if "_FULL_CONTEXT.md" in filename:
        return filename.replace("_FULL_CONTEXT.md", "")
    return None


def process_file(
    file_path: Path,
    output_dir: Path,
    scout: GeneticDataScout,
) -> ManifestEntry:
    """
    Process a single FULL_CONTEXT.md file through Data Scout.

    Args:
        file_path: Path to the input markdown file
        output_dir: Directory for output files
        scout: Configured GeneticDataScout instance

    Returns:
        ManifestEntry with processing result
    """
    pmid = extract_pmid_from_filename(file_path.name)
    if not pmid:
        pmid = file_path.stem  # Fallback to filename without extension

    logger.info(f"Processing PMID {pmid}...")

    try:
        # Read the full context markdown
        text = file_path.read_text(encoding="utf-8")

        if not text.strip():
            return ManifestEntry(
                pmid=pmid,
                status=Status.FAILED,
                error_message="Empty file",
            )

        # Run scout analysis
        report = scout.scan(text, pmid=pmid)

        # Generate output files
        files_created = []

        # Condensed markdown with only kept zones
        md_output = output_dir / f"{pmid}_DATA_ZONES.md"
        md_content = scout.format_markdown(report, text)
        md_output.write_text(md_content, encoding="utf-8")
        files_created.append(str(md_output))

        # Full JSON report
        json_output = output_dir / f"{pmid}_DATA_ZONES.json"
        json_content = scout.to_full_json(report)
        json_output.write_text(json_content, encoding="utf-8")
        files_created.append(str(json_output))

        logger.info(
            f"  Found {report.total_zones_found} zones, kept {report.zones_kept} "
            f"({report.compression_ratio:.1%} compression)"
        )

        return ManifestEntry(
            pmid=pmid,
            status=Status.SUCCESS,
            files_created=files_created,
        )

    except Exception as e:
        logger.error(f"  Error processing {pmid}: {e}")
        return ManifestEntry(
            pmid=pmid,
            status=Status.FAILED,
            error_message=str(e),
        )


def run_scout(
    input_path: Path,
    output_dir: Path,
    gene: str,
    manifest_out: Optional[Path] = None,
    min_relevance: float = 0.1,
    max_zones: int = 30,
) -> Manifest:
    """
    Run Data Scout on input files.

    Args:
        input_path: Directory or manifest.json path
        output_dir: Output directory for results
        gene: Gene symbol to search for
        manifest_out: Optional custom manifest output path
        min_relevance: Minimum relevance score threshold
        max_zones: Maximum zones per paper

    Returns:
        Manifest with processing results
    """
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine input mode
    input_manifest = None
    if input_path.is_file() and input_path.name.endswith(".json"):
        # Load manifest and use its directory
        input_manifest = Manifest.load(input_path)
        input_dir = input_path.parent
        logger.info(f"Loaded manifest with {len(input_manifest)} entries")
    elif input_path.is_dir():
        input_dir = input_path
        # Check for manifest.json in directory
        manifest_path = input_dir / "manifest.json"
        if manifest_path.exists():
            input_manifest = Manifest.load(manifest_path)
            logger.info(f"Found manifest with {len(input_manifest)} entries")
    else:
        raise ValueError(f"Input path must be a directory or manifest.json: {input_path}")

    # Find files to process
    files = find_input_files(input_dir, input_manifest)
    if not files:
        logger.warning("No *_FULL_CONTEXT.md files found to process")
        return Manifest(stage=Stage.SCOUT, gene=gene)

    logger.info(f"Found {len(files)} files to process for gene {gene}")

    # Initialize scout
    scout = GeneticDataScout(
        gene_symbol=gene,
        min_relevance_score=min_relevance,
        max_zones=max_zones,
    )

    # Create output manifest
    output_manifest = Manifest(stage=Stage.SCOUT, gene=gene)

    # Process each file
    for file_path in files:
        entry = process_file(file_path, output_dir, scout)
        output_manifest.add_entry(entry)

    # Save manifest
    if manifest_out is None:
        manifest_out = output_dir / "scout_manifest.json"
    output_manifest.save(manifest_out)
    logger.info(f"Saved manifest to {manifest_out}")

    # Print summary
    summary = output_manifest.summary()
    logger.info(f"Scout complete: {summary}")

    return output_manifest


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Data Scout: Identify high-value data zones in downloaded papers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process a directory of papers
    gvf scout --input downloads/ --output scout_results/ --gene KCNQ1

    # Use a download manifest
    gvf scout --input downloads/manifest.json --output results/ --gene BRCA1

    # Customize thresholds
    gvf scout --input papers/ --output results/ --gene SCN5A \\
        --min-relevance 0.2 --max-zones 20
""",
    )

    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input directory or manifest.json file",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output directory for DATA_ZONES files",
    )
    parser.add_argument(
        "--gene", "-g",
        type=str,
        required=True,
        help="Gene symbol to search for (e.g., KCNQ1, BRCA1)",
    )
    parser.add_argument(
        "--manifest-out",
        type=Path,
        default=None,
        help="Custom path for output manifest (default: OUTPUT/scout_manifest.json)",
    )
    parser.add_argument(
        "--min-relevance",
        type=float,
        default=0.1,
        help="Minimum relevance score for TEXT zones (0.0-1.0, default: 0.1)",
    )
    parser.add_argument(
        "--max-zones",
        type=int,
        default=30,
        help="Maximum zones per paper (default: 30)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        manifest = run_scout(
            input_path=args.input,
            output_dir=args.output,
            gene=args.gene,
            manifest_out=args.manifest_out,
            min_relevance=args.min_relevance,
            max_zones=args.max_zones,
        )

        # Exit with error if all failed
        if all(e.status != Status.SUCCESS for e in manifest.entries):
            sys.exit(1)

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
