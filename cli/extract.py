#!/usr/bin/env python3
"""
Standalone CLI for the Extraction stage.

Processes DATA_ZONES.md files (from Data Scout) or FULL_CONTEXT.md files
through ExpertExtractor to extract structured genetic variant data.

Usage:
    # From a directory of scout outputs (looks for scout_manifest.json)
    python -m cli.extract --input scout_results/ --output extractions/ --gene KCNQ1

    # From a scout manifest directly
    python -m cli.extract --input scout_results/scout_manifest.json --output extractions/ --gene BRCA1

    # With custom output manifest path
    python -m cli.extract --input scout_results/ --output extractions/ --gene SCN5A \
        --manifest-out extractions/extraction_manifest.json

    # Force full-text mode (skip DATA_ZONES, use FULL_CONTEXT.md)
    python -m cli.extract --input downloads/ --output extractions/ --gene KCNH2 --full-text
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Optional

from pipeline.extraction import ExpertExtractor
from utils.manifest import Manifest, ManifestEntry, Stage, Status
from utils.models import Paper

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# =============================================================================
# INPUT VALIDATION
# =============================================================================

# Gene symbol pattern: 1-10 uppercase letters/numbers, optionally followed by orf suffix
GENE_SYMBOL_PATTERN = re.compile(r"^[A-Z][A-Z0-9]{0,9}(?:orf\d+)?$", re.IGNORECASE)


class ValidationError(Exception):
    """Raised when input validation fails."""

    pass


def validate_gene_symbol(gene: str) -> None:
    """
    Validate gene symbol format.

    Args:
        gene: Gene symbol to validate

    Raises:
        ValidationError: If gene symbol is invalid
    """
    if not gene:
        raise ValidationError("Gene symbol is required but was not provided")

    if not GENE_SYMBOL_PATTERN.match(gene):
        raise ValidationError(
            f"Invalid gene symbol format: '{gene}'. "
            f"Expected format: 1-10 alphanumeric characters starting with a letter "
            f"(e.g., BRCA1, TP53, SCN5A)"
        )


def validate_input_path(input_path: Path) -> None:
    """
    Validate that input path exists and is accessible.

    Args:
        input_path: Path to input directory or manifest file

    Raises:
        ValidationError: If input path is invalid
    """
    if not input_path.exists():
        raise ValidationError(f"Input path does not exist: {input_path}")

    if input_path.is_file():
        if not input_path.name.endswith(".json"):
            raise ValidationError(f"Input file must be a JSON manifest: {input_path}")
        # Check file is readable
        if not os.access(input_path, os.R_OK):
            raise ValidationError(f"Input file is not readable: {input_path}")
    elif input_path.is_dir():
        # Check directory is readable
        if not os.access(input_path, os.R_OK):
            raise ValidationError(f"Input directory is not readable: {input_path}")
    else:
        raise ValidationError(f"Input path must be a file or directory: {input_path}")


def validate_output_directory(output_dir: Path) -> None:
    """
    Validate that output directory is writable.

    Args:
        output_dir: Path to output directory

    Raises:
        ValidationError: If output directory is not writable
    """
    if output_dir.exists():
        if not output_dir.is_dir():
            raise ValidationError(
                f"Output path exists but is not a directory: {output_dir}"
            )
        if not os.access(output_dir, os.W_OK):
            raise ValidationError(f"Output directory is not writable: {output_dir}")
    else:
        # Check if parent directory is writable (so we can create output_dir)
        parent = output_dir.parent
        if not parent.exists():
            raise ValidationError(
                f"Parent directory does not exist: {parent}. "
                f"Cannot create output directory: {output_dir}"
            )
        if not os.access(parent, os.W_OK):
            raise ValidationError(
                f"Cannot create output directory (parent not writable): {output_dir}"
            )


def validate_manifest_has_success(manifest: Manifest, manifest_path: Path) -> None:
    """
    Validate that a manifest has at least one SUCCESS entry.

    Args:
        manifest: Loaded manifest object
        manifest_path: Path to manifest (for error message)

    Raises:
        ValidationError: If manifest has no SUCCESS entries
    """
    successful = manifest.get_successful()
    if not successful:
        failed_count = len(manifest.get_failed())
        total = len(manifest.entries)
        raise ValidationError(
            f"Manifest has no successful entries to process: {manifest_path}. "
            f"({total} total entries, {failed_count} failed)"
        )


def validate_extract_inputs(
    input_path: Path,
    output_dir: Path,
    gene: str,
) -> None:
    """
    Validate all inputs for the extract stage.

    Args:
        input_path: Path to input directory or manifest
        output_dir: Path to output directory
        gene: Gene symbol

    Raises:
        ValidationError: If any validation fails
    """
    validate_gene_symbol(gene)
    validate_input_path(input_path)
    validate_output_directory(output_dir)


def find_input_files(
    input_path: Path,
    manifest: Optional[Manifest] = None,
    use_full_text: bool = False,
) -> list[tuple[Path, str]]:
    """
    Find input files to process.

    If manifest provided, uses file paths from successful entries.
    Otherwise, globs for *_DATA_ZONES.md or *_FULL_CONTEXT.md files.

    Args:
        input_path: Directory containing scout outputs or downloaded papers
        manifest: Optional scout manifest to get file paths from
        use_full_text: If True, look for FULL_CONTEXT.md instead of DATA_ZONES.md

    Returns:
        List of (path, pmid) tuples
    """
    results = []

    if manifest:
        # Get successful entries from manifest
        successful = manifest.get_successful()
        for entry in successful:
            # Use files_created from manifest entry if available
            if entry.files_created:
                for file_path_str in entry.files_created:
                    file_path = Path(file_path_str)
                    # Prefer DATA_ZONES.md files unless full_text mode
                    if use_full_text:
                        if "_FULL_CONTEXT.md" in file_path.name:
                            results.append((file_path, entry.pmid))
                            break
                    else:
                        if "_DATA_ZONES.md" in file_path.name:
                            results.append((file_path, entry.pmid))
                            break
                else:
                    # Fallback: construct expected path
                    if use_full_text:
                        pattern = f"{entry.pmid}_FULL_CONTEXT.md"
                    else:
                        pattern = f"{entry.pmid}_DATA_ZONES.md"
                    matches = list(input_path.glob(pattern))
                    if matches:
                        results.append((matches[0], entry.pmid))
                    else:
                        logger.warning(f"PMID {entry.pmid}: No input file found")
            else:
                # No files_created, try to find by PMID
                if use_full_text:
                    pattern = f"{entry.pmid}_FULL_CONTEXT.md"
                else:
                    pattern = f"{entry.pmid}_DATA_ZONES.md"
                matches = list(input_path.glob(pattern))
                if matches:
                    results.append((matches[0], entry.pmid))
                else:
                    logger.warning(f"PMID {entry.pmid}: No input file found")
        return results
    else:
        # Glob for files
        if use_full_text:
            pattern = "*_FULL_CONTEXT.md"
        else:
            pattern = "*_DATA_ZONES.md"

        for file_path in input_path.glob(pattern):
            pmid = extract_pmid_from_filename(file_path.name, use_full_text)
            if pmid:
                results.append((file_path, pmid))

        return results


def extract_pmid_from_filename(
    filename: str, use_full_text: bool = False
) -> Optional[str]:
    """Extract PMID from filename like '12345678_DATA_ZONES.md' or '12345678_FULL_CONTEXT.md'."""
    if use_full_text:
        if "_FULL_CONTEXT.md" in filename:
            return filename.replace("_FULL_CONTEXT.md", "")
    else:
        if "_DATA_ZONES.md" in filename:
            return filename.replace("_DATA_ZONES.md", "")
    return None


def process_file(
    file_path: Path,
    pmid: str,
    output_dir: Path,
    gene_symbol: str,
    extractor: ExpertExtractor,
    use_full_text: bool = False,
) -> tuple[ManifestEntry, Optional[dict]]:
    """
    Process a single file through ExpertExtractor.

    Args:
        file_path: Path to the input markdown file
        pmid: PubMed ID for this paper
        output_dir: Directory for output files
        gene_symbol: Target gene symbol
        extractor: Configured ExpertExtractor instance
        use_full_text: Whether input is FULL_CONTEXT vs DATA_ZONES

    Returns:
        Tuple of (ManifestEntry, extracted_data or None)
    """
    logger.info(f"Processing PMID {pmid}...")

    try:
        # Read the input file
        text = file_path.read_text(encoding="utf-8")

        if not text.strip():
            return ManifestEntry(
                pmid=pmid,
                status=Status.FAILED,
                error_message="Empty file",
            ), None

        # Create Paper object for extraction
        paper = Paper(
            pmid=pmid,
            gene_symbol=gene_symbol,
            title=f"Paper {pmid}",  # We don't have title metadata here
            full_text=text,
        )

        # Run extraction
        result = extractor.extract(paper)

        if not result.success:
            logger.warning(f"  Extraction failed: {result.error}")
            return ManifestEntry(
                pmid=pmid,
                status=Status.FAILED,
                error_message=result.error or "Extraction failed",
            ), None

        # Count variants
        extracted_data = result.extracted_data or {}
        variants = extracted_data.get("variants", [])
        num_variants = len(variants)

        # Generate output file
        files_created = []
        output_file = output_dir / f"{pmid}_extraction.json"

        # Add model info to metadata
        if "extraction_metadata" not in extracted_data:
            extracted_data["extraction_metadata"] = {}
        extracted_data["extraction_metadata"]["model_used"] = result.model_used
        extracted_data["extraction_metadata"]["source_file"] = str(file_path)

        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(extracted_data, f, indent=2)
        files_created.append(str(output_file))

        logger.info(f"  ✓ Extracted {num_variants} variants using {result.model_used}")

        # Create manifest entry with extended info
        entry = ManifestEntry(
            pmid=pmid,
            status=Status.SUCCESS,
            files_created=files_created,
        )

        return entry, extracted_data

    except Exception as e:
        logger.error(f"  Error processing {pmid}: {e}")
        return ManifestEntry(
            pmid=pmid,
            status=Status.FAILED,
            error_message=str(e),
        ), None


def run_extraction(
    input_path: Path,
    output_dir: Path,
    gene: str,
    manifest_out: Optional[Path] = None,
    models: Optional[list[str]] = None,
    tier_threshold: int = 1,
    use_full_text: bool = False,
) -> Manifest:
    """
    Run extraction on input files.

    Args:
        input_path: Directory or manifest.json path
        output_dir: Output directory for extraction results
        gene: Gene symbol being analyzed
        manifest_out: Optional custom manifest output path
        models: Optional list of models to use (overrides config)
        tier_threshold: Variant threshold for trying next model
        use_full_text: If True, use FULL_CONTEXT.md instead of DATA_ZONES.md

    Returns:
        Manifest with processing results

    Raises:
        ValidationError: If input validation fails
    """
    # Validate inputs before processing
    validate_extract_inputs(input_path, output_dir, gene)

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine input mode
    input_manifest = None
    if input_path.is_file() and input_path.name.endswith(".json"):
        # Load manifest and use its directory
        input_manifest = Manifest.load(input_path)
        input_dir = input_path.parent
        logger.info(
            f"Loaded manifest with {len(input_manifest)} entries (stage: {input_manifest.stage.value})"
        )
        # Validate manifest has successful entries
        validate_manifest_has_success(input_manifest, input_path)
    elif input_path.is_dir():
        input_dir = input_path
        # Check for scout_manifest.json or manifest.json in directory
        scout_manifest_path = input_dir / "scout_manifest.json"
        manifest_path = input_dir / "manifest.json"

        if scout_manifest_path.exists():
            input_manifest = Manifest.load(scout_manifest_path)
            logger.info(f"Found scout_manifest.json with {len(input_manifest)} entries")
            # Validate manifest has successful entries
            validate_manifest_has_success(input_manifest, scout_manifest_path)
        elif manifest_path.exists():
            input_manifest = Manifest.load(manifest_path)
            logger.info(f"Found manifest.json with {len(input_manifest)} entries")
            # Validate manifest has successful entries
            validate_manifest_has_success(input_manifest, manifest_path)
    else:
        raise ValidationError(
            f"Input path must be a directory or manifest.json: {input_path}"
        )

    # Find files to process
    files = find_input_files(input_dir, input_manifest, use_full_text)
    if not files:
        file_type = "FULL_CONTEXT.md" if use_full_text else "DATA_ZONES.md"
        logger.warning(f"No *_{file_type} files found to process")
        return Manifest(stage=Stage.EXTRACT, gene=gene)

    logger.info(f"Found {len(files)} files to process for gene {gene}")

    # Initialize extractor
    extractor = ExpertExtractor(
        models=models,
        tier_threshold=tier_threshold,
        fulltext_dir=str(input_dir),
    )

    # Create output manifest
    output_manifest = Manifest(stage=Stage.EXTRACT, gene=gene)

    # Track variant statistics
    total_variants = 0
    variants_by_pmid = {}

    # Process each file
    for file_path, pmid in files:
        entry, extracted_data = process_file(
            file_path=file_path,
            pmid=pmid,
            output_dir=output_dir,
            gene_symbol=gene,
            extractor=extractor,
            use_full_text=use_full_text,
        )
        output_manifest.add_entry(entry)

        if extracted_data:
            num_variants = len(extracted_data.get("variants", []))
            total_variants += num_variants
            variants_by_pmid[pmid] = num_variants

    # Save manifest with extended metadata
    if manifest_out is None:
        manifest_out = output_dir / "extraction_manifest.json"

    # Save the manifest
    output_manifest.save(manifest_out)
    logger.info(f"Saved manifest to {manifest_out}")

    # Also save a summary JSON with variant counts
    summary = {
        "gene": gene,
        "total_papers_processed": len(files),
        "successful_extractions": len(output_manifest.get_successful()),
        "failed_extractions": len(output_manifest.get_failed()),
        "total_variants_extracted": total_variants,
        "variants_by_pmid": variants_by_pmid,
    }
    summary_file = output_dir / "extraction_summary.json"
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Saved summary to {summary_file}")

    # Print summary
    manifest_summary = output_manifest.summary()
    logger.info(f"Extraction complete: {manifest_summary}")
    logger.info(f"Total variants extracted: {total_variants}")

    return output_manifest


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Extract structured variant data from papers using AI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process scout outputs
    gvf extract --input scout_results/ --output extractions/ --gene KCNQ1

    # Use a scout manifest
    gvf extract --input scout_results/scout_manifest.json --output extractions/ --gene BRCA1

    # Force full-text mode (bypass DATA_ZONES)
    gvf extract --input downloads/ --output extractions/ --gene SCN5A --full-text

    # Customize model selection
    gvf extract --input scout_results/ --output extractions/ --gene KCNH2 \\
        --model openai/gpt-4o --tier-threshold 5
""",
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Input directory or manifest.json file (scout_manifest.json or manifest.json)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output directory for extraction JSON files",
    )
    parser.add_argument(
        "--gene",
        "-g",
        type=str,
        required=True,
        help="Gene symbol being analyzed (e.g., KCNQ1, BRCA1)",
    )
    parser.add_argument(
        "--manifest-out",
        type=Path,
        default=None,
        help="Custom path for output manifest (default: OUTPUT/extraction_manifest.json)",
    )
    parser.add_argument(
        "--model",
        type=str,
        action="append",
        dest="models",
        default=None,
        help="Model to use for extraction (can be repeated for tiered fallback)",
    )
    parser.add_argument(
        "--tier-threshold",
        type=int,
        default=1,
        help="If first model finds fewer variants, try next model (default: 1)",
    )
    parser.add_argument(
        "--full-text",
        action="store_true",
        help="Use FULL_CONTEXT.md files instead of DATA_ZONES.md",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        manifest = run_extraction(
            input_path=args.input,
            output_dir=args.output,
            gene=args.gene,
            manifest_out=args.manifest_out,
            models=args.models,
            tier_threshold=args.tier_threshold,
            use_full_text=args.full_text,
        )

        # Exit with error if all failed
        if manifest.entries and all(
            e.status != Status.SUCCESS for e in manifest.entries
        ):
            sys.exit(1)

    except ValidationError as e:
        logger.error(f"Input validation failed: {e}")
        sys.exit(2)
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
