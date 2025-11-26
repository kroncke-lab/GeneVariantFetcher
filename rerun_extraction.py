#!/usr/bin/env python3
"""
Re-run LLM extraction on already-downloaded papers.

This script processes existing markdown files from a previous workflow run
and re-extracts variant data using the current extraction pipeline.

Usage:
    python3.12 rerun_extraction.py <workflow_output_folder> <gene_symbol>

Example:
    python3.12 rerun_extraction.py automated_output/TTR/20251125_114028 TTR
"""

import os
import sys
import json
import logging
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def rerun_extraction_on_folder(
    workflow_folder: str,
    gene_symbol: str,
    tier_threshold: int = 1,
    max_workers: int = 5
):
    """
    Re-run extraction on all markdown files in a workflow output folder.

    Args:
        workflow_folder: Path to the workflow output folder (e.g., automated_output/TTR/20251125_114028)
        gene_symbol: Gene symbol for filtering variants (e.g., "TTR")
        tier_threshold: If the first model finds fewer variants than this, the next model is tried.
        max_workers: Number of parallel workers for processing papers
    """
    from models import Paper
    from pipeline.extraction import ExpertExtractor
    from pipeline.aggregation import aggregate_penetrance

    workflow_path = Path(workflow_folder)
    if not workflow_path.exists():
        logger.error(f"Folder not found: {workflow_folder}")
        return {"success": False, "error": "Folder not found"}

    logger.info("="*80)
    logger.info(f"RE-RUNNING EXTRACTION FOR GENE: {gene_symbol}")
    logger.info(f"Source folder: {workflow_folder}")
    logger.info("="*80)

    # Find markdown files in pmc_fulltext subfolder
    fulltext_dir = workflow_path / "pmc_fulltext"
    if not fulltext_dir.exists():
        logger.error(f"Full-text directory not found: {fulltext_dir}")
        return {"success": False, "error": "Full-text directory not found"}

    markdown_files = list(fulltext_dir.glob("*_FULL_CONTEXT.md"))
    logger.info(f"Found {len(markdown_files)} markdown files to process")

    if not markdown_files:
        logger.warning("No markdown files found to process")
        return {"success": False, "error": "No markdown files found"}

    # Create new extraction output directory
    extraction_dir = workflow_path / "extractions_rerun" / datetime.now().strftime("%Y%m%d_%H%M%S")
    extraction_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving new extractions to: {extraction_dir}")

    # Initialize extractor
    extractor = ExpertExtractor(tier_threshold=tier_threshold)

    def process_paper_file(md_file):
        """Process a single paper file (for parallel execution)"""
        # Extract PMID from filename (format: PMID_12345678_FULL_CONTEXT.md)
        filename_parts = md_file.stem.split('_')
        pmid_match = filename_parts[1] if len(filename_parts) > 1 and filename_parts[0] == "PMID" else filename_parts[0]

        if not pmid_match:
            logger.warning(f"Could not extract PMID from filename: {md_file.name}")
            return None

        logger.info(f"Processing PMID {pmid_match}...")

        # Read the markdown content
        with open(md_file, 'r', encoding='utf-8') as f:
            full_text = f.read()

        # Extract title from markdown if available (first line after "# ")
        title = f"Paper {pmid_match}"
        for line in full_text.split('\n')[:20]:
            if line.startswith('# '):
                title = line[2:].strip()
                break

        # Create paper object with gene_symbol
        paper = Paper(
            pmid=pmid_match,
            title=title,
            full_text=full_text,
            gene_symbol=gene_symbol
        )

        # Extract variant data
        result = extractor.extract(paper)

        if result.success:
            # Save extraction result
            output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid_match}.json"
            with open(output_file, 'w') as f:
                json.dump(result.extracted_data, f, indent=2)
            logger.info(f"✓ PMID {pmid_match} - Saved extraction to {output_file.name}")
            return result
        else:
            logger.warning(f"✗ PMID {pmid_match} - Extraction failed: {result.error}")
            return None

    # Process papers in parallel
    extractions = []
    successful_count = 0
    failed_count = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {executor.submit(process_paper_file, md_file): md_file for md_file in markdown_files}

        for future in as_completed(future_to_file):
            result = future.result()
            if result:
                extractions.append(result)
                successful_count += 1
            else:
                failed_count += 1

    logger.info(f"\n✓ Extraction complete: {successful_count} successful, {failed_count} failed")

    # ============================================================================
    # STEP 4: Aggregate Penetrance Data
    # ============================================================================
    logger.info("\n📊 STEP 4: Aggregating penetrance data...")

    # Collect all extraction JSON files
    extraction_files = list(extraction_dir.glob(f"{gene_symbol}_PMID_*.json"))
    logger.info(f"Found {len(extraction_files)} extraction files to aggregate")

    all_extracted_data = []
    for json_file in extraction_files:
        with open(json_file, 'r') as f:
            data = json.load(f)
            all_extracted_data.append(data)

    # Aggregate penetrance data
    penetrance_summary = aggregate_penetrance(
        extracted_papers=all_extracted_data,
        gene_symbol=gene_symbol
    )

    # Save aggregated results
    summary_file = workflow_path / f"{gene_symbol}_penetrance_summary_rerun_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(summary_file, 'w') as f:
        json.dump(penetrance_summary, f, indent=2)

    logger.info(f"✓ Penetrance summary saved to: {summary_file}")

    # Print summary statistics
    total_variants = penetrance_summary.get('total_variants', 0)
    variants_with_penetrance = sum(
        1 for v in penetrance_summary.get('variants', [])
        if v.get('aggregated_penetrance', {}).get('total_carriers', 0) > 0
    )
    total_carriers = sum(
        v.get('aggregated_penetrance', {}).get('total_carriers', 0)
        for v in penetrance_summary.get('variants', [])
    )
    total_affected = sum(
        v.get('aggregated_penetrance', {}).get('affected', 0)
        for v in penetrance_summary.get('variants', [])
    )

    logger.info("\n" + "="*80)
    logger.info("RE-EXTRACTION COMPLETE!")
    logger.info("="*80)
    logger.info(f"Gene: {gene_symbol}")
    logger.info(f"Papers processed: {len(markdown_files)}")
    logger.info(f"Papers with extractions: {successful_count}")
    logger.info(f"Total variants found: {total_variants}")
    logger.info(f"Variants with penetrance data: {variants_with_penetrance}")
    logger.info(f"Total carriers observed: {total_carriers}")
    logger.info(f"Total affected carriers: {total_affected}")
    logger.info(f"\nNew extractions saved to: {extraction_dir}")
    logger.info(f"Penetrance summary: {summary_file}")
    logger.info("="*80)

    return {
        "success": True,
        "papers_processed": len(markdown_files),
        "successful_extractions": successful_count,
        "failed_extractions": failed_count,
        "total_variants": total_variants,
        "variants_with_penetrance": variants_with_penetrance,
        "total_carriers": total_carriers,
        "total_affected": total_affected,
        "extraction_dir": str(extraction_dir),
        "summary_file": str(summary_file)
    }


def main():
    """Command-line entrypoint."""
    if len(sys.argv) < 3:
        print("Usage: python3.12 rerun_extraction.py <workflow_output_folder> <gene_symbol>")
        print("\nExample:")
        print("  python3.12 rerun_extraction.py automated_output/TTR/20251125_114028 TTR")
        sys.exit(1)

    workflow_folder = sys.argv[1]
    gene_symbol = sys.argv[2]

    result = rerun_extraction_on_folder(
        workflow_folder=workflow_folder,
        gene_symbol=gene_symbol,
        tier_threshold=1,
        max_workers=5
    )

    if not result["success"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
