#!/usr/bin/env python3
"""
Automated Workflow Example: From Gene to Variant Data

This script demonstrates the complete automated workflow:
1. Fetch relevant PMIDs from PubMind for a gene
2. Download full-text articles from PubMed Central
3. Extract individual-level variant and patient data
4. Save structured results to JSON

This is the "hands-free" approach - just provide a gene symbol and get
comprehensive variant/patient data from all publicly available papers.
"""

import os
import sys
import json
import logging
from pathlib import Path
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def automated_variant_extraction_workflow(
    gene_symbol: str,
    email: str,
    max_pmids: int = 100,
    max_papers_to_download: int = 50,
    output_dir: str = "automated_output"
):
    """
    Complete automated workflow from gene symbol to extracted variant data.

    Args:
        gene_symbol: Gene to search for (e.g., "BRCA1", "SCN5A")
        email: Your email for NCBI E-utilities (required)
        max_pmids: Maximum PMIDs to fetch from PubMind/PubMed
        max_papers_to_download: Maximum papers to download full-text
        output_dir: Directory to save all outputs

    Returns:
        Dictionary with workflow results and statistics
    """
    from pubmind_fetcher import fetch_pmids_for_gene
    from harvest_pmc_fulltext import PMCHarvester

    output_path = Path(output_dir) / gene_symbol / datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("="*80)
    logger.info(f"AUTOMATED WORKFLOW FOR GENE: {gene_symbol}")
    logger.info("="*80)

    # ============================================================================
    # STEP 1: Fetch PMIDs from PubMind
    # ============================================================================
    logger.info("\n📚 STEP 1: Discovering relevant papers from PubMind...")

    pmids_file = output_path / f"{gene_symbol}_pmids.txt"
    pmids = fetch_pmids_for_gene(
        gene_symbol=gene_symbol,
        email=email,
        max_results=max_pmids,
        output_file=pmids_file
    )

    logger.info(f"✓ Found {len(pmids)} PMIDs for {gene_symbol}")
    logger.info(f"✓ Saved PMID list to: {pmids_file}")

    if not pmids:
        logger.warning(f"No PMIDs found for {gene_symbol}. Workflow terminated.")
        return {"success": False, "error": "No PMIDs found"}

    # ============================================================================
    # STEP 2: Download Full-Text Papers from PMC
    # ============================================================================
    logger.info("\n📥 STEP 2: Downloading full-text papers from PubMed Central...")

    harvest_dir = output_path / "pmc_fulltext"
    harvester = PMCHarvester(output_dir=str(harvest_dir))

    # Limit downloads to avoid excessive processing time
    pmids_to_download = pmids[:max_papers_to_download]
    logger.info(f"Downloading up to {len(pmids_to_download)} papers (out of {len(pmids)} total PMIDs)")

    harvester.harvest(pmids_to_download, delay=2.0)

    # Check how many were successfully downloaded
    success_log = harvest_dir / "successful_downloads.csv"
    if success_log.exists():
        import pandas as pd
        successful_downloads = pd.read_csv(success_log)
        num_downloaded = len(successful_downloads)
        logger.info(f"✓ Successfully downloaded {num_downloaded} full-text papers")
    else:
        num_downloaded = 0
        logger.warning("No papers were successfully downloaded from PMC")

    # ============================================================================
    # STEP 3: Extract Variant and Patient Data
    # ============================================================================
    logger.info("\n🧬 STEP 3: Extracting variant and patient data using AI...")

    extraction_dir = output_path / "extractions"
    extraction_dir.mkdir(exist_ok=True)

    # Get list of downloaded markdown files
    markdown_files = list(harvest_dir.glob("*_FULL_CONTEXT.md"))
    logger.info(f"Found {len(markdown_files)} markdown files to process")

    # Process each paper
    from models import Paper
    from extractor import ExpertExtractor

    extractor = ExpertExtractor(model="gpt-4o")
    extractions = []

    for md_file in markdown_files:
        # Extract PMID from filename (format: PMID_12345678_FULL_CONTEXT.md)
        pmid_match = md_file.stem.split('_')[1] if '_' in md_file.stem else None

        if not pmid_match:
            logger.warning(f"Could not extract PMID from filename: {md_file.name}")
            continue

        logger.info(f"\nProcessing PMID {pmid_match}...")

        # Read the markdown content
        with open(md_file, 'r', encoding='utf-8') as f:
            full_text = f.read()

        # Create paper object
        paper = Paper(
            pmid=pmid_match,
            title=f"Paper {pmid_match}",  # Title would be in the markdown
            full_text=full_text,
            gene_symbol=gene_symbol
        )

        # Extract data
        try:
            extraction_result = extractor.extract(paper)

            if extraction_result.success:
                extractions.append(extraction_result)

                # Save individual extraction
                output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid_match}.json"
                with open(output_file, 'w') as f:
                    json.dump(extraction_result.extracted_data, f, indent=2)

                num_variants = extraction_result.extracted_data.get(
                    'extraction_metadata', {}
                ).get('total_variants_found', 0)

                logger.info(f"✓ Extracted {num_variants} variants from PMID {pmid_match}")
                logger.info(f"✓ Saved to: {output_file}")
            else:
                logger.warning(f"✗ Extraction failed for PMID {pmid_match}: {extraction_result.error}")

        except Exception as e:
            logger.error(f"Error processing PMID {pmid_match}: {e}")
            continue

    # ============================================================================
    # STEP 4: Compile Results and Statistics
    # ============================================================================
    logger.info("\n📊 STEP 4: Compiling results and statistics...")

    total_variants = 0
    for extraction in extractions:
        if extraction.extracted_data:
            total_variants += extraction.extracted_data.get(
                'extraction_metadata', {}
            ).get('total_variants_found', 0)

    # Create summary report
    summary = {
        "gene_symbol": gene_symbol,
        "workflow_timestamp": datetime.now().isoformat(),
        "statistics": {
            "pmids_discovered": len(pmids),
            "papers_downloaded": num_downloaded,
            "papers_extracted": len(extractions),
            "total_variants_found": total_variants,
            "success_rate": f"{len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%"
        },
        "output_locations": {
            "pmid_list": str(pmids_file),
            "full_text_papers": str(harvest_dir),
            "extractions": str(extraction_dir)
        }
    }

    summary_file = output_path / f"{gene_symbol}_workflow_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info("\n" + "="*80)
    logger.info("WORKFLOW COMPLETE!")
    logger.info("="*80)
    logger.info(f"Gene: {gene_symbol}")
    logger.info(f"PMIDs discovered: {len(pmids)}")
    logger.info(f"Papers downloaded: {num_downloaded}")
    logger.info(f"Papers with extractions: {len(extractions)}")
    logger.info(f"Total variants found: {total_variants}")
    logger.info(f"Success rate: {len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%")
    logger.info(f"\nAll outputs saved to: {output_path}")
    logger.info(f"Summary report: {summary_file}")
    logger.info("="*80)

    return summary


def main():
    """
    Main entry point for automated workflow.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Automated variant extraction workflow from gene symbol to structured data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract all variant data for BRCA1
  python example_automated_workflow.py BRCA1 --email your@email.com

  # Extract data for SCN5A with custom limits
  python example_automated_workflow.py SCN5A --email your@email.com --max-pmids 200 --max-downloads 100

  # Quick test with small dataset
  python example_automated_workflow.py TP53 --email your@email.com --max-pmids 10 --max-downloads 5
        """
    )

    parser.add_argument("gene", help="Gene symbol (e.g., BRCA1, SCN5A, TP53)")
    parser.add_argument("--email", "-e", required=True, help="Your email for NCBI E-utilities")
    parser.add_argument("--max-pmids", type=int, default=100,
                       help="Maximum PMIDs to fetch (default: 100)")
    parser.add_argument("--max-downloads", type=int, default=50,
                       help="Maximum papers to download (default: 50)")
    parser.add_argument("--output", "-o", default="automated_output",
                       help="Output directory (default: automated_output)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Check for API keys
    if not (os.getenv("AI_INTEGRATIONS_OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY")):
        logger.error("⚠️  ERROR: OpenAI API key not found!")
        logger.error("Please set AI_INTEGRATIONS_OPENAI_API_KEY or OPENAI_API_KEY")
        logger.error("Example: export AI_INTEGRATIONS_OPENAI_API_KEY='your-key-here'")
        sys.exit(1)

    # Run workflow
    try:
        automated_variant_extraction_workflow(
            gene_symbol=args.gene,
            email=args.email,
            max_pmids=args.max_pmids,
            max_papers_to_download=args.max_downloads,
            output_dir=args.output
        )

        # Exit with success code
        sys.exit(0)

    except KeyboardInterrupt:
        logger.warning("\n⚠️  Workflow interrupted by user")
        sys.exit(1)

    except Exception as e:
        logger.error(f"\n❌ Workflow failed with error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
