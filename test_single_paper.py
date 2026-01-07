#!/usr/bin/env python3
"""
Test script to extract data from a single paper by PMID.

Usage:
    python test_single_paper.py --pmid 16116052 --gene KCNH2 --email your@email.com
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure logging
from utils.logging_utils import setup_logging, get_logger
setup_logging(level=logging.INFO)
logger = get_logger(__name__)


def test_single_paper(pmid: str, gene_symbol: str, email: str, output_dir: str = None):
    """
    Test extraction workflow for a single PMID.

    Args:
        pmid: PubMed ID to test
        gene_symbol: Gene symbol (e.g., "KCNH2", "SCN5A")
        email: Your email for NCBI E-utilities
        output_dir: Optional output directory (defaults to ./test_output_{pmid})
    """
    from harvesting import PMCHarvester
    from harvesting.abstracts import fetch_and_save_abstracts
    from pipeline.extraction import ExpertExtractor
    from utils.models import Paper

    # Setup output directory
    if output_dir is None:
        output_dir = f"./test_output_{pmid}"

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("="*80)
    logger.info(f"TESTING SINGLE PAPER EXTRACTION")
    logger.info("="*80)
    logger.info(f"PMID: {pmid}")
    logger.info(f"Gene: {gene_symbol}")
    logger.info(f"Output: {output_path}")
    logger.info("="*80)

    # ============================================================================
    # STEP 1: Fetch Abstract and Metadata
    # ============================================================================
    logger.info("\n📝 STEP 1: Fetching abstract and metadata...")

    abstract_dir = output_path / "abstract_json"
    abstract_records = fetch_and_save_abstracts(
        pmids=[pmid],
        output_dir=str(abstract_dir),
        email=email,
    )

    if pmid not in abstract_records:
        logger.error(f"Failed to fetch abstract for PMID {pmid}")
        return {"success": False, "error": "Abstract fetch failed"}

    logger.info(f"✓ Saved abstract to: {abstract_records[pmid]}")

    # Load abstract for display
    with open(abstract_records[pmid], 'r', encoding='utf-8') as f:
        abstract_data = json.load(f)

    logger.info(f"\nPaper Title: {abstract_data['metadata']['title']}")
    logger.info(f"Authors: {', '.join(abstract_data['metadata']['authors'][:3])} et al.")
    logger.info(f"Journal: {abstract_data['metadata']['journal']}")
    logger.info(f"Year: {abstract_data['metadata']['year']}")

    # ============================================================================
    # STEP 2: Download Full-Text from PMC
    # ============================================================================
    logger.info("\n📥 STEP 2: Downloading full-text from PubMed Central...")

    harvest_dir = output_path / "pmc_fulltext"
    harvester = PMCHarvester(output_dir=str(harvest_dir), gene_symbol=gene_symbol)

    harvester.harvest([pmid], delay=1.0)

    # Check if download was successful
    success_log = harvest_dir / "successful_downloads.csv"
    if not success_log.exists():
        logger.warning("⚠️  Paper may not be available in PMC or download failed")
        logger.info("Checking for downloaded files...")

    # Look for the downloaded markdown file
    data_zones_file = harvest_dir / f"{pmid}_DATA_ZONES.md"
    full_context_file = harvest_dir / f"{pmid}_FULL_CONTEXT.md"

    markdown_file = None
    if data_zones_file.exists():
        markdown_file = data_zones_file
        logger.info(f"✓ Found DATA_ZONES file: {markdown_file}")
    elif full_context_file.exists():
        markdown_file = full_context_file
        logger.info(f"✓ Found FULL_CONTEXT file: {markdown_file}")
    else:
        logger.error("❌ No markdown file found. Paper may not be in PMC.")
        logger.info("\nPossible reasons:")
        logger.info("1. Paper is not open access in PubMed Central")
        logger.info("2. PMID is invalid or not in PMC database")
        logger.info("3. Network/API error during download")
        logger.info("\nYou can manually download the paper and place it as:")
        logger.info(f"   {full_context_file}")
        return {"success": False, "error": "Full-text download failed"}

    # ============================================================================
    # STEP 3: Extract Variant Data
    # ============================================================================
    logger.info("\n🧬 STEP 3: Extracting variant and patient data...")

    extraction_dir = output_path / "extractions"
    extraction_dir.mkdir(exist_ok=True)

    # Read the markdown content
    with open(markdown_file, 'r', encoding='utf-8') as f:
        full_text = f.read()

    logger.info(f"Full-text length: {len(full_text):,} characters")

    # Create paper object
    paper = Paper(
        pmid=pmid,
        title=abstract_data['metadata']['title'],
        abstract=abstract_data['abstract'],
        full_text=full_text,
        gene_symbol=gene_symbol,
        authors=abstract_data['metadata']['authors'],
        journal=abstract_data['metadata']['journal'],
        publication_date=abstract_data['metadata']['year'],
    )

    # Extract data
    extractor = ExpertExtractor(tier_threshold=1, fulltext_dir=str(harvest_dir))

    try:
        extraction_result = extractor.extract(paper)

        if extraction_result.success:
            # Save extraction
            output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid}.json"
            with open(output_file, 'w') as f:
                json.dump(extraction_result.extracted_data, f, indent=2)

            # Display summary
            num_variants = extraction_result.extracted_data.get(
                'extraction_metadata', {}
            ).get('total_variants_found', 0)

            logger.info(f"\n✓ EXTRACTION SUCCESSFUL!")
            logger.info(f"✓ Found {num_variants} variant(s)")
            logger.info(f"✓ Saved to: {output_file}")

            # Display sample of extracted data
            logger.info("\n" + "="*80)
            logger.info("EXTRACTION SUMMARY")
            logger.info("="*80)

            variants = extraction_result.extracted_data.get('variants', [])
            logger.info(f"\nTotal variants: {len(variants)}")

            if variants:
                logger.info("\nFirst few variants:")
                for i, variant in enumerate(variants[:3], 1):
                    logger.info(f"\n{i}. Variant: {variant.get('variant_name', 'Unknown')}")
                    logger.info(f"   Protein change: {variant.get('protein_change', 'N/A')}")
                    logger.info(f"   DNA change: {variant.get('dna_change', 'N/A')}")

                    carriers = variant.get('carriers', [])
                    logger.info(f"   Carriers: {len(carriers)}")
                    if carriers:
                        affected = sum(1 for c in carriers if c.get('clinical_status') == 'Affected')
                        unaffected = sum(1 for c in carriers if c.get('clinical_status') == 'Unaffected')
                        logger.info(f"      Affected: {affected}, Unaffected: {unaffected}")

            logger.info("\n" + "="*80)
            logger.info(f"Full extraction saved to: {output_file}")
            logger.info("="*80)

            return {
                "success": True,
                "pmid": pmid,
                "variants_found": num_variants,
                "output_file": str(output_file),
                "extraction_data": extraction_result.extracted_data
            }

        else:
            logger.error(f"❌ Extraction failed: {extraction_result.error}")
            return {
                "success": False,
                "error": extraction_result.error
            }

    except Exception as e:
        logger.error(f"❌ Extraction error: {e}", exc_info=True)
        return {
            "success": False,
            "error": str(e)
        }


def main():
    parser = argparse.ArgumentParser(
        description="Test extraction for a single paper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test KCNH2 paper (PMID 16116052)
  python test_single_paper.py --pmid 16116052 --gene KCNH2 --email your@email.com

  # Test with custom output directory
  python test_single_paper.py --pmid 16116052 --gene KCNH2 --email your@email.com --output ./my_test
        """
    )

    parser.add_argument("--pmid", required=True, help="PubMed ID to test")
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g., KCNH2, SCN5A)")
    parser.add_argument("--email", "-e", required=True, help="Your email for NCBI E-utilities")
    parser.add_argument("--output", "-o", help="Output directory (default: ./test_output_{PMID})")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Check for API keys
    if not (os.getenv("AI_INTEGRATIONS_OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY")):
        logger.error("⚠️  ERROR: OpenAI API key not found!")
        logger.error("Please set AI_INTEGRATIONS_OPENAI_API_KEY or OPENAI_API_KEY")
        sys.exit(1)

    # Run test
    try:
        result = test_single_paper(
            pmid=args.pmid,
            gene_symbol=args.gene,
            email=args.email,
            output_dir=args.output,
        )

        if result["success"]:
            logger.info("\n✅ Test completed successfully!")
            sys.exit(0)
        else:
            logger.error(f"\n❌ Test failed: {result.get('error', 'Unknown error')}")
            sys.exit(1)

    except KeyboardInterrupt:
        logger.warning("\n⚠️  Test interrupted by user")
        sys.exit(1)

    except Exception as e:
        logger.error(f"\n❌ Test failed with error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
