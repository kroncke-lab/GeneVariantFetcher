#!/usr/bin/env python3
"""
Minimal test to verify GVF extraction pipeline works on a single PMID.
Uses an existing FULL_CONTEXT.md file to test extraction.
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Ensure project root is in path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Load environment
from dotenv import load_dotenv

load_dotenv()


def test_extraction():
    """Test extraction on a single PMID."""

    # Use existing FULL_CONTEXT.md file
    test_file = Path("gvf_output/browser_downloads/PMID_16414944_FULL_CONTEXT.md")

    if not test_file.exists():
        logger.error(f"Test file not found: {test_file}")
        return False, "Test file not found"

    logger.info(f"Using test file: {test_file}")

    # Read the content
    content = test_file.read_text()
    content_length = len(content)
    logger.info(f"Content length: {content_length} characters")

    # Import extraction module
    from pipeline.extraction import ExpertExtractor
    from utils.models import Paper

    # Create a Paper object (gene_symbol is passed in the Paper object)
    paper = Paper(
        pmid="16414944",
        title="Genetic Testing in the Long QT Syndrome",
        abstract="In long QT syndrome (LQTS), disease severity and response to therapy vary...",
        full_text=content,
        gene_symbol="KCNH2",  # Gene symbol is part of Paper
    )

    # Initialize extractor
    logger.info("Initializing ExpertExtractor...")
    extractor = ExpertExtractor(tier_threshold=1)

    # Run extraction
    logger.info("Running extraction (this may take a moment)...")
    start_time = datetime.now()

    try:
        result = extractor.extract(paper)
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        logger.info(f"Extraction completed in {duration:.1f} seconds")

        if result.success:
            variant_count = len(result.variants)
            logger.info(f"SUCCESS: Extracted {variant_count} variants")

            # Show first few variants
            if result.variants:
                logger.info("Sample variants:")
                for v in result.variants[:5]:
                    logger.info(f"  - {v.variant_notation or v.variant_id}")

            return True, {
                "pmid": "16414944",
                "variants_found": variant_count,
                "duration_seconds": duration,
                "model_used": result.model_used
                if hasattr(result, "model_used")
                else "unknown",
            }
        else:
            logger.error(f"Extraction failed: {result.error}")
            return False, f"Extraction error: {result.error}"

    except Exception as e:
        logger.error(f"Exception during extraction: {e}")
        import traceback

        traceback.print_exc()
        return False, str(e)


def main():
    logger.info("=" * 60)
    logger.info("GVF Single PMID Test")
    logger.info("=" * 60)

    success, result = test_extraction()

    # Write results
    results_path = Path("TEST_RESULTS.md")

    with open(results_path, "w") as f:
        f.write("# GVF Pipeline Test Results\n\n")
        f.write(f"**Test Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        if success:
            f.write("## ✅ Test PASSED\n\n")
            f.write(f"- **PMID:** {result['pmid']}\n")
            f.write(f"- **Variants Found:** {result['variants_found']}\n")
            f.write(f"- **Duration:** {result['duration_seconds']:.1f} seconds\n")
            f.write(f"- **Model Used:** {result.get('model_used', 'unknown')}\n")
        else:
            f.write("## ❌ Test FAILED\n\n")
            f.write(f"- **Error:** {result}\n")

        f.write("\n## Test Details\n\n")
        f.write(
            "- **Test File:** `gvf_output/browser_downloads/PMID_16414944_FULL_CONTEXT.md`\n"
        )
        f.write("- **Gene Symbol:** KCNH2\n")
        f.write("- **Test Type:** Full-text extraction using ExpertExtractor\n")

    logger.info(f"Results written to: {results_path}")

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
