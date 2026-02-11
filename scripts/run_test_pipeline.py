#!/usr/bin/env python3
"""Run the full pipeline on test PMIDs."""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from dotenv import load_dotenv

load_dotenv()

import logging

from gui.worker import PipelineWorker, create_job

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def parse_pmids(filepath: Path) -> list[str]:
    """Parse PMIDs from test file."""
    pmids = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                # Extract just the PMID (first field)
                pmid = line.split()[0]
                if pmid.isdigit():
                    pmids.append(pmid)
    return pmids


def main():
    # Parse PMIDs
    pmid_file = Path(__file__).parent / "test_pmids.txt"
    pmids = parse_pmids(pmid_file)
    logger.info(f"Loaded {len(pmids)} test PMIDs")

    # Create job with specific PMIDs (skips discovery)
    job = create_job(
        gene_symbol="KCNH2",
        email="brett.kroncke@gmail.com",
        output_dir="tests/test_run_output",
        max_pmids=100,  # Not used when specific_pmids is set
        max_papers_to_download=100,  # Download all
        tier_threshold=1,
        specific_pmids=pmids,  # Skip discovery, use these PMIDs
        enable_tier1=False,  # Skip keyword filtering for test
        enable_tier2=False,  # Skip LLM filtering for test
        scout_enabled=True,  # Run data scout
    )

    logger.info(f"Created job: {job.job_id}")
    logger.info(f"Output dir: {job.output_path}")

    # Run the pipeline
    worker = PipelineWorker()
    try:
        result = worker.run(job)
        logger.info("=" * 80)
        logger.info("PIPELINE COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Discovered PMIDs: {len(result.discovered_pmids)}")
        logger.info(f"Downloaded PMIDs: {len(result.downloaded_pmids)}")
        logger.info(f"Extracted PMIDs: {len(result.extracted_pmids)}")
        logger.info(f"Output directory: {result.output_path}")
    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
