#!/usr/bin/env python3
"""Run extraction on existing folder with downloaded papers."""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from dotenv import load_dotenv

load_dotenv()

import logging

from gui.worker import PipelineWorker, create_folder_job

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    # Use the existing folder with downloaded papers
    folder_path = "tests/test_run_output/KCNH2/20260125_113957"

    # Resolve to absolute path
    folder_path = str(Path(folder_path).resolve())

    logger.info(f"Processing existing folder: {folder_path}")

    # Check what files exist
    pmc_dir = Path(folder_path) / "pmc_fulltext"
    full_context_files = list(pmc_dir.glob("*_FULL_CONTEXT.md"))
    data_zones_files = list(pmc_dir.glob("*_DATA_ZONES.md"))

    logger.info(f"Found {len(full_context_files)} FULL_CONTEXT files")
    logger.info(f"Found {len(data_zones_files)} DATA_ZONES files")

    # Get PMIDs from files
    full_context_pmids = [
        f.name.replace("_FULL_CONTEXT.md", "") for f in full_context_files
    ]
    data_zones_pmids = [f.name.replace("_DATA_ZONES.md", "") for f in data_zones_files]

    # Check for existing extractions
    extraction_dir = Path(folder_path) / "extractions"
    existing_extractions = (
        list(extraction_dir.glob("*.json")) if extraction_dir.exists() else []
    )
    extraction_pmids = [
        f.name.split("_PMID_")[1].replace(".json", "")
        for f in existing_extractions
        if "_PMID_" in f.name
    ]

    logger.info(f"Found {len(extraction_pmids)} existing extractions")

    # Create folder job (skips discovery/download, goes to scout/extract)
    job = create_folder_job(
        folder_path=folder_path,
        gene_symbol="KCNH2",
        email="brett.kroncke@gmail.com",
        run_scout=False,  # DATA_ZONES already exist from download step
        skip_already_extracted=True,
        tier_threshold=1,
        scout_enabled=True,  # Prefer DATA_ZONES for extraction
        full_context_pmids=full_context_pmids,
        data_zones_pmids=data_zones_pmids,
        extraction_pmids=extraction_pmids,
    )

    logger.info(f"Created folder job: {job.job_id}")
    logger.info(f"Starting step: {job.current_step}")

    # Run the pipeline (will start at extraction step)
    worker = PipelineWorker()
    try:
        result = worker.run(job)
        logger.info("=" * 80)
        logger.info("PIPELINE COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Extracted PMIDs: {len(result.extracted_pmids)}")
        logger.info(f"Output directory: {result.output_path}")

        # Show extraction summary
        extraction_dir = Path(folder_path) / "extractions"
        if extraction_dir.exists():
            extraction_files = list(extraction_dir.glob("*.json"))
            logger.info(f"Extraction files created: {len(extraction_files)}")

    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
