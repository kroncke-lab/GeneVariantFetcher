#!/usr/bin/env python3
"""
Run KCNH2 full download pipeline WITHOUT extraction.
Discovery → Synonyms → PMIDs → Abstracts → Filter → Download → Scout only.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

# Add project root
sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv
load_dotenv()

from gui.worker import PipelineWorker, create_job, LoggingProgressCallback
from gui.checkpoint import PipelineStep
from utils.logging_utils import setup_logging

setup_logging(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    output_dir = "/mnt/temp2/kronckbm/gvf_output"

    job = create_job(
        gene_symbol="KCNH2",
        email="brett.kroncke@gmail.com",
        output_dir=output_dir,
        max_pmids=10000,
        max_papers_to_download=1000,
        auto_synonyms=True,
        enable_tier1=True,
        enable_tier2=True,
        tier2_confidence_threshold=0.5,
        scout_enabled=True,
        use_pubmed=True,
        use_europepmc=False,
    )

    logger.info(f"Job created: {job.job_id}")
    logger.info(f"Gene: {job.gene_symbol}")
    logger.info(f"Auto synonyms: {job.auto_synonyms}")
    logger.info(f"Output: {job.output_dir}")

    worker = PipelineWorker()
    callback = LoggingProgressCallback()
    worker.progress_callback = callback

    # Monkey-patch: skip extraction, aggregation, migration steps
    original_run = worker.run

    def run_download_only(checkpoint, resume=False):
        """Run pipeline but stop before extraction."""
        # We'll run the full pipeline but override _step_extract_variants to be a no-op
        return original_run(checkpoint, resume=resume)

    # Make extraction a no-op
    def _noop_extract(checkpoint):
        logger.info("=== SKIPPING EXTRACTION (download-only mode) ===")
        checkpoint.step_stats[PipelineStep.EXTRACTING_VARIANTS.value] = {
            "skipped": True,
            "reason": "download-only test"
        }
        return checkpoint

    def _noop_aggregate(checkpoint):
        logger.info("=== SKIPPING AGGREGATION (download-only mode) ===")
        return checkpoint

    def _noop_migrate(checkpoint):
        logger.info("=== SKIPPING MIGRATION (download-only mode) ===")
        return checkpoint

    worker._step_extract_variants = _noop_extract
    worker._step_aggregate_data = _noop_aggregate
    worker._step_migrate_database = _noop_migrate

    logger.info("Starting download-only pipeline...")
    start = datetime.now()

    try:
        result = worker.run(job)
        elapsed = (datetime.now() - start).total_seconds()
        logger.info(f"\n{'='*60}")
        logger.info(f"PIPELINE COMPLETE in {elapsed:.1f}s")
        logger.info(f"Job ID: {result.job_id}")
        logger.info(f"Final step: {result.current_step.value}")
        logger.info(f"Output: {result.output_path}")

        # Print stats
        for step_name, stats in result.step_stats.items():
            logger.info(f"\n--- {step_name} ---")
            for k, v in stats.items():
                logger.info(f"  {k}: {v}")

    except Exception as e:
        elapsed = (datetime.now() - start).total_seconds()
        logger.error(f"Pipeline failed after {elapsed:.1f}s: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
