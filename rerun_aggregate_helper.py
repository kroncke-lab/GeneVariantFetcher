#!/usr/bin/env python3
"""
Helper script for aggregating penetrance data from rerun extraction outputs.

Usage:
    python rerun_aggregate_helper.py \
        --extraction-dir /path/to/workflow/extractions_rerun/20251126_161413 \
        --gene KCNH2

The script writes a summary JSON next to the workflow folder unless you override
the output file via --output-file.
"""

from argparse import ArgumentParser
from datetime import datetime
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = ArgumentParser(description="Aggregate penetrance from rerun extractions")
    parser.add_argument(
        "--extraction-dir",
        "-e",
        required=True,
        help="Path to the rerun extraction directory containing *_PMID_*.json files",
    )
    parser.add_argument(
        "--gene",
        "-g",
        required=True,
        help="Gene symbol used when the extractions were run (e.g., KCNH2)",
    )
    parser.add_argument(
        "--output-file",
        "-o",
        help="Optional path to write the penetrance summary JSON",
    )
    return parser.parse_args()


def main():
    from pipeline.aggregation import aggregate_penetrance

    args = parse_args()
    extraction_dir = Path(args.extraction_dir)
    if not extraction_dir.exists():
        logger.error("Extraction directory not found: %s", extraction_dir)
        raise SystemExit(1)

    if args.output_file:
        summary_file = Path(args.output_file)
    else:
        parent = extraction_dir.parent
        workflow_dir = parent.parent if parent.name == "extractions_rerun" else parent
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        summary_file = workflow_dir / f"{args.gene}_penetrance_summary_rerun_{timestamp}.json"

    summary_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Aggregating penetrance data from %s", extraction_dir)
    result = aggregate_penetrance(
        extraction_dir=extraction_dir,
        gene_symbol=args.gene,
        output_file=summary_file,
    )

    logger.info("Saved penetrance summary (%d variants) to %s", result.get("total_variants", 0), summary_file)


if __name__ == "__main__":
    main()

