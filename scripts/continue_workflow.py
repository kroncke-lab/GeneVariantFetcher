#!/usr/bin/env python3
"""
Continue an interrupted GVF workflow from filtering step.

This script:
1. Loads already-fetched abstracts
2. Extracts already-filtered PMIDs from the log
3. Continues filtering remaining PMIDs
4. Downloads and extracts papers that passed
"""

import json
import logging
import os
import re
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv
load_dotenv()

from config.settings import get_settings
from pipeline.steps import (
    filter_papers,
    download_fulltext,
    extract_variants,
    aggregate_data,
)
from utils.logging_utils import setup_logging, get_logger

setup_logging(level=logging.INFO)
logger = get_logger(__name__)


def extract_filtered_pmids_from_log(log_path: Path) -> tuple[set, set]:
    """Extract passed and failed PMIDs from workflow log."""
    passed = set()
    failed = set()
    
    with open(log_path) as f:
        for line in f:
            if "Intern filter: pass" in line:
                m = re.search(r"PMID (\d+)", line)
                if m:
                    passed.add(m.group(1))
            elif "Intern filter: fail" in line:
                m = re.search(r"PMID (\d+)", line)
                if m:
                    failed.add(m.group(1))
    
    return passed, failed


def load_abstracts(abstract_dir: Path) -> dict:
    """Load abstract record paths (pmid -> file path mapping)."""
    records = {}
    for f in abstract_dir.glob("*.json"):
        try:
            with open(f) as fp:
                data = json.load(fp)
                pmid = data.get("pmid") or f.stem.split("_")[-1]
                # filter_papers expects pmid -> path, not pmid -> data
                records[str(pmid)] = str(f)
        except Exception as e:
            logger.warning(f"Failed to load {f}: {e}")
    return records


def continue_workflow(
    run_dir: str,
    gene_symbol: str,
    max_downloads: int = 500,
    skip_remaining_filter: bool = False,
):
    """Continue interrupted workflow from run directory."""
    run_path = Path(run_dir)
    settings = get_settings()
    
    # Find key files
    log_path = run_path / f"{gene_symbol}_workflow.log"
    pmids_path = run_path / f"{gene_symbol}_pmids.txt"
    abstract_dir = run_path / "abstract_json"
    
    if not log_path.exists():
        logger.error(f"Log file not found: {log_path}")
        return
    
    # Add file handler to log
    file_handler = logging.FileHandler(log_path, mode='a')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)
    
    logger.info("=" * 80)
    logger.info(f"CONTINUING WORKFLOW FOR {gene_symbol}")
    logger.info("=" * 80)
    
    # Load all PMIDs
    with open(pmids_path) as f:
        all_pmids = [line.strip() for line in f if line.strip()]
    logger.info(f"Total PMIDs in run: {len(all_pmids)}")
    
    # Extract already-filtered PMIDs
    already_passed, already_failed = extract_filtered_pmids_from_log(log_path)
    already_filtered = already_passed | already_failed
    logger.info(f"Already filtered: {len(already_filtered)} ({len(already_passed)} passed, {len(already_failed)} failed)")
    
    # Find PMIDs still needing filtering
    remaining_pmids = [p for p in all_pmids if p not in already_filtered]
    logger.info(f"Remaining to filter: {len(remaining_pmids)}")
    
    # Load abstracts
    logger.info("Loading abstract records...")
    abstract_records = load_abstracts(abstract_dir)
    logger.info(f"Loaded {len(abstract_records)} abstracts")
    
    # Filter remaining PMIDs (unless skipping)
    if remaining_pmids and not skip_remaining_filter:
        logger.info(f"\n🧹 Filtering {len(remaining_pmids)} remaining PMIDs...")
        
        filter_result = filter_papers(
            pmids=remaining_pmids,
            abstract_records=abstract_records,
            gene_symbol=gene_symbol,
            output_path=run_path,
            enable_tier1=settings.enable_tier1,
            enable_tier2=settings.enable_tier2,
            use_clinical_triage=False,
            tier1_min_keywords=settings.tier1_min_keywords,
            tier2_confidence_threshold=settings.tier2_confidence_threshold,
        )
        
        new_passed = set(filter_result.data.get("filtered_pmids", []))
        new_failed = set(filter_result.data.get("dropped_pmids", []))
        
        already_passed.update(new_passed)
        already_failed.update(new_failed)
        
        logger.info(f"Filtering complete. New totals: {len(already_passed)} passed, {len(already_failed)} failed")
    
    # Now we have all passed PMIDs
    filtered_pmids = list(already_passed)
    logger.info(f"\n📥 Downloading {min(len(filtered_pmids), max_downloads)} papers...")
    
    # Download full-text papers
    download_result = download_fulltext(
        pmids=filtered_pmids,
        output_path=run_path,
        gene_symbol=gene_symbol,
        max_papers=max_downloads,
        delay=2.0,
    )
    
    harvest_dir = download_result.data.get("harvest_dir")
    downloaded_pmids = download_result.data.get("downloaded_pmids", [])
    abstract_only_pmids = download_result.data.get("abstract_only_pmids", [])
    
    logger.info(f"✓ Downloaded {len(downloaded_pmids)} full-text papers")
    logger.info(f"✓ {len(abstract_only_pmids)} papers will use abstract-only")
    
    # Extract variants
    logger.info("\n🧬 Extracting variants...")
    extraction_dir = run_path / "extractions"
    
    extract_result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol=gene_symbol,
        tier_threshold=1,
    )
    
    extracted_count = extract_result.stats.get("papers_extracted", 0)
    variants_found = extract_result.stats.get("variants_extracted", 0)
    
    logger.info(f"✓ Extracted from {extracted_count} papers")
    logger.info(f"✓ Found {variants_found} variants")
    
    # Aggregate
    logger.info("\n📊 Aggregating results...")
    aggregate_result = aggregate_data(
        extraction_dir=extraction_dir,
        output_path=run_path,
        gene_symbol=gene_symbol,
    )
    
    logger.info(f"✓ Aggregation complete: {aggregate_result.stats}")
    
    logger.info("\n" + "=" * 80)
    logger.info("WORKFLOW CONTINUATION COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Continue interrupted GVF workflow")
    parser.add_argument("run_dir", help="Path to the run directory (e.g., /path/to/KCNH2/20260202_173749)")
    parser.add_argument("--gene", "-g", required=True, help="Gene symbol (e.g., KCNH2)")
    parser.add_argument("--max-downloads", "-m", type=int, default=500, help="Max papers to download")
    parser.add_argument("--skip-filter", action="store_true", help="Skip filtering remaining PMIDs")
    
    args = parser.parse_args()
    
    continue_workflow(
        run_dir=args.run_dir,
        gene_symbol=args.gene,
        max_downloads=args.max_downloads,
        skip_remaining_filter=args.skip_filter,
    )
