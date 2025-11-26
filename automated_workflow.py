#!/usr/bin/env python3
"""
Automated Workflow Entrypoint: From Gene to Variant Data

This script provides the production-ready, end-to-end automated workflow:
1. Fetch relevant PMIDs from PubMind for a gene
2. Download full-text articles from PubMed Central
3. Extract individual-level variant and patient data
4. Save structured results to JSON and aggregate penetrance metrics

Run this script directly for the canonical workflow. The example script
(`example_automated_workflow.py`) is now a thin wrapper around this module
for demonstration purposes.
"""

import os
import sys
import json
import logging
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv
from concurrent.futures import ThreadPoolExecutor, as_completed

# Load environment variables from .env file
load_dotenv()

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
    output_dir: str = "automated_output",
    tier_threshold: int = 1,
    skip_extraction: bool = False,
    model: str = None,
    max_workers: int = None,
    max_workers_cap: int = 8,
    rate_limit_delay: float = 0.0,
    migrate_to_sqlite: bool = False,
    db_path: str = None,
    cleanup_after_migrate: bool = False,
):
    """
    Complete automated workflow from gene symbol to extracted variant data.

    Args:
        gene_symbol: Gene to search for (e.g., "BRCA1", "SCN5A")
        email: Your email for NCBI E-utilities (required)
        max_pmids: Maximum PMIDs to fetch from PubMind/PubMed
        max_papers_to_download: Maximum papers to download full-text
        output_dir: Directory to save all outputs
        tier_threshold: If the first model finds fewer variants than this, the next model is tried.
        skip_extraction: If True, skip LLM extraction (harvest only)
        model: LLM model to use (e.g., 'gpt-4o', 'claude-3-opus'). If None, uses config default.
        max_workers: Number of parallel workers. If None, auto-detects based on CPU count.
        max_workers_cap: Maximum cap for auto-detected workers (default: 8)
        rate_limit_delay: Delay in seconds between API requests to avoid rate limits
        migrate_to_sqlite: If True, migrate extraction results to SQLite database
        db_path: SQLite database path. If None, auto-generates as {gene_symbol}.db
        cleanup_after_migrate: If True, archive PMC files after successful migration
    """
    from pubmind_fetcher import fetch_pmids_for_gene
    from harvesting import PMCHarvester

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
    # STEP 3: Extract Variant and Patient Data (Optional)
    # ============================================================================
    if not skip_extraction:
        logger.info("\n🧬 STEP 3: Extracting variant and patient data using AI...")

        extraction_dir = output_path / "extractions"
        extraction_dir.mkdir(exist_ok=True)

        # Get list of downloaded markdown files
        markdown_files = list(harvest_dir.glob("*_FULL_CONTEXT.md"))
        logger.info(f"Found {len(markdown_files)} markdown files to process")

        # Process papers in parallel (OPTIMIZED for 3-5x speedup!)
        from models import Paper
        from pipeline.extraction import ExpertExtractor

        # Override model if specified
        if model:
            extractor = ExpertExtractor(models=[model], tier_threshold=tier_threshold)
            logger.info(f"Using custom model: {model}")
        else:
            extractor = ExpertExtractor(tier_threshold=tier_threshold)

        extractions = []

        def process_paper_file(md_file):
            """Process a single paper file (for parallel execution)"""
            import threading
            import time

            thread_id = threading.current_thread().name
            start_time = time.time()

            # Extract PMID from filename (format: PMID_12345678_FULL_CONTEXT.md)
            pmid_match = md_file.stem.split('_')[1] if '_' in md_file.stem else None

            if not pmid_match:
                logger.warning(f"[{thread_id}] Could not extract PMID from filename: {md_file.name}")
                return None

            logger.info(f"[{thread_id}] Processing PMID {pmid_match}...")

            # Rate limiting
            if rate_limit_delay > 0:
                time.sleep(rate_limit_delay)

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
                    # Save individual extraction (thread-safe: each writes to unique file)
                    output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid_match}.json"
                    with open(output_file, 'w') as f:
                        json.dump(extraction_result.extracted_data, f, indent=2)

                    num_variants = extraction_result.extracted_data.get(
                        'extraction_metadata', {}
                    ).get('total_variants_found', 0)

                    elapsed = time.time() - start_time
                    logger.info(f"[{thread_id}] ✓ Completed PMID {pmid_match} in {elapsed:.1f}s ({num_variants} variants)")
                    logger.info(f"[{thread_id}] ✓ Saved to: {output_file}")

                    return extraction_result
                else:
                    elapsed = time.time() - start_time
                    logger.warning(f"[{thread_id}] ✗ Extraction failed for PMID {pmid_match} after {elapsed:.1f}s: {extraction_result.error}")
                    return None

            except Exception as e:
                elapsed = time.time() - start_time
                logger.error(f"[{thread_id}] Error processing PMID {pmid_match} after {elapsed:.1f}s: {e}")
                return None

        # Process papers in parallel with ThreadPoolExecutor
        # LLM calls are I/O-bound, so threading provides excellent speedup

        # Determine worker count
        if max_workers:
            num_workers = max_workers
        else:
            # Auto-detect based on CPUs
            cpu_count = os.cpu_count() or 4
            num_workers = min(max_workers_cap, cpu_count)

        num_workers = min(num_workers, len(markdown_files)) if markdown_files else 1

        logger.info(f"🔧 Parallelism configured: {num_workers} workers (cap: {max_workers_cap}, rate limit: {rate_limit_delay}s)")

        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            # Submit all papers for processing
            future_to_file = {
                executor.submit(process_paper_file, md_file): md_file
                for md_file in markdown_files
            }

            # Collect results as they complete
            completed = 0
            total = len(markdown_files)

            for future in as_completed(future_to_file):
                md_file = future_to_file[future]
                completed += 1

                try:
                    result = future.result()
                    if result:
                        extractions.append(result)

                    # Progress indicator
                    logger.info(f"📊 Progress: {completed}/{total} papers ({completed/total*100:.1f}%)")
                except Exception as e:
                    logger.error(f"⚠ Failed to process {md_file.name}: {e}")

    else:
        # Skip extraction - harvest only mode
        logger.info("\n⏭️  STEP 3: Skipping extraction (--skip-extraction flag set)")
        logger.info("Papers have been downloaded and are ready for later extraction")
        extractions = []
        extraction_dir = None

    # ============================================================================
    # STEP 4: Aggregate Penetrance Data (Optional)
    # ============================================================================
    if not skip_extraction and extraction_dir:
        logger.info("\n📊 STEP 4: Aggregating penetrance data across papers...")

        from pipeline.aggregation import aggregate_penetrance

        penetrance_summary_file = output_path / f"{gene_symbol}_penetrance_summary.json"
        penetrance_summary = aggregate_penetrance(
            extraction_dir=extraction_dir,
            gene_symbol=gene_symbol,
            output_file=penetrance_summary_file
        )

        logger.info(f"✓ Aggregated penetrance data for {penetrance_summary['total_variants']} variants")
        logger.info(f"✓ Saved penetrance summary to: {penetrance_summary_file}")
    else:
        logger.info("\n⏭️  STEP 4: Skipping aggregation (no extraction performed)")
        penetrance_summary = {"total_variants": 0, "variants": []}

    # ============================================================================
    # STEP 5: Migrate to SQLite (Optional)
    # ============================================================================
    migration_stats = None

    if migrate_to_sqlite and extraction_dir:
        logger.info("\n🗄️  STEP 5: Migrating to SQLite database...")

        try:
            from migrate_to_sqlite import (
                create_database_schema,
                migrate_extraction_directory as migrate_dir,
                cleanup_data_directory
            )

            # Determine database path
            if db_path:
                db_file = db_path
            else:
                db_file = f"{gene_symbol}.db"

            logger.info(f"Database: {db_file}")

            # Create schema and migrate
            conn = create_database_schema(db_file)
            migration_stats = migrate_dir(conn, extraction_dir)

            logger.info(f"✓ Migrated {migration_stats['successful']}/{migration_stats['total_files']} papers to {db_file}")

            # Optional cleanup
            if cleanup_after_migrate:
                logger.info("🧹 Running cleanup and archival...")
                cleanup_results = cleanup_data_directory(
                    output_path,
                    delete_empty_dirs=True,
                    archive_pmc=True,
                    delete_pmc_after_archive=True,
                    dry_run=False
                )
                logger.info(f"✓ Archived {len(cleanup_results['archives_created'])} directories")

            conn.close()

        except Exception as e:
            logger.error(f"⚠️  Migration failed: {e}", exc_info=True)
            migration_stats = {"error": str(e)}

    elif migrate_to_sqlite and not extraction_dir:
        logger.info("\n⏭️  STEP 5: Skipping SQLite migration (no extraction performed)")
    else:
        logger.info("\n⏭️  STEP 5: Skipping SQLite migration (use --migrate-to-sqlite to enable)")

    # ============================================================================
    # STEP 6: Compile Results and Statistics
    # ============================================================================
    logger.info("\n📊 STEP 6: Compiling results and statistics...")

    total_variants = 0
    for extraction in extractions:
        if extraction.extracted_data:
            total_variants += extraction.extracted_data.get(
                'extraction_metadata', {}
            ).get('total_variants_found', 0)

    # Calculate total carriers and affected from penetrance summary
    total_carriers = sum(
        v.get("aggregated_penetrance", {}).get("total_carriers", 0) or 0
        for v in penetrance_summary.get("variants", [])
    )
    total_affected = sum(
        v.get("aggregated_penetrance", {}).get("affected", 0) or 0
        for v in penetrance_summary.get("variants", [])
    )

    # Create summary report
    summary = {
        "gene_symbol": gene_symbol,
        "workflow_timestamp": datetime.now().isoformat(),
        "statistics": {
            "pmids_discovered": len(pmids),
            "papers_downloaded": num_downloaded,
            "papers_extracted": len(extractions),
            "total_variants_found": total_variants,
            "variants_with_penetrance_data": penetrance_summary["total_variants"],
            "total_carriers_observed": total_carriers,
            "total_affected_carriers": total_affected,
            "success_rate": f"{len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%"
        },
        "extraction": {
            "skip_extraction": skip_extraction,
            "model_used": model or "default (from config)",
            "parallel_workers": num_workers if not skip_extraction else 0,
            "rate_limit_delay": rate_limit_delay
        },
        "output_locations": {
            "pmid_list": str(pmids_file),
            "full_text_papers": str(harvest_dir),
            "extractions": str(extraction_dir) if extraction_dir else None,
            "penetrance_summary": str(penetrance_summary_file) if not skip_extraction else None
        },
        "penetrance_validation": {
            "errors": penetrance_summary.get("validation", {}).get("error_count", 0),
            "warnings": penetrance_summary.get("validation", {}).get("warning_count", 0)
        }
    }

    # Add migration stats if performed
    if migration_stats:
        summary["migration"] = {
            "database_path": db_file if migrate_to_sqlite and extraction_dir else None,
            "papers_migrated": migration_stats.get("successful", 0),
            "migration_errors": migration_stats.get("failed", 0),
            "cleanup_performed": cleanup_after_migrate
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
    logger.info(f"Variants with penetrance data: {penetrance_summary['total_variants']}")
    logger.info(f"Total carriers observed: {total_carriers}")
    logger.info(f"Total affected carriers: {total_affected}")
    logger.info(f"Success rate: {len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%")
    logger.info(f"\nAll outputs saved to: {output_path}")
    logger.info(f"Summary report: {summary_file}")
    logger.info(f"Penetrance summary: {penetrance_summary_file}")
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
  python automated_workflow.py BRCA1 --email your@email.com

  # Extract data for SCN5A with custom limits
  python automated_workflow.py SCN5A --email your@email.com --max-pmids 200 --max-downloads 100

  # Quick test with small dataset
  python automated_workflow.py TP53 --email your@email.com --max-pmids 10 --max-downloads 5
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
    parser.add_argument("--tier-threshold", type=int, default=None,
                       help="If the first model finds fewer variants than this, the next model is tried (default: from .env TIER3_THRESHOLD or 1). Set to 0 to only use first model.")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose logging")
    parser.add_argument("--skip-extraction", action="store_true",
                       help="Skip LLM extraction (harvest papers only, no API key required)")
    parser.add_argument("--model", type=str, default=None,
                       help="LLM model to use (e.g., 'gpt-4o', 'claude-3-opus', 'anthropic/claude-3-5-sonnet-20241022')")
    parser.add_argument("--max-workers", type=int, default=None,
                       help="Number of parallel workers (default: auto-detect based on CPU)")
    parser.add_argument("--max-workers-cap", type=int, default=8,
                       help="Maximum cap for auto-detected workers (default: 8)")
    parser.add_argument("--rate-limit-delay", type=float, default=0.0,
                       help="Delay in seconds between API requests to avoid rate limits (default: 0)")
    parser.add_argument("--migrate-to-sqlite", action="store_true",
                       help="Migrate extraction results to SQLite database after workflow")
    parser.add_argument("--db-path", type=str, default=None,
                       help="SQLite database path (default: auto-detect as {GENE}.db)")
    parser.add_argument("--cleanup-after-migrate", action="store_true",
                       help="Archive PMC files after successful migration")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Check for API keys (soft warning with escape hatch)
    has_api_key = bool(os.getenv("AI_INTEGRATIONS_OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY"))

    if not has_api_key:
        if args.skip_extraction:
            logger.warning("⚠️  No OpenAI API key found - extraction will be skipped")
        else:
            logger.error("⚠️  ERROR: OpenAI API key not found!")
            logger.error("Please set AI_INTEGRATIONS_OPENAI_API_KEY or OPENAI_API_KEY")
            logger.error("Example: export AI_INTEGRATIONS_OPENAI_API_KEY='your-key-here'")
            logger.error("Or use --skip-extraction to run harvest-only workflow")
            sys.exit(1)

    # Get tier threshold from settings if not provided via CLI
    from config.settings import get_settings
    tier_threshold = args.tier_threshold
    if tier_threshold is None:
        settings = get_settings()
        tier_threshold = settings.tier3_threshold

    # Run workflow
    try:
        automated_variant_extraction_workflow(
            gene_symbol=args.gene,
            email=args.email,
            max_pmids=args.max_pmids,
            max_papers_to_download=args.max_downloads,
            output_dir=args.output,
            tier_threshold=tier_threshold,
            skip_extraction=args.skip_extraction,
            model=args.model,
            max_workers=args.max_workers,
            max_workers_cap=args.max_workers_cap,
            rate_limit_delay=args.rate_limit_delay,
            migrate_to_sqlite=args.migrate_to_sqlite,
            db_path=args.db_path,
            cleanup_after_migrate=args.cleanup_after_migrate,
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
