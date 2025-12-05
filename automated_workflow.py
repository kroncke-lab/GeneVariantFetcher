#!/usr/bin/env python3
"""
Automated Workflow Entrypoint: From Gene to Variant Data

This script provides the production-ready, end-to-end automated workflow:
1. Fetch relevant PMIDs from PubMind, PubMed, and Europe PMC for a gene
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

# Configure logging using centralized utility
from utils.logging_utils import setup_logging, get_logger
setup_logging(level=logging.INFO)
logger = get_logger(__name__)


def automated_variant_extraction_workflow(
    gene_symbol: str,
    email: str,
    output_dir: str,
    max_pmids: Path = 100,
    max_papers_to_download: Path = 50,
    tier_threshold: Path = 1,
    use_clinical_triage: os = False,
):
    """
    Complete automated workflow from gene symbol to extracted variant data.

    Args:
        gene_symbol: Gene to search for (e.g., "BRCA1", "SCN5A")
        email: Your email for NCBI E-utilities (required)
        output_dir: Directory to save all outputs (required)
        max_pmids: Maximum PMIDs to fetch from active sources (integer)
        max_papers_to_download: Maximum papers to download full-text (integer)
        tier_threshold: If the first model finds fewer variants than this, the next model is tried (integer).
        use_clinical_triage: Use ClinicalDataTriageFilter for Tier 2 instead of InternFilter.
    """
    from gene_literature.discovery import discover_pmids_for_gene
    from harvesting import PMCHarvester
    from config.settings import get_settings

    output_path = Path(output_dir) / gene_symbol / datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path.mkdir(parents=True, exist_ok=True)

    settings = get_settings()

    logger.info("="*80)
    logger.info(f"AUTOMATED WORKFLOW FOR GENE: {gene_symbol}")
    logger.info("="*80)

    # ============================================================================
    # STEP 1: Fetch PMIDs from PubMind, PubMed, and Europe PMC
    # ============================================================================
    logger.info(
        "\n📚 STEP 1: Discovering relevant papers from PubMind, PubMed, and Europe PMC..."
    )

    pubmind_pmids_file = output_path / f"{gene_symbol}_pmids_pubmind.txt"
    pubmed_pmids_file = output_path / f"{gene_symbol}_pmids_pubmed.txt"
    combined_pmids_file = output_path / f"{gene_symbol}_pmids.txt"

    use_pubmind = settings.use_pubmind
    use_pubmed = settings.use_pubmed and not settings.pubmind_only

    if settings.pubmind_only and not settings.use_pubmind:
        logger.warning("PUBMIND_ONLY is set but USE_PUBMIND is false; enabling PubMind to avoid empty discovery.")
        use_pubmind = True

    if not use_pubmind and not use_pubmed:
        logger.warning("No PMID sources enabled; defaulting to PubMind to proceed.")
        use_pubmind = True

    logger.info(
        "PMID source selection -> PubMind: %s | PubMed keyword search: %s",
        "on" if use_pubmind else "off",
        "on" if use_pubmed else "off",
    )

    pmid_discovery = discover_pmids_for_gene(
        gene_symbol=gene_symbol,
        email=email,
        max_results=max_pmids,
        pubmind_output=pubmind_pmids_file,
        pubmed_output=pubmed_pmids_file,
        combined_output=combined_pmids_file,
        api_key=os.getenv("NCBI_API_KEY"),
        use_pubmind=use_pubmind,
        use_pubmed=use_pubmed,
    )
    pmids = pmid_discovery.combined_pmids

    logger.info(
        "✓ Found %d PubMind PMIDs, %d PubMed PMIDs, and %d Europe PMC PMIDs",
        len(pmid_discovery.pubmind_pmids),
        len(pmid_discovery.pubmed_pmids),
        len(pmid_discovery.europepmc_pmids),
    )
    logger.info("✓ Using %d unique PMIDs after merging sources", len(pmids))
    logger.info(f"✓ Saved combined PMID list to: {combined_pmids_file}")

    if not pmids:
        logger.warning(f"No PMIDs found for {gene_symbol}. Workflow terminated.")
        return {"success": False, "error": "No PMIDs found"}

    # ============================================================================
    # STEP 1.5: Fetch Abstracts and Metadata
    # ============================================================================
    logger.info("\n📝 STEP 1.5: Fetching abstracts and metadata for discovered PMIDs...")

    from harvesting.abstracts import fetch_and_save_abstracts

    abstract_dir = output_path / "abstract_json"
    abstract_records = fetch_and_save_abstracts(
        pmids=pmids,
        output_dir=str(abstract_dir),
        email=email,
    )

    logger.info(
        "✓ Saved abstracts for %d PMIDs to %s", len(abstract_records), abstract_dir
    )

    # Filter abstracts before attempting full-text downloads
    logger.info("\n🧹 STEP 1.6: Filtering papers by relevance before download...")

    from pipeline.filters import KeywordFilter, InternFilter, ClinicalDataTriageFilter
    from utils.models import Paper, FilterDecision, FilterResult, FilterTier

    keyword_filter = KeywordFilter(min_keyword_matches=settings.tier1_min_keywords)
    tier2_filter = (
        ClinicalDataTriageFilter()
        if use_clinical_triage
        else InternFilter(confidence_threshold=settings.tier2_confidence_threshold)
    )
    tier2_filter_name = "ClinicalDataTriageFilter" if use_clinical_triage else "InternFilter"

    filtered_pmids = []
    dropped_pmids = []

    for pmid in pmids:
        record_path = abstract_records.get(pmid)

        if not record_path:
            logger.warning("PMID %s has no saved abstract JSON; dropping from download queue", pmid)
            dropped_pmids.append((pmid, "Missing abstract JSON"))
            continue

        record_path = Path(record_path)

        if not record_path.exists():
            logger.warning("PMID %s has no saved abstract JSON; dropping from download queue", pmid)
            dropped_pmids.append((pmid, "Missing abstract JSON"))
            continue

        try:
            with record_path.open("r", encoding="utf-8") as f:
                record = json.load(f)
        except Exception as e:
            logger.error("Failed to read abstract for PMID %s: %s", pmid, e)
            dropped_pmids.append((pmid, "Abstract JSON read error"))
            continue

        metadata = record.get("metadata", {})
        paper = Paper(
            pmid=pmid,
            title=metadata.get("title"),
            abstract=record.get("abstract"),
            authors=metadata.get("authors"),
            journal=metadata.get("journal"),
            publication_date=metadata.get("year"),
            gene_symbol=gene_symbol,
            source="PubMed",
        )

        if settings.enable_tier1:
            tier1_result = keyword_filter.filter(paper)
            if tier1_result.decision is not FilterDecision.PASS:
                logger.info(
                    "PMID %s dropped at Tier 1 (KeywordFilter): %s",
                    pmid,
                    tier1_result.reason,
                )
                dropped_pmids.append((pmid, tier1_result.reason))
                continue

        if settings.enable_tier2:
            if use_clinical_triage:
                triage_result = tier2_filter.triage_paper(paper, gene_symbol)
                decision = FilterDecision.PASS if triage_result.get("decision") == "KEEP" else FilterDecision.FAIL
                confidence = triage_result.get("confidence")
                reason = triage_result.get("reason", "No reason provided")

                if (
                    decision is FilterDecision.PASS
                    and confidence is not None
                    and confidence < settings.tier2_confidence_threshold
                ):
                    decision = FilterDecision.FAIL
                    reason = (
                        f"Low confidence ({confidence:.2f} < {settings.tier2_confidence_threshold}): "
                        f"{reason}"
                    )

                tier2_result = FilterResult(
                    decision=decision,
                    tier=FilterTier.TIER_2_INTERN,
                    reason=reason,
                    pmid=pmid,
                    confidence=confidence,
                    metadata={"model": getattr(tier2_filter, "model", None)},
                )
            else:
                tier2_result = tier2_filter.filter(paper)

            if tier2_result.decision is not FilterDecision.PASS:
                logger.info(
                    "PMID %s dropped at Tier 2 (%s): %s",
                    pmid,
                    tier2_filter_name,
                    tier2_result.reason,
                )
                dropped_pmids.append((pmid, tier2_result.reason))
                continue

        filtered_pmids.append(pmid)

    logger.info(
        "Filtering complete: %d passed filters, %d dropped before download",
        len(filtered_pmids),
        len(dropped_pmids),
    )

    if dropped_pmids:
        logger.debug(
            "Dropped PMIDs and reasons: %s",
            "; ".join([f"{pmid} ({reason})" for pmid, reason in dropped_pmids]),
        )

    # ============================================================================
    # STEP 2: Download Full-Text Papers from PMC
    # ============================================================================
    logger.info("\n📥 STEP 2: Downloading full-text papers from PubMed Central...")

    harvest_dir = output_path / "pmc_fulltext"
    harvester = PMCHarvester(output_dir=str(harvest_dir))

    # Limit downloads to avoid excessive processing time
    pmids_to_download = filtered_pmids[:max_papers_to_download]
    logger.info(
        f"Downloading up to {len(pmids_to_download)} papers after filtering (out of {len(pmids)} total PMIDs)"
    )

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

    # Process papers in parallel (OPTIMIZED for 3-5x speedup!)
    from utils.models import Paper
    from pipeline.extraction import ExpertExtractor

    extractor = ExpertExtractor(tier_threshold=tier_threshold)
    extractions = []

    def process_paper_file(md_file):
        """Process a single paper file (for parallel execution)"""
        from utils.pmid_utils import extract_pmid_from_filename

        # Extract PMID from filename using shared utility
        pmid_match = extract_pmid_from_filename(md_file)

        if not pmid_match:
            logger.warning(f"Could not extract PMID from filename: {md_file.name}")
            return None

        logger.info(f"Processing PMID {pmid_match}...")

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

                logger.info(f"✓ Extracted {num_variants} variants from PMID {pmid_match}")
                logger.info(f"✓ Saved to: {output_file}")

                return extraction_result
            else:
                logger.warning(f"✗ Extraction failed for PMID {pmid_match}: {extraction_result.error}")
                return None

        except Exception as e:
            logger.error(f"Error processing PMID {pmid_match}: {e}")
            return None

    # Process papers in parallel with ThreadPoolExecutor
    # LLM calls are I/O-bound, so threading provides excellent speedup
    max_workers = min(8, len(markdown_files)) if markdown_files else 1

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all papers for processing
        future_to_file = {
            executor.submit(process_paper_file, md_file): md_file
            for md_file in markdown_files
        }

        # Collect results as they complete
        completed = 0
        for future in as_completed(future_to_file):
            md_file = future_to_file[future]
            completed += 1

            try:
                result = future.result()
                if result:
                    extractions.append(result)
                logger.info(f"✓ Completed {completed}/{len(markdown_files)} papers")
            except Exception as e:
                logger.error(f"⚠ Failed to process {md_file.name}: {e}")

    # ============================================================================
    # STEP 4: Aggregate Penetrance Data
    # ============================================================================
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

    # ============================================================================
    # STEP 5: Migrate to SQLite Database
    # ============================================================================
    logger.info("\n💾 STEP 5: Migrating data to SQLite database...")

    from harvesting.migrate_to_sqlite import create_database_schema, migrate_extraction_directory

    # Create gene-specific database
    db_path = output_path / f"{gene_symbol}.db"
    conn = create_database_schema(str(db_path))

    # Migrate all extraction files
    migration_stats = migrate_extraction_directory(conn, extraction_dir)
    conn.close()

    logger.info(f"✓ Migrated {migration_stats['successful']}/{migration_stats['total_files']} extractions to SQLite")
    logger.info(f"✓ Database saved to: {db_path}")

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
        "output_locations": {
            "pmid_list": str(combined_pmids_file),
            "full_text_papers": str(harvest_dir),
            "extractions": str(extraction_dir),
            "penetrance_summary": str(penetrance_summary_file),
            "sqlite_database": str(db_path)
        },
        "database_migration": {
            "successful": migration_stats["successful"],
            "failed": migration_stats["failed"],
            "total_files": migration_stats["total_files"]
        },
        "penetrance_validation": {
            "errors": penetrance_summary.get("validation", {}).get("error_count", 0),
            "warnings": penetrance_summary.get("validation", {}).get("warning_count", 0)
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
    logger.info(f"Variants with penetrance data: {penetrance_summary['total_variants']}")
    logger.info(f"Total carriers observed: {total_carriers}")
    logger.info(f"Total affected carriers: {total_affected}")
    logger.info(f"Success rate: {len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%")
    logger.info(f"\n💾 Database migrated: {migration_stats['successful']}/{migration_stats['total_files']} extractions")
    logger.info(f"\nAll outputs saved to: {output_path}")
    logger.info(f"Summary report: {summary_file}")
    logger.info(f"Penetrance summary: {penetrance_summary_file}")
    logger.info(f"SQLite database: {db_path}")
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
  python automated_workflow.py BRCA1 --email your@email.com --output /path/to/output

  # Extract data for SCN5A with custom limits
  python automated_workflow.py SCN5A --email your@email.com --output ./results --max-pmids 200 --max-downloads 100

  # Quick test with small dataset
  python automated_workflow.py TP53 --email your@email.com --output ./test_output --max-pmids 10 --max-downloads 5
        """
    )

    parser.add_argument("gene", help="Gene symbol (e.g., BRCA1, SCN5A, TP53)")
    parser.add_argument("--email", "-e", required=True, help="Your email for NCBI E-utilities")
    parser.add_argument("--output", "-o", required=True,
                       help="Output directory for all data and analyses (required)")
    parser.add_argument("--max-pmids", type=Path, default=100,
                       help="Maximum PMIDs to fetch (default: 100)")
    parser.add_argument("--max-downloads", type=Path, default=50,
                       help="Maximum papers to download (default: 50)")
    parser.add_argument("--tier-threshold", type=Path, default=None,
                       help="If the first model finds fewer variants than this, the next model is tried (default: from .env TIER3_THRESHOLD or 1). Set to 0 to only use first model.")
    parser.add_argument("--clinical-triage", action="store_true",
                       help="Use ClinicalDataTriageFilter for Tier 2 filtering instead of InternFilter")
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
            output_dir=args.output,
            max_pmids=args.max_pmids,
            max_papers_to_download=args.max_downloads,
            tier_threshold=tier_threshold,
            use_clinical_triage=args.clinical_triage,
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
