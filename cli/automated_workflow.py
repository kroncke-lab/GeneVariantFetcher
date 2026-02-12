#!/usr/bin/env python3
"""
Automated Workflow Entrypoint: From Gene to Variant Data

This script provides the production-ready, end-to-end automated workflow:
1. Fetch relevant PMIDs from PubMind, PubMed, and Europe PMC for a gene
2. Download full-text articles from PubMed Central
3. Extract individual-level variant and patient data
4. Save structured results to JSON and aggregate penetrance metrics

This module uses the shared step implementations from pipeline/steps.py,
ensuring consistency with the GUI workflow.
"""

import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from dotenv import load_dotenv

from cli.scout import run_scout

# Load environment variables from .env file
load_dotenv()

# Configure logging using centralized utility
from gui.checkpoint import CheckpointManager, JobCheckpoint, PipelineStep
from utils.logging_utils import get_logger, setup_logging
from utils.run_manifest import RunManifestManager

setup_logging(level=logging.INFO)
logger = get_logger(__name__)


def automated_variant_extraction_workflow(
    gene_symbol: str,
    email: str,
    output_dir: str,
    max_pmids: int = 100,
    max_papers_to_download: int = 50,
    tier_threshold: int = 1,
    use_clinical_triage: bool = False,
    auto_synonyms: bool = True,
    synonyms: list[str] | None = None,
    scout_first: bool = False,
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
        auto_synonyms: Automatically discover and use gene synonyms from NCBI Gene database (default: True).
        synonyms: List of manually specified gene synonyms to include in searches.
        scout_first: Run Data Scout before extraction to identify high-value data zones.
    """
    from config.settings import get_settings
    from pipeline.steps import (
        aggregate_data,
        discover_synonyms,
        download_fulltext,
        extract_variants,
        fetch_abstracts,
        fetch_pmids,
        filter_papers,
        migrate_to_sqlite,
    )

    # Setup output directory
    output_path = (
        Path(output_dir) / gene_symbol / datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    output_path.mkdir(parents=True, exist_ok=True)

    # Set up file logging
    log_file = output_path / f"{gene_symbol}_workflow.log"
    setup_logging(level=logging.INFO, log_file=log_file)
    logger.info(f"Logging to file: {log_file}")

    # Initialize run manifest
    run_manifest = RunManifestManager.create_for_workflow(gene_symbol, str(output_path))
    run_manifest.set_config(
        max_pmids=max_pmids,
        max_papers_to_download=max_papers_to_download,
        tier_threshold=tier_threshold,
        use_clinical_triage=use_clinical_triage,
        auto_synonyms=auto_synonyms,
        synonyms=synonyms or [],
        scout_first=scout_first,
    )

    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager()

    # Create unique job ID for this run
    job_id = run_manifest.run_id

    # Create initial checkpoint
    checkpoint = JobCheckpoint(
        job_id=job_id,
        gene_symbol=gene_symbol,
        email=email,
        output_dir=str(output_path),
        max_pmids=max_pmids,
        max_papers_to_download=max_papers_to_download,
        tier_threshold=tier_threshold,
        use_clinical_triage=use_clinical_triage,
        auto_synonyms=auto_synonyms,
        synonyms=synonyms or [],
    )

    # Save initial checkpoint and manifest
    checkpoint_manager.save(checkpoint)
    run_manifest.save()

    settings = get_settings()

    logger.info("=" * 80)
    logger.info(f"AUTOMATED WORKFLOW FOR GENE: {gene_symbol}")
    logger.info("=" * 80)

    # Track statistics for final summary
    workflow_stats = {
        "gene_symbol": gene_symbol,
        "workflow_timestamp": datetime.now().isoformat(),
    }

    # =========================================================================
    # STEP 0: Discover Gene Synonyms (if enabled)
    # =========================================================================
    all_synonyms: list[str] = list(synonyms) if synonyms else []

    checkpoint.update_step(PipelineStep.DISCOVERING_SYNONYMS)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("discovering_synonyms")
    run_manifest.save()

    if auto_synonyms:
        logger.info("\n🔍 STEP 0: Discovering gene synonyms from NCBI Gene database...")

        synonym_result = discover_synonyms(
            gene_symbol=gene_symbol,
            email=email,
            existing_synonyms=all_synonyms,
            api_key=os.getenv("NCBI_API_KEY"),
        )

        if synonym_result.success:
            all_synonyms = synonym_result.data.get("synonyms", [])
            logger.info(
                f"✓ Discovered {synonym_result.stats.get('synonyms_found', 0)} gene synonyms"
            )
            if all_synonyms:
                logger.info(f"✓ Total synonyms to use: {', '.join(all_synonyms)}")

            # Update checkpoint with discovered synonyms
            checkpoint.synonyms = all_synonyms
            checkpoint.step_progress["did_synonyms"] = True
            checkpoint_manager.save(checkpoint)

        else:
            logger.warning(f"Failed to discover synonyms: {synonym_result.error}")
            logger.warning("Continuing without synonym expansion")
            run_manifest.add_warning(
                f"Synonym discovery failed: {synonym_result.error}"
            )

    elif all_synonyms:
        logger.info(
            f"\n🔍 Using {len(all_synonyms)} manually specified synonyms: {', '.join(all_synonyms)}"
        )
        checkpoint.synonyms = all_synonyms
        checkpoint_manager.save(checkpoint)

    # =========================================================================
    # STEP 1: Fetch PMIDs from PubMind, PubMed, and Europe PMC
    # =========================================================================
    checkpoint.update_step(PipelineStep.FETCHING_PMIDS)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("fetching_pmids")
    run_manifest.save()

    logger.info(
        "\n📚 STEP 1: Discovering relevant papers from PubMind, PubMed, and Europe PMC..."
    )

    pmid_result = fetch_pmids(
        gene_symbol=gene_symbol,
        email=email,
        output_path=output_path,
        max_results=max_pmids,
        synonyms=all_synonyms if all_synonyms else None,
        use_pubmind=settings.use_pubmind,
        use_pubmed=settings.use_pubmed and not settings.pubmind_only,
        use_europepmc=settings.use_europepmc,
        api_key=os.getenv("NCBI_API_KEY"),
    )

    if not pmid_result.success:
        logger.error(f"PMID discovery failed: {pmid_result.error}")
        checkpoint.mark_failed(pmid_result.error, PipelineStep.FETCHING_PMIDS)
        checkpoint_manager.save(checkpoint)
        run_manifest.add_error(f"PMID discovery failed: {pmid_result.error}")
        run_manifest.finalize(success=False)
        return {"success": False, "error": pmid_result.error}

    pmids = pmid_result.data.get("pmids", [])
    logger.info(
        f"✓ Found {pmid_result.stats.get('pubmind_count', 0)} PubMind PMIDs, "
        f"{pmid_result.stats.get('pubmed_count', 0)} PubMed PMIDs, "
        f"and {pmid_result.stats.get('europepmc_count', 0)} Europe PMC PMIDs"
    )
    logger.info(f"✓ Using {len(pmids)} unique PMIDs after merging sources")

    if not pmids:
        logger.warning(f"No PMIDs found for {gene_symbol}. Workflow terminated.")
        checkpoint.mark_failed("No PMIDs found", PipelineStep.FETCHING_PMIDS)
        checkpoint_manager.save(checkpoint)
        run_manifest.add_warning("No PMIDs found for gene")
        run_manifest.finalize(success=False)
        return {"success": False, "error": "No PMIDs found"}

    workflow_stats["pmids_discovered"] = len(pmids)

    # Update checkpoint and manifest with discovered PMIDs
    checkpoint.discovered_pmids = pmids
    checkpoint_manager.save(checkpoint)
    run_manifest.update_statistics(pmids_discovered=len(pmids))

    # =========================================================================
    # STEP 1.5: Fetch Abstracts and Metadata
    # =========================================================================
    checkpoint.update_step(PipelineStep.FETCHING_ABSTRACTS)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("fetching_abstracts")
    run_manifest.save()

    logger.info(
        "\n📝 STEP 1.5: Fetching abstracts and metadata for discovered PMIDs..."
    )

    abstract_result = fetch_abstracts(
        pmids=pmids,
        output_path=output_path,
        email=email,
    )

    if not abstract_result.success:
        logger.error(f"Abstract fetch failed: {abstract_result.error}")
        checkpoint.mark_failed(abstract_result.error, PipelineStep.FETCHING_ABSTRACTS)
        checkpoint_manager.save(checkpoint)
        run_manifest.add_error(f"Abstract fetch failed: {abstract_result.error}")
        run_manifest.finalize(success=False)
        return {"success": False, "error": abstract_result.error}

    abstract_records = abstract_result.data.get("abstract_records", {})
    abstract_dir = abstract_result.data.get("abstract_dir")
    logger.info(
        f"✓ Saved abstracts for {len(abstract_records)} PMIDs to {abstract_dir}"
    )

    # =========================================================================
    # STEP 1.6: Filter Papers by Relevance
    # =========================================================================
    checkpoint.update_step(PipelineStep.FILTERING_PAPERS)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("filtering_papers")
    run_manifest.save()

    logger.info("\n🧹 STEP 1.6: Filtering papers by relevance before download...")

    filter_result = filter_papers(
        pmids=pmids,
        abstract_records=abstract_records,
        gene_symbol=gene_symbol,
        output_path=output_path,
        enable_tier1=settings.enable_tier1,
        enable_tier2=settings.enable_tier2,
        use_clinical_triage=use_clinical_triage,
        tier1_min_keywords=settings.tier1_min_keywords,
        tier2_confidence_threshold=settings.tier2_confidence_threshold,
    )

    filtered_pmids = filter_result.data.get("filtered_pmids", [])
    dropped_pmids = filter_result.data.get("dropped_pmids", [])

    logger.info(
        f"Filtering complete: {len(filtered_pmids)} passed filters, "
        f"{len(dropped_pmids)} dropped before download"
    )

    workflow_stats["pmids_filtered_out"] = len(dropped_pmids)
    workflow_stats["pmids_passed_filters"] = len(filtered_pmids)

    # Update checkpoint and manifest with filtered results
    checkpoint.filtered_pmids = filtered_pmids
    checkpoint_manager.save(checkpoint)
    run_manifest.update_statistics(
        pmids_filtered_out=len(dropped_pmids), pmids_passed_filters=len(filtered_pmids)
    )

    # =========================================================================
    # STEP 2: Download Full-Text Papers from PMC
    # =========================================================================
    checkpoint.update_step(PipelineStep.DOWNLOADING_FULLTEXT)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("downloading_fulltext")
    run_manifest.save()

    logger.info("\n📥 STEP 2: Downloading full-text papers from PubMed Central...")

    download_result = download_fulltext(
        pmids=filtered_pmids,
        output_path=output_path,
        gene_symbol=gene_symbol,
        max_papers=max_papers_to_download,
        delay=2.0,
    )

    harvest_dir = download_result.data.get("harvest_dir")
    downloaded_pmids = download_result.data.get("downloaded_pmids", [])
    abstract_only_pmids = download_result.data.get("abstract_only_pmids", [])

    recovered_pmids = download_result.data.get("recovered_pmids", [])
    logger.info(f"✓ Successfully obtained {len(downloaded_pmids)} full-text papers")
    if recovered_pmids:
        logger.info(
            f"  ↳ {len(recovered_pmids)} recovered from prior runs, "
            f"{len(downloaded_pmids) - len(recovered_pmids)} freshly downloaded"
        )
    logger.info(
        f"✓ {len(abstract_only_pmids)} papers will use abstract-only extraction"
    )

    workflow_stats["papers_downloaded"] = len(downloaded_pmids)
    workflow_stats["papers_download_failed"] = len(abstract_only_pmids)

    # Update checkpoint and manifest
    checkpoint.downloaded_pmids = downloaded_pmids
    checkpoint_manager.save(checkpoint)
    run_manifest.update_statistics(
        papers_downloaded=len(downloaded_pmids),
        papers_download_failed=len(abstract_only_pmids),
    )
    run_manifest.update_output_locations(harvest_dir=str(harvest_dir))

    # =========================================================================
    # STEP 2.5 (Optional): Run Data Scout for Better Context
    # =========================================================================
    scout_manifest_path = None
    if scout_first:
        logger.info(
            "\n🔍 STEP 2.5: Running Data Scout to identify high-value data zones..."
        )

        scout_output_dir = output_path / "scout_output"
        scout_output_dir.mkdir(parents=True, exist_ok=True)

        try:
            scout_result = run_scout(
                input_path=str(harvest_dir),
                output_dir=str(scout_output_dir),
                gene=gene_symbol,
                min_relevance=0.1,
                max_zones=30,
            )

            scout_manifest_path = scout_output_dir / "scout_manifest.json"
            if scout_manifest_path.exists():
                logger.info(
                    f"✓ Data Scout completed. Manifest saved to: {scout_manifest_path}"
                )
                run_manifest.update_output_locations(
                    scout_manifest=str(scout_manifest_path)
                )
            else:
                logger.warning("Data Scout completed but no manifest generated")

        except Exception as e:
            logger.warning(
                f"Data Scout failed: {e}. Continuing with standard extraction."
            )
            run_manifest.add_warning(f"Data Scout failed: {e}")

    # =========================================================================
    # STEP 3: Extract Variant and Patient Data
    # =========================================================================
    checkpoint.update_step(PipelineStep.EXTRACTING_VARIANTS)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("extracting_variants")
    run_manifest.save()

    logger.info("\n🧬 STEP 3: Extracting variant and patient data using AI...")

    extraction_dir = output_path / "extractions"

    # If scout_first was used, pass the scout manifest for better context
    extract_result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol=gene_symbol,
        abstract_records=abstract_records,
        abstract_only_pmids=abstract_only_pmids,
        tier_threshold=tier_threshold,
        max_workers=8,
    )

    extractions = extract_result.data.get("extractions", [])
    extraction_failures = extract_result.data.get("failures", [])

    logger.info(f"✓ Extracted data from {len(extractions)} papers")
    if extraction_failures:
        logger.warning(f"✗ {len(extraction_failures)} extraction failures")

    workflow_stats["papers_extracted"] = len(extractions)
    workflow_stats["extraction_failures"] = len(extraction_failures)
    workflow_stats["total_variants_found"] = extract_result.stats.get(
        "total_variants", 0
    )

    # Update checkpoint and manifest
    checkpoint.extracted_pmids = [e.pmid for e in extractions if hasattr(e, "pmid")]
    checkpoint_manager.save(checkpoint)

    # Count abstract vs fulltext extractions
    abstract_extraction_count = sum(
        1
        for e in extractions
        if e.extracted_data
        and e.extracted_data.get("extraction_metadata", {}).get("abstract_only")
    )

    run_manifest.update_statistics(
        papers_extracted=len(extractions),
        papers_extraction_failed=len(extraction_failures),
        total_variants_found=extract_result.stats.get("total_variants", 0),
        papers_from_fulltext=len(extractions) - abstract_extraction_count,
        papers_from_abstract_only=abstract_extraction_count,
    )
    run_manifest.update_output_locations(extractions_dir=str(extraction_dir))

    # =========================================================================
    # STEP 4: Aggregate Penetrance Data
    # =========================================================================
    checkpoint.update_step(PipelineStep.AGGREGATING_DATA)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("aggregating_data")
    run_manifest.save()

    logger.info("\n📊 STEP 4: Aggregating penetrance data across papers...")

    aggregate_result = aggregate_data(
        extraction_dir=extraction_dir,
        gene_symbol=gene_symbol,
        output_path=output_path,
    )

    penetrance_summary = aggregate_result.data.get("summary", {})
    penetrance_summary_file = output_path / f"{gene_symbol}_penetrance_summary.json"

    logger.info(
        f"✓ Aggregated penetrance data for {aggregate_result.stats.get('variants_aggregated', 0)} variants"
    )
    logger.info(f"✓ Saved penetrance summary to: {penetrance_summary_file}")

    workflow_stats["variants_with_penetrance"] = aggregate_result.stats.get(
        "variants_aggregated", 0
    )

    # Update checkpoint and manifest
    checkpoint_manager.save(checkpoint)
    run_manifest.update_statistics(
        variants_with_penetrance_data=aggregate_result.stats.get(
            "variants_aggregated", 0
        )
    )
    run_manifest.update_output_locations(
        penetrance_summary=str(penetrance_summary_file)
    )

    # =========================================================================
    # STEP 5: Migrate to SQLite Database
    # =========================================================================
    checkpoint.update_step(PipelineStep.MIGRATING_DATABASE)
    checkpoint_manager.save(checkpoint)
    run_manifest.update_status("migrating_database")
    run_manifest.save()

    logger.info("\n💾 STEP 5: Migrating data to SQLite database...")

    db_path = output_path / f"{gene_symbol}.db"

    migrate_result = migrate_to_sqlite(
        extraction_dir=extraction_dir,
        db_path=db_path,
    )

    logger.info(
        f"✓ Migrated {migrate_result.stats.get('successful', 0)}/"
        f"{migrate_result.stats.get('total_files', 0)} extractions to SQLite"
    )
    logger.info(f"✓ Database saved to: {db_path}")

    workflow_stats["migration_successful"] = migrate_result.stats.get("successful", 0)
    workflow_stats["migration_failed"] = migrate_result.stats.get("failed", 0)

    # Update manifest with database location
    run_manifest.update_output_locations(sqlite_database=str(db_path))

    # =========================================================================
    # STEP 6: Compile Results and Statistics
    # =========================================================================
    logger.info("\n📊 STEP 6: Compiling results and statistics...")

    # Calculate totals from penetrance summary
    total_carriers = sum(
        v.get("aggregated_penetrance", {}).get("total_carriers", 0) or 0
        for v in penetrance_summary.get("variants", [])
    )
    total_affected = sum(
        v.get("aggregated_penetrance", {}).get("affected", 0) or 0
        for v in penetrance_summary.get("variants", [])
    )

    # Count abstract-only extractions
    abstract_extraction_count = sum(
        1
        for e in extractions
        if e.extracted_data
        and e.extracted_data.get("extraction_metadata", {}).get("abstract_only")
    )

    # Create final summary
    summary = {
        "gene_symbol": gene_symbol,
        "workflow_timestamp": workflow_stats["workflow_timestamp"],
        "statistics": {
            "pmids_discovered": workflow_stats.get("pmids_discovered", 0),
            "pmids_filtered_out": workflow_stats.get("pmids_filtered_out", 0),
            "pmids_passed_filters": workflow_stats.get("pmids_passed_filters", 0),
            "papers_downloaded": workflow_stats.get("papers_downloaded", 0),
            "papers_download_failed": workflow_stats.get("papers_download_failed", 0),
            "papers_extracted": workflow_stats.get("papers_extracted", 0),
            "papers_from_fulltext": len(extractions) - abstract_extraction_count,
            "papers_from_abstract_only": abstract_extraction_count,
            "papers_extraction_failed": workflow_stats.get("extraction_failures", 0),
            "total_variants_found": workflow_stats.get("total_variants_found", 0),
            "variants_with_penetrance_data": workflow_stats.get(
                "variants_with_penetrance", 0
            ),
            "total_carriers_observed": total_carriers,
            "total_affected_carriers": total_affected,
            "success_rate": f"{len(extractions) / len(pmids) * 100:.1f}%"
            if pmids
            else "0%",
        },
        "output_locations": {
            "pmid_list": str(output_path / f"{gene_symbol}_pmids.txt"),
            "pmid_status": str(output_path / "pmid_status"),
            "full_text_papers": str(harvest_dir),
            "extractions": str(extraction_dir),
            "penetrance_summary": str(penetrance_summary_file),
            "sqlite_database": str(db_path),
            "workflow_log": str(log_file),
        },
        "database_migration": {
            "successful": migrate_result.stats.get("successful", 0),
            "failed": migrate_result.stats.get("failed", 0),
            "total_files": migrate_result.stats.get("total_files", 0),
        },
        "penetrance_validation": {
            "errors": penetrance_summary.get("validation", {}).get("error_count", 0),
            "warnings": penetrance_summary.get("validation", {}).get(
                "warning_count", 0
            ),
        },
    }

    # =========================================================================
    # STEP 6: Compile Results, Save Manifest, and Update Checkpoints
    # =========================================================================
    logger.info("\n📊 STEP 6: Finalizing results and saving manifests...")

    # Calculate success rate
    success_rate = f"{len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%"

    summary_file = output_path / f"{gene_symbol}_workflow_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    # Finalize checkpoint and run manifest
    checkpoint.mark_completed()
    checkpoint_manager.save(checkpoint)

    # Final manifest updates
    run_manifest.update_statistics(
        **workflow_stats,
        success_rate=success_rate,
        total_carriers_observed=total_carriers,
        total_affected_carriers=total_affected,
    )
    run_manifest.update_output_locations(
        workflow_summary=str(summary_file), workflow_log=str(log_file)
    )
    manifest_file = run_manifest.finalize(success=True)

    # Print final summary
    logger.info("\n" + "=" * 80)
    logger.info("WORKFLOW COMPLETE!")
    logger.info("=" * 80)
    logger.info(f"Gene: {gene_symbol}")
    logger.info(f"Run ID: {job_id}")
    logger.info(f"PMIDs discovered: {workflow_stats.get('pmids_discovered', 0)}")
    logger.info(f"  - Filtered out: {workflow_stats.get('pmids_filtered_out', 0)}")
    logger.info(f"  - Passed filters: {workflow_stats.get('pmids_passed_filters', 0)}")
    logger.info(f"Papers downloaded: {workflow_stats.get('papers_downloaded', 0)}")
    logger.info(
        f"  - Download failures: {workflow_stats.get('papers_download_failed', 0)}"
    )
    logger.info(f"Papers with extractions: {len(extractions)}")
    logger.info(f"  - From full-text: {len(extractions) - abstract_extraction_count}")
    logger.info(f"  - From abstract only: {abstract_extraction_count}")
    logger.info(
        f"  - Extraction failures: {workflow_stats.get('extraction_failures', 0)}"
    )
    logger.info(
        f"Total variants found: {workflow_stats.get('total_variants_found', 0)}"
    )
    logger.info(
        f"Variants with penetrance data: {workflow_stats.get('variants_with_penetrance', 0)}"
    )
    logger.info(f"Total carriers observed: {total_carriers}")
    logger.info(f"Total affected carriers: {total_affected}")
    logger.info(f"Success rate: {success_rate}")
    logger.info(
        f"💾 Database migrated: {migrate_result.stats.get('successful', 0)}/"
        f"{migrate_result.stats.get('total_files', 0)} extractions"
    )
    logger.info(f"\nAll outputs saved to: {output_path}")
    logger.info(f"Summary report: {summary_file}")
    logger.info(f"Run manifest: {manifest_file}")
    logger.info(f"Penetrance summary: {penetrance_summary_file}")
    logger.info(f"SQLite database: {db_path}")
    logger.info(f"Workflow log: {log_file}")
    logger.info("=" * 80)

    # Clean up checkpoint since workflow completed successfully
    checkpoint_manager.delete(job_id)

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
        """,
    )

    parser.add_argument("gene", help="Gene symbol (e.g., BRCA1, SCN5A, TP53)")
    parser.add_argument(
        "--email", "-e", required=True, help="Your email for NCBI E-utilities"
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output directory for all data and analyses (required)",
    )
    parser.add_argument(
        "--max-pmids",
        type=int,
        default=100,
        help="Maximum PMIDs to fetch (default: 100)",
    )
    parser.add_argument(
        "--max-downloads",
        type=int,
        default=50,
        help="Maximum papers to download (default: 50)",
    )
    parser.add_argument(
        "--tier-threshold",
        type=int,
        default=None,
        help="If the first model finds fewer variants than this, the next model is tried (default: from .env TIER3_THRESHOLD or 1). Set to 0 to only use first model.",
    )
    parser.add_argument(
        "--clinical-triage",
        action="store_true",
        help="Use ClinicalDataTriageFilter for Tier 2 filtering instead of InternFilter",
    )
    parser.add_argument(
        "--auto-synonyms",
        action="store_true",
        help="Automatically discover and use gene synonyms from NCBI Gene database",
    )
    parser.add_argument(
        "--synonym",
        action="append",
        dest="synonyms",
        default=None,
        help="Manually specify gene synonym (can be used multiple times)",
    )
    parser.add_argument(
        "--scout-first",
        action="store_true",
        help="Run Data Scout before extraction to identify high-value data zones for better context",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Check for API keys
    if not os.getenv("OPENAI_API_KEY"):
        logger.error("⚠️  ERROR: OPENAI_API_KEY not found in environment!")
        logger.error("Please set OPENAI_API_KEY in your .env file")
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
            auto_synonyms=args.auto_synonyms,
            synonyms=args.synonyms,
            scout_first=args.scout_first,
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
