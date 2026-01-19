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
import csv
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
    max_pmids: int = 100,
    max_papers_to_download: int = 50,
    tier_threshold: int = 1,
    use_clinical_triage: bool = False,
    auto_synonyms: bool = False,
    synonyms: list[str] | None = None,
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
        auto_synonyms: Automatically discover and use gene synonyms from NCBI Gene database.
        synonyms: List of manually specified gene synonyms to include in searches.
    """
    from gene_literature.discovery import discover_pmids_for_gene
    from harvesting import PMCHarvester
    from config.settings import get_settings

    output_path = Path(output_dir) / gene_symbol / datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path.mkdir(parents=True, exist_ok=True)

    # Set up file logging in the output directory for troubleshooting
    log_file = output_path / "workflow.log"
    setup_logging(level=logging.INFO, log_file=log_file)
    logger.info(f"Logging to file: {log_file}")

    settings = get_settings()

    logger.info("="*80)
    logger.info(f"AUTOMATED WORKFLOW FOR GENE: {gene_symbol}")
    logger.info("="*80)

    # ============================================================================
    # STEP 0: Discover Gene Synonyms (if enabled)
    # ============================================================================
    all_synonyms: list[str] = list(synonyms) if synonyms else []

    if auto_synonyms:
        logger.info("\n🔍 STEP 0: Discovering gene synonyms from NCBI Gene database...")

        from gene_literature.synonym_finder import SynonymFinder, automatic_synonym_selection

        synonym_finder = SynonymFinder(
            email=email,
            api_key=os.getenv("NCBI_API_KEY"),
        )

        try:
            found_synonyms = synonym_finder.find_gene_synonyms(
                gene_symbol,
                include_other_designations=False,  # Skip verbose designations
            )

            # Use automatic selection for batch mode
            auto_selected = automatic_synonym_selection(
                gene_symbol,
                found_synonyms,
                include_official=True,
                include_aliases=True,
                include_other_designations=False,
                only_relevant=False,  # No LLM checking in automated mode
            )

            # Merge with manually provided synonyms (avoid duplicates)
            existing_set = set(s.lower() for s in all_synonyms)
            for syn in auto_selected:
                if syn.lower() not in existing_set and syn.lower() != gene_symbol.lower():
                    all_synonyms.append(syn)
                    existing_set.add(syn.lower())

            logger.info("✓ Discovered %d gene synonyms", len(auto_selected))
            if all_synonyms:
                logger.info("✓ Total synonyms to use: %s", ", ".join(all_synonyms))

        except Exception as e:
            logger.warning("Failed to discover synonyms: %s", e)
            logger.warning("Continuing without synonym expansion")

    elif all_synonyms:
        logger.info("\n🔍 Using %d manually specified synonyms: %s", len(all_synonyms), ", ".join(all_synonyms))

    # ============================================================================
    # STEP 1: Fetch PMIDs from PubMind, PubMed, and Europe PMC
    # ============================================================================
    logger.info(
        "\n📚 STEP 1: Discovering relevant papers from PubMind, PubMed, and Europe PMC..."
    )

    pubmind_pmids_file = output_path / f"{gene_symbol}_pmids_pubmind.txt"
    pubmed_pmids_file = output_path / f"{gene_symbol}_pmids_pubmed.txt"
    combined_pmids_file = output_path / f"{gene_symbol}_pmids.txt"

    # Log which sources are enabled (determined by settings)
    use_pubmind = settings.use_pubmind
    use_pubmed = settings.use_pubmed and not settings.pubmind_only

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
        settings=settings,
        synonyms=all_synonyms if all_synonyms else None,
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

    # Persist dropped PMIDs to disk for debugging and analysis
    pmid_status_dir = output_path / "pmid_status"
    pmid_status_dir.mkdir(parents=True, exist_ok=True)
    filtered_out_file = pmid_status_dir / "filtered_out.csv"
    with open(filtered_out_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["PMID", "Reason"])
        writer.writerows(dropped_pmids)
    logger.info(f"✓ Wrote {len(dropped_pmids)} filtered-out PMIDs to {filtered_out_file}")

    # NOTE: Abstract extraction is now integrated into the main extraction pipeline (Step 3)
    # Papers that pass filters but fail download will have their abstracts extracted
    # using the same ExpertExtractor, ensuring all data flows into the final SQLite database.

    # ============================================================================
    # STEP 2: Download Full-Text Papers from PMC
    # ============================================================================
    logger.info("\n📥 STEP 2: Downloading full-text papers from PubMed Central...")

    harvest_dir = output_path / "pmc_fulltext"
    harvester = PMCHarvester(output_dir=str(harvest_dir), gene_symbol=gene_symbol)

    # Limit downloads to avoid excessive processing time
    pmids_to_download = filtered_pmids[:max_papers_to_download]
    logger.info(
        f"Downloading up to {len(pmids_to_download)} papers after filtering (out of {len(pmids)} total PMIDs)"
    )

    harvester.harvest(pmids_to_download, delay=2.0)

    # Check how many were successfully downloaded
    success_log = harvest_dir / "successful_downloads.csv"
    successfully_downloaded_pmids = set()
    if success_log.exists():
        import pandas as pd
        successful_downloads = pd.read_csv(success_log)
        num_downloaded = len(successful_downloads)
        successfully_downloaded_pmids = set(str(p) for p in successful_downloads["PMID"].tolist())
        logger.info(f"✓ Successfully downloaded {num_downloaded} full-text papers")
    else:
        num_downloaded = 0
        logger.warning("No papers were successfully downloaded from PMC")

    # ============================================================================
    # STEP 2.5: Identify Papers for Abstract-Only Extraction
    # ============================================================================
    # Papers that passed the relevance filter but failed to download will be
    # extracted using their abstracts in Step 3.
    logger.info("\n📌 STEP 2.5: Identifying papers that need abstract-only extraction...")

    # Load paywalled/missing log if it exists
    paywalled_log = harvest_dir / "paywalled_missing.csv"
    paywalled_pmids: dict[str, str] = {}  # PMID -> Reason
    if paywalled_log.exists():
        try:
            paywalled_df = pd.read_csv(paywalled_log)
            for _, row in paywalled_df.iterrows():
                paywalled_pmids[str(row["PMID"])] = row.get("Reason", "Unknown")
        except Exception as e:
            logger.warning(f"Could not read paywalled_missing.csv: {e}")

    # Identify papers that passed filters but need abstract-only extraction
    abstract_only_pmids = []
    for pmid in pmids_to_download:
        pmid_str = str(pmid)
        if pmid_str not in successfully_downloaded_pmids:
            abstract_only_pmids.append(pmid_str)

    logger.info(
        f"✓ {len(successfully_downloaded_pmids)} papers downloaded successfully, "
        f"{len(abstract_only_pmids)} will use abstract-only extraction"
    )

    # ============================================================================
    # STEP 3: Extract Variant and Patient Data
    # ============================================================================
    logger.info("\n🧬 STEP 3: Extracting variant and patient data using AI...")

    extraction_dir = output_path / "extractions"
    extraction_dir.mkdir(exist_ok=True)

    # Get list of downloaded markdown files - prefer DATA_ZONES.md over FULL_CONTEXT.md
    # First, find all available PMIDs from either file type
    data_zones_files = {f.name.replace("_DATA_ZONES.md", ""): f for f in harvest_dir.glob("*_DATA_ZONES.md")}
    full_context_files = {f.name.replace("_FULL_CONTEXT.md", ""): f for f in harvest_dir.glob("*_FULL_CONTEXT.md")}

    # Merge: prefer DATA_ZONES.md when available, fall back to FULL_CONTEXT.md
    all_pmids = set(data_zones_files.keys()) | set(full_context_files.keys())
    markdown_files = []
    for pmid in all_pmids:
        if pmid in data_zones_files:
            markdown_files.append(data_zones_files[pmid])
        elif pmid in full_context_files:
            markdown_files.append(full_context_files[pmid])

    logger.info(f"Found {len(data_zones_files)} DATA_ZONES.md and {len(full_context_files)} FULL_CONTEXT.md files")
    logger.info(f"Processing {len(markdown_files)} unique papers with full-text (preferring DATA_ZONES.md)")

    # Also prepare papers that need abstract-only extraction
    abstract_only_papers = []
    for pmid in abstract_only_pmids:
        record_path = abstract_records.get(pmid)
        if record_path and Path(record_path).exists():
            abstract_only_papers.append((pmid, record_path))

    if abstract_only_papers:
        logger.info(f"Also processing {len(abstract_only_papers)} papers using abstract-only extraction")

    # Process papers in parallel (OPTIMIZED for 3-5x speedup!)
    from utils.models import Paper
    from pipeline.extraction import ExpertExtractor

    extractor = ExpertExtractor(tier_threshold=tier_threshold, fulltext_dir=str(harvest_dir))
    extractions = []

    def process_paper_file(md_file):
        """Process a single paper file (for parallel execution).

        Returns:
            tuple: (extraction_result or None, failure_info or None)
                   where failure_info is (pmid, error_message)
        """
        from utils.pmid_utils import extract_pmid_from_filename

        # Extract PMID from filename using shared utility
        pmid_match = extract_pmid_from_filename(md_file)

        if not pmid_match:
            logger.warning(f"Could not extract PMID from filename: {md_file.name}")
            return (None, (md_file.name, "Could not extract PMID from filename"))

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

                return (extraction_result, None)
            else:
                logger.warning(f"✗ Extraction failed for PMID {pmid_match}: {extraction_result.error}")
                return (None, (pmid_match, extraction_result.error))

        except Exception as e:
            logger.error(f"Error processing PMID {pmid_match}: {e}")
            return (None, (pmid_match, str(e)))

    # Function to process abstract-only papers
    def process_abstract_paper(pmid: str, record_path: str):
        """Process a paper using only its abstract (for papers that failed download).

        Returns:
            tuple: (extraction_result or None, failure_info or None)
                   where failure_info is (pmid, error_message)
        """
        logger.info(f"Processing PMID {pmid} (ABSTRACT ONLY)...")

        try:
            with open(record_path, "r", encoding="utf-8") as f:
                record = json.load(f)

            metadata = record.get("metadata", {})
            abstract_text = record.get("abstract")

            if not abstract_text:
                return (None, (pmid, "No abstract available"))

            # Create paper object with abstract as the text source
            paper = Paper(
                pmid=pmid,
                title=metadata.get("title", f"Paper {pmid}"),
                abstract=abstract_text,
                full_text=None,  # No full text - extraction will use abstract
                gene_symbol=gene_symbol
            )

            # Extract data using abstract
            extraction_result = extractor.extract(paper)

            if extraction_result.success:
                # Mark as abstract-only extraction
                if extraction_result.extracted_data:
                    if "extraction_metadata" not in extraction_result.extracted_data:
                        extraction_result.extracted_data["extraction_metadata"] = {}
                    extraction_result.extracted_data["extraction_metadata"]["source_type"] = "abstract_only"
                    extraction_result.extracted_data["extraction_metadata"]["abstract_only"] = True

                # Save individual extraction
                output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid}.json"
                with open(output_file, 'w') as f:
                    json.dump(extraction_result.extracted_data, f, indent=2)

                num_variants = extraction_result.extracted_data.get(
                    'extraction_metadata', {}
                ).get('total_variants_found', 0)

                logger.info(f"✓ Extracted {num_variants} variants from PMID {pmid} (ABSTRACT ONLY)")
                logger.info(f"✓ Saved to: {output_file}")

                return (extraction_result, None)
            else:
                logger.warning(f"✗ Abstract extraction failed for PMID {pmid}: {extraction_result.error}")
                return (None, (pmid, f"Abstract extraction: {extraction_result.error}"))

        except Exception as e:
            logger.error(f"Error processing abstract for PMID {pmid}: {e}")
            return (None, (pmid, str(e)))

    # Process papers in parallel with ThreadPoolExecutor
    # LLM calls are I/O-bound, so threading provides excellent speedup
    total_papers = len(markdown_files) + len(abstract_only_papers)
    max_workers = min(8, total_papers) if total_papers > 0 else 1

    extraction_failures = []
    abstract_extraction_count = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all full-text papers for processing
        future_to_source = {}
        for md_file in markdown_files:
            future = executor.submit(process_paper_file, md_file)
            future_to_source[future] = ("fulltext", md_file)

        # Submit all abstract-only papers for processing
        for pmid, record_path in abstract_only_papers:
            future = executor.submit(process_abstract_paper, pmid, record_path)
            future_to_source[future] = ("abstract", pmid)

        # Collect results as they complete
        completed = 0
        for future in as_completed(future_to_source):
            source_type, source_info = future_to_source[future]
            completed += 1

            try:
                result, failure_info = future.result()
                if result:
                    extractions.append(result)
                    if source_type == "abstract":
                        abstract_extraction_count += 1
                if failure_info:
                    extraction_failures.append(failure_info)
                logger.info(f"✓ Completed {completed}/{total_papers} papers")
            except Exception as e:
                if source_type == "fulltext":
                    logger.error(f"⚠ Failed to process {source_info.name}: {e}")
                    from utils.pmid_utils import extract_pmid_from_filename
                    pmid = extract_pmid_from_filename(source_info) or source_info.name
                else:
                    logger.error(f"⚠ Failed to process abstract for {source_info}: {e}")
                    pmid = source_info
                extraction_failures.append((pmid, f"Executor error: {e}"))

    # Persist extraction failures to disk
    extraction_failures_file = pmid_status_dir / "extraction_failures.csv"
    with open(extraction_failures_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["PMID", "Error"])
        writer.writerows(extraction_failures)
    logger.info(f"✓ Wrote {len(extraction_failures)} extraction failures to {extraction_failures_file}")

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

    # Count download failures from paywalled_missing.csv if it exists
    paywalled_log = harvest_dir / "paywalled_missing.csv"
    num_download_failures = 0
    if paywalled_log.exists():
        import pandas as pd
        try:
            paywalled_df = pd.read_csv(paywalled_log)
            num_download_failures = len(paywalled_df)
        except Exception:
            pass

    # Create summary report
    summary = {
        "gene_symbol": gene_symbol,
        "workflow_timestamp": datetime.now().isoformat(),
        "statistics": {
            "pmids_discovered": len(pmids),
            "pmids_filtered_out": len(dropped_pmids),
            "pmids_passed_filters": len(filtered_pmids),
            "papers_downloaded": num_downloaded,
            "papers_download_failed": num_download_failures,
            "papers_extracted": len(extractions),
            "papers_from_fulltext": len(extractions) - abstract_extraction_count,
            "papers_from_abstract_only": abstract_extraction_count,
            "papers_extraction_failed": len(extraction_failures),
            "total_variants_found": total_variants,
            "variants_with_penetrance_data": penetrance_summary["total_variants"],
            "total_carriers_observed": total_carriers,
            "total_affected_carriers": total_affected,
            "success_rate": f"{len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%"
        },
        "output_locations": {
            "pmid_list": str(combined_pmids_file),
            "pmid_status": str(pmid_status_dir),
            "full_text_papers": str(harvest_dir),
            "extractions": str(extraction_dir),
            "penetrance_summary": str(penetrance_summary_file),
            "sqlite_database": str(db_path),
            "workflow_log": str(log_file)
        },
        "pmid_status": {
            "filtered_out_file": str(filtered_out_file),
            "filtered_out_count": len(dropped_pmids),
            "extraction_failures_file": str(extraction_failures_file),
            "extraction_failures_count": len(extraction_failures),
            "download_failures_file": str(paywalled_log) if paywalled_log.exists() else None,
            "download_failures_count": num_download_failures,
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
    logger.info(f"  - Filtered out: {len(dropped_pmids)}")
    logger.info(f"  - Passed filters: {len(filtered_pmids)}")
    logger.info(f"Papers downloaded: {num_downloaded}")
    logger.info(f"  - Download failures: {num_download_failures}")
    logger.info(f"Papers with extractions: {len(extractions)}")
    logger.info(f"  - From full-text: {len(extractions) - abstract_extraction_count}")
    logger.info(f"  - From abstract only: {abstract_extraction_count}")
    logger.info(f"  - Extraction failures: {len(extraction_failures)}")
    logger.info(f"Total variants found: {total_variants}")
    logger.info(f"Variants with penetrance data: {penetrance_summary['total_variants']}")
    logger.info(f"Total carriers observed: {total_carriers}")
    logger.info(f"Total affected carriers: {total_affected}")
    logger.info(f"Success rate: {len(extractions) / len(pmids) * 100:.1f}%" if pmids else "0%")
    logger.info(f"💾 Database migrated: {migration_stats['successful']}/{migration_stats['total_files']} extractions")
    logger.info(f"\n📋 PMID status files: {pmid_status_dir}")
    logger.info(f"\nAll outputs saved to: {output_path}")
    logger.info(f"Summary report: {summary_file}")
    logger.info(f"Penetrance summary: {penetrance_summary_file}")
    logger.info(f"SQLite database: {db_path}")
    logger.info(f"Workflow log: {log_file}")
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
    parser.add_argument("--max-pmids", type=int, default=100,
                       help="Maximum PMIDs to fetch (default: 100)")
    parser.add_argument("--max-downloads", type=int, default=50,
                       help="Maximum papers to download (default: 50)")
    parser.add_argument("--tier-threshold", type=int, default=None,
                       help="If the first model finds fewer variants than this, the next model is tried (default: from .env TIER3_THRESHOLD or 1). Set to 0 to only use first model.")
    parser.add_argument("--clinical-triage", action="store_true",
                       help="Use ClinicalDataTriageFilter for Tier 2 filtering instead of InternFilter")
    parser.add_argument("--auto-synonyms", action="store_true",
                       help="Automatically discover and use gene synonyms from NCBI Gene database")
    parser.add_argument("--synonym", action="append", dest="synonyms", default=None,
                       help="Manually specify gene synonym (can be used multiple times)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose logging")

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
