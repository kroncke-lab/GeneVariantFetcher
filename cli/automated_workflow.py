#!/usr/bin/env python3
"""
Automated Workflow Entrypoint: From Gene to Variant Data

This script provides the production-ready, end-to-end automated workflow:
1. Fetch relevant PMIDs from PubMind, PubMed, and Europe PMC for a gene
2. Download full-text articles from PubMed Central
3. Extract individual-level variant and patient data
4. Save structured results to JSON and aggregate penetrance metrics

This module uses the shared step implementations from pipeline/steps.py,
matching the `gvf-run` orchestration path.
"""

import csv
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from cli.scout import run_scout

# Configure logging using centralized utility
from pipeline.checkpoint import CheckpointManager, JobCheckpoint, PipelineStep
from utils.bootstrap import has_llm_provider_key, initialize_runtime
from utils.logging_utils import get_logger, setup_logging
from utils.run_manifest import RunManifestManager

logger = get_logger(__name__)


def _apply_paper_scope_gate(
    *,
    pmids: list[str],
    abstract_records: dict[str, str],
    gene_symbol: str,
    output_path: Path,
) -> tuple[list[str], list[tuple[str, str]]]:
    """Reject obvious non-human target-gene papers before any model call.

    Explicit PMID manifests intentionally bypass the recall-oriented Tier 1 and
    Tier 2 filters, but they must never bypass deterministic clinical scope.
    Missing/unreadable abstracts pass through so full-text extraction can apply
    the same gate later.
    """

    from pipeline.filters import names_nonhuman_ortholog

    kept: list[str] = []
    dropped: list[tuple[str, str]] = []
    for pmid in pmids:
        record_path = abstract_records.get(pmid)
        record: dict = {}
        if record_path:
            try:
                record = json.loads(Path(record_path).read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                record = {}
        metadata = record.get("metadata") or {}
        title = metadata.get("title")
        abstract = record.get("abstract")
        # Prefer the title: an abstract can legitimately mention an animal
        # ortholog as background in an otherwise human study.
        scope_text = title if isinstance(title, str) and title.strip() else abstract
        if names_nonhuman_ortholog(scope_text, gene_symbol):
            dropped.append(
                (pmid, "paper studies a non-human ortholog of the target gene")
            )
        else:
            kept.append(pmid)

    status_dir = output_path / "pmid_status"
    status_dir.mkdir(parents=True, exist_ok=True)
    exclusion_path = status_dir / "paper_scope_exclusions.csv"
    with exclusion_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("pmid", "reason"))
        writer.writeheader()
        writer.writerows({"pmid": pmid, "reason": reason} for pmid, reason in dropped)
    return kept, dropped


def automated_variant_extraction_workflow(
    gene_symbol: str,
    email: str,
    output_dir: str,
    max_pmids: int = 1500,
    max_papers_to_download: int | None = None,
    tier_threshold: int = 1,
    use_clinical_triage: bool = False,
    auto_synonyms: bool = True,
    synonyms: list[str] | None = None,
    scout_first: bool = False,
    disease: str | None = None,
    include_all_clinigen_phenotypes: bool | None = None,
    pmids: list[str] | None = None,
    extraction_top_n: int | None = None,
    extraction_priority_offset: int = 0,
    extraction_triage_mode: str | None = None,
    extraction_triage_model: str | None = None,
    extraction_triage_include_defer: bool = False,
    extraction_triage_max_llm: int | None = None,
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
        disease: Optional disease term (e.g. "atrial fibrillation"). When set,
            disease clause is appended to PubMed queries and Tier-2 filter prompts
            prioritize original patient/functional data. When None (default),
            behavior is unchanged.
        include_all_clinigen_phenotypes: Include every ClinGen disease label
            for the gene in the run's disease context, even when disease is set.
        pmids: Optional explicit PMID list. When provided, skips synonym discovery,
            PMID search, AND recall-oriented Tier 1/Tier 2 relevance filtering.
            The deterministic paper-scope gate still runs before harvest so a
            non-human target-gene paper cannot enter a clinical dataset.
        extraction_top_n: Optional pre-LLM priority cap. When set, full-text and
            abstract candidates are scored for original variant/count evidence
            and only the top N are submitted for extraction first.
        extraction_priority_offset: Number of priority-ranked candidates to skip
            before applying extraction_top_n. Use 500 to process the next band
            after a top-500 pilot.
        extraction_triage_mode: Optional cheap triage pass after prioritization
            ("deterministic", "hybrid", or "llm"). When set, only triage-selected
            papers are submitted for expensive extraction.
        extraction_triage_model: Optional LLM model for triage; defaults to Tier-2.
        extraction_triage_include_defer: Extract triage "defer" papers as well.
        extraction_triage_max_llm: Optional cap on LLM triage calls.
    """
    initialize_runtime()
    setup_logging(level=logging.INFO)
    if extraction_top_n is None:
        env_top_n = os.environ.get("GVF_EXTRACTION_TOP_N", "").strip()
        if env_top_n:
            try:
                extraction_top_n = int(env_top_n)
            except ValueError:
                logger.warning(
                    "Ignoring invalid GVF_EXTRACTION_TOP_N=%r; expected integer",
                    env_top_n,
                )
                extraction_top_n = None
    if extraction_triage_mode is None:
        extraction_triage_mode = (
            os.environ.get("GVF_EXTRACTION_TRIAGE_MODE", "").strip() or None
        )
    if extraction_triage_model is None:
        extraction_triage_model = (
            os.environ.get("GVF_EXTRACTION_TRIAGE_MODEL", "").strip() or None
        )
    if extraction_triage_max_llm is None:
        env_triage_max = os.environ.get("GVF_EXTRACTION_TRIAGE_MAX_LLM", "").strip()
        if env_triage_max:
            try:
                extraction_triage_max_llm = int(env_triage_max)
            except ValueError:
                logger.warning(
                    "Ignoring invalid GVF_EXTRACTION_TRIAGE_MAX_LLM=%r; expected integer",
                    env_triage_max,
                )
    if include_all_clinigen_phenotypes is None:
        include_all_clinigen_phenotypes = os.environ.get(
            "GVF_INCLUDE_ALL_CLINGEN_PHENOTYPES", ""
        ).strip().lower() in {"1", "true", "yes", "on", "all"}

    # Track whether the caller supplied an explicit PMID list. We need this
    # later (after STEP 1 has reassigned `pmids`) to skip Tier 1/2 filtering.
    explicit_pmids_provided = bool(pmids)

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
        preprocess_papers,
    )

    # Setup output directory.  GVF_RESUME_DIR lets you point a re-run at an
    # existing timestamped run dir so cached abstracts / pmid_status / harvested
    # papers are reused (resume semantics) instead of starting fresh.
    import os as _os

    resume_dir = _os.environ.get("GVF_RESUME_DIR")
    if resume_dir:
        output_path = Path(resume_dir)
    else:
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

    # Every normal workflow gets a durable per-call LLM trace unless explicitly
    # disabled.  An explicit GVF_LLM_TRACE_DIR can place traces on encrypted or
    # high-capacity storage; otherwise they live with the run they explain.
    from utils.llm_trace import (
        TRACE_INDEX_NAME,
        configure_llm_tracing,
        resolve_trace_location,
        tracing_enabled_by_environment,
    )

    trace_enabled = tracing_enabled_by_environment()
    # GVF_LLM_TRACE_DIR is a storage BASE: this run gets its own child so
    # sequential runs never mix records (which would make every later manifest
    # rebuild raise TraceRunMismatchError). The selection is exported so
    # gvf-run's post-extraction stages continue in the SAME tree under the SAME
    # id instead of resolving a second one.
    trace_location = resolve_trace_location(
        run_manifest.run_id, default_root=output_path / "llm_traces"
    )
    trace_root = trace_location.root
    trace_run_id = trace_location.run_id
    configure_llm_tracing(trace_root, run_id=trace_run_id, enabled=trace_enabled)
    if trace_enabled:
        os.environ["GVF_LLM_TRACE_DIR"] = str(trace_root)
        os.environ["GVF_LLM_TRACE_RUN_ID"] = trace_run_id
        logger.info("LLM traces -> %s (run %s)", trace_root, trace_run_id)
    run_manifest.update_output_locations(
        llm_traces=str(trace_root),
        llm_trace_index=str(trace_root / TRACE_INDEX_NAME),
    )
    run_manifest.set_config(
        llm_tracing_enabled=trace_enabled, llm_trace_run_id=trace_run_id
    )

    from gene_literature.disease_context import build_gene_disease_context

    disease_context = build_gene_disease_context(
        gene_symbol,
        disease,
        include_all_clinigen_phenotypes=include_all_clinigen_phenotypes,
    )
    disease_context_path = output_path / "gene_disease_context.json"
    disease_context.save(disease_context_path)
    effective_disease = disease_context.prompt_disease or disease
    disease_terms = disease_context.disease_terms
    if disease_terms:
        logger.info(
            "Disease context for %s: %s",
            gene_symbol,
            "; ".join(disease_terms),
        )
    if disease_context.clinigen_curations:
        logger.info(
            "ClinGen gene-validity curations for %s: %d total, %d selected for this run",
            gene_symbol,
            len(disease_context.clinigen_curations),
            len(disease_context.selected_clinigen_curations),
        )
    for warning in disease_context.warnings:
        run_manifest.add_warning(warning)
    run_manifest.update_output_locations(gene_disease_context=str(disease_context_path))
    run_manifest.set_config(
        max_pmids=max_pmids,
        max_papers_to_download=max_papers_to_download,
        tier_threshold=tier_threshold,
        use_clinical_triage=use_clinical_triage,
        auto_synonyms=auto_synonyms,
        synonyms=synonyms or [],
        scout_first=scout_first,
        disease=effective_disease,
        requested_disease=disease,
        include_all_clinigen_phenotypes=include_all_clinigen_phenotypes,
        disease_terms=disease_terms,
        extraction_top_n=extraction_top_n,
        extraction_priority_offset=extraction_priority_offset,
        extraction_triage_mode=extraction_triage_mode,
        extraction_triage_model=extraction_triage_model,
        extraction_triage_include_defer=extraction_triage_include_defer,
        extraction_triage_max_llm=extraction_triage_max_llm,
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

    if pmids:
        logger.info(
            "\n⏭️  STEP 0: Skipping synonym discovery — explicit PMID list provided."
        )
        auto_synonyms = False

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

    if pmids:
        logger.info(
            "\n📚 STEP 1: Skipping discovery — using %d explicit PMIDs.",
            len(pmids),
        )
        # Persist the list to the same combined-output path that fetch_pmids would
        # have produced, so downstream tools that look at <gene>_pmids.txt still work.
        combined_file = output_path / f"{gene_symbol}_pmids.txt"
        combined_file.write_text("\n".join(pmids) + "\n", encoding="utf-8")
    else:
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
            disease=disease,
            disease_terms=disease_terms,
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

    scope_pmids, scope_dropped = _apply_paper_scope_gate(
        pmids=pmids,
        abstract_records=abstract_records,
        gene_symbol=gene_symbol,
        output_path=output_path,
    )
    if scope_dropped:
        logger.warning(
            "Paper-scope gate rejected %d PMID(s) before harvest: %s",
            len(scope_dropped),
            ", ".join(pmid for pmid, _ in scope_dropped),
        )

    if explicit_pmids_provided:
        logger.info(
            "\n⏭️  STEP 1.6: Skipping Tier 1/Tier 2 filtering — explicit PMID "
            "list provided. %d scope-eligible PMIDs go to harvest + extraction.",
            len(scope_pmids),
        )
        filtered_pmids = list(scope_pmids)
        dropped_pmids: list = list(scope_dropped)
    else:
        logger.info("\n🧹 STEP 1.6: Filtering papers by relevance before download...")

        filter_result = filter_papers(
            pmids=scope_pmids,
            abstract_records=abstract_records,
            gene_symbol=gene_symbol,
            output_path=output_path,
            enable_tier1=settings.enable_tier1,
            enable_tier2=settings.enable_tier2,
            use_clinical_triage=use_clinical_triage,
            tier1_min_keywords=settings.tier1_min_keywords,
            tier2_confidence_threshold=settings.tier2_confidence_threshold,
            disease=effective_disease,
        )

        filtered_pmids = filter_result.data.get("filtered_pmids", [])
        dropped_pmids = scope_dropped + filter_result.data.get("dropped_pmids", [])

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
    # STEP 2.2: Deterministic Preprocessing (Non-destructive)
    # =========================================================================
    logger.info("\n🧹 STEP 2.2: Preprocessing full-text into *_CLEANED.md files...")

    preprocess_result = preprocess_papers(
        harvest_dir=harvest_dir,
        gene_symbol=gene_symbol,
    )
    pre_stats = preprocess_result.stats
    if not pre_stats.get("skipped"):
        logger.info(
            "✓ Preprocessed %s papers (%s cleaned files)",
            pre_stats.get("processed", 0),
            pre_stats.get("cleaned_files_written", 0),
        )
        if pre_stats.get("token_savings_pct", 0) > 0:
            logger.info(
                "  ↳ token reduction estimate: %s%%",
                pre_stats.get("token_savings_pct"),
            )

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
                input_path=Path(harvest_dir),
                output_dir=Path(scout_output_dir),
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
    if extraction_top_n and extraction_top_n > 0:
        logger.info(
            "Priority gate enabled: extracting %s ranked papers after offset %s",
            extraction_top_n,
            extraction_priority_offset,
        )
    if extraction_triage_mode:
        logger.info(
            "Cheap extraction triage enabled: mode=%s model=%s include_defer=%s",
            extraction_triage_mode,
            extraction_triage_model or "(tier-2 default)",
            extraction_triage_include_defer,
        )

    extraction_dir = output_path / "extractions"
    priority_report_dir = output_path / "extraction_priority"
    triage_report_dir = output_path / "extraction_triage"

    # If scout_first was used, pass the scout manifest for better context.
    # max_workers is intentionally omitted so extract_variants picks the
    # provider-aware default (Anthropic 10, Azure 3) from Settings.
    extract_result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol=gene_symbol,
        disease=effective_disease,
        abstract_records=abstract_records,
        abstract_only_pmids=abstract_only_pmids,
        candidate_pmids=filtered_pmids,
        tier_threshold=tier_threshold,
        priority_top_n=extraction_top_n,
        priority_offset=extraction_priority_offset,
        priority_report_dir=priority_report_dir,
        priority_disease_terms=disease_terms,
        triage_mode=extraction_triage_mode,
        triage_model=extraction_triage_model,
        triage_report_dir=triage_report_dir,
        triage_include_defer=extraction_triage_include_defer,
        triage_max_llm_candidates=extraction_triage_max_llm,
    )

    extractions = extract_result.data.get("extractions", [])
    extraction_failures = extract_result.data.get("failures", [])

    logger.info(f"✓ Extracted data from {len(extractions)} papers")
    if extraction_failures:
        logger.warning(f"✗ {len(extraction_failures)} extraction failures")

    workflow_stats["papers_extracted"] = len(extractions)
    workflow_stats["extraction_failures"] = len(extraction_failures)
    workflow_stats["extraction_initial_failures"] = extract_result.stats.get(
        "initial_failures", 0
    )
    workflow_stats["extraction_retry_attempted"] = extract_result.stats.get(
        "extraction_retry_attempted", 0
    )
    workflow_stats["extraction_retry_succeeded"] = extract_result.stats.get(
        "extraction_retry_succeeded", 0
    )
    workflow_stats["extraction_retry_failed"] = extract_result.stats.get(
        "extraction_retry_failed", 0
    )
    workflow_stats["extraction_retry_skipped"] = extract_result.stats.get(
        "extraction_retry_skipped", 0
    )
    workflow_stats["total_variants_found"] = extract_result.stats.get(
        "total_variants", 0
    )
    if extract_result.stats.get("priority_top_n"):
        workflow_stats["extraction_priority_candidates"] = extract_result.stats.get(
            "priority_candidates", 0
        )
        workflow_stats["extraction_priority_selected"] = extract_result.stats.get(
            "priority_selected", 0
        )
    if extract_result.stats.get("triage_total"):
        workflow_stats["extraction_triage_total"] = extract_result.stats.get(
            "triage_total", 0
        )
        workflow_stats["extraction_triage_selected"] = extract_result.stats.get(
            "triage_selected_for_extraction", 0
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
        papers_extraction_initial_failures=extract_result.stats.get(
            "initial_failures", 0
        ),
        extraction_retry_attempted=extract_result.stats.get(
            "extraction_retry_attempted", 0
        ),
        extraction_retry_succeeded=extract_result.stats.get(
            "extraction_retry_succeeded", 0
        ),
        extraction_retry_failed=extract_result.stats.get("extraction_retry_failed", 0),
        extraction_retry_skipped=extract_result.stats.get(
            "extraction_retry_skipped", 0
        ),
        dense_table_overflow_records=extract_result.stats.get(
            "dense_table_overflow_records", 0
        ),
        scanner_cap_trips=extract_result.stats.get("scanner_cap_trips", 0),
        table_merge_cap_trips=extract_result.stats.get("table_merge_cap_trips", 0),
        table_candidates_omitted_after_dedupe=extract_result.stats.get(
            "table_candidates_omitted_after_dedupe", 0
        ),
        missing_supplement_ref_pmids=extract_result.stats.get(
            "missing_supplement_ref_pmids", 0
        ),
        total_variants_found=extract_result.stats.get("total_variants", 0),
        papers_from_fulltext=len(extractions) - abstract_extraction_count,
        papers_from_abstract_only=abstract_extraction_count,
        extraction_priority_candidates=extract_result.stats.get(
            "priority_candidates", 0
        ),
        extraction_priority_selected=extract_result.stats.get("priority_selected", 0),
        extraction_priority_fulltext_selected=extract_result.stats.get(
            "priority_fulltext_selected", 0
        ),
        extraction_priority_abstract_selected=extract_result.stats.get(
            "priority_abstract_selected", 0
        ),
        extraction_triage_total=extract_result.stats.get("triage_total", 0),
        extraction_triage_extract_now=extract_result.stats.get("triage_extract_now", 0),
        extraction_triage_defer=extract_result.stats.get("triage_defer", 0),
        extraction_triage_skip=extract_result.stats.get("triage_skip", 0),
        extraction_triage_selected=extract_result.stats.get(
            "triage_selected_for_extraction", 0
        ),
    )
    run_manifest.update_output_locations(extractions_dir=str(extraction_dir))
    if extract_result.data.get("priority_report_dir"):
        run_manifest.update_output_locations(
            extraction_priority_dir=extract_result.data["priority_report_dir"]
        )
    if extract_result.data.get("triage_report_dir"):
        run_manifest.update_output_locations(
            extraction_triage_dir=extract_result.data["triage_report_dir"]
        )
    if extract_result.data.get("retry_report_path"):
        run_manifest.update_output_locations(
            extraction_retry_summary=extract_result.data["retry_report_path"]
        )
    if extract_result.data.get("dense_table_overflow_report"):
        run_manifest.update_output_locations(
            dense_table_overflow=extract_result.data["dense_table_overflow_report"]
        )

    # =========================================================================
    # STEP 3.5: Source-Completeness Report & Zero-Variant QA
    # =========================================================================
    logger.info("\n🔎 STEP 3.5: Generating source-completeness report...")

    # --- Build completeness metrics from the FINALIZED extraction records ---
    # The denominator is the list of papers extraction actually produced a
    # record for, not what the download step intended to fetch. Reading step
    # intent is what made a resume run report 0/50 full text while its own
    # per-paper records said 50/50 (see pipeline/source_ledger.py).
    from pipeline.source_ledger import build_source_ledger

    ledger = build_source_ledger(
        [ext.extracted_data for ext in extractions if ext.extracted_data],
        requested_pmids=[str(p) for p in filtered_pmids],
        search_dirs=[Path(harvest_dir)] if harvest_dir else [],
    )
    completeness = {
        "total_pmids_discovered": workflow_stats.get("pmids_discovered", 0),
        "pmids_filtered_out": workflow_stats.get("pmids_filtered_out", 0),
        "pmids_passed_filters": workflow_stats.get("pmids_passed_filters", 0),
    }
    completeness.update(ledger.as_completeness_dict())

    # Acquisition-step intent is retained beside the ledger, clearly labelled,
    # because "the fetcher failed but the corpus cache served it" is useful
    # operational signal. It is no longer the denominator.
    completeness["acquisition_step"] = {
        "downloaded_pmids": sorted(str(p) for p in downloaded_pmids),
        "download_failed_pmids": sorted(str(p) for p in abstract_only_pmids),
    }

    # Count stub/empty files (FULL_CONTEXT < 500 bytes)
    stub_pmids = []
    if harvest_dir and Path(harvest_dir).exists():
        for fc_file in Path(harvest_dir).glob("*_FULL_CONTEXT.md"):
            if fc_file.stat().st_size < 500:
                stub_pmids.append(fc_file.name.replace("_FULL_CONTEXT.md", ""))
    completeness["stub_files_count"] = len(stub_pmids)
    completeness["stub_pmids"] = sorted(stub_pmids)

    supp_dirs = list(Path(harvest_dir).glob("*_supplements")) if harvest_dir else []
    completeness["papers_with_supplements"] = len(supp_dirs)

    zero_variant_pmids = completeness["zero_variant_pmids"]

    # Save completeness report
    completeness_file = output_path / "source_completeness.json"
    with open(completeness_file, "w") as f:
        json.dump(completeness, f, indent=2)

    logger.info(f"✓ Source completeness report: {completeness_file}")
    logger.info(f"  Full-text coverage: {completeness['fulltext_coverage_pct']}%")
    logger.info(f"  Abstract-only: {completeness['papers_abstract_only']}")
    logger.info(f"  Source unverified: {completeness['papers_source_unverified']}")
    logger.info(f"  Stub/empty files: {completeness['stub_files_count']}")
    logger.info(f"  Zero-variant papers: {completeness['zero_variant_papers']}")
    logger.info(
        f"  Single-carrier-only papers: {completeness['single_carrier_papers']}"
    )
    if completeness["source_class_discrepancies"]:
        logger.warning(
            "⚠ %d paper(s) where the declared source class disagrees with the "
            "recorded source bytes; see source_class_discrepancies",
            len(completeness["source_class_discrepancies"]),
        )
    if completeness["missing_extraction_pmids"]:
        logger.warning(
            "⚠ %d requested paper(s) produced no extraction record",
            len(completeness["missing_extraction_pmids"]),
        )

    if zero_variant_pmids:
        logger.warning(
            f"⚠ {len(zero_variant_pmids)} papers passed relevance filters but "
            f"yielded ZERO variants — potential extraction failures"
        )

    run_manifest.update_output_locations(source_completeness=str(completeness_file))

    # --- QA re-extraction for zero-variant papers ---
    # Papers that passed Tier 2 relevance but yielded 0 variants are
    # suspicious. Re-extract them with the QA model (typically Opus) to
    # check if variants were missed.
    qa_model = os.environ.get("GVF_QA_MODEL", "").strip()
    qa_max = int(os.environ.get("GVF_QA_MAX_PAPERS", "20"))

    if zero_variant_pmids and qa_model:
        qa_targets = zero_variant_pmids[:qa_max]
        logger.info(
            f"\n🔍 STEP 3.6: QA re-extraction for {len(qa_targets)} zero-variant "
            f"papers with {qa_model}..."
        )

        from pipeline.extraction import ExpertExtractor
        from utils.models import Paper as PaperModel

        qa_extractor = ExpertExtractor(
            models=[qa_model],
            fulltext_dir=str(harvest_dir),
            tier_threshold=0,  # accept any result
        )

        qa_recovered = 0
        for zv_pmid in qa_targets:
            # Load the paper text
            fc_path = Path(harvest_dir) / f"{zv_pmid}_FULL_CONTEXT.md"
            if not fc_path.exists():
                continue
            text = fc_path.read_text(encoding="utf-8", errors="replace")
            if len(text) < 200:
                continue

            paper_obj = PaperModel(
                pmid=zv_pmid,
                full_text=text,
                gene_symbol=gene_symbol,
            )

            try:
                qa_result = qa_extractor.extract(paper_obj)
                if qa_result and qa_result.extracted_data:
                    qa_variants = qa_result.extracted_data.get("variants", [])
                    if qa_variants:
                        qa_recovered += 1
                        # Save QA audit output outside the migrated top-level
                        # extraction glob, then replace the canonical JSON.
                        qa_dir = extraction_dir / "qa_reextractions"
                        qa_dir.mkdir(exist_ok=True)
                        qa_output = qa_dir / f"{gene_symbol}_PMID_{zv_pmid}_qa.json"
                        with open(qa_output, "w") as qf:
                            json.dump(qa_result.extracted_data, qf, indent=2)
                        # Also overwrite the original extraction file
                        orig_output = (
                            extraction_dir / f"{gene_symbol}_PMID_{zv_pmid}.json"
                        )
                        with open(orig_output, "w") as of:
                            json.dump(qa_result.extracted_data, of, indent=2)
                        logger.info(
                            f"  ✓ QA recovered {len(qa_variants)} variants from PMID {zv_pmid}"
                        )
            except Exception as e:
                logger.warning(f"  QA extraction failed for {zv_pmid}: {e}")

        completeness["qa_reextracted"] = len(qa_targets)
        completeness["qa_recovered"] = qa_recovered
        # Re-save completeness with QA data
        with open(completeness_file, "w") as f:
            json.dump(completeness, f, indent=2)
        logger.info(
            f"✓ QA re-extraction: recovered variants from {qa_recovered}/{len(qa_targets)} papers"
        )
    elif zero_variant_pmids:
        logger.info(
            "ℹ Set GVF_QA_MODEL=anthropic/claude-opus-4-7 to re-extract "
            f"{len(zero_variant_pmids)} zero-variant papers with a stronger model"
        )

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

    if not migrate_result.success:
        # A complete migration crash returns empty stats, so migration_failed is
        # 0 and the partial-failure branch below would miss it entirely.
        run_manifest.add_warning(
            f"SQLite migration failed completely: "
            f"{migrate_result.error or 'unknown error'}"
        )
    elif workflow_stats["migration_failed"]:
        # Record on the run manifest so a partial migration is visible in the
        # durable run status, not just the transient log.
        run_manifest.add_warning(
            f"SQLite migration: {workflow_stats['migration_failed']}/"
            f"{migrate_result.stats.get('total_files', 0)} extraction files "
            f"failed to migrate (rolled back; no rows written for them)"
        )

    # Update manifest with database location
    run_manifest.update_output_locations(sqlite_database=str(db_path))

    # =========================================================================
    # STEP 6: Compile Results and Statistics
    # =========================================================================
    logger.info("\n📊 STEP 6: Compiling results and statistics...")

    # Derive headline counts from the SQLite DB (canonical sink) instead of the
    # in-memory aggregator output. The aggregator pass and the migrator pass
    # are independent — when the aggregator returns 0 (e.g. variant-key
    # mismatch) but the migrator wrote 478 penetrance rows, we want the
    # workflow_summary to reflect what landed in the DB.
    db_penetrance_variants = 0
    db_total_carriers = 0
    db_total_affected = 0
    db_unique_pmids = 0
    db_variant_paper_links = 0
    try:
        import sqlite3 as _sqlite3

        with _sqlite3.connect(str(db_path)) as _dbconn:
            _c = _dbconn.cursor()
            _c.execute("SELECT COUNT(DISTINCT variant_id) FROM penetrance_data")
            db_penetrance_variants = _c.fetchone()[0] or 0
            _c.execute(
                "SELECT COALESCE(SUM(total_carriers_observed), 0), "
                "COALESCE(SUM(affected_count), 0) FROM penetrance_data"
            )
            row = _c.fetchone()
            db_total_carriers, db_total_affected = int(row[0]), int(row[1])
            _c.execute("SELECT COUNT(DISTINCT pmid) FROM variant_papers")
            db_unique_pmids = _c.fetchone()[0] or 0
            _c.execute("SELECT COUNT(*) FROM variant_papers")
            db_variant_paper_links = _c.fetchone()[0] or 0
    except Exception as _db_err:
        logger.warning(
            f"Could not derive headline counts from DB ({db_path}): {_db_err}. "
            f"Falling back to in-memory aggregator totals."
        )
        db_penetrance_variants = sum(
            1
            for v in penetrance_summary.get("variants", [])
            if v.get("aggregated_penetrance", {}).get("total_carriers")
        )
        db_total_carriers = sum(
            v.get("aggregated_penetrance", {}).get("total_carriers", 0) or 0
            for v in penetrance_summary.get("variants", [])
        )
        db_total_affected = sum(
            v.get("aggregated_penetrance", {}).get("affected", 0) or 0
            for v in penetrance_summary.get("variants", [])
        )

    total_carriers = db_total_carriers
    total_affected = db_total_affected

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
            "variants_with_penetrance_data": db_penetrance_variants,
            "total_carriers_observed": total_carriers,
            "total_affected_carriers": total_affected,
            "unique_pmids_with_variants": db_unique_pmids,
            "variant_paper_links": db_variant_paper_links,
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
            "llm_traces": str(trace_root) if trace_enabled else None,
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
    trace_manifest_path = None
    trace_report_path = None
    if trace_enabled:
        try:
            from utils.llm_trace import TRACE_MANIFEST_NAME, build_trace_manifest
            from utils.llm_trace_html import (
                TRACE_REPORT_NAME,
                build_trace_html_report,
            )

            candidate_manifest_path = trace_root / TRACE_MANIFEST_NAME
            build_trace_manifest(
                trace_root,
                output_path=candidate_manifest_path,
                run_id=trace_run_id,
            )
            trace_manifest_path = candidate_manifest_path
            candidate_report_path = output_path / TRACE_REPORT_NAME
            report_data = build_trace_html_report(
                trace_root,
                output_path=candidate_report_path,
                run_dir=output_path,
                title=f"{gene_symbol} · LLM curation trace review",
                run_id=trace_run_id,
            )
            trace_report_path = candidate_report_path
            summary["files"]["llm_trace_report"] = str(trace_report_path)
            summary["llm_trace"] = {
                "integrity_level": report_data["integrity"]["level"],
                "omissions": len(report_data.get("omissions") or []),
                "sharded": bool(report_data.get("sharded")),
                "missing_decision_links": report_data["coverage"][
                    "missing_decision_links"
                ],
            }
            with open(summary_file, "w") as f:
                json.dump(summary, f, indent=2)
        except Exception as trace_error:  # tracing remains best-effort
            logger.warning("Could not finalize LLM trace artifacts: %s", trace_error)
    run_manifest.update_statistics(
        **workflow_stats,
        success_rate=success_rate,
        total_carriers_observed=total_carriers,
        total_affected_carriers=total_affected,
    )
    run_manifest.update_output_locations(
        workflow_summary=str(summary_file),
        workflow_log=str(log_file),
        llm_trace_manifest=(
            str(trace_manifest_path) if trace_manifest_path is not None else None
        ),
        llm_trace_report=(
            str(trace_report_path) if trace_report_path is not None else None
        ),
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
    initialize_runtime()
    setup_logging(level=logging.INFO)

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
        default=1500,
        help="Maximum PMIDs to fetch (default: 1500)",
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

    # Require at least one LLM provider key. The pipeline routes via
    # TIER*_MODEL strings and Settings resolves the current provider.
    if not has_llm_provider_key():
        logger.error("⚠️  ERROR: No LLM provider API key found in environment!")
        logger.error(
            "Set OPENAI_API_KEY, AZURE_AI_API_KEY, or ANTHROPIC_API_KEY in your .env file."
        )
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
