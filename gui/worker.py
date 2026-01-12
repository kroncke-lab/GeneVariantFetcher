"""
Background worker for running pipeline jobs with checkpointing.

This module wraps the automated workflow to:
1. Save checkpoints after each major step
2. Support resume from any checkpoint
3. Stream progress to connected clients
"""

import json
import logging
import os
import csv
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from dotenv import load_dotenv

load_dotenv()

from gui.checkpoint import (
    CheckpointManager,
    JobCheckpoint,
    PipelineStep,
)

logger = logging.getLogger(__name__)


class ProgressCallback:
    """Callback interface for progress updates."""

    def on_step_start(self, step: PipelineStep, message: str):
        """Called when a step starts."""
        pass

    def on_step_complete(self, step: PipelineStep, message: str, stats: Dict[str, Any]):
        """Called when a step completes."""
        pass

    def on_log(self, message: str):
        """Called for log messages."""
        pass

    def on_error(self, step: PipelineStep, error: str):
        """Called on error."""
        pass


class LoggingProgressCallback(ProgressCallback):
    """Progress callback that logs to standard logger."""

    def on_step_start(self, step: PipelineStep, message: str):
        logger.info(f"[{step.value}] Starting: {message}")

    def on_step_complete(self, step: PipelineStep, message: str, stats: Dict[str, Any]):
        logger.info(f"[{step.value}] Complete: {message}")
        if stats:
            logger.info(f"  Stats: {stats}")

    def on_log(self, message: str):
        logger.info(message)

    def on_error(self, step: PipelineStep, error: str):
        logger.error(f"[{step.value}] Error: {error}")


class PipelineWorker:
    """
    Runs the variant extraction pipeline with checkpoint support.

    This worker breaks the pipeline into discrete, resumable steps.
    Each step is saved to disk, allowing recovery after interruption.
    """

    def __init__(
        self,
        checkpoint_manager: Optional[CheckpointManager] = None,
        progress_callback: Optional[ProgressCallback] = None,
    ):
        self.checkpoint_manager = checkpoint_manager or CheckpointManager()
        self.progress = progress_callback or LoggingProgressCallback()
        self._stop_requested = False

    def request_stop(self):
        """Request graceful stop after current step."""
        self._stop_requested = True

    def _should_stop(self) -> bool:
        """Check if stop was requested."""
        return self._stop_requested

    def _save_checkpoint(self, checkpoint: JobCheckpoint):
        """Save checkpoint and notify progress."""
        self.checkpoint_manager.save(checkpoint)

    def _log(self, checkpoint: JobCheckpoint, message: str):
        """Log message and add to checkpoint."""
        checkpoint.add_log(message)
        self.progress.on_log(message)

    def run(self, checkpoint: JobCheckpoint, resume: bool = False) -> JobCheckpoint:
        """
        Run the pipeline from current checkpoint state.

        Args:
            checkpoint: Job checkpoint with configuration
            resume: If True, skip already-completed steps

        Returns:
            Updated checkpoint with final state
        """
        self._stop_requested = False

        # Determine starting step
        if resume and checkpoint.current_step != PipelineStep.PENDING:
            start_step = checkpoint.current_step
            self._log(checkpoint, f"Resuming from step: {start_step.value}")
        else:
            start_step = PipelineStep.PENDING
            checkpoint.started_at = datetime.now().isoformat()

        # Create output directory structure
        output_path = checkpoint.output_path
        output_path.mkdir(parents=True, exist_ok=True)

        try:
            # Execute steps in order, skipping completed ones if resuming
            steps = [
                (PipelineStep.DISCOVERING_SYNONYMS, self._step_discover_synonyms),
                (PipelineStep.FETCHING_PMIDS, self._step_fetch_pmids),
                (PipelineStep.FETCHING_ABSTRACTS, self._step_fetch_abstracts),
                (PipelineStep.FILTERING_PAPERS, self._step_filter_papers),
                (PipelineStep.DOWNLOADING_FULLTEXT, self._step_download_fulltext),
                (PipelineStep.EXTRACTING_VARIANTS, self._step_extract_variants),
                (PipelineStep.AGGREGATING_DATA, self._step_aggregate_data),
                (PipelineStep.MIGRATING_DATABASE, self._step_migrate_database),
            ]

            started = not resume
            for step, step_func in steps:
                # Skip until we reach the resume point
                if not started:
                    if step == start_step:
                        started = True
                    else:
                        continue

                # Check for stop request
                if self._should_stop():
                    self._log(checkpoint, "Stop requested, saving checkpoint...")
                    self._save_checkpoint(checkpoint)
                    return checkpoint

                # Run step
                self.progress.on_step_start(
                    step, PipelineStep.get_display_name(step)
                )
                checkpoint.update_step(step)
                self._save_checkpoint(checkpoint)

                try:
                    stats = step_func(checkpoint)
                    self.progress.on_step_complete(
                        step, f"Completed {step.value}", stats or {}
                    )
                except Exception as e:
                    checkpoint.mark_failed(str(e), step)
                    self._save_checkpoint(checkpoint)
                    self.progress.on_error(step, str(e))
                    raise

            # Mark completed
            checkpoint.mark_completed()
            self._save_checkpoint(checkpoint)
            self._log(checkpoint, "Pipeline completed successfully!")

            return checkpoint

        except Exception as e:
            logger.exception(f"Pipeline failed: {e}")
            if checkpoint.current_step != PipelineStep.FAILED:
                checkpoint.mark_failed(str(e))
                self._save_checkpoint(checkpoint)
            raise

    def _step_discover_synonyms(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 0: Discover gene synonyms (optional)."""
        if not checkpoint.auto_synonyms:
            self._log(checkpoint, "Synonym discovery disabled, skipping...")
            return {"skipped": True}

        self._log(checkpoint, f"Discovering synonyms for {checkpoint.gene_symbol}...")

        from gene_literature.synonym_finder import SynonymFinder, automatic_synonym_selection

        synonym_finder = SynonymFinder(
            email=checkpoint.email,
            api_key=os.getenv("NCBI_API_KEY"),
        )

        found_synonyms = synonym_finder.find_gene_synonyms(
            checkpoint.gene_symbol,
            include_other_designations=False,
        )

        auto_selected = automatic_synonym_selection(
            checkpoint.gene_symbol,
            found_synonyms,
            include_official=True,
            include_aliases=True,
            include_other_designations=False,
            only_relevant=False,
        )

        # Merge with manually provided synonyms
        existing = set(s.lower() for s in checkpoint.synonyms)
        for syn in auto_selected:
            if syn.lower() not in existing and syn.lower() != checkpoint.gene_symbol.lower():
                checkpoint.synonyms.append(syn)
                existing.add(syn.lower())

        self._log(checkpoint, f"Found {len(auto_selected)} synonyms: {', '.join(checkpoint.synonyms)}")

        return {"synonyms_found": len(auto_selected), "total_synonyms": len(checkpoint.synonyms)}

    def _step_fetch_pmids(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 1: Fetch PMIDs from literature sources."""
        self._log(checkpoint, "Fetching PMIDs from PubMind, PubMed, and Europe PMC...")

        from gene_literature.discovery import discover_pmids_for_gene
        from config.settings import get_settings

        output_path = checkpoint.output_path
        settings = get_settings()

        # Override settings with checkpoint config
        settings.use_pubmind = checkpoint.use_pubmind
        settings.use_pubmed = checkpoint.use_pubmed
        settings.use_europepmc = checkpoint.use_europepmc

        pubmind_file = output_path / f"{checkpoint.gene_symbol}_pmids_pubmind.txt"
        pubmed_file = output_path / f"{checkpoint.gene_symbol}_pmids_pubmed.txt"
        combined_file = output_path / f"{checkpoint.gene_symbol}_pmids.txt"

        pmid_discovery = discover_pmids_for_gene(
            gene_symbol=checkpoint.gene_symbol,
            email=checkpoint.email,
            max_results=checkpoint.max_pmids,
            pubmind_output=pubmind_file,
            pubmed_output=pubmed_file,
            combined_output=combined_file,
            api_key=os.getenv("NCBI_API_KEY"),
            settings=settings,
            synonyms=checkpoint.synonyms if checkpoint.synonyms else None,
        )

        checkpoint.discovered_pmids = list(pmid_discovery.combined_pmids)
        self._log(checkpoint, f"Discovered {len(checkpoint.discovered_pmids)} unique PMIDs")

        return {
            "pubmind_count": len(pmid_discovery.pubmind_pmids),
            "pubmed_count": len(pmid_discovery.pubmed_pmids),
            "europepmc_count": len(pmid_discovery.europepmc_pmids),
            "total_unique": len(checkpoint.discovered_pmids),
        }

    def _step_fetch_abstracts(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 1.5: Fetch abstracts and metadata."""
        self._log(checkpoint, "Fetching abstracts and metadata...")

        from harvesting.abstracts import fetch_and_save_abstracts

        output_path = checkpoint.output_path
        abstract_dir = output_path / "abstract_json"

        abstract_records = fetch_and_save_abstracts(
            pmids=checkpoint.discovered_pmids,
            output_dir=str(abstract_dir),
            email=checkpoint.email,
        )

        # Store the paths for filtering step
        checkpoint.step_progress["abstract_records"] = {
            pmid: str(path) for pmid, path in abstract_records.items()
        }

        self._log(checkpoint, f"Fetched abstracts for {len(abstract_records)} PMIDs")

        return {"abstracts_fetched": len(abstract_records)}

    def _step_filter_papers(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 1.6: Filter papers by relevance."""
        self._log(checkpoint, "Filtering papers by relevance...")

        from pipeline.filters import KeywordFilter, InternFilter, ClinicalDataTriageFilter
        from utils.models import Paper, FilterDecision
        from config.settings import get_settings

        settings = get_settings()
        output_path = checkpoint.output_path

        # Get abstract records from previous step
        abstract_records = checkpoint.step_progress.get("abstract_records", {})

        keyword_filter = KeywordFilter(min_keyword_matches=settings.tier1_min_keywords)
        tier2_filter = (
            ClinicalDataTriageFilter()
            if checkpoint.use_clinical_triage
            else InternFilter(confidence_threshold=checkpoint.tier2_confidence_threshold)
        )

        filtered_pmids = []
        dropped_pmids = []

        for pmid in checkpoint.discovered_pmids:
            record_path = abstract_records.get(pmid)
            if not record_path or not Path(record_path).exists():
                dropped_pmids.append((pmid, "Missing abstract"))
                continue

            try:
                with open(record_path, "r", encoding="utf-8") as f:
                    record = json.load(f)
            except Exception as e:
                dropped_pmids.append((pmid, f"Read error: {e}"))
                continue

            metadata = record.get("metadata", {})
            paper = Paper(
                pmid=pmid,
                title=metadata.get("title"),
                abstract=record.get("abstract"),
                authors=metadata.get("authors"),
                journal=metadata.get("journal"),
                publication_date=metadata.get("year"),
                gene_symbol=checkpoint.gene_symbol,
                source="PubMed",
            )

            # Tier 1 filtering
            if checkpoint.enable_tier1:
                tier1_result = keyword_filter.filter(paper)
                if tier1_result.decision is not FilterDecision.PASS:
                    dropped_pmids.append((pmid, tier1_result.reason))
                    continue

            # Tier 2 filtering
            if checkpoint.enable_tier2:
                if checkpoint.use_clinical_triage:
                    triage_result = tier2_filter.triage_paper(paper, checkpoint.gene_symbol)
                    if triage_result.get("decision") != "KEEP":
                        dropped_pmids.append((pmid, triage_result.get("reason", "Failed triage")))
                        continue
                else:
                    tier2_result = tier2_filter.filter(paper)
                    if tier2_result.decision is not FilterDecision.PASS:
                        dropped_pmids.append((pmid, tier2_result.reason))
                        continue

            filtered_pmids.append(pmid)

        checkpoint.filtered_pmids = filtered_pmids

        # Save dropped PMIDs for reference
        pmid_status_dir = output_path / "pmid_status"
        pmid_status_dir.mkdir(parents=True, exist_ok=True)
        filtered_out_file = pmid_status_dir / "filtered_out.csv"
        with open(filtered_out_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["PMID", "Reason"])
            writer.writerows(dropped_pmids)

        self._log(checkpoint, f"Passed filters: {len(filtered_pmids)}, Dropped: {len(dropped_pmids)}")

        return {
            "passed_filters": len(filtered_pmids),
            "dropped": len(dropped_pmids),
        }

    def _step_download_fulltext(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 2: Download full-text papers."""
        self._log(checkpoint, "Downloading full-text papers from PubMed Central...")

        from harvesting import PMCHarvester

        output_path = checkpoint.output_path
        harvest_dir = output_path / "pmc_fulltext"

        harvester = PMCHarvester(
            output_dir=str(harvest_dir),
            gene_symbol=checkpoint.gene_symbol,
        )

        # Limit downloads
        pmids_to_download = checkpoint.filtered_pmids[:checkpoint.max_papers_to_download]
        self._log(checkpoint, f"Downloading up to {len(pmids_to_download)} papers...")

        harvester.harvest(pmids_to_download, delay=2.0)

        # Check results
        success_log = harvest_dir / "successful_downloads.csv"
        num_downloaded = 0
        if success_log.exists():
            import pandas as pd
            successful_downloads = pd.read_csv(success_log)
            num_downloaded = len(successful_downloads)
            checkpoint.downloaded_pmids = successful_downloads["PMID"].astype(str).tolist()

        self._log(checkpoint, f"Successfully downloaded {num_downloaded} papers")

        return {"downloaded": num_downloaded, "attempted": len(pmids_to_download)}

    def _step_extract_variants(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 3: Extract variant data using LLM."""
        self._log(checkpoint, "Extracting variant data using AI...")

        from utils.models import Paper
        from pipeline.extraction import ExpertExtractor
        from utils.pmid_utils import extract_pmid_from_filename

        output_path = checkpoint.output_path
        harvest_dir = output_path / "pmc_fulltext"
        extraction_dir = output_path / "extractions"
        extraction_dir.mkdir(exist_ok=True)

        # Find markdown files (prefer DATA_ZONES.md)
        data_zones = {f.name.replace("_DATA_ZONES.md", ""): f for f in harvest_dir.glob("*_DATA_ZONES.md")}
        full_context = {f.name.replace("_FULL_CONTEXT.md", ""): f for f in harvest_dir.glob("*_FULL_CONTEXT.md")}

        all_pmids = set(data_zones.keys()) | set(full_context.keys())
        markdown_files = []
        for pmid in all_pmids:
            if pmid in data_zones:
                markdown_files.append(data_zones[pmid])
            elif pmid in full_context:
                markdown_files.append(full_context[pmid])

        self._log(checkpoint, f"Processing {len(markdown_files)} papers for extraction...")

        extractor = ExpertExtractor(
            tier_threshold=checkpoint.tier_threshold,
            fulltext_dir=str(harvest_dir),
        )

        extractions = []
        extraction_failures = []

        def process_paper(md_file):
            pmid_match = extract_pmid_from_filename(md_file)
            if not pmid_match:
                return (None, (md_file.name, "Could not extract PMID"))

            with open(md_file, "r", encoding="utf-8") as f:
                full_text = f.read()

            paper = Paper(
                pmid=pmid_match,
                title=f"Paper {pmid_match}",
                full_text=full_text,
                gene_symbol=checkpoint.gene_symbol,
            )

            try:
                result = extractor.extract(paper)
                if result.success:
                    output_file = extraction_dir / f"{checkpoint.gene_symbol}_PMID_{pmid_match}.json"
                    with open(output_file, "w") as f:
                        json.dump(result.extracted_data, f, indent=2)
                    return (result, None)
                else:
                    return (None, (pmid_match, result.error))
            except Exception as e:
                return (None, (pmid_match, str(e)))

        # Process in parallel
        max_workers = min(8, len(markdown_files)) if markdown_files else 1
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_paper, f): f for f in markdown_files}
            completed = 0
            for future in as_completed(futures):
                completed += 1
                result, failure = future.result()
                if result:
                    extractions.append(result)
                if failure:
                    extraction_failures.append(failure)
                if completed % 5 == 0:
                    self._log(checkpoint, f"Processed {completed}/{len(markdown_files)} papers...")

        checkpoint.extracted_pmids = [str(e.extracted_data.get("pmid", "")) for e in extractions if e.extracted_data]

        # Save failures
        pmid_status_dir = output_path / "pmid_status"
        pmid_status_dir.mkdir(parents=True, exist_ok=True)
        failures_file = pmid_status_dir / "extraction_failures.csv"
        with open(failures_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["PMID", "Error"])
            writer.writerows(extraction_failures)

        # Count variants
        total_variants = sum(
            e.extracted_data.get("extraction_metadata", {}).get("total_variants_found", 0)
            for e in extractions if e.extracted_data
        )

        self._log(checkpoint, f"Extracted {total_variants} variants from {len(extractions)} papers")

        return {
            "papers_extracted": len(extractions),
            "failures": len(extraction_failures),
            "total_variants": total_variants,
        }

    def _step_aggregate_data(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 4: Aggregate penetrance data."""
        self._log(checkpoint, "Aggregating penetrance data...")

        from pipeline.aggregation import aggregate_penetrance

        output_path = checkpoint.output_path
        extraction_dir = output_path / "extractions"
        summary_file = output_path / f"{checkpoint.gene_symbol}_penetrance_summary.json"

        penetrance_summary = aggregate_penetrance(
            extraction_dir=extraction_dir,
            gene_symbol=checkpoint.gene_symbol,
            output_file=summary_file,
        )

        total_variants = penetrance_summary.get("total_variants", 0)
        self._log(checkpoint, f"Aggregated data for {total_variants} variants")

        return {"variants_aggregated": total_variants}

    def _step_migrate_database(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 5: Migrate to SQLite database."""
        self._log(checkpoint, "Creating SQLite database...")

        from harvesting.migrate_to_sqlite import create_database_schema, migrate_extraction_directory

        output_path = checkpoint.output_path
        extraction_dir = output_path / "extractions"
        db_path = output_path / f"{checkpoint.gene_symbol}.db"

        conn = create_database_schema(str(db_path))
        migration_stats = migrate_extraction_directory(conn, extraction_dir)
        conn.close()

        self._log(
            checkpoint,
            f"Migrated {migration_stats['successful']}/{migration_stats['total_files']} extractions to database",
        )

        # Write final summary
        self._write_workflow_summary(checkpoint, migration_stats)

        return migration_stats

    def _write_workflow_summary(self, checkpoint: JobCheckpoint, migration_stats: Dict[str, Any]):
        """Write final workflow summary JSON."""
        output_path = checkpoint.output_path

        summary = {
            "job_id": checkpoint.job_id,
            "gene_symbol": checkpoint.gene_symbol,
            "workflow_timestamp": checkpoint.created_at,
            "completed_at": checkpoint.completed_at,
            "statistics": {
                "pmids_discovered": len(checkpoint.discovered_pmids),
                "pmids_passed_filters": len(checkpoint.filtered_pmids),
                "papers_downloaded": len(checkpoint.downloaded_pmids),
                "papers_extracted": len(checkpoint.extracted_pmids),
            },
            "database_migration": migration_stats,
            "output_locations": {
                "base_dir": str(output_path),
                "sqlite_database": str(output_path / f"{checkpoint.gene_symbol}.db"),
                "extractions": str(output_path / "extractions"),
                "pmid_status": str(output_path / "pmid_status"),
            },
        }

        summary_file = output_path / f"{checkpoint.gene_symbol}_workflow_summary.json"
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2)


def create_job(
    gene_symbol: str,
    email: str,
    output_dir: str,
    **kwargs,
) -> JobCheckpoint:
    """
    Create a new pipeline job.

    Args:
        gene_symbol: Gene to analyze
        email: Email for NCBI
        output_dir: Base output directory
        **kwargs: Additional configuration options

    Returns:
        New JobCheckpoint ready to run
    """
    import uuid

    job_id = f"{gene_symbol.lower()}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"

    checkpoint = JobCheckpoint(
        job_id=job_id,
        gene_symbol=gene_symbol.upper(),
        email=email,
        output_dir=output_dir,
        max_pmids=kwargs.get("max_pmids", 100),
        max_papers_to_download=kwargs.get("max_papers_to_download", 50),
        tier_threshold=kwargs.get("tier_threshold", 1),
        use_clinical_triage=kwargs.get("use_clinical_triage", False),
        auto_synonyms=kwargs.get("auto_synonyms", False),
        synonyms=kwargs.get("synonyms", []),
        enable_tier1=kwargs.get("enable_tier1", True),
        enable_tier2=kwargs.get("enable_tier2", True),
        use_pubmind=kwargs.get("use_pubmind", True),
        use_pubmed=kwargs.get("use_pubmed", True),
        use_europepmc=kwargs.get("use_europepmc", False),
        tier2_confidence_threshold=kwargs.get("tier2_confidence_threshold", 0.5),
        scout_enabled=kwargs.get("scout_enabled", True),
    )

    return checkpoint
