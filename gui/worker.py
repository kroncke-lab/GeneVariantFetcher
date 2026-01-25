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
from typing import Any, Callable, Dict, List, Optional, Set

from dotenv import load_dotenv

load_dotenv()

from gui.checkpoint import (
    CheckpointManager,
    JobCheckpoint,
    PipelineStep,
)
from utils.logging_utils import setup_logging, get_logger

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

        # Log job configuration for debugging
        self._log(checkpoint, f"=== PIPELINE START ===")
        self._log(checkpoint, f"Job ID: {checkpoint.job_id}")
        self._log(checkpoint, f"Gene: {checkpoint.gene_symbol}")
        self._log(checkpoint, f"is_folder_job: {checkpoint.is_folder_job}")
        self._log(checkpoint, f"is_resume_job: {checkpoint.is_resume_job}")
        self._log(checkpoint, f"resume_stage: {checkpoint.resume_stage}")
        self._log(checkpoint, f"folder_path: {checkpoint.folder_path}")
        self._log(checkpoint, f"run_scout_on_folder: {checkpoint.run_scout_on_folder}")
        self._log(checkpoint, f"scout_enabled: {checkpoint.scout_enabled}")
        self._log(
            checkpoint, f"skip_already_extracted: {checkpoint.skip_already_extracted}"
        )
        self._log(checkpoint, f"resume parameter: {resume}")
        self._log(checkpoint, f"current_step: {checkpoint.current_step.value}")

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

        # Set up file logging in the output directory for troubleshooting
        log_file = output_path / "workflow.log"
        setup_logging(level=logging.INFO, log_file=log_file)
        self._log(checkpoint, f"Logging to file: {log_file}")

        try:
            # Determine which steps to run based on job type
            if checkpoint.is_resume_job and checkpoint.resume_stage == "downloading":
                # Resume from downloading stage - includes download, scout, extract, etc.
                steps = [
                    (PipelineStep.DOWNLOADING_FULLTEXT, self._step_download_fulltext),
                ]
                if checkpoint.run_scout_on_folder:
                    steps.append((PipelineStep.SCOUTING_DATA, self._step_scout_data))
                steps.extend(
                    [
                        (PipelineStep.EXTRACTING_VARIANTS, self._step_extract_variants),
                        (PipelineStep.AGGREGATING_DATA, self._step_aggregate_data),
                        (PipelineStep.MIGRATING_DATABASE, self._step_migrate_database),
                    ]
                )
            elif checkpoint.is_folder_job:
                # Folder jobs skip discovery/download, go straight to scouting/extraction
                steps = []
                if checkpoint.run_scout_on_folder:
                    steps.append((PipelineStep.SCOUTING_DATA, self._step_scout_data))
                steps.extend(
                    [
                        (PipelineStep.EXTRACTING_VARIANTS, self._step_extract_variants),
                        (PipelineStep.AGGREGATING_DATA, self._step_aggregate_data),
                        (PipelineStep.MIGRATING_DATABASE, self._step_migrate_database),
                    ]
                )
            else:
                # Full pipeline
                steps = [
                    (PipelineStep.DISCOVERING_SYNONYMS, self._step_discover_synonyms),
                    (PipelineStep.FETCHING_PMIDS, self._step_fetch_pmids),
                    (PipelineStep.FETCHING_ABSTRACTS, self._step_fetch_abstracts),
                    (PipelineStep.FILTERING_PAPERS, self._step_filter_papers),
                    (PipelineStep.DOWNLOADING_FULLTEXT, self._step_download_fulltext),
                    (PipelineStep.SCOUTING_DATA, self._step_scout_data),
                    (PipelineStep.EXTRACTING_VARIANTS, self._step_extract_variants),
                    (PipelineStep.AGGREGATING_DATA, self._step_aggregate_data),
                    (PipelineStep.MIGRATING_DATABASE, self._step_migrate_database),
                ]

            # Log which steps will be executed
            step_names = [s[0].value for s in steps]
            self._log(checkpoint, f"=== STEPS TO EXECUTE: {step_names} ===")
            self._log(
                checkpoint, f"start_step: {start_step.value}, resume param: {resume}"
            )

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
                self._log(checkpoint, f">>> Entering step: {step.value}")
                self.progress.on_step_start(step, PipelineStep.get_display_name(step))
                checkpoint.update_step(step)
                self._save_checkpoint(checkpoint)

                try:
                    stats = step_func(checkpoint)
                    # Store step stats on checkpoint for later access
                    if stats:
                        checkpoint.step_progress[step.value] = stats
                        self._save_checkpoint(checkpoint)
                    self.progress.on_step_complete(
                        step, f"Completed {step.value}", stats or {}
                    )
                    self._log(
                        checkpoint,
                        f">>> Step completed: {step.value} with stats={stats}",
                    )
                except Exception as e:
                    import traceback

                    error_tb = traceback.format_exc()
                    self._log(checkpoint, f">>> Step FAILED: {step.value}")
                    self._log(checkpoint, f">>> Error: {e}")
                    self._log(checkpoint, f">>> Traceback: {error_tb}")
                    checkpoint.mark_failed(str(e), step)
                    self._save_checkpoint(checkpoint)
                    self.progress.on_error(step, str(e))
                    raise

            # Mark completed
            checkpoint.mark_completed()
            self._save_checkpoint(checkpoint)
            self._log(checkpoint, "=== PIPELINE COMPLETED SUCCESSFULLY ===")
            self._log(
                checkpoint,
                f"Final stats: extracted_pmids={len(checkpoint.extracted_pmids)}",
            )

            return checkpoint

        except Exception as e:
            import traceback

            error_traceback = traceback.format_exc()
            logger.exception(f"Pipeline failed: {e}")
            self._log(checkpoint, f"=== PIPELINE FAILED ===")
            self._log(checkpoint, f"Error: {e}")
            self._log(checkpoint, f"Traceback: {error_traceback}")
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

        from gene_literature.synonym_finder import (
            SynonymFinder,
            automatic_synonym_selection,
        )

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
            if (
                syn.lower() not in existing
                and syn.lower() != checkpoint.gene_symbol.lower()
            ):
                checkpoint.synonyms.append(syn)
                existing.add(syn.lower())

        self._log(
            checkpoint,
            f"Found {len(auto_selected)} synonyms: {', '.join(checkpoint.synonyms)}",
        )

        return {
            "synonyms_found": len(auto_selected),
            "total_synonyms": len(checkpoint.synonyms),
        }

    def _step_fetch_pmids(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 1: Fetch PMIDs from literature sources."""

        # If specific PMIDs provided, use those directly (for testing)
        if checkpoint.specific_pmids:
            self._log(
                checkpoint,
                f"Using {len(checkpoint.specific_pmids)} specific PMID(s) for testing...",
            )
            checkpoint.discovered_pmids = list(checkpoint.specific_pmids)

            # Save to combined file for consistency
            output_path = checkpoint.output_path
            combined_file = output_path / f"{checkpoint.gene_symbol}_pmids.txt"
            with open(combined_file, "w") as f:
                for pmid in checkpoint.specific_pmids:
                    f.write(f"{pmid}\n")

            self._log(
                checkpoint,
                f"Using specific PMIDs: {', '.join(checkpoint.specific_pmids)}",
            )

            return {
                "pubmind_count": 0,
                "pubmed_count": 0,
                "europepmc_count": 0,
                "specific_pmids": len(checkpoint.specific_pmids),
                "total_unique": len(checkpoint.discovered_pmids),
            }

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
        self._log(
            checkpoint, f"Discovered {len(checkpoint.discovered_pmids)} unique PMIDs"
        )

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

        from pipeline.filters import (
            KeywordFilter,
            InternFilter,
            ClinicalDataTriageFilter,
        )
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
            else InternFilter(
                confidence_threshold=checkpoint.tier2_confidence_threshold
            )
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
                    triage_result = tier2_filter.triage_paper(
                        paper, checkpoint.gene_symbol
                    )
                    if triage_result.get("decision") != "KEEP":
                        dropped_pmids.append(
                            (pmid, triage_result.get("reason", "Failed triage"))
                        )
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

        self._log(
            checkpoint,
            f"Passed filters: {len(filtered_pmids)}, Dropped: {len(dropped_pmids)}",
        )

        return {
            "passed_filters": len(filtered_pmids),
            "dropped": len(dropped_pmids),
        }

    def _step_download_fulltext(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 2: Download full-text papers."""
        self._log(checkpoint, "Downloading full-text papers from PubMed Central...")

        from harvesting import PMCHarvester

        # Determine output directory based on job type
        if checkpoint.is_resume_job and checkpoint.folder_path:
            # Resume job: use the existing folder structure
            folder_path = Path(checkpoint.folder_path)
            harvest_dir = folder_path / "pmc_fulltext"
        else:
            output_path = checkpoint.output_path
            harvest_dir = output_path / "pmc_fulltext"

        harvester = PMCHarvester(
            output_dir=str(harvest_dir),
            gene_symbol=checkpoint.gene_symbol,
        )

        # Determine which PMIDs to download
        if checkpoint.is_resume_job and checkpoint.pmids_to_download:
            # Resume job: download the specific missing PMIDs
            pmids_to_download = checkpoint.pmids_to_download[
                : checkpoint.max_papers_to_download
            ]
            self._log(
                checkpoint,
                f"Resuming download for {len(pmids_to_download)} missing papers...",
            )
        else:
            # Normal job: use filtered PMIDs
            pmids_to_download = checkpoint.filtered_pmids[
                : checkpoint.max_papers_to_download
            ]
            self._log(
                checkpoint, f"Downloading up to {len(pmids_to_download)} papers..."
            )

        harvester.harvest(pmids_to_download, delay=2.0)

        # Check results - for resume jobs, merge with existing downloads
        success_log = harvest_dir / "successful_downloads.csv"
        num_downloaded = 0
        newly_downloaded = []

        if success_log.exists():
            import pandas as pd

            successful_downloads = pd.read_csv(success_log)
            num_downloaded = len(successful_downloads)
            newly_downloaded = successful_downloads["PMID"].astype(str).tolist()

        # Merge with existing downloaded PMIDs (for resume jobs)
        if checkpoint.is_resume_job:
            existing_pmids = set(checkpoint.downloaded_pmids)
            new_pmids = set(newly_downloaded)
            checkpoint.downloaded_pmids = list(existing_pmids | new_pmids)
            self._log(
                checkpoint,
                f"Resume complete. Total downloaded: {len(checkpoint.downloaded_pmids)} papers ({len(new_pmids - existing_pmids)} new)",
            )
        else:
            checkpoint.downloaded_pmids = newly_downloaded
            self._log(checkpoint, f"Successfully downloaded {num_downloaded} papers")

        return {"downloaded": num_downloaded, "attempted": len(pmids_to_download)}

    def _step_scout_data(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 2.5: Run Data Scout to create DATA_ZONES files."""
        self._log(checkpoint, "=== ZONE FINDING START ===")

        if not checkpoint.scout_enabled:
            self._log(
                checkpoint, "Data Scout disabled (scout_enabled=False), skipping..."
            )
            self._log(checkpoint, "=== ZONE FINDING END (skipped) ===")
            return {"skipped": True}

        self._log(checkpoint, "Running Data Scout to identify high-value data zones...")

        from pipeline.data_scout import GeneticDataScout

        # Determine the fulltext directory
        if checkpoint.is_folder_job and checkpoint.folder_path:
            folder_path = Path(checkpoint.folder_path)
            fulltext_dir = folder_path / "pmc_fulltext"
            if not fulltext_dir.exists():
                fulltext_dir = folder_path
        else:
            fulltext_dir = checkpoint.output_path / "pmc_fulltext"

        if not fulltext_dir.exists():
            self._log(
                checkpoint,
                f"No fulltext directory found at {fulltext_dir}, skipping scout...",
            )
            self._log(checkpoint, "=== ZONE FINDING END (no fulltext dir) ===")
            return {"skipped": True, "reason": "no_fulltext_dir"}

        # Find FULL_CONTEXT files that don't have corresponding DATA_ZONES
        full_context_files = list(fulltext_dir.glob("*_FULL_CONTEXT.md"))
        existing_zones = {
            f.name.replace("_DATA_ZONES.md", "")
            for f in fulltext_dir.glob("*_DATA_ZONES.md")
        }

        files_to_scout = [
            f
            for f in full_context_files
            if f.name.replace("_FULL_CONTEXT.md", "") not in existing_zones
        ]

        if not files_to_scout:
            self._log(
                checkpoint,
                f"All {len(existing_zones)} papers already scouted, skipping...",
            )
            self._log(checkpoint, "=== ZONE FINDING END (already complete) ===")
            return {"skipped": True, "already_scouted": len(existing_zones)}

        self._log(checkpoint, f"Scouting {len(files_to_scout)} papers...")

        scout = GeneticDataScout(
            gene_symbol=checkpoint.gene_symbol,
            min_relevance_score=checkpoint.scout_min_relevance,
            max_zones=30,
        )

        scouted = 0
        errors = 0

        for md_file in files_to_scout:
            try:
                with open(md_file, "r", encoding="utf-8") as f:
                    full_text = f.read()

                # Extract PMID from filename
                pmid_match = md_file.name.replace("_FULL_CONTEXT.md", "")

                # Run scout
                report = scout.scan(full_text, pmid=pmid_match)

                # Write DATA_ZONES.md
                zones_file = md_file.parent / f"{pmid_match}_DATA_ZONES.md"
                zones_markdown = scout.format_markdown(report, full_text)
                with open(zones_file, "w", encoding="utf-8") as f:
                    f.write(zones_markdown)

                # Write zones JSON for reference
                zones_json_file = md_file.parent / f"{pmid_match}_DATA_ZONES.json"
                with open(zones_json_file, "w", encoding="utf-8") as f:
                    f.write(scout.to_full_json(report))

                scouted += 1
                if scouted % 5 == 0:
                    self._log(
                        checkpoint, f"Scouted {scouted}/{len(files_to_scout)} papers..."
                    )

            except Exception as e:
                errors += 1
                self._log(checkpoint, f"Error scouting {md_file.name}: {e}")

        self._log(checkpoint, f"Scouted {scouted} papers ({errors} errors)")
        self._log(checkpoint, "=== ZONE FINDING END (success) ===")

        return {
            "scouted": scouted,
            "errors": errors,
            "already_scouted": len(existing_zones),
        }

    def _step_extract_variants(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 3: Extract variant data using LLM."""
        self._log(checkpoint, "=== VARIANT EXTRACTION START ===")
        self._log(checkpoint, "Extracting variant data using AI...")

        from utils.models import Paper
        from pipeline.extraction import ExpertExtractor
        from utils.pmid_utils import extract_pmid_from_filename

        # Determine directories based on job type
        if checkpoint.is_folder_job and checkpoint.folder_path:
            base_path = Path(checkpoint.folder_path)
            harvest_dir = base_path / "pmc_fulltext"
            if not harvest_dir.exists():
                harvest_dir = base_path
            extraction_dir = base_path / "extractions"
            self._log(checkpoint, f"Using folder job paths:")
        else:
            base_path = checkpoint.output_path
            harvest_dir = base_path / "pmc_fulltext"
            extraction_dir = base_path / "extractions"
            self._log(checkpoint, f"Using standard paths:")

        self._log(checkpoint, f"  base_path: {base_path}")
        self._log(checkpoint, f"  harvest_dir: {harvest_dir}")
        self._log(checkpoint, f"  extraction_dir: {extraction_dir}")

        extraction_dir.mkdir(exist_ok=True)

        # Find markdown files (prefer DATA_ZONES.md)
        data_zones = {
            f.name.replace("_DATA_ZONES.md", ""): f
            for f in harvest_dir.glob("*_DATA_ZONES.md")
        }
        full_context = {
            f.name.replace("_FULL_CONTEXT.md", ""): f
            for f in harvest_dir.glob("*_FULL_CONTEXT.md")
        }

        self._log(
            checkpoint,
            f"Found {len(data_zones)} DATA_ZONES files, {len(full_context)} FULL_CONTEXT files",
        )

        all_pmids = set(data_zones.keys()) | set(full_context.keys())
        self._log(checkpoint, f"Total unique PMIDs to consider: {len(all_pmids)}")

        # Skip already extracted PMIDs if configured (for folder jobs)
        if checkpoint.skip_already_extracted and checkpoint.extracted_pmids:
            already_extracted = set(checkpoint.extracted_pmids)
            skipped_count = len(all_pmids & already_extracted)
            all_pmids = all_pmids - already_extracted
            self._log(
                checkpoint,
                f"skip_already_extracted=True: {skipped_count} in extracted_pmids list, {len(all_pmids)} remaining",
            )
        else:
            self._log(
                checkpoint,
                f"skip_already_extracted={checkpoint.skip_already_extracted}, extracted_pmids count={len(checkpoint.extracted_pmids)}",
            )

        markdown_files = []
        for pmid in all_pmids:
            if pmid in data_zones:
                markdown_files.append(data_zones[pmid])
            elif pmid in full_context:
                markdown_files.append(full_context[pmid])

        self._log(
            checkpoint, f"Processing {len(markdown_files)} papers for extraction..."
        )

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
                    output_file = (
                        extraction_dir
                        / f"{checkpoint.gene_symbol}_PMID_{pmid_match}.json"
                    )
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
                    self._log(
                        checkpoint,
                        f"Processed {completed}/{len(markdown_files)} papers...",
                    )

        checkpoint.extracted_pmids = [
            str(e.extracted_data.get("pmid", ""))
            for e in extractions
            if e.extracted_data
        ]

        # Save failures
        pmid_status_dir = base_path / "pmid_status"
        pmid_status_dir.mkdir(parents=True, exist_ok=True)
        failures_file = pmid_status_dir / "extraction_failures.csv"
        with open(failures_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["PMID", "Error"])
            writer.writerows(extraction_failures)

        # Count variants
        total_variants = sum(
            e.extracted_data.get("extraction_metadata", {}).get(
                "total_variants_found", 0
            )
            for e in extractions
            if e.extracted_data
        )

        self._log(
            checkpoint,
            f"Extracted {total_variants} variants from {len(extractions)} papers",
        )
        self._log(checkpoint, f"Extraction failures: {len(extraction_failures)}")
        self._log(checkpoint, f"Failures file: {failures_file}")
        self._log(checkpoint, "=== VARIANT EXTRACTION END ===")

        return {
            "papers_extracted": len(extractions),
            "failures": len(extraction_failures),
            "total_variants": total_variants,
        }

    def _step_aggregate_data(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 4: Aggregate penetrance data."""
        self._log(checkpoint, "=== AGGREGATION START ===")
        self._log(checkpoint, "Aggregating penetrance data...")

        from pipeline.aggregation import aggregate_penetrance

        # Determine directories based on job type
        if checkpoint.is_folder_job and checkpoint.folder_path:
            base_path = Path(checkpoint.folder_path)
        else:
            base_path = checkpoint.output_path

        extraction_dir = base_path / "extractions"
        summary_file = base_path / f"{checkpoint.gene_symbol}_penetrance_summary.json"

        self._log(checkpoint, f"  base_path: {base_path}")
        self._log(checkpoint, f"  extraction_dir: {extraction_dir}")
        self._log(checkpoint, f"  summary_file: {summary_file}")

        penetrance_summary = aggregate_penetrance(
            extraction_dir=extraction_dir,
            gene_symbol=checkpoint.gene_symbol,
            output_file=summary_file,
        )

        total_variants = penetrance_summary.get("total_variants", 0)
        self._log(checkpoint, f"Aggregated data for {total_variants} variants")
        self._log(checkpoint, "=== AGGREGATION END ===")

        return {"variants_aggregated": total_variants}

    def _step_migrate_database(self, checkpoint: JobCheckpoint) -> Dict[str, Any]:
        """Step 5: Migrate to SQLite database."""
        self._log(checkpoint, "=== DATABASE MIGRATION START ===")
        self._log(checkpoint, "Creating SQLite database...")

        from harvesting.migrate_to_sqlite import (
            create_database_schema,
            migrate_extraction_directory,
        )

        # Determine directories based on job type
        if checkpoint.is_folder_job and checkpoint.folder_path:
            base_path = Path(checkpoint.folder_path)
        else:
            base_path = checkpoint.output_path

        extraction_dir = base_path / "extractions"
        db_path = base_path / f"{checkpoint.gene_symbol}.db"

        self._log(checkpoint, f"  base_path: {base_path}")
        self._log(checkpoint, f"  extraction_dir: {extraction_dir}")
        self._log(checkpoint, f"  DB PATH: {db_path}")

        # Check if extraction directory exists and has files
        if not extraction_dir.exists():
            self._log(
                checkpoint, f"WARNING: extraction_dir does not exist: {extraction_dir}"
            )
        else:
            extraction_files = list(extraction_dir.glob("*.json"))
            self._log(
                checkpoint, f"Found {len(extraction_files)} extraction JSON files"
            )

        conn = create_database_schema(str(db_path))
        migration_stats = migrate_extraction_directory(conn, extraction_dir)
        conn.close()

        self._log(
            checkpoint,
            f"Migrated {migration_stats['successful']}/{migration_stats['total_files']} extractions to database",
        )
        self._log(checkpoint, f"Migration stats: {migration_stats}")
        self._log(checkpoint, f"Database written to: {db_path}")
        self._log(checkpoint, "=== DATABASE MIGRATION END ===")

        # Write final summary
        self._write_workflow_summary(checkpoint, migration_stats)

        return migration_stats

    def _write_workflow_summary(
        self, checkpoint: JobCheckpoint, migration_stats: Dict[str, Any]
    ):
        """Write final workflow summary JSON."""
        # Determine base path based on job type
        if checkpoint.is_folder_job and checkpoint.folder_path:
            output_path = Path(checkpoint.folder_path)
        else:
            output_path = checkpoint.output_path

        # Get total variants from extraction step stats or penetrance summary
        total_variants = 0

        # First try to get from step_progress (extraction step)
        extraction_stats = checkpoint.step_progress.get(
            PipelineStep.EXTRACTING_VARIANTS.value, {}
        )
        if "total_variants" in extraction_stats:
            total_variants = extraction_stats["total_variants"]

        # If not found, try aggregation step
        if total_variants == 0:
            aggregation_stats = checkpoint.step_progress.get(
                PipelineStep.AGGREGATING_DATA.value, {}
            )
            if "variants_aggregated" in aggregation_stats:
                total_variants = aggregation_stats["variants_aggregated"]

        # As fallback, read from penetrance summary file
        if total_variants == 0:
            penetrance_file = (
                output_path / f"{checkpoint.gene_symbol}_penetrance_summary.json"
            )
            if penetrance_file.exists():
                try:
                    with open(penetrance_file, "r") as f:
                        penetrance_data = json.load(f)
                        total_variants = penetrance_data.get("total_variants", 0)
                except Exception:
                    pass

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
                "total_variants_found": total_variants,
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
        specific_pmids=kwargs.get("specific_pmids", []),
        enable_tier1=kwargs.get("enable_tier1", True),
        enable_tier2=kwargs.get("enable_tier2", True),
        use_pubmind=kwargs.get("use_pubmind", True),
        use_pubmed=kwargs.get("use_pubmed", True),
        use_europepmc=kwargs.get("use_europepmc", False),
        tier2_confidence_threshold=kwargs.get("tier2_confidence_threshold", 0.5),
        scout_enabled=kwargs.get("scout_enabled", True),
    )

    return checkpoint


def create_folder_job(
    folder_path: str,
    gene_symbol: str,
    email: str,
    run_scout: bool = False,
    skip_already_extracted: bool = True,
    tier_threshold: int = 1,
    scout_enabled: bool = True,
    scout_min_relevance: float = 0.3,
    full_context_pmids: List[str] = None,
    data_zones_pmids: List[str] = None,
    extraction_pmids: List[str] = None,
) -> JobCheckpoint:
    """
    Create a job that processes an existing folder with downloaded papers.

    This job skips discovery, filtering, and download steps, going directly
    to extraction (or scouting first if enabled).

    Args:
        folder_path: Path to folder containing pmc_fulltext/ with markdown files
        gene_symbol: Gene symbol for extraction
        email: Email for NCBI (may be needed for some operations)
        run_scout: Whether to run Data Scout on unscouted files first
        skip_already_extracted: Skip PMIDs that already have extractions
        tier_threshold: Model cascade threshold for extraction
        scout_enabled: Enable DATA_ZONES preference during extraction
        scout_min_relevance: Minimum relevance score for scout zones
        full_context_pmids: List of PMIDs with FULL_CONTEXT files
        data_zones_pmids: List of PMIDs with DATA_ZONES files
        extraction_pmids: List of PMIDs already extracted (to skip)

    Returns:
        JobCheckpoint configured for folder-based extraction
    """
    import uuid

    job_id = f"{gene_symbol.lower()}_folder_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"

    # Determine starting step
    if run_scout:
        starting_step = PipelineStep.SCOUTING_DATA
    else:
        starting_step = PipelineStep.EXTRACTING_VARIANTS

    checkpoint = JobCheckpoint(
        job_id=job_id,
        gene_symbol=gene_symbol.upper(),
        email=email,
        output_dir=folder_path,  # Use folder as output dir for path calculations
        tier_threshold=tier_threshold,
        scout_enabled=scout_enabled,
        scout_min_relevance=scout_min_relevance,
        # Folder job specific settings
        is_folder_job=True,
        folder_path=folder_path,
        run_scout_on_folder=run_scout,
        skip_already_extracted=skip_already_extracted,
        # Pre-populate PMID lists from analysis
        downloaded_pmids=list(
            set((full_context_pmids or []) + (data_zones_pmids or []))
        ),
        extracted_pmids=extraction_pmids or [],
        # Start at appropriate step
        current_step=starting_step,
    )

    return checkpoint


def create_resume_job(
    folder_path: str,
    gene_symbol: str,
    email: str,
    resume_stage: str,
    pmids_to_download: List[str] = None,
    max_papers_to_download: int = 500,
    run_scout: bool = True,
    skip_already_extracted: bool = True,
    tier_threshold: int = 1,
    scout_enabled: bool = True,
    scout_min_relevance: float = 0.3,
    full_context_pmids: List[str] = None,
    data_zones_pmids: List[str] = None,
    extraction_pmids: List[str] = None,
) -> JobCheckpoint:
    """
    Create a job that resumes an interrupted pipeline from a specific stage.

    This job handles two main scenarios:
    1. Resume downloading - when the pipeline was interrupted during download
    2. Resume extraction - when downloading is complete but extraction was interrupted

    Args:
        folder_path: Path to folder with partial pipeline results
        gene_symbol: Gene symbol for the pipeline
        email: Email for NCBI E-utilities
        resume_stage: Stage to resume from ('downloading' or 'extraction')
        pmids_to_download: PMIDs to download (for downloading stage)
        max_papers_to_download: Maximum papers to download
        run_scout: Whether to run Data Scout after downloading
        skip_already_extracted: Skip PMIDs that already have extractions
        tier_threshold: Model cascade threshold for extraction
        scout_enabled: Enable DATA_ZONES preference during extraction
        scout_min_relevance: Minimum relevance score for scout zones
        full_context_pmids: List of PMIDs with FULL_CONTEXT files
        data_zones_pmids: List of PMIDs with DATA_ZONES files
        extraction_pmids: List of PMIDs already extracted (to skip)

    Returns:
        JobCheckpoint configured for resuming the pipeline
    """
    import uuid

    job_id = f"{gene_symbol.lower()}_resume_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"

    # Determine starting step based on resume stage
    if resume_stage == "downloading":
        starting_step = PipelineStep.DOWNLOADING_FULLTEXT
    else:
        starting_step = (
            PipelineStep.SCOUTING_DATA
            if run_scout
            else PipelineStep.EXTRACTING_VARIANTS
        )

    checkpoint = JobCheckpoint(
        job_id=job_id,
        gene_symbol=gene_symbol.upper(),
        email=email,
        output_dir=folder_path,
        max_papers_to_download=max_papers_to_download,
        tier_threshold=tier_threshold,
        scout_enabled=scout_enabled,
        scout_min_relevance=scout_min_relevance,
        # Resume job specific settings
        is_resume_job=True,
        resume_stage=resume_stage,
        pmids_to_download=pmids_to_download or [],
        # Also set folder job flags for path handling
        is_folder_job=True,
        folder_path=folder_path,
        run_scout_on_folder=run_scout,
        skip_already_extracted=skip_already_extracted,
        # Pre-populate PMID lists from existing files
        downloaded_pmids=list(
            set((full_context_pmids or []) + (data_zones_pmids or []))
        ),
        extracted_pmids=extraction_pmids or [],
        # For resume, also pre-populate filtered_pmids to include both downloaded and to-download
        filtered_pmids=list(
            set(
                (full_context_pmids or [])
                + (data_zones_pmids or [])
                + (pmids_to_download or [])
            )
        ),
        # Start at appropriate step
        current_step=starting_step,
    )

    return checkpoint
