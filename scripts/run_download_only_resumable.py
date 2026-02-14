#!/usr/bin/env python3
"""
Resumable KCNH2 download-only pipeline.

Supports checkpointing and resumption after OOM/kill. Each stage:
1. Discover PMIDs (PubMed + PubMind) - skips if pmids.txt exists
2. Fetch abstracts (~6000) - skips already-fetched abstracts
3. Tier 2 LLM filter (Grok, ~1/sec) - resumes via filter_progress.jsonl
4. Download full-text papers - skips already-downloaded papers

Usage:
    python scripts/run_download_only_resumable.py                # Auto-detect resume
    python scripts/run_download_only_resumable.py --resume       # Force resume mode
    python scripts/run_download_only_resumable.py --fresh        # Force fresh start
    python scripts/run_download_only_resumable.py --status       # Show current state
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set

# Add project root
sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv
load_dotenv()

from gui.checkpoint import CheckpointManager, JobCheckpoint, PipelineStep
from utils.logging_utils import setup_logging

setup_logging(level=logging.INFO)
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

CONFIG = {
    "gene_symbol": "KCNH2",
    "email": "brett.kroncke@gmail.com",
    "output_dir": "/mnt/temp2/kronckbm/gvf_output",
    "max_pmids": 10000,
    "max_papers_to_download": 1000,
    "tier2_confidence_threshold": 0.5,
}


# =============================================================================
# Resume Detection & State Recovery
# =============================================================================

def find_latest_run_dir(output_dir: str, gene_symbol: str) -> Optional[Path]:
    """Find the most recent run directory for this gene that has actual content."""
    gene_dir = Path(output_dir) / gene_symbol
    if not gene_dir.exists():
        return None
    
    run_dirs = [d for d in gene_dir.iterdir() if d.is_dir() and d.name[0].isdigit()]
    if not run_dirs:
        return None
    
    # Sort by timestamp (format: YYYYMMDD_HHMMSS), newest first
    run_dirs.sort(key=lambda d: d.name, reverse=True)
    
    # Find the first directory with actual content (PMIDs file or abstracts)
    for run_dir in run_dirs:
        pmids_file = run_dir / f"{gene_symbol}_pmids.txt"
        abstract_dir = run_dir / "abstract_json"
        if pmids_file.exists() or abstract_dir.exists():
            return run_dir
    
    # Fall back to newest if none have content
    return run_dirs[0]


def scan_existing_state(run_dir: Path, gene_symbol: str) -> Dict:
    """Scan a run directory to determine current state."""
    state = {
        "run_dir": run_dir,
        "pmids_file": None,
        "discovered_pmids": [],
        "abstract_dir": None,
        "existing_abstracts": set(),
        "filter_progress_file": None,
        "filtered_pmids": [],
        "filter_completed": False,
        "fulltext_dir": None,
        "downloaded_pmids": set(),
    }
    
    # Check for PMIDs file
    pmids_file = run_dir / f"{gene_symbol}_pmids.txt"
    if pmids_file.exists():
        state["pmids_file"] = pmids_file
        state["discovered_pmids"] = [
            line.strip() for line in pmids_file.read_text().strip().split("\n")
            if line.strip()
        ]
        logger.info(f"Found {len(state['discovered_pmids'])} discovered PMIDs")
    
    # Check for abstracts
    abstract_dir = run_dir / "abstract_json"
    if abstract_dir.exists():
        state["abstract_dir"] = abstract_dir
        state["existing_abstracts"] = {
            f.stem for f in abstract_dir.glob("*.json")
        }
        logger.info(f"Found {len(state['existing_abstracts'])} existing abstracts")
    
    # Check filter progress
    pmid_status_dir = run_dir / "pmid_status"
    filter_progress = pmid_status_dir / "filter_progress.jsonl"
    if filter_progress.exists():
        state["filter_progress_file"] = filter_progress
        # Count processed and passed
        processed = 0
        passed = []
        with open(filter_progress, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                    processed += 1
                    if rec.get("final_decision", "").upper() == "PASS":
                        passed.append(rec.get("pmid"))
                except:
                    continue
        state["filter_processed_count"] = processed
        state["filtered_pmids"] = passed
        # Check if filtering is complete
        if state["discovered_pmids"]:
            state["filter_completed"] = processed >= len(state["discovered_pmids"])
        logger.info(f"Filter progress: {processed} processed, {len(passed)} passed")
    
    # Check for filtered_pmids.txt (final filter output)
    filtered_pmids_file = pmid_status_dir / "filtered_pmids.txt"
    if filtered_pmids_file.exists():
        state["filtered_pmids"] = [
            line.strip() for line in filtered_pmids_file.read_text().strip().split("\n")
            if line.strip()
        ]
        state["filter_completed"] = True
        logger.info(f"Found completed filter output: {len(state['filtered_pmids'])} PMIDs")
    
    # Check for downloaded papers
    fulltext_dir = run_dir / "pmc_fulltext"
    if fulltext_dir.exists():
        state["fulltext_dir"] = fulltext_dir
        state["downloaded_pmids"] = {
            f.name.replace("_FULL_CONTEXT.md", "")
            for f in fulltext_dir.glob("*_FULL_CONTEXT.md")
        }
        logger.info(f"Found {len(state['downloaded_pmids'])} downloaded papers")
    
    return state


def determine_resume_stage(state: Dict) -> Optional[PipelineStep]:
    """Determine which stage to resume from based on existing state."""
    
    # No PMIDs discovered yet
    if not state["discovered_pmids"]:
        return PipelineStep.FETCHING_PMIDS
    
    # Have PMIDs but not all abstracts fetched
    if state["existing_abstracts"]:
        missing_abstracts = set(state["discovered_pmids"]) - state["existing_abstracts"]
        if missing_abstracts:
            logger.info(f"Missing {len(missing_abstracts)} abstracts")
            return PipelineStep.FETCHING_ABSTRACTS
    else:
        return PipelineStep.FETCHING_ABSTRACTS
    
    # Have abstracts but filtering not complete
    if not state["filter_completed"]:
        return PipelineStep.FILTERING_PAPERS
    
    # Filtering done but not all papers downloaded
    if state["filtered_pmids"]:
        downloaded = state["downloaded_pmids"]
        to_download = set(state["filtered_pmids"][:CONFIG["max_papers_to_download"]])
        missing_downloads = to_download - downloaded
        if missing_downloads:
            logger.info(f"Missing {len(missing_downloads)} paper downloads")
            return PipelineStep.DOWNLOADING_FULLTEXT
    
    # Everything done
    return PipelineStep.COMPLETED


def print_status(state: Dict):
    """Print current pipeline status."""
    print("\n" + "="*60)
    print("PIPELINE STATUS")
    print("="*60)
    
    print(f"\nRun directory: {state['run_dir']}")
    print(f"\nStage Progress:")
    print(f"  1. PMID Discovery:    {len(state['discovered_pmids']):,} PMIDs")
    print(f"  2. Abstract Fetch:    {len(state['existing_abstracts']):,} / {len(state['discovered_pmids']):,}")
    
    filter_processed = state.get("filter_processed_count", 0)
    filter_passed = len(state["filtered_pmids"])
    filter_status = "COMPLETE" if state["filter_completed"] else f"{filter_processed:,} processed"
    print(f"  3. Tier 2 Filter:     {filter_passed:,} passed ({filter_status})")
    
    downloaded = len(state["downloaded_pmids"])
    to_download = min(filter_passed, CONFIG["max_papers_to_download"])
    print(f"  4. Full-text DL:      {downloaded:,} / {to_download:,}")
    
    resume_stage = determine_resume_stage(state)
    if resume_stage == PipelineStep.COMPLETED:
        print(f"\n✓ Pipeline COMPLETE")
    else:
        print(f"\n→ Resume from: {resume_stage.value}")
    
    print("="*60 + "\n")


# =============================================================================
# Resumable Pipeline Steps
# =============================================================================

def step_fetch_pmids_resumable(checkpoint: JobCheckpoint, state: Dict) -> Dict:
    """Fetch PMIDs, skipping if already done."""
    from pipeline import steps as workflow_steps
    import os
    
    # Check if PMIDs already exist
    if state["discovered_pmids"]:
        logger.info(f"PMIDs already discovered: {len(state['discovered_pmids'])}")
        checkpoint.discovered_pmids = state["discovered_pmids"]
        return {"skipped": True, "total_unique": len(state["discovered_pmids"])}
    
    logger.info("Discovering PMIDs from PubMed + PubMind...")
    
    result = workflow_steps.fetch_pmids(
        gene_symbol=checkpoint.gene_symbol,
        email=checkpoint.email,
        output_path=checkpoint.output_path,
        max_results=checkpoint.max_pmids,
        synonyms=checkpoint.synonyms if checkpoint.synonyms else None,
        use_pubmind=True,
        use_pubmed=True,
        use_europepmc=False,
        api_key=os.getenv("NCBI_API_KEY"),
    )
    
    if result.success:
        checkpoint.discovered_pmids = result.data.get("pmids", [])
        logger.info(f"Discovered {len(checkpoint.discovered_pmids)} unique PMIDs")
    else:
        logger.error(f"PMID discovery failed: {result.error}")
        checkpoint.discovered_pmids = []
    
    return result.stats


def step_fetch_abstracts_resumable(checkpoint: JobCheckpoint, state: Dict) -> Dict:
    """Fetch abstracts, skipping already-fetched ones."""
    from harvesting.abstracts import fetch_and_save_abstracts
    from utils.pubmed_utils import batch_fetch_metadata, fetch_paper_abstract
    import json
    
    abstract_dir = checkpoint.output_path / "abstract_json"
    abstract_dir.mkdir(parents=True, exist_ok=True)
    
    # Find which abstracts we still need
    all_pmids = set(checkpoint.discovered_pmids)
    existing = {f.stem for f in abstract_dir.glob("*.json")}
    missing = all_pmids - existing
    
    logger.info(f"Abstracts: {len(existing)} existing, {len(missing)} to fetch")
    
    if not missing:
        # All abstracts exist, build the records dict
        abstract_records = {
            pmid: str(abstract_dir / f"{pmid}.json")
            for pmid in existing
        }
        checkpoint.step_progress["abstract_records"] = abstract_records
        return {"skipped": True, "abstracts_fetched": len(existing)}
    
    # Fetch missing abstracts in batches
    missing_list = list(missing)
    batch_size = 200  # Fetch metadata in batches
    fetched = 0
    
    for i in range(0, len(missing_list), batch_size):
        batch = missing_list[i:i+batch_size]
        logger.info(f"Fetching abstracts batch {i//batch_size + 1} ({len(batch)} PMIDs)...")
        
        try:
            # Batch fetch metadata
            metadata_dict = batch_fetch_metadata(batch, email=checkpoint.email)
            
            for pmid in batch:
                output_file = abstract_dir / f"{pmid}.json"
                
                # Skip if file was created by another process
                if output_file.exists():
                    fetched += 1
                    continue
                
                metadata = metadata_dict.get(pmid, {})
                
                try:
                    abstract = fetch_paper_abstract(pmid, email=checkpoint.email)
                except Exception as e:
                    logger.warning(f"Failed to fetch abstract for {pmid}: {e}")
                    abstract = None
                
                # Build author list
                raw_authors = metadata.get("AuthorList", [])
                authors = []
                for author in raw_authors:
                    if isinstance(author, dict):
                        name = author.get("Name") or " ".join(
                            part for part in [author.get("ForeName"), author.get("LastName")]
                            if part
                        )
                        if name:
                            authors.append(str(name))
                    else:
                        authors.append(str(author))
                
                # Extract year
                year = None
                for field in ("PubDate", "EPubDate", "SortPubDate"):
                    date_value = metadata.get(field)
                    if date_value:
                        year_str = str(date_value).strip().split(" ")[0]
                        if year_str and year_str[:4].isdigit():
                            year = year_str[:4]
                            break
                
                record = {
                    "metadata": {
                        "pmid": pmid,
                        "title": metadata.get("Title"),
                        "authors": authors,
                        "journal": metadata.get("FullJournalName") or metadata.get("Source"),
                        "year": year,
                    },
                    "abstract": abstract,
                }
                
                with open(output_file, "w", encoding="utf-8") as f:
                    json.dump(record, f, indent=2)
                
                fetched += 1
                
        except Exception as e:
            logger.error(f"Batch fetch error: {e}")
            continue
        
        logger.info(f"Progress: {fetched + len(existing)} / {len(all_pmids)} abstracts")
    
    # Build final records dict
    abstract_records = {
        pmid: str(abstract_dir / f"{pmid}.json")
        for pmid in all_pmids
        if (abstract_dir / f"{pmid}.json").exists()
    }
    checkpoint.step_progress["abstract_records"] = abstract_records
    
    return {"abstracts_fetched": len(abstract_records), "newly_fetched": fetched}


def step_filter_papers_resumable(checkpoint: JobCheckpoint, state: Dict) -> Dict:
    """Filter papers - automatically resumes via filter_progress.jsonl."""
    from pipeline import steps as workflow_steps
    
    abstract_records = checkpoint.step_progress.get("abstract_records", {})
    
    # If no abstract records in checkpoint, rebuild from disk
    if not abstract_records:
        abstract_dir = checkpoint.output_path / "abstract_json"
        if abstract_dir.exists():
            abstract_records = {
                f.stem: str(f) for f in abstract_dir.glob("*.json")
            }
            checkpoint.step_progress["abstract_records"] = abstract_records
    
    logger.info(f"Filtering {len(checkpoint.discovered_pmids)} PMIDs...")
    logger.info("(Resume support: filter_progress.jsonl tracks progress)")
    
    # This function already has resume support via filter_progress.jsonl
    result = workflow_steps.filter_papers(
        pmids=checkpoint.discovered_pmids,
        abstract_records=abstract_records,
        gene_symbol=checkpoint.gene_symbol,
        output_path=checkpoint.output_path,
        enable_tier1=checkpoint.enable_tier1,
        enable_tier2=checkpoint.enable_tier2,
        use_clinical_triage=checkpoint.use_clinical_triage,
        tier1_min_keywords=1,
        tier2_confidence_threshold=checkpoint.tier2_confidence_threshold,
    )
    
    if result.success:
        checkpoint.filtered_pmids = result.data.get("filtered_pmids", [])
        logger.info(f"Filter complete: {len(checkpoint.filtered_pmids)} passed")
    else:
        logger.error(f"Filtering failed: {result.error}")
    
    return result.stats


def step_download_fulltext_resumable(checkpoint: JobCheckpoint, state: Dict) -> Dict:
    """Download full-text papers, skipping already-downloaded."""
    from pipeline import steps as workflow_steps
    
    pmids_to_download = checkpoint.filtered_pmids[:checkpoint.max_papers_to_download]
    
    # Check what's already downloaded
    fulltext_dir = checkpoint.output_path / "pmc_fulltext"
    already_downloaded = set()
    if fulltext_dir.exists():
        already_downloaded = {
            f.name.replace("_FULL_CONTEXT.md", "")
            for f in fulltext_dir.glob("*_FULL_CONTEXT.md")
        }
    
    remaining = [p for p in pmids_to_download if p not in already_downloaded]
    
    logger.info(f"Downloads: {len(already_downloaded)} existing, {len(remaining)} remaining")
    
    if not remaining:
        checkpoint.downloaded_pmids = list(already_downloaded)
        return {
            "skipped": True,
            "downloaded": len(already_downloaded),
            "attempted": len(pmids_to_download),
        }
    
    # Download remaining papers
    result = workflow_steps.download_fulltext(
        pmids=remaining,
        output_path=checkpoint.output_path,
        gene_symbol=checkpoint.gene_symbol,
        max_papers=len(remaining),
        delay=2.0,
    )
    
    # Merge with already downloaded
    newly_downloaded = set(result.data.get("downloaded_pmids", []))
    all_downloaded = already_downloaded | newly_downloaded
    checkpoint.downloaded_pmids = list(all_downloaded)
    
    result.stats["previously_downloaded"] = len(already_downloaded)
    result.stats["total_downloaded"] = len(all_downloaded)
    
    logger.info(f"Download complete: {len(all_downloaded)} total papers")
    
    return result.stats


def step_scout_data_resumable(checkpoint: JobCheckpoint, state: Dict) -> Dict:
    """Run Data Scout, skipping already-scouted papers."""
    from pipeline import steps as workflow_steps
    
    fulltext_dir = checkpoint.output_path / "pmc_fulltext"
    
    if not fulltext_dir.exists():
        return {"skipped": True, "reason": "no_fulltext_dir"}
    
    result = workflow_steps.run_data_scout(
        harvest_dir=fulltext_dir,
        gene_symbol=checkpoint.gene_symbol,
        min_relevance=0.3,
    )
    
    return result.stats


# =============================================================================
# Main Pipeline Runner
# =============================================================================

def run_resumable_pipeline(
    resume: bool = True,
    fresh: bool = False,
) -> JobCheckpoint:
    """Run the download-only pipeline with resume support."""
    from gui.worker import LoggingProgressCallback
    import uuid
    
    cm = CheckpointManager()
    
    # Determine run mode
    existing_run = find_latest_run_dir(CONFIG["output_dir"], CONFIG["gene_symbol"])
    
    if fresh or not existing_run:
        # Fresh start
        logger.info("Starting fresh pipeline run...")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        job_id = f"{CONFIG['gene_symbol'].lower()}_{timestamp}_{uuid.uuid4().hex[:6]}"
        
        checkpoint = JobCheckpoint(
            job_id=job_id,
            gene_symbol=CONFIG["gene_symbol"],
            email=CONFIG["email"],
            output_dir=CONFIG["output_dir"],
            max_pmids=CONFIG["max_pmids"],
            max_papers_to_download=CONFIG["max_papers_to_download"],
            auto_synonyms=True,
            enable_tier1=True,
            enable_tier2=True,
            tier2_confidence_threshold=CONFIG["tier2_confidence_threshold"],
            scout_enabled=True,
            use_pubmed=True,
            use_pubmind=True,
            use_europepmc=False,
        )
        state = {
            "run_dir": checkpoint.output_path,
            "discovered_pmids": [],
            "existing_abstracts": set(),
            "filter_completed": False,
            "filtered_pmids": [],
            "downloaded_pmids": set(),
        }
    else:
        # Resume existing run
        logger.info(f"Resuming from: {existing_run}")
        state = scan_existing_state(existing_run, CONFIG["gene_symbol"])
        
        resume_stage = determine_resume_stage(state)
        if resume_stage == PipelineStep.COMPLETED:
            logger.info("Pipeline already complete!")
            print_status(state)
            return None
        
        # Create checkpoint from existing state
        timestamp = existing_run.name
        job_id = f"{CONFIG['gene_symbol'].lower()}_{timestamp}_resume"
        
        checkpoint = JobCheckpoint(
            job_id=job_id,
            gene_symbol=CONFIG["gene_symbol"],
            email=CONFIG["email"],
            output_dir=CONFIG["output_dir"],
            max_pmids=CONFIG["max_pmids"],
            max_papers_to_download=CONFIG["max_papers_to_download"],
            auto_synonyms=True,
            enable_tier1=True,
            enable_tier2=True,
            tier2_confidence_threshold=CONFIG["tier2_confidence_threshold"],
            scout_enabled=True,
            use_pubmed=True,
            use_pubmind=True,
            use_europepmc=False,
            # Pre-populate from existing state
            discovered_pmids=state["discovered_pmids"],
            filtered_pmids=state["filtered_pmids"],
            downloaded_pmids=list(state["downloaded_pmids"]),
            current_step=resume_stage,
        )
        
        # Override output_path to use existing run dir
        checkpoint._existing_output_path = existing_run
    
    # Create output directory
    output_path = getattr(checkpoint, "_existing_output_path", None) or checkpoint.output_path
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Override output_path property for this run
    if hasattr(checkpoint, "_existing_output_path"):
        # Monkey-patch to use existing path
        JobCheckpoint.output_path = property(lambda self: self._existing_output_path)
        checkpoint._existing_output_path = output_path
    
    logger.info(f"Output directory: {output_path}")
    logger.info(f"Starting from step: {checkpoint.current_step.value}")
    
    start_time = datetime.now()
    
    # Define pipeline steps (download-only, no extraction)
    steps = [
        (PipelineStep.FETCHING_PMIDS, step_fetch_pmids_resumable),
        (PipelineStep.FETCHING_ABSTRACTS, step_fetch_abstracts_resumable),
        (PipelineStep.FILTERING_PAPERS, step_filter_papers_resumable),
        (PipelineStep.DOWNLOADING_FULLTEXT, step_download_fulltext_resumable),
        (PipelineStep.SCOUTING_DATA, step_scout_data_resumable),
    ]
    
    # Track step order for resume logic
    step_order = [s[0] for s in steps]
    current_idx = 0
    if checkpoint.current_step in step_order:
        current_idx = step_order.index(checkpoint.current_step)
    
    step_stats = {}
    
    try:
        for idx, (step, step_func) in enumerate(steps):
            # Skip steps before resume point
            if idx < current_idx:
                logger.info(f"Skipping {step.value} (already complete)")
                continue
            
            logger.info(f"\n{'='*60}")
            logger.info(f"STEP: {step.value}")
            logger.info(f"{'='*60}")
            
            checkpoint.current_step = step
            cm.save(checkpoint)
            
            stats = step_func(checkpoint, state)
            step_stats[step.value] = stats
            
            logger.info(f"Stats: {stats}")
            cm.save(checkpoint)
        
        # Mark complete
        checkpoint.current_step = PipelineStep.COMPLETED
        checkpoint.completed_at = datetime.now().isoformat()
        cm.save(checkpoint)
        
        elapsed = (datetime.now() - start_time).total_seconds()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"PIPELINE COMPLETE in {elapsed:.1f}s")
        logger.info(f"{'='*60}")
        logger.info(f"PMIDs discovered: {len(checkpoint.discovered_pmids)}")
        logger.info(f"PMIDs passed filter: {len(checkpoint.filtered_pmids)}")
        logger.info(f"Papers downloaded: {len(checkpoint.downloaded_pmids)}")
        logger.info(f"Output: {output_path}")
        
        return checkpoint
        
    except KeyboardInterrupt:
        logger.info("\nInterrupted! Progress saved. Run again to resume.")
        cm.save(checkpoint)
        raise
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        checkpoint.mark_failed(str(e), checkpoint.current_step)
        cm.save(checkpoint)
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Resumable KCNH2 download-only pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from last checkpoint (default behavior)"
    )
    parser.add_argument(
        "--fresh", action="store_true",
        help="Start fresh, ignoring any existing state"
    )
    parser.add_argument(
        "--status", action="store_true",
        help="Show current pipeline status and exit"
    )
    
    args = parser.parse_args()
    
    if args.status:
        existing_run = find_latest_run_dir(CONFIG["output_dir"], CONFIG["gene_symbol"])
        if existing_run:
            state = scan_existing_state(existing_run, CONFIG["gene_symbol"])
            print_status(state)
        else:
            print("No existing pipeline run found.")
        return
    
    try:
        run_resumable_pipeline(resume=not args.fresh, fresh=args.fresh)
    except KeyboardInterrupt:
        sys.exit(1)
    except Exception:
        sys.exit(1)


if __name__ == "__main__":
    main()
