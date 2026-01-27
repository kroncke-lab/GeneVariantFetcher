"""
FastAPI server for GeneVariantFetcher GUI.

Provides:
- REST API for job management
- WebSocket for real-time progress updates
- Static file serving for frontend
"""

import asyncio
import json
import logging
import os
import threading
from contextlib import asynccontextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from fastapi import (
    FastAPI,
    HTTPException,
    WebSocket,
    WebSocketDisconnect,
    BackgroundTasks,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles

from gui.checkpoint import CheckpointManager, JobCheckpoint, PipelineStep
from gui.models import (
    BrowserFetchRequest,
    BrowserFetchStatus,
    FolderAnalysis,
    FolderJobRequest,
    JobCreateRequest,
    JobListResponse,
    JobResponse,
    PipelineStageBreakdown,
    PmidFilteredEntry,
    PmidPaywalledEntry,
    ResumeJobRequest,
)
from gui.worker import PipelineWorker, ProgressCallback, create_job
from gui.routes import settings_router, browser_router

# Import post-processing functions for browser fetch
try:
    from cli.fetch_manager import convert_to_markdown, run_data_scout

    FETCH_MANAGER_AVAILABLE = True
except ImportError:
    FETCH_MANAGER_AVAILABLE = False
    convert_to_markdown = None
    run_data_scout = None

# Path to .env file
ENV_FILE_PATH = Path(__file__).resolve().parent.parent / ".env"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# Track running browser fetch tasks
_browser_fetch_tasks: Dict[str, BrowserFetchStatus] = {}
# Store full logs separately (can be large)
_browser_fetch_logs: Dict[str, List[str]] = {}


# =============================================================================
# WebSocket Connection Manager
# =============================================================================


class ConnectionManager:
    """Manages WebSocket connections for progress updates."""

    def __init__(self):
        self.active_connections: Dict[str, Set[WebSocket]] = {}
        self._lock = threading.Lock()

    async def connect(self, websocket: WebSocket, job_id: str):
        """Accept and register a WebSocket connection."""
        await websocket.accept()
        with self._lock:
            if job_id not in self.active_connections:
                self.active_connections[job_id] = set()
            self.active_connections[job_id].add(websocket)
        logger.info(f"WebSocket connected for job {job_id}")

    def disconnect(self, websocket: WebSocket, job_id: str):
        """Remove a WebSocket connection."""
        with self._lock:
            if job_id in self.active_connections:
                self.active_connections[job_id].discard(websocket)
                if not self.active_connections[job_id]:
                    del self.active_connections[job_id]
        logger.info(f"WebSocket disconnected for job {job_id}")

    async def broadcast(self, job_id: str, message: Dict[str, Any]):
        """Broadcast message to all connections for a job."""
        with self._lock:
            connections = list(self.active_connections.get(job_id, set()))

        for websocket in connections:
            try:
                await websocket.send_json(message)
            except Exception:
                self.disconnect(websocket, job_id)


# Global instances
checkpoint_manager = CheckpointManager()
connection_manager = ConnectionManager()
active_workers: Dict[str, PipelineWorker] = {}


# =============================================================================
# Progress Callback for WebSocket
# =============================================================================


class WebSocketProgressCallback(ProgressCallback):
    """Progress callback that broadcasts to WebSocket clients."""

    def __init__(self, job_id: str, loop: asyncio.AbstractEventLoop):
        self.job_id = job_id
        self.loop = loop

    def _broadcast(self, message: Dict[str, Any]):
        """Send message to all connected clients."""
        asyncio.run_coroutine_threadsafe(
            connection_manager.broadcast(self.job_id, message),
            self.loop,
        )

    def on_step_start(self, step: PipelineStep, message: str):
        self._broadcast(
            {
                "type": "step_start",
                "step": step.value,
                "step_display": PipelineStep.get_display_name(step),
                "message": message,
                "timestamp": datetime.now().isoformat(),
            }
        )

    def on_step_complete(self, step: PipelineStep, message: str, stats: Dict[str, Any]):
        self._broadcast(
            {
                "type": "step_complete",
                "step": step.value,
                "step_display": PipelineStep.get_display_name(step),
                "message": message,
                "stats": stats,
                "timestamp": datetime.now().isoformat(),
            }
        )

    def on_log(self, message: str):
        self._broadcast(
            {
                "type": "log",
                "message": message,
                "timestamp": datetime.now().isoformat(),
            }
        )

    def on_error(self, step: PipelineStep, error: str):
        self._broadcast(
            {
                "type": "error",
                "step": step.value,
                "error": error,
                "timestamp": datetime.now().isoformat(),
            }
        )


# =============================================================================
# Background Job Runner
# =============================================================================


def run_job_in_background(job_id: str, resume: bool, loop: asyncio.AbstractEventLoop):
    """Run a pipeline job in a background thread."""
    try:
        checkpoint = checkpoint_manager.load(job_id)
        if not checkpoint:
            logger.error(f"Job {job_id} not found")
            return

        # Create worker with WebSocket progress callback
        callback = WebSocketProgressCallback(job_id, loop)
        worker = PipelineWorker(
            checkpoint_manager=checkpoint_manager,
            progress_callback=callback,
        )

        # Store worker reference for stop requests
        active_workers[job_id] = worker

        try:
            worker.run(checkpoint, resume=resume)
        finally:
            active_workers.pop(job_id, None)

        # Broadcast completion
        final_checkpoint = checkpoint_manager.load(job_id)
        asyncio.run_coroutine_threadsafe(
            connection_manager.broadcast(
                job_id,
                {
                    "type": "job_complete",
                    "status": final_checkpoint.current_step.value
                    if final_checkpoint
                    else "unknown",
                    "timestamp": datetime.now().isoformat(),
                },
            ),
            loop,
        )

    except Exception as e:
        logger.exception(f"Job {job_id} failed: {e}")
        asyncio.run_coroutine_threadsafe(
            connection_manager.broadcast(
                job_id,
                {
                    "type": "error",
                    "error": str(e),
                    "timestamp": datetime.now().isoformat(),
                },
            ),
            loop,
        )


# =============================================================================
# FastAPI App
# =============================================================================


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown events."""
    # Check for incomplete jobs on startup
    incomplete = checkpoint_manager.get_incomplete_jobs()
    if incomplete:
        logger.warning(f"Found {len(incomplete)} incomplete jobs that can be resumed")
        for job in incomplete:
            logger.warning(
                f"  - {job.job_id}: {job.gene_symbol} at step {job.current_step.value}"
            )
    yield
    # Cleanup on shutdown
    for job_id, worker in list(active_workers.items()):
        logger.info(f"Requesting stop for job {job_id}")
        worker.request_stop()


app = FastAPI(
    title="GeneVariantFetcher GUI",
    description="Web interface for extracting genetic variant data from literature",
    version="1.0.0",
    lifespan=lifespan,
)

# CORS for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve static files
static_dir = Path(__file__).parent / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

# Include routers from gui.routes
app.include_router(settings_router)
app.include_router(browser_router)


# =============================================================================
# Helper Functions
# =============================================================================


def checkpoint_to_response(checkpoint: JobCheckpoint) -> JobResponse:
    """Convert checkpoint to API response."""
    stats = {
        "pmids_discovered": len(checkpoint.discovered_pmids),
        "pmids_filtered": len(checkpoint.filtered_pmids),
        "papers_downloaded": len(checkpoint.downloaded_pmids),
        "papers_extracted": len(checkpoint.extracted_pmids),
    }

    # Include total_variants from extraction or aggregation step
    extraction_stats = checkpoint.step_progress.get("extracting_variants", {})
    aggregation_stats = checkpoint.step_progress.get("aggregating_data", {})
    if "total_variants" in extraction_stats:
        stats["total_variants"] = extraction_stats["total_variants"]
    elif "variants_aggregated" in aggregation_stats:
        stats["total_variants"] = aggregation_stats["variants_aggregated"]

    # Merge current step progress
    stats.update(checkpoint.step_progress.get(checkpoint.current_step.value, {}))

    return JobResponse(
        job_id=checkpoint.job_id,
        gene_symbol=checkpoint.gene_symbol,
        current_step=checkpoint.current_step.value,
        current_step_display=PipelineStep.get_display_name(checkpoint.current_step),
        created_at=checkpoint.created_at,
        updated_at=checkpoint.updated_at,
        started_at=checkpoint.started_at,
        completed_at=checkpoint.completed_at,
        error_message=checkpoint.error_message,
        is_resumable=checkpoint.is_resumable,
        output_dir=str(checkpoint.output_path),
        statistics=stats,
    )


# =============================================================================
# API Routes
# =============================================================================


@app.get("/", response_class=HTMLResponse)
async def root():
    """Serve the main UI."""
    index_path = static_dir / "index.html"
    if index_path.exists():
        return FileResponse(index_path)
    return HTMLResponse("<h1>GeneVariantFetcher GUI</h1><p>Static files not found.</p>")


@app.get("/api/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "ok", "timestamp": datetime.now().isoformat()}


@app.get("/api/jobs", response_model=JobListResponse)
async def list_jobs():
    """List all jobs."""
    jobs = checkpoint_manager.list_jobs()
    incomplete = [j for j in jobs if j.is_resumable]

    return JobListResponse(
        jobs=[checkpoint_to_response(j) for j in jobs],
        incomplete_count=len(incomplete),
    )


@app.get("/api/jobs/incomplete", response_model=JobListResponse)
async def list_incomplete_jobs():
    """List jobs that can be resumed."""
    jobs = checkpoint_manager.get_incomplete_jobs()
    return JobListResponse(
        jobs=[checkpoint_to_response(j) for j in jobs],
        incomplete_count=len(jobs),
    )


@app.get("/api/jobs/{job_id}", response_model=JobResponse)
async def get_job(job_id: str):
    """Get job details."""
    checkpoint = checkpoint_manager.load(job_id)
    if not checkpoint:
        raise HTTPException(status_code=404, detail="Job not found")
    return checkpoint_to_response(checkpoint)


@app.post("/api/jobs", response_model=JobResponse)
async def create_new_job(request: JobCreateRequest, background_tasks: BackgroundTasks):
    """Create and start a new pipeline job."""
    # Create checkpoint
    checkpoint = create_job(
        gene_symbol=request.gene_symbol,
        email=request.email,
        output_dir=request.output_dir,
        max_pmids=request.max_pmids,
        max_papers_to_download=request.max_papers_to_download,
        tier_threshold=request.tier_threshold,
        use_clinical_triage=request.use_clinical_triage,
        auto_synonyms=request.auto_synonyms,
        synonyms=request.synonyms,
        specific_pmids=request.specific_pmids,
        enable_tier1=request.enable_tier1,
        enable_tier2=request.enable_tier2,
        use_pubmind=request.use_pubmind,
        use_pubmed=request.use_pubmed,
        use_europepmc=request.use_europepmc,
        tier2_confidence_threshold=request.tier2_confidence_threshold,
        scout_enabled=request.scout_enabled,
    )

    # Save initial checkpoint
    checkpoint_manager.save(checkpoint)

    # Start job in background
    loop = asyncio.get_event_loop()
    background_tasks.add_task(
        lambda: threading.Thread(
            target=run_job_in_background,
            args=(checkpoint.job_id, False, loop),
            daemon=True,
        ).start()
    )

    return checkpoint_to_response(checkpoint)


@app.post("/api/jobs/{job_id}/resume", response_model=JobResponse)
async def resume_job(job_id: str, background_tasks: BackgroundTasks):
    """Resume an interrupted job."""
    checkpoint = checkpoint_manager.load(job_id)
    if not checkpoint:
        raise HTTPException(status_code=404, detail="Job not found")

    if not checkpoint.is_resumable:
        raise HTTPException(
            status_code=400,
            detail=f"Job cannot be resumed (status: {checkpoint.current_step.value})",
        )

    if job_id in active_workers:
        raise HTTPException(status_code=400, detail="Job is already running")

    # Start job in background
    loop = asyncio.get_event_loop()
    background_tasks.add_task(
        lambda: threading.Thread(
            target=run_job_in_background,
            args=(job_id, True, loop),
            daemon=True,
        ).start()
    )

    return checkpoint_to_response(checkpoint)


@app.post("/api/jobs/{job_id}/stop")
async def stop_job(job_id: str):
    """Request graceful stop of a running job."""
    if job_id not in active_workers:
        raise HTTPException(status_code=400, detail="Job is not running")

    active_workers[job_id].request_stop()
    return {"status": "stop_requested", "job_id": job_id}


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job and its checkpoint."""
    if job_id in active_workers:
        raise HTTPException(status_code=400, detail="Cannot delete running job")

    if checkpoint_manager.delete(job_id):
        return {"status": "deleted", "job_id": job_id}
    else:
        raise HTTPException(status_code=404, detail="Job not found")


@app.get("/api/jobs/{job_id}/logs")
async def get_job_logs(job_id: str, limit: int = 100):
    """Get recent log lines for a job."""
    checkpoint = checkpoint_manager.load(job_id)
    if not checkpoint:
        raise HTTPException(status_code=404, detail="Job not found")

    return {
        "job_id": job_id,
        "logs": checkpoint.log_lines[-limit:],
    }


@app.get("/api/jobs/{job_id}/results")
async def get_job_results(job_id: str):
    """Get results for a completed job."""
    checkpoint = checkpoint_manager.load(job_id)
    if not checkpoint:
        raise HTTPException(status_code=404, detail="Job not found")

    if checkpoint.current_step != PipelineStep.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed")

    output_path = checkpoint.output_path
    summary_file = output_path / f"{checkpoint.gene_symbol}_workflow_summary.json"
    penetrance_file = output_path / f"{checkpoint.gene_symbol}_penetrance_summary.json"

    results = {
        "job_id": job_id,
        "output_path": str(output_path),
        "files": {},
    }

    if summary_file.exists():
        with open(summary_file) as f:
            results["workflow_summary"] = json.load(f)

    if penetrance_file.exists():
        with open(penetrance_file) as f:
            results["penetrance_summary"] = json.load(f)

    # List output files
    results["files"] = {
        "database": str(output_path / f"{checkpoint.gene_symbol}.db"),
        "extractions": [str(f) for f in (output_path / "extractions").glob("*.json")]
        if (output_path / "extractions").exists()
        else [],
        "pmid_status": [str(f) for f in (output_path / "pmid_status").glob("*.csv")]
        if (output_path / "pmid_status").exists()
        else [],
    }

    return results


# =============================================================================
# Folder Analysis API (for importing existing downloads)
# =============================================================================


def _extract_pmid_from_filename(filename: str) -> Optional[str]:
    """Extract PMID from various filename formats."""
    import re

    # Patterns: 12345678_FULL_CONTEXT.md, PMID_12345678_FULL_CONTEXT.md, etc.
    patterns = [
        r"^(\d{6,10})_(?:FULL_CONTEXT|DATA_ZONES)\.md$",
        r"^PMID_(\d{6,10})_(?:FULL_CONTEXT|DATA_ZONES)\.md$",
        r"^[A-Z]+_PMID_(\d{6,10})\.json$",  # SCN5A_PMID_12345678.json
    ]
    for pattern in patterns:
        match = re.match(pattern, filename, re.IGNORECASE)
        if match:
            return match.group(1)
    return None


def _extract_gene_from_extraction_filename(filename: str) -> Optional[str]:
    """Extract gene symbol from extraction filename."""
    import re

    # Pattern: GENE_PMID_12345678.json
    match = re.match(r"^([A-Z0-9]+)_PMID_\d+\.json$", filename, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return None


def _build_pipeline_breakdown(
    folder_path: Path,
    detected_gene: Optional[str],
    actual_fulltext_pmids: Optional[List[str]] = None,
) -> PipelineStageBreakdown:
    """
    Build a detailed breakdown of PMID status through each pipeline stage.

    This helps identify where the pipeline may have been interrupted vs where
    PMIDs were deliberately filtered out.

    Args:
        folder_path: Path to the job folder
        detected_gene: Gene symbol detected from folder structure
        actual_fulltext_pmids: PMIDs from actual *_FULL_CONTEXT.md files (source of truth)
    """
    import csv

    breakdown = PipelineStageBreakdown()

    # 1. Discovery stage: Read {GENE}_pmids.txt
    discovered_pmids = set()
    pmids_file = None

    # Try to find pmids file with detected gene first
    if detected_gene:
        pmids_file = folder_path / f"{detected_gene}_pmids.txt"
        if not pmids_file.exists():
            pmids_file = None

    # Fall back to glob pattern
    if pmids_file is None:
        pmid_files = list(folder_path.glob("*_pmids.txt"))
        # Filter out source-specific files (pubmind, pubmed)
        pmid_files = [
            f
            for f in pmid_files
            if "_pubmind" not in f.name and "_pubmed" not in f.name
        ]
        if pmid_files:
            pmids_file = pmid_files[0]

    if pmids_file and pmids_file.exists():
        try:
            with open(pmids_file, "r") as f:
                for line in f:
                    pmid = line.strip()
                    if pmid and pmid.isdigit():
                        discovered_pmids.add(pmid)
            breakdown.has_discovery_data = True
            breakdown.discovered_count = len(discovered_pmids)
            breakdown.discovered_pmids = sorted(discovered_pmids, key=int)
        except Exception:
            pass

    # 2. Harvesting stage: Read abstract_json/ directory
    harvested_pmids = set()
    abstract_dir = folder_path / "abstract_json"
    if abstract_dir.exists():
        for json_file in abstract_dir.glob("*.json"):
            pmid = json_file.stem
            if pmid.isdigit():
                harvested_pmids.add(pmid)
        breakdown.has_harvest_data = True
        breakdown.harvested_count = len(harvested_pmids)
        breakdown.harvested_pmids = sorted(harvested_pmids, key=int)

    # 3. Filtering stage: Read pmid_status/filtered_out.csv
    filtered_out_pmids = set()
    filtered_entries = []
    filtered_out_file = folder_path / "pmid_status" / "filtered_out.csv"
    if filtered_out_file.exists():
        try:
            with open(filtered_out_file, "r", newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pmid = row.get("PMID", "").strip()
                    reason = row.get("Reason", "Unknown").strip()
                    if pmid:
                        filtered_out_pmids.add(pmid)
                        filtered_entries.append(
                            PmidFilteredEntry(pmid=pmid, reason=reason)
                        )
            breakdown.has_filter_data = True
            breakdown.filtered_out_count = len(filtered_out_pmids)
            breakdown.filtered_out_entries = filtered_entries
        except Exception:
            pass

    # 4. Download stage: Use actual files as source of truth, fall back to CSV
    fulltext_dir = folder_path / "pmc_fulltext"

    # Prefer actual FULL_CONTEXT files over CSV (CSV may be incomplete/stale)
    downloaded_pmids = set()
    if actual_fulltext_pmids:
        # Use actual files found on disk - this is the authoritative source
        downloaded_pmids = set(actual_fulltext_pmids)
        breakdown.has_download_data = True
        breakdown.downloaded_count = len(downloaded_pmids)
        breakdown.downloaded_pmids = sorted(
            downloaded_pmids, key=lambda x: int(x) if x.isdigit() else 0
        )
    else:
        # Fall back to CSV if no actual files provided
        success_log = fulltext_dir / "successful_downloads.csv"
        if success_log.exists():
            try:
                with open(success_log, "r", newline="", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        pmid = row.get("PMID", "").strip()
                        if pmid:
                            downloaded_pmids.add(pmid)
                breakdown.has_download_data = True
                breakdown.downloaded_count = len(downloaded_pmids)
                breakdown.downloaded_pmids = sorted(
                    downloaded_pmids, key=lambda x: int(x) if x.isdigit() else 0
                )
            except Exception:
                pass

    # Paywalled/missing - exclude PMIDs that have actual downloaded files
    paywalled_pmids = set()
    paywalled_entries = []
    paywalled_log = fulltext_dir / "paywalled_missing.csv"
    if paywalled_log.exists():
        try:
            with open(paywalled_log, "r", newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pmid = row.get("PMID", "").strip()
                    reason = row.get("Reason", "Unknown").strip()
                    url = row.get("URL", "").strip() or None
                    # Skip if this PMID actually has a downloaded file
                    if pmid and pmid not in downloaded_pmids:
                        paywalled_pmids.add(pmid)
                        paywalled_entries.append(
                            PmidPaywalledEntry(pmid=pmid, reason=reason, url=url)
                        )
            breakdown.paywalled_count = len(paywalled_pmids)
            breakdown.paywalled_entries = paywalled_entries
        except Exception:
            pass

    # 5. Calculate gaps (PMIDs lost at each stage that weren't deliberately filtered)

    # Harvest gap: discovered but not harvested (and not filtered)
    if breakdown.has_discovery_data and breakdown.has_harvest_data:
        expected_harvested = discovered_pmids - filtered_out_pmids
        harvest_gap = expected_harvested - harvested_pmids
        breakdown.harvest_gap_count = len(harvest_gap)
        breakdown.harvest_gap_pmids = sorted(
            harvest_gap, key=lambda x: int(x) if x.isdigit() else 0
        )

    # Download gap: harvested but not downloaded (and not paywalled, not filtered)
    if breakdown.has_harvest_data:
        # PMIDs that passed filters = harvested - filtered_out
        pmids_to_download = harvested_pmids - filtered_out_pmids
        # Expected: should be either downloaded or paywalled
        expected_accounted = downloaded_pmids | paywalled_pmids
        download_gap = pmids_to_download - expected_accounted
        breakdown.download_gap_count = len(download_gap)
        breakdown.download_gap_pmids = sorted(
            download_gap, key=lambda x: int(x) if x.isdigit() else 0
        )

    # 6. Detect likely interruption
    if breakdown.has_discovery_data:
        total_discovered = breakdown.discovered_count

        # Check harvest stage
        if breakdown.harvest_gap_count > 0:
            gap_pct = (breakdown.harvest_gap_count / total_discovered) * 100
            if gap_pct > 5:  # More than 5% lost = likely interruption
                breakdown.likely_interrupted = True
                breakdown.interruption_stage = "harvesting"
                breakdown.interruption_details = (
                    f"{breakdown.harvest_gap_count} PMIDs ({gap_pct:.1f}%) discovered but not harvested. "
                    f"Process may have been interrupted during abstract fetching. "
                    f"Last harvested: {breakdown.harvested_pmids[-1] if breakdown.harvested_pmids else 'N/A'}"
                )

        # Check download stage
        elif breakdown.download_gap_count > 0:
            gap_pct = (breakdown.download_gap_count / total_discovered) * 100
            if gap_pct > 5:
                breakdown.likely_interrupted = True
                breakdown.interruption_stage = "downloading"
                breakdown.interruption_details = (
                    f"{breakdown.download_gap_count} PMIDs ({gap_pct:.1f}%) passed filters but not downloaded. "
                    f"Process may have been interrupted during full-text download."
                )

    return breakdown


@app.get("/api/analyze-folder", response_model=FolderAnalysis)
async def analyze_folder(path: str):
    """
    Analyze an existing folder to detect downloaded papers, scouting status, and extractions.

    This endpoint helps users understand what files exist in a folder and what processing
    steps are needed to complete extraction.
    """
    folder_path = Path(path).expanduser().resolve()

    # Security check
    allowed_roots = [Path.home(), Path("/tmp"), Path.cwd()]
    is_allowed = any(
        folder_path == root or root in folder_path.parents for root in allowed_roots
    )
    if not is_allowed:
        raise HTTPException(
            status_code=403, detail="Access to this directory is not allowed"
        )

    if not folder_path.exists():
        return FolderAnalysis(
            path=str(folder_path),
            is_valid=False,
            status_message="Directory does not exist",
        )

    if not folder_path.is_dir():
        return FolderAnalysis(
            path=str(folder_path),
            is_valid=False,
            status_message="Path is not a directory",
        )

    # Look for pmc_fulltext subdirectory or files directly in the folder
    fulltext_dir = folder_path / "pmc_fulltext"
    if fulltext_dir.exists():
        search_dir = fulltext_dir
    else:
        search_dir = folder_path

    # Find FULL_CONTEXT and DATA_ZONES files
    full_context_files = list(search_dir.glob("*_FULL_CONTEXT.md"))
    data_zones_files = list(search_dir.glob("*_DATA_ZONES.md"))

    full_context_pmids = []
    for f in full_context_files:
        pmid = _extract_pmid_from_filename(f.name)
        if pmid:
            full_context_pmids.append(pmid)

    data_zones_pmids = []
    for f in data_zones_files:
        pmid = _extract_pmid_from_filename(f.name)
        if pmid:
            data_zones_pmids.append(pmid)

    # Look for extractions directory
    extraction_dir = folder_path / "extractions"
    extraction_pmids = []
    detected_genes = set()

    if extraction_dir.exists():
        extraction_files = list(extraction_dir.glob("*_PMID_*.json"))
        for f in extraction_files:
            gene = _extract_gene_from_extraction_filename(f.name)
            if gene:
                detected_genes.add(gene)
            # Extract PMID from filename
            import re

            match = re.search(r"PMID_(\d+)", f.name)
            if match:
                extraction_pmids.append(match.group(1))

    # Also try to detect gene from folder name (e.g., output/SCN5A/20240101_120000)
    folder_gene = None
    if (
        folder_path.parent.name
        and folder_path.parent.name.isupper()
        and len(folder_path.parent.name) <= 10
    ):
        folder_gene = folder_path.parent.name
        detected_genes.add(folder_gene)

    # Determine status
    has_fulltext = len(full_context_pmids) > 0 or len(data_zones_pmids) > 0
    has_scouted = len(data_zones_pmids) > 0
    all_pmids = set(full_context_pmids) | set(data_zones_pmids)
    scouting_complete = has_scouted and set(data_zones_pmids) >= set(full_context_pmids)
    has_extractions = len(extraction_pmids) > 0

    # Build recommendations
    recommendations = []
    if not has_fulltext:
        recommendations.append(
            "No downloaded papers found. Run a full pipeline or download papers first."
        )
    elif not has_scouted:
        recommendations.append(
            "Papers not yet scouted. Enable 'Run Data Scout' to create condensed DATA_ZONES files for faster extraction."
        )
    elif not scouting_complete:
        unscouted = set(full_context_pmids) - set(data_zones_pmids)
        recommendations.append(
            f"{len(unscouted)} papers have not been scouted yet. Enable 'Run Data Scout' to process them."
        )

    if has_fulltext and not has_extractions:
        recommendations.append(
            "Ready for extraction. Start extraction to process papers."
        )
    elif has_extractions:
        unextracted = all_pmids - set(extraction_pmids)
        if unextracted:
            recommendations.append(
                f"{len(unextracted)} papers have not been extracted yet."
            )
        else:
            recommendations.append(
                "All papers have been extracted. You can re-run extraction if needed."
            )

    # Build status message
    status_parts = []
    if full_context_pmids:
        status_parts.append(f"{len(full_context_pmids)} full-text papers")
    if data_zones_pmids:
        status_parts.append(f"{len(data_zones_pmids)} scouted")
    if extraction_pmids:
        status_parts.append(f"{len(extraction_pmids)} extracted")

    status_message = (
        "Found: " + ", ".join(status_parts) if status_parts else "No papers found"
    )

    # Build pipeline breakdown for debugging
    detected_gene = folder_gene or (
        list(detected_genes)[0] if len(detected_genes) == 1 else None
    )
    pipeline_breakdown = _build_pipeline_breakdown(
        folder_path, detected_gene, actual_fulltext_pmids=full_context_pmids
    )

    # Add recommendations based on pipeline breakdown
    if pipeline_breakdown.likely_interrupted:
        recommendations.insert(
            0,
            f"Pipeline may have been interrupted at {pipeline_breakdown.interruption_stage} stage. {pipeline_breakdown.interruption_details}",
        )
    elif (
        pipeline_breakdown.has_discovery_data
        and pipeline_breakdown.harvest_gap_count > 0
    ):
        recommendations.append(
            f"{pipeline_breakdown.harvest_gap_count} PMIDs were discovered but not harvested."
        )
    elif pipeline_breakdown.download_gap_count > 0:
        recommendations.append(
            f"{pipeline_breakdown.download_gap_count} PMIDs passed filters but were not downloaded or marked as paywalled."
        )

    return FolderAnalysis(
        path=str(folder_path),
        is_valid=has_fulltext,
        gene_symbol=detected_gene,
        detected_gene_symbols=sorted(detected_genes),
        full_context_count=len(full_context_pmids),
        data_zones_count=len(data_zones_pmids),
        extraction_count=len(extraction_pmids),
        full_context_pmids=full_context_pmids,
        data_zones_pmids=data_zones_pmids,
        extraction_pmids=extraction_pmids,
        has_fulltext=has_fulltext,
        has_scouted=has_scouted,
        scouting_complete=scouting_complete,
        has_extractions=has_extractions,
        status_message=status_message,
        recommendations=recommendations,
        pipeline_breakdown=pipeline_breakdown,
    )


@app.post("/api/jobs/from-folder", response_model=JobResponse)
async def create_job_from_folder(
    request: FolderJobRequest, background_tasks: BackgroundTasks
):
    """
    Create a job that processes an existing folder with downloaded papers.

    This skips the discovery, abstract fetching, filtering, and download steps,
    starting directly at extraction (or optionally scouting first).
    """
    from gui.worker import create_folder_job

    folder_path = Path(request.folder_path).expanduser().resolve()

    # Validate folder
    analysis = await analyze_folder(str(folder_path))
    if not analysis.is_valid:
        raise HTTPException(
            status_code=400, detail=f"Invalid folder: {analysis.status_message}"
        )

    # Create checkpoint for folder job
    checkpoint = create_folder_job(
        folder_path=str(folder_path),
        gene_symbol=request.gene_symbol,
        email=request.email,
        run_scout=request.run_scout,
        skip_already_extracted=request.skip_already_extracted,
        tier_threshold=request.tier_threshold,
        scout_enabled=request.scout_enabled,
        scout_min_relevance=request.scout_min_relevance,
        full_context_pmids=analysis.full_context_pmids,
        data_zones_pmids=analysis.data_zones_pmids,
        extraction_pmids=analysis.extraction_pmids
        if request.skip_already_extracted
        else [],
    )

    # Save initial checkpoint
    checkpoint_manager.save(checkpoint)

    # Start job in background
    loop = asyncio.get_event_loop()
    background_tasks.add_task(
        lambda: threading.Thread(
            target=run_job_in_background,
            args=(checkpoint.job_id, False, loop),
            daemon=True,
        ).start()
    )

    return checkpoint_to_response(checkpoint)


@app.post("/api/jobs/resume-folder", response_model=JobResponse)
async def resume_job_from_folder(
    request: ResumeJobRequest, background_tasks: BackgroundTasks
):
    """
    Resume an interrupted pipeline from a folder.

    This endpoint handles two main scenarios:
    1. Resume downloading - when the pipeline was interrupted during download stage
    2. Resume extraction - when downloading is complete but extraction was interrupted

    The job will continue from the specified stage and run through completion.
    """
    from gui.worker import create_resume_job

    folder_path = Path(request.folder_path).expanduser().resolve()

    # Security check
    allowed_roots = [Path.home(), Path("/tmp"), Path.cwd()]
    is_allowed = any(
        folder_path == root or root in folder_path.parents for root in allowed_roots
    )
    if not is_allowed:
        raise HTTPException(
            status_code=403, detail="Access to this directory is not allowed"
        )

    if not folder_path.exists():
        raise HTTPException(status_code=400, detail="Folder does not exist")

    # Get current folder analysis for context
    analysis = await analyze_folder(str(folder_path))

    # Validate resume stage
    valid_stages = ["downloading", "extraction"]
    if request.resume_stage not in valid_stages:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid resume stage. Must be one of: {valid_stages}",
        )

    # For downloading stage, we need PMIDs to download
    if request.resume_stage == "downloading" and not request.pmids_to_download:
        raise HTTPException(
            status_code=400,
            detail="PMIDs to download must be provided when resuming from downloading stage",
        )

    # Create checkpoint for resume job
    checkpoint = create_resume_job(
        folder_path=str(folder_path),
        gene_symbol=request.gene_symbol,
        email=request.email,
        resume_stage=request.resume_stage,
        pmids_to_download=request.pmids_to_download,
        max_papers_to_download=request.max_papers_to_download,
        run_scout=request.run_scout,
        skip_already_extracted=request.skip_already_extracted,
        tier_threshold=request.tier_threshold,
        scout_enabled=request.scout_enabled,
        scout_min_relevance=request.scout_min_relevance,
        full_context_pmids=analysis.full_context_pmids,
        data_zones_pmids=analysis.data_zones_pmids,
        extraction_pmids=analysis.extraction_pmids
        if request.skip_already_extracted
        else [],
    )

    # Save initial checkpoint
    checkpoint_manager.save(checkpoint)

    # Start job in background
    loop = asyncio.get_event_loop()
    background_tasks.add_task(
        lambda: threading.Thread(
            target=run_job_in_background,
            args=(checkpoint.job_id, False, loop),
            daemon=True,
        ).start()
    )

    return checkpoint_to_response(checkpoint)


# =============================================================================
# Browser Fetch for Paywalled Papers
# =============================================================================


def _run_browser_fetch(
    task_id: str,
    csv_path: str,
    headless: bool,
    max_papers: Optional[int],
    use_claude: bool,
    wait_for_captcha: bool = False,
    retry_failures_only: bool = False,
    gene_symbol: Optional[str] = None,
    auto_process: bool = True,
    fallback_to_manual: bool = True,
):
    """Run browser fetch in a background thread with real-time log streaming."""
    import subprocess
    import sys
    import re

    # Initialize logs storage
    _browser_fetch_logs[task_id] = []
    max_recent_logs = 20  # Keep last N lines in status for quick polling

    def add_log(line: str):
        """Add a log line and update status."""
        line = line.strip()
        if not line:
            return

        # Store in full logs
        _browser_fetch_logs[task_id].append(line)

        # Update recent logs (keep last N)
        task = _browser_fetch_tasks[task_id]
        recent = list(task.recent_logs)
        recent.append(line)
        if len(recent) > max_recent_logs:
            recent = recent[-max_recent_logs:]
        task.recent_logs = recent

        # Parse line for progress info
        # Match: [1/3] Processing PMID 12345678
        progress_match = re.search(r"\[(\d+)/(\d+)\]\s*Processing\s*PMID\s*(\d+)", line)
        if progress_match:
            task.progress = int(progress_match.group(1))
            task.total = int(progress_match.group(2))
            task.current_pmid = progress_match.group(3)

        # Match: >>> STEP 1: Direct PDF URL detected
        step_match = re.search(r">>>\s*(STEP\s*\d+[^:]*:?\s*[^\n]*)", line)
        if step_match:
            task.current_step = step_match.group(1).strip()

        # Match: Fetching PMID 12345678
        pmid_match = re.search(r"Fetching\s*PMID\s*(\d+)", line)
        if pmid_match:
            task.current_pmid = pmid_match.group(1)

        # Match: SUCCESS: Downloaded N file(s)
        if "SUCCESS:" in line and "Downloaded" in line:
            task.successful = (task.successful or 0) + 1

        # Match: FAILED:
        if "FAILED:" in line or "WARNING" in line and "FAILED" in line:
            task.failed = (task.failed or 0) + 1

    try:
        _browser_fetch_tasks[task_id].status = "running"
        add_log(f"Starting browser fetch for {csv_path}")

        # Build command
        cmd = [
            sys.executable,
            "-u",
            "browser_fetch.py",
            csv_path,
        ]  # -u for unbuffered output
        if headless:
            cmd.append("--headless")
        if max_papers:
            cmd.extend(["--max", str(max_papers)])
        if use_claude:
            cmd.append("--use-claude")
        if wait_for_captcha:
            cmd.append("--wait-for-captcha")
        if retry_failures_only:
            cmd.append("--retry-failures")
        if not fallback_to_manual:
            cmd.append("--no-fallback-manual")

        add_log(f"Command: {' '.join(cmd)}")

        # Run with Popen for real-time output
        # Set PYTHONUNBUFFERED to ensure subprocess doesn't buffer
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Combine stderr into stdout
            text=True,
            bufsize=1,  # Line buffered
            cwd=Path(__file__).parent.parent,
            env=env,
        )

        # Read output line by line using readline() for immediate reads
        try:
            while True:
                line = process.stdout.readline()
                if not line and process.poll() is not None:
                    break
                if line:
                    add_log(line)
                    logger.info(f"[BrowserFetch {task_id}] {line.strip()}")
        except Exception as e:
            add_log(f"Error reading output: {e}")

        # Wait for process to complete
        return_code = process.wait(timeout=3600)

        # Parse final stats from logs
        full_output = "\n".join(_browser_fetch_logs[task_id])
        success_match = re.search(r"Successful:\s*(\d+)/(\d+)", full_output)
        if success_match:
            _browser_fetch_tasks[task_id].successful = int(success_match.group(1))
            _browser_fetch_tasks[task_id].total = int(success_match.group(2))
            _browser_fetch_tasks[task_id].failed = (
                _browser_fetch_tasks[task_id].total
                - _browser_fetch_tasks[task_id].successful
            )

        if return_code == 0:
            add_log("Browser fetch completed successfully")

            # Post-processing: convert to markdown and run scout
            if auto_process and FETCH_MANAGER_AVAILABLE and gene_symbol:
                add_log(
                    "Starting post-processing: converting to markdown and running scout..."
                )
                _browser_fetch_tasks[task_id].current_step = "Post-processing"

                target_dir = Path(csv_path).parent
                converted = 0
                scouted = 0

                # Find PMIDs with downloaded files (Main_Text.pdf) but no FULL_CONTEXT.md
                main_pdfs = list(target_dir.glob("*_Main_Text.pdf"))
                for pdf_path in main_pdfs:
                    pmid = pdf_path.name.split("_")[0]
                    full_context = target_dir / f"{pmid}_FULL_CONTEXT.md"

                    if not full_context.exists():
                        add_log(f"Converting PMID {pmid} to markdown...")
                        try:
                            result = convert_to_markdown(pmid, target_dir)
                            if result:
                                converted += 1
                                add_log(f"  Created {pmid}_FULL_CONTEXT.md")

                                # Run scout
                                add_log(f"Running scout for PMID {pmid}...")
                                scout_result = run_data_scout(
                                    pmid, target_dir, gene_symbol
                                )
                                if scout_result:
                                    scouted += 1
                                    add_log(f"  Created {pmid}_DATA_ZONES.md")
                        except Exception as e:
                            add_log(f"  Error processing {pmid}: {e}")

                add_log(
                    f"Post-processing complete: {converted} converted, {scouted} scouted"
                )
            elif auto_process and not FETCH_MANAGER_AVAILABLE:
                add_log(
                    "Warning: Post-processing skipped - fetch_manager not available"
                )
            elif auto_process and not gene_symbol:
                add_log("Warning: Post-processing skipped - gene_symbol not provided")

            _browser_fetch_tasks[task_id].status = "completed"
        else:
            _browser_fetch_tasks[task_id].status = "failed"
            _browser_fetch_tasks[
                task_id
            ].error = f"Process exited with code {return_code}"
            add_log(f"Browser fetch failed with exit code {return_code}")

    except subprocess.TimeoutExpired:
        _browser_fetch_tasks[task_id].status = "failed"
        _browser_fetch_tasks[task_id].error = "Browser fetch timed out after 1 hour"
        add_log("ERROR: Browser fetch timed out after 1 hour")
        if "process" in locals():
            process.kill()
    except Exception as e:
        _browser_fetch_tasks[task_id].status = "failed"
        _browser_fetch_tasks[task_id].error = str(e)
        add_log(f"ERROR: {e}")


@app.post("/api/browser-fetch/start", response_model=BrowserFetchStatus)
async def start_browser_fetch(
    request: BrowserFetchRequest, background_tasks: BackgroundTasks
):
    """
    Start browser-based fetching of paywalled papers.

    This runs browser_fetch.py in the background to download papers
    that couldn't be fetched automatically.
    """
    import uuid

    csv_path = Path(request.csv_path).expanduser().resolve()

    if not csv_path.exists():
        raise HTTPException(status_code=404, detail=f"CSV file not found: {csv_path}")

    if not csv_path.name.endswith(".csv"):
        raise HTTPException(status_code=400, detail="File must be a CSV")

    # Create task
    task_id = str(uuid.uuid4())[:8]
    output_dir = str(csv_path.parent)

    _browser_fetch_tasks[task_id] = BrowserFetchStatus(
        task_id=task_id,
        status="starting",
        output_dir=output_dir,
    )

    # Start in background
    background_tasks.add_task(
        lambda: threading.Thread(
            target=_run_browser_fetch,
            args=(
                task_id,
                str(csv_path),
                request.headless,
                request.max_papers,
                request.use_claude,
                request.wait_for_captcha,
                request.retry_failures_only,
                request.gene_symbol,
                request.auto_process,
                request.fallback_to_manual,
            ),
            daemon=True,
        ).start()
    )

    return _browser_fetch_tasks[task_id]


@app.get("/api/browser-fetch/status/{task_id}", response_model=BrowserFetchStatus)
async def get_browser_fetch_status(task_id: str):
    """Get the status of a browser fetch task."""
    if task_id not in _browser_fetch_tasks:
        raise HTTPException(status_code=404, detail=f"Task not found: {task_id}")

    return _browser_fetch_tasks[task_id]


@app.get("/api/browser-fetch/tasks")
async def list_browser_fetch_tasks():
    """List all browser fetch tasks."""
    return {"tasks": list(_browser_fetch_tasks.values())}


@app.get("/api/browser-fetch/logs/{task_id}")
async def get_browser_fetch_logs(task_id: str, tail: int = 100):
    """
    Get full logs for a browser fetch task.

    Args:
        task_id: The task ID
        tail: Number of lines from the end (default 100, 0 for all)
    """
    if task_id not in _browser_fetch_logs:
        raise HTTPException(
            status_code=404, detail=f"Logs not found for task: {task_id}"
        )

    logs = _browser_fetch_logs[task_id]
    if tail > 0 and len(logs) > tail:
        logs = logs[-tail:]

    return {
        "task_id": task_id,
        "logs": logs,
        "total_lines": len(_browser_fetch_logs[task_id]),
    }


# =============================================================================
# WebSocket for Real-time Updates
# =============================================================================


@app.websocket("/ws/{job_id}")
async def websocket_endpoint(websocket: WebSocket, job_id: str):
    """WebSocket endpoint for real-time job progress."""
    await connection_manager.connect(websocket, job_id)

    # Send current state
    checkpoint = checkpoint_manager.load(job_id)
    if checkpoint:
        await websocket.send_json(
            {
                "type": "state",
                "job": checkpoint_to_response(checkpoint).model_dump(),
                "logs": checkpoint.log_lines[-50:],
            }
        )

    try:
        while True:
            # Keep connection alive and handle client messages
            data = await websocket.receive_text()
            if data == "ping":
                await websocket.send_json({"type": "pong"})
    except WebSocketDisconnect:
        connection_manager.disconnect(websocket, job_id)


# =============================================================================
# Main Entry Point
# =============================================================================


def main():
    """Run the server."""
    import argparse
    import uvicorn

    parser = argparse.ArgumentParser(description="GeneVariantFetcher GUI Server")
    parser.add_argument("--host", default="localhost", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload")

    args = parser.parse_args()

    print(f"\n{'='*60}")
    print("GeneVariantFetcher GUI Server")
    print(f"{'='*60}")
    print(f"Starting server at http://{args.host}:{args.port}")
    print("Open this URL in your browser to access the GUI")
    print(f"{'='*60}\n")

    uvicorn.run(
        "gui.server:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
