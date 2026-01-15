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

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from gui.checkpoint import CheckpointManager, JobCheckpoint, PipelineStep
from gui.worker import PipelineWorker, ProgressCallback, create_job

# Path to .env file
ENV_FILE_PATH = Path(__file__).resolve().parent.parent / ".env"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# =============================================================================
# Request/Response Models
# =============================================================================


class JobCreateRequest(BaseModel):
    """Request to create a new pipeline job."""

    gene_symbol: str = Field(..., description="Gene symbol (e.g., BRCA1, SCN5A)")
    email: str = Field(..., description="Email for NCBI E-utilities")
    output_dir: str = Field(..., description="Output directory for results")

    # Pipeline limits
    max_pmids: int = Field(100, ge=1, le=20000, description="Maximum PMIDs to fetch")
    max_papers_to_download: int = Field(50, ge=1, le=500, description="Maximum papers to download")

    # Synonym settings
    auto_synonyms: bool = Field(False, description="Auto-discover gene synonyms")
    synonyms: List[str] = Field(default_factory=list, description="Manual synonyms")

    # Specific PMIDs for testing (skip discovery if provided)
    specific_pmids: List[str] = Field(default_factory=list, description="Specific PMIDs for testing (skips discovery)")

    # Advanced settings
    tier_threshold: int = Field(1, ge=0, le=10, description="Model cascade threshold")
    use_clinical_triage: bool = Field(False, description="Use clinical triage filter")
    enable_tier1: bool = Field(True, description="Enable keyword filtering")
    enable_tier2: bool = Field(True, description="Enable LLM filtering")
    use_pubmind: bool = Field(True, description="Use PubMind source")
    use_pubmed: bool = Field(True, description="Use PubMed source")
    use_europepmc: bool = Field(False, description="Use Europe PMC source")
    tier2_confidence_threshold: float = Field(0.5, ge=0.0, le=1.0)
    scout_enabled: bool = Field(True, description="Enable DATA_ZONES condensation")


class JobResponse(BaseModel):
    """Response with job information."""

    job_id: str
    gene_symbol: str
    current_step: str
    current_step_display: str
    created_at: str
    updated_at: str
    started_at: Optional[str]
    completed_at: Optional[str]
    error_message: Optional[str]
    is_resumable: bool
    output_dir: str
    statistics: Dict[str, Any] = Field(default_factory=dict)


class JobListResponse(BaseModel):
    """Response with list of jobs."""

    jobs: List[JobResponse]
    incomplete_count: int


class EnvSettingsResponse(BaseModel):
    """Response with environment settings."""

    settings: Dict[str, Optional[str]]
    env_file_exists: bool
    env_file_path: str


class EnvSettingsUpdateRequest(BaseModel):
    """Request to update environment settings."""

    settings: Dict[str, str]


class DirectoryEntry(BaseModel):
    """A directory entry for the file browser."""

    name: str
    path: str
    is_dir: bool
    is_readable: bool = True


class DirectoryListResponse(BaseModel):
    """Response with directory contents."""

    current_path: str
    parent_path: Optional[str]
    entries: List[DirectoryEntry]
    can_create: bool = True


class FolderAnalysis(BaseModel):
    """Analysis of an existing folder with downloaded papers."""

    path: str
    is_valid: bool
    gene_symbol: Optional[str] = None
    detected_gene_symbols: List[str] = Field(default_factory=list)

    # File counts
    full_context_count: int = 0
    data_zones_count: int = 0
    extraction_count: int = 0

    # Lists of PMIDs found
    full_context_pmids: List[str] = Field(default_factory=list)
    data_zones_pmids: List[str] = Field(default_factory=list)
    extraction_pmids: List[str] = Field(default_factory=list)

    # Status flags
    has_fulltext: bool = False
    has_scouted: bool = False
    scouting_complete: bool = False
    has_extractions: bool = False

    # Messages
    status_message: str = ""
    recommendations: List[str] = Field(default_factory=list)


class FolderJobRequest(BaseModel):
    """Request to create a job from an existing folder."""

    folder_path: str = Field(..., description="Path to folder with downloaded papers")
    gene_symbol: str = Field(..., description="Gene symbol for extraction")
    email: str = Field(..., description="Email for NCBI (if needed)")

    # Options
    run_scout: bool = Field(False, description="Run Data Scout on unscouted files")
    skip_already_extracted: bool = Field(True, description="Skip PMIDs that already have extractions")

    # Advanced settings
    tier_threshold: int = Field(1, ge=0, le=10)
    scout_enabled: bool = Field(True)
    scout_min_relevance: float = Field(0.3, ge=0.0, le=1.0)


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
        self._broadcast({
            "type": "step_start",
            "step": step.value,
            "step_display": PipelineStep.get_display_name(step),
            "message": message,
            "timestamp": datetime.now().isoformat(),
        })

    def on_step_complete(self, step: PipelineStep, message: str, stats: Dict[str, Any]):
        self._broadcast({
            "type": "step_complete",
            "step": step.value,
            "step_display": PipelineStep.get_display_name(step),
            "message": message,
            "stats": stats,
            "timestamp": datetime.now().isoformat(),
        })

    def on_log(self, message: str):
        self._broadcast({
            "type": "log",
            "message": message,
            "timestamp": datetime.now().isoformat(),
        })

    def on_error(self, step: PipelineStep, error: str):
        self._broadcast({
            "type": "error",
            "step": step.value,
            "error": error,
            "timestamp": datetime.now().isoformat(),
        })


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
            connection_manager.broadcast(job_id, {
                "type": "job_complete",
                "status": final_checkpoint.current_step.value if final_checkpoint else "unknown",
                "timestamp": datetime.now().isoformat(),
            }),
            loop,
        )

    except Exception as e:
        logger.exception(f"Job {job_id} failed: {e}")
        asyncio.run_coroutine_threadsafe(
            connection_manager.broadcast(job_id, {
                "type": "error",
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            }),
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
            logger.warning(f"  - {job.job_id}: {job.gene_symbol} at step {job.current_step.value}")
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
        "extractions": [str(f) for f in (output_path / "extractions").glob("*.json")] if (output_path / "extractions").exists() else [],
        "pmid_status": [str(f) for f in (output_path / "pmid_status").glob("*.csv")] if (output_path / "pmid_status").exists() else [],
    }

    return results


# =============================================================================
# Settings Management API
# =============================================================================

# Define which settings are configurable via the GUI
CONFIGURABLE_SETTINGS = {
    # API Keys (required)
    "OPENAI_API_KEY": {"label": "OpenAI API Key", "type": "password", "required": True, "group": "API Keys"},
    "ANTHROPIC_API_KEY": {"label": "Anthropic API Key", "type": "password", "required": False, "group": "API Keys"},
    "NCBI_EMAIL": {"label": "NCBI Email", "type": "email", "required": True, "group": "API Keys"},
    "NCBI_API_KEY": {"label": "NCBI API Key (optional)", "type": "password", "required": False, "group": "API Keys"},

    # Model Configuration
    "TIER2_MODEL": {"label": "Tier 2 Model", "type": "text", "default": "gpt-4o-mini", "group": "Models"},
    "TIER3_MODELS": {"label": "Tier 3 Models (comma-separated)", "type": "text", "default": "gpt-4o-mini,gpt-4o", "group": "Models"},

    # Pipeline Defaults
    "ENABLE_TIER1": {"label": "Enable Tier 1 (Keyword Filter)", "type": "checkbox", "default": "true", "group": "Pipeline"},
    "ENABLE_TIER2": {"label": "Enable Tier 2 (LLM Filter)", "type": "checkbox", "default": "true", "group": "Pipeline"},
    "TIER2_CONFIDENCE_THRESHOLD": {"label": "Filter Confidence Threshold", "type": "number", "default": "0.5", "group": "Pipeline"},

    # Literature Sources
    "USE_PUBMIND": {"label": "Use PubMind", "type": "checkbox", "default": "true", "group": "Sources"},
    "USE_PUBMED": {"label": "Use PubMed", "type": "checkbox", "default": "true", "group": "Sources"},
    "USE_EUROPEPMC": {"label": "Use Europe PMC", "type": "checkbox", "default": "false", "group": "Sources"},

    # Scout Configuration
    "SCOUT_ENABLED": {"label": "Enable DATA_ZONES (faster extraction)", "type": "checkbox", "default": "true", "group": "Scout"},
    "SCOUT_MIN_RELEVANCE": {"label": "Scout Min Relevance", "type": "number", "default": "0.3", "group": "Scout"},

    # Output Defaults
    "DEFAULT_OUTPUT_DIR": {"label": "Default Output Directory", "type": "directory", "default": "./output", "group": "Output"},
}


def parse_env_file(path: Path) -> Dict[str, str]:
    """Parse a .env file and return key-value pairs."""
    settings = {}
    if not path.exists():
        return settings

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue
            # Parse key=value
            if "=" in line:
                key, _, value = line.partition("=")
                key = key.strip()
                value = value.strip()
                # Remove quotes if present
                if (value.startswith('"') and value.endswith('"')) or \
                   (value.startswith("'") and value.endswith("'")):
                    value = value[1:-1]
                settings[key] = value
    return settings


def write_env_file(path: Path, settings: Dict[str, str]):
    """Write settings to a .env file, preserving comments and structure."""
    existing_lines = []
    existing_keys = set()

    # Read existing file to preserve structure and comments
    if path.exists():
        with open(path, "r") as f:
            for line in f:
                original_line = line
                line = line.strip()
                if not line or line.startswith("#"):
                    existing_lines.append(original_line)
                elif "=" in line:
                    key, _, _ = line.partition("=")
                    key = key.strip()
                    existing_keys.add(key)
                    if key in settings:
                        # Update existing key
                        value = settings[key]
                        # Quote value if it contains spaces or special chars
                        if " " in value or "#" in value:
                            value = f'"{value}"'
                        existing_lines.append(f"{key}={value}\n")
                    else:
                        # Keep original line
                        existing_lines.append(original_line)

    # Add any new keys that weren't in the file
    new_keys = set(settings.keys()) - existing_keys
    if new_keys:
        if existing_lines and not existing_lines[-1].endswith("\n"):
            existing_lines.append("\n")
        existing_lines.append("\n# Added via GUI\n")
        for key in sorted(new_keys):
            value = settings[key]
            if " " in value or "#" in value:
                value = f'"{value}"'
            existing_lines.append(f"{key}={value}\n")

    # Write file
    with open(path, "w") as f:
        f.writelines(existing_lines)


@app.get("/api/settings", response_model=EnvSettingsResponse)
async def get_settings():
    """Get current environment settings."""
    # Load from .env file
    file_settings = parse_env_file(ENV_FILE_PATH)

    # Also check current environment (might have been set externally)
    settings = {}
    for key, config in CONFIGURABLE_SETTINGS.items():
        # Priority: .env file > environment variable > default
        if key in file_settings:
            settings[key] = file_settings[key]
        elif key in os.environ:
            settings[key] = os.environ[key]
        else:
            settings[key] = config.get("default", "")

    return EnvSettingsResponse(
        settings=settings,
        env_file_exists=ENV_FILE_PATH.exists(),
        env_file_path=str(ENV_FILE_PATH),
    )


@app.get("/api/settings/schema")
async def get_settings_schema():
    """Get the schema for configurable settings."""
    return {"settings": CONFIGURABLE_SETTINGS}


@app.post("/api/settings")
async def update_settings(request: EnvSettingsUpdateRequest):
    """Update environment settings."""
    # Validate that only known settings are being updated
    unknown = set(request.settings.keys()) - set(CONFIGURABLE_SETTINGS.keys())
    if unknown:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown settings: {', '.join(unknown)}"
        )

    # Load existing settings
    existing = parse_env_file(ENV_FILE_PATH)

    # Merge with new settings (only update non-empty values)
    for key, value in request.settings.items():
        if value or key not in existing:
            existing[key] = value

    # Write back to file
    write_env_file(ENV_FILE_PATH, existing)

    # Reload dotenv to update current process
    from dotenv import load_dotenv
    load_dotenv(ENV_FILE_PATH, override=True)

    logger.info(f"Settings updated: {list(request.settings.keys())}")

    return {"status": "ok", "updated": list(request.settings.keys())}


@app.get("/api/settings/validate")
async def validate_settings():
    """Validate current settings and return any issues."""
    file_settings = parse_env_file(ENV_FILE_PATH)

    issues = []
    warnings = []

    # Check required settings
    if not file_settings.get("OPENAI_API_KEY") and not os.environ.get("OPENAI_API_KEY"):
        issues.append("OPENAI_API_KEY is required for variant extraction")

    if not file_settings.get("NCBI_EMAIL") and not os.environ.get("NCBI_EMAIL"):
        issues.append("NCBI_EMAIL is required for literature fetching")

    # Check for placeholder values
    placeholders = {"your_api_key", "changeme", "your-openai-api-key"}
    for key, value in file_settings.items():
        if value.lower() in placeholders:
            issues.append(f"{key} appears to be a placeholder value")

    # Warnings
    if not file_settings.get("NCBI_API_KEY") and not os.environ.get("NCBI_API_KEY"):
        warnings.append("NCBI_API_KEY not set - API requests may be rate-limited")

    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
    }


# =============================================================================
# Folder Analysis API (for importing existing downloads)
# =============================================================================


def _extract_pmid_from_filename(filename: str) -> Optional[str]:
    """Extract PMID from various filename formats."""
    import re
    # Patterns: 12345678_FULL_CONTEXT.md, PMID_12345678_FULL_CONTEXT.md, etc.
    patterns = [
        r'^(\d{6,10})_(?:FULL_CONTEXT|DATA_ZONES)\.md$',
        r'^PMID_(\d{6,10})_(?:FULL_CONTEXT|DATA_ZONES)\.md$',
        r'^[A-Z]+_PMID_(\d{6,10})\.json$',  # SCN5A_PMID_12345678.json
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
    match = re.match(r'^([A-Z0-9]+)_PMID_\d+\.json$', filename, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return None


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
        folder_path == root or root in folder_path.parents
        for root in allowed_roots
    )
    if not is_allowed:
        raise HTTPException(status_code=403, detail="Access to this directory is not allowed")

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
            match = re.search(r'PMID_(\d+)', f.name)
            if match:
                extraction_pmids.append(match.group(1))

    # Also try to detect gene from folder name (e.g., output/SCN5A/20240101_120000)
    folder_gene = None
    if folder_path.parent.name and folder_path.parent.name.isupper() and len(folder_path.parent.name) <= 10:
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
        recommendations.append("No downloaded papers found. Run a full pipeline or download papers first.")
    elif not has_scouted:
        recommendations.append("Papers not yet scouted. Enable 'Run Data Scout' to create condensed DATA_ZONES files for faster extraction.")
    elif not scouting_complete:
        unscouted = set(full_context_pmids) - set(data_zones_pmids)
        recommendations.append(f"{len(unscouted)} papers have not been scouted yet. Enable 'Run Data Scout' to process them.")

    if has_fulltext and not has_extractions:
        recommendations.append("Ready for extraction. Start extraction to process papers.")
    elif has_extractions:
        unextracted = all_pmids - set(extraction_pmids)
        if unextracted:
            recommendations.append(f"{len(unextracted)} papers have not been extracted yet.")
        else:
            recommendations.append("All papers have been extracted. You can re-run extraction if needed.")

    # Build status message
    status_parts = []
    if full_context_pmids:
        status_parts.append(f"{len(full_context_pmids)} full-text papers")
    if data_zones_pmids:
        status_parts.append(f"{len(data_zones_pmids)} scouted")
    if extraction_pmids:
        status_parts.append(f"{len(extraction_pmids)} extracted")

    status_message = "Found: " + ", ".join(status_parts) if status_parts else "No papers found"

    return FolderAnalysis(
        path=str(folder_path),
        is_valid=has_fulltext,
        gene_symbol=folder_gene or (list(detected_genes)[0] if len(detected_genes) == 1 else None),
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
    )


@app.post("/api/jobs/from-folder", response_model=JobResponse)
async def create_job_from_folder(request: FolderJobRequest, background_tasks: BackgroundTasks):
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
            status_code=400,
            detail=f"Invalid folder: {analysis.status_message}"
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
        extraction_pmids=analysis.extraction_pmids if request.skip_already_extracted else [],
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
# Directory Browser API
# =============================================================================


@app.get("/api/browse")
async def browse_directory(path: str = "~") -> DirectoryListResponse:
    """Browse a directory and return its contents."""
    # Expand user home directory
    browse_path = Path(path).expanduser().resolve()

    # Security: Don't allow browsing outside of reasonable locations
    # Allow home directory and subdirectories, /tmp, and relative paths from cwd
    allowed_roots = [
        Path.home(),
        Path("/tmp"),
        Path.cwd(),
    ]

    is_allowed = any(
        browse_path == root or root in browse_path.parents
        for root in allowed_roots
    )

    if not is_allowed:
        raise HTTPException(
            status_code=403,
            detail="Access to this directory is not allowed"
        )

    if not browse_path.exists():
        raise HTTPException(status_code=404, detail="Directory not found")

    if not browse_path.is_dir():
        raise HTTPException(status_code=400, detail="Path is not a directory")

    entries = []
    try:
        for item in sorted(browse_path.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
            # Skip hidden files (optional)
            if item.name.startswith("."):
                continue

            try:
                is_readable = os.access(item, os.R_OK)
                entries.append(DirectoryEntry(
                    name=item.name,
                    path=str(item),
                    is_dir=item.is_dir(),
                    is_readable=is_readable,
                ))
            except PermissionError:
                entries.append(DirectoryEntry(
                    name=item.name,
                    path=str(item),
                    is_dir=item.is_dir(),
                    is_readable=False,
                ))
    except PermissionError:
        raise HTTPException(status_code=403, detail="Permission denied")

    # Calculate parent path
    parent_path = None
    if browse_path != Path.home() and browse_path.parent != browse_path:
        parent_path = str(browse_path.parent)

    # Check if we can create directories here
    can_create = os.access(browse_path, os.W_OK)

    return DirectoryListResponse(
        current_path=str(browse_path),
        parent_path=parent_path,
        entries=entries,
        can_create=can_create,
    )


@app.post("/api/browse/create")
async def create_directory(path: str):
    """Create a new directory."""
    new_path = Path(path).expanduser().resolve()

    # Security check
    allowed_roots = [Path.home(), Path("/tmp"), Path.cwd()]
    is_allowed = any(
        root in new_path.parents or new_path == root
        for root in allowed_roots
    )

    if not is_allowed:
        raise HTTPException(
            status_code=403,
            detail="Cannot create directory in this location"
        )

    if new_path.exists():
        raise HTTPException(status_code=400, detail="Path already exists")

    try:
        new_path.mkdir(parents=True, exist_ok=True)
        return {"status": "ok", "path": str(new_path)}
    except PermissionError:
        raise HTTPException(status_code=403, detail="Permission denied")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/browse/shortcuts")
async def get_directory_shortcuts():
    """Get common directory shortcuts."""
    shortcuts = []

    # Home directory
    home = Path.home()
    shortcuts.append({"name": "Home", "path": str(home), "icon": "home"})

    # Desktop (if exists)
    desktop = home / "Desktop"
    if desktop.exists():
        shortcuts.append({"name": "Desktop", "path": str(desktop), "icon": "desktop"})

    # Documents (if exists)
    documents = home / "Documents"
    if documents.exists():
        shortcuts.append({"name": "Documents", "path": str(documents), "icon": "folder"})

    # Current working directory
    cwd = Path.cwd()
    shortcuts.append({"name": "Project Root", "path": str(cwd), "icon": "code"})

    # Default output directory
    default_output = cwd / "output"
    shortcuts.append({"name": "Default Output", "path": str(default_output), "icon": "output"})

    return {"shortcuts": shortcuts}


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
        await websocket.send_json({
            "type": "state",
            "job": checkpoint_to_response(checkpoint).model_dump(),
            "logs": checkpoint.log_lines[-50:],
        })

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
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload")

    args = parser.parse_args()

    print(f"\n{'='*60}")
    print("GeneVariantFetcher GUI Server")
    print(f"{'='*60}")
    print(f"Starting server at http://{args.host}:{args.port}")
    print(f"Open this URL in your browser to access the GUI")
    print(f"{'='*60}\n")

    uvicorn.run(
        "gui.server:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
