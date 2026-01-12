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
    max_pmids: int = Field(100, ge=1, le=1000, description="Maximum PMIDs to fetch")
    max_papers_to_download: int = Field(50, ge=1, le=500, description="Maximum papers to download")

    # Synonym settings
    auto_synonyms: bool = Field(False, description="Auto-discover gene synonyms")
    synonyms: List[str] = Field(default_factory=list, description="Manual synonyms")

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
