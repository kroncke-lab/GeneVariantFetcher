"""
Checkpoint system for resumable pipeline jobs.

Saves job state to disk after each major pipeline step, allowing resume
after interruption (computer shutdown, crash, etc.).
"""

import json
import logging
from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional
import fcntl
import os

logger = logging.getLogger(__name__)


class PipelineStep(str, Enum):
    """Pipeline steps in execution order."""
    PENDING = "pending"
    DISCOVERING_SYNONYMS = "discovering_synonyms"
    FETCHING_PMIDS = "fetching_pmids"
    FETCHING_ABSTRACTS = "fetching_abstracts"
    FILTERING_PAPERS = "filtering_papers"
    DOWNLOADING_FULLTEXT = "downloading_fulltext"
    SCOUTING_DATA = "scouting_data"  # New step for running Data Scout
    EXTRACTING_VARIANTS = "extracting_variants"
    AGGREGATING_DATA = "aggregating_data"
    MIGRATING_DATABASE = "migrating_database"
    COMPLETED = "completed"
    FAILED = "failed"

    @classmethod
    def get_display_name(cls, step: "PipelineStep") -> str:
        """Get human-readable name for step."""
        names = {
            cls.PENDING: "Pending",
            cls.DISCOVERING_SYNONYMS: "Discovering Gene Synonyms",
            cls.FETCHING_PMIDS: "Fetching PMIDs",
            cls.FETCHING_ABSTRACTS: "Fetching Abstracts",
            cls.FILTERING_PAPERS: "Filtering Papers",
            cls.DOWNLOADING_FULLTEXT: "Downloading Full-Text",
            cls.SCOUTING_DATA: "Scouting Data Zones",
            cls.EXTRACTING_VARIANTS: "Extracting Variants",
            cls.AGGREGATING_DATA: "Aggregating Data",
            cls.MIGRATING_DATABASE: "Creating Database",
            cls.COMPLETED: "Completed",
            cls.FAILED: "Failed",
        }
        return names.get(step, str(step))


@dataclass
class JobCheckpoint:
    """Persistent job state that survives restarts."""

    # Job identification
    job_id: str
    gene_symbol: str
    email: str
    output_dir: str

    # Configuration
    max_pmids: int = 100
    max_papers_to_download: int = 50
    tier_threshold: int = 1
    use_clinical_triage: bool = False
    auto_synonyms: bool = False
    synonyms: List[str] = field(default_factory=list)

    # Advanced settings
    enable_tier1: bool = True
    enable_tier2: bool = True
    use_pubmind: bool = True
    use_pubmed: bool = True
    use_europepmc: bool = False
    tier2_confidence_threshold: float = 0.5
    scout_enabled: bool = True
    scout_min_relevance: float = 0.3

    # Folder job mode (skip discovery/download, start from extraction)
    is_folder_job: bool = False
    folder_path: Optional[str] = None
    run_scout_on_folder: bool = False
    skip_already_extracted: bool = True

    # State tracking
    current_step: PipelineStep = PipelineStep.PENDING
    step_progress: Dict[str, Any] = field(default_factory=dict)

    # Timestamps
    created_at: str = ""
    updated_at: str = ""
    started_at: Optional[str] = None
    completed_at: Optional[str] = None

    # Results from completed steps (for resume)
    discovered_pmids: List[str] = field(default_factory=list)
    filtered_pmids: List[str] = field(default_factory=list)
    downloaded_pmids: List[str] = field(default_factory=list)
    extracted_pmids: List[str] = field(default_factory=list)

    # Error tracking
    error_message: Optional[str] = None
    error_step: Optional[str] = None

    # Log buffer for UI
    log_lines: List[str] = field(default_factory=list)

    def __post_init__(self):
        if not self.created_at:
            self.created_at = datetime.now().isoformat()
        self.updated_at = datetime.now().isoformat()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        data = asdict(self)
        data["current_step"] = self.current_step.value
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "JobCheckpoint":
        """Create from dictionary."""
        data = data.copy()
        if "current_step" in data:
            data["current_step"] = PipelineStep(data["current_step"])
        return cls(**data)

    def update_step(self, step: PipelineStep, progress: Optional[Dict[str, Any]] = None):
        """Update current step and optional progress info."""
        self.current_step = step
        self.updated_at = datetime.now().isoformat()
        if progress:
            self.step_progress[step.value] = progress

    def add_log(self, message: str):
        """Add a log line (keeps last 500 lines)."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_lines.append(f"[{timestamp}] {message}")
        if len(self.log_lines) > 500:
            self.log_lines = self.log_lines[-500:]

    def mark_failed(self, error: str, step: Optional[PipelineStep] = None):
        """Mark job as failed with error info."""
        self.current_step = PipelineStep.FAILED
        self.error_message = error
        self.error_step = step.value if step else None
        self.updated_at = datetime.now().isoformat()

    def mark_completed(self):
        """Mark job as completed."""
        self.current_step = PipelineStep.COMPLETED
        self.completed_at = datetime.now().isoformat()
        self.updated_at = datetime.now().isoformat()

    @property
    def is_resumable(self) -> bool:
        """Check if job can be resumed."""
        return self.current_step not in (
            PipelineStep.PENDING,
            PipelineStep.COMPLETED,
            PipelineStep.FAILED,
        )

    @property
    def output_path(self) -> Path:
        """Get the output path for this job."""
        # Use timestamp from created_at for consistent directory naming
        timestamp = datetime.fromisoformat(self.created_at).strftime("%Y%m%d_%H%M%S")
        return Path(self.output_dir) / self.gene_symbol / timestamp


class CheckpointManager:
    """Manages job checkpoints on disk."""

    CHECKPOINT_DIR = ".gvf_jobs"
    CHECKPOINT_FILE = "checkpoint.json"

    def __init__(self, base_dir: Optional[Path] = None):
        """
        Initialize checkpoint manager.

        Args:
            base_dir: Base directory for storing checkpoints.
                      Defaults to ~/.gvf_jobs
        """
        if base_dir is None:
            base_dir = Path.home() / self.CHECKPOINT_DIR
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def _job_dir(self, job_id: str) -> Path:
        """Get directory for a specific job."""
        return self.base_dir / job_id

    def _checkpoint_path(self, job_id: str) -> Path:
        """Get checkpoint file path for a job."""
        return self._job_dir(job_id) / self.CHECKPOINT_FILE

    def save(self, checkpoint: JobCheckpoint) -> None:
        """
        Save checkpoint to disk atomically.

        Uses file locking and atomic write to prevent corruption.
        """
        job_dir = self._job_dir(checkpoint.job_id)
        job_dir.mkdir(parents=True, exist_ok=True)

        checkpoint_path = self._checkpoint_path(checkpoint.job_id)
        temp_path = checkpoint_path.with_suffix(".tmp")

        checkpoint.updated_at = datetime.now().isoformat()

        try:
            # Write to temp file first
            with open(temp_path, "w", encoding="utf-8") as f:
                # Get exclusive lock
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                try:
                    json.dump(checkpoint.to_dict(), f, indent=2)
                    f.flush()
                    os.fsync(f.fileno())
                finally:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

            # Atomic rename
            temp_path.rename(checkpoint_path)
            logger.debug(f"Saved checkpoint for job {checkpoint.job_id}")

        except Exception as e:
            logger.error(f"Failed to save checkpoint: {e}")
            if temp_path.exists():
                temp_path.unlink()
            raise

    def load(self, job_id: str) -> Optional[JobCheckpoint]:
        """Load checkpoint from disk."""
        checkpoint_path = self._checkpoint_path(job_id)

        if not checkpoint_path.exists():
            return None

        try:
            with open(checkpoint_path, "r", encoding="utf-8") as f:
                fcntl.flock(f.fileno(), fcntl.LOCK_SH)
                try:
                    data = json.load(f)
                finally:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

            return JobCheckpoint.from_dict(data)

        except Exception as e:
            logger.error(f"Failed to load checkpoint for {job_id}: {e}")
            return None

    def delete(self, job_id: str) -> bool:
        """Delete a job checkpoint."""
        job_dir = self._job_dir(job_id)
        if job_dir.exists():
            import shutil
            shutil.rmtree(job_dir)
            return True
        return False

    def list_jobs(self) -> List[JobCheckpoint]:
        """List all saved job checkpoints."""
        jobs = []
        if not self.base_dir.exists():
            return jobs

        for job_dir in self.base_dir.iterdir():
            if job_dir.is_dir():
                checkpoint = self.load(job_dir.name)
                if checkpoint:
                    jobs.append(checkpoint)

        # Sort by creation time, newest first
        jobs.sort(key=lambda j: j.created_at, reverse=True)
        return jobs

    def get_incomplete_jobs(self) -> List[JobCheckpoint]:
        """Get jobs that were interrupted and can be resumed."""
        return [j for j in self.list_jobs() if j.is_resumable]

    def get_active_job(self) -> Optional[JobCheckpoint]:
        """Get currently running job (if any)."""
        for job in self.list_jobs():
            if job.current_step not in (
                PipelineStep.PENDING,
                PipelineStep.COMPLETED,
                PipelineStep.FAILED,
            ):
                return job
        return None

    def cleanup_old_jobs(self, max_age_days: int = 30) -> int:
        """Remove completed/failed jobs older than max_age_days."""
        from datetime import timedelta

        cutoff = datetime.now() - timedelta(days=max_age_days)
        removed = 0

        for job in self.list_jobs():
            if job.current_step in (PipelineStep.COMPLETED, PipelineStep.FAILED):
                job_date = datetime.fromisoformat(job.updated_at)
                if job_date < cutoff:
                    self.delete(job.job_id)
                    removed += 1

        return removed
