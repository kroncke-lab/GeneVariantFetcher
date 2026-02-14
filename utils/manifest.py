"""
Manifest utilities for tracking pipeline stage results.

Provides structured tracking of processing outcomes across pipeline stages,
enabling resumability, audit trails, and stage handoffs.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Optional

SCHEMA_VERSION = "1.0"


class Status(str, Enum):
    """Processing outcome status codes."""

    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"
    PAYWALL = "PAYWALL"
    CAPTCHA = "CAPTCHA"
    TIMEOUT = "TIMEOUT"


class Stage(str, Enum):
    """Pipeline stages."""

    DOWNLOAD = "download"
    SCOUT = "scout"
    EXTRACT = "extract"
    MIGRATE = "migrate"
    REFETCH = "refetch"


@dataclass
class ManifestEntry:
    """A single processing result entry."""

    pmid: str
    status: Status
    error_message: Optional[str] = None
    files_created: list[str] = field(default_factory=list)
    timestamp: Optional[str] = None

    def __post_init__(self):
        # Auto-set timestamp if not provided
        if self.timestamp is None:
            self.timestamp = datetime.now(timezone.utc).isoformat()

        # Convert string status to enum if needed
        if isinstance(self.status, str):
            self.status = Status(self.status)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        d = {
            "pmid": self.pmid,
            "status": self.status.value,
            "timestamp": self.timestamp,
            "files_created": self.files_created,
        }
        if self.error_message:
            d["error_message"] = self.error_message
        return d

    @classmethod
    def from_dict(cls, data: dict) -> ManifestEntry:
        """Create entry from dictionary."""
        return cls(
            pmid=data["pmid"],
            status=Status(data["status"]),
            error_message=data.get("error_message"),
            files_created=data.get("files_created", []),
            timestamp=data.get("timestamp"),
        )


@dataclass
class Manifest:
    """
    Pipeline stage manifest for tracking processing results.

    Attributes:
        stage: Pipeline stage that produced this manifest
        gene: Gene symbol being processed (optional)
        created_at: ISO timestamp when manifest was created
        entries: List of processing result entries
        schema_version: Schema version for compatibility
    """

    stage: Stage
    gene: Optional[str] = None
    created_at: Optional[str] = None
    entries: list[ManifestEntry] = field(default_factory=list)
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self):
        # Auto-set created_at if not provided
        if self.created_at is None:
            self.created_at = datetime.now(timezone.utc).isoformat()

        # Convert string stage to enum if needed
        if isinstance(self.stage, str):
            self.stage = Stage(self.stage)

    def add_entry(self, entry: ManifestEntry) -> None:
        """Add a processing result entry."""
        self.entries.append(entry)

    def get_by_status(self, status: Status) -> list[ManifestEntry]:
        """Get all entries with a specific status."""
        return [e for e in self.entries if e.status == status]

    def get_successful(self) -> list[ManifestEntry]:
        """Get all successful entries."""
        return self.get_by_status(Status.SUCCESS)

    def get_failed(self) -> list[ManifestEntry]:
        """Get all non-successful entries."""
        return [e for e in self.entries if e.status != Status.SUCCESS]

    def get_pmids_by_status(self, status: Status) -> list[str]:
        """Get PMIDs with a specific status."""
        return [e.pmid for e in self.get_by_status(status)]

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        d = {
            "schema_version": self.schema_version,
            "stage": self.stage.value,
            "created_at": self.created_at,
            "entries": [e.to_dict() for e in self.entries],
        }
        if self.gene:
            d["gene"] = self.gene
        return d

    def save(self, path: str | Path) -> None:
        """
        Save manifest to file atomically.

        Uses write-to-temp-then-rename pattern to prevent corruption
        from interrupted writes.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        # Write to temp file in same directory (for atomic rename)
        fd, tmp_path = tempfile.mkstemp(
            dir=path.parent, prefix=".manifest_", suffix=".tmp"
        )
        try:
            with os.fdopen(fd, "w") as f:
                json.dump(self.to_dict(), f, indent=2)
            # Atomic rename
            os.replace(tmp_path, path)
        except Exception:
            # Clean up temp file on failure
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            raise

    @classmethod
    def load(cls, path: str | Path) -> Manifest:
        """Load manifest from file."""
        path = Path(path)
        with open(path) as f:
            data = json.load(f)

        return cls(
            schema_version=data.get("schema_version", SCHEMA_VERSION),
            stage=Stage(data["stage"]),
            gene=data.get("gene"),
            created_at=data.get("created_at"),
            entries=[ManifestEntry.from_dict(e) for e in data.get("entries", [])],
        )

    def __len__(self) -> int:
        """Return number of entries."""
        return len(self.entries)

    def summary(self) -> dict[str, int]:
        """Get count of entries by status."""
        counts = {}
        for entry in self.entries:
            status = entry.status.value
            counts[status] = counts.get(status, 0) + 1
        return counts
