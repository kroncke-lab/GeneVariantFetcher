"""Tests for manifest utilities."""

import json
import os
import tempfile
from pathlib import Path

import pytest

from utils.manifest import (
    Manifest,
    ManifestEntry,
    Status,
    Stage,
    SCHEMA_VERSION,
)


class TestManifestEntry:
    """Tests for ManifestEntry dataclass."""

    def test_create_entry_basic(self):
        """Test creating a basic entry."""
        entry = ManifestEntry(pmid="12345678", status=Status.SUCCESS)
        assert entry.pmid == "12345678"
        assert entry.status == Status.SUCCESS
        assert entry.timestamp is not None
        assert entry.error_message is None
        assert entry.files_created == []

    def test_create_entry_with_all_fields(self):
        """Test creating an entry with all fields."""
        entry = ManifestEntry(
            pmid="12345678",
            status=Status.FAILED,
            error_message="Connection refused",
            files_created=["path/to/file.pdf"],
            timestamp="2025-01-30T14:30:00Z",
        )
        assert entry.pmid == "12345678"
        assert entry.status == Status.FAILED
        assert entry.error_message == "Connection refused"
        assert entry.files_created == ["path/to/file.pdf"]
        assert entry.timestamp == "2025-01-30T14:30:00Z"

    def test_status_from_string(self):
        """Test that string status is converted to enum."""
        entry = ManifestEntry(pmid="12345678", status="SUCCESS")
        assert entry.status == Status.SUCCESS
        assert isinstance(entry.status, Status)

    def test_to_dict(self):
        """Test dictionary serialization."""
        entry = ManifestEntry(
            pmid="12345678",
            status=Status.PAYWALL,
            error_message="Requires subscription",
            files_created=["meta.json"],
            timestamp="2025-01-30T14:30:00Z",
        )
        d = entry.to_dict()
        assert d["pmid"] == "12345678"
        assert d["status"] == "PAYWALL"
        assert d["error_message"] == "Requires subscription"
        assert d["files_created"] == ["meta.json"]

    def test_to_dict_omits_empty_error(self):
        """Test that None error_message is omitted from dict."""
        entry = ManifestEntry(pmid="12345678", status=Status.SUCCESS)
        d = entry.to_dict()
        assert "error_message" not in d

    def test_from_dict(self):
        """Test creating entry from dictionary."""
        data = {
            "pmid": "87654321",
            "status": "TIMEOUT",
            "error_message": "Timed out",
            "files_created": ["a.pdf", "b.txt"],
            "timestamp": "2025-01-30T15:00:00Z",
        }
        entry = ManifestEntry.from_dict(data)
        assert entry.pmid == "87654321"
        assert entry.status == Status.TIMEOUT
        assert entry.error_message == "Timed out"
        assert entry.files_created == ["a.pdf", "b.txt"]


class TestManifest:
    """Tests for Manifest class."""

    def test_create_manifest_basic(self):
        """Test creating a basic manifest."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        assert manifest.stage == Stage.DOWNLOAD
        assert manifest.gene is None
        assert manifest.created_at is not None
        assert manifest.entries == []
        assert manifest.schema_version == SCHEMA_VERSION

    def test_create_manifest_with_gene(self):
        """Test creating manifest with gene specified."""
        manifest = Manifest(stage=Stage.SCOUT, gene="BRCA1")
        assert manifest.gene == "BRCA1"
        assert manifest.stage == Stage.SCOUT

    def test_stage_from_string(self):
        """Test that string stage is converted to enum."""
        manifest = Manifest(stage="extract")
        assert manifest.stage == Stage.EXTRACT
        assert isinstance(manifest.stage, Stage)

    def test_add_entry(self):
        """Test adding entries."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        entry = ManifestEntry(pmid="12345678", status=Status.SUCCESS)
        manifest.add_entry(entry)
        assert len(manifest) == 1
        assert manifest.entries[0].pmid == "12345678"

    def test_get_by_status(self):
        """Test filtering by status."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        manifest.add_entry(ManifestEntry(pmid="1", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="2", status=Status.FAILED))
        manifest.add_entry(ManifestEntry(pmid="3", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="4", status=Status.PAYWALL))

        success = manifest.get_by_status(Status.SUCCESS)
        assert len(success) == 2
        assert {e.pmid for e in success} == {"1", "3"}

    def test_get_successful(self):
        """Test getting successful entries."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        manifest.add_entry(ManifestEntry(pmid="1", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="2", status=Status.FAILED))

        success = manifest.get_successful()
        assert len(success) == 1
        assert success[0].pmid == "1"

    def test_get_failed(self):
        """Test getting non-successful entries."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        manifest.add_entry(ManifestEntry(pmid="1", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="2", status=Status.FAILED))
        manifest.add_entry(ManifestEntry(pmid="3", status=Status.PAYWALL))

        failed = manifest.get_failed()
        assert len(failed) == 2
        assert {e.pmid for e in failed} == {"2", "3"}

    def test_get_pmids_by_status(self):
        """Test getting PMIDs by status."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        manifest.add_entry(ManifestEntry(pmid="111", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="222", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="333", status=Status.TIMEOUT))

        pmids = manifest.get_pmids_by_status(Status.SUCCESS)
        assert pmids == ["111", "222"]

    def test_summary(self):
        """Test summary counts."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        manifest.add_entry(ManifestEntry(pmid="1", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="2", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="3", status=Status.FAILED))
        manifest.add_entry(ManifestEntry(pmid="4", status=Status.PAYWALL))

        summary = manifest.summary()
        assert summary == {"SUCCESS": 2, "FAILED": 1, "PAYWALL": 1}

    def test_to_dict(self):
        """Test dictionary serialization."""
        manifest = Manifest(
            stage=Stage.EXTRACT,
            gene="TP53",
            created_at="2025-01-30T12:00:00Z",
        )
        manifest.add_entry(ManifestEntry(
            pmid="12345",
            status=Status.SUCCESS,
            timestamp="2025-01-30T12:01:00Z",
        ))

        d = manifest.to_dict()
        assert d["schema_version"] == SCHEMA_VERSION
        assert d["stage"] == "extract"
        assert d["gene"] == "TP53"
        assert d["created_at"] == "2025-01-30T12:00:00Z"
        assert len(d["entries"]) == 1

    def test_to_dict_omits_empty_gene(self):
        """Test that None gene is omitted from dict."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        d = manifest.to_dict()
        assert "gene" not in d


class TestManifestPersistence:
    """Tests for save/load functionality."""

    def test_save_and_load(self, tmp_path):
        """Test round-trip save and load."""
        manifest = Manifest(stage=Stage.DOWNLOAD, gene="BRCA1")
        manifest.add_entry(ManifestEntry(
            pmid="12345678",
            status=Status.SUCCESS,
            files_created=["file1.pdf", "file2.json"],
        ))
        manifest.add_entry(ManifestEntry(
            pmid="87654321",
            status=Status.PAYWALL,
            error_message="Subscription required",
        ))

        path = tmp_path / "manifest.json"
        manifest.save(path)

        # Verify file exists and is valid JSON
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert data["stage"] == "download"
        assert data["gene"] == "BRCA1"
        assert len(data["entries"]) == 2

        # Load and verify
        loaded = Manifest.load(path)
        assert loaded.stage == Stage.DOWNLOAD
        assert loaded.gene == "BRCA1"
        assert len(loaded) == 2
        assert loaded.entries[0].pmid == "12345678"
        assert loaded.entries[0].status == Status.SUCCESS
        assert loaded.entries[1].status == Status.PAYWALL
        assert loaded.entries[1].error_message == "Subscription required"

    def test_save_creates_parent_dirs(self, tmp_path):
        """Test that save creates parent directories."""
        manifest = Manifest(stage=Stage.SCOUT)
        path = tmp_path / "nested" / "dirs" / "manifest.json"
        manifest.save(path)
        assert path.exists()

    def test_atomic_write(self, tmp_path):
        """Test that save uses atomic write pattern."""
        manifest = Manifest(stage=Stage.DOWNLOAD)
        manifest.add_entry(ManifestEntry(pmid="123", status=Status.SUCCESS))
        
        path = tmp_path / "manifest.json"
        manifest.save(path)
        
        # No temp files should remain
        temp_files = list(tmp_path.glob(".manifest_*.tmp"))
        assert len(temp_files) == 0

    def test_overwrite_existing(self, tmp_path):
        """Test overwriting an existing manifest."""
        path = tmp_path / "manifest.json"
        
        # Create initial manifest
        m1 = Manifest(stage=Stage.DOWNLOAD, gene="GENE1")
        m1.add_entry(ManifestEntry(pmid="111", status=Status.SUCCESS))
        m1.save(path)
        
        # Overwrite with new manifest
        m2 = Manifest(stage=Stage.SCOUT, gene="GENE2")
        m2.add_entry(ManifestEntry(pmid="222", status=Status.FAILED))
        m2.add_entry(ManifestEntry(pmid="333", status=Status.SUCCESS))
        m2.save(path)
        
        # Load and verify it's the new one
        loaded = Manifest.load(path)
        assert loaded.stage == Stage.SCOUT
        assert loaded.gene == "GENE2"
        assert len(loaded) == 2


class TestStatusEnum:
    """Tests for Status enum."""

    def test_all_statuses_exist(self):
        """Verify all expected status codes exist."""
        expected = {"SUCCESS", "FAILED", "SKIPPED", "PAYWALL", "CAPTCHA", "TIMEOUT"}
        actual = {s.value for s in Status}
        assert actual == expected

    def test_status_is_string(self):
        """Verify Status values are strings for JSON compatibility."""
        for status in Status:
            assert isinstance(status.value, str)


class TestStageEnum:
    """Tests for Stage enum."""

    def test_all_stages_exist(self):
        """Verify all expected stages exist."""
        expected = {"download", "scout", "extract", "migrate"}
        actual = {s.value for s in Stage}
        assert actual == expected
