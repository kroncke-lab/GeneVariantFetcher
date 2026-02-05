"""Tests for the extract CLI module."""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from cli.extract import (
    find_input_files,
    extract_pmid_from_filename,
    run_extraction,
)
from utils.manifest import Manifest, ManifestEntry, Status, Stage


class TestExtractPmidFromFilename:
    """Tests for extract_pmid_from_filename function."""

    def test_extract_from_data_zones(self):
        """Test extracting PMID from DATA_ZONES filename."""
        assert extract_pmid_from_filename("12345678_DATA_ZONES.md") == "12345678"
        assert extract_pmid_from_filename("PMC12345_DATA_ZONES.md") == "PMC12345"

    def test_extract_from_full_context(self):
        """Test extracting PMID from FULL_CONTEXT filename."""
        assert (
            extract_pmid_from_filename("12345678_FULL_CONTEXT.md", use_full_text=True)
            == "12345678"
        )
        assert (
            extract_pmid_from_filename("PMC12345_FULL_CONTEXT.md", use_full_text=True)
            == "PMC12345"
        )

    def test_wrong_suffix_returns_none(self):
        """Test that wrong suffix returns None."""
        assert (
            extract_pmid_from_filename("12345678_FULL_CONTEXT.md", use_full_text=False)
            is None
        )
        assert (
            extract_pmid_from_filename("12345678_DATA_ZONES.md", use_full_text=True)
            is None
        )


class TestFindInputFiles:
    """Tests for find_input_files function."""

    def test_find_files_from_manifest(self, tmp_path):
        """Test finding files using manifest entries."""
        # Create test files
        zones_file = tmp_path / "12345678_DATA_ZONES.md"
        zones_file.write_text("# Data Zones for 12345678\n\nSome content")

        # Create manifest with file reference
        manifest = Manifest(stage=Stage.SCOUT, gene="BRCA1")
        manifest.add_entry(
            ManifestEntry(
                pmid="12345678",
                status=Status.SUCCESS,
                files_created=[str(zones_file)],
            )
        )

        files = find_input_files(tmp_path, manifest)
        assert len(files) == 1
        assert files[0] == (zones_file, "12345678")

    def test_find_files_from_manifest_fallback(self, tmp_path):
        """Test finding files when manifest entry doesn't have files_created."""
        # Create test files
        zones_file = tmp_path / "12345678_DATA_ZONES.md"
        zones_file.write_text("# Data Zones")

        # Create manifest without files_created
        manifest = Manifest(stage=Stage.SCOUT, gene="BRCA1")
        manifest.add_entry(
            ManifestEntry(
                pmid="12345678",
                status=Status.SUCCESS,
                files_created=[],  # Empty
            )
        )

        files = find_input_files(tmp_path, manifest)
        assert len(files) == 1
        assert files[0][1] == "12345678"

    def test_find_files_glob_mode(self, tmp_path):
        """Test finding files by globbing when no manifest provided."""
        # Create test files
        (tmp_path / "11111111_DATA_ZONES.md").write_text("content1")
        (tmp_path / "22222222_DATA_ZONES.md").write_text("content2")
        (tmp_path / "other_file.txt").write_text("ignore me")

        files = find_input_files(tmp_path, manifest=None)
        assert len(files) == 2
        pmids = {f[1] for f in files}
        assert pmids == {"11111111", "22222222"}

    def test_find_full_text_files(self, tmp_path):
        """Test finding FULL_CONTEXT.md files."""
        # Create test files
        (tmp_path / "12345678_FULL_CONTEXT.md").write_text("full content")

        files = find_input_files(tmp_path, manifest=None, use_full_text=True)
        assert len(files) == 1
        assert files[0][1] == "12345678"

    def test_manifest_skips_failed_entries(self, tmp_path):
        """Test that failed manifest entries are skipped."""
        # Create file for both PMIDs
        (tmp_path / "11111111_DATA_ZONES.md").write_text("content1")
        (tmp_path / "22222222_DATA_ZONES.md").write_text("content2")

        manifest = Manifest(stage=Stage.SCOUT, gene="TEST")
        manifest.add_entry(ManifestEntry(pmid="11111111", status=Status.SUCCESS))
        manifest.add_entry(ManifestEntry(pmid="22222222", status=Status.FAILED))

        files = find_input_files(tmp_path, manifest)
        # Only successful entry should be returned
        assert len(files) == 1
        assert files[0][1] == "11111111"


class TestRunExtraction:
    """Integration tests for run_extraction function."""

    @patch("cli.extract.ExpertExtractor")
    def test_run_extraction_with_manifest(self, mock_extractor_class, tmp_path):
        """Test running extraction with a manifest input."""
        input_dir = tmp_path / "scout_output"
        input_dir.mkdir()
        output_dir = tmp_path / "extractions"

        # Create scout manifest
        scout_manifest = Manifest(stage=Stage.SCOUT, gene="BRCA1")
        scout_manifest.add_entry(
            ManifestEntry(
                pmid="12345678",
                status=Status.SUCCESS,
                files_created=[str(input_dir / "12345678_DATA_ZONES.md")],
            )
        )
        scout_manifest.save(input_dir / "scout_manifest.json")

        # Create the DATA_ZONES file
        zones_file = input_dir / "12345678_DATA_ZONES.md"
        zones_file.write_text(
            "# PMID 12345678\n## Data Zone 1\n| Variant | Effect |\n| p.R1234H | Pathogenic |"
        )

        # Mock the extractor
        mock_extractor = Mock()
        mock_result = Mock()
        mock_result.success = True
        mock_result.extracted_data = {
            "variants": [{"protein_notation": "p.R1234H"}],
            "extraction_metadata": {"total_variants_found": 1},
        }
        mock_result.model_used = "test-model"
        mock_result.error = None
        mock_extractor.extract.return_value = mock_result
        mock_extractor_class.return_value = mock_extractor

        # Run extraction
        manifest = run_extraction(
            input_path=input_dir / "scout_manifest.json",
            output_dir=output_dir,
            gene="BRCA1",
        )

        # Verify results
        assert len(manifest) == 1
        assert manifest.entries[0].status == Status.SUCCESS
        assert manifest.entries[0].pmid == "12345678"
        assert len(manifest.entries[0].files_created) > 0

        # Verify output files
        assert (output_dir / "extraction_manifest.json").exists()
        assert (output_dir / "extraction_summary.json").exists()
        assert (output_dir / "12345678_extraction.json").exists()

    @patch("cli.extract.ExpertExtractor")
    def test_run_extraction_directory_with_auto_manifest(
        self, mock_extractor_class, tmp_path
    ):
        """Test running extraction when directory contains scout_manifest.json."""
        input_dir = tmp_path / "scout_output"
        input_dir.mkdir()
        output_dir = tmp_path / "extractions"

        # Create scout manifest in directory
        scout_manifest = Manifest(stage=Stage.SCOUT, gene="TP53")
        scout_manifest.add_entry(
            ManifestEntry(
                pmid="99999999",
                status=Status.SUCCESS,
            )
        )
        scout_manifest.save(input_dir / "scout_manifest.json")

        # Create DATA_ZONES file
        (input_dir / "99999999_DATA_ZONES.md").write_text("# Test content")

        # Mock
        mock_extractor = Mock()
        mock_result = Mock()
        mock_result.success = True
        mock_result.extracted_data = {"variants": []}
        mock_result.model_used = "test"
        mock_result.error = None
        mock_extractor.extract.return_value = mock_result
        mock_extractor_class.return_value = mock_extractor

        # Run with directory (not manifest file directly)
        manifest = run_extraction(
            input_path=input_dir,  # Directory, not manifest
            output_dir=output_dir,
            gene="TP53",
        )

        # Should still find and use the manifest
        assert len(manifest) == 1
        assert manifest.entries[0].pmid == "99999999"

    def test_run_extraction_empty_directory(self, tmp_path):
        """Test running extraction on empty directory."""
        input_dir = tmp_path / "empty"
        input_dir.mkdir()
        output_dir = tmp_path / "output"

        manifest = run_extraction(
            input_path=input_dir,
            output_dir=output_dir,
            gene="TEST",
        )

        # Should return empty manifest without error
        assert len(manifest) == 0
        assert manifest.stage == Stage.EXTRACT


class TestManifestStageValue:
    """Test that extraction manifest uses correct stage."""

    @patch("cli.extract.ExpertExtractor")
    def test_manifest_has_extract_stage(self, mock_extractor_class, tmp_path):
        """Verify output manifest has EXTRACT stage."""
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        output_dir = tmp_path / "output"

        (input_dir / "12345_DATA_ZONES.md").write_text("test")

        mock_extractor = Mock()
        mock_result = Mock()
        mock_result.success = True
        mock_result.extracted_data = {"variants": []}
        mock_result.model_used = "test"
        mock_result.error = None
        mock_extractor.extract.return_value = mock_result
        mock_extractor_class.return_value = mock_extractor

        manifest = run_extraction(
            input_path=input_dir,
            output_dir=output_dir,
            gene="TEST",
        )

        assert manifest.stage == Stage.EXTRACT

        # Also check saved manifest
        loaded = Manifest.load(output_dir / "extraction_manifest.json")
        assert loaded.stage == Stage.EXTRACT
