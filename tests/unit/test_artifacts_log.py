"""Tests for harvesting.artifacts_log."""

from __future__ import annotations

import json
from pathlib import Path

from harvesting.artifacts_log import (
    ArtifactsLog,
    FigureArtifact,
    SupplementArtifact,
)


def test_artifacts_log_writes_summary(tmp_path: Path):
    log = ArtifactsLog(
        pmid="12345",
        output_dir=tmp_path,
        pmcid="PMC999",
        doi="10.1/test",
        gene_symbol="KCNH2",
    )
    log.record_main_text(
        source="pmc_xml",
        chars=1234,
        figure_captions=2,
        table_captions=1,
        supplement_descriptions=3,
    )
    fig = tmp_path / "fig1.png"
    fig.write_bytes(b"\x89PNGfake")
    log.record_figure(
        ArtifactsLog.figure_artifact_from_path(
            fig,
            source="pmc_html",
            caption_label="Figure 1",
            caption_text="A pedigree carrying KCNH2 G604S",
        )
    )
    log.record_supplement_dict(
        filename="mmc1.xlsx",
        path=str(tmp_path / "mmc1.xlsx"),
        url="https://example.com/mmc1.xlsx",
        size_bytes=2048,
        source="elsevier_api",
        description="Supplementary Table S1: 86 missense variants",
        converted=True,
        converted_chars=8000,
        figures_extracted=0,
    )
    log.add_note("first run")

    out = log.save()
    assert out.exists()
    data = json.loads(out.read_text())
    assert data["pmid"] == "12345"
    assert data["pmcid"] == "PMC999"
    assert data["main_text"]["chars"] == 1234
    assert data["main_text"]["figure_captions_count"] == 2
    assert len(data["figures"]) == 1
    assert data["figures"][0]["caption_label"] == "Figure 1"
    assert data["figures"][0]["size_bytes"] == len(b"\x89PNGfake")
    assert len(data["supplements"]) == 1
    assert data["supplements"][0]["source"] == "elsevier_api"
    assert data["supplements"][0]["converted_chars"] == 8000
    assert data["summary"]["figure_count"] == 1
    assert data["summary"]["figures_with_captions"] == 1
    assert data["summary"]["supplement_count"] == 1
    assert data["summary"]["supplements_converted"] == 1
    assert data["summary"]["supplements_total_chars"] == 8000
    assert data["notes"] == ["first run"]


def test_supplement_artifact_from_missing_path_records_zero_size(tmp_path: Path):
    missing = tmp_path / "does-not-exist.zip"
    artifact = ArtifactsLog.supplement_artifact_from_path(
        filepath=missing,
        url="https://example.com/x.zip",
        source="scraper",
    )
    assert isinstance(artifact, SupplementArtifact)
    assert artifact.size_bytes == 0
    assert not artifact.converted


def test_figure_artifact_records_size(tmp_path: Path):
    p = tmp_path / "fig.png"
    p.write_bytes(b"x" * 17)
    fig = ArtifactsLog.figure_artifact_from_path(p, source="pmc_html")
    assert isinstance(fig, FigureArtifact)
    assert fig.size_bytes == 17
    assert fig.filename == "fig.png"
