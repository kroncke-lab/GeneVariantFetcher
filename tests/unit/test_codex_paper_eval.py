"""Regression coverage for the extraction-blinded Codex paper evaluation."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from benchmarks.codex_paper_eval.build_report_artifact import build_payload
from benchmarks.codex_paper_eval.run_eval import (
    command_extract,
    digest,
    matches,
    material_digest_errors,
    selection_metadata,
    write_json,
    write_markdown_report,
)


COUNT_FIELDS = ("carriers", "affected", "unaffected")
GENES = ("SCN5A", "KCNH2", "KCNQ1", "RYR2")


def _count_metric(asserted: int = 1, predicted: int = 1) -> dict:
    return {
        "gold_asserted": asserted,
        "predicted": predicted,
        "recall": predicted / asserted if asserted else None,
        "mae": 0.0 if predicted else None,
        "rmse": 0.0 if predicted else None,
    }


def _report_fixture() -> dict:
    papers = []
    by_gene = {}
    for index, gene in enumerate(GENES, 1):
        count = {field: _count_metric() for field in COUNT_FIELDS}
        paper = {
            "gene": gene,
            "pmid": str(index),
            "tool": "text",
            "tool_rationale": "Running text contained the evidence.",
            "source_completeness": "full_text",
            "elapsed_seconds": 1.0,
            "token_usage": {
                "input_tokens": 20,
                "output_tokens": 5,
                "total_tokens": 25,
            },
            "tp": 1,
            "fp": 0,
            "fn": 0,
            "precision": 1.0,
            "recall": 1.0,
            "f1": 1.0,
            "count": count,
            "matched_variants": [{"predicted": "A1V", "gold": "A1V"}],
            "missed_gold": [],
            "extra_predictions": [],
            "count_errors": [],
        }
        papers.append(paper)
        by_gene[gene] = {
            "papers": 1,
            "tp": 1,
            "fp": 0,
            "fn": 0,
            "precision": 1.0,
            "recall": 1.0,
            "f1": 1.0,
            "elapsed_seconds": 1.0,
            "token_usage": {
                "input_tokens": 20,
                "output_tokens": 5,
                "total_tokens": 25,
            },
            "count": count,
        }
    overall = {
        "papers": 4,
        "tp": 4,
        "fp": 0,
        "fn": 0,
        "precision": 1.0,
        "recall": 1.0,
        "f1": 1.0,
        "elapsed_seconds": 4.0,
        "token_usage": {
            "input_tokens": 80,
            "output_tokens": 20,
            "total_tokens": 100,
        },
        "count": {field: _count_metric(4, 4) for field in COUNT_FIELDS},
    }
    return {
        "run_id": "fixture",
        "seed": 7,
        "locked_at": "2026-07-24T00:00:00+00:00",
        "scored_at": "2026-07-24T00:01:00+00:00",
        "overall": overall,
        "by_gene": by_gene,
        "papers": papers,
        "selection": {
            "mode": "manifest",
            "population": "fixed manifest `fixture.tsv` (4 papers)",
            "description": (
                "Paper selection used the fixed manifest `fixture.tsv` (4 papers). "
                "Routing and extraction were gold-value-blind."
            ),
        },
        "tools_used": {"text": 4},
        "token_usage": {
            "input_tokens": 100,
            "output_tokens": 23,
            "total_tokens": 123,
        },
        "timing": {
            "wall_seconds": 60.0,
            "summed_paper_seconds": 4.0,
            "started_at": "2026-07-24T00:00:00+00:00",
            "completed_at": "2026-07-24T00:01:00+00:00",
        },
    }


def test_deletion_range_does_not_match_single_residue_deletion():
    assert not matches("K1505del", "K1505_Q1507del", "SCN5A")
    assert not matches("K1505_Q1507del", "K1505del", "SCN5A")
    assert matches("K1505_Q1507del", "K1505_Q1507del", "SCN5A")


def test_material_digests_cover_source_artifact_pdf_and_figure(tmp_path: Path):
    source = tmp_path / "source.md"
    artifact = tmp_path / "artifact.json"
    pdf = tmp_path / "paper.pdf"
    figure = tmp_path / "figure.png"
    source.write_text("source")
    artifact.write_text("{}")
    pdf.write_bytes(b"%PDF fixture")
    figure.write_bytes(b"PNG fixture")
    paper = {
        "gene": "SCN5A",
        "pmid": "1",
        "source": str(source),
        "source_sha256": digest(source),
        "artifacts": str(artifact),
        "artifacts_sha256": digest(artifact),
        "pdfs": [str(pdf)],
        "pdf_sha256": {str(pdf): digest(pdf)},
        "figures": [str(figure)],
        "figure_sha256": {str(figure): digest(figure)},
    }

    assert material_digest_errors(paper) == []

    figure.write_bytes(b"changed")
    assert any(
        "figure changed after selection" in error
        for error in material_digest_errors(paper)
    )


def test_api_usage_is_checkpointed_before_response_parsing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source = tmp_path / "source.md"
    source.write_text("SCN5A A1V carrier evidence")
    selection = {
        "papers": [
            {
                "gene": "SCN5A",
                "pmid": "1",
                "source": str(source),
                "source_sha256": digest(source),
                "artifacts": None,
                "artifacts_sha256": None,
                "pdfs": [],
                "pdf_sha256": {},
                "figures": [],
                "figure_sha256": {},
            }
        ]
    }
    predictions = {
        "token_usage": None,
        "papers": [
            {
                "gene": "SCN5A",
                "pmid": "1",
                "tool": None,
                "token_usage": None,
                "variants": [],
            }
        ],
    }
    write_json(tmp_path / "selection.json", selection)
    write_json(tmp_path / "predictions.json", predictions)

    responses = iter(
        [
            SimpleNamespace(
                output_text=json.dumps(
                    {
                        "tool": "text",
                        "tool_rationale": "Only text was available.",
                        "source_completeness": "full_text",
                    }
                ),
                usage=SimpleNamespace(input_tokens=5, output_tokens=1),
            ),
            SimpleNamespace(
                output_text="not valid JSON",
                usage=SimpleNamespace(input_tokens=7, output_tokens=2),
            ),
        ]
    )
    fake_openai = ModuleType("openai")
    fake_openai.OpenAI = lambda **_kwargs: SimpleNamespace(
        responses=SimpleNamespace(create=lambda **_kwargs: next(responses))
    )
    monkeypatch.setitem(sys.modules, "openai", fake_openai)
    monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.invalid")
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-key")
    args = SimpleNamespace(
        run_dir=tmp_path,
        timeout=1,
        model="test-model",
        force=False,
        max_artifact_chars=100,
        max_source_chars=1000,
        route_preview_chars=500,
        max_ocr_images=1,
        route_reasoning_effort="medium",
        reasoning_effort="high",
        max_output_tokens=100,
    )

    with pytest.raises(json.JSONDecodeError):
        command_extract(args)

    checkpoint = json.loads((tmp_path / "predictions.json").read_text())
    assert checkpoint["token_usage"]["total_tokens"] == 15
    assert checkpoint["papers"][0]["token_usage"]["total_tokens"] == 15


def test_selection_metadata_describes_manifest_and_random_modes(tmp_path: Path):
    manifest = tmp_path / "comparison.tsv"
    manifest.write_text("SCN5A\t1\n")
    fixed = selection_metadata(
        {
            "paper_manifest": str(manifest),
            "papers": [{"gene": "SCN5A", "pmid": "1"}],
            "seed": 7,
            "per_gene": 5,
        }
    )
    random = selection_metadata(
        {
            "paper_manifest": None,
            "papers": [{"gene": "SCN5A", "pmid": "1"}] * 4,
            "seed": 11,
            "per_gene": 1,
        }
    )

    assert fixed["mode"] == "manifest"
    assert "comparison.tsv" in fixed["description"]
    assert random["mode"] == "random"
    assert "seed 11" in random["description"]
    assert "high-carrier" not in fixed["description"] + random["description"]


def test_markdown_and_artifact_narratives_are_derived_from_report(tmp_path: Path):
    report = _report_fixture()
    write_json(tmp_path / "report.json", report)

    markdown_path = tmp_path / "report.md"
    write_markdown_report(report, markdown_path)
    markdown = markdown_path.read_text()
    assert "fixture.tsv" in markdown
    assert "gold-ranked" not in markdown
    assert "high-carrier" not in markdown

    payload = build_payload(tmp_path)
    manifest = payload["manifest"]
    executive = next(
        block["body"]
        for block in manifest["blocks"]
        if block["id"] == "executive_summary"
    )
    assert manifest["title"] == "Codex 4-Paper Extraction Evaluation"
    assert "100.0% precision" in executive
    assert "4 gold variant rows" in executive
    assert "123 exact API tokens" in executive
    assert "93.1%" not in json.dumps(payload)
