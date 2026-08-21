"""Production evaluation source binding must reflect what extraction read."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from benchmarks.codex_paper_eval import rebind_production_sources as rebind


def _digest(path: Path) -> str:
    return rebind.digest(path)


def _write_extraction(record: Path, source: Path) -> None:
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(
        json.dumps(
            {
                "extraction_metadata": {
                    "source_file": str(source),
                    "source_sha256": _digest(source),
                    "source_size_bytes": source.stat().st_size,
                },
                "variants": [],
            }
        )
    )


def _write_completed_status(gene_run: Path, gene: str) -> None:
    (gene_run / f"{gene}.db").touch()
    (gene_run / "RUN_STATUS.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "exit_code": 0,
                "stage_failures": [],
                "active_db": f"{gene}.db",
            }
        )
    )


def test_extraction_source_uses_and_verifies_persisted_source_metadata(tmp_path: Path):
    gene_run = tmp_path / "production" / "SCN5A" / "run"
    pmc = gene_run / "pmc_fulltext"
    pmc.mkdir(parents=True)
    full = pmc / "123_FULL_CONTEXT.md"
    cleaned = pmc / "123_CLEANED.md"
    full.write_text("full")
    cleaned.write_text("cleaned source actually read")
    record = gene_run / "extractions" / "SCN5A_PMID_123.json"
    _write_extraction(record, cleaned)
    _write_completed_status(gene_run, "SCN5A")

    source, actual_record, status = rebind.extraction_source(
        gene_run, "SCN5A", "123", full
    )

    assert source == cleaned.resolve()
    assert actual_record == record.resolve()
    assert status == "extraction_metadata_verified"


def test_extraction_source_fails_closed_on_changed_input(tmp_path: Path):
    gene_run = tmp_path / "production" / "SCN5A" / "run"
    source = gene_run / "pmc_fulltext" / "123_CLEANED.md"
    source.parent.mkdir(parents=True)
    source.write_text("original")
    record = gene_run / "extractions" / "SCN5A_PMID_123.json"
    _write_extraction(record, source)
    source.write_text("changed after extraction")

    with pytest.raises(SystemExit, match="source SHA mismatch"):
        rebind.extraction_source(gene_run, "SCN5A", "123", source)


def test_rebind_main_locks_all_text_representations_and_exact_primary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run_dir = tmp_path / "eval"
    production_root = run_dir / "production_runs"
    gene_run = production_root / "SCN5A" / "20260821_000000"
    pmc = gene_run / "pmc_fulltext"
    pmc.mkdir(parents=True)
    full = pmc / "123_FULL_CONTEXT.md"
    cleaned = pmc / "123_CLEANED.md"
    zones = pmc / "123_DATA_ZONES.md"
    full.write_text("full representation")
    cleaned.write_text("cleaned representation actually read")
    zones.write_text("condensed representation")
    record = gene_run / "extractions" / "SCN5A_PMID_123.json"
    _write_extraction(record, cleaned)
    _write_completed_status(gene_run, "SCN5A")
    original = tmp_path / "eligibility.md"
    original.write_text("eligibility representation")
    run_dir.mkdir(exist_ok=True)
    (run_dir / "selection.json").write_text(
        json.dumps(
            {
                "papers": [
                    {
                        "gene": "SCN5A",
                        "pmid": "123",
                        "source": str(original),
                        "source_sha256": _digest(original),
                    }
                ]
            }
        )
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rebind_production_sources.py",
            "--run-dir",
            str(run_dir),
            "--production-root",
            str(production_root),
        ],
    )

    assert rebind.main() == 0

    selection = json.loads((run_dir / "selection.json").read_text())
    paper = selection["papers"][0]
    audit = selection["production_source_rebinding"]
    assert paper["source"] == str(cleaned.resolve())
    assert paper["source_selection"]["policy"] == "extraction_metadata_verified"
    assert paper["production_extraction_record"] == str(record.resolve())
    assert paper["representations"] == [
        str(cleaned.resolve()),
        str(zones.resolve()),
        str(full.resolve()),
    ]
    assert audit["binding_status_counts"] == {"extraction_metadata_verified": 1}
    assert audit["selected_representation_counts"] == {"CLEANED.md": 1}
