"""Offline regression net for the gvf-run cold-path orchestrator.

`run_gvf_pipeline` is the turnkey driver (doctor -> extract -> source-qc ->
source-recovery -> layers -> report). The expensive bits (NCBI/LLM extraction
and the recovery subprocesses) are network-bound, so this module monkeypatches
*only* those boundaries and asserts the chain itself wires together: that it
runs to completion, threads `--disease` into extraction, defaults
source-recovery ON, and emits a RUN_REPORT.md. This catches chain regressions
in seconds without a multi-hour live run (the previous gap: the cold path had
no non-mocked completion test at all).
"""

from __future__ import annotations

import json
from pathlib import Path

import cli.gvf_run as gvf_run


def _ok_doctor() -> dict:
    return {
        "ok": True,
        "required": {"NCBI_EMAIL": True},
        "recommended": {},
        "llm_providers": {"ANTHROPIC_API_KEY": True},
        "unlocks": {"ELSEVIER_INSTTOKEN": True},
        "ncbi_reachable": True,
    }


def _fake_extract_factory(captured: dict):
    """Return a step_extract stand-in that records its kwargs and seeds a DB."""

    def fake_extract(
        gene, email, output_dir, pmid_file, max_pmids, resume_dir, disease=None
    ):
        captured["disease"] = disease
        captured["gene"] = gene
        run_dir = Path(output_dir) / gene / "run1"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / f"{gene}.db").write_bytes(b"sqlite")
        return run_dir

    return fake_extract


def test_pipeline_completes_and_writes_report(tmp_path: Path, monkeypatch):
    """Full chain runs to a 0 exit code and writes RUN_REPORT.md (no gold)."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="testgene",
        email="x@example.com",
        output=output,
        source_recovery=False,
        skip=["layers", "source-qc"],
    )

    assert rc == 0
    # gene is upper-cased before the run dir is created
    report = output / "TESTGENE" / "run1" / "RUN_REPORT.md"
    assert report.exists(), "RUN_REPORT.md was not produced by the cold path"
    assert (output / "TESTGENE_RUN_REPORT.md").exists(), "report not copied to root"
    body = report.read_text()
    assert "# GVF Run Report — TESTGENE" in body


def test_source_recovery_on_by_default(tmp_path: Path, monkeypatch):
    """With no flag passed, the turnkey driver runs source-qc + source-recovery."""
    captured: dict = {}
    calls: list[str] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    def fake_qc(gene, run_dir, outdir):
        calls.append("source-qc")
        outdir.mkdir(parents=True, exist_ok=True)
        summary = outdir / "summary.json"
        summary.write_text(json.dumps({"pmid_coverage": {}}), encoding="utf-8")
        return summary

    def fake_recovery(
        gene, run_dir, source_qc_dir, gold, run_recovery_layers, timeout_s
    ):
        calls.append("source-recovery")
        return None

    monkeypatch.setattr(gvf_run, "step_source_qc", fake_qc)
    monkeypatch.setattr(gvf_run, "step_source_recovery", fake_recovery)

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        # NOTE: source_recovery intentionally omitted -> exercises the default
        skip=["layers"],
    )

    assert rc == 0
    assert "source-qc" in calls, "source-qc did not run under the default"
    assert "source-recovery" in calls, "source-recovery is not ON by default"


def test_corpus_sync_runs_by_default_and_is_skippable(tmp_path: Path, monkeypatch):
    """The write-back corpus-sync step runs by default; --no-corpus-sync skips it."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    calls: list = []
    monkeypatch.setattr(
        gvf_run, "step_corpus_sync", lambda run_dir: calls.append(run_dir)
    )

    # default: corpus_sync on
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "a",
        source_recovery=False,
        skip=["layers", "source-qc"],
    )
    assert rc == 0 and len(calls) == 1, "corpus sync should run by default"

    # explicit off
    calls.clear()
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "b",
        source_recovery=False,
        skip=["layers", "source-qc"],
        corpus_sync=False,
    )
    assert rc == 0 and calls == [], "corpus sync should be skipped when off"


def test_disease_threads_into_extraction(tmp_path: Path, monkeypatch):
    """--disease reaches step_extract (the gene-disease-pair path)."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
        disease="long QT syndrome",
    )

    assert rc == 0
    assert captured["disease"] == "long QT syndrome"
