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
        gene,
        email,
        output_dir,
        pmid_file,
        max_pmids,
        resume_dir,
        disease=None,
        **kwargs,
    ):
        captured["disease"] = disease
        captured["gene"] = gene
        captured["extract_kwargs"] = kwargs
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

    status_path = report.parent / "RUN_STATUS.json"
    status = json.loads(status_path.read_text())
    assert status["status"] == "completed"
    assert status["severity"] == "ok"
    assert status["active_db"] == "TESTGENE.db"
    assert (status_path.parent / status["active_db"]).resolve() == (
        report.parent / "TESTGENE.db"
    ).resolve()


def test_source_recovery_on_by_default(tmp_path: Path, monkeypatch):
    """With no flag passed, the turnkey driver runs source-qc + source-recovery."""
    captured: dict = {}
    calls: list[str] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    def fake_qc(gene, run_dir, outdir, stage_failures=None):
        calls.append("source-qc")
        outdir.mkdir(parents=True, exist_ok=True)
        summary = outdir / "summary.json"
        summary.write_text(json.dumps({"pmid_coverage": {}}), encoding="utf-8")
        return summary

    def fake_recovery(
        gene,
        run_dir,
        source_qc_dir,
        gold,
        run_recovery_layers,
        timeout_s,
        stage_failures=None,
    ):
        calls.append("source-recovery")
        refreshed_db = run_dir / f"{gene}.refresh.db"
        refreshed_db.write_bytes(b"refreshed sqlite")
        return gvf_run.SourceRecoveryResult(
            fetch_dir=run_dir / "source_qc" / "fetch",
            outcome_summary=run_dir / "source_qc" / "outcome.json",
            fetched_source_override=run_dir / "source_qc" / "override.csv",
            active_db=refreshed_db,
        )

    monkeypatch.setattr(gvf_run, "step_source_qc", fake_qc)
    monkeypatch.setattr(gvf_run, "step_source_recovery", fake_recovery)

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        # NOTE: source_recovery intentionally omitted -> exercises the default
        skip=["layers"],
    )

    assert rc == 0
    assert "source-qc" in calls, "source-qc did not run under the default"
    assert "source-recovery" in calls, "source-recovery is not ON by default"
    status_path = output / "TESTGENE" / "run1" / "RUN_STATUS.json"
    status = json.loads(status_path.read_text())
    assert status["active_db"] == "TESTGENE.refresh.db"
    assert (
        status_path.parent / status["active_db"]
    ).read_bytes() == b"refreshed sqlite"


def test_corpus_sync_runs_by_default_and_is_skippable(tmp_path: Path, monkeypatch):
    """The write-back corpus-sync step runs by default; --no-corpus-sync skips it."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    calls: list = []
    monkeypatch.setattr(
        gvf_run,
        "step_corpus_sync",
        lambda run_dir, stage_warnings=None: calls.append(run_dir),
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


def test_full_coverage_is_opt_in(tmp_path: Path, monkeypatch):
    """Default gvf-run keeps the bounded/non-full-coverage path untouched."""
    captured: dict = {}
    calls: list[str] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "step_full_coverage_walk",
        lambda **kwargs: calls.append("walk") or {},
    )
    monkeypatch.setattr(
        gvf_run, "step_carrier_guard", lambda **kwargs: calls.append("guard") or {}
    )
    monkeypatch.setattr(
        gvf_run, "step_vf_enrich", lambda **kwargs: calls.append("vf") or {}
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
    )

    assert rc == 0
    assert calls == []
    assert captured["extract_kwargs"].get("extraction_top_n") is None


def test_full_coverage_seeds_and_walks_from_next_offset(tmp_path: Path, monkeypatch):
    """--full-coverage seeds a bounded priority batch, then continues walking."""
    captured: dict = {}
    calls: list[tuple[str, dict]] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    def fake_walk(**kwargs):
        calls.append(("walk", kwargs))
        return {"walked": True}

    monkeypatch.setattr(gvf_run, "step_full_coverage_walk", fake_walk)
    monkeypatch.setattr(
        gvf_run,
        "step_carrier_guard",
        lambda **kwargs: calls.append(("guard", kwargs)) or {},
    )
    monkeypatch.setattr(
        gvf_run, "step_vf_enrich", lambda **kwargs: calls.append(("vf", kwargs)) or {}
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
        full_coverage=True,
        extraction_model="azure_ai/gpt-5.4",
        extraction_workers=7,
        taper_min_variants=11,
    )

    assert rc == 0
    assert captured["extract_kwargs"]["extraction_top_n"] == 1000
    assert calls[0][0] == "walk"
    assert calls[0][1]["model"] == "azure_ai/gpt-5.4"
    assert calls[0][1]["max_workers"] == 7
    assert calls[0][1]["start_offset"] == 1000
    assert calls[0][1]["min_new_variants"] == 11
    assert [name for name, _ in calls] == ["walk", "guard", "vf"]


def test_full_coverage_quality_steps_are_toggleable(tmp_path: Path, monkeypatch):
    captured: dict = {}
    calls: list[str] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "step_full_coverage_walk",
        lambda **kwargs: calls.append("walk") or {},
    )
    monkeypatch.setattr(
        gvf_run, "step_carrier_guard", lambda **kwargs: calls.append("guard") or {}
    )
    monkeypatch.setattr(
        gvf_run, "step_vf_enrich", lambda **kwargs: calls.append("vf") or {}
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
        full_coverage=True,
        carrier_guard=False,
        vf_enrich=False,
    )

    assert rc == 0
    assert calls == ["walk"]


def test_stage_failure_sets_nonzero_exit_and_status(tmp_path: Path, monkeypatch):
    """A completeness-stage failure -> non-zero exit (EXIT_STAGE_WARNINGS) plus a
    machine-readable RUN_STATUS.json, so a fleet detects the incomplete run
    instead of reading a clean exit 0. Also guards the no-double-run fix."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run, "step_corpus_sync", lambda run_dir, stage_warnings=None: None
    )
    monkeypatch.setattr(gvf_run, "step_backfill_metadata", lambda **kwargs: None)

    qc_calls: list[str] = []

    def failing_qc(gene, run_dir, outdir, stage_failures=None):
        qc_calls.append("qc")
        if stage_failures is not None:
            stage_failures.append("source-qc (source_acquisition_audit.py) exited 1")
        return None  # QC failed -> no summary

    monkeypatch.setattr(gvf_run, "step_source_qc", failing_qc)

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=True,
        skip=["layers"],
    )

    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    # A failed Step-3 QC must not be re-run by the Step-4 fallback.
    assert qc_calls == ["qc"], "source-qc ran twice after a failure"

    status = json.loads((output / "TESTGENE" / "run1" / "RUN_STATUS.json").read_text())
    assert status["status"] == "completed_with_warnings"
    assert status["severity"] == "warning"
    assert status["exit_code"] == gvf_run.EXIT_STAGE_WARNINGS
    assert status["active_db"] == "TESTGENE.db"
    assert any("source-qc" in f for f in status["stage_failures"])


def test_fatal_extract_failure_does_not_write_completed_status(
    tmp_path: Path, monkeypatch
):
    """Fatal exit 3 must not leave a RUN_STATUS that claims completion."""
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)

    def failing_extract(*, gene, output_dir, **kwargs):
        # Model an extractor that creates its run directory before raising.
        (Path(output_dir) / gene / "run1").mkdir(parents=True)
        raise RuntimeError("fatal extraction failure")

    monkeypatch.setattr(gvf_run, "step_extract", failing_extract)

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
    )

    assert rc == 3
    assert not (output / "TESTGENE" / "run1" / "RUN_STATUS.json").exists()


def test_trust_gate_runs_by_default(tmp_path: Path, monkeypatch):
    """The per-fact trust gate runs by default — the primary quality control."""
    captured: dict = {}
    calls: list = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "step_trust_gate",
        lambda db: calls.append(db) or {"trusted": 0, "quarantine": 0},
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
    )
    assert rc == 0
    assert len(calls) == 1, "trust gate should run by default"


def test_trust_gate_is_skippable(tmp_path: Path, monkeypatch):
    captured: dict = {}
    calls: list = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: calls.append(db) or {})

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc", "trust-gate"],
    )
    assert rc == 0
    assert calls == [], "trust gate should be skippable via skip=['trust-gate']"


def test_paper_final_check_runs_by_default(tmp_path: Path, monkeypatch):
    """The per-paper final check (sniff test) runs by default — Step 3.8."""
    captured: dict = {}
    calls: list = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check",
        lambda db, run_dir, gene: calls.append(db) or {"papers": 0, "checked": 0},
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
    )
    assert rc == 0
    assert len(calls) == 1, "paper final check should run by default"


def test_paper_final_check_is_skippable(tmp_path: Path, monkeypatch):
    captured: dict = {}
    calls: list = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check",
        lambda db, run_dir, gene: calls.append(db) or {},
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc", "paper-final-check"],
    )
    assert rc == 0
    assert calls == [], "final check should be skippable via skip=['paper-final-check']"


def test_paper_final_check_errors_warn_and_reach_run_outputs(
    tmp_path: Path, monkeypatch, caplog
):
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {})
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check",
        lambda db, run_dir, gene: {
            "papers": 2,
            "checked": 2,
            "skipped": 0,
            "source_grounded": 0,
            "flagged_facts": 0,
            "missing_carriers": 0,
            "error": 2,
        },
    )

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        corpus_sync=False,
        skip=["layers", "source-qc", "metadata-backfill"],
    )

    warning = "paper final check failed for all 2 checked paper(s)"
    assert rc == 0  # soft QC warning does not discard an otherwise usable run
    assert warning in caplog.text
    run_dir = output / "TESTGENE" / "run1"
    report = (run_dir / "RUN_REPORT.md").read_text()
    status = json.loads((run_dir / "RUN_STATUS.json").read_text())
    assert warning in report
    assert "## Paper Final Check" in report
    assert warning in status["stage_warnings"]
    assert status["exit_code"] == 0
