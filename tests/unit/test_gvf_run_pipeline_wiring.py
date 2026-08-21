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
import os
from pathlib import Path

import pytest

import cli.gvf_run as gvf_run


@pytest.fixture(autouse=True)
def _skip_institutional_preflight(monkeypatch):
    """These tests exercise pipeline WIRING in a hermetic, mocked environment with
    no live EZproxy. Opt out of the institutional-access preflight (covered
    directly in test_institutional_preflight.py, and by
    ``test_full_run_blocks_when_institutional_access_degraded`` below) the same way
    CI/automation does, so a full-recovery run does not halt at the live probe.
    """
    monkeypatch.setenv("GVF_PREFLIGHT_SKIP", "1")
    # The wiring fixture writes a placeholder, not a real SQLite database. Keep
    # the deterministic composer inert unless a test explicitly replaces it.
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check_gate",
        lambda db: {"skipped": "wiring fixture"},
    )
    # Keep the wiring suite hermetic. The dedicated corpus-sync test below
    # replaces this stub with a recorder and verifies default/disabled behavior.
    # Calling the real step would scan and rewrite the user's local multi-GB
    # corpus once per test.
    monkeypatch.setattr(
        gvf_run,
        "step_corpus_sync",
        lambda run_dir, stage_warnings=None, gene=None: None,
    )


def _ok_doctor() -> dict:
    return {
        "ok": True,
        "required": {"NCBI_EMAIL": True},
        "recommended": {},
        "llm_providers": {"ANTHROPIC_API_KEY": True},
        "unlocks": {"ELSEVIER_INSTTOKEN": True},
        "ncbi_reachable": True,
    }


def _fake_extract_factory(captured: dict, name: str = "run1"):
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
        run_dir = Path(output_dir) / gene / name
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / f"{gene}.db").write_bytes(b"sqlite")
        return run_dir

    return fake_extract


def test_full_run_blocks_when_institutional_access_degraded(
    tmp_path: Path, monkeypatch
):
    """A full-dataset run (source_recovery on, no --pmid-file) must HALT at the
    institutional preflight when EZproxy is unconfigured — before any extraction —
    rather than silently harvesting abstract-only. This is the guard's whole point.
    """
    # Run the real guard, but mock the probe to report degraded access. The live
    # probe is non-deterministic (network) and run_gvf_pipeline reloads .env at
    # startup (initialize_runtime -> load_dotenv), so env-manipulation cannot
    # reliably simulate "no access". The probe's own logic is covered in
    # test_institutional_preflight.py; here we assert a should_block report HALTS
    # the run before extraction.
    monkeypatch.delenv("GVF_PREFLIGHT_SKIP", raising=False)
    import cli.institutional_preflight as ip

    monkeypatch.setattr(
        ip,
        "probe_institutional_access",
        lambda **kw: ip.AccessReport(
            viable=False,
            should_block=True,
            reason="degraded (test)",
            lines=["(test) EZproxy configured: NO"],
        ),
    )
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)

    def _extract_must_not_run(*args, **kwargs):
        raise AssertionError("extraction ran despite a blocked preflight")

    monkeypatch.setattr(gvf_run, "step_extract", _extract_must_not_run)

    rc = gvf_run.run_gvf_pipeline(
        gene="SCN5A",
        email="t@example.org",
        output=tmp_path,
        source_recovery=True,
        corpus_sync=False,
    )
    assert rc == gvf_run.EXIT_INSTITUTIONAL_BLOCK


def test_allow_degraded_institutional_overrides_block(tmp_path: Path, monkeypatch):
    """--allow-degraded-institutional lets a degraded full run PROCEED past the
    preflight, but the run is flagged degraded (recorded as a stage failure in
    RUN_STATUS.json) so it is never a silent success.
    """
    monkeypatch.delenv("GVF_PREFLIGHT_SKIP", raising=False)
    import cli.institutional_preflight as ip

    monkeypatch.setattr(
        ip,
        "probe_institutional_access",
        lambda **kw: ip.AccessReport(
            viable=False,
            should_block=True,
            reason="degraded (test)",
            lines=["(test) EZproxy configured: NO"],
        ),
    )
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    def fake_qc(gene, run_dir, outdir, stage_failures=None):
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
        return gvf_run.SourceRecoveryResult(
            fetch_dir=run_dir / "source_qc" / "fetch",
            outcome_summary=run_dir / "source_qc" / "outcome.json",
            fetched_source_override=run_dir / "source_qc" / "override.csv",
            active_db=None,
        )

    monkeypatch.setattr(gvf_run, "step_source_qc", fake_qc)
    monkeypatch.setattr(gvf_run, "step_source_recovery", fake_recovery)

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=True,  # full run -> the guard runs
        allow_degraded_institutional=True,  # ...but the operator overrides the block
        corpus_sync=False,
        skip=["layers"],
    )

    # The override let extraction + recovery proceed past the dead-access guard...
    assert captured.get("gene") == "TESTGENE"
    # ...and the run is flagged degraded (stage failure -> EXIT_STAGE_WARNINGS),
    # so an overridden run leaves a machine-readable trail, not a clean success.
    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    status = json.loads((output / "TESTGENE" / "run1" / "RUN_STATUS.json").read_text())
    assert any("institutional access degraded" in f for f in status["stage_failures"])


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
    assert status["gold_access"] == {
        "disabled": False,
        "mode": "auto_discovery_allowed",
    }
    assert (status_path.parent / status["active_db"]).resolve() == (
        report.parent / "TESTGENE.db"
    ).resolve()


def test_gold_free_run_never_discovers_or_passes_gold_to_layers(
    tmp_path: Path, monkeypatch
):
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "_find_gold",
        lambda gene: pytest.fail("gold discovery ran during --gold-free-run"),
    )

    def fake_layers(**kwargs):
        captured["layer_gold"] = kwargs["gold"]
        return None

    monkeypatch.setattr(gvf_run, "step_layers", fake_layers)

    gvf_run.run_gvf_pipeline(
        gene="KCNH2",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        corpus_sync=False,
        gold_free_run=True,
        skip=[
            "source-qc",
            "count-repair",
            "phenotype-count-guard",
            "claim-verification",
            "count-recovery",
            "vf-enrich",
            "trust-gate",
            "paper-final-check",
            "paper-final-check-gate",
            "metadata-backfill",
        ],
    )

    assert captured["layer_gold"] is None
    status = json.loads(next((tmp_path / "out").rglob("RUN_STATUS.json")).read_text())
    assert status["gold_access"] == {"disabled": True, "mode": "disabled"}


def test_count_recovery_runs_only_when_enabled_and_records_partial_failures(
    tmp_path: Path, monkeypatch
):
    captured: dict = {}
    calls: list[tuple[str, Path, Path]] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {"trusted": 0})
    monkeypatch.setenv("COUNT_RECOVERY_ENABLED", "true")
    monkeypatch.setattr(
        gvf_run,
        "step_count_recovery",
        lambda gene, db, run_dir: (
            calls.append((gene, db, run_dir))
            or {"papers_failed": 2, "counts_accepted": 3, "counts_written": 3}
        ),
    )

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        corpus_sync=False,
        skip=[
            "layers",
            "source-qc",
            "paper-final-check",
            "paper-final-check-gate",
            "metadata-backfill",
        ],
    )

    assert rc == 0
    assert len(calls) == 1
    gene, db, run_dir = calls[0]
    assert gene == "TESTGENE"
    assert db == run_dir / "TESTGENE.db"
    status = json.loads((run_dir / "RUN_STATUS.json").read_text())
    assert (
        "count recovery had model/parse failures for 2 paper(s)"
        in status["stage_warnings"]
    )
    # Stats reach RUN_STATUS.json instead of being logged and dropped.
    assert status["count_recovery"]["status"] == "ran"
    assert status["count_recovery"]["counts_written"] == 3


def test_count_recovery_refuses_when_trust_gate_is_skipped(tmp_path: Path, monkeypatch):
    """Recovered counts land as quarantine, so 3.7 must run to promote them."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setenv("COUNT_RECOVERY_ENABLED", "true")
    monkeypatch.setattr(
        gvf_run,
        "step_count_recovery",
        lambda **kwargs: pytest.fail("recovery ran with the trust gate skipped"),
    )

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        corpus_sync=False,
        skip=[
            "layers",
            "source-qc",
            "trust-gate",
            "paper-final-check",
            "paper-final-check-gate",
            "metadata-backfill",
        ],
    )

    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    run_dir = output / "TESTGENE" / "run1"
    status = json.loads((run_dir / "RUN_STATUS.json").read_text())
    assert status["count_recovery"]["status"] == "refused"
    assert any("count recovery refused" in f for f in status["stage_failures"])


def test_total_count_recovery_failure_is_a_stage_failure(tmp_path: Path, monkeypatch):
    """Every paper failing is a failure, not a silently-successful no-op."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {"trusted": 0})
    monkeypatch.setenv("COUNT_RECOVERY_ENABLED", "true")
    monkeypatch.setattr(
        gvf_run,
        "step_count_recovery",
        lambda gene, db, run_dir: {
            "papers_attempted": 4,
            "papers_failed": 4,
            "batch_failures": 4,
            "counts_accepted": 0,
            "counts_written": 0,
            "all_batches_failed": True,
        },
    )

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        corpus_sync=False,
        skip=[
            "layers",
            "source-qc",
            "paper-final-check",
            "paper-final-check-gate",
            "metadata-backfill",
        ],
    )

    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    status = json.loads((output / "TESTGENE" / "run1" / "RUN_STATUS.json").read_text())
    assert status["severity"] == "warning"
    assert status["count_recovery"]["status"] == "failed"
    assert any("failed on every attempted paper" in f for f in status["stage_failures"])


def test_count_recovery_explicit_skip_wins_over_enabled_flag(
    tmp_path: Path, monkeypatch
):
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setenv("COUNT_RECOVERY_ENABLED", "true")
    monkeypatch.setattr(
        gvf_run,
        "step_count_recovery",
        lambda **kwargs: pytest.fail("explicitly skipped count recovery ran"),
    )

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        corpus_sync=False,
        skip=[
            "layers",
            "source-qc",
            "count-recovery",
            "trust-gate",
            "paper-final-check",
            "paper-final-check-gate",
            "metadata-backfill",
        ],
    )
    assert rc == 0


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
        lambda run_dir, stage_warnings=None, gene=None: calls.append((run_dir, gene)),
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
    assert calls[0][1] == "TESTGENE", (
        "the run's gene must reach corpus sync as the --assume-gene hint"
    )

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
        gvf_run,
        "step_corpus_sync",
        lambda run_dir, stage_warnings=None, gene=None: None,
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
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check_gate",
        lambda db: {"applied_facts": 0, "unresolved_actionable": 0},
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
    gate_calls: list = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check",
        lambda db, run_dir, gene: calls.append(db) or {},
    )
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check_gate",
        lambda db: (
            gate_calls.append(db) or {"applied_facts": 0, "unresolved_actionable": 0}
        ),
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
    assert len(gate_calls) == 1, "stored grounded findings must still be recomposed"
    report = (tmp_path / "out" / "TESTGENE" / "run1" / "RUN_REPORT.md").read_text(
        encoding="utf-8"
    )
    assert "live reviewer was skipped" in report
    assert "Gate-applied facts" in report


def test_disabled_live_reviewer_still_reports_durable_gate(tmp_path: Path, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {})
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check",
        lambda db, run_dir, gene: {"skipped": "disabled by configuration"},
    )
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check_gate",
        lambda db: {
            "papers": 1,
            "applied_facts": 2,
            "applied_bindings": 2,
            "unresolved_actionable": 0,
            "unverified_actionable": 0,
            "missing_carriers": 0,
            "ungrounded_actionable": 0,
            "stale_papers": 0,
            "advisory_reason_flags": 0,
            "trusted": 8,
            "quarantine": 2,
        },
    )

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        skip=["layers", "source-qc"],
    )

    assert rc == 0
    report = (output / "TESTGENE" / "run1" / "RUN_REPORT.md").read_text()
    assert "_Skipped: disabled by configuration_" in report
    assert "Gate-applied facts" in report
    assert "| 2 | 2 | 0 | 0 | 0 | 0 | 8 | 2 |" in report


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
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check_gate",
        lambda db: {
            "applied_facts": 0,
            "unresolved_actionable": 0,
            "missing_carriers": 0,
            "ungrounded_actionable": 0,
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


def test_paper_final_check_gate_makes_missing_or_unbound_evidence_nonzero(
    tmp_path: Path, monkeypatch
):
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {})
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check",
        lambda db, run_dir, gene: {
            "papers": 1,
            "checked": 1,
            "skipped": 0,
            "source_grounded": 1,
            "flagged_facts": 1,
            "missing_carriers": 1,
            "error": 0,
        },
    )
    monkeypatch.setattr(
        gvf_run,
        "step_paper_final_check_gate",
        lambda db: {
            "applied_facts": 0,
            "unresolved_actionable": 1,
            "missing_carriers": 1,
            "ungrounded_actionable": 0,
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

    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    status = json.loads((output / "TESTGENE" / "run1" / "RUN_STATUS.json").read_text())
    assert any("binding" in value for value in status["stage_failures"])
    assert any("re-extraction required" in value for value in status["stage_failures"])


# ---------------------------------------------------------------------------
# Per-run trace isolation under GVF_LLM_TRACE_DIR. Each of these fails when the
# override is treated as one exact directory rather than a storage base.
# ---------------------------------------------------------------------------


def _trace_children(base: Path) -> list[str]:
    return sorted(p.name for p in base.iterdir() if p.is_dir())


def _run_twice(tmp_path: Path, monkeypatch, *, skip_extract_second: bool) -> Path:
    """Run the pipeline twice in ONE process and return the trace base."""
    from utils.llm_trace import reset_llm_tracing

    base = tmp_path / "vol"
    base.mkdir()
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", str(base))
    monkeypatch.delenv("GVF_LLM_TRACE_RUN_ID", raising=False)
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {"trusted": 0})

    output = tmp_path / "out"
    common = dict(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        corpus_sync=False,
        skip=[
            "layers",
            "source-qc",
            "paper-final-check",
            "paper-final-check-gate",
            "metadata-backfill",
        ],
    )

    captured: dict = {}
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    assert gvf_run.run_gvf_pipeline(**common) == 0
    reset_llm_tracing()

    if skip_extract_second:
        monkeypatch.setattr(
            gvf_run,
            "step_extract",
            lambda **kwargs: pytest.fail("extract ran on the skip-extract pass"),
        )
        second = dict(common, skip=[*common["skip"], "extract"])
    else:
        captured2: dict = {}
        monkeypatch.setattr(
            gvf_run, "step_extract", _fake_extract_factory(captured2, name="run2")
        )
        second = common
    assert gvf_run.run_gvf_pipeline(**second) == 0
    reset_llm_tracing()
    return base


def test_trace_env_identity_is_restored_after_a_run(tmp_path: Path, monkeypatch):
    """A leaked run id makes the NEXT run adopt this run's identity."""
    from utils.llm_trace import reset_llm_tracing

    base = tmp_path / "vol"
    base.mkdir()
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", str(base))
    monkeypatch.delenv("GVF_LLM_TRACE_RUN_ID", raising=False)
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_trust_gate", lambda db: {"trusted": 0})
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))

    assert (
        gvf_run.run_gvf_pipeline(
            gene="TESTGENE",
            email="x@example.com",
            output=tmp_path / "out",
            source_recovery=False,
            corpus_sync=False,
            skip=[
                "layers",
                "source-qc",
                "paper-final-check",
                "paper-final-check-gate",
                "metadata-backfill",
            ],
        )
        == 0
    )
    reset_llm_tracing()

    assert os.environ["GVF_LLM_TRACE_DIR"] == str(base)
    assert "GVF_LLM_TRACE_RUN_ID" not in os.environ


def test_two_normal_runs_under_one_base_get_separate_trace_children(
    tmp_path: Path, monkeypatch
):
    base = _run_twice(tmp_path, monkeypatch, skip_extract_second=False)

    children = _trace_children(base)
    assert len(children) == 2, f"runs shared a trace directory: {children}"
    assert all(name.startswith("gvfrun-") for name in children), children
    # Each child's manifest is attributed to exactly one run, so a later rebuild
    # never raises TraceRunMismatchError. (Extraction is mocked here, so the
    # record sets are empty; the point is that the trees are separate and each
    # manifest names one run.)
    from utils.llm_trace import build_trace_manifest

    manifest_run_ids = set()
    for name in children:
        manifest = build_trace_manifest(base / name, run_id=name)
        assert set(manifest["record_run_ids"]) <= {name}
        manifest_run_ids.add(manifest["run_id"])
    assert manifest_run_ids == set(children)


def test_skip_extract_second_run_does_not_join_the_first_runs_tree(
    tmp_path: Path, monkeypatch
):
    """Skip-extract reuses the run DIR; it must not reuse the run's trace tree."""
    base = _run_twice(tmp_path, monkeypatch, skip_extract_second=True)

    children = _trace_children(base)
    assert len(children) == 2, f"skip-extract joined the earlier tree: {children}"


def test_standalone_recovery_gets_a_per_run_child_of_the_base(
    tmp_path: Path, monkeypatch
):
    """Two sequential standalone recoveries must not share a trace directory."""
    import sqlite3

    from utils.llm_trace import reset_llm_tracing
    from scripts.recover_counts import main as recover_main

    db = tmp_path / "T.db"
    con = sqlite3.connect(db)
    con.executescript(
        """
        CREATE TABLE variants(variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT);
        CREATE TABLE variant_papers(variant_id INTEGER, pmid TEXT,
            source_location TEXT, key_quotes TEXT, source_layer TEXT);
        CREATE TABLE penetrance_data(penetrance_id INTEGER PRIMARY KEY,
            variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
            affected_count INTEGER, unaffected_count INTEGER);
        """
    )
    con.commit()
    con.close()

    base = tmp_path / "vol"
    base.mkdir()
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", str(base))
    monkeypatch.delenv("GVF_LLM_TRACE_RUN_ID", raising=False)

    for _ in range(2):
        assert (
            recover_main(
                ["--db", str(db), "--gene", "TESTGENE", "--dry-run", "--no-backup"]
            )
            == 0
        )
        reset_llm_tracing()

    children = _trace_children(base)
    assert len(children) == 2, f"sequential recoveries shared a directory: {children}"
    assert all(name.startswith("recover-counts-") for name in children), children


def test_vf_enrich_runs_outside_full_coverage_when_warehouse_available(
    tmp_path: Path, monkeypatch
):
    """The wrong-gene FP quarantine must not be full-coverage-only.

    It is the only check that validates a variant's gene against an independent
    reference, and it was reachable only in full-coverage discovery mode — so
    every calibrated `--pmid-file` run, which is how the collaborator-facing
    review sets are built, silently skipped it. Even the "SKIPPED" log line sat
    inside the same branch, so nothing reported that it had not run.
    """
    captured: dict = {}
    calls: list[str] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        gvf_run, "step_vf_enrich", lambda **kwargs: calls.append("vf") or {}
    )
    warehouse = tmp_path / "variants.db"
    warehouse.write_text("")
    monkeypatch.setattr(
        "utils.gene_metadata.default_variantfeatures_db_path", lambda: warehouse
    )
    monkeypatch.setattr("utils.env_utils.local_data_discovery_disabled", lambda: False)

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
    )
    assert rc == 0
    assert calls == ["vf"]


def test_vf_enrich_unavailable_is_a_warning_not_silence(tmp_path: Path, monkeypatch):
    """An absent warehouse must be loud. Silently skipping the only independent
    gene check is how wrong-gene rows reached collaborator staging."""
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        "utils.gene_metadata.default_variantfeatures_db_path", lambda: None
    )
    monkeypatch.setattr("utils.env_utils.local_data_discovery_disabled", lambda: False)

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        source_recovery=False,
        skip=["layers", "source-qc"],
    )
    assert rc == 0
    statuses = list((tmp_path / "out").rglob("RUN_STATUS.json"))
    assert statuses, "no RUN_STATUS.json written"
    status = json.loads(statuses[0].read_text())
    warnings = " ".join(status.get("stage_warnings") or [])
    assert "vf-enrich" in warnings and "VARIANTFEATURES_DB" in warnings


def test_required_vf_enrich_unavailable_fails_acceptance(tmp_path: Path, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        "utils.gene_metadata.default_variantfeatures_db_path", lambda: None
    )
    monkeypatch.setattr("utils.env_utils.local_data_discovery_disabled", lambda: False)

    output = tmp_path / "out"
    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=output,
        source_recovery=False,
        skip=["layers", "source-qc"],
        require_vf_enrich=True,
    )
    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    status = json.loads((output / "TESTGENE" / "run1" / "RUN_STATUS.json").read_text())
    failures = " ".join(status["stage_failures"])
    assert "vf-enrich" in failures and "VARIANTFEATURES_DB" in failures


def test_publish_refused_when_any_required_stage_failed(tmp_path: Path, monkeypatch):
    captured: dict = {}
    publish_calls: list[dict] = []
    monkeypatch.setattr(gvf_run, "doctor", _ok_doctor)
    monkeypatch.setattr(gvf_run, "step_extract", _fake_extract_factory(captured))
    monkeypatch.setattr(
        "utils.gene_metadata.default_variantfeatures_db_path", lambda: None
    )
    monkeypatch.setattr("utils.env_utils.local_data_discovery_disabled", lambda: False)
    monkeypatch.setattr(
        gvf_run,
        "step_publish_review",
        lambda **kwargs: publish_calls.append(kwargs) or True,
    )
    manifest = tmp_path / "TESTGENE.txt"
    manifest.write_text("12345\n")

    rc = gvf_run.run_gvf_pipeline(
        gene="TESTGENE",
        email="x@example.com",
        output=tmp_path / "out",
        pmid_file=manifest,
        source_recovery=False,
        publish_review=True,
        skip=["layers", "source-qc"],
    )
    assert rc == gvf_run.EXIT_STAGE_WARNINGS
    assert publish_calls == []
