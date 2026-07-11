"""RUN_REPORT.md must surface non-fatal recovery-stage failures.

Regression guard: recovery subprocesses (run_all_layers.py, fetch_paywalled.py,
etc.) log-and-continue on a nonzero exit. Those failures must reach the durable
run status instead of a bare "✅ Done" hiding a swallowed error.
"""

from cli.gvf_run import step_report


def _report(
    tmp_path,
    stage_failures,
    *,
    stage_warnings=None,
    paper_final_check=None,
):
    out = tmp_path / "RUN_REPORT.md"
    step_report(
        gene="LDLR",
        db=tmp_path / "LDLR.db",
        run_dir=tmp_path,
        layer_outdir=None,
        source_qc_summary=None,
        source_recovery=None,
        doctor_status={},
        started_at=0.0,
        duration_s=60.0,
        out_path=out,
        stage_failures=stage_failures,
        stage_warnings=stage_warnings,
        paper_final_check=paper_final_check,
    )
    return out.read_text()


def test_report_surfaces_stage_failures(tmp_path):
    text = _report(tmp_path, ["recovery layers (run_all_layers.py) exited 1"])
    assert "## Stage Warnings" in text
    assert "run_all_layers.py) exited 1" in text
    assert "completed with warnings" in text


def test_report_clean_when_no_stage_failures(tmp_path):
    text = _report(tmp_path, [])
    assert "## Stage Warnings" not in text
    assert "all stages ok" in text


def test_report_surfaces_paper_final_check_errors_as_best_effort_warning(tmp_path):
    warning = "paper final check failed for all 2 checked paper(s)"
    text = _report(
        tmp_path,
        [],
        stage_warnings=[warning],
        paper_final_check={
            "papers": 3,
            "checked": 2,
            "skipped": 1,
            "skipped_empty_no_source": 1,
            "source_grounded": 1,
            "flagged_facts": 0,
            "error": 2,
            "missing_carriers": 0,
        },
    )
    assert "## Best-effort Warnings" in text
    assert warning in text
    assert "## Paper Final Check" in text
    assert "| 3 | 2 | 1 | 1 | 0 | 2 | 0 |" in text
    assert "explicitly skipped" in text
    assert "best-effort warnings recorded" in text


def test_step_layers_returns_none_on_failure(tmp_path, monkeypatch):
    """A failed recovery-layers subprocess must not hand back a bogus
    progression.csv path; it returns None and records the stage failure."""
    import cli.gvf_run as gvf_run

    class _Result:
        returncode = 1
        stderr = "boom"

    monkeypatch.setattr(gvf_run.subprocess, "run", lambda *a, **k: _Result())
    failures: list[str] = []
    out = gvf_run.step_layers(
        gene="LDLR",
        db=tmp_path / "LDLR.db",
        run_dir=tmp_path,
        gold=None,
        outdir=tmp_path / "layers",
        with_v12=None,
        stage_failures=failures,
    )
    assert out is None
    assert any("run_all_layers.py" in f for f in failures)
