"""RUN_REPORT.md must surface non-fatal recovery-stage failures.

Regression guard: recovery subprocesses (run_all_layers.py, fetch_paywalled.py,
etc.) log-and-continue on a nonzero exit. Those failures must reach the durable
run status instead of a bare "✅ Done" hiding a swallowed error.
"""

from cli.gvf_run import step_report


def _report(tmp_path, stage_failures):
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
