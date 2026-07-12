"""Focused tests for curated benchmark extraction-run status handling."""

import json
from pathlib import Path
from types import SimpleNamespace

from benchmarks.curated_extraction_eval import run_benchmark as benchmark


def _run_dir_from_cmd(cmd: list[str]) -> Path:
    output = Path(cmd[cmd.index("--output") + 1])
    gene = cmd[cmd.index("gvf-run") + 1]
    run_dir = output / gene / "run1"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def _extract(tmp_path: Path, *, fast: bool = False) -> benchmark.ExtractionBatch:
    return benchmark.do_extract(
        "test@example.com",
        tmp_path / "extract",
        ["KCNH2"],
        fast=fast,
        source_recovery_timeout_s=120,
        use_local_source_corpus=False,
    )


def test_exit_three_warning_status_is_scoreable_and_degraded(tmp_path, monkeypatch):
    def fake_run(cmd, **kwargs):
        run_dir = _run_dir_from_cmd(cmd)
        db = run_dir / "KCNH2.final.db"
        db.write_bytes(b"sqlite")
        # Deliberately omit severity: pre-schema statuses must still work.
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "gene": "KCNH2",
                    "status": "completed_with_warnings",
                    "exit_code": 3,
                    "stage_failures": ["source recovery failed"],
                    "active_db": db.name,
                }
            )
        )
        return SimpleNamespace(returncode=3)

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)

    result = _extract(tmp_path)

    assert result.dbs["KCNH2"].name == "KCNH2.final.db"
    assert result.degraded == {"KCNH2": "source recovery failed"}
    assert result.failures == {}


def test_nonzero_exit_never_consumes_stale_status_or_db(tmp_path, monkeypatch):
    run_dir = tmp_path / "extract" / "KCNH2" / "KCNH2" / "old"
    run_dir.mkdir(parents=True)
    (run_dir / "KCNH2.db").write_bytes(b"old")
    (run_dir / "RUN_STATUS.json").write_text(
        json.dumps(
            {
                "gene": "KCNH2",
                "status": "completed_with_warnings",
                "exit_code": 3,
            }
        )
    )
    monkeypatch.setattr(
        benchmark.subprocess,
        "run",
        lambda cmd, **kwargs: SimpleNamespace(returncode=3),
    )

    result = _extract(tmp_path)

    assert result.dbs == {}
    assert "without a fresh RUN_STATUS.json" in result.failures["KCNH2"]


def test_fatal_fresh_status_is_not_scoreable(tmp_path, monkeypatch):
    def fake_run(cmd, **kwargs):
        run_dir = _run_dir_from_cmd(cmd)
        (run_dir / "KCNH2.db").write_bytes(b"sqlite")
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "gene": "KCNH2",
                    "status": "failed",
                    "severity": "fatal",
                    "exit_code": 1,
                    "active_db": "KCNH2.db",
                }
            )
        )
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)

    result = _extract(tmp_path)

    assert result.dbs == {}
    assert "unscoreable" in result.failures["KCNH2"]


def test_exit_zero_without_status_uses_only_fresh_db_and_marks_degraded(
    tmp_path, monkeypatch
):
    def fake_run(cmd, **kwargs):
        run_dir = _run_dir_from_cmd(cmd)
        (run_dir / "KCNH2.db").write_bytes(b"sqlite")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)

    result = _extract(tmp_path)

    assert result.dbs["KCNH2"].name == "KCNH2.db"
    assert "no fresh RUN_STATUS.json" in result.degraded["KCNH2"]


def test_active_db_path_takes_precedence_over_other_fresh_dbs(tmp_path, monkeypatch):
    def fake_run(cmd, **kwargs):
        run_dir = _run_dir_from_cmd(cmd)
        active = run_dir / "chosen" / "KCNH2.active.db"
        active.parent.mkdir()
        active.write_bytes(b"active")
        (run_dir / "KCNH2.newer.db").write_bytes(b"other")
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "gene": "KCNH2",
                    "status": "completed",
                    "severity": "ok",
                    "exit_code": 0,
                    "active_db": "chosen/KCNH2.active.db",
                }
            )
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)

    result = _extract(tmp_path)

    assert result.dbs["KCNH2"].name == "KCNH2.active.db"
    assert result.degraded == {}
    assert result.failures == {}


def test_status_exit_must_match_subprocess_exit():
    scoreable, degraded, failure = benchmark._classify_status(
        {
            "status": "completed_with_warnings",
            "severity": "warning",
            "exit_code": 3,
        },
        actual_exit=0,
    )

    assert scoreable is False
    assert degraded is None
    assert "disagrees" in failure


def test_fresh_status_requires_explicit_active_db(tmp_path, monkeypatch):
    def fake_run(cmd, **kwargs):
        run_dir = _run_dir_from_cmd(cmd)
        (run_dir / "KCNH2.intermediate.db").write_bytes(b"not final")
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "gene": "KCNH2",
                    "status": "completed_with_warnings",
                    "severity": "warning",
                    "exit_code": 3,
                    "stage_failures": ["source recovery failed"],
                }
            )
        )
        return SimpleNamespace(returncode=3)

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)

    result = _extract(tmp_path)

    assert result.dbs == {}
    assert "no active_db" in result.failures["KCNH2"]


def test_fast_extract_is_always_degraded(tmp_path, monkeypatch):
    def fake_run(cmd, **kwargs):
        run_dir = _run_dir_from_cmd(cmd)
        db = run_dir / "KCNH2.db"
        db.write_bytes(b"sqlite")
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "gene": "KCNH2",
                    "status": "completed",
                    "severity": "ok",
                    "exit_code": 0,
                    "stage_failures": [],
                    "active_db": db.name,
                }
            )
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)

    result = _extract(tmp_path, fast=True)

    assert "--fast disables source recovery" in result.degraded["KCNH2"]
