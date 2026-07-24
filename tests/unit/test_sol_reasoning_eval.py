from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from types import ModuleType

import pytest

from benchmarks.curated_extraction_eval import run_sol_reasoning_eval as runner


MANIFEST_COLUMNS = [
    "gene",
    "pmid",
    "strategy",
    "gold_variant_rows",
    "gold_carriers",
    "gold_affected",
    "gold_unaffected",
    "gold_source",
    "title",
    "why_selected",
    "corpus_source",
    "source_status",
    "pubmed_url",
]


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {column: row.get(column, "") for column in MANIFEST_COLUMNS}
            )


def _fake_sol_module(tmp_path: Path, calls: list[dict] | None = None) -> ModuleType:
    module_file = tmp_path / "fake_sol_extractor.py"
    module_file.write_text("# deterministic fake module\n", encoding="utf-8")
    module = ModuleType("fake_sol_extractor")
    module.__file__ = str(module_file)
    module.SOL_PROTOCOL_VERSION = "test-v1"

    def strip_irrelevant_context(text: str) -> str:
        return text.replace("## References\nignore me", "")

    def extract_paper(**kwargs):
        if calls is not None:
            calls.append(kwargs)
        cleaned = strip_irrelevant_context(kwargs["source_text"])
        return (
            {
                "paper_metadata": {
                    "pmid": kwargs["pmid"],
                    "title": kwargs.get("title") or "Synthetic paper",
                    "gene_symbol": kwargs["gene"],
                },
                "variants": [],
                "extraction_metadata": {
                    "total_variants_found": 0,
                    "model_used": kwargs["model"],
                },
            },
            {
                "model": kwargs["model"],
                "requested_reasoning_effort": kwargs["reasoning_effort"],
                "effective_reasoning_effort": kwargs["reasoning_effort"],
                "source_sha256": kwargs["source_sha256"],
                "cleaned_source_sha256": hashlib.sha256(
                    cleaned.encode("utf-8")
                ).hexdigest(),
                "usage": {"input_tokens": 11, "output_tokens": 7},
                "tool_calls": 2,
            },
        )

    module.strip_irrelevant_context = strip_irrelevant_context
    module.extract_paper = extract_paper
    return module


def _synthetic_fixture(tmp_path: Path) -> tuple[Path, Path, dict[str, dict[str, str]]]:
    manifest_path = tmp_path / "manifest.csv"
    sources_dir = tmp_path / "sources"
    source = sources_dir / "TEST" / "12345678" / "12345678_FULL_CONTEXT.md"
    source.parent.mkdir(parents=True)
    source.write_text(
        "# Synthetic study\nOne human carrier.\n## References\nignore me",
        encoding="utf-8",
    )
    row = {
        "gene": "TEST",
        "pmid": "12345678",
        "strategy": "text",
        "title": "Synthetic paper",
        "corpus_source": "missing/fallback.md",
    }
    _write_manifest(manifest_path, [row])
    return manifest_path, sources_dir, {"TEST:12345678": row}


def test_parse_efforts_keeps_max_literal_and_rejects_aliases() -> None:
    assert runner.parse_efforts("none,xhigh,max") == ("none", "xhigh", "max")
    with pytest.raises(runner.SolEvalError, match="unsupported"):
        runner.parse_efforts("minimal")
    with pytest.raises(runner.SolEvalError, match="duplicate"):
        runner.parse_efforts("max,max")


def test_resolve_source_prefers_frozen_source_then_manifest_fallback(
    tmp_path: Path,
) -> None:
    sources = tmp_path / "sources"
    frozen = sources / "GENE" / "111" / "111_FULL_CONTEXT.md"
    frozen.parent.mkdir(parents=True)
    frozen.write_text("frozen", encoding="utf-8")
    fallback = tmp_path / "fallback.md"
    fallback.write_text("fallback", encoding="utf-8")
    row = {"gene": "GENE", "pmid": "111", "corpus_source": "fallback.md"}

    assert runner.resolve_source(row, sources_dir=sources, repo=tmp_path) == frozen
    frozen.unlink()
    assert runner.resolve_source(row, sources_dir=sources, repo=tmp_path) == fallback


def test_interleaved_schedule_rotates_effort_order() -> None:
    def paper(key: str) -> runner.PaperInput:
        gene, pmid = key.split(":")
        return runner.PaperInput(
            key=key,
            gene=gene,
            pmid=pmid,
            title="",
            strategy="text",
            source_path="/tmp/source",
            source_bytes=1,
            source_chars=1,
            source_sha256="a" * 64,
            cleaned_chars=1,
            cleaned_source_sha256="b" * 64,
        )

    schedule = runner.interleaved_schedule(
        [paper("A:1"), paper("B:2")], ("none", "high", "max")
    )
    assert [(task.paper.key, task.effort) for task in schedule] == [
        ("A:1", "none"),
        ("A:1", "high"),
        ("A:1", "max"),
        ("B:2", "high"),
        ("B:2", "max"),
        ("B:2", "none"),
    ]


def test_success_receipt_resumes_only_when_all_hashes_match(tmp_path: Path) -> None:
    _manifest, sources, selected = _synthetic_fixture(tmp_path)
    calls: list[dict] = []
    sol_module = _fake_sol_module(tmp_path, calls)
    paper = runner.build_source_manifest(
        selected, sources_dir=sources, sol_module=sol_module, repo=tmp_path
    )[0]
    task = runner.TaskSpec(ordinal=0, paper=paper, effort="max")
    run_dir = tmp_path / "run"

    first = runner._execute_task(
        task,
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )
    second = runner._execute_task(
        task,
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )
    assert first.status == second.status == "success"
    assert second.reused is True
    assert len(calls) == 1

    extraction_path = runner._task_paths(run_dir, task)["extraction"]
    extraction_path.write_text("{}\n", encoding="utf-8")
    with pytest.raises(runner.ResumeIdentityError, match="receipt hash mismatch"):
        runner._execute_task(
            task,
            run_dir=run_dir,
            model=runner.DEFAULT_MODEL,
            config_sha="config-hash",
            sol_module=sol_module,
        )


def test_effective_max_cannot_be_reported_as_xhigh(tmp_path: Path) -> None:
    _manifest, sources, selected = _synthetic_fixture(tmp_path)
    sol_module = _fake_sol_module(tmp_path)
    original_extract = sol_module.extract_paper

    def aliased_extract(**kwargs):
        extraction, telemetry = original_extract(**kwargs)
        telemetry["effective_reasoning_effort"] = "xhigh"
        return extraction, telemetry

    sol_module.extract_paper = aliased_extract
    paper = runner.build_source_manifest(
        selected, sources_dir=sources, sol_module=sol_module, repo=tmp_path
    )[0]
    task = runner.TaskSpec(ordinal=0, paper=paper, effort="max")
    result = runner._execute_task(
        task,
        run_dir=tmp_path / "run",
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )

    assert result.status == "failed"
    assert "effective_reasoning_effort mismatch" in (result.error or "")
    assert not runner._task_paths(tmp_path / "run", task)["receipt"].exists()


def test_nested_usage_keeps_reasoning_and_cache_dimensions() -> None:
    usage = runner._token_usage(
        {
            "responses": {
                "usage": {
                    "input_tokens": 100,
                    "cached_input_tokens": 40,
                    "cache_write_input_tokens": 12,
                    "output_tokens": 30,
                    "reasoning_tokens": 20,
                    "total_tokens": 130,
                }
            }
        }
    )
    assert usage == {
        "input_tokens": 100,
        "cached_input_tokens": 40,
        "cache_write_input_tokens": 12,
        "output_tokens": 30,
        "reasoning_tokens": 20,
        "total_tokens": 130,
    }
    # Reasoning is diagnostic and already included in output billing.
    assert runner._estimated_cost(
        usage,
        runner.PriceConfig(
            input_per_million=1,
            cached_input_per_million=0.5,
            output_per_million=10,
        ),
    ) == pytest.approx((60 + 20 + 300) / 1_000_000)


def test_model_loop_failure_is_separate_from_migration_success(tmp_path: Path) -> None:
    _manifest, sources, selected = _synthetic_fixture(tmp_path)
    sol_module = _fake_sol_module(tmp_path)
    original_extract = sol_module.extract_paper

    def empty_failure(**kwargs):
        extraction, telemetry = original_extract(**kwargs)
        telemetry.pop("tool_calls", None)
        extraction["extraction_metadata"]["responses"] = {
            "ok": False,
            "status": "http_error",
            "error": "provider rejected request",
            "telemetry": {"tool_calls_executed": 3, "tool_errors": 1},
        }
        return extraction, telemetry

    sol_module.extract_paper = empty_failure
    paper = runner.build_source_manifest(
        selected, sources_dir=sources, sol_module=sol_module, repo=tmp_path
    )[0]
    task = runner.TaskSpec(ordinal=0, paper=paper, effort="none")
    run_dir = tmp_path / "run"
    result = runner._execute_task(
        task,
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )
    rows = runner._effort_telemetry_rows("none", [paper], run_dir)

    assert result.status == "failed"
    assert rows[0]["status"] == "failed"
    assert rows[0]["model_tool_loop_failed"] is True
    assert rows[0]["model_status"] == "http_error"
    assert rows[0]["tool_calls"] == 3
    assert rows[0]["input_tokens"] == 11
    assert rows[0]["output_tokens"] == 7
    paths = runner._task_paths(run_dir, task)
    assert paths["extraction"].is_file()
    assert paths["telemetry"].is_file()
    assert not paths["receipt"].exists()
    summary = runner._summarize_effort(
        "none",
        rows,
        prices=runner.PriceConfig(),
        score_summary=None,
        scoring_reason="test",
    )
    assert summary["papers_successful"] == 0
    assert summary["papers_failed"] == 1
    assert summary["model_tool_loop_failures"] == 1


def test_model_loop_failure_without_receipt_is_retried(tmp_path: Path) -> None:
    _manifest, sources, selected = _synthetic_fixture(tmp_path)
    calls: list[dict] = []
    sol_module = _fake_sol_module(tmp_path, calls)
    original_extract = sol_module.extract_paper

    def fail_once(**kwargs):
        extraction, telemetry = original_extract(**kwargs)
        if len(calls) == 1:
            extraction.setdefault("extraction_metadata", {})["responses"] = {
                "ok": False,
                "status": "incomplete",
                "error": "temporary incomplete response",
            }
        return extraction, telemetry

    sol_module.extract_paper = fail_once
    paper = runner.build_source_manifest(
        selected, sources_dir=sources, sol_module=sol_module, repo=tmp_path
    )[0]
    task = runner.TaskSpec(ordinal=0, paper=paper, effort="high")
    run_dir = tmp_path / "run"

    first = runner._execute_task(
        task,
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )
    second = runner._execute_task(
        task,
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )

    assert first.status == "failed"
    assert second.status == "success"
    assert second.reused is False
    assert len(calls) == 2
    assert runner._task_paths(run_dir, task)["receipt"].is_file()


def test_custom_manifest_run_never_scores_and_uses_no_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, sources, _selected = _synthetic_fixture(tmp_path)
    sol_module = _fake_sol_module(tmp_path)
    monkeypatch.setattr(runner, "_load_sol_module", lambda: sol_module)

    def scorer_must_not_run(*args, **kwargs):
        raise AssertionError("subset/custom manifest must not be scored")

    monkeypatch.setattr(runner, "score_complete_effort", scorer_must_not_run)
    outdir = tmp_path / "out"
    exit_code = runner.main(
        [
            "--manifest",
            str(manifest),
            "--sources-dir",
            str(sources),
            "--efforts",
            "none",
            "--run-id",
            "synthetic",
            "--outdir",
            str(outdir),
        ]
    )

    assert exit_code == 0
    comparison = json.loads(
        (outdir / "synthetic" / "comparison.json").read_text(encoding="utf-8")
    )
    result = comparison["efforts"]["none"]
    assert result["papers_total"] == 1
    assert result["papers_successful"] == 1
    assert result["scored"] is False
    assert "custom manifest" in result["scoring_skipped_reason"]


def test_dry_run_writes_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, sources, _selected = _synthetic_fixture(tmp_path)
    sol_module = _fake_sol_module(tmp_path)
    monkeypatch.setattr(runner, "_load_sol_module", lambda: sol_module)
    outdir = tmp_path / "never-created"

    assert (
        runner.main(
            [
                "--manifest",
                str(manifest),
                "--sources-dir",
                str(sources),
                "--efforts",
                "none,max",
                "--run-id",
                "dry",
                "--outdir",
                str(outdir),
                "--dry-run",
            ]
        )
        == 0
    )
    assert not outdir.exists()


def test_resume_rejects_config_change(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    kwargs = {
        "run_dir": run_dir,
        "run_id": "identity",
        "config": {"model": "sol", "efforts": ["none"]},
        "source_records": [{"key": "G:1", "source_sha256": "a"}],
        "schedule_records": [{"ordinal": 0, "key": "G:1", "effort": "none"}],
    }
    runner._prepare_run(**kwargs)
    runner._prepare_run(**kwargs)
    with pytest.raises(runner.ResumeIdentityError, match="config_manifest"):
        runner._prepare_run(**{**kwargs, "config": {"model": "other"}})


def test_database_migration_is_complete_and_integral(tmp_path: Path) -> None:
    _manifest, sources, selected = _synthetic_fixture(tmp_path)
    sol_module = _fake_sol_module(tmp_path)
    paper = runner.build_source_manifest(
        selected, sources_dir=sources, sol_module=sol_module, repo=tmp_path
    )[0]
    task = runner.TaskSpec(ordinal=0, paper=paper, effort="none")
    run_dir = tmp_path / "run"
    result = runner._execute_task(
        task,
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
        sol_module=sol_module,
    )
    assert result.status == "success"

    dbs, manifest = runner.build_gene_databases(
        "none",
        [paper],
        run_dir=run_dir,
        model=runner.DEFAULT_MODEL,
        config_sha="config-hash",
    )
    assert dbs["TEST"].is_file()
    assert manifest["databases"]["TEST"]["papers_expected"] == 1
    assert manifest["databases"]["TEST"]["papers_observed"] == 1
    assert manifest["databases"]["TEST"]["integrity_check"] == ["ok"]
