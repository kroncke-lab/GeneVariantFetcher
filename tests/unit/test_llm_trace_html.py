from __future__ import annotations

import json
import re
from pathlib import Path
from types import SimpleNamespace

import pytest

from utils.llm_trace_html import (
    build_trace_html_report,
    collect_trace_report_data,
    render_trace_report,
)
from utils.llm_trace import (
    TRACE_MANIFEST_NAME,
    build_trace_manifest,
    capture_llm_call,
    configure_llm_tracing,
    llm_trace_scope,
    record_trace_event,
    reset_llm_tracing,
)


@pytest.fixture(autouse=True)
def _reset_trace_configuration():
    reset_llm_tracing()
    yield
    reset_llm_tracing()


def _traced_run(tmp_path: Path) -> tuple[Path, Path]:
    run_dir = tmp_path / "run"
    trace_root = configure_llm_tracing(
        run_dir / "llm_traces", run_id="trace-report-fixture"
    )
    prompt = "Target gene: KCNQ1\nPMID: 12345678\nCheck </script> safely."
    response = SimpleNamespace(
        output_text=json.dumps(
            {
                "variants": [{"variant": "c.1A>T"}],
                "curation_rationale": "The table reports two carriers.",
            }
        ),
        usage=SimpleNamespace(input_tokens=14, output_tokens=9, total_tokens=23),
        model_dump=lambda mode=None: {
            "id": "response-1",
            "output_text": '{"variants":[{"variant":"c.1A>T"}]}',
            "usage": {"input_tokens": 14, "output_tokens": 9, "total_tokens": 23},
        },
    )
    with llm_trace_scope(
        gene="KCNQ1",
        pmid="12345678",
        stage="paper_curation",
        component="fixture",
    ):
        _, call_ref = capture_llm_call(
            provider="fixture",
            requested_model="model-a",
            resolved_model="provider/model-a",
            request={"model": "model-a", "input": prompt},
            call=lambda: response,
        )
        record_trace_event(
            "paper_curation_decision",
            {
                "selected_model": "model-a",
                "variant_count": 1,
                "curation_rationale": "The table reports two carriers.",
                "accepted_response_trace_id": call_ref["trace_id"],
                "selection_policy": "Use the final complete JSON response.",
            },
        )
    build_trace_manifest(
        trace_root,
        output_path=trace_root / TRACE_MANIFEST_NAME,
        run_id="trace-report-fixture",
    )
    (run_dir / "selection.json").write_text(
        json.dumps(
            {
                "papers": [
                    {
                        "gene": "KCNQ1",
                        "pmid": "12345678",
                        "title": "A traced channelopathy paper",
                        "source_selection": {
                            "selected": "/corpus/KCNQ1/12345678_FULL_CONTEXT.md",
                            "candidates": [
                                {"path": "/corpus/KCNQ1/12345678_FULL_CONTEXT.md"},
                                {"path": "/corpus/KCNQ1/12345678_CLEANED.md"},
                            ],
                            "rationale": "Selected the richer complete rendering.",
                        },
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    return run_dir, trace_root


def test_report_groups_call_and_decision_by_paper(tmp_path: Path):
    run_dir, trace_root = _traced_run(tmp_path)

    data = collect_trace_report_data(trace_root, run_dir=run_dir)

    assert data["integrity"]["valid"] is True
    assert data["summary"] == {
        "paper_count": 1,
        "group_count": 1,
        "trace_count": 2,
        "llm_call_count": 1,
        "decision_event_count": 1,
        "failure_count": 0,
        "token_count": 23,
    }
    (paper,) = data["papers"]
    assert paper["gene"] == "KCNQ1"
    assert paper["pmid"] == "12345678"
    assert paper["title"] == "A traced channelopathy paper"
    assert paper["metadata"]["source_selection"]["rationale"].startswith(
        "Selected the richer"
    )
    assert [record["record_type"] for record in paper["records"]] == [
        "llm_call",
        "decision_event",
    ]


def test_report_is_self_contained_parseable_and_script_safe(tmp_path: Path):
    run_dir, trace_root = _traced_run(tmp_path)
    data = collect_trace_report_data(trace_root, run_dir=run_dir)

    rendered = render_trace_report(data)

    assert "Search gene, PMID, title" in rendered
    assert "Sent to model" in rendered
    assert "Why this outcome" in rendered
    assert "Source choice" in rendered
    assert "Open accepted model call" in rendered
    assert "Exact trace record" in rendered
    assert "Check </script> safely." not in rendered
    assert "Check <\\/script> safely." in rendered
    # The report is normally opened from disk. Restricted file-origin clipboard
    # or storage permissions must not abort initialization or leave Copy inert.
    assert (
        "navigator.clipboard.writeText(text).then(done).catch(copyWithTextarea)"
        in rendered
    )
    assert 'safeStorageGet("gvf-trace-theme")' in rendered
    assert 'safeStorageSet("gvf-trace-theme", next)' in rendered
    match = re.search(r"const DATA = (.*?);\n    const state", rendered, re.S)
    assert match
    embedded = json.loads(match.group(1))
    assert embedded["papers"][0]["records"][0]["request"]["payload"]["input"].endswith(
        "safely."
    )

    output = run_dir / "llm_trace_report.html"
    built = build_trace_html_report(
        trace_root,
        output_path=output,
        run_dir=run_dir,
        title="Fixture trace review",
    )
    assert output.is_file()
    assert "Fixture trace review" in output.read_text(encoding="utf-8")
    assert built["summary"]["trace_count"] == 2


def test_report_surfaces_manifest_tampering(tmp_path: Path):
    run_dir, trace_root = _traced_run(tmp_path)
    trace_path = next(
        path for path in trace_root.rglob("*.json") if path.name != TRACE_MANIFEST_NAME
    )
    trace_path.write_text(trace_path.read_text(encoding="utf-8") + " ")

    data = collect_trace_report_data(trace_root, run_dir=run_dir)

    assert data["integrity"]["valid"] is False
    assert any(
        "changed after manifest" in error for error in data["integrity"]["errors"]
    )


def test_report_names_its_integrity_level_and_run_id(tmp_path: Path):
    """A manifest generated by the report is 'generated_now', never 'verified'."""
    run_dir, trace_root = _traced_run(tmp_path)

    verified = collect_trace_report_data(
        trace_root, run_dir=run_dir, run_id="trace-report-fixture"
    )
    assert verified["integrity"]["level"] == "write_time_verified"
    assert verified["integrity"]["manifest_generated_by_this_report"] is False
    assert verified["run_id"] == "trace-report-fixture"
    rendered = render_trace_report(verified)
    assert "Verified against write-time digests." in rendered
    assert "run ${DATA.run_id}" in rendered or "DATA.run_id" in rendered

    # Remove the manifest so the report has to generate one itself.
    (trace_root / TRACE_MANIFEST_NAME).unlink()
    generated = collect_trace_report_data(trace_root, run_dir=run_dir)
    assert generated["integrity"]["level"] == "generated_now"
    assert generated["integrity"]["manifest_generated_by_this_report"] is True
    assert "not a tamper guarantee" in render_trace_report(generated)


def test_report_lists_stage_coverage_and_missing_decision_links(tmp_path: Path):
    run_dir = tmp_path / "run"
    trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="coverage")
    # A stage with model calls and NO decision event is a real adjudicability gap.
    with llm_trace_scope(
        gene="KCNH2", pmid="111", stage="paper_final_check", component="fixture"
    ):
        capture_llm_call(
            provider="fixture",
            requested_model="gpt-fixture",
            resolved_model="gpt-fixture",
            request={"input": "check this paper"},
            call=lambda: SimpleNamespace(output_text="{}"),
        )

    data = collect_trace_report_data(trace_root, run_dir=run_dir)

    stages = {entry["stage"]: entry for entry in data["coverage"]["stages"]}
    assert stages["paper_final_check"]["llm_calls"] == 1
    assert stages["paper_final_check"]["decisions"] == 0
    (gap,) = data["coverage"]["missing_decision_links"]
    assert gap["stage"] == "paper_final_check"
    assert gap["expected_event"] == "paper_final_check_decision"
    assert "Route coverage" in render_trace_report(data)


def test_long_bodies_are_bounded_and_the_omission_is_visible(tmp_path: Path):
    run_dir = tmp_path / "run"
    trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="bounded")
    huge = "PROMPT-" + ("x" * 60_000) + "-TAIL"
    with llm_trace_scope(gene="KCNH2", pmid="222", stage="paper_variant_extraction"):
        capture_llm_call(
            provider="fixture",
            requested_model="model-a",
            resolved_model="model-a",
            request={"input": huge},
            call=lambda: SimpleNamespace(output_text="ok"),
        )

    data = collect_trace_report_data(trace_root, run_dir=run_dir, max_field_chars=1000)

    embedded = data["papers"][0]["records"][0]["request"]["payload"]["input"]
    assert len(embedded) < 2000
    assert embedded.startswith("PROMPT-")
    # Never silent: the marker states the extent and names the record on disk.
    assert "truncated" in embedded
    assert data["papers"][0]["records"][0]["_report_path"] in embedded
    (omission,) = [o for o in data["omissions"] if o["kind"] == "body_truncated"]
    assert omission["characters_total"] > 60_000
    assert omission["characters_embedded"] == 1000
    assert omission["sha256"]
    rendered = render_trace_report(data)
    assert "Omitted from this page" in rendered
    assert "-TAIL" not in rendered  # the dropped tail really is not embedded


def test_large_runs_shard_per_paper_with_a_linked_index(tmp_path: Path):
    run_dir = tmp_path / "run"
    trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="sharded")
    for index in range(5):
        with llm_trace_scope(
            gene="KCNH2", pmid=f"90000{index}", stage="paper_variant_extraction"
        ):
            capture_llm_call(
                provider="fixture",
                requested_model="model-a",
                resolved_model="model-a",
                request={"input": f"paper {index}"},
                call=lambda: SimpleNamespace(output_text="ok"),
            )

    output = run_dir / "llm_trace_report.html"
    data = build_trace_html_report(
        trace_root, output_path=output, run_dir=run_dir, max_papers_per_file=2
    )

    assert data["sharded"] is True
    shard_dir = Path(data["shard_dir"])
    shards = sorted(shard_dir.glob("*.html"))
    assert len(shards) == 5
    index_html = output.read_text(encoding="utf-8")
    assert '"mode":"index"' in index_html
    for shard in shards:
        assert f"{shard_dir.name}/{shard.name}" in index_html
        body = shard.read_text(encoding="utf-8")
        assert '"mode":"shard"' in body
        # Each shard links back and embeds exactly one paper.
        assert "all papers in this run" in body


def test_accepted_and_discarded_calls_are_labelled(tmp_path: Path):
    run_dir = tmp_path / "run"
    trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="attempts")
    ids = []
    with llm_trace_scope(gene="KCNH2", pmid="333", stage="paper_variant_extraction"):
        for attempt, role in ((1, "primary"), (2, "json_repair")):
            with llm_trace_scope(attempt=attempt, operation=role):
                _, summary = capture_llm_call(
                    provider="fixture",
                    requested_model="model-a",
                    resolved_model="model-a",
                    request={"input": f"attempt {attempt}"},
                    call=lambda: SimpleNamespace(output_text="ok"),
                )
            ids.append(summary["trace_id"])
        record_trace_event(
            "paper_extraction_selection",
            {
                "accepted_response_trace_id": ids[1],
                "discarded_trace_ids": [ids[0]],
                "repaired": True,
                "attempt_trace_links": [
                    {
                        "attempt": 1,
                        "role": "primary",
                        "trace_id": ids[0],
                        "outcome": "parse_failed",
                    },
                    {
                        "attempt": 2,
                        "role": "json_repair",
                        "trace_id": ids[1],
                        "outcome": "accepted",
                    },
                ],
            },
        )

    data = collect_trace_report_data(trace_root, run_dir=run_dir)

    (paper,) = data["papers"]
    assert paper["accepted_trace_ids"] == [ids[1]]
    assert paper["discarded_trace_ids"] == [ids[0]]
    by_id = {record.get("trace_id"): record for record in paper["records"]}
    assert by_id[ids[0]]["_is_discarded"] is True
    assert by_id[ids[0]]["_is_accepted"] is False
    assert by_id[ids[1]]["_is_accepted"] is True
    assert by_id[ids[1]]["_attempt_role"] == "json_repair"
    rendered = render_trace_report(data)
    assert "Retries &amp; repairs" in rendered
    assert "Open accepted model call" in rendered


def test_production_source_completeness_supplies_paper_metadata(tmp_path: Path):
    """docs' adjudication step 2 must apply to production, not only benchmarks."""
    run_dir = tmp_path / "run"
    trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="prod")
    with llm_trace_scope(gene="KCNH2", pmid="444", stage="paper_variant_extraction"):
        capture_llm_call(
            provider="fixture",
            requested_model="model-a",
            resolved_model="model-a",
            request={"input": "x"},
            call=lambda: SimpleNamespace(output_text="ok"),
        )
    (run_dir / "source_completeness.json").write_text(
        json.dumps(
            {
                "papers": [
                    {
                        "gene": "KCNH2",
                        "pmid": "444",
                        "title": "A production paper",
                        "status": "usable_fulltext",
                        "source_file": "pmc_fulltext/444_FULL_CONTEXT.md",
                        "chars": 51783,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    data = collect_trace_report_data(trace_root, run_dir=run_dir)

    (paper,) = data["papers"]
    assert paper["title"] == "A production paper"
    assert paper["metadata"]["source_status"] == "usable_fulltext"
    assert paper["metadata"]["source_chars"] == 51783
    assert "Source read" in render_trace_report(data)
