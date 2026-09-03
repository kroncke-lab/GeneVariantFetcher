"""Regression coverage for the extraction-blinded Codex paper evaluation."""

from __future__ import annotations

import importlib.util
import json
import sqlite3
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from benchmarks.codex_paper_eval.build_report_artifact import build_payload
from pipeline.prompts import TABLE_ATTRIBUTION_GUIDANCE
from utils.llm_trace import (
    build_trace_manifest,
    capture_llm_call,
    configure_llm_tracing,
    llm_trace_scope,
    reset_llm_tracing,
)
import benchmarks.codex_paper_eval.run_eval as run_eval_module
from benchmarks.codex_paper_eval.run_eval import (
    CARDIAC_GENES,
    DEFAULT_GOLD,
    EXTRACTION_INSTRUCTIONS,
    aggregate,
    choose_source,
    compare_paired_reports,
    command_extract,
    command_lock,
    command_prepare,
    digest,
    digest_gold_root,
    effective_effort,
    gold_count_eligible_pmids,
    gold_variant_eligible_pmids,
    gold_csv_path,
    load_gold,
    matches,
    material_digest_errors,
    production_run_status_lock_entries,
    production_trace_lock_entries,
    read_paper_manifest,
    record_tranche_consumption,
    reasoning_params,
    score_one,
    score_prediction_lane,
    looks_truncated_json,
    selection_metadata,
    supports_images,
    usable_sources,
    validated_score_gold_root,
    validate_predictions,
    write_json,
    write_markdown_report,
)


COUNT_FIELDS = ("carriers", "affected", "unaffected")
GENES = ("SCN5A", "KCNH2", "KCNQ1", "RYR2")


@pytest.fixture(autouse=True)
def _reset_trace_configuration():
    reset_llm_tracing()
    yield
    reset_llm_tracing()


def _count_metric(asserted: int = 1, predicted: int = 1) -> dict:
    return {
        "gold_asserted": asserted,
        "predicted": predicted,
        "recall": predicted / asserted if asserted else None,
        "mae": 0.0 if predicted else None,
        "rmse": 0.0 if predicted else None,
    }


def test_load_gold_prefers_adjudicated_v2_counts_and_explicit_nulls(tmp_path):
    gold = tmp_path / "SCN5A_recall_input.csv"
    gold.write_text(
        "variant,pmid,carriers,affected,unaffected,gold_v2_carriers,"
        "gold_v2_affected,gold_v2_unaffected,gold_v2_status\n"
        "S1103Y,20470418,85,39,46,26,17,9,adjudicated_variant_carrier_count\n"
        "A1V,20470418,2,1,1,1,1,,adjudicated_null_unaffected\n"
        "G2D,20470418,112,112,0,,,,excluded_duplicate_current_cohort\n"
        "R2W,20470418,3,2,1,,,,\n"
        "V3A,20470418,4,3,1\n"
    )

    assert load_gold(tmp_path, "SCN5A", "20470418") == [
        {"variant": "S1103Y", "carriers": 26, "affected": 17, "unaffected": 9},
        {"variant": "A1V", "carriers": 1, "affected": 1, "unaffected": None},
        {"variant": "R2W", "carriers": 3, "affected": 2, "unaffected": 1},
        {"variant": "V3A", "carriers": 4, "affected": 3, "unaffected": 1},
    ]


def test_counted_precision_metrics_keep_the_two_denominators_distinct():
    prediction = {
        "variants": [
            {"variant": "A1V", "carriers": None},
            {"variant": "A2V", "carriers": 1},
            {"variant": "A3V", "carriers": None},
            {"variant": "A4V", "carriers": 1},
        ]
    }
    gold = [
        {"variant": "A1V", "carriers": 1},
        {"variant": "A2V", "carriers": 1},
    ]

    score = score_one("KCNH2", "1", prediction, gold)

    assert score["counted_precision"] == {
        "matched_rows": 2,
        "counted_extra_rows": 1,
        "precision_vs_counted_gold_pmids": pytest.approx(2 / 3),
        "count_bearing_matched_rows": 1,
        "precision_among_count_bearing_predictions": pytest.approx(1 / 2),
    }
    combined = aggregate([score])
    assert combined["counted_precision"]["precision_vs_counted_gold_pmids"] == (
        pytest.approx(2 / 3)
    )
    assert combined["counted_precision"][
        "precision_among_count_bearing_predictions"
    ] == pytest.approx(1 / 2)


def test_load_gold_rejects_unknown_v2_status(tmp_path):
    gold = tmp_path / "SCN5A_recall_input.csv"
    gold.write_text(
        "variant,pmid,carriers,affected,unaffected,gold_v2_carriers,"
        "gold_v2_affected,gold_v2_unaffected,gold_v2_status\n"
        "A1V,1,2,1,1,,,,needs_review\n"
    )

    with pytest.raises(ValueError, match="Unknown gold_v2_status"):
        load_gold(tmp_path, "SCN5A", "1")


def test_gold_count_eligibility_respects_v2_nulls_and_exclusions(tmp_path):
    gold = tmp_path / "SCN5A_recall_input.csv"
    gold.write_text(
        "variant,pmid,carriers,affected,unaffected,gold_v2_carriers,"
        "gold_v2_affected,gold_v2_unaffected,gold_v2_status\n"
        "A1V,1,2,1,1,2,1,,adjudicated_null_unaffected\n"
        "A2V,2,2,1,1,2,1,1,excluded_duplicate_current_cohort\n"
        "A3V,3,2,1,1,,,,\n"
        "A4V,4,9,9,0,3,2,1,adjudicated_variant_carrier_count\n"
    )

    assert gold_count_eligible_pmids(tmp_path, "SCN5A") == {"3", "4"}
    assert gold_variant_eligible_pmids(tmp_path, "SCN5A") == {"1", "3", "4"}


def test_provenance_lanes_score_paper_discovery_separately_from_linkage(tmp_path):
    (tmp_path / "KCNH2_recall_input.csv").write_text(
        "variant,pmid,carriers,affected,unaffected\nA1V,1,1,1,0\nA2V,1,1,1,0\n"
    )
    selection = {"papers": [{"gene": "KCNH2", "pmid": "1"}]}
    paper_row = {
        "variant": "A1V",
        "carriers": 1,
        "affected": 1,
        "unaffected": 0,
    }
    linkage_row = {
        "variant": "A2V",
        "carriers": None,
        "affected": None,
        "unaffected": None,
    }
    predictions = {
        "primary_score_lane": "paper_derived",
        "papers": [
            {
                "gene": "KCNH2",
                "pmid": "1",
                "variants": [paper_row],
                "comparison_variants": {"linkage_assisted": [paper_row, linkage_row]},
            }
        ],
    }

    paper_scores, _, _ = score_prediction_lane(
        selection, predictions, tmp_path, "paper_derived"
    )
    linkage_scores, _, _ = score_prediction_lane(
        selection, predictions, tmp_path, "linkage_assisted"
    )

    assert aggregate(paper_scores)["recall"] == pytest.approx(0.5)
    assert aggregate(linkage_scores)["recall"] == pytest.approx(1.0)


def test_scoring_burn_ledger_is_append_only_and_idempotent(tmp_path: Path):
    registry = tmp_path / "suite" / "registry.json"
    registry.parent.mkdir()
    registry.write_text("{}")
    log = registry.parent / "consumption.jsonl"
    log.write_text("")
    write_json(
        tmp_path / "setup.json",
        {
            "run_id": "baseline-run",
            "cohort": {
                "id": "mixed_01",
                "registry": str(registry),
                "registry_sha256": "registry-digest",
                "comparison_arm": "baseline",
                "consumption_log": {
                    "path": log.name,
                    "schema_version": 1,
                },
            },
            "repository": {"runtime_source": {"sha256": "runtime"}},
        },
    )
    lock = {"selection_sha256": "selection", "predictions_sha256": "prediction"}

    first = record_tranche_consumption(tmp_path, lock)
    second = record_tranche_consumption(tmp_path, lock)

    assert first == second
    assert first["registry_sha256"] == "registry-digest"
    assert len(log.read_text().splitlines()) == 1
    with pytest.raises(SystemExit, match="consumed arm"):
        record_tranche_consumption(
            tmp_path, {**lock, "predictions_sha256": "different"}
        )


def test_registered_score_requires_exact_pinned_composite_gold(tmp_path: Path):
    gold = tmp_path / "answer_key"
    gold.mkdir()
    (gold / "KCNH2_recall_input.csv").write_text("pmid,variant\n1,A1V\n")
    (gold / "provenance.tsv").write_text(
        "gene\tpmid\tgold_provenance\nKCNH2\t1\ttest\n"
    )
    write_json(
        tmp_path / "setup.json",
        {
            "cohort": {
                "gold_root": str(gold),
                "gold_root_sha256": digest_gold_root(gold),
                "consumption_log": {"path": "consumption.jsonl"},
            }
        },
    )

    assert validated_score_gold_root(tmp_path, gold) == gold.resolve()
    wrong = tmp_path / "wrong"
    wrong.mkdir()
    with pytest.raises(SystemExit, match="differs from the setup-pinned"):
        validated_score_gold_root(tmp_path, wrong)

    (gold / "KCNH2_recall_input.csv").write_text("pmid,variant\n1,A2V\n")
    with pytest.raises(SystemExit, match="answer-key digest mismatch"):
        validated_score_gold_root(tmp_path, gold)


def test_registered_tranche_cannot_be_relocked_without_setup(tmp_path: Path):
    suite = tmp_path / "suite"
    suite.mkdir()
    manifest = suite / "tranche.tsv"
    manifest.write_text("KCNH2\t1\n")
    write_json(
        suite / "registry.json",
        {
            "consumption_log": {"path": "consumption.jsonl"},
            "tiers": [{"id": "mixed_01", "manifest": manifest.name}],
        },
    )
    write_json(
        tmp_path / "selection.json",
        {"paper_manifest": str(manifest), "papers": []},
    )
    write_json(tmp_path / "predictions.json", {"papers": []})

    with pytest.raises(SystemExit, match="burn contract"):
        command_lock(SimpleNamespace(run_dir=tmp_path))


def test_production_run_status_is_digest_bound_and_gold_free(tmp_path: Path):
    status_path = tmp_path / "RUN_STATUS.json"
    write_json(
        status_path,
        {
            "status": "completed",
            "exit_code": 0,
            "stage_failures": [],
            "gold_access": {
                "disabled": True,
                "gold_derived_alias_files_disabled": True,
            },
        },
    )
    predictions = {
        "strategy": "production_gvf_run",
        "primary_score_lane": "paper_derived",
        "papers": [{"gene": "KCNH2", "pmid": "1"}],
        "production_run_statuses": [
            {"gene": "KCNH2", "status": str(status_path), "sha256": digest(status_path)}
        ],
    }

    locked, errors = production_run_status_lock_entries(predictions)
    assert errors == []
    assert locked == predictions["production_run_statuses"]

    status_path.write_text("{}")
    _, errors = production_run_status_lock_entries(predictions)
    assert errors == ["production_run_statuses[0]: RUN_STATUS digest mismatch"]


def test_paired_report_comparison_applies_registered_pmid_cluster_rule():
    def report(arm: str, improved: bool) -> dict:
        return {
            "primary_score_lane": "paper_derived",
            "tranche_consumption": {
                "tier_id": "mixed_01",
                "comparison_arm": arm,
                "registry_sha256": "registry-digest",
            },
            "papers": [
                {
                    "gene": "KCNH2",
                    "pmid": str(pmid),
                    "tp": 2 if improved else 1,
                    "fp": 0,
                    "fn": 0 if improved else 1,
                }
                for pmid in (1, 2, 3)
            ],
        }

    design = {
        "decision_rule": {
            "primary": {
                "minimum_observed_delta": 0.01,
                "noninferiority_margin": -0.01,
                "one_sided_confidence_level": 0.95,
            },
            "precision_guardrail": {
                "noninferiority_margin": -0.02,
                "one_sided_confidence_level": 0.95,
            },
            "confidence_interval": {
                "method": "paired_cluster_bootstrap_nearest_rank",
                "cluster_unit": "PMID",
                "resamples": 1000,
                "seed": 7,
            },
        }
    }

    result = compare_paired_reports(
        report("baseline", False),
        report("candidate", True),
        design,
        tier_id="mixed_01",
        phase="discovery",
        registry_sha256="registry-digest",
    )

    assert result["cluster_count"] == 3
    assert result["metrics"]["recall"]["delta_candidate_minus_baseline"] == 0.5
    assert result["criteria"] == {
        "recall_pass": True,
        "precision_noninferiority_pass": True,
    }
    assert result["decision"] == "advance_to_confirmation"

    unchanged = compare_paired_reports(
        report("baseline", False),
        report("candidate", False),
        design,
        tier_id="mixed_01",
        phase="discovery",
        registry_sha256="registry-digest",
    )
    assert unchanged["passed"] is False
    assert unchanged["decision"] == "reject_or_revise_candidate"


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
                "telemetry_available": True,
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
            "telemetry_available": True,
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
        "integrity": {"llm_trace_manifest_sha256": "fixture-digest"},
    }


def test_deletion_range_does_not_match_single_residue_deletion():
    assert not matches("K1505del", "K1505_Q1507del", "SCN5A")
    assert not matches("K1505_Q1507del", "K1505del", "SCN5A")
    assert matches("K1505_Q1507del", "K1505_Q1507del", "SCN5A")


def test_material_digests_cover_production_record_and_every_representation(
    tmp_path: Path,
):
    source = tmp_path / "source.md"
    full_context = tmp_path / "full_context.md"
    extraction_record = tmp_path / "extraction.json"
    artifact = tmp_path / "artifact.json"
    pdf = tmp_path / "paper.pdf"
    figure = tmp_path / "figure.png"
    source.write_text("source")
    full_context.write_text("full source")
    extraction_record.write_text("{}")
    artifact.write_text("{}")
    pdf.write_bytes(b"%PDF fixture")
    figure.write_bytes(b"PNG fixture")
    paper = {
        "gene": "SCN5A",
        "pmid": "1",
        "source": str(source),
        "source_sha256": digest(source),
        "production_extraction_record": str(extraction_record),
        "production_extraction_record_sha256": digest(extraction_record),
        "representations": [str(source), str(full_context)],
        "representation_sha256": {
            str(source): digest(source),
            str(full_context): digest(full_context),
        },
        "artifacts": str(artifact),
        "artifacts_sha256": digest(artifact),
        "pdfs": [str(pdf)],
        "pdf_sha256": {str(pdf): digest(pdf)},
        "figures": [str(figure)],
        "figure_sha256": {str(figure): digest(figure)},
    }

    assert material_digest_errors(paper) == []

    full_context.write_text("changed")
    assert any(
        "text representation changed after selection" in error
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
        route_max_output_tokens=1600,
        max_output_tokens_ceiling=100000,
        legacy_table_material=False,
    )

    with pytest.raises(json.JSONDecodeError):
        command_extract(args)

    checkpoint = json.loads((tmp_path / "predictions.json").read_text())
    assert checkpoint["token_usage"]["total_tokens"] == 15
    assert checkpoint["papers"][0]["token_usage"]["total_tokens"] == 15
    refs = checkpoint["papers"][0]["llm_trace_refs"]
    assert [(ref["context"]["stage"], ref["success"]) for ref in refs] == [
        ("representation_route", True),
        ("representation_route_decision", None),
        ("paper_curation", True),
    ]
    trace_files = list((tmp_path / "llm_traces" / "SCN5A" / "1").glob("*.json"))
    assert len(trace_files) == 3


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


def test_markdown_does_not_invent_legacy_telemetry_or_traces(tmp_path: Path):
    report = _report_fixture()
    report["integrity"] = {
        "llm_trace_manifest_sha256": None,
        "llm_trace_report_sha256": None,
    }
    for paper in report["papers"]:
        paper["elapsed_seconds"] = 0.0
        paper["token_usage"] = {
            "telemetry_available": False,
            "total_tokens": None,
        }
    report["token_usage"] = {"telemetry_available": False}
    report["timing"] = {
        "wall_seconds": 0.0,
        "started_at": None,
        "completed_at": None,
    }

    markdown_path = tmp_path / "legacy_report.md"
    write_markdown_report(report, markdown_path)
    markdown = markdown_path.read_text()

    assert "token and timing telemetry was not captured" in markdown
    assert "zero placeholders must not be interpreted as zero cost" in markdown
    assert "Exact per-call LLM traces are not attached" in markdown
    assert "`llm_traces/<GENE>/<PMID>/`" not in markdown
    assert "| n/a | n/a |" in markdown


def test_markdown_discloses_production_projection_prelock_scoring(tmp_path: Path):
    report = _report_fixture()
    report["prelock_gold_usage"] = {
        "read_only_layer_scoring_possible": True,
        "scores_feed_back_into_predictions": False,
    }

    markdown_path = tmp_path / "projection_report.md"
    write_markdown_report(report, markdown_path)
    markdown = markdown_path.read_text()

    assert "may have read registered gold for read-only layer scorecards" in markdown
    assert "not the stricter native lock-before-any-gold-read protocol" in markdown


def test_markdown_discloses_explicit_gold_free_production_projection(tmp_path: Path):
    report = _report_fixture()
    report["prelock_gold_usage"] = {
        "read_only_layer_scoring_possible": False,
        "scores_feed_back_into_predictions": False,
        "provenance": "all_production_run_statuses_record_gold_access_disabled",
    }

    markdown_path = tmp_path / "gold_free_projection_report.md"
    write_markdown_report(report, markdown_path)
    markdown = markdown_path.read_text()

    assert "gold was used only for PMID eligibility" in markdown
    assert "may have read registered gold" not in markdown


def test_markdown_discloses_locked_production_trace_manifests(tmp_path: Path):
    report = _report_fixture()
    report["integrity"] = {
        "llm_trace_manifest_sha256": None,
        "production_trace_manifests": [{"gene": "SCN5A", "sha256": "abc"}],
    }

    markdown_path = tmp_path / "projection_trace_report.md"
    write_markdown_report(report, markdown_path)
    markdown = markdown_path.read_text()

    assert "production `gvf-run` trace manifest for every gene" in markdown
    assert "evaluation projection does not copy or relabel" in markdown
    assert "legacy runs require a rerun" not in markdown


def test_prediction_validation_accepts_explicit_deterministic_zero_tokens():
    paper = {
        "gene": "KCNH2",
        "pmid": "1",
        "tool": "table",
        "tool_rationale": "deterministic fixed-width parser",
        "elapsed_seconds": 0.0,
        "source_completeness": "full_text",
        "curation_rationale": "No model call was necessary.",
        "token_usage": {"telemetry_available": True, "total_tokens": 0},
        "llm_trace_refs": [
            {"context": {"stage": stage}}
            for stage in (
                "representation_route",
                "representation_route_decision",
                "paper_curation",
                "paper_curation_decision",
            )
        ],
        "variants": [],
    }
    selection = {"papers": [{"gene": "KCNH2", "pmid": "1"}]}

    assert (
        validate_predictions(selection, {"schema_version": 2, "papers": [paper]}) == []
    )


def _rendering(path: Path, rows: int, variants: int, padding: int = 0) -> Path:
    """Write a candidate rendering with a known number of table rows/variants."""
    lines = [f"| {chr(65 + i % 26)}{i + 1}V | {i} | carrier |" for i in range(rows)]
    lines += [f"prose mentioning p.Ala{i + 1}Val here" for i in range(variants)]
    lines.append("x" * padding)
    path.write_text("\n".join(lines))
    return path


def test_choose_source_prefers_a_strictly_richer_rendering(tmp_path: Path):
    full_context = _rendering(tmp_path / "a_FULL_CONTEXT.md", rows=0, variants=1)
    cleaned = _rendering(tmp_path / "a_CLEANED.md", rows=12, variants=9)

    assert choose_source([full_context, cleaned]) == cleaned


def test_choose_source_keeps_priority_when_neither_rendering_dominates(tmp_path: Path):
    """More table rows but less prose is a trade, not an improvement."""
    full_context = _rendering(
        tmp_path / "b_FULL_CONTEXT.md", rows=8, variants=9, padding=5000
    )
    cleaned = _rendering(tmp_path / "b_CLEANED.md", rows=12, variants=9)

    assert choose_source([full_context, cleaned]) == full_context


def test_usable_sources_selects_the_richer_rendering_on_disk(tmp_path: Path):
    paper_dir = tmp_path / "KCNQ1" / "17470695"
    paper_dir.mkdir(parents=True)
    _rendering(paper_dir / "17470695_FULL_CONTEXT.md", rows=0, variants=1, padding=3000)
    cleaned = _rendering(
        paper_dir / "17470695_CLEANED.md", rows=40, variants=30, padding=3000
    )

    papers = usable_sources(tmp_path, "KCNQ1", minimum_chars=100)

    assert [Path(p["source"]) for p in papers] == [cleaned.resolve()]
    audit = papers[0]["source_selection"]
    assert audit["policy"] == "pareto_richness"
    assert audit["selected"] == str(cleaned.resolve())
    assert len(audit["candidates"]) == 2
    assert "Pareto-dominated" in audit["rationale"]


def test_usable_sources_can_limit_expensive_hashing_to_eligible_pmids(tmp_path: Path):
    for pmid in ("17470695", "99999999"):
        paper_dir = tmp_path / "KCNQ1" / pmid
        paper_dir.mkdir(parents=True)
        _rendering(
            paper_dir / f"{pmid}_FULL_CONTEXT.md",
            rows=1,
            variants=1,
            padding=3000,
        )

    papers = usable_sources(
        tmp_path,
        "KCNQ1",
        minimum_chars=100,
        include_pmids={"17470695"},
    )

    assert [paper["pmid"] for paper in papers] == ["17470695"]


def _extract_captured_prompts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    source_text: str,
    artifact: str,
    route_tool: str = "table",
) -> tuple[dict, list]:
    source = tmp_path / "source.md"
    source.write_text(source_text)
    artifacts = tmp_path / "artifacts.json"
    artifacts.write_text(artifact)
    write_json(
        tmp_path / "selection.json",
        {
            "papers": [
                {
                    "gene": "KCNQ1",
                    "pmid": "17470695",
                    "source": str(source),
                    "source_sha256": digest(source),
                    "artifacts": str(artifacts),
                    "artifacts_sha256": digest(artifacts),
                    "pdfs": [],
                    "pdf_sha256": {},
                    "figures": [],
                    "figure_sha256": {},
                }
            ]
        },
    )
    write_json(
        tmp_path / "predictions.json",
        {
            "schema_version": 2,
            "token_usage": None,
            "papers": [
                {
                    "gene": "KCNQ1",
                    "pmid": "17470695",
                    "tool": None,
                    "token_usage": None,
                    "variants": [],
                }
            ],
        },
    )

    prompts: list = []
    responses = iter(
        [
            SimpleNamespace(
                output_text=json.dumps(
                    {
                        "tool": route_tool,
                        "tool_rationale": "Captions imply variant tables.",
                        "source_completeness": "partial_text",
                    }
                ),
                usage=SimpleNamespace(input_tokens=5, output_tokens=1),
            ),
            SimpleNamespace(
                output_text=json.dumps(
                    {
                        "notes": "No in-scope rows.",
                        "curation_rationale": (
                            "The supplied material had no qualifying human "
                            "variant/count row."
                        ),
                        "variants": [],
                    }
                ),
                usage=SimpleNamespace(input_tokens=7, output_tokens=2),
            ),
        ]
    )

    def create(**kwargs):
        prompts.append(kwargs["input"])
        return next(responses)

    fake_openai = ModuleType("openai")
    fake_openai.OpenAI = lambda **_kwargs: SimpleNamespace(
        responses=SimpleNamespace(create=create)
    )
    monkeypatch.setitem(sys.modules, "openai", fake_openai)
    monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.invalid")
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-key")
    command_extract(
        SimpleNamespace(
            run_dir=tmp_path,
            timeout=1,
            model="test-model",
            force=False,
            max_artifact_chars=1000,
            max_source_chars=10000,
            route_preview_chars=5000,
            max_ocr_images=1,
            route_reasoning_effort="medium",
            reasoning_effort="high",
            max_output_tokens=100,
            route_max_output_tokens=1600,
            max_output_tokens_ceiling=100000,
            legacy_table_material=False,
        )
    )
    return json.loads((tmp_path / "predictions.json").read_text()), prompts


def test_successful_eval_has_four_hash_lockable_trace_stages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    checkpoint, _ = _extract_captured_prompts(
        tmp_path,
        monkeypatch,
        source_text=(
            "| Variant | carriers |\n|---|---|\n| A1V | 2 |\nKCNQ1 patient evidence."
        ),
        artifact="{}",
        route_tool="table",
    )

    refs = checkpoint["papers"][0]["llm_trace_refs"]
    assert [ref["context"]["stage"] for ref in refs] == [
        "representation_route",
        "representation_route_decision",
        "paper_curation",
        "paper_curation_decision",
    ]
    setup_path = tmp_path / "setup.json"
    write_json(setup_path, {"cohort": {"comparison_arm": "baseline"}})
    command_lock(SimpleNamespace(run_dir=tmp_path))
    lock = json.loads((tmp_path / "LOCK.json").read_text())
    assert lock["setup_sha256"] == digest(setup_path)
    assert lock["llm_trace_manifest_sha256"]
    assert lock["llm_trace_report_sha256"]
    manifest = json.loads((tmp_path / "llm_traces" / "trace_manifest.json").read_text())
    assert manifest["llm_call_count"] == 2
    assert manifest["decision_event_count"] == 2
    report = tmp_path / "llm_trace_report.html"
    assert report.is_file()
    assert "Sent to model" in report.read_text(encoding="utf-8")
    assert report.stat().st_mode & 0o222 == 0
    assert setup_path.stat().st_mode & 0o222 == 0


def test_table_route_is_not_offered_without_real_table_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """An artifact full of captions must not advertise a table representation.

    Regression for KCNQ1 17470695: the artifact wrapper made the table string
    non-empty, the router picked ``table``, and the payload carried captions plus a
    keyword preview instead of any table rows.
    """
    checkpoint, prompts = _extract_captured_prompts(
        tmp_path,
        monkeypatch,
        source_text="KCNQ1 running text with no table rows at all.",
        artifact='{"table_captions_count": 5, "main_text": "TABLE 2. Mutations"}',
    )
    paper = checkpoint["papers"][0]

    assert paper["representations_available"] == ["text"]
    assert paper["tool"] == "text"
    assert "which was unavailable" in paper["tool_rationale"]
    assert "## TABLE PREVIEW\n\n" in prompts[0]
    assert "### Full/partial running text" in prompts[1]
    assert "no table rows at all" in prompts[1]


def test_parsed_artifact_reaches_the_model_on_a_non_table_route(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """The artifact used to ride along only on the table route."""
    _, prompts = _extract_captured_prompts(
        tmp_path,
        monkeypatch,
        source_text="KCNQ1 running text with no table rows at all.",
        artifact='{"supplements": ["Supplemental_Material.docx"]}',
    )

    assert "### Parsed artifact data" in prompts[1]
    assert "Supplemental_Material.docx" in prompts[1]
    assert "## PARSED ARTIFACT (attached to every route)" in prompts[0]


def test_grok_gets_no_reasoning_block_and_is_labelled_honestly():
    """Grok-4-class deployments 400 on reasoning.effort; they reason by default."""
    assert reasoning_params("gpt-5.6-sol", "high") == {"reasoning": {"effort": "high"}}
    assert effective_effort("gpt-5.6-sol", "high") == "high"

    assert reasoning_params("grok-4.3", "high") == {}
    assert effective_effort("grok-4.3", "high") == "model_default"


def test_ocr_route_is_withheld_from_a_model_without_vision():
    """grok-4.3 rejects image input; offering ocr would fail the paper outright."""
    assert supports_images("gpt-5.6-sol")
    assert not supports_images("grok-4.3")


def test_truncated_json_is_distinguished_from_garbage():
    """Azure reports status=completed even when output is cut off mid-string."""
    assert looks_truncated_json('{"variants": [{"variant": "A1V", "evidence": "row 1')
    assert looks_truncated_json('```json\n{"variants": [{"variant": "A1V"')
    assert not looks_truncated_json('{"variants": []}')
    assert not looks_truncated_json("not valid JSON")


def test_legacy_source_selection_ablation_restores_first_candidate(tmp_path: Path):
    """The ablation must reproduce the old behaviour exactly, or it proves nothing."""
    paper_dir = tmp_path / "KCNQ1" / "17470695"
    paper_dir.mkdir(parents=True)
    full = _rendering(
        paper_dir / "17470695_FULL_CONTEXT.md", rows=0, variants=1, padding=3000
    )
    _rendering(paper_dir / "17470695_CLEANED.md", rows=40, variants=30, padding=3000)

    legacy = usable_sources(
        tmp_path, "KCNQ1", minimum_chars=100, legacy_source_selection=True
    )
    fixed = usable_sources(tmp_path, "KCNQ1", minimum_chars=100)

    assert Path(legacy[0]["source"]) == full.resolve()
    assert Path(fixed[0]["source"]) != full.resolve()


def test_extraction_prompt_carries_table_attribution_guidance():
    """Compilation tables caused the 10973849 over-attribution regression."""
    # Pinned to the canonical copy so the harness cannot silently drift from it,
    # and so the guidance survives this harness being replaced.
    assert TABLE_ATTRIBUTION_GUIDANCE in EXTRACTION_INSTRUCTIONS
    assert "compilation citing other" in TABLE_ATTRIBUTION_GUIDANCE
    assert (
        "families, individuals, alleles, probands, cases" in TABLE_ATTRIBUTION_GUIDANCE
    )
    assert "Count only what this study observed" in TABLE_ATTRIBUTION_GUIDANCE
    assert "{" not in TABLE_ATTRIBUTION_GUIDANCE, "must be str.format-safe"
    # The measured dominant failure was omitting uncountable variants entirely.
    assert (
        "always\nemit the variant even when all three counts are null"
        in EXTRACTION_INSTRUCTIONS
    )


# ---------------------------------------------------------------------------
# BRCA2 arm: manifest support + gold fallback to the adjudicated overrides
# ---------------------------------------------------------------------------

GOLD_CSV_HEADER = "variant,pmid,carriers,affected,unaffected\n"


def _write_gold_csv(path: Path, pmid: str = "1") -> None:
    path.write_text(GOLD_CSV_HEADER + f"A1V,{pmid},2,1,1\n")


def _paper_dir(corpus: Path, gene: str, pmid: str) -> None:
    paper = corpus / gene / pmid
    paper.mkdir(parents=True)
    _rendering(paper / f"{pmid}_FULL_CONTEXT.md", rows=10, variants=5, padding=500)


def _prepare_args(tmp_path: Path, **overrides) -> SimpleNamespace:
    defaults = dict(
        seed=1,
        per_gene=1,
        minimum_chars=100,
        corpus_root=tmp_path / "corpus",
        gold_root=tmp_path / "gold",
        paper_manifest=None,
        runs_dir=tmp_path / "runs",
        run_id="run",
        legacy_source_selection=False,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_gold_csv_path_prefers_root_then_adjudicated_overrides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = tmp_path / "gold"
    overrides = tmp_path / "overrides"
    root.mkdir()
    overrides.mkdir()
    _write_gold_csv(root / "KCNH2_recall_input.csv")
    _write_gold_csv(overrides / "KCNH2_recall_input.csv")
    _write_gold_csv(overrides / "BRCA2_recall_input.csv")
    monkeypatch.setattr(run_eval_module, "GOLD_OVERRIDES", overrides)

    # An explicit root always wins for genes it carries; only absent genes fall
    # back, and a gene in neither location fails loudly instead of scoring empty.
    assert gold_csv_path(root, "KCNH2") == root / "KCNH2_recall_input.csv"
    assert gold_csv_path(root, "BRCA2") == overrides / "BRCA2_recall_input.csv"
    with pytest.raises(SystemExit, match="no gold CSV for APOE"):
        gold_csv_path(root, "APOE")


def test_paper_manifest_accepts_brca2_and_rejects_unregistered_genes(tmp_path: Path):
    good = tmp_path / "good.tsv"
    good.write_text("SCN5A\t123\nBRCA2\t26848529\n")
    assert read_paper_manifest(good) == [("SCN5A", "123"), ("BRCA2", "26848529")]

    bad = tmp_path / "bad.tsv"
    bad.write_text("BMPR2\t123\n")
    with pytest.raises(SystemExit):
        read_paper_manifest(bad)


@pytest.mark.parametrize(
    "manifest_name, expected_total, expected_brca2",
    [
        ("highcarrier48_plus_brca2_collaborator2_20260811.tsv", 50, 2),
        ("brca2_2_collaborator_reviewed_20260811.tsv", 2, 2),
        # Historical frozen manifests remain test-covered for reproducibility.
        ("highcarrier48_plus_brca2_20260810.tsv", 56, 8),
        ("brca2_8_papers_20260810.tsv", 8, 8),
    ],
)
def test_shipped_manifests_are_fully_gold_count_eligible(
    manifest_name: str, expected_total: int, expected_brca2: int
):
    """Every row of the shipped manifests must clear prepare's eligibility rule.

    Pins active and historical manifests against drift in either answer key:
    cardiac rows against the manual gold standard and BRCA2 rows against the
    fallback override. Active BRCA2 membership is separately restricted to the
    two lead-approved Variant Browser papers.
    """
    manifest = Path(run_eval_module.__file__).parent / manifest_name
    papers = read_paper_manifest(manifest)
    assert len(papers) == expected_total
    assert sum(1 for gene, _ in papers if gene == "BRCA2") == expected_brca2
    for gene in dict.fromkeys(gene for gene, _ in papers):
        eligible = gold_count_eligible_pmids(DEFAULT_GOLD, gene)
        missing = [pmid for g, pmid in papers if g == gene and pmid not in eligible]
        assert not missing, f"{gene}: not gold-count-eligible: {missing}"


def test_random_prepare_samples_cardiac_genes_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Seeded random runs must not require BRCA2 gold, source, or overrides."""
    gold = tmp_path / "gold"
    gold.mkdir()
    for gene in CARDIAC_GENES:
        _paper_dir(tmp_path / "corpus", gene, "1")
        _write_gold_csv(gold / f"{gene}_recall_input.csv")
    monkeypatch.setattr(run_eval_module, "GOLD_OVERRIDES", tmp_path / "absent")

    command_prepare(_prepare_args(tmp_path))

    selection = json.loads((tmp_path / "runs" / "run" / "selection.json").read_text())
    assert sorted(p["gene"] for p in selection["papers"]) == sorted(CARDIAC_GENES)
    assert set(selection["eligible_counts"]) == set(CARDIAC_GENES)
    assert set(selection["gold_sources"]) == set(CARDIAC_GENES)


def test_manifest_prepare_scores_brca2_against_override_gold(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    gold = tmp_path / "gold"
    gold.mkdir()
    overrides = tmp_path / "overrides"
    overrides.mkdir()
    _paper_dir(tmp_path / "corpus", "BRCA2", "99")
    _write_gold_csv(overrides / "BRCA2_recall_input.csv", pmid="99")
    monkeypatch.setattr(run_eval_module, "GOLD_OVERRIDES", overrides)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("BRCA2\t99\n")

    command_prepare(_prepare_args(tmp_path, paper_manifest=manifest))

    selection = json.loads((tmp_path / "runs" / "run" / "selection.json").read_text())
    assert [(p["gene"], p["pmid"]) for p in selection["papers"]] == [("BRCA2", "99")]
    assert selection["gold_sources"] == {
        "BRCA2": str(overrides / "BRCA2_recall_input.csv")
    }
    # A cardiac-only manifest run in the same layout must not touch BRCA2 pools.
    assert "SCN5A" not in selection["eligible_counts"]


def test_schema1_production_import_locks_without_call_telemetry():
    """Schema 1 is the external-import contract (production gvf-run projection).

    gvf-run does not aggregate per-paper wall time or exact token usage, so
    those checks must bind only to schema >= 2 (harness-native extraction) —
    otherwise the production baseline can never be locked and scored.
    """
    from benchmarks.codex_paper_eval.run_eval import validate_predictions

    paper = {
        "gene": "BRCA2",
        "pmid": "99",
        "tool": "text",
        "tool_rationale": "Production gvf-run strategy.",
        "source_completeness": "corpus_as_locked",
        "elapsed_seconds": None,
        "token_usage": {"telemetry_available": False, "total_tokens": None},
        "variants": [
            {
                "variant": "c.1T>A",
                "carriers": 1,
                "affected": 1,
                "unaffected": 0,
                "evidence": "Table 1 row",
                "source_location": "Table 1",
            }
        ],
    }
    selection = {"papers": [{"gene": "BRCA2", "pmid": "99"}]}

    assert (
        validate_predictions({**selection}, {"schema_version": 1, "papers": [paper]})
        == []
    )
    native_errors = validate_predictions(
        {**selection}, {"schema_version": 2, "papers": [dict(paper)]}
    )
    assert any("elapsed_seconds" in e for e in native_errors)
    assert any("token telemetry" in e for e in native_errors)


def test_production_projection_applies_field_level_trust_masks():
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location("db_to_predictions", converter_path)
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)

    mask = json.dumps(
        {
            "total_carriers": "quarantine",
            "affected": "trusted",
            "unaffected": "trusted",
        }
    )
    assert (
        converter.trusted_count(
            7,
            field="total_carriers",
            trust_tier="quarantine",
            field_trust=mask,
        )
        is None
    )
    assert (
        converter.trusted_count(
            5,
            field="affected",
            trust_tier="quarantine",
            field_trust=mask,
        )
        == 5
    )
    assert (
        converter.trusted_count(
            2,
            field="unaffected",
            trust_tier="quarantine",
            field_trust=None,
        )
        is None
    )


def test_production_projection_holds_ambiguous_identity_classes(tmp_path: Path):
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location(
        "db_to_predictions_identity", converter_path
    )
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)
    db = tmp_path / "BRCA2.db"
    con = sqlite3.connect(db)
    con.executescript(
        """
        CREATE TABLE variants (
          variant_id INTEGER PRIMARY KEY, protein_notation TEXT, cdna_notation TEXT
        );
        CREATE TABLE variant_papers (
          variant_id INTEGER, pmid TEXT, source_location TEXT, key_quotes TEXT,
          source_layer TEXT
        );
        CREATE TABLE penetrance_data (
          variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
          affected_count INTEGER, unaffected_count INTEGER, trust_tier TEXT,
          field_trust TEXT
        );
        CREATE TABLE vf_enrichment (
          variant_id INTEGER PRIMARY KEY, matched INTEGER, fp_class TEXT
        );
        INSERT INTO variants VALUES (1, 'p.Ala1Val', 'c.1C>T');
        INSERT INTO variants VALUES (2, 'p.Ala2Val', 'c.2C>T');
        INSERT INTO variants VALUES (3, 'p.Ala3Val', 'c.3C>T');
        INSERT INTO variants VALUES (4, 'p.Ala4Val', 'c.4C>T');
        INSERT INTO variants VALUES (5, 'p.Ser227del', 'c.828_830del');
        INSERT INTO variants VALUES (6, 'p.Ala6Val', NULL);
        INSERT INTO variants VALUES (7, 'p.Ala7Val', 'c.?');
        INSERT INTO variant_papers VALUES (1, '99', 'Table 1', '[]', 'llm_table');
        INSERT INTO variant_papers VALUES (2, '99', 'Table 1', '[]', 'llm_table');
        INSERT INTO variant_papers VALUES (3, '99', 'Table 1', '[]', 'llm_table');
        INSERT INTO variant_papers VALUES (4, '99', 'Table 1', '[]', 'llm_table');
        INSERT INTO variant_papers VALUES (5, '99', 'Results', '[]', 'llm_text');
        INSERT INTO variant_papers VALUES (6, '99', 'Table 1', '[]', 'llm_table');
        INSERT INTO variant_papers VALUES (7, '99', 'Table 1', '[]', 'llm_table');
        INSERT INTO penetrance_data VALUES (1, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO penetrance_data VALUES (2, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO penetrance_data VALUES (3, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO penetrance_data VALUES (4, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO penetrance_data VALUES (5, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO penetrance_data VALUES (6, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO penetrance_data VALUES (7, '99', 1, 1, 0, 'trusted', '{}');
        INSERT INTO vf_enrichment VALUES (1, 1, NULL);
        INSERT INTO vf_enrichment VALUES (2, 0, 'novel_in_range');
        INSERT INTO vf_enrichment VALUES (3, 0, 'residue_unverified');
        INSERT INTO vf_enrichment VALUES (4, 0, 'known_isoform_offset');
        INSERT INTO vf_enrichment VALUES (5, 0, 'residue_offset_suspect');
        INSERT INTO vf_enrichment VALUES (6, 0, 'residue_offset_suspect');
        INSERT INTO vf_enrichment VALUES (7, 0, 'residue_offset_suspect');
        """
    )
    con.commit()
    con.close()

    rows, dropped = converter.rows_for_gene(
        db,
        {"99"},
        set(),
        trust_mode="trusted",
        identity_mode="trusted",
    )

    assert [row["variant"] for row in rows["99"]] == [
        "p.Ala1Val c.1C>T",
        "p.Ala2Val c.2C>T",
        "p.Ala4Val c.4C>T",
        "p.Ser227del c.828_830del",
    ]
    assert dict(dropped) == {"99": 3}


def test_production_projection_keeps_valid_structural_only_identity(tmp_path: Path):
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location(
        "db_to_predictions_structural", converter_path
    )
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)
    db = tmp_path / "RYR2.db"
    con = sqlite3.connect(db)
    con.executescript(
        """
        CREATE TABLE variants (
          variant_id INTEGER PRIMARY KEY, protein_notation TEXT, cdna_notation TEXT,
          variant_class TEXT, structural_description TEXT
        );
        CREATE TABLE variant_papers (
          variant_id INTEGER, pmid TEXT, source_location TEXT, key_quotes TEXT,
          source_layer TEXT
        );
        CREATE TABLE penetrance_data (
          variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
          affected_count INTEGER, unaffected_count INTEGER, trust_tier TEXT,
          field_trust TEXT
        );
        CREATE TABLE vf_enrichment (
          variant_id INTEGER PRIMARY KEY, matched INTEGER, fp_class TEXT
        );
        INSERT INTO variants VALUES (
          1, NULL, NULL, 'exon_deletion',
          'deletion of exon 3 resulting in p.Asn57-Gly91del'
        );
        INSERT INTO variants VALUES (
          2, NULL, NULL, 'exon_deletion', 'deletion of exon 4'
        );
        INSERT INTO variants VALUES (
          3, NULL, NULL, 'exon_deletion', 'large event near the N terminus'
        );
        INSERT INTO variant_papers VALUES (
          1, '19216760', 'Results', '[]', 'llm_text'
        );
        INSERT INTO variant_papers VALUES (
          2, '19216760', 'Results', '[]', 'llm_text'
        );
        INSERT INTO variant_papers VALUES (
          3, '19216760', 'Results', '[]', 'llm_text'
        );
        INSERT INTO penetrance_data VALUES (
          1, '19216760', 4, NULL, NULL, 'trusted', '{}'
        );
        INSERT INTO penetrance_data VALUES (
          2, '19216760', 1, NULL, NULL, 'trusted', '{}'
        );
        INSERT INTO penetrance_data VALUES (
          3, '19216760', 1, NULL, NULL, 'trusted', '{}'
        );
        INSERT INTO vf_enrichment VALUES (1, 0, 'no_notation_suspect');
        INSERT INTO vf_enrichment VALUES (2, 0, 'wrong_gene_residue_mismatch');
        INSERT INTO vf_enrichment VALUES (3, 0, 'no_notation_suspect');
        """
    )
    con.commit()
    con.close()

    rows, dropped = converter.rows_for_gene(
        db,
        {"19216760"},
        set(),
        trust_mode="trusted",
        identity_mode="trusted",
    )

    assert rows["19216760"] == [
        {
            "variant": "deletion of exon 3 resulting in p.Asn57-Gly91del",
            "carriers": 4,
            "affected": None,
            "unaffected": None,
            "evidence": "no quote captured; source_layer=llm_text",
            "source_location": "Results",
            "source_layer": "llm_text",
        }
    ]
    assert dict(dropped) == {"19216760": 2}


def test_production_projection_does_not_merge_distinct_structural_breakpoints():
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location(
        "db_to_predictions_structural_twins", converter_path
    )
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)
    rows = [
        {
            "variant": (
                "c.168-301_c.273+722del1128 deletion of exon 3 resulting in "
                "p.Asn57-Gly91del"
            ),
            "carriers": 4,
            "affected": None,
            "unaffected": None,
            "source_layer": "llm_text",
        },
        {
            "variant": (
                "c.168-228_c.273+793del1126 deletion of exon 3 resulting in "
                "p.Asn57-Gly91del"
            ),
            "carriers": 2,
            "affected": None,
            "unaffected": None,
            "source_layer": "llm_text",
        },
    ]

    merged = converter.merge_same_variant(rows, "RYR2")
    assert len(merged) == 2
    assert {row["variant"] for row in merged} == {row["variant"] for row in rows}


def test_production_projection_keeps_exact_and_generic_structural_rows_distinct():
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location(
        "db_to_predictions_structural_generic", converter_path
    )
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)
    rows = [
        {
            "variant": "c.168-301_c.273+722del1128 deletion of exon 3",
            "carriers": 4,
            "affected": None,
            "unaffected": None,
            "source_layer": "llm_text",
        },
        {
            "variant": "deletion of exon 3",
            "carriers": 3,
            "affected": None,
            "unaffected": None,
            "source_layer": "llm_text",
        },
    ]

    merged = converter.merge_same_variant(rows, "RYR2")
    assert len(merged) == 2
    assert {row["variant"] for row in merged} == {row["variant"] for row in rows}


def test_production_projection_binds_and_revalidates_external_trace_manifest(
    tmp_path: Path,
):
    trace_root = configure_llm_tracing(
        tmp_path / "production" / "llm_traces", run_id="production-run"
    )
    with llm_trace_scope(gene="KCNH2", pmid="1", stage="fixture"):
        capture_llm_call(
            provider="fixture",
            requested_model="fixture-model",
            resolved_model="fixture-model",
            request={"input": "fixture"},
            call=lambda: SimpleNamespace(
                output_text="fixture",
                usage=SimpleNamespace(input_tokens=1, output_tokens=1, total_tokens=2),
            ),
        )
    manifest_path = trace_root / "trace_manifest.json"
    manifest = build_trace_manifest(
        trace_root, output_path=manifest_path, run_id="production-run"
    )
    predictions = {
        "strategy": "production_gvf_run",
        "papers": [{"gene": "KCNH2", "pmid": "1"}],
        "production_trace_manifests": [
            {
                "gene": "KCNH2",
                "manifest": str(manifest_path),
                "sha256": digest(manifest_path),
                "run_id": manifest["run_id"],
                "llm_call_count": manifest["llm_call_count"],
                "decision_event_count": manifest["decision_event_count"],
                "integrity_level": manifest["verification"]["level"],
            }
        ],
    }

    locked, roots, errors = production_trace_lock_entries(predictions)
    assert errors == []
    assert locked[0]["sha256"] == digest(manifest_path)
    assert roots == [trace_root]

    manifest_path.write_text(manifest_path.read_text() + "\n")
    _locked, _roots, errors = production_trace_lock_entries(predictions)
    assert any("digest mismatch" in error for error in errors)


def test_new_production_projection_cannot_omit_trace_manifests():
    _locked, _roots, errors = production_trace_lock_entries(
        {
            "strategy": "production_gvf_run",
            "papers": [{"gene": "KCNH2", "pmid": "1"}],
        }
    )

    assert errors == [
        "production_gvf_run predictions require production_trace_manifests"
    ]


def test_markdown_report_renders_only_genes_present_in_run():
    report = _report_fixture()
    report["by_gene"]["BRCA2"] = report["by_gene"]["SCN5A"]

    lines = []

    class _Sink:
        def write_text(self, text):
            lines.append(text)

    write_markdown_report(report, _Sink())
    assert "| BRCA2 |" in lines[0]

    del report["by_gene"]["BRCA2"]
    write_markdown_report(report, _Sink())
    assert "| BRCA2 |" not in lines[1]


def test_table_route_carries_full_text_not_a_keyword_preview(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Provenance needs the caption/footnotes; a 6k excerpt strips them.

    Regression for KCNH2 10973849, where a compilation table read as first-party data.
    """
    marker = "SENTINEL_PROSE_FAR_FROM_ANY_VARIANT_MENTION"
    source = (
        "| Mutation | n | Type |\n| A1V | 2 | Missense |\n| G2R | 3 | Missense |\n"
        # Long enough to fall outside targeted_preview's reach, short enough to
        # survive max_source_chars, so the assertion tests the route not the cap.
        + ("filler prose. " * 400)
        + marker
    )
    _, prompts = _extract_captured_prompts(
        tmp_path, monkeypatch, source_text=source, artifact="{}", route_tool="table"
    )

    assert "### Structured table rows" in prompts[1]
    assert marker in prompts[1], "full running text must reach the table route"


def test_matches_legacy_stop_star_and_frameshift_notations():
    # HGVS '*' stop against legacy X, inside compound protein+cDNA strings.
    assert matches("p.R518* c.1522C>T", "R518X", "KCNQ1")
    assert matches("p.Q530* c.1588C>T", "Q530X", "KCNQ1")
    # Legacy inline-stop-distance frameshifts and 3-letter fs without new AA.
    assert matches("p.Arg192fs c.573_577del", "R192CFS91X", "KCNQ1")
    assert matches("P400fs", "P400FS/62X", "KCNQ1")
    # Distinct missense at the same codon must never match.
    assert not matches("D609N", "D609G", "KCNH2")
    assert not matches("p.Asp609Asn c.1825G>A", "D609G", "KCNH2")


def test_matches_reference_less_protein_range_duplication():
    """Legacy tables can omit residues while preserving an exact range/op."""
    assert matches(
        "p.360_361dupKQ c.1071_1076dupGAAGCA",
        "360_361DUPKQ",
        "KCNQ1",
    )
    assert not matches(
        "c.1071_1076dupGAAGCA",
        "360_361DUPKQ",
        "KCNQ1",
    )
    assert not any(
        candidate.lower().startswith("c.")
        for candidate in run_eval_module.variant_candidates("p.360_361dupKQ", "KCNQ1")
    )


def test_matches_normalizes_legacy_arrow_spellings():
    assert matches("c.2592+1G->A", "c.2592+1G>A", "KCNH2")
    assert matches("c.2592+1G→A", "c.2592+1G>A", "KCNH2")


def test_matches_splice_shorthand_tolerates_separators():
    assert matches("A344/SP", "A344SP", "KCNQ1")
    assert matches("A344A/SPLICE", "A344SP", "KCNQ1")


def test_splice_bridge_requires_intronic_offset():
    # Curators encode splice-site variants as a terminal event at the
    # flanking codon: c.477+1G>A is M159X.
    assert matches("c.477+1G>A", "M159X", "KCNQ1")
    # An exonic substitution at the same codon must NOT codon-bridge.
    assert not matches("c.477G>A", "M159X", "KCNQ1")
    assert not matches("c.940G>A", "G314X", "KCNQ1")


def test_structural_alleles_are_gene_map_not_hardcoded():
    assert matches("EXON 3 DELETION", "p.Asn57_Gly91del", "RYR2")
    assert matches("p.N57_G91del", "EXON 3 DELETION", "RYR2")
    assert not matches("EXON 3 DELETION", "p.Asn57_Gly91del", "KCNH2")
    assert matches("ΔKPQ", "p.K1505_Q1507del", "SCN5A")
    assert not matches("ΔKPQ", "p.K1505_Q1507del", "KCNH2")


def test_translation_bridge_is_nucleotide_verified():
    # c.1127G>A turns an Arg codon (CGC) into His (CAC): same allele.
    assert matches("c.1127G>A", "R376H", "SCN5A")
    assert matches("c.1826A>G", "D609G", "KCNH2")
    # The same cDNA change can never be the other missense at that codon.
    assert not matches("c.1826A>G", "D609N", "KCNH2")
    # Codon-index agreement alone is not enough.
    assert not matches("c.1130G>A", "R376H", "SCN5A")


def test_merge_notation_twins_bounds():
    from benchmarks.codex_paper_eval.run_eval import merge_notation_twins

    rows = [
        {
            "variant": "p.Arg376His",
            "carriers": 12,
            "affected": None,
            "unaffected": None,
        },
        {"variant": "c.1127G>A", "carriers": 12, "affected": 3, "unaffected": None},
        {"variant": "D609G", "carriers": 2, "affected": 2, "unaffected": None},
        {"variant": "D609N", "carriers": 1, "affected": 0, "unaffected": None},
    ]
    merged, twins = merge_notation_twins(rows, "SCN5A")
    assert twins == 1 and len(merged) == 3
    # The kept row already carries counts, so the complementary twin fields
    # are REFUSED (a union of two partially-counted rows can invent a count
    # vector no single source row asserted).
    assert merged[0]["carriers"] == 12 and merged[0]["affected"] is None

    # Conflicting non-null counts refuse the merge (different patients).
    conflicted = [
        {
            "variant": "p.Arg376His",
            "carriers": 12,
            "affected": None,
            "unaffected": None,
        },
        {"variant": "c.1127G>A", "carriers": 4, "affected": None, "unaffected": None},
    ]
    merged, twins = merge_notation_twins(conflicted, "SCN5A")
    assert twins == 0 and len(merged) == 2


def test_score_one_merges_twins_before_scoring():
    prediction = {
        "variants": [
            {
                "variant": "R376H",
                "carriers": None,
                "affected": None,
                "unaffected": None,
            },
            {
                "variant": "c.1127G>A",
                "carriers": 12,
                "affected": None,
                "unaffected": None,
            },
        ]
    }
    gold = [{"variant": "R376H", "carriers": 12, "affected": None, "unaffected": None}]
    score = score_one("SCN5A", "1", prediction, gold)
    # One TP, no stranded cDNA FP, and the twin's count reaches the match.
    assert (score["tp"], score["fp"], score["fn"]) == (1, 0, 0)
    assert score["merged_notation_twin_rows"] == 1
    assert score["count"]["carriers"]["predicted"] == 1


def test_score_one_reports_zero_stratified_count_recall():
    prediction = {
        "variants": [
            {"variant": "A1V", "carriers": 3, "affected": None, "unaffected": 0},
            {"variant": "A2V", "carriers": None, "affected": None, "unaffected": None},
        ]
    }
    gold = [
        {"variant": "A1V", "carriers": 3, "affected": 3, "unaffected": 0},
        {"variant": "A2V", "carriers": 2, "affected": 1, "unaffected": 0},
    ]
    score = score_one("KCNH2", "1", prediction, gold)
    carriers = score["count"]["carriers"]
    assert carriers["gold_asserted_nonzero"] == 2
    assert carriers["gold_asserted_zero"] == 0
    assert carriers["recall_nonzero_gold"] == pytest.approx(1 / 2)
    unaffected = score["count"]["unaffected"]
    assert unaffected["gold_asserted_zero"] == 2
    assert unaffected["recall_zero_gold"] == pytest.approx(1 / 2)
    assert unaffected["recall_nonzero_gold"] is None
    combined = aggregate([score])
    assert combined["count"]["unaffected"]["gold_asserted_zero"] == 2
    assert combined["count"]["carriers"]["recall_nonzero_gold"] == pytest.approx(1 / 2)


def test_linkage_codon_shadows_are_excluded_from_the_projection(tmp_path: Path):
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location("db_to_predictions", converter_path)
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)

    source = tmp_path / "PMID_FULL_CONTEXT.md"
    source.write_text(
        "The proband carried R376H. A relative carried p.Ala561Thr as well."
    )
    rows = [
        # Independently-extracted anchors at codons 376 and 561.
        {
            "variant": "R376H",
            "source_layer": "llm_text",
            "carriers": 2,
            "affected": 1,
            "unaffected": None,
        },
        {
            "variant": "A561V",
            "source_layer": "llm_table",
            "carriers": None,
            "affected": None,
            "unaffected": None,
        },
        # Ungrounded ClinVar neighbor at the SAME codon: enumeration artifact.
        {
            "variant": "p.Arg376Cys",
            "source_layer": "clinvar",
            "carriers": None,
            "affected": None,
            "unaffected": None,
        },
        # ClinVar row at an anchored codon that IS in the text stays.
        {
            "variant": "p.Ala561Thr",
            "source_layer": "clinvar",
            "carriers": None,
            "affected": None,
            "unaffected": None,
        },
        # Ungrounded ClinVar row at an UNANCHORED codon stays: on papers whose
        # tables never reached disk, linkage is the only recall signal.
        {
            "variant": "p.Ser277del",
            "source_layer": "clinvar",
            "carriers": None,
            "affected": None,
            "unaffected": None,
        },
    ]
    kept, dropped = converter.drop_linkage_shadows(rows, "KCNQ1", str(source))
    kept_variants = [r["variant"] for r in kept]
    assert dropped == 1
    assert "p.Arg376Cys" not in kept_variants
    assert {"R376H", "A561V", "p.Ala561Thr", "p.Ser277del"} <= set(kept_variants)

    # Missing or unreadable source keeps everything.
    kept_all, dropped_none = converter.drop_linkage_shadows(
        rows, "KCNQ1", str(tmp_path / "missing.md")
    )
    assert dropped_none == 0 and len(kept_all) == len(rows)


def test_production_projection_classifies_origin_not_later_witnesses():
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location(
        "db_to_predictions_lanes", converter_path
    )
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)

    assert converter.provenance_lane("llm_table,clinvar") == "paper_derived"
    assert converter.provenance_lane("mixed") == "unclassified"
    assert converter.provenance_lane("manual_or_legacy") == "unclassified"
    assert converter.provenance_lane("clinvar,pubtator") == "external_linkage"
    assert converter.provenance_lane("pubtator") == "external_linkage"
    assert converter.provenance_lane(None) == "unclassified"


def test_production_projection_filters_strict_provenance_lanes(tmp_path: Path):
    converter_path = (
        Path(__file__).parents[2]
        / "benchmarks/codex_paper_eval"
        / "db_to_predictions.py"
    )
    spec = importlib.util.spec_from_file_location(
        "db_to_predictions_lane_filter", converter_path
    )
    assert spec and spec.loader
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)

    db = tmp_path / "KCNH2.db"
    con = sqlite3.connect(db)
    con.executescript(
        """
        CREATE TABLE variants (
          variant_id INTEGER PRIMARY KEY, protein_notation TEXT, cdna_notation TEXT
        );
        CREATE TABLE variant_papers (
          variant_id INTEGER, pmid TEXT, source_location TEXT, key_quotes TEXT,
          source_layer TEXT
        );
        CREATE TABLE penetrance_data (
          variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
          affected_count INTEGER, unaffected_count INTEGER
        );
        INSERT INTO variants VALUES (1, 'p.Ala1Val', NULL);
        INSERT INTO variants VALUES (2, 'p.Ala2Val', NULL);
        INSERT INTO variants VALUES (3, 'p.Ala3Val', NULL);
        INSERT INTO variant_papers VALUES (1, '99', 'Table 1', 'paper', 'llm_table,clinvar');
        INSERT INTO variant_papers VALUES (2, '99', 'ClinVar', 'link', 'clinvar');
        INSERT INTO variant_papers VALUES (3, '99', 'legacy', 'unknown', 'mystery');
        """
    )
    con.commit()
    con.close()

    kwargs = {"trust_mode": "all", "identity_mode": "all"}
    paper, _ = converter.rows_for_gene(
        db, {"99"}, set(), provenance="paper_derived", **kwargs
    )
    linkage, _ = converter.rows_for_gene(
        db, {"99"}, set(), provenance="external_linkage", **kwargs
    )
    scored_union, _ = converter.rows_for_gene(
        db, {"99"}, set(), provenance="scored_union", **kwargs
    )

    assert [row["variant"] for row in paper["99"]] == ["p.Ala1Val"]
    assert [row["variant"] for row in linkage["99"]] == ["p.Ala2Val"]
    assert {row["variant"] for row in scored_union["99"]} == {
        "p.Ala1Val",
        "p.Ala2Val",
    }


def test_twin_merge_never_fuses_distinct_splice_alleles():
    from benchmarks.codex_paper_eval.run_eval import merge_notation_twins

    # matches() may splice-bridge both intronic alleles to M159X for SCORING,
    # but identity must not: protein-first order used to fuse all three.
    rows = [
        {"variant": "M159X", "carriers": None, "affected": None, "unaffected": None},
        {"variant": "c.477+1G>A", "carriers": 2, "affected": 1, "unaffected": None},
        {"variant": "c.477+2T>C", "carriers": 1, "affected": 1, "unaffected": None},
    ]
    merged, twins = merge_notation_twins(rows, "KCNQ1")
    assert len(merged) == 3 and twins == 0
    # Scoring still bridges: matches() keeps the recall win.
    assert matches("c.477+1G>A", "M159X", "KCNQ1")


def test_deletion_span_bridge_matches_endpoint_spellings():
    """``c.693delCA`` and ``c.692_693delCA`` are the same two-base event.

    Curated gold only ever uses the span form (0 of 6971 rows use the
    single-coordinate spelling), so without this bridge the paper's spelling
    is scored as a false positive next to the gold row it is identical to.
    """
    assert run_eval_module.matches("c.693delCA", "c.692_693delCA", "SCN5A")
    assert run_eval_module.matches("c.692delCA", "c.692_693delCA", "SCN5A")


def test_deletion_span_bridge_refuses_unrelated_events():
    match = run_eval_module.cdna_deletion_endpoint_match
    # Different deleted bases are different alleles.
    assert not match("c.693delCA", "c.692_693delCT")
    # A coordinate that is not an endpoint of the span never bridges.
    assert not match("c.694delCA", "c.692_693delCA")
    # Two open spellings, or two closed spellings, are compared as text.
    assert not match("c.693delCA", "c.694delCA")
    assert not match("c.692_693delCA", "c.693_694delCA")
    # Single-base deletions are already unambiguous and must not bridge to a
    # range.
    assert not match("c.123delA", "c.122_123delA")


def test_deletion_span_parse_keeps_the_open_end_open():
    parse = run_eval_module.cdna_deletion_span
    assert parse("c.692_693delCA") == (692, 693, "CA")
    assert parse("c.693delCA") == (693, None, "CA")
    assert parse("c.1234A>G") is None


def test_deletion_span_bridge_is_allele_identity_for_twin_merge():
    """The bridge is identity, not scoring fuzz, so twin-merge may use it.

    ``merge_notation_twins`` refuses a row that is identical to more than one
    kept row, which is what keeps the residual ``c.693_694delCA`` ambiguity
    from fusing distinct alleles.
    """
    assert run_eval_module.twin_identical("c.693delCA", "c.692_693delCA", "SCN5A")
    rows = [
        {"variant": "c.693delCA", "carriers": 2, "affected": None, "unaffected": None},
        {
            "variant": "c.692_693delCA",
            "carriers": None,
            "affected": None,
            "unaffected": None,
        },
    ]
    merged, twins = run_eval_module.merge_notation_twins(rows, "SCN5A")
    assert twins == 1
    assert len(merged) == 1
    assert merged[0]["carriers"] == 2


def test_twin_merge_refuses_an_ambiguous_deletion_endpoint():
    rows = [
        {
            "variant": "c.692_693delCA",
            "carriers": 1,
            "affected": None,
            "unaffected": None,
        },
        {
            "variant": "c.693_694delCA",
            "carriers": 1,
            "affected": None,
            "unaffected": None,
        },
        {"variant": "c.693delCA", "carriers": 2, "affected": None, "unaffected": None},
    ]
    merged, twins = run_eval_module.merge_notation_twins(rows, "SCN5A")
    assert twins == 0
    assert len(merged) == 3
