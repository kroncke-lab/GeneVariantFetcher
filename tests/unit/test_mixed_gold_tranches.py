"""Integrity checks for the complete mixed-gold protocol-regression suite."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from benchmarks.codex_paper_eval.db_to_predictions import trace_usage
from benchmarks.evaluation_tiers.build_mixed_tranches import (
    digest_answer_key,
    excluded_pmids,
    sha256_file,
)


REPO = Path(__file__).parents[2]
SUITE = REPO / "benchmarks" / "evaluation_tiers" / "mixed_gold"


def _manifest_rows(path: Path) -> list[tuple[str, str]]:
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            gene, pmid = line.split()
            rows.append((gene, pmid))
    return rows


def test_exclusion_manifests_remove_whole_articles_and_bind_their_bytes(
    tmp_path: Path,
):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    first.write_text("# consumed\nKCNH2\t123\nKCNQ1\t123\n")
    second.write_text("SCN5A\t456\n")

    pmids, metadata = excluded_pmids([first, second])

    assert pmids == {"123", "456"}
    assert [record["sha256"] for record in metadata] == [
        sha256_file(first),
        sha256_file(second),
    ]


def test_every_runnable_attempt_is_assigned_once_and_articles_are_atomic():
    registry = json.loads((SUITE / "registry.json").read_text())
    with (SUITE / "inventory.tsv").open(newline="") as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    included = {
        (row["gene"], row["pmid"]) for row in inventory if row["status"] == "included"
    }

    assigned: list[tuple[str, str]] = []
    pmid_tranches: dict[str, set[str]] = defaultdict(set)
    with (SUITE / "answer_key" / "provenance.tsv").open(newline="") as handle:
        provenance = {
            (row["gene"], row["pmid"]): row["gold_provenance"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    for tier in registry["tiers"]:
        manifest = SUITE / tier["manifest"]
        rows = _manifest_rows(manifest)
        assert sha256_file(manifest) == tier["sha256"]
        assert len(rows) == tier["attempt_count"]
        assert len({pmid for _, pmid in rows}) == tier["unique_pmid_count"]
        assert Counter(gene for gene, _ in rows) == Counter(tier["gene_attempt_counts"])
        assert Counter(provenance[row] for row in rows) == Counter(
            tier["gold_provenance_attempt_counts"]
        )
        assert 1 <= len(rows) <= registry["target_tranche_size"]
        for gene, pmid in rows:
            assigned.append((gene, pmid))
            pmid_tranches[pmid].add(tier["id"])

    assert len(assigned) == len(set(assigned))
    assert set(assigned) == included
    assert all(len(tranches) == 1 for tranches in pmid_tranches.values())
    assert len(assigned) == registry["inventory"]["source_available_attempts"]


def test_inventory_is_complete_and_costs_reconcile():
    registry = json.loads((SUITE / "registry.json").read_text())
    with (SUITE / "inventory.tsv").open(newline="") as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    statuses = Counter(row["status"] for row in inventory)

    assert sha256_file(SUITE / "inventory.tsv") == registry["inventory"]["sha256"]
    assert len(inventory) == registry["inventory"]["gold_attempts"] == 1534
    assert statuses == {
        "included": 1422,
        "source_unavailable": 111,
        "quarantined": 1,
    }
    assert registry["primary_score_lane"] == "paper_derived"
    assert registry["comparison_score_lanes"] == ["linkage_assisted"]
    assert registry["evaluation_design"]["primary_endpoint"] == (
        "paper_derived_micro_variant_identity_recall"
    )
    assert registry["evaluation_design"]["cluster_unit"] == "PMID"
    rule = registry["evaluation_design"]["decision_rule"]
    assert rule["delta_definition"] == "candidate_minus_baseline_on_the_same_tranche"
    assert rule["primary"]["minimum_observed_delta"] == 0.01
    assert rule["primary"]["noninferiority_margin"] == -0.01
    assert rule["precision_guardrail"]["noninferiority_margin"] == -0.02
    assert rule["confidence_interval"] == {
        "method": "paired_cluster_bootstrap_nearest_rank",
        "cluster_unit": "PMID",
        "resamples": 10000,
        "seed": registry["selection_seed"],
    }

    estimated = sum(tier["estimated_cost_usd"] for tier in registry["tiers"])
    budget = sum(tier["budget_with_headroom_usd"] for tier in registry["tiers"])
    assert estimated == pytest.approx(
        registry["cost_model"]["estimated_suite_cost_usd"], abs=1e-4
    )
    assert budget == pytest.approx(
        registry["cost_model"]["budget_with_headroom_usd"], abs=1e-3
    )
    assert registry["cost_model"]["paired_estimated_suite_cost_usd"] == (
        pytest.approx(estimated * 2, abs=1e-4)
    )
    assert registry["cost_model"]["paired_budget_with_headroom_usd"] == (
        pytest.approx(registry["cost_model"]["budget_with_headroom_usd"] * 2)
    )
    assert all(
        tier["primary_score_lane"] == "paper_derived"
        and tier["comparison_score_lanes"] == ["linkage_assisted"]
        and tier["eligibility_mode"] == "variant"
        and tier["paired_estimated_cost_usd"]
        == pytest.approx(tier["estimated_cost_usd"] * 2, abs=2e-4)
        for tier in registry["tiers"]
    )
    assert digest_answer_key(SUITE / "answer_key") == registry["answer_key"]["sha256"]
    for source in registry["gold_inputs"]:
        path = REPO / source["path"]
        assert path.is_file()
        assert sha256_file(path) == source["sha256"]
    calibration = REPO / registry["cost_model"]["calibration"]
    assert sha256_file(calibration) == registry["cost_model"]["calibration_sha256"]
    cost_profile = json.loads(calibration.read_text())
    for observed in cost_profile["calibrations"].values():
        source = REPO / observed["source"]
        assert source.is_file()
        assert sha256_file(source) == observed["source_sha256"]


def test_cost_calibration_matches_exact_locked_trace_usage():
    profile = json.loads(
        (REPO / "benchmarks" / "evaluation_tiers" / "cost_calibration.json").read_text()
    )
    for gene, calibration in profile["calibrations"].items():
        source = REPO / calibration["source"]
        if source.name == "predictions.json":
            predictions = json.loads(source.read_text())
            actual: dict[str, Counter] = defaultdict(Counter)
            papers = [paper for paper in predictions["papers"] if paper["gene"] == gene]
            assert len(papers) == calibration["attempts"]
            for paper in papers:
                for model, usage in paper["token_usage"]["models"].items():
                    actual[model].update(
                        {
                            "calls": usage["llm_calls"],
                            "input_tokens": usage["input_tokens"],
                            "output_tokens": usage["output_tokens"],
                        }
                    )
        else:
            usage, _ = trace_usage(source.parent)
            actual = {
                model: Counter(
                    {
                        "calls": values["llm_calls"],
                        "input_tokens": values["input_tokens"],
                        "output_tokens": values["output_tokens"],
                    }
                )
                for model, values in usage["models"].items()
            }
        assert {model: dict(values) for model, values in actual.items()} == calibration[
            "models"
        ]
