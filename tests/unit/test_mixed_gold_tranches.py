"""Integrity checks for the complete mixed-gold protocol-regression suite."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from utils.gold_standard import gold_digest_lineage

from benchmarks.evaluation_tiers.build_mixed_tranches import (
    digest_answer_key,
    excluded_pmids,
    sha256_file,
)


REPO = Path(__file__).parents[2]
SUITE = REPO / "benchmarks" / "evaluation_tiers" / "mixed_gold"


# Usage-only receipts are exported from the original ignored operator traces.
# Keep their bytes pinned without rewriting historical cost/registry digests.
USAGE_RECEIPTS = REPO / "tests/fixtures/cost_calibration_usage.json"
USAGE_RECEIPTS_SHA256 = (
    "3617fc19464fa32734c8b933a3a5f0c058faa6790962f01e53e5062f6c588b97"
)


def _usage_receipts():
    assert sha256_file(USAGE_RECEIPTS) == USAGE_RECEIPTS_SHA256
    receipts = json.loads(USAGE_RECEIPTS.read_text())
    assert receipts["calibration_sha256"] == sha256_file(
        REPO / "benchmarks/evaluation_tiers/cost_calibration.json"
    )
    return receipts["calibrations"]


def _manifest_rows(path: Path) -> list[tuple[str, str]]:
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            gene, pmid = line.split()
            rows.append((gene, pmid))
    return rows


def test_usage_export_ignores_unrelated_json_but_rejects_malformed_calls(
    tmp_path, monkeypatch
):
    from benchmarks.evaluation_tiers import export_cost_usage_receipts as exporter

    monkeypatch.setattr(exporter, "ROOT", tmp_path)
    traces = tmp_path / "traces"
    traces.mkdir()
    manifest = traces / "manifest.json"
    manifest.write_text('{"schema_version": 1}')
    (traces / "unrelated.json").write_text("[1, 2, 3]")
    (traces / "call.json").write_text(
        json.dumps(
            {
                "record_type": "llm_call",
                "context": {"model": "offline"},
                "response": {
                    "text": "private response",
                    "usage": {
                        "prompt_tokens": 100,
                        "total_tokens": 140,
                    },
                },
            }
        )
    )
    profile = tmp_path / "profile.json"
    profile.write_text(
        json.dumps(
            {
                "calibrations": {
                    "GENE": {
                        "source": "traces/manifest.json",
                        "source_sha256": sha256_file(manifest),
                        "attempts": 1,
                        "models": {
                            "offline": {
                                "calls": 1,
                                "input_tokens": 100,
                                "output_tokens": 40,
                            }
                        },
                    }
                }
            }
        )
    )
    original = profile.read_bytes()
    destination = tmp_path / "receipt.json"
    receipt = exporter.export(profile, destination)
    calls = receipt["calibrations"]["GENE"]["calls"]
    assert len(calls) == 1
    assert calls[0]["input_tokens"] == 100 and calls[0]["output_tokens"] == 40
    assert "private response" not in destination.read_text()
    assert profile.read_bytes() == original

    # Ignore unrelated records, but never silently omit a malformed model call
    # from the cost ledger or write a new receipt from incomplete usage.
    valid_call = json.loads((traces / "call.json").read_text())
    destination.unlink()
    for field, malformed in (
        ("response", "error text"),
        ("context", []),
        ("response", {"usage": [100, 140]}),
    ):
        (traces / "call.json").write_text(json.dumps({**valid_call, field: malformed}))
        with pytest.raises(ValueError, match="Malformed llm_call"):
            exporter.export(profile, destination)
        assert not destination.exists()
        assert profile.read_bytes() == original


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
        # The registry pins the gold it was built from. An approved adjudication
        # may move the live file forward, but only through the append-only
        # revision log; an unrecorded edit still fails here.
        assert source["sha256"] in gold_digest_lineage(path), (
            f"{source['path']}: pinned digest is neither the live file nor a "
            "recorded revision in gene_variant_fetcher_gold_standard/gold_revisions.jsonl"
        )
    calibration = REPO / registry["cost_model"]["calibration"]
    assert sha256_file(calibration) == registry["cost_model"]["calibration_sha256"]
    cost_profile = json.loads(calibration.read_text())
    receipts = _usage_receipts()
    for gene, observed in cost_profile["calibrations"].items():
        source = REPO / observed["source"]
        if source.name == "predictions.json":
            assert sha256_file(source) == observed["source_sha256"]
        else:
            # Full traces are deliberately gitignored. A fresh checkout validates
            # the exported source identity and usage receipts, not local files.
            receipt = receipts[gene]
            assert receipt["source"] == observed["source"]
            assert receipt["source_sha256"] == observed["source_sha256"]
            assert receipt["attempts"] == observed["attempts"]


def test_cost_calibration_matches_locked_predictions_and_trace_usage_receipts():
    profile = json.loads(
        (REPO / "benchmarks" / "evaluation_tiers" / "cost_calibration.json").read_text()
    )
    receipts = _usage_receipts()
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
            actual = defaultdict(Counter)
            for call in receipts[gene]["calls"]:
                assert len(call["trace_sha256"]) == 64
                actual[call["model"]].update(
                    {
                        "calls": 1,
                        "input_tokens": call["input_tokens"],
                        "output_tokens": call["output_tokens"],
                    }
                )
        assert {model: dict(values) for model, values in actual.items()} == calibration[
            "models"
        ]
