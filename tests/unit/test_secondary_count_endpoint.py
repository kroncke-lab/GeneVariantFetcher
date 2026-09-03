"""Pin the preregistered secondary count endpoint.

Registered 2026-09-03 before the tranche 02 candidate lock
(docs/evidence/mixed_gold_count_endpoint_preregistration_20260903.md). The
identity rule stays primary; this endpoint may never pass while identity
non-inferiority fails, and an identity miss or an abstention counts as the
full carrier error.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.codex_paper_eval.secondary_count_endpoint import (
    evaluate_secondary_count_endpoint,
    per_attempt_terms,
)

RULE = json.loads(
    (
        Path(__file__).resolve().parents[2]
        / "benchmarks/evaluation_tiers/mixed_gold/secondary_endpoints.json"
    ).read_text()
)


def _gold_root(
    tmp_path: Path, rows: dict[str, list[tuple[str, str, int | None]]]
) -> Path:
    root = tmp_path / "gold"
    root.mkdir()
    for gene, gold_rows in rows.items():
        with (root / f"{gene}_recall_input.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                [
                    "variant",
                    "pmid",
                    "carriers",
                    "affected",
                    "unaffected",
                    "gold_v2_carriers",
                    "gold_v2_affected",
                    "gold_v2_unaffected",
                    "gold_v2_status",
                    "gold_v2_note",
                    "gold_v2_source",
                ]
            )
            for variant, pmid, carriers in gold_rows:
                w.writerow(
                    [
                        variant,
                        pmid,
                        "" if carriers is None else carriers,
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                    ]
                )
    return root


def _report(pairs_by_paper: dict[tuple[str, str], list[tuple[str, str]]]) -> dict:
    return {
        "papers": [
            {
                "gene": gene,
                "pmid": pmid,
                "matched_variants": [{"gold": g, "predicted": p} for g, p in pairs],
            }
            for (gene, pmid), pairs in pairs_by_paper.items()
        ]
    }


def _predictions(values_by_paper: dict[tuple[str, str], dict[str, int | None]]) -> dict:
    return {
        "papers": [
            {
                "gene": gene,
                "pmid": pmid,
                "variants": [
                    {"variant": name, "carriers": value}
                    for name, value in values.items()
                ],
            }
            for (gene, pmid), values in values_by_paper.items()
        ]
    }


GOLD = {
    "KCNQ1": [
        ("A341V", "1", 10),
        ("R231C", "1", 4),
        ("G168R", "1", None),
        ("R518X", "2", 3),
    ]
}
IDENTITY_OK = {
    "recall": {
        "one_sided_lower_confidence_bound": 0.0,
        "delta_candidate_minus_baseline": 0.0,
    },
    "precision": {
        "one_sided_lower_confidence_bound": 0.0,
        "delta_candidate_minus_baseline": 0.0,
    },
}


def test_missed_or_abstained_rows_carry_the_full_gold_error(tmp_path):
    gold_root = _gold_root(tmp_path, GOLD)
    report = _report({("KCNQ1", "1"): [("A341V", "p.Ala341Val")], ("KCNQ1", "2"): []})
    predictions = _predictions(
        {("KCNQ1", "1"): {"p.Ala341Val": None}, ("KCNQ1", "2"): {}}
    )
    terms = per_attempt_terms(report, predictions, gold_root)
    # A341V matched but no value (error 10), R231C missed (4), G168R has no
    # asserted value (skipped), paper 2's R518X missed (3).
    assert terms[("KCNQ1", "1")] == {
        "abs_error": 14,
        "asserted": 2,
        "matched": 1,
        "matched_supplied": 0,
    }
    assert terms[("KCNQ1", "2")] == {
        "abs_error": 3,
        "asserted": 1,
        "matched": 0,
        "matched_supplied": 0,
    }


def test_candidate_that_supplies_correct_counts_passes_and_abstainer_fails(tmp_path):
    gold_root = _gold_root(tmp_path, GOLD)
    pairs = {
        ("KCNQ1", "1"): [("A341V", "p.Ala341Val"), ("R231C", "p.Arg231Cys")],
        ("KCNQ1", "2"): [("R518X", "p.Arg518*")],
    }
    baseline = _report(pairs)
    candidate = _report(pairs)
    base_pred = _predictions(
        {
            ("KCNQ1", "1"): {"p.Ala341Val": None, "p.Arg231Cys": None},
            ("KCNQ1", "2"): {"p.Arg518*": None},
        }
    )
    cand_pred = _predictions(
        {
            ("KCNQ1", "1"): {"p.Ala341Val": 10, "p.Arg231Cys": 4},
            ("KCNQ1", "2"): {"p.Arg518*": 3},
        }
    )
    result = evaluate_secondary_count_endpoint(
        baseline_report=baseline,
        candidate_report=candidate,
        baseline_predictions=base_pred,
        candidate_predictions=cand_pred,
        gold_root=gold_root,
        rule=RULE,
        identity_metrics=IDENTITY_OK,
    )
    assert result["observed"]["baseline"]["end_to_end_mae"] == pytest.approx(17 / 3)
    assert result["observed"]["candidate"]["end_to_end_mae"] == 0.0
    assert result["criteria"] == {
        "identity_guard_pass": True,
        "end_to_end_mae_pass": True,
        "coverage_pass": True,
    }
    assert result["passed"] is True
    # The reverse direction must fail on both the observed margin and the bound.
    reverse = evaluate_secondary_count_endpoint(
        baseline_report=baseline,
        candidate_report=candidate,
        baseline_predictions=cand_pred,
        candidate_predictions=base_pred,
        gold_root=gold_root,
        rule=RULE,
        identity_metrics=IDENTITY_OK,
    )
    assert reverse["passed"] is False
    assert reverse["criteria"]["end_to_end_mae_pass"] is False
    assert reverse["criteria"]["coverage_pass"] is False


def test_identity_guard_blocks_a_count_pass(tmp_path):
    gold_root = _gold_root(tmp_path, GOLD)
    pairs = {
        ("KCNQ1", "1"): [("A341V", "p.Ala341Val"), ("R231C", "p.Arg231Cys")],
        ("KCNQ1", "2"): [("R518X", "p.Arg518*")],
    }
    result = evaluate_secondary_count_endpoint(
        baseline_report=_report(pairs),
        candidate_report=_report(pairs),
        baseline_predictions=_predictions(
            {
                ("KCNQ1", "1"): {"p.Ala341Val": None, "p.Arg231Cys": None},
                ("KCNQ1", "2"): {"p.Arg518*": None},
            }
        ),
        candidate_predictions=_predictions(
            {
                ("KCNQ1", "1"): {"p.Ala341Val": 10, "p.Arg231Cys": 4},
                ("KCNQ1", "2"): {"p.Arg518*": 3},
            }
        ),
        gold_root=gold_root,
        rule=RULE,
        identity_metrics={
            "recall": {
                "one_sided_lower_confidence_bound": -0.03,
                "delta_candidate_minus_baseline": -0.01,
            },
            "precision": {
                "one_sided_lower_confidence_bound": 0.0,
                "delta_candidate_minus_baseline": 0.0,
            },
        },
    )
    assert result["criteria"]["end_to_end_mae_pass"] is True
    assert result["criteria"]["identity_guard_pass"] is False
    assert result["passed"] is False
    # A negative observed recall delta alone blocks the pass, even with a clean bound.
    dropped = evaluate_secondary_count_endpoint(
        baseline_report=_report(pairs),
        candidate_report=_report(pairs),
        baseline_predictions=_predictions(
            {
                ("KCNQ1", "1"): {"p.Ala341Val": None, "p.Arg231Cys": None},
                ("KCNQ1", "2"): {"p.Arg518*": None},
            }
        ),
        candidate_predictions=_predictions(
            {
                ("KCNQ1", "1"): {"p.Ala341Val": 10, "p.Arg231Cys": 4},
                ("KCNQ1", "2"): {"p.Arg518*": 3},
            }
        ),
        gold_root=gold_root,
        rule=RULE,
        identity_metrics={
            "recall": {
                "one_sided_lower_confidence_bound": -0.005,
                "delta_candidate_minus_baseline": -0.004,
            },
            "precision": {
                "one_sided_lower_confidence_bound": 0.0,
                "delta_candidate_minus_baseline": 0.0,
            },
        },
    )
    assert dropped["criteria"]["identity_guard_pass"] is False


def test_rule_file_matches_the_preregistration_text():
    assert RULE["end_to_end_mae"]["maximum_observed_delta"] == -0.05
    assert RULE["end_to_end_mae"]["upper_bound_must_be_below"] == 0.0
    assert RULE["coverage_on_matched"]["noninferiority_margin"] == -0.05
    assert RULE["identity_guard"]["recall_noninferiority_margin"] == -0.01
    assert RULE["identity_guard"]["precision_noninferiority_margin"] == -0.02
    assert RULE["identity_guard"]["observed_recall_delta_minimum"] == 0.0
    assert RULE["confidence_interval"]["seed"] == 2026090301
    assert RULE["confidence_interval"]["resamples"] == 10000
