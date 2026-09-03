"""Secondary, preregistered count endpoint for paired mixed-gold comparisons.

Registered text: ``docs/evidence/mixed_gold_count_endpoint_preregistration_20260903.md``
(committed before the tranche 02 candidate was locked). Rule parameters live in
``benchmarks/evaluation_tiers/mixed_gold/secondary_endpoints.json`` so the
registry digest bound into already-scored arms stays byte-identical.

Everything here is computed from artifacts that exist once both arms are
scored: each arm's locked ``report.json`` (which gold string paired with which
prediction string, per attempt), each arm's locked ``predictions.json`` (the
carrier value each prediction supplied), and the tranche's answer key (the
gold carrier value per row). No score is rewritten.

End-to-end carrier error, per gold row with an asserted carrier value:
``|gold − predicted|`` where ``predicted`` is the paired prediction's supplied
carrier value, or **0** when the row was not paired or the paired prediction
supplied no value. This mirrors ``cli/compare_variants.compute_end_to_end_count_error``.
"""

from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

COUNT_FIELD = "carriers"


def _load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _paper_key(paper: dict) -> tuple[str, str]:
    return str(paper.get("gene")), str(paper.get("pmid"))


def _prediction_values(
    predictions: dict,
) -> dict[tuple[str, str], dict[str, int | None]]:
    """Paper-lane prediction string → supplied carrier value (None when absent)."""
    out: dict[tuple[str, str], dict[str, int | None]] = {}
    for paper in predictions.get("papers") or []:
        key = _paper_key(paper)
        values: dict[str, int | None] = {}
        for variant in paper.get("variants") or []:
            name = str(variant.get("variant") or "")
            raw = variant.get(COUNT_FIELD)
            value: int | None
            try:
                value = int(raw) if raw is not None else None
            except (TypeError, ValueError):
                value = None
            # Twin rows merged before scoring keep one string; keep the first
            # supplied value and never let a later None erase it.
            if name not in values or (values[name] is None and value is not None):
                values[name] = value
        out[key] = values
    return out


def _gold_rows(gold_root: Path, gene: str, pmid: str) -> list[dict]:
    from benchmarks.codex_paper_eval.run_eval import load_gold

    return load_gold(gold_root, gene, pmid)


def per_attempt_terms(
    report: dict,
    predictions: dict,
    gold_root: Path,
) -> dict[tuple[str, str], dict[str, int]]:
    """Per attempt: sum of end-to-end carrier errors, asserted rows, matched
    rows, matched rows with a supplied carrier value."""
    pred_values = _prediction_values(predictions)
    terms: dict[tuple[str, str], dict[str, int]] = {}
    for paper in report.get("papers") or []:
        key = _paper_key(paper)
        pairs: dict[str, list[str]] = defaultdict(list)
        for pair in paper.get("matched_variants") or []:
            pairs[str(pair.get("gold"))].append(str(pair.get("predicted")))
        supplied_by_pred = pred_values.get(key, {})
        term = {"abs_error": 0, "asserted": 0, "matched": 0, "matched_supplied": 0}
        for gold_row in _gold_rows(gold_root, key[0], key[1]):
            gold_name = str(gold_row.get("variant") or "").strip()
            predicted_name = pairs[gold_name].pop(0) if pairs.get(gold_name) else None
            gold_value = gold_row.get(COUNT_FIELD)
            if predicted_name is not None:
                term["matched"] += 1
            supplied = (
                supplied_by_pred.get(predicted_name)
                if predicted_name is not None
                else None
            )
            if predicted_name is not None and supplied is not None:
                term["matched_supplied"] += 1
            if not isinstance(gold_value, int) or isinstance(gold_value, bool):
                continue
            term["asserted"] += 1
            term["abs_error"] += abs(
                gold_value - (supplied if supplied is not None else 0)
            )
        terms[key] = term
    return terms


def _rate(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def _percentile(values: list[float], level: float, *, upper: bool) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    rank = level if upper else (1.0 - level)
    index = min(len(ordered) - 1, max(0, int(round(rank * len(ordered))) - 1))
    return ordered[index]


def evaluate_secondary_count_endpoint(
    *,
    baseline_report: dict,
    candidate_report: dict,
    baseline_predictions: dict,
    candidate_predictions: dict,
    gold_root: Path,
    rule: dict,
    identity_metrics: dict,
) -> dict[str, Any]:
    """Apply the preregistered secondary count rule; never mutates inputs."""
    base_terms = per_attempt_terms(baseline_report, baseline_predictions, gold_root)
    cand_terms = per_attempt_terms(candidate_report, candidate_predictions, gold_root)
    if set(base_terms) != set(cand_terms):
        raise SystemExit("secondary endpoint: paired reports differ in attempts")
    for key in base_terms:
        if base_terms[key]["asserted"] != cand_terms[key]["asserted"]:
            raise SystemExit(f"secondary endpoint: asserted gold rows differ for {key}")

    clusters: dict[str, dict[str, dict[str, int]]] = defaultdict(
        lambda: {
            "baseline": {
                "abs_error": 0,
                "asserted": 0,
                "matched": 0,
                "matched_supplied": 0,
            },
            "candidate": {
                "abs_error": 0,
                "asserted": 0,
                "matched": 0,
                "matched_supplied": 0,
            },
        }
    )
    for key in base_terms:
        for arm, terms in (("baseline", base_terms), ("candidate", cand_terms)):
            for field, value in terms[key].items():
                clusters[key[1]][arm][field] += value

    def aggregate(cluster_ids: list[str]) -> dict[str, dict[str, float]]:
        sums = {
            arm: {"abs_error": 0, "asserted": 0, "matched": 0, "matched_supplied": 0}
            for arm in ("baseline", "candidate")
        }
        for cluster_id in cluster_ids:
            for arm in sums:
                for field, value in clusters[cluster_id][arm].items():
                    sums[arm][field] += value
        out: dict[str, dict[str, float]] = {}
        for arm, s in sums.items():
            out[arm] = {
                "end_to_end_mae": _rate(s["abs_error"], s["asserted"]),
                "coverage_on_matched": _rate(s["matched_supplied"], s["matched"]),
                "asserted_rows": s["asserted"],
                "matched_rows": s["matched"],
                "matched_supplied_rows": s["matched_supplied"],
            }
        return out

    ids = sorted(clusters)
    observed = aggregate(ids)
    mae_delta = (
        observed["candidate"]["end_to_end_mae"] - observed["baseline"]["end_to_end_mae"]
    )
    cov_delta = (
        observed["candidate"]["coverage_on_matched"]
        - observed["baseline"]["coverage_on_matched"]
    )

    boot = rule["confidence_interval"]
    rng = random.Random(int(boot["seed"]))
    mae_deltas: list[float] = []
    cov_deltas: list[float] = []
    for _ in range(int(boot["resamples"])):
        sampled = [rng.choice(ids) for _ in ids]
        agg = aggregate(sampled)
        mae_deltas.append(
            agg["candidate"]["end_to_end_mae"] - agg["baseline"]["end_to_end_mae"]
        )
        cov_deltas.append(
            agg["candidate"]["coverage_on_matched"]
            - agg["baseline"]["coverage_on_matched"]
        )
    level = float(rule["one_sided_confidence_level"])
    mae_upper = _percentile(mae_deltas, level, upper=True)
    cov_lower = _percentile(cov_deltas, level, upper=False)

    guard = rule["identity_guard"]
    recall_lb = float(identity_metrics["recall"]["one_sided_lower_confidence_bound"])
    precision_lb = float(
        identity_metrics["precision"]["one_sided_lower_confidence_bound"]
    )
    recall_observed = float(
        identity_metrics["recall"]["delta_candidate_minus_baseline"]
    )
    identity_guard_pass = (
        recall_lb >= float(guard["recall_noninferiority_margin"])
        and precision_lb >= float(guard["precision_noninferiority_margin"])
        and recall_observed >= float(guard.get("observed_recall_delta_minimum", 0.0))
    )
    mae_pass = mae_delta <= float(
        rule["end_to_end_mae"]["maximum_observed_delta"]
    ) and (mae_upper < float(rule["end_to_end_mae"]["upper_bound_must_be_below"]))
    coverage_pass = cov_delta >= float(
        rule["coverage_on_matched"]["minimum_observed_delta"]
    ) and (cov_lower >= float(rule["coverage_on_matched"]["noninferiority_margin"]))
    passed = identity_guard_pass and mae_pass and coverage_pass
    return {
        "endpoint": "paper_derived_end_to_end_carrier_mae",
        "role": "secondary_preregistered_2026-09-03",
        "count_field": COUNT_FIELD,
        "cluster_unit": "PMID",
        "cluster_count": len(ids),
        "observed": observed,
        "deltas": {
            "end_to_end_mae": {
                "delta_candidate_minus_baseline": mae_delta,
                "one_sided_upper_confidence_bound": mae_upper,
            },
            "coverage_on_matched": {
                "delta_candidate_minus_baseline": cov_delta,
                "one_sided_lower_confidence_bound": cov_lower,
            },
        },
        "criteria": {
            "identity_guard_pass": identity_guard_pass,
            "end_to_end_mae_pass": mae_pass,
            "coverage_pass": coverage_pass,
        },
        "passed": passed,
        "registered_rule": rule,
    }
