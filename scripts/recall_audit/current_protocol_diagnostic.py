#!/usr/bin/env python3
"""Build the current two-cohort recall and count diagnostic.

The diagnostic is deliberately read-only with respect to scored runs.  It
combines their locked reports with ``fn_root_cause.py`` outputs, recomputes the
preregistered end-to-end carrier endpoint, and (optionally) replays the current
trust gate on temporary database copies.  It never rewrites a lock or score.

Example
-------
python scripts/recall_audit/current_protocol_diagnostic.py \
  --out-dir docs/evidence/current_protocol_diagnostic_20260904 \
  --generated-at 2026-09-04T16:45:00Z \
  --replay-trust-gate
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks.codex_paper_eval.run_eval import load_gold  # noqa: E402
from benchmarks.codex_paper_eval.secondary_count_endpoint import (  # noqa: E402
    _prediction_values,
)
from pipeline.trust_gate import apply_trust_gate  # noqa: E402


@dataclass(frozen=True)
class Cohort:
    cohort_id: str
    label: str
    score_lane: str
    run_dir: Path
    gold_root: Path
    root_cause_tsv: Path


COHORTS = (
    Cohort(
        "gold_120_current",
        "gold_120 current lock (118 attempts)",
        "trusted legacy projection (paper plus ClinVar/PubTator linkage)",
        REPO / "benchmarks/codex_paper_eval/runs/20260902_false_zero_recovery_gold118",
        REPO / "gene_variant_fetcher_gold_standard/normalized",
        Path("gold118_fn_root_cause/fn_root_cause.tsv"),
    ),
    Cohort(
        "mixed_continuation_120",
        "mixed continuation tranche 01 (120 attempts)",
        "paper-derived primary",
        REPO
        / "benchmarks/codex_paper_eval/runs/20260903_protocol_cont120_01_candidate",
        REPO / "benchmarks/evaluation_tiers/mixed_gold_continuation_120/answer_key",
        Path("cont120_fn_root_cause/fn_root_cause.tsv"),
    ),
)


CAUSE_GROUPS = {
    "acquisition": "No usable source / row absent from acquired text",
    "unknown_notation": "Figure-only or notation not string-searchable",
    "source_selection": "Source staged incorrectly",
    "condensing": "Source lost during condensing",
    "model_missed": "Model omitted visible row",
    "parser_dropped": "Parser dropped model output",
    "postprocess_dropped": "Post-processing dropped row",
    "postprocess_dropped(quarantined)": "Post-processing quarantined row",
    "paper_row_lost_to_linkage_origin": "Paper row re-attributed to linkage",
    "projection_dropped": "Trusted projection dropped row",
    "matcher": "Scorer did not pair projected row",
}

SOURCE_CLASS_LABELS = {
    "source_absent": "no acquired source",
    "text_absent_stub_body": "abstract/stub only",
    "text_absent_substitution": "substitution absent from acquired bytes",
    "text_absent_figures_present": "likely figure-only",
    "text_absent_notation_inconclusive": "notation not searchable",
    "present_in_body": "present in body; later-stage loss",
}

COUNT_CAUSE_LABELS = {
    "identity_missed": "Identity missed",
    "matched_count_missing": "Matched identity, count omitted",
    "supplied_wrong": "Count supplied, wrong value",
    "supplied_exact": "Count supplied, exact",
}


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty dataset: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _load_root_causes(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _paper_key(paper: dict[str, Any]) -> tuple[str, str]:
    return str(paper.get("gene")), str(paper.get("pmid"))


def _count_diagnostic(
    report: dict[str, Any], predictions: dict[str, Any], gold_root: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, float | int]]:
    """Return error decomposition, ranked papers, and carrier endpoint stats."""
    predicted_values = _prediction_values(predictions)
    aggregate: dict[str, dict[str, int]] = defaultdict(
        lambda: {"rows": 0, "abs_error": 0}
    )
    by_paper: dict[tuple[str, str], dict[str, dict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: {"rows": 0, "abs_error": 0})
    )
    strata = {
        "all": {"rows": 0, "abs_error": 0},
        "zero": {"rows": 0, "abs_error": 0},
        "nonzero": {"rows": 0, "abs_error": 0},
    }
    for paper in report.get("papers") or []:
        key = _paper_key(paper)
        pairs: dict[str, list[str]] = defaultdict(list)
        for pair in paper.get("matched_variants") or []:
            pairs[str(pair.get("gold"))].append(str(pair.get("predicted")))
        for gold_row in load_gold(gold_root, *key):
            gold_value = gold_row.get("carriers")
            if not isinstance(gold_value, int) or isinstance(gold_value, bool):
                continue
            gold_name = str(gold_row.get("variant") or "").strip()
            predicted_name = pairs[gold_name].pop(0) if pairs.get(gold_name) else None
            predicted_value = (
                predicted_values.get(key, {}).get(predicted_name)
                if predicted_name is not None
                else None
            )
            if predicted_name is None:
                cause = "identity_missed"
            elif predicted_value is None:
                cause = "matched_count_missing"
            elif predicted_value == gold_value:
                cause = "supplied_exact"
            else:
                cause = "supplied_wrong"
            error = abs(gold_value - (predicted_value or 0))
            stratum = "zero" if gold_value == 0 else "nonzero"
            strata["all"]["rows"] += 1
            strata["all"]["abs_error"] += error
            strata[stratum]["rows"] += 1
            strata[stratum]["abs_error"] += error
            aggregate[cause]["rows"] += 1
            aggregate[cause]["abs_error"] += error
            by_paper[key][cause]["rows"] += 1
            by_paper[key][cause]["abs_error"] += error

    total_error = sum(x["abs_error"] for x in aggregate.values())
    total_rows = sum(x["rows"] for x in aggregate.values())
    decomposition = []
    for cause in (
        "identity_missed",
        "matched_count_missing",
        "supplied_wrong",
        "supplied_exact",
    ):
        values = aggregate[cause]
        decomposition.append(
            {
                "cause": cause,
                "cause_label": COUNT_CAUSE_LABELS[cause],
                "rows": values["rows"],
                "absolute_error_units": values["abs_error"],
                "share_of_error": (
                    values["abs_error"] / total_error if total_error else 0.0
                ),
            }
        )

    ranked = []
    for (gene, pmid), causes in by_paper.items():
        paper_error = sum(x["abs_error"] for x in causes.values())
        if not paper_error:
            continue
        dominant = max(causes, key=lambda cause: causes[cause]["abs_error"])
        ranked.append(
            {
                "gene": gene,
                "pmid": pmid,
                "absolute_error_units": paper_error,
                "asserted_rows": sum(x["rows"] for x in causes.values()),
                "dominant_cause": COUNT_CAUSE_LABELS[dominant],
                "identity_miss_error": causes["identity_missed"]["abs_error"],
                "omitted_count_error": causes["matched_count_missing"]["abs_error"],
                "wrong_count_error": causes["supplied_wrong"]["abs_error"],
            }
        )
    ranked.sort(
        key=lambda row: (-row["absolute_error_units"], row["gene"], row["pmid"])
    )
    endpoint = {
        "end_to_end_mae": total_error / total_rows if total_rows else 0.0,
        "end_to_end_mae_nonzero_gold": (
            strata["nonzero"]["abs_error"] / strata["nonzero"]["rows"]
            if strata["nonzero"]["rows"]
            else 0.0
        ),
        "end_to_end_mae_zero_gold": (
            strata["zero"]["abs_error"] / strata["zero"]["rows"]
            if strata["zero"]["rows"]
            else 0.0
        ),
        "gold_zero_rows": strata["zero"]["rows"],
        "gold_nonzero_rows": strata["nonzero"]["rows"],
        "absolute_error_units": total_error,
    }
    return decomposition, ranked, endpoint


def _format_breakdown(counter: Counter[str], labels: dict[str, str]) -> str:
    return "; ".join(
        f"{count} {labels.get(key, key.replace('_', ' '))}"
        for key, count in counter.most_common()
    )


def _trust_gate_replay(run_dir: Path) -> dict[str, int]:
    """Replay the gate on temporary DB copies and count material changes."""
    totals = {
        "databases": 0,
        "source_backed_outliers": 0,
        "materially_changed_rows": 0,
    }
    for status_path in sorted(
        (run_dir / "production_runs").glob("*/*/RUN_STATUS.json")
    ):
        status = _json(status_path)
        active_db = status.get("active_db")
        if not active_db:
            continue
        source = status_path.parent / str(active_db)
        if not source.is_file():
            continue
        totals["databases"] += 1
        with tempfile.TemporaryDirectory(prefix="gvf_current_protocol_") as temp_dir:
            copied = Path(temp_dir) / source.name
            shutil.copy2(source, copied)
            con = sqlite3.connect(copied)
            con.row_factory = sqlite3.Row
            columns = {
                row[1] for row in con.execute("pragma table_info(penetrance_data)")
            }
            material_columns = [
                column
                for column in (
                    "penetrance_id",
                    "trust_tier",
                    "trust_reasons",
                    "field_trust",
                    "total_carriers_observed",
                    "affected_count",
                    "unaffected_count",
                )
                if column in columns
            ]
            before = {
                row["penetrance_id"]: tuple(
                    row[column] for column in material_columns[1:]
                )
                for row in con.execute(
                    f"select {','.join(material_columns)} from penetrance_data"
                )
            }
            con.close()
            stats = apply_trust_gate(copied)
            totals["source_backed_outliers"] += int(
                stats.get("paper_outlier_source_supported", 0)
            )
            con = sqlite3.connect(copied)
            con.row_factory = sqlite3.Row
            after = {
                row["penetrance_id"]: tuple(
                    row[column] for column in material_columns[1:]
                )
                for row in con.execute(
                    f"select {','.join(material_columns)} from penetrance_data"
                )
            }
            con.close()
            totals["materially_changed_rows"] += sum(
                before.get(row_id) != values for row_id, values in after.items()
            )
    return totals


def _solution_test(
    name: str, cohort: str, comparison: dict[str, Any]
) -> dict[str, Any]:
    count = comparison["secondary_count_endpoint"]
    return {
        "test": name,
        "evaluation_set": cohort,
        "attempts": comparison["attempt_count"],
        "recall_baseline": comparison["metrics"]["recall"]["baseline"],
        "recall_candidate": comparison["metrics"]["recall"]["candidate"],
        "recall_delta_pp": 100
        * comparison["metrics"]["recall"]["delta_candidate_minus_baseline"],
        "carrier_e2e_mae_baseline": count["observed"]["baseline"]["end_to_end_mae"],
        "carrier_e2e_mae_candidate": count["observed"]["candidate"]["end_to_end_mae"],
        "carrier_e2e_mae_delta": count["deltas"]["end_to_end_mae"][
            "delta_candidate_minus_baseline"
        ],
        "carrier_coverage_delta_pp": 100
        * count["deltas"]["coverage_on_matched"]["delta_candidate_minus_baseline"],
        "identity_decision": comparison["decision"],
        "count_endpoint_passed": count["passed"],
    }


def _gold_120b_status() -> dict[str, Any]:
    registry = _json(REPO / "benchmarks/evaluation_tiers/registry.json")
    tier = next(row for row in registry["tiers"] if row["id"] == "gold_120b")
    scored_runs = []
    for report_path in (REPO / "benchmarks/codex_paper_eval/runs").glob(
        "*/report.json"
    ):
        report = _json(report_path)
        manifest = str((report.get("selection") or {}).get("paper_manifest") or "")
        if Path(manifest).name == tier["manifest"]:
            scored_runs.append(report_path.parent.name)
    return {
        "tier": "gold_120b",
        "registered_attempts": tier["attempt_count"],
        "registered_unique_pmids": tier["unique_pmid_count"],
        "matching_scored_runs": scored_runs,
        "status": "scored" if scored_runs else "registered but no locked report found",
    }


def build(
    out_dir: Path, *, generated_at: str, replay_trust_gate: bool
) -> dict[str, Any]:
    cohort_metrics: list[dict[str, Any]] = []
    gene_metrics: list[dict[str, Any]] = []
    fn_causes: list[dict[str, Any]] = []
    problem_papers: list[dict[str, Any]] = []
    count_decomposition: list[dict[str, Any]] = []
    mae_problem_papers: list[dict[str, Any]] = []
    counterfactuals: list[dict[str, Any]] = []
    trust_replays: list[dict[str, Any]] = []
    sources: dict[str, str] = {}

    for cohort in COHORTS:
        report_path = cohort.run_dir / "report.json"
        predictions_path = cohort.run_dir / "predictions.json"
        root_path = out_dir / cohort.root_cause_tsv
        report = _json(report_path)
        predictions = _json(predictions_path)
        root_rows = _load_root_causes(root_path)
        paper_lookup = {_paper_key(paper): paper for paper in report["papers"]}
        by_leaf = Counter(row["leaf"] for row in root_rows)
        by_paper: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in root_rows:
            by_paper[(row["gene"], row["pmid"])].append(row)

        decomposition, mae_papers, count_endpoint = _count_diagnostic(
            report, predictions, cohort.gold_root
        )
        overall = report["overall"]
        carrier = overall["count"]["carriers"]
        cohort_metrics.append(
            {
                "cohort_id": cohort.cohort_id,
                "cohort": cohort.label,
                "run_id": report["run_id"],
                "score_lane": cohort.score_lane,
                "attempts": overall["papers"],
                "unique_pmids": len({str(p["pmid"]) for p in report["papers"]}),
                "gold_rows": overall["tp"] + overall["fn"],
                "tp": overall["tp"],
                "fp": overall["fp"],
                "fn": overall["fn"],
                "precision": overall["precision"],
                "recall": overall["recall"],
                "f1": overall["f1"],
                "carrier_values_supplied": carrier["predicted"],
                "carrier_supply_rate": carrier["recall"],
                "carrier_conditional_mae": carrier["mae"],
                "carrier_end_to_end_mae": count_endpoint["end_to_end_mae"],
                "carrier_end_to_end_mae_nonzero_gold": count_endpoint[
                    "end_to_end_mae_nonzero_gold"
                ],
                "carrier_end_to_end_mae_zero_gold": count_endpoint[
                    "end_to_end_mae_zero_gold"
                ],
                "carrier_gold_zero_rows": count_endpoint["gold_zero_rows"],
                "carrier_gold_nonzero_rows": count_endpoint["gold_nonzero_rows"],
                "affected_supply_rate": overall["count"]["affected"]["recall"],
                "unaffected_supply_rate": overall["count"]["unaffected"]["recall"],
            }
        )
        for gene, values in sorted(report["by_gene"].items()):
            gene_metrics.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort": cohort.label,
                    "gene": gene,
                    "tp": values["tp"],
                    "fp": values["fp"],
                    "fn": values["fn"],
                    "precision": values["precision"],
                    "recall": values["recall"],
                    "carrier_supply_rate": values["count"]["carriers"]["recall"],
                    "carrier_conditional_mae": values["count"]["carriers"]["mae"],
                }
            )
        for leaf, count in by_leaf.most_common():
            fn_causes.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort": cohort.label,
                    "cause": leaf,
                    "cause_label": CAUSE_GROUPS.get(leaf, leaf.replace("_", " ")),
                    "missed_rows": count,
                    "share_of_fn": count / overall["fn"],
                    "gold_rows": overall["tp"] + overall["fn"],
                }
            )

        ranked_papers = sorted(
            by_paper.items(), key=lambda item: (-len(item[1]), item[0])
        )
        for rank, ((gene, pmid), rows) in enumerate(ranked_papers[:12], start=1):
            leaves = Counter(row["leaf"] for row in rows)
            source_classes = Counter(row["sweep_class"] for row in rows)
            paper = paper_lookup[(gene, pmid)]
            problem_papers.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort": cohort.label,
                    "rank_in_cohort": rank,
                    "gene": gene,
                    "pmid": pmid,
                    "false_negatives": len(rows),
                    "true_positives": paper["tp"],
                    "false_positives": paper["fp"],
                    "primary_bottleneck": CAUSE_GROUPS.get(
                        leaves.most_common(1)[0][0], leaves.most_common(1)[0][0]
                    ),
                    "stage_breakdown": _format_breakdown(leaves, CAUSE_GROUPS),
                    "source_breakdown": _format_breakdown(
                        source_classes, SOURCE_CLASS_LABELS
                    ),
                }
            )
        for row in decomposition:
            count_decomposition.append(
                {"cohort_id": cohort.cohort_id, "cohort": cohort.label, **row}
            )
        for rank, row in enumerate(mae_papers[:12], start=1):
            mae_problem_papers.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort": cohort.label,
                    "rank_in_cohort": rank,
                    **row,
                }
            )

        gold_rows = overall["tp"] + overall["fn"]
        reachable = overall["fn"] - by_leaf["acquisition"] - by_leaf["unknown_notation"]
        scenarios = (
            ("Current", 0),
            ("Perfect recovery of reachable pipeline misses", reachable),
            ("Perfect recovery of hard acquisition misses", by_leaf["acquisition"]),
            (
                "Perfect recovery of hard + figure/notation source misses",
                by_leaf["acquisition"] + by_leaf["unknown_notation"],
            ),
        )
        for scenario, recovered in scenarios:
            counterfactuals.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort": cohort.label,
                    "scenario": scenario,
                    "rows_recovered": recovered,
                    "recall_upper_bound": (overall["tp"] + recovered) / gold_rows,
                    "note": "Diagnostic upper bound; not a measured candidate result",
                }
            )
        if replay_trust_gate:
            trust_replays.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort": cohort.label,
                    **_trust_gate_replay(cohort.run_dir),
                }
            )
        sources[str(report_path.relative_to(REPO))] = _sha256(report_path)
        sources[str(predictions_path.relative_to(REPO))] = _sha256(predictions_path)
        sources[str(root_path.relative_to(REPO))] = _sha256(root_path)

    comparison_paths = (
        (
            "v1 paper-agnostic reading/count bundle",
            "mixed tranche 01 (49 attempts)",
            REPO
            / "benchmarks/codex_paper_eval/runs/20260903_protocol_mixed01_candidate/compare_discovery_with_secondary.json",
        ),
        (
            "v2 paper-agnostic reading/count bundle",
            "mixed tranche 02 (49 attempts)",
            REPO
            / "benchmarks/codex_paper_eval/runs/20260903_protocol_mixed02_candidate/compare_discovery_with_secondary.json",
        ),
        (
            "v2 paper-agnostic reading/count bundle",
            "mixed continuation tranche 01 (120 attempts)",
            REPO
            / "benchmarks/codex_paper_eval/runs/20260903_protocol_cont120_01_candidate/compare_discovery_with_secondary.json",
        ),
    )
    solution_tests = []
    for name, cohort_label, path in comparison_paths:
        solution_tests.append(_solution_test(name, cohort_label, _json(path)))
        sources[str(path.relative_to(REPO))] = _sha256(path)
    for replay in trust_replays:
        metric = next(
            row for row in cohort_metrics if row["cohort_id"] == replay["cohort_id"]
        )
        solution_tests.append(
            {
                "test": "tg6 source-backed table-outlier exception",
                "evaluation_set": replay["cohort"],
                "attempts": metric["attempts"],
                "recall_baseline": metric["recall"],
                "recall_candidate": metric["recall"],
                "recall_delta_pp": 0.0,
                "carrier_e2e_mae_baseline": metric["carrier_end_to_end_mae"],
                "carrier_e2e_mae_candidate": metric["carrier_end_to_end_mae"],
                "carrier_e2e_mae_delta": 0.0,
                "carrier_coverage_delta_pp": 0.0,
                "identity_decision": (
                    "no aggregate effect: "
                    f"{replay['source_backed_outliers']} eligible rows, "
                    f"{replay['materially_changed_rows']} material DB changes"
                ),
                "count_endpoint_passed": None,
            }
        )

    scope_status = _gold_120b_status()
    registry_path = REPO / "benchmarks/evaluation_tiers/registry.json"
    sources[str(registry_path.relative_to(REPO))] = _sha256(registry_path)
    payload = {
        "schema_version": 1,
        "generated_at": generated_at,
        "scope": {
            "audited_cohorts": [cohort.cohort_id for cohort in COHORTS],
            "gold_120b": scope_status,
            "note": (
                "The two audited cohorts are the latest locked, scored cohorts near "
                "120 attempts. The registered gold_120b tier is reported separately "
                "because no matching locked report exists."
            ),
        },
        "cohort_metrics": cohort_metrics,
        "gene_metrics": gene_metrics,
        "fn_causes": fn_causes,
        "problem_papers": problem_papers,
        "count_error_decomposition": count_decomposition,
        "mae_problem_papers": mae_problem_papers,
        "recall_counterfactuals": counterfactuals,
        "solution_tests": solution_tests,
        "trust_gate_replays": trust_replays,
        "sources": [
            {"path": path, "sha256": digest} for path, digest in sorted(sources.items())
        ],
    }
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--generated-at",
        default=datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    )
    parser.add_argument("--replay-trust-gate", action="store_true")
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = build(
        out_dir,
        generated_at=args.generated_at,
        replay_trust_gate=args.replay_trust_gate,
    )
    (out_dir / "diagnostic.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    for key in (
        "cohort_metrics",
        "gene_metrics",
        "fn_causes",
        "problem_papers",
        "count_error_decomposition",
        "mae_problem_papers",
        "recall_counterfactuals",
        "solution_tests",
        "trust_gate_replays",
    ):
        rows = payload[key]
        if rows:
            _write_csv(out_dir / f"{key}.csv", rows)
    print(
        json.dumps(
            {
                "out_dir": str(out_dir),
                "rows": {k: len(v) for k, v in payload.items() if isinstance(v, list)},
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
