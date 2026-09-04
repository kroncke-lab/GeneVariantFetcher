#!/usr/bin/env python3
"""Companion difference view of the phenotype-count recovery figure family.

The canonical per-protocol-change figure is the stratified phenotype-count
recovery figure (``scripts/build_stratified_phenotype_count_recovery.py`` on top
of the base plotting module ``scripts/build_phenotype_count_recovery.py`` and the
opened-tranche combiner ``scripts/build_combined_phenotype_count_recovery.py``;
``run_eval.py score`` regenerates it after every scored candidate arm). That
figure is an agreement plot over identity-matched rows. This script renders the
same evaluation as a difference plot for one locked run, reusing the base
module's palette, text and marker-size conventions, and adds what the
agreement plot cannot show: the identity misses, and the frozen baseline arm on
the identical gold rows.

For every gold row with an asserted count, the figure plots the automated count
minus the reference count against the reference count, one panel per measure
(carriers, affected, unaffected). A gold row the protocol never matched
(identity miss) and a matched row with no automated count (abstention) are both
evaluated as zero, so they sit on the lower boundary ``difference = -reference``
in their own marker styles; supplied rows sit on the zero line when exact. This
is the end-to-end definition preregistered for the secondary count endpoint
(``benchmarks/codex_paper_eval/secondary_count_endpoint.py``) applied to all
three count measures, so the figure shows how well the whole process works
against the gold standard rather than agreement conditional on count supply.

When the run is the candidate arm of a paired tranche, the frozen baseline arm
is drawn on the identical gold rows above it and the paired deltas (with the
registered bootstrap bounds when a ``run_eval.py compare`` JSON exists) are
printed under the panels, so a protocol change is read as the difference
between two difference plots.

No plotting package is required. The script writes an editable SVG, a tidy CSV
and a JSON summary straight from the locked artifacts, then best-effort PNG/PDF
copies. ``run_eval.py score`` renders it for every scored run and
``run_eval.py compare`` re-renders the candidate with the registered bounds;
run it by hand to regenerate a figure or to copy one into ``docs/figures``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import build_phenotype_count_recovery as base  # noqa: E402

RUNS_DIR = REPO / "benchmarks" / "codex_paper_eval" / "runs"
DEFAULT_RUN = RUNS_DIR / "20260824_aha_table_sourcebound_gold118"
DEFAULT_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
FIGURE_STEM = "gold_difference"
FIELDS = ("carriers", "affected", "unaffected")
FIELD_TITLES = {
    "carriers": "Carriers",
    "affected": "Affected individuals",
    "unaffected": "Unaffected individuals",
}
STATUSES = ("identity_miss", "abstained", "supplied")
STATUS_LABELS = {
    "supplied": "automated count supplied",
    "abstained": "variant matched, no automated count (scored 0)",
    "identity_miss": "variant not found (scored 0)",
}
EVALUATION_RULE = (
    "difference = automated count minus reference count for every gold row "
    "with an asserted value; an identity miss or a matched row without an "
    "automated count is evaluated as 0 (the full reference count is the error)"
)

# Palette, text element and marker-size rule come from the base plotting module
# so the companion view and the canonical figure keep one look.
COLORS = dict(base.COLORS)
text_element = base.text_element
bubble_radius = base.bubble_radius
STATUS_STYLE = {
    "supplied": {
        "fill": COLORS["blue"],
        "stroke": COLORS["navy"],
        "stroke_width": 2.5,
        "opacity": 0.9,
        "dash": "",
    },
    "abstained": {
        "fill": COLORS["coral_pale"],
        "stroke": COLORS["coral"],
        "stroke_width": 3.0,
        "opacity": 0.96,
        "dash": "",
    },
    "identity_miss": {
        "fill": COLORS["pale"],
        "stroke": COLORS["slate"],
        "stroke_width": 3.0,
        "opacity": 0.96,
        "dash": "7 5",
    },
}

WIDTH = 2400
PNG_MAX_WIDTH = 2160
PANEL = 620.0
PANEL_X = (170.0, 940.0, 1710.0)
INSET = 46.0
FIRST_ROW_TOP = 400.0
ROW_PITCH = 1010.0
TICK_VALUES = (0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
MARKER_SIZE_KEY_COUNTS = (1, 5, 25)


# --------------------------------------------------------------------------- data
def _run_eval():
    """Import the benchmark scorer lazily; it is the matching/gold authority."""
    from benchmarks.codex_paper_eval import run_eval

    return run_eval


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def to_int(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(str(value)))
        except (TypeError, ValueError):
            return None


def verify_lock(run_dir: Path) -> dict[str, Any]:
    """Refuse to plot a run whose selection or predictions moved after locking."""
    lock_path = run_dir / "LOCK.json"
    if not lock_path.is_file():
        raise SystemExit(f"{run_dir} has no LOCK.json; only locked runs are plotted")
    lock = read_json(lock_path)
    checked = {
        "selection_sha256": sha256(run_dir / "selection.json"),
        "predictions_sha256": sha256(run_dir / "predictions.json"),
    }
    for key, observed in checked.items():
        if observed != lock.get(key):
            raise SystemExit(f"lock mismatch for {key}: {observed} != {lock.get(key)}")
    return {**lock, "verified_hashes": checked}


def resolve_gold_sources(
    report: dict[str, Any], gold_root: Path | None
) -> dict[str, Path]:
    """Gold CSV per gene: an explicit root wins, then the paths the score recorded."""
    run_eval = _run_eval()
    genes = sorted({str(paper["gene"]) for paper in report["papers"]})
    recorded = report.get("gold_sources") or {}
    sources: dict[str, Path] = {}
    for gene in genes:
        if gold_root is not None:
            sources[gene] = run_eval.gold_csv_path(gold_root, gene)
        elif recorded.get(gene):
            sources[gene] = Path(str(recorded[gene]))
        else:
            sources[gene] = run_eval.gold_csv_path(DEFAULT_GOLD, gene)
    missing = [str(path) for path in sources.values() if not path.is_file()]
    if missing:
        raise SystemExit("gold source missing: " + ", ".join(missing))
    return sources


def load_gold_rows(path: Path, pmid: str) -> list[dict[str, Any]]:
    """The scorer's gold rows for one paper, plus which value set each came from."""
    run_eval = _run_eval()
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for source in csv.DictReader(handle):
            if str(source.get("pmid", "")).strip() != pmid:
                continue
            if run_eval.gold_row_excluded(source):
                continue
            status = str(source.get("gold_v2_status") or "").strip()
            rows.append(
                {
                    "variant": str(source.get("variant", "")).strip(),
                    **{
                        field: run_eval.authoritative_gold_count(
                            source, field, parser=run_eval.to_int
                        )
                        for field in FIELDS
                    },
                    "gold_value_set": "gold_v2" if status else "baseline",
                    "gold_v2_status": status,
                }
            )
    return rows


def predicted_count_rows(
    paper_prediction: dict[str, Any], gene: str
) -> dict[str, dict[str, Any]]:
    """Primary-lane prediction string -> the counts the scorer saw for it.

    The scorer merges notation twins before matching, which can move a count
    from a twin row onto the kept row, so the same merge is reproduced here.
    """
    run_eval = _run_eval()
    merged, _twins = run_eval.merge_notation_twins(
        list(paper_prediction.get("variants") or []), gene
    )
    by_name: dict[str, dict[str, Any]] = {}
    for row in merged:
        name = str(row.get("variant") or "")
        home = by_name.setdefault(
            name, {"row": row, **{field: None for field in FIELDS}}
        )
        for field in FIELDS:
            value = to_int(row.get(field))
            if home[field] is None and value is not None:
                home[field] = value
    return by_name


def build_rows(
    arm: str,
    report: dict[str, Any],
    predictions: dict[str, Any],
    gold_sources: dict[str, Path],
) -> list[dict[str, Any]]:
    """One row per (gold row, measure) with an asserted reference count."""
    prediction_map = {
        (str(paper["gene"]), str(paper["pmid"])): paper
        for paper in predictions["papers"]
    }
    rows: list[dict[str, Any]] = []
    matched_total = 0
    missed_total = 0
    for paper in report["papers"]:
        gene, pmid = str(paper["gene"]), str(paper["pmid"])
        predicted = prediction_map.get((gene, pmid))
        if predicted is None:
            raise SystemExit(f"locked predictions lack {gene} PMID {pmid}")
        if "matched_variants" not in paper:
            raise SystemExit(
                f"report for {gene} PMID {pmid} has no matched_variants; "
                "re-score the run with the current scorer first"
            )
        gold = load_gold_rows(gold_sources[gene], pmid)
        pairs: dict[str, list[str]] = defaultdict(list)
        for pair in paper.get("matched_variants") or []:
            pairs[str(pair.get("gold"))].append(str(pair.get("predicted")))
        if sum(len(names) for names in pairs.values()) != int(paper["tp"]):
            raise SystemExit(f"{gene} PMID {pmid}: matched pairs do not equal TP")
        counts = predicted_count_rows(predicted, gene)
        matched_here = 0
        for index, gold_row in enumerate(gold):
            queue = pairs.get(gold_row["variant"])
            predicted_name = queue.pop(0) if queue else None
            supplied: dict[str, Any] | None = None
            if predicted_name is not None:
                matched_here += 1
                supplied = counts.get(predicted_name)
                if supplied is None:
                    raise SystemExit(
                        f"{gene} PMID {pmid}: matched prediction {predicted_name!r} "
                        "is absent from the locked predictions"
                    )
            for field in FIELDS:
                gold_value = gold_row[field]
                if gold_value is None:
                    continue
                if supplied is None:
                    status, raw = "identity_miss", None
                else:
                    raw = supplied[field]
                    status = "supplied" if raw is not None else "abstained"
                evaluated = raw if raw is not None else 0
                difference = evaluated - int(gold_value)
                prediction_row = supplied["row"] if supplied else {}
                rows.append(
                    {
                        "analysis_run": str(report["run_id"]),
                        "arm": arm,
                        "gene": gene,
                        "pmid": pmid,
                        "gold_row_index": index,
                        "gold_variant": gold_row["variant"],
                        "predicted_variant": predicted_name or "",
                        "measure": field,
                        "gold_count": int(gold_value),
                        "automated_count_raw": raw,
                        "automated_count_evaluated": evaluated,
                        "status": status,
                        "difference": difference,
                        "absolute_difference": abs(difference),
                        "exact": int(difference == 0),
                        "gold_value_set": gold_row["gold_value_set"],
                        "gold_v2_status": gold_row["gold_v2_status"],
                        "gold_provenance": str(paper.get("gold_provenance") or ""),
                        "prediction_source_layer": str(
                            prediction_row.get("source_layer") or ""
                        ),
                        "prediction_source_location": str(
                            prediction_row.get("source_location") or ""
                        ),
                    }
                )
        if any(pairs.values()):
            leftovers = [name for name, names in pairs.items() if names]
            raise SystemExit(
                f"{gene} PMID {pmid}: matched gold strings not in gold: {leftovers}"
            )
        missed_here = len(gold) - matched_here
        if matched_here != int(paper["tp"]) or missed_here != int(paper["fn"]):
            raise SystemExit(f"{gene} PMID {pmid}: TP/FN do not reproduce")
        matched_total += matched_here
        missed_total += missed_here
    overall = report["overall"]
    if matched_total != int(overall["tp"]) or missed_total != int(overall["fn"]):
        raise SystemExit("identity totals do not reproduce the locked report")
    return rows


def _percentile(values: list[float], level: float) -> float:
    ordered = sorted(values)
    rank = min(len(ordered), max(1, math.ceil(level * len(ordered))))
    return ordered[rank - 1]


def summarize_field(
    rows: list[dict[str, Any]], report_metric: dict[str, Any] | None
) -> dict[str, Any]:
    """End-to-end difference summary for one measure, checked against the report."""
    if not rows:
        return {"gold_rows": 0}
    differences = [int(row["difference"]) for row in rows]
    absolute = [abs(value) for value in differences]
    by_status = Counter(str(row["status"]) for row in rows)
    matched = by_status["supplied"] + by_status["abstained"]
    supplied = [row for row in rows if row["status"] == "supplied"]
    supplied_abs = [int(row["absolute_difference"]) for row in supplied]
    exact = sum(value == 0 for value in differences)
    zero_baseline_exact = sum(
        1 for row in rows if row["gold_count"] == 0 and row["status"] != "supplied"
    )
    mean_difference = statistics.fmean(differences)
    sd_difference = statistics.pstdev(differences) if len(differences) > 1 else 0.0
    record: dict[str, Any] = {
        "gold_rows": len(rows),
        "gold_rows_zero": sum(row["gold_count"] == 0 for row in rows),
        "gold_max": max(row["gold_count"] for row in rows),
        "identity_matched": matched,
        "identity_miss": by_status["identity_miss"],
        "identity_coverage": matched / len(rows),
        "supplied": by_status["supplied"],
        "abstained": by_status["abstained"],
        "supplied_fraction": by_status["supplied"] / len(rows),
        "coverage_on_matched": (by_status["supplied"] / matched) if matched else None,
        "exact": exact,
        "exact_fraction": exact / len(rows),
        "zero_baseline_exact": zero_baseline_exact,
        "zero_baseline_fraction": zero_baseline_exact / len(rows),
        "supplied_exact": sum(value == 0 for value in supplied_abs),
        "supplied_exact_fraction": (
            sum(value == 0 for value in supplied_abs) / len(supplied)
            if supplied
            else None
        ),
        "within_one": sum(value <= 1 for value in absolute),
        "within_one_fraction": sum(value <= 1 for value in absolute) / len(rows),
        "greater_than_one": sum(value > 1 for value in absolute),
        "bias": mean_difference,
        "end_to_end_mae": statistics.fmean(absolute),
        "rmse": math.sqrt(statistics.fmean(value**2 for value in absolute)),
        "median_absolute_difference": statistics.median(absolute),
        "conditional_bias": (
            statistics.fmean(int(row["difference"]) for row in supplied)
            if supplied
            else None
        ),
        "conditional_mae": statistics.fmean(supplied_abs) if supplied else None,
        "difference_min": min(differences),
        "difference_max": max(differences),
        "empirical_limits": {
            "lower_2_5": _percentile(differences, 0.025),
            "upper_97_5": _percentile(differences, 0.975),
        },
        "bland_altman": {
            "mean": mean_difference,
            "sd": sd_difference,
            "lower": mean_difference - 1.96 * sd_difference,
            "upper": mean_difference + 1.96 * sd_difference,
            "note": "reported for completeness; the differences are zero-inflated",
        },
    }
    if report_metric:
        asserted = int(report_metric.get("gold_asserted") or 0)
        predicted = int(report_metric.get("predicted") or 0)
        if asserted != record["gold_rows"]:
            raise SystemExit("asserted gold rows do not reproduce the locked report")
        if predicted != record["supplied"]:
            raise SystemExit("supplied counts do not reproduce the locked report")
        reported_mae = report_metric.get("mae")
        mine = record["conditional_mae"]
        if (reported_mae is None) != (mine is None) or (
            reported_mae is not None and abs(float(reported_mae) - mine) > 1e-9
        ):
            raise SystemExit("conditional MAE does not reproduce the locked report")
    return record


def identity_summary(report: dict[str, Any]) -> dict[str, Any]:
    overall = report["overall"]
    tp, fp, fn = (int(overall[key]) for key in ("tp", "fp", "fn"))
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "recall": tp / (tp + fn) if tp + fn else None,
        "precision": tp / (tp + fp) if tp + fp else None,
        "primary_score_lane": report.get("primary_score_lane") or "primary",
    }


def load_arm(label: str, run_dir: Path, gold_root: Path | None) -> dict[str, Any]:
    """Everything the figure needs about one locked arm."""
    run_dir = run_dir.resolve()
    lock = verify_lock(run_dir)
    report = read_json(run_dir / "report.json")
    predictions = read_json(run_dir / "predictions.json")
    gold_sources = resolve_gold_sources(report, gold_root)
    rows = build_rows(label, report, predictions, gold_sources)
    count_metrics = (report.get("overall") or {}).get("count") or {}
    metrics = {
        field: summarize_field(
            [row for row in rows if row["measure"] == field],
            count_metrics.get(field),
        )
        for field in FIELDS
    }
    setup_path = run_dir / "setup.json"
    fingerprint = None
    if setup_path.is_file():
        setup = read_json(setup_path)
        fingerprint = ((setup.get("repository") or {}).get("runtime_source") or {}).get(
            "sha256"
        )
    consumption = report.get("tranche_consumption") or {}
    return {
        "label": label,
        "run_dir": run_dir,
        "run_id": str(report["run_id"]),
        "locked_at": str(lock.get("locked_at") or report.get("locked_at") or ""),
        "integrity": lock["verified_hashes"],
        "runtime_source_sha256": fingerprint,
        "tier_id": consumption.get("tier_id"),
        "comparison_arm": consumption.get("comparison_arm"),
        "attempts": len(report["papers"]),
        "unique_pmids": len({str(paper["pmid"]) for paper in report["papers"]}),
        "gold_sources": {
            gene: {"path": str(path), "sha256": sha256(path)}
            for gene, path in gold_sources.items()
        },
        "identity": identity_summary(report),
        "metrics": metrics,
        "rows": rows,
    }


def discover_baseline(run_dir: Path, report: dict[str, Any]) -> Path | None:
    """For a scored candidate arm, find its tranche's scored baseline arm."""
    consumption = report.get("tranche_consumption") or {}
    if consumption.get("comparison_arm") != "candidate":
        return None
    setup_path = run_dir / "setup.json"
    if not setup_path.is_file():
        return None
    cohort = read_json(setup_path).get("cohort") or {}
    registry = str(cohort.get("registry") or "").strip()
    relative = str((cohort.get("consumption_log") or {}).get("path") or "").strip()
    if not registry or not relative:
        return None
    log_path = Path(registry).resolve().parent / relative
    if not log_path.is_file():
        return None
    baselines = []
    for raw in log_path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        entry = json.loads(raw)
        if (
            entry.get("tier_id") == consumption.get("tier_id")
            and entry.get("comparison_arm") == "baseline"
        ):
            baselines.append(entry)
    if not baselines:
        return None
    latest = sorted(baselines, key=lambda entry: str(entry.get("scored_at") or ""))[-1]
    baseline_dir = run_dir.parent / str(latest.get("run_id") or "")
    return baseline_dir if (baseline_dir / "report.json").is_file() else None


def discover_comparison(run_dir: Path) -> Path | None:
    """The newest ``compare`` JSON in the run that names this run as candidate."""
    found: list[tuple[int, float, Path]] = []
    for path in sorted(run_dir.glob("compare*.json")):
        try:
            data = read_json(path)
        except (OSError, ValueError):
            continue
        candidate = ((data.get("integrity") or {}).get("candidate_report")) or ""
        if not candidate or Path(candidate).resolve().parent != run_dir.resolve():
            continue
        found.append(
            (
                1 if data.get("secondary_count_endpoint") else 0,
                path.stat().st_mtime,
                path,
            )
        )
    return sorted(found)[-1][2] if found else None


def load_comparison(
    path: Path, candidate: dict[str, Any], baseline: dict[str, Any] | None
) -> dict[str, Any]:
    """Registered paired decision for these two arms, checked against our rows."""
    data = read_json(path)
    integrity = data.get("integrity") or {}
    candidate_report = str(integrity.get("candidate_report") or "")
    if Path(candidate_report).resolve().parent != candidate["run_dir"]:
        raise SystemExit(f"{path} does not compare this run as the candidate arm")
    baseline_report = str(integrity.get("baseline_report") or "")
    if (
        baseline is not None
        and Path(baseline_report).resolve().parent != (baseline["run_dir"])
    ):
        raise SystemExit(f"{path} does not compare {baseline['run_dir']} as baseline")
    secondary = data.get("secondary_count_endpoint") or {}
    observed = secondary.get("observed") or {}
    if baseline is not None and observed:
        for arm in (baseline, candidate):
            mine = arm["metrics"].get(secondary.get("count_field") or "carriers") or {}
            theirs = observed.get(arm["label"]) or {}
            for key in ("end_to_end_mae", "coverage_on_matched"):
                if key in theirs and mine.get(key) is not None:
                    if abs(float(theirs[key]) - float(mine[key])) > 1e-9:
                        raise SystemExit(
                            f"{key} for the {arm['label']} arm does not reproduce "
                            f"{path.name}"
                        )
    return {
        "path": str(path),
        "sha256": sha256(path),
        "phase": data.get("phase"),
        "decision": data.get("decision"),
        "passed": data.get("passed"),
        "metrics": data.get("metrics") or {},
        "criteria": data.get("criteria") or {},
        "secondary_count_endpoint": (
            {
                "count_field": secondary.get("count_field"),
                "observed": observed,
                "deltas": secondary.get("deltas") or {},
                "criteria": secondary.get("criteria") or {},
                "passed": secondary.get("passed"),
            }
            if secondary
            else None
        ),
    }


def pair_arms(candidate: dict[str, Any], baseline: dict[str, Any]) -> dict[str, Any]:
    """Row-paired deltas: the same gold row under two protocols."""

    def key(row: dict[str, Any]) -> tuple[str, str, int, str]:
        return (row["gene"], row["pmid"], int(row["gold_row_index"]), row["measure"])

    candidate_rows = {key(row): row for row in candidate["rows"]}
    baseline_rows = {key(row): row for row in baseline["rows"]}
    if set(candidate_rows) != set(baseline_rows):
        raise SystemExit("paired arms do not evaluate the same gold rows")
    for item, row in candidate_rows.items():
        if row["gold_count"] != baseline_rows[item]["gold_count"]:
            raise SystemExit(f"gold count differs between arms for {item}")
    per_field: dict[str, Any] = {}
    for field in FIELDS:
        keys = [item for item in candidate_rows if item[3] == field]
        improved = worsened = unchanged = 0
        transitions: Counter[str] = Counter()
        for item in keys:
            before = abs(int(baseline_rows[item]["difference"]))
            after = abs(int(candidate_rows[item]["difference"]))
            if after < before:
                improved += 1
            elif after > before:
                worsened += 1
            else:
                unchanged += 1
            transitions[
                f"{baseline_rows[item]['status']}->{candidate_rows[item]['status']}"
            ] += 1
        deltas = {}
        for metric in (
            "end_to_end_mae",
            "bias",
            "rmse",
            "coverage_on_matched",
            "supplied_fraction",
            "exact_fraction",
            "identity_coverage",
            "supplied",
            "exact",
            "identity_matched",
        ):
            before_value = baseline["metrics"][field].get(metric)
            after_value = candidate["metrics"][field].get(metric)
            deltas[metric] = (
                after_value - before_value
                if before_value is not None and after_value is not None
                else None
            )
        per_field[field] = {
            "rows": len(keys),
            "rows_improved": improved,
            "rows_worsened": worsened,
            "rows_unchanged": unchanged,
            "status_transitions": dict(sorted(transitions.items())),
            "delta_candidate_minus_baseline": deltas,
        }
    identity = {}
    for metric in ("recall", "precision"):
        before_value = baseline["identity"].get(metric)
        after_value = candidate["identity"].get(metric)
        identity[metric] = {
            "baseline": before_value,
            "candidate": after_value,
            "delta_candidate_minus_baseline": (
                after_value - before_value
                if before_value is not None and after_value is not None
                else None
            ),
        }
    return {
        "definition": "candidate minus baseline on the same gold rows",
        "gold_rows": len(candidate_rows),
        "identity": identity,
        "by_measure": per_field,
    }


# ------------------------------------------------------------------------ drawing
def nice_ceiling(value: float) -> int:
    for candidate in (20, 50, 110, 220, 550, 1100, 2200):
        if value <= candidate * 0.97:
            return candidate
    return int(math.ceil(value * 1.1))


def x_coord(value: float, x0: float, xmax: float) -> float:
    usable = PANEL - 2 * INSET
    return x0 + INSET + math.sqrt(max(value, 0) / xmax) * usable


def y_coord(difference: float, top: float, ymax: float) -> float:
    half = (PANEL - 2 * INSET) / 2
    center = top + PANEL / 2
    if difference == 0:
        return center
    return center - math.copysign(math.sqrt(abs(difference) / ymax), difference) * half


def fmt_pct(value: float | None, digits: int = 1) -> str:
    return "n/a" if value is None else f"{100 * value:.{digits}f}%"


def fmt_num(value: float | None, digits: int = 3) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}".replace("-", "−")


def fmt_signed(value: float | None, digits: int = 3, suffix: str = "") -> str:
    if value is None:
        return "n/a"
    sign = "+" if value > 0 else ("−" if value < 0 else "±")
    return f"{sign}{abs(value):.{digits}f}{suffix}"


def fmt_int(value: int) -> str:
    return f"{value:+d}".replace("-", "−") if value else "0"


def fmt_pp(value: float | None) -> str:
    return "n/a" if value is None else fmt_signed(100 * value, 2, " pp")


def grouped_markers(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One marker per (reference, automated) coordinate; wedges split it by status."""
    groups: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(int(row["gold_count"]), int(row["automated_count_evaluated"]))].append(
            row
        )
    result = []
    for (gold, evaluated), members in groups.items():
        by_status = Counter(str(row["status"]) for row in members)
        result.append(
            {
                "gold": gold,
                "evaluated": evaluated,
                "difference": evaluated - gold,
                "n": len(members),
                "members": members,
                "by_status": {
                    status: by_status[status]
                    for status in STATUSES
                    if by_status[status]
                },
            }
        )
    result.sort(key=lambda group: (-group["n"], group["gold"], group["evaluated"]))
    return result


def _wedge_path(cx: float, cy: float, radius: float, start: float, end: float) -> str:
    x_start = cx + radius * math.cos(start)
    y_start = cy + radius * math.sin(start)
    x_end = cx + radius * math.cos(end)
    y_end = cy + radius * math.sin(end)
    large_arc = 1 if end - start > math.pi else 0
    return (
        f"M {cx:.1f} {cy:.1f} L {x_start:.1f} {y_start:.1f} "
        f"A {radius:.1f} {radius:.1f} 0 {large_arc} 1 {x_end:.1f} {y_end:.1f} Z"
    )


def draw_marker(
    group: dict[str, Any], x0: float, top: float, xmax: float, ymax: float
) -> str:
    cx = x_coord(group["gold"], x0, xmax)
    cy = y_coord(group["difference"], top, ymax)
    radius = bubble_radius(group["n"])
    names = "; ".join(
        f"{row['gene']} {row['gold_variant']} PMID {row['pmid']} ({row['status']})"
        for row in group["members"][:12]
    )
    if group["n"] > 12:
        names += f"; … {group['n'] - 12} more"
    shares = ", ".join(
        f"{count} {STATUS_LABELS[status]}"
        for status, count in group["by_status"].items()
    )
    title = escape(
        f"reference {group['gold']}, automated {group['evaluated']}, difference "
        f"{group['difference']:+d}; {group['n']} row(s): {shares}: {names}"
    )
    if len(group["by_status"]) == 1:
        status = next(iter(group["by_status"]))
        style = STATUS_STYLE[status]
        dash = f' stroke-dasharray="{style["dash"]}"' if style["dash"] else ""
        return (
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{radius:.1f}" '
            f'fill="{style["fill"]}" stroke="{style["stroke"]}" '
            f'stroke-width="{style["stroke_width"]}" fill-opacity="{style["opacity"]}"'
            f"{dash}><title>{title}</title></circle>"
        )
    # A shared coordinate (typically misses and abstentions on the floor) keeps one
    # marker whose area is the row count; wedges show the status shares.
    parts = [f"<g><title>{title}</title>"]
    angle = -math.pi / 2
    for status, count in group["by_status"].items():
        style = STATUS_STYLE[status]
        sweep = 2 * math.pi * count / group["n"]
        dash = f' stroke-dasharray="{style["dash"]}"' if style["dash"] else ""
        parts.append(
            f'<path d="{_wedge_path(cx, cy, radius, angle, angle + sweep)}" '
            f'fill="{style["fill"]}" stroke="{style["stroke"]}" stroke-width="2.2" '
            f'fill-opacity="{style["opacity"]}"{dash}/>'
        )
        angle += sweep
    parts.append("</g>")
    return "".join(parts)


def choose_labels(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Largest supplied disagreements plus the largest count lost to a miss."""
    supplied = sorted(
        (row for row in rows if row["status"] == "supplied" and row["difference"]),
        key=lambda row: (-abs(int(row["difference"])), row["gene"], row["pmid"]),
    )
    chosen = supplied[:3]
    lost = sorted(
        (row for row in rows if row["status"] != "supplied" and row["gold_count"]),
        key=lambda row: (-int(row["gold_count"]), row["gene"], row["pmid"]),
    )
    if lost:
        chosen.append(lost[0])
    return chosen


def _boxes_overlap(a: tuple[float, ...], b: tuple[float, ...]) -> bool:
    return not (
        a[0] + a[2] < b[0]
        or b[0] + b[2] < a[0]
        or a[1] + a[3] < b[1]
        or b[1] + b[3] < a[1]
    )


def _overlap_area(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    width = min(a[0] + a[2], b[0] + b[2]) - max(a[0], b[0])
    height = min(a[1] + a[3], b[1] + b[3]) - max(a[1], b[1])
    return width * height if width > 0 and height > 0 else 0.0


LABEL_OFFSETS = (0, -62, 62, -124, 124, -186, 186, -248, 248, -310, 310, -372, 372)


def draw_labels(
    rows: list[dict[str, Any]],
    groups: list[dict[str, Any]],
    x0: float,
    top: float,
    xmax: float,
    ymax: float,
    blocked: list[tuple[float, float, float, float]] | None = None,
    soft: list[tuple[float, float, float, float]] | None = None,
) -> list[str]:
    """Label the chosen rows beside their markers without overlapping each other.

    Each label tries both sides of its marker and a widening ladder of vertical
    offsets, takes the first spot that overlaps nothing already placed (other
    labels and the reference-line captions), and otherwise the least-overlapping
    spot, so a crowded panel degrades gracefully instead of oscillating.
    """
    size_of = {(group["gold"], group["evaluated"]): group["n"] for group in groups}
    placed: list[tuple[float, float, float, float]] = list(blocked or [])
    output: list[str] = []
    for row in choose_labels(rows):
        gold = int(row["gold_count"])
        evaluated = int(row["automated_count_evaluated"])
        mx = x_coord(gold, x0, xmax)
        my = y_coord(evaluated - gold, top, ymax)
        radius = bubble_radius(size_of.get((gold, evaluated), 1))
        title = f"{row['gene']} {row['gold_variant']}"
        if row["status"] == "supplied":
            detail = f"PMID {row['pmid']} · ref {gold} · auto {evaluated}"
        elif row["status"] == "abstained":
            detail = f"PMID {row['pmid']} · ref {gold} · no automated count"
        else:
            detail = f"PMID {row['pmid']} · ref {gold} · variant not found"
        width = max(len(title) * 0.56 * 26, len(detail) * 0.52 * 22)
        room_right = (x0 + PANEL - 12) - (mx + radius + 14)
        room_left = (mx - radius - 14) - (x0 + 12)
        sides = ("start", "end") if room_right >= room_left else ("end", "start")
        best: tuple[float, str, float, float, tuple[float, ...]] | None = None
        for side in sides:
            lx = mx + radius + 14 if side == "start" else mx - radius - 14
            box_x = lx if side == "start" else lx - width
            outside = max(0.0, x0 + 6 - box_x) + max(
                0.0, box_x + width - (x0 + PANEL - 6)
            )
            for offset in LABEL_OFFSETS:
                ly = min(max(my - 8 + offset, top + 44), top + PANEL - 40)
                box = (box_x, ly - 26, width, 56)
                cost = sum(_overlap_area(box, other) for other in placed)
                cost += 0.35 * sum(_overlap_area(box, other) for other in soft or ())
                cost += outside * 56 + abs(offset) * 0.01
                if best is None or cost < best[0]:
                    best = (cost, side, lx, ly, box)
                if cost < 1.0:
                    break
            if best is not None and best[0] < 1.0:
                break
        assert best is not None
        _cost, anchor, lx, ly, box = best
        placed.append(box)
        if abs((ly + 10) - my) > radius + 16:
            start_x = mx + radius if anchor == "start" else mx - radius
            end_x = lx - 4 if anchor == "start" else lx + 4
            output.append(
                f'<line x1="{start_x:.1f}" y1="{my:.1f}" x2="{end_x:.1f}" '
                f'y2="{ly + 10:.1f}" class="leader"/>'
            )
        output.append(text_element(lx, ly, title, "label-title", anchor=anchor))
        output.append(text_element(lx, ly + 26, detail, "label-detail", anchor=anchor))
    return output


def draw_panel(
    rows: list[dict[str, Any]],
    metrics: dict[str, Any],
    field: str,
    x0: float,
    top: float,
    xmax: float,
    ymax: float,
    *,
    clip_id: str,
    y_labels: bool,
) -> list[str]:
    lines: list[str] = [
        text_element(x0, top - 28, FIELD_TITLES[field], "facet-title"),
        (
            f'<defs><clipPath id="{clip_id}"><rect x="{x0}" y="{top}" '
            f'width="{PANEL}" height="{PANEL}"/></clipPath></defs>'
        ),
        (
            f'<rect x="{x0}" y="{top}" width="{PANEL}" height="{PANEL}" '
            f'fill="{COLORS["white"]}"/>'
        ),
    ]
    x_ticks = [tick for tick in TICK_VALUES if tick <= xmax]
    y_ticks = [tick for tick in TICK_VALUES if 0 < tick <= ymax]
    last_label_x: float | None = None
    for tick in x_ticks:
        x = x_coord(tick, x0, xmax)
        lines.append(
            f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + PANEL}" class="grid"/>'
        )
        if last_label_x is None or x - last_label_x >= 30:
            lines.append(
                text_element(x, top + PANEL + 37, str(tick), "tick", anchor="middle")
            )
            last_label_x = x
    labeled_y: list[float] = []
    signed_ticks = sorted(
        [-tick for tick in y_ticks] + [0] + y_ticks, key=lambda tick: (abs(tick), tick)
    )
    for tick in signed_ticks:
        y = y_coord(tick, top, ymax)
        lines.append(
            f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + PANEL}" y2="{y:.1f}" class="grid"/>'
        )
        if y_labels and all(abs(y - other) >= 28 for other in labeled_y):
            lines.append(
                text_element(x0 - 18, y + 10, fmt_int(tick), "tick", anchor="end")
            )
            labeled_y.append(y)
    # Lower boundary: a missing or abstained count sits exactly on difference = -ref.
    floor_points = " ".join(
        f"{x_coord(v, x0, xmax):.1f},{y_coord(-v, top, ymax):.1f}"
        for v in [xmax * (i / 80) ** 2 for i in range(81)]
    )
    lines.append(
        f'<polyline points="{floor_points}" class="floor-line" '
        f'clip-path="url(#{clip_id})"/>'
    )
    zero_y = y_coord(0, top, ymax)
    lines.append(
        f'<line x1="{x0}" y1="{zero_y:.1f}" x2="{x0 + PANEL}" y2="{zero_y:.1f}" '
        'class="zero-line"/>'
    )
    groups = grouped_markers(rows)
    marker_boxes: list[tuple[float, float, float, float]] = []
    for group in groups:
        cx = x_coord(group["gold"], x0, xmax)
        cy = y_coord(group["difference"], top, ymax)
        radius = bubble_radius(group["n"])
        marker_boxes.append((cx - radius, cy - radius, 2 * radius, 2 * radius))
    reference_lines = []
    if metrics.get("gold_rows"):
        limits = metrics["empirical_limits"]
        reference_lines = [
            (
                metrics["bias"],
                "bias-line",
                "ref-label-bias",
                f"bias {fmt_signed(metrics['bias'], 2)}",
            ),
            (
                limits["upper_97_5"],
                "limit-line",
                "ref-label",
                f"97.5%: {fmt_int(limits['upper_97_5'])}",
            ),
            (
                limits["lower_2_5"],
                "limit-line",
                "ref-label",
                f"2.5%: {fmt_int(limits['lower_2_5'])}",
            ),
        ]
    label_ys: list[float] = []
    blocked: list[tuple[float, float, float, float]] = []
    for value, class_name, label_class, text in reference_lines:
        if abs(value) > ymax:
            continue
        y = y_coord(value, top, ymax)
        lines.append(
            f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + PANEL}" y2="{y:.1f}" '
            f'class="{class_name}"/>'
        )
        label_y = y - 8
        while any(abs(label_y - other) < 26 for other in label_ys):
            label_y -= 26
        label_ys.append(label_y)
        width = len(text) * 0.55 * 22
        # Captions sit on whichever edge of the panel has fewer markers under it.
        right_box = (x0 + PANEL - 10 - width, label_y - 22, width, 28)
        left_box = (x0 + 10, label_y - 22, width, 28)
        right_cost = sum(_overlap_area(right_box, other) for other in marker_boxes)
        left_cost = sum(_overlap_area(left_box, other) for other in marker_boxes)
        if left_cost < right_cost:
            blocked.append(left_box)
            lines.append(text_element(x0 + 10, label_y, text, label_class))
        else:
            blocked.append(right_box)
            lines.append(
                text_element(x0 + PANEL - 10, label_y, text, label_class, anchor="end")
            )
    lines.append(f'<g clip-path="url(#{clip_id})">')
    lines.extend(draw_marker(group, x0, top, xmax, ymax) for group in groups)
    lines.append("</g>")
    lines.extend(draw_labels(rows, groups, x0, top, xmax, ymax, blocked, marker_boxes))
    lines.append(
        f'<rect x="{x0}" y="{top}" width="{PANEL}" height="{PANEL}" class="axis"/>'
    )
    lines.append(
        text_element(
            x0 + PANEL / 2,
            top + PANEL + 92,
            "Reference (gold) count",
            "axis-label",
            anchor="middle",
        )
    )
    if metrics.get("gold_rows"):
        stats_top = top + PANEL + 140
        limits = metrics["empirical_limits"]
        stats = [
            (
                f"{metrics['gold_rows']} gold rows · matched "
                f"{fmt_pct(metrics['identity_coverage'])} · supplied "
                f"{fmt_pct(metrics['supplied_fraction'])}"
            ),
            (
                f"exact {fmt_pct(metrics['exact_fraction'])} (zero-baseline "
                f"{fmt_pct(metrics['zero_baseline_fraction'])}) · within ±1 "
                f"{fmt_pct(metrics['within_one_fraction'])}"
            ),
            (
                f"bias {fmt_signed(metrics['bias'], 2)} · e2e MAE "
                f"{fmt_num(metrics['end_to_end_mae'])} · cond. MAE "
                f"{fmt_num(metrics['conditional_mae'])} (n {metrics['supplied']})"
            ),
            (
                f"empirical 95% limits {fmt_int(limits['lower_2_5'])} to "
                f"{fmt_int(limits['upper_97_5'])} · RMSE {fmt_num(metrics['rmse'], 2)}"
            ),
        ]
        for index, text in enumerate(stats):
            lines.append(text_element(x0, stats_top + 34 * index, text, "stats"))
    else:
        lines.append(
            text_element(x0, top + PANEL + 140, "no asserted gold rows", "stats")
        )
    return lines


def draw_legend(y: float) -> list[str]:
    lines: list[str] = []
    x = 80.0
    for status in ("supplied", "abstained", "identity_miss"):
        style = STATUS_STYLE[status]
        dash = f' stroke-dasharray="{style["dash"]}"' if style["dash"] else ""
        lines.append(
            f'<circle cx="{x + 14:.1f}" cy="{y - 9:.1f}" r="14" fill="{style["fill"]}" '
            f'stroke="{style["stroke"]}" stroke-width="{style["stroke_width"]}"{dash}/>'
        )
        label = STATUS_LABELS[status]
        lines.append(text_element(x + 38, y, label, "legend"))
        x += 38 + len(label) * 0.5 * 27 + 46
    key_text = "Marker area · gold rows at one coordinate (wedges split it by status):"
    key_y = y + 48
    lines.append(text_element(80, key_y, key_text, "legend"))
    cx = 80 + len(key_text) * 0.5 * 27 + 44
    for count in MARKER_SIZE_KEY_COUNTS:
        radius = bubble_radius(count)
        lines.append(
            f'<circle cx="{cx:.1f}" cy="{key_y - 9:.1f}" r="{radius:.1f}" class="size-key"/>'
        )
        lines.append(text_element(cx + radius + 8, key_y, str(count), "legend-value"))
        cx += radius * 2 + 82
    return lines


def arm_heading(arm: dict[str, Any], paired: bool) -> str:
    parts = []
    if paired:
        parts.append(f"{arm['label'].capitalize()} arm")
    else:
        parts.append("Run")
    parts.append(arm["run_id"])
    if arm.get("runtime_source_sha256"):
        parts.append(f"runtime {str(arm['runtime_source_sha256'])[:8]}")
    if arm.get("locked_at"):
        parts.append(f"locked {str(arm['locked_at'])[:10]}")
    parts.append(f"{arm['attempts']} gene–paper attempts")
    return " · ".join(parts)


def draw_delta_strip(
    top: float,
    candidate: dict[str, Any],
    baseline: dict[str, Any],
    paired: dict[str, Any],
    comparison: dict[str, Any] | None,
) -> list[str]:
    lines = [
        f'<line x1="80" y1="{top - 40:.1f}" x2="{WIDTH - 80}" y2="{top - 40:.1f}" '
        'class="divider"/>',
        text_element(
            80,
            top,
            (
                f"Candidate − baseline on the same {paired['gold_rows']} gold rows "
                "(paired by gold row; deltas are candidate minus baseline)"
            ),
            "arm-title",
        ),
    ]
    secondary = (comparison or {}).get("secondary_count_endpoint") or None
    for x0, field in zip(PANEL_X, FIELDS, strict=True):
        before = baseline["metrics"][field]
        after = candidate["metrics"][field]
        pair = paired["by_measure"][field]
        if not before.get("gold_rows"):
            lines.append(
                text_element(x0, top + 50, FIELD_TITLES[field], "stats-strong")
            )
            lines.append(text_element(x0, top + 84, "no asserted gold rows", "stats"))
            continue
        delta = pair["delta_candidate_minus_baseline"]
        mae_text = (
            f"e2e MAE {fmt_num(before['end_to_end_mae'])} → "
            f"{fmt_num(after['end_to_end_mae'])} ({fmt_signed(delta['end_to_end_mae'])})"
        )
        coverage_text = (
            f"coverage on matched {fmt_pct(before['coverage_on_matched'])} → "
            f"{fmt_pct(after['coverage_on_matched'])} "
            f"({fmt_pp(delta['coverage_on_matched'])})"
        )
        if secondary and secondary.get("count_field") == field:
            deltas = secondary.get("deltas") or {}
            mae_bound = (deltas.get("end_to_end_mae") or {}).get(
                "one_sided_upper_confidence_bound"
            )
            coverage_bound = (deltas.get("coverage_on_matched") or {}).get(
                "one_sided_lower_confidence_bound"
            )
            if mae_bound is not None:
                mae_text += f" · registered UB {fmt_signed(mae_bound)}"
            if coverage_bound is not None:
                coverage_text += f" · LB {fmt_pp(coverage_bound)}"
        texts = [
            mae_text,
            (
                f"supplied {before['supplied']} → {after['supplied']} · exact "
                f"{before['exact']} → {after['exact']} · bias "
                f"{fmt_signed(before['bias'], 2)} → {fmt_signed(after['bias'], 2)}"
            ),
            (
                f"rows better {pair['rows_improved']} · worse {pair['rows_worsened']} · "
                f"unchanged {pair['rows_unchanged']}"
            ),
            coverage_text,
        ]
        lines.append(text_element(x0, top + 50, FIELD_TITLES[field], "stats-strong"))
        for index, text in enumerate(texts):
            lines.append(text_element(x0, top + 84 + 32 * index, text, "stats"))
    identity = paired["identity"]
    identity_parts = []
    for metric in ("recall", "precision"):
        item = identity[metric]
        text = (
            f"{metric} {fmt_pct(item['baseline'], 2)} → {fmt_pct(item['candidate'], 2)} "
            f"({fmt_pp(item['delta_candidate_minus_baseline'])}"
        )
        bound = ((comparison or {}).get("metrics") or {}).get(metric) or {}
        if bound.get("one_sided_lower_confidence_bound") is not None:
            text += f"; LB {fmt_pp(bound['one_sided_lower_confidence_bound'])}"
        identity_parts.append(text + ")")
    if comparison:
        identity_parts.append(f"registered decision: {comparison.get('decision')}")
        if secondary:
            verdict = "passed" if secondary.get("passed") else "not passed"
            identity_parts.append(f"secondary count rule: {verdict}")
    else:
        identity_parts.append("no registered comparison found; bounds not shown")
    lines.append(
        text_element(
            80,
            top + 236,
            "Identity (primary endpoint): " + " · ".join(identity_parts),
            "stats-strong",
        )
    )
    return lines


def render_svg(
    arms: list[dict[str, Any]],
    *,
    paired: dict[str, Any] | None,
    comparison: dict[str, Any] | None,
    title: str,
) -> tuple[str, dict[str, Any]]:
    """Render every arm on shared axes; returns the SVG and the axis contract."""
    all_rows = [row for arm in arms for row in arm["rows"]]
    max_gold = max((int(row["gold_count"]) for row in all_rows), default=0)
    max_diff = max((abs(int(row["difference"])) for row in all_rows), default=0)
    xmax = nice_ceiling(max(max_gold, 10))
    ymax = nice_ceiling(max(max_diff, 10))
    is_paired = paired is not None and len(arms) == 2
    strip_height = 300 if is_paired else 0
    height = int(
        FIRST_ROW_TOP + ROW_PITCH * (len(arms) - 1) + PANEL + 300 + strip_height
    )
    primary = arms[-1]
    subtitle_parts = []
    if primary.get("tier_id"):
        subtitle_parts.append(f"tranche {primary['tier_id']}")
    if is_paired:
        subtitle_parts.append(
            "frozen baseline versus candidate protocol on the same manifest and gold"
        )
    subtitle_parts.append(
        "misses and abstentions scored as 0 · square-root axes · every asserted gold row"
    )
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH / 300:.2f}in" '
            f'height="{height / 300:.2f}in" viewBox="0 0 {WIDTH} {height}" role="img" '
            'aria-labelledby="figure-title figure-desc">'
        ),
        f'<title id="figure-title">{escape(title)}</title>',
        (
            '<desc id="figure-desc">Difference plot of automated minus reference '
            "counts against the reference count for every gold row with an asserted "
            "value. Identity misses and matched rows without an automated count are "
            "evaluated as zero and drawn in their own styles.</desc>"
        ),
        "<metadata>",
        *(
            f"Arm {escape(arm['label'])}: {escape(arm['run_id'])}, locked "
            f"{escape(str(arm['locked_at']))}"
            for arm in arms
        ),
        f"Evaluation rule: {EVALUATION_RULE}",
        f"Axes: reference 0..{xmax} (square root); difference ±{ymax} (signed square root)",
        "Marker radius is 8 + 2*sqrt(rows), capped at 28 viewBox units",
        "</metadata>",
        "<defs><style><![CDATA[",
        (
            'text { font-family: "Avenir Next Condensed", "Source Sans 3", '
            '"Helvetica Neue", Arial, sans-serif; fill: #17324D; }'
        ),
        ".title { font-size: 56px; font-weight: 600; letter-spacing: -0.3px; }",
        ".subtitle { font-size: 30px; font-weight: 500; fill: #687584; }",
        ".arm-title { font-size: 36px; font-weight: 600; }",
        ".facet-title { font-size: 40px; font-weight: 600; }",
        ".tick { font-size: 27px; font-weight: 400; fill: #687584; }",
        ".axis-label { font-size: 30px; font-weight: 500; }",
        ".stats { font-size: 25px; font-weight: 500; }",
        ".stats-strong { font-size: 27px; font-weight: 600; }",
        ".legend { font-size: 27px; font-weight: 500; }",
        ".legend-value { font-size: 23px; font-weight: 600; fill: #687584; }",
        ".size-key { fill: #FFFFFF; stroke: #17324D; stroke-width: 2.5; }",
        ".ref-label { font-size: 22px; font-weight: 600; fill: #687584; }",
        ".ref-label-bias { font-size: 22px; font-weight: 600; fill: #C85A63; }",
        ".label-title { font-size: 26px; font-weight: 600; }",
        ".label-detail { font-size: 22px; font-weight: 500; fill: #687584; }",
        ".leader { stroke: #687584; stroke-width: 2.2; fill: none; }",
        ".grid { stroke: #DCE3E8; stroke-width: 2; }",
        ".divider { stroke: #DCE3E8; stroke-width: 3; }",
        ".axis { stroke: #17324D; stroke-width: 3.5; fill: none; }",
        ".zero-line { stroke: #17324D; stroke-width: 3.5; }",
        ".bias-line { stroke: #C85A63; stroke-width: 3.5; stroke-dasharray: 14 10; }",
        ".limit-line { stroke: #687584; stroke-width: 2.6; stroke-dasharray: 5 7; }",
        ".floor-line { stroke: #9AA6B2; stroke-width: 2.4; stroke-dasharray: 3 9; fill: none; }",
        "]]></style></defs>",
        f'<rect width="{WIDTH}" height="{height}" fill="{COLORS["off_white"]}"/>',
        text_element(80, 82, title, "title"),
        text_element(80, 138, " · ".join(subtitle_parts), "subtitle"),
        *draw_legend(204),
    ]
    for arm_index, arm in enumerate(arms):
        row_top = FIRST_ROW_TOP + ROW_PITCH * arm_index
        lines.append(
            text_element(80, row_top - 82, arm_heading(arm, is_paired), "arm-title")
        )
        for x0, field in zip(PANEL_X, FIELDS, strict=True):
            lines.extend(
                draw_panel(
                    [row for row in arm["rows"] if row["measure"] == field],
                    arm["metrics"][field],
                    field,
                    x0,
                    row_top,
                    xmax,
                    ymax,
                    clip_id=f"clip-{arm_index}-{field}",
                    y_labels=(x0 == PANEL_X[0]),
                )
            )
        lines.append(
            f'<text x="60" y="{row_top + PANEL / 2:.1f}" class="axis-label" '
            f'text-anchor="middle" transform="rotate(-90 60 {row_top + PANEL / 2:.1f})">'
            "Automated − reference</text>"
        )
    if is_paired and paired is not None:
        strip_top = FIRST_ROW_TOP + ROW_PITCH * (len(arms) - 1) + PANEL + 320
        lines.extend(draw_delta_strip(strip_top, arms[1], arms[0], paired, comparison))
    lines.append("</svg>")
    contract = {
        "x_axis": f"reference (gold) count, square-root scale, 0 to {xmax}",
        "y_axis": f"automated minus reference, signed square-root scale, ±{ymax}",
        "x_max": xmax,
        "y_max": ymax,
        "width": WIDTH,
        "height": height,
    }
    return "\n".join(lines) + "\n", contract


# ------------------------------------------------------------------------ exports
def square_page(svg: str, width: int, height: int, side: int) -> str:
    """Return the SVG with a square viewBox/page so rasterizers keep every pixel."""
    root = svg.find("<svg")
    if root < 0:
        raise ValueError("not an SVG document")
    close = svg.index(">", root)
    tag = svg[root:close]
    tag = re.sub(r'viewBox="0 0 [0-9.]+ [0-9.]+"', f'viewBox="0 0 {side} {side}"', tag)
    tag = re.sub(r'\bwidth="[0-9.]+in"', f'width="{side / 300:.2f}in"', tag)
    tag = re.sub(r'\bheight="[0-9.]+in"', f'height="{side / 300:.2f}in"', tag)
    return svg[:root] + tag + svg[close:]


def export_svg_copies(
    svg_path: Path, png_path: Path, pdf_path: Path, *, width: int, height: int
) -> list[Path]:
    """Best-effort PNG/PDF copies without adding a Python dependency."""
    png_width = min(width, PNG_MAX_WIDTH)
    png_height = round(height * png_width / width)
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run(
            [
                rsvg,
                "--format=png",
                f"--width={png_width}",
                "--output",
                str(png_path),
                str(svg_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [rsvg, "--format=pdf", "--output", str(pdf_path), str(svg_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return [png_path, pdf_path]
    qlmanage = shutil.which("qlmanage")
    ffmpeg = shutil.which("ffmpeg")
    sips = shutil.which("sips")
    if qlmanage and ffmpeg and sips:
        # Quick Look renders into an s-by-s canvas anchored at the top-left and
        # scales the page width to s, so a portrait page would lose its bottom.
        # Rendering a square copy of the page (blank padding on the right or the
        # bottom) and cropping back to the viewBox keeps every pixel.
        page_side = max(width, height)
        raster_side = max(png_width, png_height)
        with tempfile.TemporaryDirectory(prefix="gvf-gold-difference-") as tmp:
            tmp_dir = Path(tmp)
            square = tmp_dir / svg_path.name
            square.write_text(
                square_page(
                    svg_path.read_text(encoding="utf-8"), width, height, page_side
                ),
                encoding="utf-8",
            )
            subprocess.run(
                [
                    qlmanage,
                    "-t",
                    "-s",
                    str(raster_side),
                    "-o",
                    str(tmp_dir),
                    str(square),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            rendered = tmp_dir / f"{square.name}.png"
            subprocess.run(
                [
                    ffmpeg,
                    "-y",
                    "-loglevel",
                    "error",
                    "-i",
                    str(rendered),
                    "-vf",
                    f"crop={png_width}:{png_height}:0:0",
                    "-frames:v",
                    "1",
                    str(png_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sips, "-s", "format", "pdf", str(svg_path), "--out", str(pdf_path)],
                check=True,
                capture_output=True,
                text=True,
            )
        return [png_path, pdf_path]
    if sips:
        with tempfile.TemporaryDirectory(prefix="gvf-gold-difference-") as tmp:
            intermediate = Path(tmp) / "figure.png"
            subprocess.run(
                [
                    sips,
                    "-s",
                    "format",
                    "png",
                    str(svg_path),
                    "--out",
                    str(intermediate),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [
                    sips,
                    "--resampleWidth",
                    str(png_width),
                    str(intermediate),
                    "--out",
                    str(png_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sips, "-s", "format", "pdf", str(svg_path), "--out", str(pdf_path)],
                check=True,
                capture_output=True,
                text=True,
            )
        return [png_path, pdf_path]
    print(
        "warning: no rsvg-convert or sips executable; retained editable SVG only",
        file=sys.stderr,
    )
    return []


CSV_COLUMNS = (
    "analysis_run",
    "arm",
    "gene",
    "pmid",
    "gold_row_index",
    "gold_variant",
    "predicted_variant",
    "measure",
    "gold_count",
    "automated_count_raw",
    "automated_count_evaluated",
    "status",
    "difference",
    "absolute_difference",
    "exact",
    "gold_value_set",
    "gold_v2_status",
    "gold_provenance",
    "prediction_source_layer",
    "prediction_source_location",
)


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CSV_COLUMNS))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    column: ("" if row.get(column) is None else row.get(column))
                    for column in CSV_COLUMNS
                }
            )


def arm_record(arm: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value for key, value in arm.items() if key not in {"rows", "run_dir"}
    } | {"run_dir": str(arm["run_dir"])}


def build_figure(
    run_dir: Path,
    *,
    gold_root: Path | None = None,
    baseline_run_dir: Path | None = None,
    comparison_path: Path | None = None,
    auto_baseline: bool = True,
    title: str | None = None,
) -> dict[str, Any]:
    """Load the arm(s), pair them, render, and return everything for the outputs."""
    candidate = load_arm("run", run_dir, gold_root)
    report = read_json(candidate["run_dir"] / "report.json")
    if baseline_run_dir is None and auto_baseline:
        baseline_run_dir = discover_baseline(candidate["run_dir"], report)
    baseline = None
    paired = None
    if baseline_run_dir is not None:
        candidate["label"] = "candidate"
        for row in candidate["rows"]:
            row["arm"] = "candidate"
        baseline = load_arm("baseline", baseline_run_dir, gold_root)
        paired = pair_arms(candidate, baseline)
    if comparison_path is None and baseline is not None:
        comparison_path = discover_comparison(candidate["run_dir"])
    comparison = (
        load_comparison(comparison_path, candidate, baseline)
        if comparison_path
        else None
    )
    arms = [baseline, candidate] if baseline is not None else [candidate]
    if title is None:
        title = (
            "Automated minus reference counts on the gold standard: protocol change"
            if baseline is not None
            else "Automated minus reference counts on the gold standard"
        )
    svg, contract = render_svg(arms, paired=paired, comparison=comparison, title=title)
    rows = [row for arm in arms for row in arm["rows"]]
    summary = {
        "figure": FIGURE_STEM,
        "title": title,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "evaluation_rule": {
            "text": EVALUATION_RULE,
            "identity_miss_evaluated_as": 0,
            "abstention_evaluated_as": 0,
            "stored_missing_value": None,
            "mirrors": (
                "benchmarks/codex_paper_eval/secondary_count_endpoint.py "
                "(end-to-end carrier error) extended to affected and unaffected"
            ),
        },
        "row_grain": "gold row with an asserted value, per count measure",
        "figure_contract": {
            **contract,
            "chart_family": (
                "faceted difference plots (automated minus reference against "
                "reference) with grouped markers sized by row count"
            ),
            "reference_lines": [
                "zero difference",
                "mean difference (bias)",
                "empirical 2.5th and 97.5th percentiles of the difference",
                "difference = -reference (the boundary where a missing count sits)",
            ],
            "marker_radius": "min(28, 8 + 2*sqrt(rows at the coordinate))",
            "labels": (
                "the three largest supplied disagreements and the largest reference "
                "count lost to a miss or abstention, placed automatically"
            ),
            "renderer": "editable SVG; vector PDF and PNG copies downstream",
        },
        "arms": {arm["label"]: arm_record(arm) for arm in arms},
        "paired": paired,
        "comparison": comparison,
    }
    return {"svg": svg, "rows": rows, "summary": summary, "contract": contract}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    parser.add_argument(
        "--gold-root",
        type=Path,
        default=None,
        help="override the gold CSV root; default is the root the score recorded",
    )
    parser.add_argument(
        "--baseline-run-dir",
        type=Path,
        default=None,
        help="paired frozen-baseline arm; discovered from the consumption log by default",
    )
    parser.add_argument(
        "--no-baseline",
        action="store_true",
        help="plot the run alone even when a paired baseline arm exists",
    )
    parser.add_argument(
        "--comparison-json",
        type=Path,
        default=None,
        help="run_eval.py compare output; the newest one in the run dir is used by default",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="output directory (default <run-dir>/figures); data goes to <out-dir>/data",
    )
    parser.add_argument("--stem", default=FIGURE_STEM)
    parser.add_argument("--title", default=None)
    parser.add_argument(
        "--copy-to",
        type=Path,
        default=None,
        help="also copy the outputs here, named by run id (for docs/figures)",
    )
    parser.add_argument("--svg-only", action="store_true", help="skip PNG/PDF export")
    args = parser.parse_args()

    result = build_figure(
        args.run_dir,
        gold_root=args.gold_root,
        baseline_run_dir=args.baseline_run_dir,
        comparison_path=args.comparison_json,
        auto_baseline=not args.no_baseline,
        title=args.title,
    )
    out_dir = (args.out_dir or (args.run_dir.resolve() / "figures")).resolve()
    data_dir = out_dir / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    svg_path = out_dir / f"{args.stem}.svg"
    png_path = out_dir / f"{args.stem}.png"
    pdf_path = out_dir / f"{args.stem}.pdf"
    csv_path = data_dir / f"{args.stem}.csv"
    json_path = data_dir / f"{args.stem}.json"
    svg_path.write_text(result["svg"], encoding="utf-8")
    write_csv(result["rows"], csv_path)
    exported: list[Path] = []
    if not args.svg_only:
        exported = export_svg_copies(
            svg_path,
            png_path,
            pdf_path,
            width=result["contract"]["width"],
            height=result["contract"]["height"],
        )
    outputs = {"svg": svg_path, "csv": csv_path, "json": json_path}
    for path in exported:
        outputs[path.suffix.lstrip(".")] = path
    result["summary"]["outputs"] = {key: str(path) for key, path in outputs.items()}
    json_path.write_text(
        json.dumps(result["summary"], indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(result['rows'])} plotted rows to {csv_path}")
    print(f"wrote summary to {json_path}")
    print(f"wrote editable figure to {svg_path}")
    for path in exported:
        print(f"wrote figure copy to {path}")
    if args.copy_to is not None:
        run_id = result["summary"]["arms"][
            "candidate" if "candidate" in result["summary"]["arms"] else "run"
        ]["run_id"]
        copy_dir = args.copy_to.resolve()
        (copy_dir / "data").mkdir(parents=True, exist_ok=True)
        for key, path in outputs.items():
            target = (copy_dir / "data" if key in {"csv", "json"} else copy_dir) / (
                f"{run_id}.{key}"
            )
            shutil.copyfile(path, target)
            print(f"copied {key} to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
