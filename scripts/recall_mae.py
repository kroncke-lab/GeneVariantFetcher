#!/usr/bin/env python3
"""Mean Absolute Error reporting for GVF recall metrics.

Two modes:

* ``rows`` — per-row absolute count error on matched gold variants.
  Consumes one or more per-gene ``discrepancies.csv`` files emitted by
  ``cli/compare_variants.py`` (via ``scripts/run_recall_suite.py``).
  For each gold↔SQLite match where both sides report a count, the script
  collects ``|excel - sqlite|`` for total carriers, affected, and
  unaffected. It also reports recall gap (``1 - recall``) for the
  set-membership metrics where MAE does not apply: PMIDs, variant rows,
  and unique variants.

* ``runs`` — absolute difference between two run summaries. Consumes two
  ``summary.json`` files produced by ``scripts/run_recall_suite.py`` and
  reports ``|after.recall - before.recall|`` for each metric, per gene
  and aggregated.

Outputs (under --outdir):

* ``mae_per_gene.csv`` — one row per gene per metric
* ``mae_summary.json`` — full structured report
* ``report.md`` — human-readable summary

Example::

    python scripts/recall_mae.py rows \
        --summary recall_metrics/post_insttoken_20260521/summary.json \
        --outdir recall_metrics/mae/baseline_20260521

    python scripts/recall_mae.py runs \
        --before recall_metrics/post_insttoken_20260521/summary.json \
        --after  recall_metrics/post_reextract_20260521/summary.json \
        --outdir recall_metrics/mae/delta_20260521
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

BASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_OUTROOT = BASE_DIR / "recall_metrics" / "mae"

COUNT_METRICS = ("carriers", "affected", "unaffected")
SET_METRICS = ("pmids", "variant_rows", "unique_variants")

EXCEL_KEY = {
    "carriers": "excel_carriers_total",
    "affected": "excel_affected",
    "unaffected": "excel_unaffected",
}
SQLITE_KEY = {
    "carriers": "sqlite_carriers_total",
    "affected": "sqlite_affected",
    "unaffected": "sqlite_unaffected",
}
DIFF_KEY = {
    "carriers": "carriers_diff",
    "affected": "affected_diff",
    "unaffected": "unaffected_diff",
}


@dataclass
class CountStats:
    """Distribution of absolute errors for one count field."""

    n: int = 0
    mae: float | None = None
    median_abs_err: float | None = None
    max_abs_err: float | None = None
    sum_abs_err: float = 0.0
    sum_excel: float = 0.0
    sum_sqlite: float = 0.0
    n_zero_error: int = 0


@dataclass
class GeneRowReport:
    """Per-gene per-row MAE results."""

    gene: str
    discrepancies_path: str | None
    n_matched_rows: int = 0
    counts: dict[str, CountStats] = field(default_factory=dict)
    set_metrics: dict[str, dict[str, float | int | None]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# rows mode
# ---------------------------------------------------------------------------


def _to_float(value: str | None) -> float | None:
    if value is None:
        return None
    s = value.strip()
    if not s:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if math.isnan(v):
        return None
    return v


def _is_matched_row(row: dict[str, str]) -> bool:
    """A row is a true gold↔SQLite match when neither side is missing."""
    missing_sqlite = str(row.get("missing_in_sqlite", "")).strip().lower() == "true"
    missing_excel = str(row.get("missing_in_excel", "")).strip().lower() == "true"
    return not (missing_sqlite or missing_excel)


def _gather_count_errors(
    rows: Iterable[dict[str, str]],
) -> tuple[int, dict[str, list[tuple[float, float, float]]]]:
    """Return (n_matched_rows, {metric: [(excel, sqlite, abs_err), ...]})."""
    samples: dict[str, list[tuple[float, float, float]]] = {
        m: [] for m in COUNT_METRICS
    }
    n_matched = 0
    for row in rows:
        if not _is_matched_row(row):
            continue
        n_matched += 1
        for metric in COUNT_METRICS:
            excel = _to_float(row.get(EXCEL_KEY[metric]))
            sqlite = _to_float(row.get(SQLITE_KEY[metric]))
            if excel is None or sqlite is None:
                continue
            samples[metric].append((excel, sqlite, abs(excel - sqlite)))
    return n_matched, samples


def _count_stats(samples: list[tuple[float, float, float]]) -> CountStats:
    if not samples:
        return CountStats()
    errs = [s[2] for s in samples]
    return CountStats(
        n=len(samples),
        mae=sum(errs) / len(errs),
        median_abs_err=statistics.median(errs),
        max_abs_err=max(errs),
        sum_abs_err=sum(errs),
        sum_excel=sum(s[0] for s in samples),
        sum_sqlite=sum(s[1] for s in samples),
        n_zero_error=sum(1 for e in errs if e == 0),
    )


def _set_metric_block(
    recall_entry: dict[str, Any] | None,
) -> dict[str, float | int | None]:
    if not recall_entry:
        return {
            "matched": None,
            "gold": None,
            "recall": None,
            "miss_count": None,
            "miss_rate": None,
        }
    matched = recall_entry.get("matched")
    gold = recall_entry.get("gold")
    recall = recall_entry.get("recall")
    try:
        matched_i = int(matched) if matched is not None else None
        gold_i = int(gold) if gold is not None else None
    except (TypeError, ValueError):
        matched_i, gold_i = None, None
    miss_count = (
        gold_i - matched_i if (matched_i is not None and gold_i is not None) else None
    )
    miss_rate = (1.0 - float(recall)) if isinstance(recall, (int, float)) else None
    return {
        "matched": matched_i,
        "gold": gold_i,
        "recall": recall,
        "miss_count": miss_count,
        "miss_rate": miss_rate,
    }


def _read_discrepancies(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _gene_report_from_outdir(
    gene: str,
    outdir: Path,
    recall_entry: dict[str, Any] | None,
) -> GeneRowReport:
    discrepancies = outdir / "discrepancies.csv"
    if not discrepancies.exists():
        return GeneRowReport(
            gene=gene,
            discrepancies_path=None,
            counts={m: CountStats() for m in COUNT_METRICS},
            set_metrics={m: _set_metric_block(None) for m in SET_METRICS},
        )

    n_matched, samples = _gather_count_errors(_read_discrepancies(discrepancies))
    counts = {m: _count_stats(samples[m]) for m in COUNT_METRICS}

    recall = (recall_entry or {}).get("recall", {}) if recall_entry else {}
    set_metrics = {m: _set_metric_block(recall.get(m)) for m in SET_METRICS}

    return GeneRowReport(
        gene=gene,
        discrepancies_path=str(discrepancies),
        n_matched_rows=n_matched,
        counts=counts,
        set_metrics=set_metrics,
    )


def _aggregate_counts(reports: list[GeneRowReport]) -> dict[str, CountStats]:
    """Pool absolute errors across genes."""
    aggregate: dict[str, CountStats] = {}
    for metric in COUNT_METRICS:
        n = sum(r.counts[metric].n for r in reports)
        sum_err = sum(r.counts[metric].sum_abs_err for r in reports)
        sum_excel = sum(r.counts[metric].sum_excel for r in reports)
        sum_sqlite = sum(r.counts[metric].sum_sqlite for r in reports)
        n_zero = sum(r.counts[metric].n_zero_error for r in reports)
        per_gene_maxes = [
            r.counts[metric].max_abs_err
            for r in reports
            if r.counts[metric].max_abs_err is not None
        ]
        aggregate[metric] = CountStats(
            n=n,
            mae=(sum_err / n) if n else None,
            # Pooled median is undefined without raw samples; leave None
            median_abs_err=None,
            max_abs_err=max(per_gene_maxes) if per_gene_maxes else None,
            sum_abs_err=sum_err,
            sum_excel=sum_excel,
            sum_sqlite=sum_sqlite,
            n_zero_error=n_zero,
        )
    return aggregate


def _aggregate_set_metrics(
    aggregate_recall: dict[str, dict[str, Any]] | None,
) -> dict[str, dict[str, float | int | None]]:
    if not aggregate_recall:
        return {m: _set_metric_block(None) for m in SET_METRICS}
    return {m: _set_metric_block(aggregate_recall.get(m)) for m in SET_METRICS}


def _count_stats_to_dict(stats: CountStats) -> dict[str, Any]:
    return {
        "n": stats.n,
        "mae": stats.mae,
        "median_abs_err": stats.median_abs_err,
        "max_abs_err": stats.max_abs_err,
        "sum_abs_err": stats.sum_abs_err,
        "sum_excel": stats.sum_excel,
        "sum_sqlite": stats.sum_sqlite,
        "n_zero_error": stats.n_zero_error,
        "fraction_zero_error": (stats.n_zero_error / stats.n) if stats.n else None,
        "mean_excel": (stats.sum_excel / stats.n) if stats.n else None,
        "mean_sqlite": (stats.sum_sqlite / stats.n) if stats.n else None,
    }


def run_rows_mode(summary_path: Path, outdir: Path) -> dict[str, Any]:
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_path}")

    summary = json.loads(summary_path.read_text())
    summary_dir = summary_path.parent

    gene_reports: list[GeneRowReport] = []
    for entry in summary.get("gene_results", []):
        if entry.get("status") != "scored":
            gene_reports.append(
                GeneRowReport(
                    gene=entry.get("gene", "?"),
                    discrepancies_path=None,
                    counts={m: CountStats() for m in COUNT_METRICS},
                    set_metrics={m: _set_metric_block(None) for m in SET_METRICS},
                )
            )
            continue

        gene = entry.get("gene", "?")
        # outdir from summary may be relative to repo root; also accept absolute.
        gene_outdir_raw = entry.get("outdir") or (
            (entry.get("summary") or {}).get("input_files", {}).get("excel")
        )
        candidates = []
        if gene_outdir_raw:
            p = Path(gene_outdir_raw)
            candidates.append(p if p.is_absolute() else BASE_DIR / p)
        candidates.append(summary_dir / gene)
        chosen = next(
            (c for c in candidates if (c / "discrepancies.csv").exists()), None
        )
        if chosen is None:
            chosen = candidates[-1]

        report = _gene_report_from_outdir(
            gene=gene,
            outdir=chosen,
            recall_entry=entry.get("summary"),
        )
        gene_reports.append(report)

    aggregate_count_stats = _aggregate_counts(gene_reports)
    aggregate_set_stats = _aggregate_set_metrics(summary.get("aggregate_recall"))

    result: dict[str, Any] = {
        "mode": "rows",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_summary": str(summary_path),
        "target_recall": summary.get("target_recall"),
        "genes": [],
        "aggregate": {
            "n_matched_rows": sum(r.n_matched_rows for r in gene_reports),
            "counts": {
                m: _count_stats_to_dict(aggregate_count_stats[m]) for m in COUNT_METRICS
            },
            "set_metrics": aggregate_set_stats,
        },
    }
    for r in gene_reports:
        result["genes"].append(
            {
                "gene": r.gene,
                "discrepancies_path": r.discrepancies_path,
                "n_matched_rows": r.n_matched_rows,
                "counts": {m: _count_stats_to_dict(r.counts[m]) for m in COUNT_METRICS},
                "set_metrics": r.set_metrics,
            }
        )

    outdir.mkdir(parents=True, exist_ok=True)
    _write_rows_outputs(result, outdir)
    return result


def _write_rows_outputs(result: dict[str, Any], outdir: Path) -> None:
    (outdir / "mae_summary.json").write_text(json.dumps(result, indent=2))

    csv_path = outdir / "mae_per_gene.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "gene",
                "metric_kind",
                "metric",
                "n",
                "mae",
                "median_abs_err",
                "max_abs_err",
                "sum_abs_err",
                "fraction_zero_error",
                "mean_excel",
                "mean_sqlite",
                "matched",
                "gold",
                "recall",
                "miss_count",
                "miss_rate",
            ]
        )
        for entry in result["genes"] + [
            {
                "gene": "AGGREGATE",
                **result["aggregate"],
                "n_matched_rows": result["aggregate"]["n_matched_rows"],
            }
        ]:
            gene = entry["gene"]
            counts = entry["counts"]
            set_metrics = entry["set_metrics"]
            for metric in COUNT_METRICS:
                c = counts[metric]
                writer.writerow(
                    [
                        gene,
                        "count",
                        metric,
                        c["n"],
                        c["mae"],
                        c["median_abs_err"],
                        c["max_abs_err"],
                        c["sum_abs_err"],
                        c["fraction_zero_error"],
                        c["mean_excel"],
                        c["mean_sqlite"],
                        "",
                        "",
                        "",
                        "",
                        "",
                    ]
                )
            for metric in SET_METRICS:
                s = set_metrics[metric]
                writer.writerow(
                    [
                        gene,
                        "set",
                        metric,
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        s["matched"],
                        s["gold"],
                        s["recall"],
                        s["miss_count"],
                        s["miss_rate"],
                    ]
                )

    (outdir / "report.md").write_text(_rows_markdown(result))


def _fmt_num(v: Any, places: int = 2) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{places}f}"
    return str(v)


def _fmt_pct(v: Any) -> str:
    if v is None:
        return "—"
    if isinstance(v, (int, float)):
        return f"{float(v) * 100:.1f}%"
    return str(v)


def _rows_markdown(result: dict[str, Any]) -> str:
    lines: list[str] = [
        "# Recall MAE Report (rows mode)",
        "",
        f"- Generated: `{result['generated_at_utc']}`",
        f"- Source summary: `{result['source_summary']}`",
        f"- Target recall: `{_fmt_pct(result.get('target_recall'))}`",
        f"- Matched rows across all genes: {result['aggregate']['n_matched_rows']}",
        "",
        "## Aggregate count-field MAE (matched gold↔SQLite pairs)",
        "",
        "| Metric | n pairs | MAE | Median |abs err| | Max |abs err| | Mean gold | Mean extracted | % zero-error |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for metric in COUNT_METRICS:
        c = result["aggregate"]["counts"][metric]
        lines.append(
            "| {metric} | {n} | {mae} | {med} | {mx} | {me} | {ms} | {fz} |".format(
                metric=metric,
                n=c["n"],
                mae=_fmt_num(c["mae"], 3),
                med=_fmt_num(c["median_abs_err"], 3),
                mx=_fmt_num(c["max_abs_err"], 2),
                me=_fmt_num(c["mean_excel"], 2),
                ms=_fmt_num(c["mean_sqlite"], 2),
                fz=_fmt_pct(c["fraction_zero_error"]),
            )
        )
    lines.append("")

    lines.extend(
        [
            "## Aggregate set-membership recall (gap to perfect)",
            "",
            "| Metric | Matched / Gold | Recall | Miss rate (1 - recall) | Miss count |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for metric in SET_METRICS:
        s = result["aggregate"]["set_metrics"][metric]
        matched = s["matched"]
        gold = s["gold"]
        denom = f"{matched}/{gold}" if matched is not None and gold is not None else "—"
        lines.append(
            f"| {metric} | {denom} | {_fmt_pct(s['recall'])} | {_fmt_pct(s['miss_rate'])} | {_fmt_num(s['miss_count'])} |"
        )
    lines.append("")

    lines.extend(["## Per-gene breakdown", ""])
    for entry in result["genes"]:
        lines.append(f"### {entry['gene']}")
        lines.append("")
        if entry["discrepancies_path"]:
            lines.append(f"- discrepancies: `{entry['discrepancies_path']}`")
        lines.append(f"- matched rows: {entry['n_matched_rows']}")
        lines.append("")
        lines.append("| Metric | n pairs | MAE | Median | Max |")
        lines.append("|---|---:|---:|---:|---:|")
        for metric in COUNT_METRICS:
            c = entry["counts"][metric]
            lines.append(
                f"| {metric} | {c['n']} | {_fmt_num(c['mae'], 3)} | "
                f"{_fmt_num(c['median_abs_err'], 3)} | {_fmt_num(c['max_abs_err'], 2)} |"
            )
        lines.append("")
        lines.append("| Set metric | Matched / Gold | Recall | Miss rate |")
        lines.append("|---|---|---:|---:|")
        for metric in SET_METRICS:
            s = entry["set_metrics"][metric]
            matched = s["matched"]
            gold = s["gold"]
            denom = (
                f"{matched}/{gold}" if matched is not None and gold is not None else "—"
            )
            lines.append(
                f"| {metric} | {denom} | {_fmt_pct(s['recall'])} | {_fmt_pct(s['miss_rate'])} |"
            )
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# runs mode
# ---------------------------------------------------------------------------


METRIC_LABELS = (
    "pmids",
    "variant_rows",
    "unique_variants",
    "patients",
    "affected",
    "unaffected",
)


def _safe_float(v: Any) -> float | None:
    if isinstance(v, (int, float)) and not (isinstance(v, float) and math.isnan(v)):
        return float(v)
    return None


def _metric_delta(
    before: dict[str, Any] | None,
    after: dict[str, Any] | None,
) -> dict[str, Any]:
    before = before or {}
    after = after or {}
    rb = _safe_float(before.get("recall"))
    ra = _safe_float(after.get("recall"))
    delta_recall = (ra - rb) if (rb is not None and ra is not None) else None
    abs_delta = abs(delta_recall) if delta_recall is not None else None
    return {
        "before_matched": before.get("matched"),
        "before_gold": before.get("gold"),
        "before_recall": rb,
        "after_matched": after.get("matched"),
        "after_gold": after.get("gold"),
        "after_recall": ra,
        "delta_recall": delta_recall,
        "abs_delta_recall": abs_delta,
    }


def _genes_index(summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for entry in summary.get("gene_results", []):
        gene = (entry.get("gene") or "").upper()
        if not gene:
            continue
        recall = (entry.get("summary") or {}).get("recall") or {}
        index[gene] = recall
    return index


def run_runs_mode(before_path: Path, after_path: Path, outdir: Path) -> dict[str, Any]:
    before = json.loads(before_path.read_text())
    after = json.loads(after_path.read_text())

    before_genes = _genes_index(before)
    after_genes = _genes_index(after)
    gene_keys = sorted(set(before_genes) | set(after_genes))

    gene_results = []
    for gene in gene_keys:
        b_recall = before_genes.get(gene, {})
        a_recall = after_genes.get(gene, {})
        per_metric = {
            metric: _metric_delta(b_recall.get(metric), a_recall.get(metric))
            for metric in METRIC_LABELS
        }
        deltas = [
            d["abs_delta_recall"]
            for d in per_metric.values()
            if d["abs_delta_recall"] is not None
        ]
        mae_recall = (sum(deltas) / len(deltas)) if deltas else None
        gene_results.append(
            {
                "gene": gene,
                "metrics": per_metric,
                "mae_recall_across_metrics": mae_recall,
            }
        )

    agg_before = before.get("aggregate_recall") or {}
    agg_after = after.get("aggregate_recall") or {}
    agg_metrics = {
        metric: _metric_delta(agg_before.get(metric), agg_after.get(metric))
        for metric in METRIC_LABELS
    }
    agg_deltas = [
        d["abs_delta_recall"]
        for d in agg_metrics.values()
        if d["abs_delta_recall"] is not None
    ]
    agg_mae = (sum(agg_deltas) / len(agg_deltas)) if agg_deltas else None

    result = {
        "mode": "runs",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "before_summary": str(before_path),
        "after_summary": str(after_path),
        "target_recall": after.get("target_recall") or before.get("target_recall"),
        "aggregate": {"metrics": agg_metrics, "mae_recall_across_metrics": agg_mae},
        "genes": gene_results,
    }

    outdir.mkdir(parents=True, exist_ok=True)
    _write_runs_outputs(result, outdir)
    return result


def _write_runs_outputs(result: dict[str, Any], outdir: Path) -> None:
    (outdir / "mae_summary.json").write_text(json.dumps(result, indent=2))

    csv_path = outdir / "mae_per_gene.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "gene",
                "metric",
                "before_matched",
                "before_gold",
                "before_recall",
                "after_matched",
                "after_gold",
                "after_recall",
                "delta_recall",
                "abs_delta_recall",
            ]
        )
        for entry in result["genes"]:
            for metric, d in entry["metrics"].items():
                writer.writerow(
                    [
                        entry["gene"],
                        metric,
                        d["before_matched"],
                        d["before_gold"],
                        d["before_recall"],
                        d["after_matched"],
                        d["after_gold"],
                        d["after_recall"],
                        d["delta_recall"],
                        d["abs_delta_recall"],
                    ]
                )
        for metric, d in result["aggregate"]["metrics"].items():
            writer.writerow(
                [
                    "AGGREGATE",
                    metric,
                    d["before_matched"],
                    d["before_gold"],
                    d["before_recall"],
                    d["after_matched"],
                    d["after_gold"],
                    d["after_recall"],
                    d["delta_recall"],
                    d["abs_delta_recall"],
                ]
            )

    (outdir / "report.md").write_text(_runs_markdown(result))


def _runs_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# Recall MAE Report (runs mode)",
        "",
        f"- Generated: `{result['generated_at_utc']}`",
        f"- Before: `{result['before_summary']}`",
        f"- After: `{result['after_summary']}`",
        f"- Target recall: `{_fmt_pct(result.get('target_recall'))}`",
        f"- Aggregate MAE of |Δrecall| across metrics: "
        f"`{_fmt_num(result['aggregate']['mae_recall_across_metrics'], 4)}` "
        f"({_fmt_pct(result['aggregate']['mae_recall_across_metrics'])})",
        "",
        "## Aggregate per-metric Δrecall",
        "",
        "| Metric | Before | After | Δrecall | |Δrecall| |",
        "|---|---:|---:|---:|---:|",
    ]
    for metric in METRIC_LABELS:
        d = result["aggregate"]["metrics"][metric]
        lines.append(
            f"| {metric} | {_fmt_pct(d['before_recall'])} | {_fmt_pct(d['after_recall'])} | "
            f"{_fmt_pct(d['delta_recall'])} | {_fmt_pct(d['abs_delta_recall'])} |"
        )
    lines.extend(["", "## Per-gene Δrecall", ""])
    for entry in result["genes"]:
        lines.append(f"### {entry['gene']}")
        lines.append("")
        lines.append(
            f"- MAE of |Δrecall| across metrics: "
            f"`{_fmt_pct(entry['mae_recall_across_metrics'])}`"
        )
        lines.append("")
        lines.append("| Metric | Before | After | Δrecall |")
        lines.append("|---|---:|---:|---:|")
        for metric, d in entry["metrics"].items():
            lines.append(
                f"| {metric} | {_fmt_pct(d['before_recall'])} | {_fmt_pct(d['after_recall'])} | "
                f"{_fmt_pct(d['delta_recall'])} |"
            )
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI entry
# ---------------------------------------------------------------------------


def _default_outdir(mode: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return DEFAULT_OUTROOT / f"{mode}_{stamp}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="GVF recall MAE reporter")
    sub = parser.add_subparsers(dest="mode", required=True)

    p_rows = sub.add_parser(
        "rows", help="Per-row MAE on count fields from discrepancies.csv"
    )
    p_rows.add_argument(
        "--summary",
        type=Path,
        required=True,
        help="Path to recall-suite summary.json (rows mode looks in adjacent gene subdirs)",
    )
    p_rows.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: recall_metrics/mae/rows_<timestamp>)",
    )

    p_runs = sub.add_parser(
        "runs", help="Run-vs-run |Δrecall| across two summary.json files"
    )
    p_runs.add_argument(
        "--before", type=Path, required=True, help="Baseline summary.json"
    )
    p_runs.add_argument(
        "--after", type=Path, required=True, help="Comparison summary.json"
    )
    p_runs.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: recall_metrics/mae/runs_<timestamp>)",
    )

    args = parser.parse_args(argv)
    outdir = args.outdir or _default_outdir(args.mode)
    outdir = outdir.expanduser()

    if args.mode == "rows":
        result = run_rows_mode(args.summary.expanduser(), outdir)
    elif args.mode == "runs":
        result = run_runs_mode(
            args.before.expanduser(), args.after.expanduser(), outdir
        )
    else:  # pragma: no cover — argparse enforces required
        parser.error("mode is required")

    # Tiny stdout summary
    if args.mode == "rows":
        agg = result["aggregate"]
        print(
            f"rows MAE — matched rows: {agg['n_matched_rows']}; "
            f"carriers MAE: {agg['counts']['carriers']['mae']}, "
            f"affected MAE: {agg['counts']['affected']['mae']}, "
            f"unaffected MAE: {agg['counts']['unaffected']['mae']}"
        )
    else:
        print(
            f"runs MAE — aggregate |Δrecall|: "
            f"{result['aggregate']['mae_recall_across_metrics']}"
        )
    print(f"Wrote outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
