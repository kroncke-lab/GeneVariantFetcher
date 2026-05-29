#!/usr/bin/env python3
"""Run extraction-vs-gold recall scoring across all available GVF genes."""

from __future__ import annotations

import argparse
import json
import logging
import math
import sqlite3
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parents[1]
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from cli.compare_variants import (  # noqa: E402
    aggregate_excel_data,
    aggregate_sqlite_data,
    compare_data,
    extract_sqlite_data,
    generate_outputs,
    introspect_sqlite,
    load_excel_data,
)
from scripts.recall_audit.common import write_csv_rows  # noqa: E402
from scripts.recall_audit.paper_disagreement_report import (  # noqa: E402
    REPORT_FIELDS,
    build_paper_disagreement_rows,
    select_regression_cases,
    write_markdown as write_regression_cases_markdown,
)

DEFAULT_GOLD_DIR = BASE_DIR / "gene_variant_fetcher_gold_standard"
DEFAULT_RESULTS_DIR = BASE_DIR / "results"
METRICS_DIR = BASE_DIR / "recall_metrics"
HISTORY_FILE = METRICS_DIR / "history.jsonl"
TARGET_RECALL = 0.90

RECALL_METRICS = (
    "pmids",
    "variant_rows",
    "unique_variants",
    "patients",
    "affected",
    "unaffected",
)

logger = logging.getLogger(__name__)


def discover_gold_inputs(gold_dir: Path) -> dict[str, Path]:
    """Discover normalized per-gene recall CSVs."""
    normalized_dir = gold_dir / "normalized"
    if not normalized_dir.exists():
        return {}
    inputs: dict[str, Path] = {}
    for path in sorted(normalized_dir.glob("*_recall_input.csv")):
        gene = path.name.removesuffix("_recall_input.csv").upper()
        inputs[gene] = path
    return inputs


def load_manifest(gold_dir: Path) -> dict[str, Any]:
    """Load optional gold-standard package manifest."""
    manifest_path = gold_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    with manifest_path.open() as f:
        return json.load(f)


def parse_gene_list(value: str | None) -> list[str] | None:
    if not value:
        return None
    genes = [g.strip().upper() for g in value.split(",") if g.strip()]
    return genes or None


def parse_db_overrides(items: list[str] | None) -> dict[str, Path]:
    """Parse repeated --db GENE=/path/to/db.sqlite arguments."""
    overrides: dict[str, Path] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"--db must be GENE=/path/to/db, got {item!r}")
        gene, raw_path = item.split("=", 1)
        gene = gene.strip().upper()
        if not gene:
            raise ValueError(f"--db must include a gene name, got {item!r}")
        overrides[gene] = Path(raw_path).expanduser()
    return overrides


def _db_sort_key(path: Path) -> tuple[float, int, str]:
    """Prefer newest DB, and within a run prefer the highest vN suffix."""
    mtime = path.stat().st_mtime
    match = path.stem.rsplit("_v", 1)
    version = 0
    if len(match) == 2 and match[1].isdigit():
        version = int(match[1])
    return (mtime, version, str(path))


def find_latest_db(gene: str, results_dir: Path) -> Path | None:
    """Find the newest extraction SQLite database for a gene."""
    gene_dir = results_dir / gene
    if not gene_dir.exists():
        return None
    candidates = [p for p in gene_dir.glob(f"**/{gene}*.db") if p.is_file()]
    if not candidates:
        return None
    return max(candidates, key=_db_sort_key)


def run_gene_compare(
    *,
    gene: str,
    gold_path: Path,
    db_path: Path,
    outdir: Path,
    variant_match_mode: str = "fuzzy",
    fuzzy_threshold: float = 0.80,
) -> dict[str, Any]:
    """Compare one gene's normalized gold input against one extraction DB."""
    gene_outdir = outdir / gene
    excel_df, detected_columns = load_excel_data(gold_path, None, None)

    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        table_info = introspect_sqlite(conn)
        sqlite_df = extract_sqlite_data(conn, table_info)
    finally:
        conn.close()

    gold_aggregated = aggregate_excel_data(excel_df, detected_columns)
    sqlite_aggregated = aggregate_sqlite_data(sqlite_df)
    rows = compare_data(
        gold_aggregated,
        sqlite_aggregated,
        variant_match_mode,
        fuzzy_threshold,
    )
    summary = generate_outputs(rows, gene_outdir, gold_path, db_path)
    return {
        "gene": gene,
        "status": "scored",
        "gold_path": str(gold_path),
        "db_path": str(db_path),
        "outdir": str(gene_outdir),
        "summary": summary,
    }


def combine_recall(scored: list[dict[str, Any]]) -> dict[str, Any]:
    """Sum matched/gold counts over genes for each recall dimension."""
    aggregate: dict[str, Any] = {}
    for metric_name in RECALL_METRICS:
        matched = 0
        gold = 0
        for item in scored:
            metric = item["summary"].get("recall", {}).get(metric_name, {})
            matched += int(metric.get("matched") or 0)
            gold += int(metric.get("gold") or 0)
        aggregate[metric_name] = {
            "matched": matched,
            "gold": gold,
            "recall": matched / gold if gold else None,
        }
    return aggregate


MAE_FIELDS = ("carriers", "affected", "unaffected")


def combine_mae(scored: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate rows-mode MAE across genes by summing |err| and N."""
    aggregate: dict[str, Any] = {}
    for field_name in MAE_FIELDS:
        total = 0
        n = 0
        for item in scored:
            mae_block = item.get("summary", {}).get("mae", {}).get(field_name, {})
            total += int(mae_block.get("sum_abs_error") or 0)
            n += int(mae_block.get("n_matched") or 0)
        aggregate[field_name] = {
            "sum_abs_error": total,
            "n_matched": n,
            "mae": (total / n) if n else None,
        }
    return aggregate


def combine_precision(scored: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate the gold-PMID-restricted precision metric across genes.

    Sums the per-gene numerator (matched DB rows) and denominator components
    (matched DB rows + extra DB rows on gold PMIDs), then recomputes the
    ratio. This is an extra-rows-relative-to-gold rate / false-positive upper
    bound, NOT clean precision; see ``compute_precision_summary`` for the
    caveat. Genes whose summaries predate this metric contribute nothing.
    """
    matched_db = 0
    extra_on_gold_pmids = 0
    for item in scored:
        block = item.get("summary", {}).get("precision", {})
        matched_db += int(block.get("matched_db") or 0)
        extra_on_gold_pmids += int(block.get("extra_on_gold_pmids") or 0)
    denominator = matched_db + extra_on_gold_pmids
    return {
        "matched_db": matched_db,
        "extra_on_gold_pmids": extra_on_gold_pmids,
        "precision_vs_gold_pmids": (matched_db / denominator) if denominator else None,
        "note": (
            "Upper bound on false-positive rate, restricted to gold-curated "
            "PMIDs. Gold is a curator-selected subset, not paper-exhaustive, "
            "so 'extra' DB rows on gold PMIDs may still be true positives the "
            "curator omitted. This is an extra-rows-relative-to-gold rate, NOT "
            "clean precision."
        ),
    }


def target_gap(metric: dict[str, Any], target: float) -> int | None:
    """Rows needed to reach target recall for one metric."""
    gold = int(metric.get("gold") or 0)
    matched = int(metric.get("matched") or 0)
    if gold == 0:
        return None
    needed = math.ceil(gold * target)
    return max(0, needed - matched)


def build_run_summary(
    *,
    gene_results: list[dict[str, Any]],
    aggregate_recall: dict[str, Any],
    aggregate_mae: dict[str, Any] | None = None,
    aggregate_precision: dict[str, Any] | None = None,
    manifest: dict[str, Any],
    target: float,
) -> dict[str, Any]:
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "target_recall": target,
        "aggregate_recall": aggregate_recall,
        "aggregate_target_gaps": {
            name: target_gap(metric, target)
            for name, metric in aggregate_recall.items()
        },
        "aggregate_mae": aggregate_mae or {},
        "aggregate_precision": aggregate_precision or {},
        "gene_results": gene_results,
        "gold_manifest_generated_at_utc": manifest.get("generated_at_utc"),
    }


def write_run_summary(summary: dict[str, Any], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
    write_markdown_summary(summary, outdir / "report.md")


def write_disagreement_artifacts(
    *,
    outdir: Path,
    gene_results: list[dict[str, Any]],
    results_dir: Path,
    gold_dir: Path,
) -> dict[str, Any]:
    """Write paper-level disagreement artifacts used to choose replay fixtures."""
    db_paths = {
        str(item["gene"]).upper(): Path(item["db_path"])
        for item in gene_results
        if item.get("status") == "scored" and item.get("db_path")
    }
    rows = build_paper_disagreement_rows(
        recall_score=outdir,
        db_paths_by_gene=db_paths,
        gold_dir=gold_dir / "normalized",
        results_dir=results_dir,
    )
    report_path = outdir / "paper_disagreement_report.csv"
    selected_path = outdir / "extraction_regression_cases.csv"
    markdown_path = outdir / "extraction_regression_cases.md"
    selected = select_regression_cases(rows)

    write_csv_rows(rows, REPORT_FIELDS, report_path)
    write_csv_rows(selected, REPORT_FIELDS, selected_path)
    write_regression_cases_markdown(
        selected,
        markdown_path,
        title="Available-Source Extraction Regression Cases",
    )
    return {
        "paper_disagreement_report": str(report_path),
        "extraction_regression_cases": str(selected_path),
        "extraction_regression_cases_markdown": str(markdown_path),
        "paper_disagreement_rows": len(rows),
        "selected_regression_cases": len(selected),
    }


def _format_pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def write_markdown_summary(summary: dict[str, Any], output_path: Path) -> None:
    lines = [
        "# Multi-Gene Recall Report",
        "",
        f"- Generated: `{summary['generated_at_utc']}`",
        f"- Target recall: `{summary['target_recall']:.0%}`",
        "",
        "## Aggregate Recall",
        "",
        "| Metric | Matched / Gold | Recall | Gap to Target |",
        "|--------|----------------|--------|---------------|",
    ]
    gaps = summary.get("aggregate_target_gaps", {})
    for name, metric in summary.get("aggregate_recall", {}).items():
        lines.append(
            f"| {name} | {metric.get('matched', 0)}/{metric.get('gold', 0)} | "
            f"{_format_pct(metric.get('recall'))} | {gaps.get(name, 'n/a')} |"
        )

    mae = summary.get("aggregate_mae") or {}
    if any(block.get("n_matched") for block in mae.values()):
        lines.extend(
            [
                "",
                "## Aggregate Rows-Mode MAE",
                "",
                "Mean absolute error on matched-variant rows (lower is better; target → 0).",
                "",
                "| Field | sum \\|err\\| / N matched | MAE |",
                "|-------|---------|-----|",
            ]
        )
        for field_name in MAE_FIELDS:
            block = mae.get(field_name, {})
            mae_val = block.get("mae")
            lines.append(
                f"| {field_name} | {block.get('sum_abs_error', 0)} / {block.get('n_matched', 0)} | "
                f"{'n/a' if mae_val is None else f'{mae_val:.3f}'} |"
            )

    precision = summary.get("aggregate_precision") or {}
    if precision.get("matched_db") or precision.get("extra_on_gold_pmids"):
        pvg = precision.get("precision_vs_gold_pmids")
        lines.extend(
            [
                "",
                "## Aggregate Precision (vs gold PMIDs)",
                "",
                "Extra-rows-relative-to-gold rate, restricted to gold-curated "
                "PMIDs. Upper bound on false positives, NOT clean precision "
                "(gold is a curator-selected subset, not paper-exhaustive).",
                "",
                "| Matched DB rows | Extra DB rows on gold PMIDs | precision_vs_gold_pmids |",
                "|-----------------|-----------------------------|-------------------------|",
                f"| {precision.get('matched_db', 0)} | "
                f"{precision.get('extra_on_gold_pmids', 0)} | {_format_pct(pvg)} |",
            ]
        )

    lines.extend(["", "## Genes", ""])
    for item in summary.get("gene_results", []):
        status = item.get("status")
        gene = item.get("gene")
        if status != "scored":
            note = item.get("note") or item.get("error") or ""
            lines.append(f"- `{gene}`: `{status}` - {note}")
            continue
        recall = item.get("summary", {}).get("recall", {})
        parts = []
        for metric_name in RECALL_METRICS:
            metric = recall.get(metric_name, {})
            parts.append(f"{metric_name} {_format_pct(metric.get('recall'))}")
        lines.append(
            f"- `{gene}`: scored from `{item.get('db_path')}`; " + ", ".join(parts)
        )
    lines.append("")

    artifacts = summary.get("disagreement_artifacts")
    if artifacts:
        lines.extend(
            [
                "## Paper-Level Disagreement",
                "",
                f"- Full report: `{artifacts.get('paper_disagreement_report')}`",
                f"- Regression cases: `{artifacts.get('extraction_regression_cases')}`",
                f"- Markdown review table: `{artifacts.get('extraction_regression_cases_markdown')}`",
                f"- Selected cases: `{artifacts.get('selected_regression_cases')}`",
                "",
            ]
        )
    elif summary.get("disagreement_artifacts_error"):
        lines.extend(
            [
                "## Paper-Level Disagreement",
                "",
                f"- Report generation failed: `{summary['disagreement_artifacts_error']}`",
                "",
            ]
        )

    with output_path.open("w") as f:
        f.write("\n".join(lines))


def print_scorecard(summary: dict[str, Any]) -> None:
    print()
    print("GVF MULTI-GENE RECALL")
    print(f"Generated: {summary['generated_at_utc']}")
    print(f"Target: {summary['target_recall']:.0%}")
    print()
    print("Aggregate")
    for name, metric in summary.get("aggregate_recall", {}).items():
        gap = summary.get("aggregate_target_gaps", {}).get(name)
        gap_text = "n/a" if gap is None else str(gap)
        print(
            f"  {name:16s} "
            f"{metric.get('matched', 0):>5}/{metric.get('gold', 0):<5} "
            f"{_format_pct(metric.get('recall')):>7} "
            f"gap_to_target={gap_text}"
        )

    mae = summary.get("aggregate_mae") or {}
    if any(block.get("n_matched") for block in mae.values()):
        print()
        print("MAE (rows-mode, matched only; lower is better)")
        for field_name in MAE_FIELDS:
            block = mae.get(field_name, {})
            mae_val = block.get("mae")
            mae_text = "n/a" if mae_val is None else f"{mae_val:.3f}"
            print(
                f"  {field_name:16s} "
                f"sum|err|={block.get('sum_abs_error', 0):>5}  "
                f"n={block.get('n_matched', 0):<5} "
                f"MAE={mae_text}"
            )

    precision = summary.get("aggregate_precision") or {}
    if precision.get("matched_db") or precision.get("extra_on_gold_pmids"):
        print()
        print("Precision (vs gold PMIDs; upper bound on FPs, not clean precision)")
        print(
            f"  matched_db={precision.get('matched_db', 0):>5}  "
            f"extra_on_gold_pmids={precision.get('extra_on_gold_pmids', 0):<5} "
            f"precision_vs_gold_pmids="
            f"{_format_pct(precision.get('precision_vs_gold_pmids'))}"
        )
    print()
    print("Genes")
    for item in summary.get("gene_results", []):
        gene = item.get("gene")
        status = item.get("status")
        if status != "scored":
            print(f"  {gene:8s} {status}: {item.get('note') or item.get('error')}")
            continue
        recall = item.get("summary", {}).get("recall", {})
        row_recall = recall.get("variant_rows", {}).get("recall")
        unique_recall = recall.get("unique_variants", {}).get("recall")
        pmid_recall = recall.get("pmids", {}).get("recall")
        print(
            f"  {gene:8s} scored "
            f"pmids={_format_pct(pmid_recall)} "
            f"variant_rows={_format_pct(row_recall)} "
            f"unique_variants={_format_pct(unique_recall)}"
        )
    print()


def log_metrics(summary: dict[str, Any], notes: str = "") -> None:
    METRICS_DIR.mkdir(exist_ok=True)
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "notes": notes,
        "summary": summary,
    }
    with HISTORY_FILE.open("a") as f:
        f.write(json.dumps(entry) + "\n")
    print(f"Metrics logged to {HISTORY_FILE}")


def run_tests() -> int:
    """Run recall pytest suite; fixtures skip when external artifacts are absent."""
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/recall/", "-v", "--tb=short"],
        cwd=str(BASE_DIR),
        check=False,
    )
    return result.returncode


def selected_genes(
    requested: list[str] | None,
    gold_inputs: dict[str, Path],
    manifest: dict[str, Any],
) -> list[str]:
    if requested:
        return requested
    genes = set(gold_inputs)
    genes.update((manifest.get("genes") or {}).keys())
    return sorted(g.upper() for g in genes)


def build_gene_results(
    *,
    genes: list[str],
    gold_inputs: dict[str, Path],
    db_overrides: dict[str, Path],
    results_dir: Path,
    outdir: Path,
    manifest: dict[str, Any],
    variant_match_mode: str,
    fuzzy_threshold: float,
) -> list[dict[str, Any]]:
    gene_results: list[dict[str, Any]] = []
    manifest_genes = manifest.get("genes") or {}

    for gene in genes:
        gold_path = gold_inputs.get(gene)
        if gold_path is None:
            note = (manifest_genes.get(gene) or {}).get(
                "note", "No normalized per-PMID recall input found."
            )
            gene_results.append({"gene": gene, "status": "missing_gold", "note": note})
            continue

        db_path = db_overrides.get(gene) or find_latest_db(gene, results_dir)
        if db_path is None:
            gene_results.append(
                {
                    "gene": gene,
                    "status": "missing_db",
                    "gold_path": str(gold_path),
                    "note": f"No {gene}*.db found under {results_dir / gene}",
                }
            )
            continue
        if not db_path.exists():
            gene_results.append(
                {
                    "gene": gene,
                    "status": "missing_db",
                    "gold_path": str(gold_path),
                    "db_path": str(db_path),
                    "note": "Configured DB path does not exist.",
                }
            )
            continue

        try:
            gene_results.append(
                run_gene_compare(
                    gene=gene,
                    gold_path=gold_path,
                    db_path=db_path,
                    outdir=outdir,
                    variant_match_mode=variant_match_mode,
                    fuzzy_threshold=fuzzy_threshold,
                )
            )
        except Exception as exc:  # noqa: BLE001
            logger.exception("Recall comparison failed for %s", gene)
            gene_results.append(
                {
                    "gene": gene,
                    "status": "error",
                    "gold_path": str(gold_path),
                    "db_path": str(db_path),
                    "error": str(exc),
                }
            )

    return gene_results


def main() -> int:
    parser = argparse.ArgumentParser(description="GVF multi-gene recall suite")
    parser.add_argument(
        "--gold-dir",
        type=Path,
        default=DEFAULT_GOLD_DIR,
        help="Gold-standard package root containing normalized/*_recall_input.csv",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Root results directory to search for GENE*.db files",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=METRICS_DIR / "current",
        help="Output directory for recall reports",
    )
    parser.add_argument(
        "--genes",
        help="Comma-separated gene list. Default: genes from gold inputs and manifest.",
    )
    parser.add_argument(
        "--db",
        action="append",
        help="Override a gene DB path, e.g. --db KCNQ1=/path/KCNQ1.db",
    )
    parser.add_argument(
        "--variant-match-mode",
        choices=["exact", "fuzzy"],
        default="fuzzy",
        help="Variant matching mode passed to cli.compare_variants",
    )
    parser.add_argument(
        "--fuzzy-threshold",
        type=float,
        default=0.80,
        help="Minimum fuzzy similarity passed to cli.compare_variants",
    )
    parser.add_argument(
        "--target",
        type=float,
        default=TARGET_RECALL,
        help="Target recall used for gap calculations",
    )
    parser.add_argument("--score", action="store_true", help="Score only; skip pytest")
    parser.add_argument("--log", action="store_true", help="Append score to history")
    parser.add_argument("--notes", default="", help="Notes for --log")
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    try:
        db_overrides = parse_db_overrides(args.db)
    except ValueError as exc:
        parser.error(str(exc))

    gold_dir = args.gold_dir.expanduser()
    results_dir = args.results_dir.expanduser()
    outdir = args.outdir.expanduser()

    gold_inputs = discover_gold_inputs(gold_dir)
    manifest = load_manifest(gold_dir)
    genes = selected_genes(parse_gene_list(args.genes), gold_inputs, manifest)

    if not genes:
        parser.error(
            f"No genes discovered. Expected CSVs under {gold_dir / 'normalized'}"
        )

    gene_results = build_gene_results(
        genes=genes,
        gold_inputs=gold_inputs,
        db_overrides=db_overrides,
        results_dir=results_dir,
        outdir=outdir,
        manifest=manifest,
        variant_match_mode=args.variant_match_mode,
        fuzzy_threshold=args.fuzzy_threshold,
    )
    scored = [item for item in gene_results if item.get("status") == "scored"]
    summary = build_run_summary(
        gene_results=gene_results,
        aggregate_recall=combine_recall(scored),
        aggregate_mae=combine_mae(scored),
        aggregate_precision=combine_precision(scored),
        manifest=manifest,
        target=args.target,
    )
    try:
        summary["disagreement_artifacts"] = write_disagreement_artifacts(
            outdir=outdir,
            gene_results=gene_results,
            results_dir=results_dir,
            gold_dir=gold_dir,
        )
    except Exception as exc:  # noqa: BLE001
        logger.exception("Paper-level disagreement report generation failed")
        summary["disagreement_artifacts_error"] = str(exc)
    write_run_summary(summary, outdir)
    print_scorecard(summary)

    if args.log:
        log_metrics(summary, args.notes)

    if args.score:
        return 0

    return run_tests()


if __name__ == "__main__":
    sys.exit(main())
