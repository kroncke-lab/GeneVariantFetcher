#!/usr/bin/env python3
"""Build the bounded MCP report payload for a completed Codex paper-eval run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

GENE_ORDER = ("SCN5A", "KCNH2", "KCNQ1", "RYR2")
COUNT_FIELDS = ("carriers", "affected", "unaffected")
FAILURE_ISSUES = {
    ("KCNQ1", "17192539"): "Narrative names one mutation; remaining identities absent",
    ("KCNQ1", "17470695"): "Structured preview lacks variant-level rows",
    ("KCNQ1", "14678125"): "Only aggregate mutation-region counts preserved",
    ("SCN5A", "26746457"): "Figure 1 carrier matrix needed figure/OCR fallback",
    ("KCNH2", "14661677"): "Partial article and unusable publisher artifacts",
    ("KCNH2", "19038855"): "Variant table image not transcribed into text",
    ("KCNH2", "24667783"): "Selected table has no usable variant identifiers",
    ("RYR2", "33686871"): "Pedigree OCR recovered only one variant label",
    ("RYR2", "19926015"): "Six residual variant misses",
}


def source_path(run_dir: Path, name: str) -> str:
    return str(run_dir / name)


def table_columns(*items: tuple[str, str, str | None]) -> list[dict]:
    columns = []
    for field, label, value_type in items:
        column = {"field": field, "label": label}
        if value_type == "percent":
            column["format"] = "percent"
        elif value_type:
            column["type"] = value_type
        columns.append(column)
    return columns


def format_percent(value: float | None) -> str:
    return "n/a" if value is None else f"{100 * value:.1f}%"


def build_payload(run_dir: Path) -> dict:
    report = json.loads((run_dir / "report.json").read_text())
    model_comparison_path = run_dir / "model_comparison.csv"
    if model_comparison_path.is_file():
        with model_comparison_path.open(newline="") as fh:
            model_rows = list(csv.DictReader(fh))
    else:
        model_rows = []

    overall = report["overall"]
    paper_count = overall["papers"]
    gold_variants = overall["tp"] + overall["fn"]
    predicted_variants = overall["tp"] + overall["fp"]
    title = f"Codex {paper_count}-Paper Extraction Evaluation"
    selection = report.get("selection") or {}
    selection_description = selection.get(
        "description",
        (
            "Paper selection used the recorded evaluation set. Routing, extraction, "
            "counts, evidence, and source locations were gold-value-blind."
        ),
    )
    gene_metrics = []
    gene_prf = []
    for gene in GENE_ORDER:
        metric = report["by_gene"][gene]
        row = {
            "gene": gene,
            **{field: metric[field] for field in ("tp", "fp", "fn")},
            **{field: metric[field] for field in ("precision", "recall", "f1")},
            **{
                f"{field}_count_recall": metric["count"][field]["recall"]
                for field in COUNT_FIELDS
            },
        }
        gene_metrics.append(row)
        for label, field in (
            ("Precision", "precision"),
            ("Recall", "recall"),
            ("F1", "f1"),
        ):
            gene_prf.append(
                {
                    "gene": gene,
                    "metric": label,
                    "value": metric[field],
                    "tp": metric["tp"],
                    "fp": metric["fp"],
                    "fn": metric["fn"],
                    "gold_variants": metric["tp"] + metric["fn"],
                }
            )

    count_metrics = []
    for field in COUNT_FIELDS:
        metric = overall["count"][field]
        wrong = sum(
            1
            for paper in report["papers"]
            for error in paper["count_errors"]
            if error["field"] == field
        )
        count_metrics.append(
            {
                "field": field.title(),
                "supplied": metric["predicted"],
                "gold_asserted": metric["gold_asserted"],
                "count_recall": metric["recall"],
                "exact_rate": (
                    (metric["predicted"] - wrong) / metric["predicted"]
                    if metric["predicted"]
                    else None
                ),
                "mae": metric["mae"],
                "rmse": metric["rmse"],
            }
        )

    paper_metrics = []
    for paper in report["papers"]:
        paper_metrics.append(
            {
                "gene": paper["gene"],
                "pmid": str(paper["pmid"]),
                "tool": paper["tool"],
                **{field: paper[field] for field in ("tp", "fp", "fn")},
                **{field: paper[field] for field in ("precision", "recall", "f1")},
                **{
                    f"{field}_count_recall": paper["count"][field]["recall"]
                    for field in COUNT_FIELDS
                },
                "seconds": paper["elapsed_seconds"],
                "tokens": paper["token_usage"]["total_tokens"],
            }
        )
    ranked_failures = sorted(
        (paper for paper in report["papers"] if paper["fn"] > 0),
        key=lambda row: (-row["fn"], row["gene"], str(row["pmid"])),
    )
    failures = []
    for paper in ranked_failures[:9]:
        failures.append(
            {
                "gene": paper["gene"],
                "pmid": str(paper["pmid"]),
                "tool": paper["tool"],
                "fn": paper["fn"],
                "recall": paper["recall"],
                "dominant_issue": FAILURE_ISSUES.get(
                    (paper["gene"], str(paper["pmid"])),
                    "Residual variant misses",
                ),
            }
        )
    displayed_failure_fn = sum(row["fn"] for row in failures)
    displayed_failure_share = (
        displayed_failure_fn / overall["fn"] if overall["fn"] else 0.0
    )
    concentrated_failures = ranked_failures[:7]
    concentrated_failure_fn = sum(row["fn"] for row in concentrated_failures)
    concentrated_failure_share = (
        concentrated_failure_fn / overall["fn"] if overall["fn"] else 0.0
    )

    typed_models = []
    for row in model_rows:
        typed_models.append(
            {
                "system": row["system"],
                "score_view": row["score_view"],
                **{field: float(row[field]) for field in ("precision", "recall", "f1")},
                "extracted": int(row["extracted_variants"]),
                "total_tokens": int(row["total_tokens"]),
                "tokens_per_paper": float(row["tokens_per_paper"]),
                "wall_minutes": float(row["elapsed_minutes"]),
            }
        )
    comparator_count = sum(
        1 for row in typed_models if not row["system"].lower().startswith("codex")
    )
    if typed_models:
        comparison_body = (
            "## Model comparison\n\nThe table includes "
            f"{comparator_count} non-Codex comparator results alongside the "
            f"Codex score of {format_percent(overall['recall'])} recall. Matcher, "
            "token-accounting, and concurrency differences limit direct speed, "
            "cost, or quality ranking."
        )
    else:
        comparison_body = (
            "## Model comparison\n\nNo `model_comparison.csv` was supplied for "
            "this run, so no cross-model ranking is shown."
        )

    overview = {
        "precision": overall["precision"],
        "recall": overall["recall"],
        "f1": overall["f1"],
        **{field: overall[field] for field in ("tp", "fp", "fn")},
        "gold_variants": gold_variants,
        "predicted_variants": predicted_variants,
        **{
            f"{field}_count_recall": overall["count"][field]["recall"]
            for field in COUNT_FIELDS
        },
        "total_tokens": report["token_usage"]["total_tokens"],
        "tokens_per_paper": report["token_usage"]["total_tokens"] / overall["papers"],
        "wall_minutes": report["timing"]["wall_seconds"] / 60,
    }
    recall_rank = sorted(
        (
            (gene, report["by_gene"][gene]["recall"])
            for gene in GENE_ORDER
            if report["by_gene"][gene]["recall"] is not None
        ),
        key=lambda item: item[1],
        reverse=True,
    )
    strongest_gene, strongest_recall = recall_rank[0]
    weakest_gene, weakest_recall = recall_rank[-1]
    weakest_top_three_fn = sum(
        paper["fn"]
        for paper in sorted(
            (row for row in report["papers"] if row["gene"] == weakest_gene),
            key=lambda row: -row["fn"],
        )[:3]
    )
    audit = report.get("scoring_audit")
    if audit:
        raw = audit["pre_adjudication"]
        scorer_body = (
            "## Scorer sensitivity\n\nThe preserved pre-adjudication matcher "
            f"scored {format_percent(raw['recall'])} recall and "
            f"{format_percent(raw['f1'])} F1. A post-lock audit recovered "
            f"{audit['recovered_equivalent_matches']} notation-equivalent labels "
            "without editing predictions, producing "
            f"{format_percent(overall['recall'])} recall and "
            f"{format_percent(overall['f1'])} F1."
        )
        adjudication_rows = audit["recovered_equivalent_matches"]
        scorer_evidence = (
            "raw and audited scorer outputs, "
            f"a {adjudication_rows}-row matcher adjudication ledger"
        )
    else:
        scorer_body = (
            "## Scorer sensitivity\n\nNo separate pre-adjudication matcher result "
            "was supplied for this run; the report shows the locked scorer output."
        )
        scorer_evidence = "the locked scorer output"

    manifest = {
        "version": 1,
        "title": title,
        "surface": "report",
        "description": (
            f"Extraction-blinded evaluation of {paper_count} cardiac-gene papers with "
            "variant and phenotype-count evidence, scorer sensitivity, failure "
            "concentration, and model comparison."
        ),
        "generatedAt": report["scored_at"],
        "sources": [
            {
                "id": "report_sql",
                "label": "Reproducible DuckDB report queries",
                "path": source_path(run_dir, "report_queries.sql"),
            },
            {
                "id": "codex_report",
                "label": "Codex locked evaluation report",
                "path": source_path(run_dir, "report.json"),
            },
            {
                "id": "codex_evidence",
                "label": "Codex variant evidence ledger",
                "path": source_path(run_dir, "evidence.csv"),
            },
            {
                "id": "codex_validation",
                "label": "Independent validation notes",
                "path": source_path(run_dir, "validation_notes.md"),
            },
        ],
        "cards": [
            {
                "id": "variant_recall_card",
                "description": (
                    "Micro-averaged audited variant recall, with F1 as paired context."
                ),
                "dataset": "overview",
                "sourceId": "report_sql",
                "metrics": [
                    {"label": "Variant recall", "field": "recall", "format": "percent"},
                    {"label": "F1", "field": "f1", "format": "percent"},
                ],
            },
            {
                "id": "precision_card",
                "description": (
                    f"Share of the {predicted_variants:,} locked predictions matched "
                    "to a gold row "
                    "after notation adjudication."
                ),
                "dataset": "overview",
                "sourceId": "report_sql",
                "metrics": [
                    {
                        "label": "Variant precision",
                        "field": "precision",
                        "format": "percent",
                    },
                    {"label": "Matched variants", "field": "tp", "format": "number"},
                ],
            },
            {
                "id": "count_coverage_card",
                "description": (
                    "Share of gold count assertions supplied by the locked "
                    "extraction, reported separately for each count field."
                ),
                "dataset": "overview",
                "sourceId": "report_sql",
                "metrics": [
                    {
                        "label": "Carrier count recall",
                        "field": "carriers_count_recall",
                        "format": "percent",
                    },
                    {
                        "label": "Affected",
                        "field": "affected_count_recall",
                        "format": "percent",
                    },
                    {
                        "label": "Unaffected",
                        "field": "unaffected_count_recall",
                        "format": "percent",
                    },
                ],
            },
            {
                "id": "usage_card",
                "description": (
                    "Exact API token telemetry recorded across routing and "
                    "extraction responses."
                ),
                "dataset": "overview",
                "sourceId": "report_sql",
                "metrics": [
                    {
                        "label": "Total tokens",
                        "field": "total_tokens",
                        "format": "compact",
                    },
                    {
                        "label": "Tokens per paper",
                        "field": "tokens_per_paper",
                        "format": "compact",
                    },
                    {
                        "label": "Wall minutes",
                        "field": "wall_minutes",
                        "format": "number",
                    },
                ],
            },
        ],
        "charts": [
            {
                "id": "gene_prf_chart",
                "title": "Variant precision, recall, and F1 by gene",
                "subtitle": (
                    "KCNQ1 is the lowest-recall gene; RYR2 is the strongest "
                    "Codex result."
                ),
                "type": "bar",
                "dataset": "gene_prf",
                "sourceId": "report_sql",
                "valueFormat": "percent",
                "encodings": {
                    "x": {"field": "gene", "type": "nominal", "label": "Gene"},
                    "y": {
                        "field": "value",
                        "type": "quantitative",
                        "label": "Rate",
                        "format": "percent",
                    },
                    "color": {
                        "field": "metric",
                        "type": "nominal",
                        "label": "Metric",
                    },
                    "tooltip": [
                        {"field": "tp", "type": "quantitative", "label": "TP"},
                        {"field": "fp", "type": "quantitative", "label": "FP"},
                        {"field": "fn", "type": "quantitative", "label": "FN"},
                        {
                            "field": "gold_variants",
                            "type": "quantitative",
                            "label": "Gold variants",
                        },
                    ],
                },
                "xAxisTitle": "Gene",
                "yAxisTitle": "Rate",
            }
        ],
        "tables": [
            {
                "id": "gene_table",
                "title": "Per-gene metrics",
                "subtitle": "Micro-averaged variant and count coverage metrics.",
                "dataset": "gene_metrics",
                "sourceId": "report_sql",
                "defaultSort": {"field": "recall", "direction": "desc"},
                "columns": table_columns(
                    ("gene", "Gene", "text"),
                    ("tp", "TP", "number"),
                    ("fp", "FP", "number"),
                    ("fn", "FN", "number"),
                    ("precision", "Precision", "percent"),
                    ("recall", "Recall", "percent"),
                    ("f1", "F1", "percent"),
                    ("carriers_count_recall", "Carrier count recall", "percent"),
                    ("affected_count_recall", "Affected count recall", "percent"),
                    (
                        "unaffected_count_recall",
                        "Unaffected count recall",
                        "percent",
                    ),
                ),
            },
            {
                "id": "count_table",
                "title": "Count coverage and conditional error",
                "subtitle": (
                    "MAE and RMSE use only matched rows where both Codex and "
                    "gold supplied the field."
                ),
                "dataset": "count_metrics",
                "sourceId": "report_sql",
                "defaultSort": {"field": "count_recall", "direction": "desc"},
                "columns": table_columns(
                    ("field", "Field", "text"),
                    ("supplied", "Supplied", "number"),
                    ("gold_asserted", "Gold assertions", "number"),
                    ("count_recall", "Count recall", "percent"),
                    ("exact_rate", "Exact among supplied", "percent"),
                    ("mae", "MAE", "number"),
                    ("rmse", "RMSE", "number"),
                ),
            },
            {
                "id": "failure_table",
                "title": "Papers with the most missed gold variants",
                "subtitle": (
                    f"The displayed {len(failures)} papers account for "
                    f"{format_percent(displayed_failure_share)} of all false negatives."
                ),
                "dataset": "failure_metrics",
                "sourceId": "report_sql",
                "defaultSort": {"field": "fn", "direction": "desc"},
                "columns": table_columns(
                    ("gene", "Gene", "text"),
                    ("pmid", "PMID", "text"),
                    ("tool", "Tool", "text"),
                    ("fn", "FN", "number"),
                    ("recall", "Recall", "percent"),
                    ("dominant_issue", "Dominant issue", "text"),
                ),
            },
            {
                "id": "paper_table",
                "title": "Per-paper extraction metrics",
                "subtitle": (
                    f"All {paper_count} papers; initially sorted by missed gold variants."
                ),
                "dataset": "paper_metrics",
                "sourceId": "report_sql",
                "defaultSort": {"field": "fn", "direction": "desc"},
                "columns": table_columns(
                    ("gene", "Gene", "text"),
                    ("pmid", "PMID", "text"),
                    ("tool", "Tool", "text"),
                    ("tp", "TP", "number"),
                    ("fp", "FP", "number"),
                    ("fn", "FN", "number"),
                    ("precision", "Precision", "percent"),
                    ("recall", "Recall", "percent"),
                    ("f1", "F1", "percent"),
                    ("carriers_count_recall", "Carrier count recall", "percent"),
                    ("affected_count_recall", "Affected count recall", "percent"),
                    (
                        "unaffected_count_recall",
                        "Unaffected count recall",
                        "percent",
                    ),
                    ("seconds", "Seconds", "number"),
                    ("tokens", "Tokens", "number"),
                ),
            },
            {
                "id": "model_table",
                "title": "Aggregate model comparison",
                "subtitle": (
                    "Matcher, token, and concurrency differences limit direct ranking."
                ),
                "dataset": "model_metrics",
                "sourceId": "report_sql",
                "defaultSort": {"field": "recall", "direction": "desc"},
                "columns": table_columns(
                    ("system", "System", "text"),
                    ("score_view", "Score view", "text"),
                    ("precision", "Precision", "percent"),
                    ("recall", "Recall", "percent"),
                    ("f1", "F1", "percent"),
                    ("extracted", "Extracted", "number"),
                    ("total_tokens", "Total tokens", "number"),
                    ("tokens_per_paper", "Tokens / paper", "number"),
                    ("wall_minutes", "Wall minutes", "number"),
                ),
            },
        ],
        "blocks": [
            {
                "id": "title",
                "type": "markdown",
                "body": f"# {title}",
            },
            {
                "id": "executive_summary",
                "type": "markdown",
                "body": (
                    "## Executive Summary\n\n"
                    "- The hash-locked, extraction-blinded Codex run scored "
                    f"**{format_percent(overall['precision'])} precision, "
                    f"{format_percent(overall['recall'])} recall, and "
                    f"{format_percent(overall['f1'])} F1** over "
                    f"{gold_variants:,} gold variant rows.\n"
                    "- Count recall was "
                    f"**{format_percent(overall['count']['carriers']['recall'])} "
                    "carrier**, "
                    f"**{format_percent(overall['count']['affected']['recall'])} "
                    "affected**, and "
                    f"**{format_percent(overall['count']['unaffected']['recall'])} "
                    "unaffected**.\n"
                    f"- The top {len(concentrated_failures)} papers produced "
                    f"**{format_percent(concentrated_failure_share)} of all false "
                    "negatives**.\n"
                    f"- Codex used **{report['token_usage']['total_tokens']:,} exact "
                    "API tokens** and "
                    f"**{report['timing']['wall_seconds'] / 60:.1f} wall minutes**."
                ),
            },
            {
                "id": "headline_metrics",
                "type": "metric-strip",
                "cardIds": [
                    "variant_recall_card",
                    "precision_card",
                    "count_coverage_card",
                    "usage_card",
                ],
            },
            {
                "id": "gene_heading",
                "type": "markdown",
                "body": (
                    f"## Per-gene performance\n\n{strongest_gene} was strongest at "
                    f"{format_percent(strongest_recall)} recall. {weakest_gene} "
                    f"was weakest at {format_percent(weakest_recall)}, with its "
                    f"three highest-miss papers contributing {weakest_top_three_fn} "
                    "false negatives."
                ),
            },
            {"id": "gene_chart_block", "type": "chart", "chartId": "gene_prf_chart"},
            {"id": "gene_table_block", "type": "table", "tableId": "gene_table"},
            {
                "id": "count_heading",
                "type": "markdown",
                "body": (
                    "## Count recall and error\n\nCount recall is the share of "
                    "all gold assertions for which a matched Codex variant "
                    "supplied a value. MAE and RMSE are conditional on supplied "
                    "values, so they must not be interpreted without coverage."
                ),
            },
            {"id": "count_table_block", "type": "table", "tableId": "count_table"},
            {
                "id": "errors_heading",
                "type": "markdown",
                "body": (
                    f"## Failure concentration\n\nThe top "
                    f"{len(concentrated_failures)} papers account for "
                    f"{concentrated_failure_fn} of {overall['fn']} false negatives "
                    f"({format_percent(concentrated_failure_share)}). Review the "
                    "source-completeness and routing evidence for these papers "
                    "before attributing the misses solely to extraction."
                ),
            },
            {
                "id": "failure_table_block",
                "type": "table",
                "tableId": "failure_table",
            },
            {
                "id": "paper_heading",
                "type": "markdown",
                "body": (
                    "## Per-paper results\n\nThe table preserves the requested "
                    "paper-level precision, recall, F1, count coverage, selected "
                    "tool, elapsed time, and token use."
                ),
            },
            {"id": "paper_table_block", "type": "table", "tableId": "paper_table"},
            {
                "id": "comparison_heading",
                "type": "markdown",
                "body": comparison_body,
            },
            {"id": "model_table_block", "type": "table", "tableId": "model_table"},
            {
                "id": "method_heading",
                "type": "markdown",
                "body": (
                    f"## Blinding and scorer method\n\n{selection_description} "
                    "Gold values and row counts were absent from routing and "
                    "extraction. `predictions.json` was made "
                    "read-only and SHA-256 locked before scoring opened the gold "
                    "CSVs."
                ),
            },
            {
                "id": "scorer_heading",
                "type": "markdown",
                "body": scorer_body,
            },
            {
                "id": "limitations_heading",
                "type": "markdown",
                "body": (
                    "## Limitations and next step\n\nThe curated gold is "
                    "incomplete in some papers, unmatched predictions are "
                    "therefore soft false positives, and several count "
                    "disagreements reflect cohort or phenotype-definition scope "
                    "rather than arithmetic. Add a data-rich-paper fallback that "
                    "retries table, PDF, and figure/OCR representations when the "
                    "first route returns zero or implausibly few variants."
                ),
            },
            {
                "id": "evidence_heading",
                "type": "markdown",
                "body": (
                    "## Evidence and reproducibility\n\nThe run directory "
                    "contains the immutable prediction file, source hashes, lock "
                    f"hashes, {predicted_variants:,}-row evidence ledger, "
                    f"{paper_count}-row paper metric table, {scorer_evidence}, "
                    "validation notes, and any supplied model comparison."
                ),
            },
        ],
    }

    return {
        "surface": "report",
        "manifest": manifest,
        "snapshot": {
            "version": 1,
            "generatedAt": report["scored_at"],
            "status": "ready",
            "datasets": {
                "overview": [overview],
                "gene_prf": gene_prf,
                "gene_metrics": gene_metrics,
                "count_metrics": count_metrics,
                "failure_metrics": failures,
                "paper_metrics": paper_metrics,
                "model_metrics": typed_models,
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_payload(args.run_dir), separators=(",", ":")))


if __name__ == "__main__":
    main()
