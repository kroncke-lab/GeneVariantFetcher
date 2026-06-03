#!/usr/bin/env python3
"""Build paper-level disagreement reports from scored recall artifacts."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional

try:
    from scripts.recall_audit.common import (
        DEFAULT_RESULTS_DIR,
        canonical_variant,
        context_status,
        find_full_contexts,
        normalize_pmid,
        parse_bool,
        parse_number,
        read_csv_rows,
        repo_path,
        resolve_gold_path,
        resolve_recall_score,
        resolve_run_dir,
        source_has_unresolved_variant_supplement_refs,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        DEFAULT_RESULTS_DIR,
        canonical_variant,
        context_status,
        find_full_contexts,
        normalize_pmid,
        parse_bool,
        parse_number,
        read_csv_rows,
        repo_path,
        resolve_gold_path,
        resolve_recall_score,
        resolve_run_dir,
        source_has_unresolved_variant_supplement_refs,
        write_csv_rows,
    )


DATA_AVAILABLE_STATUSES = {
    "recovered_pmc",
    "recovered_browser",
    "recovered_supplement_only",
}

REPORT_FIELDS = [
    "gene",
    "pmid",
    "title",
    "source_status",
    "data_available",
    "source_file",
    "context_path",
    "context_bytes",
    "available_source_status",
    "available_context_path",
    "available_context_bytes",
    "source_desync",
    "source_unbound",
    "table_references",
    "table_body_markers",
    "likely_missing_table_bodies",
    "likely_missing_variant_supplement",
    "gold_rows",
    "matched_rows",
    "row_recall",
    "missing_rows",
    "missing_carriers",
    "missing_affected",
    "missing_unaffected",
    "count_mismatch_rows",
    "abs_carrier_diff",
    "abs_affected_diff",
    "abs_unaffected_diff",
    "extra_sqlite_rows",
    "disagreement_score",
    "failure_class",
    "top_missing_variants",
    "top_count_mismatch_variants",
]


def _sum_numeric(rows: Iterable[dict[str, str]], key: str) -> float:
    total = 0.0
    for row in rows:
        value = parse_number(row.get(key))
        if value is not None:
            total += value
    return total


def _abs_sum_numeric(rows: Iterable[dict[str, str]], key: str) -> float:
    total = 0.0
    for row in rows:
        value = parse_number(row.get(key))
        if value is not None:
            total += abs(value)
    return total


def _top_values(rows: Iterable[dict[str, str]], key: str, limit: int = 8) -> str:
    values = []
    seen = set()
    for row in rows:
        value = (row.get(key) or "").strip()
        if not value or value in seen:
            continue
        values.append(value)
        seen.add(value)
        if len(values) >= limit:
            break
    return ";".join(values)


def _gold_by_pmid(gold_path: Path) -> dict[str, dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    seen: set[tuple[str, str]] = set()
    for row in read_csv_rows(gold_path):
        pmid = normalize_pmid(row.get("pmid"))
        variant = canonical_variant(row.get("variant"))
        if not pmid or not variant or (pmid, variant) in seen:
            continue
        seen.add((pmid, variant))
        item = grouped.setdefault(
            pmid,
            {
                "gold_rows": 0,
                "gold_carriers": 0.0,
                "gold_affected": 0.0,
                "gold_unaffected": 0.0,
            },
        )
        item["gold_rows"] += 1
        item["gold_carriers"] += parse_number(row.get("carriers")) or 0.0
        item["gold_affected"] += parse_number(row.get("affected")) or 0.0
        item["gold_unaffected"] += parse_number(row.get("unaffected")) or 0.0
    return grouped


def _group_rows_by_pmid(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[normalize_pmid(row.get("pmid"))].append(row)
    return grouped


def _load_metadata(db_path: Path | None) -> dict[str, dict[str, Any]]:
    if not db_path or not db_path.exists():
        return {}
    sql = """
        SELECT
            p.pmid,
            p.title,
            em.source_file,
            em.source_type,
            em.abstract_only,
            em.total_variants_found,
            em.extraction_confidence,
            em.model_used,
            em.challenges,
            em.notes
        FROM papers p
        LEFT JOIN extraction_metadata em ON em.pmid = p.pmid
    """
    grouped: dict[str, dict[str, Any]] = {}
    with sqlite3.connect(f"file:{db_path}?mode=ro", uri=True) as conn:
        conn.row_factory = sqlite3.Row
        for row in conn.execute(sql):
            item = dict(row)
            pmid = normalize_pmid(item.get("pmid"))
            if not pmid:
                continue
            old = grouped.get(pmid)
            if old is None:
                grouped[pmid] = item
                continue
            old_score = int(old.get("total_variants_found") or 0)
            new_score = int(item.get("total_variants_found") or 0)
            if new_score >= old_score:
                grouped[pmid] = item
    return grouped


def _autodetect_results_dir(run_dir: Path | None) -> Path:
    if run_dir is not None:
        for candidate in (run_dir / "results", run_dir.parent / "results"):
            if candidate.exists():
                return candidate
    return DEFAULT_RESULTS_DIR


def _build_context_index(results_dir: Path) -> dict[tuple[str, str], Path]:
    index: dict[tuple[str, str], Path] = {}
    if not results_dir.exists():
        return index
    for path in results_dir.rglob("*_FULL_CONTEXT.md"):
        if not path.exists():
            continue
        pmid = normalize_pmid(path.name.removesuffix("_FULL_CONTEXT.md"))
        if not pmid:
            continue
        gene = path.relative_to(results_dir).parts[0].upper()
        key = (gene, pmid)
        current = index.get(key)
        current_size = current.stat().st_size if current and current.exists() else -1
        if path.stat().st_size > current_size:
            index[key] = path
    return index


def _is_extraction_source_path(path: Path | None) -> bool:
    if path is None:
        return False
    return path.name.endswith(
        (
            "_FULL_CONTEXT.md",
            "_CLEANED.md",
            "_DATA_ZONES.md",
        )
    )


def _metadata_mentions_abstract_only(metadata: dict[str, Any] | None) -> bool:
    if not metadata:
        return False
    fields: list[str] = []
    for key in ("source_type", "source_kind", "notes", "challenges"):
        value = metadata.get(key)
        if value:
            fields.append(str(value))
    joined = " ".join(fields).lower()
    return (
        bool(metadata.get("abstract_only"))
        or "abstract-only" in joined
        or "abstract only" in joined
        or "full text unavailable" in joined
        or "full text not available" in joined
        or "full text could not be retrieved" in joined
    )


def _table_body_metrics(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {
            "table_references": 0,
            "table_body_markers": 0,
            "likely_missing_table_bodies": False,
        }
    text = path.read_text(encoding="utf-8", errors="ignore")
    table_references = len(
        re.findall(r"\bTables?\s+\d+(?:\s*(?:and|,|-)\s*\d+)?", text)
    )
    table_body_markers = len(
        re.findall(r"(?m)^\s*\|.+\|\s*$", text)
        + re.findall(r"(?m)^\s*(?:#{1,6}\s*)?(?:Table|TABLE)\s+\d+[.:]", text)
        + re.findall(r"(?i)<table\b", text)
    )
    return {
        "table_references": table_references,
        "table_body_markers": table_body_markers,
        "likely_missing_table_bodies": (
            table_references >= 3
            and table_body_markers == 0
            and path.stat().st_size > 5_000
        ),
    }


def _sibling_full_context(source_path: Path | None, pmid: str) -> Path | None:
    if source_path is None or not source_path.exists():
        return None
    candidates: list[Path] = []
    parent = source_path.parent
    run_dir = parent.parent if parent.name in {"extractions", "pmc_fulltext"} else None
    if run_dir is not None:
        candidates.append(run_dir / "pmc_fulltext" / f"{pmid}_FULL_CONTEXT.md")
    if source_path.name.endswith("_CLEANED.md") or source_path.name.endswith(
        "_DATA_ZONES.md"
    ):
        candidates.append(parent / f"{pmid}_FULL_CONTEXT.md")
    existing = [path for path in candidates if path.exists()]
    if not existing:
        return None
    return max(existing, key=lambda path: path.stat().st_size)


def _source_info(
    *,
    gene: str,
    pmid: str,
    metadata: dict[str, Any] | None,
    results_dir: Path,
    context_index: dict[tuple[str, str], Path] | None = None,
) -> dict[str, Any]:
    source_file = (metadata or {}).get("source_file") or ""
    source_path = repo_path(source_file) if source_file else None
    context_path: Optional[Path] = None
    source_unbound = bool(
        source_path
        and source_path.exists()
        and not _is_extraction_source_path(source_path)
    )

    if source_path and source_path.exists() and _is_extraction_source_path(source_path):
        context_path = source_path

    available_context = _sibling_full_context(source_path, pmid)
    if available_context is None:
        available_context = (context_index or {}).get((gene.upper(), pmid))
    if available_context is None:
        contexts = find_full_contexts(results_dir, gene, pmid, global_search=False)
        available_context = contexts[0] if contexts else None

    status = context_status(context_path)
    metadata_abstract_only = _metadata_mentions_abstract_only(metadata)
    if metadata_abstract_only and status == "not_attempted":
        status = "abstract_only"
    table_metrics = _table_body_metrics(context_path)
    missing_variant_supplement = source_has_unresolved_variant_supplement_refs(
        context_path, gene
    )
    if status in DATA_AVAILABLE_STATUSES and (
        table_metrics["likely_missing_table_bodies"] or missing_variant_supplement
    ):
        status = "missing_table_bodies"
    context_bytes = context_path.stat().st_size if context_path else 0
    available_status = context_status(available_context)
    available_bytes = available_context.stat().st_size if available_context else 0
    source_desync = False
    if (
        available_context is not None
        and available_status in DATA_AVAILABLE_STATUSES
        and status not in DATA_AVAILABLE_STATUSES
        and metadata_abstract_only
    ):
        source_desync = True
    elif (
        available_context is not None
        and available_status in DATA_AVAILABLE_STATUSES
        and context_path is not None
        and available_context.resolve() != context_path.resolve()
        and available_bytes > context_bytes + 2_048
    ):
        source_desync = True
    return {
        "source_status": status,
        "data_available": status in DATA_AVAILABLE_STATUSES,
        "source_file": source_file,
        "context_path": str(context_path) if context_path else "",
        "context_bytes": context_bytes,
        "available_source_status": available_status,
        "available_context_path": str(available_context) if available_context else "",
        "available_context_bytes": available_bytes,
        "source_desync": source_desync,
        "source_unbound": source_unbound,
        "likely_missing_variant_supplement": missing_variant_supplement,
        **table_metrics,
    }


def _failure_class(row: dict[str, Any]) -> str:
    data_available = bool(row.get("data_available"))
    missing_rows = int(row.get("missing_rows") or 0)
    gold_rows = int(row.get("gold_rows") or 0)
    count_mismatch_rows = int(row.get("count_mismatch_rows") or 0)
    extra_rows = int(row.get("extra_sqlite_rows") or 0)
    abs_aff = float(row.get("abs_affected_diff") or 0)
    abs_unaff = float(row.get("abs_unaffected_diff") or 0)
    source_status = row.get("source_status")
    source_desync = bool(row.get("source_desync"))
    source_unbound = bool(row.get("source_unbound"))
    available_status = row.get("available_source_status")

    if source_desync and missing_rows:
        return "stale_source_desync"
    if source_unbound and available_status in DATA_AVAILABLE_STATUSES and missing_rows:
        return "source_unbound_available"
    if source_status == "missing_table_bodies" and (
        missing_rows or count_mismatch_rows
    ):
        return "source_missing_table_bodies"
    if not data_available and missing_rows:
        if source_status == "abstract_only":
            return "source_abstract_only"
        if source_status in {"paywall_stub", "not_attempted"}:
            return "source_missing_or_stub"
        return "source_not_usable"
    if missing_rows >= max(3, math.ceil(gold_rows * 0.50)):
        return "available_source_underextraction"
    if count_mismatch_rows and abs_aff > 5 and abs_unaff > 5:
        return "phenotype_or_count_direction"
    if count_mismatch_rows:
        return "count_semantics"
    if extra_rows >= max(10, gold_rows):
        return "overinclusive_extraction"
    if missing_rows:
        return "partial_underextraction"
    if extra_rows:
        return "pipeline_only_extras"
    return "matched_or_low_priority"


def _score(row: dict[str, Any]) -> float:
    score = (
        int(row.get("missing_rows") or 0) * 3
        + int(row.get("count_mismatch_rows") or 0) * 2
        + (float(row.get("abs_carrier_diff") or 0) / 10)
        + (float(row.get("abs_affected_diff") or 0) / 10)
        + (float(row.get("abs_unaffected_diff") or 0) / 10)
        + (int(row.get("extra_sqlite_rows") or 0) / 25)
    )
    if row.get("data_available"):
        score += 25
    if row.get("source_desync"):
        score += 35
    if (
        row.get("source_unbound")
        and row.get("available_source_status") in DATA_AVAILABLE_STATUSES
    ):
        score += 20
    return round(score, 3)


def build_paper_disagreement_rows(
    *,
    recall_score: Path,
    db_dir: Path | None = None,
    db_paths_by_gene: dict[str, Path] | None = None,
    gold_dir: Path | None = None,
    results_dir: Path | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    db_root = db_dir if db_dir else recall_score.parent / "dbs"
    result_root = (
        results_dir if results_dir else _autodetect_results_dir(recall_score.parent)
    )
    context_index = _build_context_index(result_root)

    for gene_dir in sorted(path for path in recall_score.iterdir() if path.is_dir()):
        gene = gene_dir.name.upper()
        gold_path = (
            gold_dir / f"{gene}_recall_input.csv"
            if gold_dir
            else resolve_gold_path(gene)
        )
        if not gold_path.exists():
            continue
        gold = _gold_by_pmid(gold_path)
        missing_by = _group_rows_by_pmid(
            read_csv_rows(gene_dir / "missing_in_sqlite.csv")
        )
        extra_rows = read_csv_rows(gene_dir / "missing_in_excel.csv")
        extra_by = _group_rows_by_pmid(
            [row for row in extra_rows if normalize_pmid(row.get("pmid")) in gold]
        )
        discrepancies = read_csv_rows(gene_dir / "discrepancies.csv")
        mismatch_by = _group_rows_by_pmid(
            [
                row
                for row in discrepancies
                if parse_bool(row.get("count_mismatch"))
                and not parse_bool(row.get("missing_in_sqlite"))
                and not parse_bool(row.get("missing_in_excel"))
            ]
        )
        db_path = (
            db_paths_by_gene.get(gene)
            if db_paths_by_gene is not None
            else db_root / f"{gene}.db"
        )
        metadata = _load_metadata(db_path)

        for pmid, gold_item in sorted(gold.items()):
            missing = missing_by.get(pmid, [])
            mismatches = mismatch_by.get(pmid, [])
            extras = extra_by.get(pmid, [])
            source = _source_info(
                gene=gene,
                pmid=pmid,
                metadata=metadata.get(pmid),
                results_dir=result_root,
                context_index=context_index,
            )
            gold_rows = int(gold_item["gold_rows"])
            missing_rows = len(missing)
            matched_rows = max(0, gold_rows - missing_rows)
            row: dict[str, Any] = {
                "gene": gene,
                "pmid": pmid,
                "title": (metadata.get(pmid) or {}).get("title") or "",
                **source,
                "gold_rows": gold_rows,
                "matched_rows": matched_rows,
                "row_recall": round(matched_rows / gold_rows, 4) if gold_rows else "",
                "missing_rows": missing_rows,
                "missing_carriers": round(
                    _sum_numeric(missing, "excel_carriers_total"), 3
                ),
                "missing_affected": round(_sum_numeric(missing, "excel_affected"), 3),
                "missing_unaffected": round(
                    _sum_numeric(missing, "excel_unaffected"), 3
                ),
                "count_mismatch_rows": len(mismatches),
                "abs_carrier_diff": round(
                    _abs_sum_numeric(mismatches, "carriers_diff"), 3
                ),
                "abs_affected_diff": round(
                    _abs_sum_numeric(mismatches, "affected_diff"), 3
                ),
                "abs_unaffected_diff": round(
                    _abs_sum_numeric(mismatches, "unaffected_diff"), 3
                ),
                "extra_sqlite_rows": len(extras),
                "top_missing_variants": _top_values(missing, "excel_variant_raw"),
                "top_count_mismatch_variants": _top_values(
                    mismatches, "excel_variant_raw"
                ),
            }
            row["failure_class"] = _failure_class(row)
            row["disagreement_score"] = _score(row)
            rows.append(row)

    return sorted(
        rows,
        key=lambda row: (
            not (
                bool(row.get("data_available"))
                or bool(row.get("source_desync"))
                or (
                    bool(row.get("source_unbound"))
                    and row.get("available_source_status") in DATA_AVAILABLE_STATUSES
                )
            ),
            -float(row.get("disagreement_score") or 0),
            row.get("gene", ""),
            row.get("pmid", ""),
        ),
    )


def select_regression_cases(
    rows: list[dict[str, Any]], *, limit: int = 36, min_per_gene: int = 6
) -> list[dict[str, Any]]:
    eligible = [
        row
        for row in rows
        if row.get("data_available")
        and row.get("failure_class") not in {"matched_or_low_priority"}
        and float(row.get("disagreement_score") or 0) > 0
    ]
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in eligible:
        by_gene[row["gene"]].append(row)

    selected: list[dict[str, Any]] = []
    selected_keys: set[tuple[str, str]] = set()
    for gene in sorted(by_gene):
        for row in by_gene[gene][:min_per_gene]:
            if len(selected) >= limit:
                break
            key = (row["gene"], row["pmid"])
            selected.append(row)
            selected_keys.add(key)

    for row in eligible:
        if len(selected) >= limit:
            break
        key = (row["gene"], row["pmid"])
        if key in selected_keys:
            continue
        selected.append(row)
        selected_keys.add(key)
    return selected


def write_markdown(rows: list[dict[str, Any]], out_path: Path, *, title: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# {title}",
        "",
        "| Gene | PMID | Class | Row Recall | Missing | Count Mismatch | Extra | Source | Top Variants |",
        "|---|---:|---|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        recall = row.get("row_recall")
        recall_text = "" if recall == "" else f"{float(recall):.1%}"
        variants = row.get("top_missing_variants") or row.get(
            "top_count_mismatch_variants"
        )
        lines.append(
            "| {gene} | {pmid} | {failure_class} | {recall} | {missing_rows} | "
            "{count_mismatch_rows} | {extra_sqlite_rows} | {source_status} | {variants} |".format(
                gene=row.get("gene", ""),
                pmid=row.get("pmid", ""),
                failure_class=row.get("failure_class", ""),
                recall=recall_text,
                missing_rows=row.get("missing_rows", 0),
                count_mismatch_rows=row.get("count_mismatch_rows", 0),
                extra_sqlite_rows=row.get("extra_sqlite_rows", 0),
                source_status=row.get("source_status", ""),
                variants=str(variants or "")[:120].replace("|", "/"),
            )
        )
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recall-score", help="Directory containing recall_score/")
    parser.add_argument("--run-dir", help="Validation run directory")
    parser.add_argument("--db-dir", help="Directory containing {GENE}.db files")
    parser.add_argument("--gold-dir", help="Directory containing *_recall_input.csv")
    parser.add_argument(
        "--results-dir", help="Run results directory with FULL_CONTEXT files"
    )
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument("--selected-out", help="Write selected regression cases CSV")
    parser.add_argument(
        "--markdown-out", help="Write selected regression cases markdown"
    )
    parser.add_argument("--limit", type=int, default=36)
    parser.add_argument("--min-per-gene", type=int, default=6)
    args = parser.parse_args()

    run_dir = resolve_run_dir(args.run_dir) if args.run_dir else None
    recall_score = resolve_recall_score(args.recall_score, run_dir)
    rows = build_paper_disagreement_rows(
        recall_score=recall_score,
        db_dir=repo_path(args.db_dir) if args.db_dir else None,
        gold_dir=repo_path(args.gold_dir) if args.gold_dir else None,
        results_dir=repo_path(args.results_dir) if args.results_dir else None,
    )
    write_csv_rows(rows, REPORT_FIELDS, repo_path(args.out))

    selected = select_regression_cases(
        rows, limit=args.limit, min_per_gene=args.min_per_gene
    )
    if args.selected_out:
        write_csv_rows(selected, REPORT_FIELDS, repo_path(args.selected_out))
    if args.markdown_out:
        write_markdown(
            selected,
            repo_path(args.markdown_out),
            title="Available-Source Extraction Regression Cases",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
