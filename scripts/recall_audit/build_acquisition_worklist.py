#!/usr/bin/env python3
"""Build a replay/refetch worklist from paper-level recall disagreements.

This is an audit helper, not a production pipeline stage. It turns
``paper_disagreement_report.csv`` plus the row-level ``missing_in_sqlite.csv``
artifact into concrete next actions:

* replay an already-available on-disk source with ``refresh_run_db.py``;
* refetch a missing/stub source with ``fetch_paywalled.py``;
* route blocked publisher cases to manual acquisition; or
* leave true extraction gaps for parser/vision work.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

try:
    from scripts.recall_audit.common import (
        canonical_variant,
        normalize_pmid,
        parse_bool,
        parse_int,
        read_csv_rows,
        repo_path,
        source_has_unresolved_variant_supplement_refs,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        canonical_variant,
        normalize_pmid,
        parse_bool,
        parse_int,
        read_csv_rows,
        repo_path,
        source_has_unresolved_variant_supplement_refs,
        write_csv_rows,
    )


DATA_AVAILABLE_STATUSES = {
    "recovered_pmc",
    "recovered_browser",
    "recovered_supplement_only",
}

SOURCE_PRESENT_STATUSES = DATA_AVAILABLE_STATUSES | {"missing_table_bodies"}
MISSING_SOURCE_STATUSES = {"abstract_only", "not_attempted", "paywall_stub"}
FETCH_ACTIONS = {"fetch"}
REPLAY_ACTIONS = {"refresh_replay"}

WORKLIST_FIELDS = [
    "gene",
    "pmid",
    "action",
    "route",
    "priority_score",
    "missing_rows",
    "missing_distinct_variants",
    "missing_carriers",
    "failure_class",
    "source_status",
    "source_file",
    "context_path",
    "context_bytes",
    "available_source_status",
    "available_context_path",
    "available_context_bytes",
    "source_desync",
    "source_unbound",
    "likely_missing_variant_supplement",
    "doi",
    "doi_source",
    "publisher_hint",
    "top_missing_variants",
    "notes",
]


DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>]+", re.IGNORECASE)


def _clean_doi(value: str) -> str:
    doi = unquote(value.strip())
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi, flags=re.IGNORECASE)
    doi = doi.strip().strip(".,;:)]}>")
    return doi


def _extract_doi(value: Any) -> str:
    if value is None:
        return ""
    text = unquote(str(value))
    match = DOI_RE.search(text)
    if not match:
        return ""
    return _clean_doi(match.group(0))


def _column(fieldnames: list[str], *candidates: str) -> str | None:
    lower_to_field = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in lower_to_field:
            return lower_to_field[candidate]
    return None


def load_doi_map(paths: list[Path]) -> dict[str, tuple[str, str]]:
    """Load PMID -> (DOI, source-path) from heterogeneous local CSV artifacts."""
    doi_by_pmid: dict[str, tuple[str, str]] = {}
    for path in paths:
        path = repo_path(path).expanduser()
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            pmid_col = _column(reader.fieldnames, "pmid", "pubmed_id", "pubmed", "id")
            if pmid_col is None:
                continue
            doi_col = _column(reader.fieldnames, "doi", "article_doi")
            scan_cols = [
                field
                for field in ("URL", "url", "Reason", "reason", "final_url", "path")
                if field in reader.fieldnames
            ]
            for row in reader:
                pmid = normalize_pmid(row.get(pmid_col))
                if not pmid or pmid in doi_by_pmid:
                    continue
                doi = _clean_doi(row.get(doi_col) or "") if doi_col else ""
                if not doi:
                    for field in scan_cols:
                        doi = _extract_doi(row.get(field))
                        if doi:
                            break
                if doi:
                    doi_by_pmid[pmid] = (doi, str(path))
    return doi_by_pmid


def load_missing_distinct_counts(path: Path | None) -> dict[str, int]:
    if path is None or not path.exists():
        return {}
    variants_by_pmid: dict[str, set[str]] = defaultdict(set)
    for row in read_csv_rows(path):
        pmid = normalize_pmid(row.get("pmid"))
        variant = canonical_variant(
            row.get("excel_variant_raw") or row.get("variant") or row.get("variant_raw")
        )
        if pmid and variant:
            variants_by_pmid[pmid].add(variant)
    return {pmid: len(variants) for pmid, variants in variants_by_pmid.items()}


def _doi_from_context_artifacts(pmid: str, *paths: str) -> tuple[str, str]:
    for raw_path in paths:
        if not raw_path:
            continue
        path = repo_path(raw_path).expanduser()
        parent = path.parent
        candidates = [
            parent / pmid / "result.json",
            parent / f"{pmid}_artifacts.json",
        ]
        for candidate in candidates:
            if not candidate.exists():
                continue
            try:
                payload = json.loads(candidate.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(payload, dict):
                continue
            doi = _clean_doi(str(payload.get("doi") or ""))
            if not doi and isinstance(payload.get("main_text"), dict):
                doi = _extract_doi(payload["main_text"])
            if doi:
                return doi, str(candidate)
    return "", ""


def default_missing_path(report: Path, gene: str) -> Path:
    return report.parent / gene.upper() / "missing_in_sqlite.csv"


def _existing_same_path(a: str, b: str) -> bool:
    if not a or not b:
        return False
    left = Path(a).expanduser()
    right = Path(b).expanduser()
    try:
        if left.exists() and right.exists():
            return left.resolve() == right.resolve()
    except OSError:
        pass
    return str(left) == str(right)


def infer_publisher(doi: str) -> str:
    lower = doi.lower()
    if lower.startswith("10.1016/"):
        return "elsevier"
    if lower.startswith(("10.1002/", "10.1111/")):
        return "wiley"
    if lower.startswith("10.1007/"):
        return "springer"
    if lower.startswith("10.1038/"):
        return "springer_nature"
    if lower.startswith("10.1089/"):
        return "liebert_sage"
    if lower.startswith("10.1161/"):
        return "aha"
    return "unknown"


def fetch_route_for_doi(doi: str) -> tuple[str, str, str]:
    publisher = infer_publisher(doi)
    if publisher == "elsevier":
        return "fetch", "fetch_elsevier_insttoken", publisher
    if publisher == "wiley":
        return "fetch", "fetch_wiley_tdm", publisher
    if publisher == "springer":
        return "fetch", "fetch_springer_api", publisher
    if publisher == "springer_nature":
        return "fetch", "fetch_springer_nature_or_browser", publisher
    if publisher == "liebert_sage":
        return "manual_or_blocked", "blocked_liebert_sage_cloudflare", publisher
    if publisher == "aha":
        return "fetch", "fetch_browser_aha", publisher
    return "fetch", "fetch_paywalled_generic", publisher


def classify_row(
    row: dict[str, str],
    *,
    gene: str,
    missing_distinct_variants: int,
    doi: str,
    doi_source: str,
) -> dict[str, Any]:
    source_status = row.get("source_status") or ""
    available_status = row.get("available_source_status") or ""
    context_path = row.get("context_path") or ""
    available_path = row.get("available_context_path") or ""
    context_bytes = parse_int(row.get("context_bytes")) or 0
    available_bytes = parse_int(row.get("available_context_bytes")) or 0
    source_desync = parse_bool(row.get("source_desync"))
    source_unbound = parse_bool(row.get("source_unbound"))
    failure_class = row.get("failure_class") or ""
    missing_rows = parse_int(row.get("missing_rows")) or 0
    missing_variant_supplement = parse_bool(
        row.get("likely_missing_variant_supplement")
    )
    source_path_for_supplement_check = available_path or context_path
    if not missing_variant_supplement and source_path_for_supplement_check:
        missing_variant_supplement = source_has_unresolved_variant_supplement_refs(
            repo_path(source_path_for_supplement_check).expanduser(),
            gene,
        )

    available_replayable = (
        available_status in DATA_AVAILABLE_STATUSES
        and bool(available_path)
        and (
            source_status not in SOURCE_PRESENT_STATUSES
            or source_desync
            or source_unbound
            or (
                not _existing_same_path(context_path, available_path)
                and available_bytes > context_bytes + 2048
            )
        )
    )

    notes: list[str] = []
    if available_replayable:
        action = "refresh_replay"
        route = "replay_available_context"
        publisher_hint = infer_publisher(doi) if doi else ""
        notes.append(
            "usable source exists on disk but current extraction is not bound to it"
        )
    elif source_status in SOURCE_PRESENT_STATUSES:
        if missing_variant_supplement:
            if doi:
                action, route, publisher_hint = fetch_route_for_doi(doi)
            else:
                action = "fetch"
                route = "doi_lookup_then_fetch"
                publisher_hint = ""
            notes.append(
                "usable main text references a missing target-gene variant supplement"
            )
        else:
            action = "extractor_or_vision"
            publisher_hint = infer_publisher(doi) if doi else ""
        if (
            not missing_variant_supplement
            and failure_class == "source_missing_table_bodies"
        ):
            route = "missing_table_body_or_supplement"
        elif failure_class == "available_source_underextraction":
            route = (
                route
                if missing_variant_supplement
                else "available_source_underextraction"
            )
        elif not missing_variant_supplement:
            route = "available_source_extraction_gap"
    elif source_status in MISSING_SOURCE_STATUSES or failure_class.startswith(
        "source_"
    ):
        if doi:
            action, route, publisher_hint = fetch_route_for_doi(doi)
            notes.append("source is missing, abstract-only, or a paywall stub")
        else:
            action = "fetch"
            route = "doi_lookup_then_fetch"
            publisher_hint = ""
            notes.append("source is missing/stub and no DOI was found locally")
    else:
        action = "inspect"
        route = "unclassified_gap"
        publisher_hint = infer_publisher(doi) if doi else ""

    priority = (
        missing_rows * 10
        + missing_distinct_variants * 5
        + (25 if action == "refresh_replay" else 0)
        + (10 if action == "fetch" else 0)
        - (20 if action == "manual_or_blocked" else 0)
    )

    return {
        "gene": row.get("gene") or "",
        "pmid": normalize_pmid(row.get("pmid")),
        "action": action,
        "route": route,
        "priority_score": priority,
        "missing_rows": missing_rows,
        "missing_distinct_variants": missing_distinct_variants,
        "missing_carriers": row.get("missing_carriers") or "",
        "failure_class": failure_class,
        "source_status": source_status,
        "source_file": row.get("source_file") or "",
        "context_path": context_path,
        "context_bytes": context_bytes,
        "available_source_status": available_status,
        "available_context_path": available_path,
        "available_context_bytes": available_bytes,
        "source_desync": source_desync,
        "source_unbound": source_unbound,
        "likely_missing_variant_supplement": missing_variant_supplement,
        "doi": doi,
        "doi_source": doi_source,
        "publisher_hint": publisher_hint,
        "top_missing_variants": row.get("top_missing_variants") or "",
        "notes": "; ".join(notes),
    }


def build_worklist(
    *,
    report: Path,
    gene: str,
    missing_path: Path | None,
    doi_sources: list[Path],
    min_missing_rows: int,
) -> list[dict[str, Any]]:
    doi_by_pmid = load_doi_map(doi_sources)
    missing_distinct = load_missing_distinct_counts(missing_path)
    rows: list[dict[str, Any]] = []
    for row in read_csv_rows(report):
        if str(row.get("gene") or "").upper() != gene.upper():
            continue
        missing_rows = parse_int(row.get("missing_rows")) or 0
        if missing_rows < min_missing_rows:
            continue
        pmid = normalize_pmid(row.get("pmid"))
        doi, doi_source = doi_by_pmid.get(pmid, ("", ""))
        if not doi:
            doi, doi_source = _doi_from_context_artifacts(
                pmid,
                row.get("available_context_path") or "",
                row.get("context_path") or "",
            )
        rows.append(
            classify_row(
                row,
                gene=gene,
                missing_distinct_variants=missing_distinct.get(pmid, 0),
                doi=doi,
                doi_source=doi_source,
            )
        )
    return sorted(
        rows,
        key=lambda row: (
            row["action"] not in {"refresh_replay", "fetch"},
            row["action"] == "manual_or_blocked",
            -int(row["priority_score"]),
            row["pmid"],
        ),
    )


def write_replay_pmids(rows: list[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pmids = [row["pmid"] for row in rows if row["action"] in REPLAY_ACTIONS]
    out_path.write_text("\n".join(pmids) + ("\n" if pmids else ""), encoding="utf-8")


def write_fetch_input(rows: list[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fetch_rows = [
        {"PMID": row["pmid"], "DOI": row["doi"], "route": row["route"]}
        for row in rows
        if row["action"] in FETCH_ACTIONS
    ]
    write_csv_rows(fetch_rows, ["PMID", "DOI", "route"], out_path)


def _recall_item(count: int, total: int) -> dict[str, Any]:
    return {
        "pmids": count,
        "gold_pmids": total,
        "recall": round(count / total, 6) if total else None,
    }


def pmid_recall_summary(
    rows: list[dict[str, Any]],
    *,
    report_rows: list[dict[str, str]],
) -> dict[str, Any]:
    """Summarize source/acquisition PMID coverage against all gold PMIDs."""
    gold_pmids = {
        normalize_pmid(row.get("pmid"))
        for row in report_rows
        if normalize_pmid(row.get("pmid"))
    }
    selected_for_fetch = {row["pmid"] for row in rows if row["action"] == "fetch"}
    selected_for_source = {
        row["pmid"]
        for row in rows
        if row["action"] in {"fetch", "refresh_replay", "manual_or_blocked"}
    }
    bound_usable = {
        normalize_pmid(row.get("pmid"))
        for row in report_rows
        if parse_bool(row.get("data_available"))
    }
    available_usable = {
        normalize_pmid(row.get("pmid"))
        for row in report_rows
        if str(row.get("available_source_status") or "") in DATA_AVAILABLE_STATUSES
    }
    total = len(gold_pmids)
    return {
        "selected_for_fetch": _recall_item(len(selected_for_fetch), total),
        "selected_for_source_acquisition_or_binding": _recall_item(
            len(selected_for_source), total
        ),
        "usable_fulltext_bound_current": _recall_item(len(bound_usable), total),
        "usable_fulltext_available_current": _recall_item(len(available_usable), total),
        "notes": {
            "selected_for_fetch": "Gold PMIDs queued for download/refetch by this worklist.",
            "selected_for_source_acquisition_or_binding": (
                "Gold PMIDs queued either for download/refetch, manual blocked "
                "acquisition, or replay from an already available source."
            ),
            "usable_fulltext_bound_current": (
                "Gold PMIDs whose currently bound extraction source is usable "
                "full text in the input paper_disagreement_report."
            ),
            "usable_fulltext_available_current": (
                "Gold PMIDs with any usable full text already present on disk, "
                "even if the current extraction is not bound to it yet."
            ),
        },
    }


def write_summary(
    rows: list[dict[str, Any]],
    out_path: Path,
    *,
    report_rows: list[dict[str, str]],
) -> None:
    by_action = Counter(row["action"] for row in rows)
    by_route = Counter(row["route"] for row in rows)
    rows_by_action = Counter()
    unique_by_action = Counter()
    for row in rows:
        rows_by_action[row["action"]] += int(row["missing_rows"] or 0)
        unique_by_action[row["action"]] += int(row["missing_distinct_variants"] or 0)
    payload = {
        "pmids": len(rows),
        "missing_rows": sum(int(row["missing_rows"] or 0) for row in rows),
        "missing_distinct_variants": sum(
            int(row["missing_distinct_variants"] or 0) for row in rows
        ),
        "by_action": dict(by_action),
        "missing_rows_by_action": dict(rows_by_action),
        "missing_distinct_variants_by_action": dict(unique_by_action),
        "by_route": dict(by_route),
        "pmid_recall": pmid_recall_summary(rows, report_rows=report_rows),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument(
        "--report", required=True, type=Path, help="paper_disagreement_report.csv"
    )
    parser.add_argument(
        "--missing",
        type=Path,
        default=None,
        help="Row-level missing_in_sqlite.csv. Default: <report-parent>/<GENE>/missing_in_sqlite.csv",
    )
    parser.add_argument(
        "--doi-source",
        action="append",
        type=Path,
        default=[],
        help="CSV to scan for PMID/DOI or PMID/URL mappings. Repeatable.",
    )
    parser.add_argument("--min-missing-rows", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True, help="Worklist CSV path")
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--replay-pmids-out", type=Path, default=None)
    parser.add_argument("--fetch-input-out", type=Path, default=None)
    args = parser.parse_args()

    report = repo_path(args.report).expanduser()
    gene = args.gene.upper()
    missing = (
        repo_path(args.missing).expanduser()
        if args.missing
        else default_missing_path(report, gene)
    )
    if not report.exists():
        sys.exit(f"Report not found: {report}")
    if not missing.exists():
        sys.exit(f"Missing rows CSV not found: {missing}")

    report_rows = [
        row
        for row in read_csv_rows(report)
        if str(row.get("gene") or "").upper() == gene
    ]
    rows = build_worklist(
        report=report,
        gene=gene,
        missing_path=missing,
        doi_sources=args.doi_source,
        min_missing_rows=args.min_missing_rows,
    )
    write_csv_rows(rows, WORKLIST_FIELDS, repo_path(args.out).expanduser())
    if args.summary_out:
        write_summary(
            rows,
            repo_path(args.summary_out).expanduser(),
            report_rows=report_rows,
        )
    if args.replay_pmids_out:
        write_replay_pmids(rows, repo_path(args.replay_pmids_out).expanduser())
    if args.fetch_input_out:
        write_fetch_input(rows, repo_path(args.fetch_input_out).expanduser())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
