#!/usr/bin/env python3
"""Retry failed extraction PMIDs in an existing GVF run directory."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.steps import aggregate_data, extract_variants, migrate_to_sqlite
from utils.bootstrap import initialize_runtime
from utils.logging_utils import setup_logging

logger = logging.getLogger("gvf.retry_failed_extractions")


def _read_pmids(path: Path) -> list[str]:
    if not path.exists():
        return []
    pmids: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        token = line.split("#", 1)[0].strip()
        if token.isdigit() and token not in seen:
            seen.add(token)
            pmids.append(token)
    return pmids


def _read_failure_pmids(path: Path) -> list[str]:
    if not path.exists():
        return []
    pmids: list[str] = []
    seen: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pmid = str(row.get("PMID") or "").strip()
            if pmid.isdigit() and pmid not in seen:
                seen.add(pmid)
                pmids.append(pmid)
    return pmids


def _read_failure_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return ["PMID", "Error"], []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or ["PMID", "Error"])
        rows = [dict(row) for row in reader]
    if "PMID" not in fieldnames:
        fieldnames.insert(0, "PMID")
    if "Error" not in fieldnames:
        fieldnames.append("Error")
    return fieldnames, rows


def _write_remaining_failure_csv(
    path: Path,
    original_fieldnames: list[str],
    original_rows: list[dict[str, str]],
    remaining_failures: list[str],
    retry_errors: dict[str, str],
) -> None:
    by_pmid = {
        str(row.get("PMID") or "").strip(): row
        for row in original_rows
        if str(row.get("PMID") or "").strip()
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=original_fieldnames)
        writer.writeheader()
        for pmid in remaining_failures:
            row = dict(by_pmid.get(pmid) or {})
            row["PMID"] = pmid
            row["Error"] = retry_errors.get(pmid) or row.get("Error") or "Retry failed"
            writer.writerow(
                {field: row.get(field, "") for field in original_fieldnames}
            )


def _abstract_records(run_dir: Path) -> dict[str, str]:
    abstract_dir = run_dir / "abstract_json"
    if not abstract_dir.is_dir():
        return {}
    return {path.stem: str(path) for path in sorted(abstract_dir.glob("*.json"))}


def _load_disease_context(run_dir: Path) -> tuple[str | None, list[str]]:
    context_path = run_dir / "gene_disease_context.json"
    if not context_path.exists():
        return None, []
    try:
        context = json.loads(context_path.read_text(encoding="utf-8"))
    except Exception as exc:
        logger.warning("Could not read disease context %s: %s", context_path, exc)
        return None, []
    disease = context.get("prompt_disease") or context.get("requested_disease")
    terms = context.get("disease_terms") or []
    return disease, [str(term) for term in terms if str(term).strip()]


def _update_manifest(run_dir: Path, summary: dict[str, Any]) -> None:
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        return
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        logger.warning("Could not read run manifest %s: %s", manifest_path, exc)
        return

    manifest["updated_at"] = datetime.now().isoformat()
    stats = manifest.setdefault("statistics", {})
    stats.update(
        {
            "papers_extraction_failed": summary["remaining_failures"],
            "targeted_retry_input_pmids": summary["input_pmids"],
            "targeted_retry_recovered_pmids": summary["recovered_pmids"],
            "targeted_retry_remaining_failures": summary["remaining_failures"],
            "targeted_retry_last_completed_at": summary["completed_at"],
        }
    )
    if "migration_successful" in summary:
        stats["migration_successful"] = summary["migration_successful"]
    if "migration_failed" in summary:
        stats["migration_failed"] = summary["migration_failed"]
    for key in (
        "dense_table_overflow_records",
        "scanner_cap_trips",
        "table_merge_cap_trips",
        "table_candidates_omitted_after_dedupe",
        "missing_supplement_ref_pmids",
    ):
        if key in summary:
            stats[key] = summary[key]
    locations = manifest.setdefault("output_locations", {})
    locations["targeted_retry_summary"] = summary["summary_path"]
    if summary.get("extraction_retry_summary"):
        locations["extraction_retry_summary"] = summary["extraction_retry_summary"]
    if summary.get("dense_table_overflow"):
        locations["dense_table_overflow"] = summary["dense_table_overflow"]
    if summary.get("dense_table_overflow_tsv"):
        locations["dense_table_overflow_tsv"] = summary["dense_table_overflow_tsv"]

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Retry failed extraction PMIDs in an existing GVF run directory."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--run-dir", required=True, type=Path, help="Existing run dir")
    parser.add_argument(
        "--failure-csv",
        type=Path,
        default=None,
        help="Failure CSV to read; defaults to <run-dir>/pmid_status/extraction_failures.csv",
    )
    parser.add_argument(
        "--pmid-file",
        type=Path,
        default=None,
        help="Optional explicit PMID list instead of a failure CSV.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Worker count for the first targeted retry pass.",
    )
    parser.add_argument(
        "--retry-attempts",
        type=int,
        default=1,
        help="Additional retry passes for failures from this targeted run.",
    )
    parser.add_argument(
        "--retry-workers",
        type=int,
        default=1,
        help="Worker count for additional retry passes.",
    )
    parser.add_argument(
        "--retry-backoff-seconds",
        type=float,
        default=30.0,
        help="Cooldown before each additional retry pass.",
    )
    parser.add_argument(
        "--no-migrate",
        action="store_true",
        help="Skip aggregation and SQLite migration after retry.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    initialize_runtime()

    gene = args.gene.upper()
    run_dir = args.run_dir.expanduser().resolve()
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    pmid_status_dir = run_dir / "pmid_status"
    pmid_status_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / f"{gene}_failed_extraction_retry.log"
    setup_logging(level=logging.INFO, log_file=log_file)

    if not run_dir.is_dir():
        logger.error("Run directory not found: %s", run_dir)
        return 2
    if not harvest_dir.is_dir():
        logger.error("Full-text directory not found: %s", harvest_dir)
        return 2

    failure_csv = args.failure_csv or (pmid_status_dir / "extraction_failures.csv")
    failure_fieldnames, original_failure_rows = _read_failure_rows(failure_csv)
    input_pmids = (
        _read_pmids(args.pmid_file)
        if args.pmid_file
        else _read_failure_pmids(failure_csv)
    )
    if not input_pmids:
        logger.info("No failed extraction PMIDs to retry for %s", gene)
        return 0

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pmid_list_path = pmid_status_dir / f"failed_extraction_retry_pmids_{timestamp}.txt"
    pmid_list_path.write_text("\n".join(input_pmids) + "\n", encoding="utf-8")

    disease, _disease_terms = _load_disease_context(run_dir)
    records = _abstract_records(run_dir)
    abstract_only_pmids = [pmid for pmid in input_pmids if pmid in records]

    logger.info(
        "Retrying %s failed extraction PMIDs for %s in %s",
        len(input_pmids),
        gene,
        run_dir,
    )
    extract_result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol=gene,
        disease=disease,
        abstract_records=records,
        abstract_only_pmids=abstract_only_pmids,
        candidate_pmids=input_pmids,
        max_workers=args.max_workers,
        retry_failed_extractions=True,
        extraction_retry_attempts=args.retry_attempts,
        extraction_retry_max_workers=args.retry_workers,
        extraction_retry_backoff_seconds=args.retry_backoff_seconds,
    )
    if not extract_result.success:
        logger.error("Targeted retry failed: %s", extract_result.error)
        return 3

    remaining_failures = [
        str(pmid) for pmid, _error in extract_result.data.get("failures", [])
    ]
    retry_errors = {
        str(pmid): str(error) for pmid, error in extract_result.data.get("failures", [])
    }
    recovered_pmids = [pmid for pmid in input_pmids if pmid not in remaining_failures]
    _write_remaining_failure_csv(
        failure_csv,
        failure_fieldnames,
        original_failure_rows,
        remaining_failures,
        retry_errors,
    )

    summary: dict[str, Any] = {
        "gene": gene,
        "run_dir": str(run_dir),
        "completed_at": datetime.now().isoformat(),
        "input_pmids": len(input_pmids),
        "recovered_pmids": len(recovered_pmids),
        "remaining_failures": len(remaining_failures),
        "pmid_list": str(pmid_list_path),
        "recovered_pmid_list": recovered_pmids,
        "remaining_failure_pmids": remaining_failures,
        "failure_csv": str(failure_csv),
        "extract_stats": extract_result.stats,
        "summary_path": "",
    }
    if extract_result.data.get("retry_report_path"):
        summary["extraction_retry_summary"] = extract_result.data["retry_report_path"]
    if extract_result.data.get("dense_table_overflow_report"):
        summary["dense_table_overflow"] = extract_result.data[
            "dense_table_overflow_report"
        ]
    if extract_result.data.get("dense_table_overflow_tsv"):
        summary["dense_table_overflow_tsv"] = extract_result.data[
            "dense_table_overflow_tsv"
        ]
    for key in (
        "dense_table_overflow_records",
        "scanner_cap_trips",
        "table_merge_cap_trips",
        "table_candidates_omitted_after_dedupe",
        "missing_supplement_ref_pmids",
    ):
        summary[key] = extract_result.stats.get(key, 0)

    if not args.no_migrate:
        aggregate_result = aggregate_data(
            extraction_dir=extraction_dir,
            gene_symbol=gene,
            output_path=run_dir,
        )
        summary["variants_with_penetrance_data"] = aggregate_result.stats.get(
            "variants_aggregated", 0
        )
        db_path = run_dir / f"{gene}.db"
        migrate_result = migrate_to_sqlite(
            extraction_dir=extraction_dir,
            db_path=db_path,
        )
        summary["migration_successful"] = migrate_result.stats.get("successful", 0)
        summary["migration_failed"] = migrate_result.stats.get("failed", 0)
        summary["sqlite_database"] = str(db_path)

    summary_path = pmid_status_dir / f"failed_extraction_retry_summary_{timestamp}.json"
    summary["summary_path"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _update_manifest(run_dir, summary)

    logger.info("Failed extraction retry summary: %s", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
