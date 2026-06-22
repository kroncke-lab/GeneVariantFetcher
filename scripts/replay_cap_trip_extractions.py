#!/usr/bin/env python3
"""Replay dense-table/scanner cap-trip PMIDs in an existing GVF run.

This is the operational path for extractor/parser fixes where an existing JSON
would otherwise be reused because the source fingerprint has not changed. It can
optionally try source recovery for PMIDs with missing supplement refs, folds any
downloaded supplements into FULL_CONTEXT, then force-reruns extraction for just
the selected PMIDs through the normal pipeline step so run-level cap QC,
failed-extraction retry, aggregation, and SQLite migration stay in sync.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from harvesting.supplement_fold import fold_supplements_into_full_context
from pipeline.steps import (
    _write_dense_table_overflow_report,
    aggregate_data,
    extract_variants,
    migrate_to_sqlite,
)
from utils.bootstrap import initialize_runtime
from utils.logging_utils import setup_logging

logger = logging.getLogger("gvf.replay_cap_trip_extractions")


def _read_pmids(path: Path) -> list[str]:
    pmids: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        token = line.split("#", 1)[0].strip()
        if token.isdigit() and token not in seen:
            seen.add(token)
            pmids.append(token)
    return pmids


def _read_cap_qc(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not path.exists():
        return {}, []
    data = json.loads(path.read_text(encoding="utf-8"))
    return data.get("summary") or {}, data.get("records") or []


def _as_int(value: Any) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _select_cap_pmids(records: list[dict[str, Any]], mode: str) -> list[str]:
    selected: list[str] = []
    seen: set[str] = set()
    for row in records:
        pmid = str(row.get("pmid") or "").strip()
        if not pmid or pmid in seen:
            continue
        scanner = bool(row.get("scanner_cap_tripped"))
        table = bool(row.get("table_merge_cap_tripped"))
        if (
            (mode == "table" and table)
            or (mode == "scanner" and scanner)
            or (mode == "any" and (table or scanner))
        ):
            seen.add(pmid)
            selected.append(pmid)
    return selected


def _pmids_with_missing_supplements(
    pmids: list[str], records: list[dict[str, Any]], harvest_dir: Path
) -> list[str]:
    by_pmid = {str(row.get("pmid") or ""): row for row in records}
    missing: list[str] = []
    for pmid in pmids:
        row = by_pmid.get(pmid) or {}
        if not row:
            artifact = harvest_dir / f"{pmid}_artifacts.json"
            if artifact.exists():
                try:
                    data = json.loads(artifact.read_text(encoding="utf-8"))
                    main_text = data.get("main_text") or {}
                    summary = data.get("summary") or {}
                    row = {
                        "supplement_refs": main_text.get(
                            "supplement_descriptions_count", 0
                        ),
                        "supplements_downloaded": summary.get(
                            "supplement_count", len(data.get("supplements") or [])
                        ),
                    }
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "Could not read artifact supplement summary for PMID %s: %s",
                        pmid,
                        exc,
                    )
        refs = _as_int(row.get("supplement_refs"))
        downloaded = _as_int(row.get("supplements_downloaded"))
        if refs > downloaded:
            missing.append(pmid)
    return missing


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
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not read disease context %s: %s", context_path, exc)
        return None, []
    disease = context.get("prompt_disease") or context.get("requested_disease")
    terms = context.get("disease_terms") or []
    return disease, [str(term) for term in terms if str(term).strip()]


def _write_fetch_input(path: Path, pmids: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["PMID"])
        writer.writeheader()
        for pmid in pmids:
            writer.writerow({"PMID": pmid})


def _run_supplement_fetch(args: argparse.Namespace, fetch_input: Path) -> int:
    command = [
        sys.executable,
        str(PROJECT_ROOT / "scripts" / "fetch_paywalled.py"),
        "--input",
        str(fetch_input),
        "--output",
        str(args.run_dir / "pmc_fulltext"),
        "--timeout-s",
        str(args.fetch_timeout_s),
        "--cookie-timeout-s",
        str(args.fetch_cookie_timeout_s),
    ]
    if args.fetch_no_cookies:
        command.append("--no-cookies")
    if args.fetch_headful:
        command.append("--no-headless")
    if args.fetch_no_chrome_channel:
        command.append("--no-chrome-channel")
    if args.browser_profile_dir:
        command.extend(["--browser-profile-dir", str(args.browser_profile_dir)])
    for allowed in args.fetch_allow or []:
        command.extend(["--allow", allowed])

    logger.info("Running supplement fetch: %s", " ".join(command))
    completed = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
    if completed.returncode:
        logger.warning(
            "Supplement fetch exited with %s; continuing with any recovered sources",
            completed.returncode,
        )
    return completed.returncode


def _fold_target_supplements(harvest_dir: Path, pmids: list[str]) -> list[str]:
    folded: list[str] = []
    for pmid in pmids:
        try:
            if fold_supplements_into_full_context(pmid, harvest_dir) is not None:
                folded.append(pmid)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Supplement fold failed for PMID %s: %s", pmid, exc)
    return folded


def _load_all_extraction_results(extraction_dir: Path, gene: str) -> list[Any]:
    from utils.models import ExtractionResult

    results: list[Any] = []
    for path in sorted(extraction_dir.glob(f"{gene}_PMID_*.json")):
        pmid = path.stem.rsplit("_PMID_", 1)[-1]
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:  # noqa: BLE001
            logger.warning("Could not read extraction JSON %s: %s", path, exc)
            continue
        metadata = data.get("extraction_metadata") if isinstance(data, dict) else {}
        results.append(
            ExtractionResult(
                pmid=pmid,
                success=True,
                extracted_data=data,
                model_used=(metadata or {}).get("model_used"),
            )
        )
    return results


def _refresh_run_level_cap_qc(
    *,
    run_dir: Path,
    harvest_dir: Path,
    extraction_dir: Path,
    gene: str,
) -> tuple[Path | None, Path | None, dict[str, Any]]:
    return _write_dense_table_overflow_report(
        extractions=_load_all_extraction_results(extraction_dir, gene),
        harvest_dir=harvest_dir,
        output_dir=run_dir,
    )


def _update_manifest(run_dir: Path, summary: dict[str, Any]) -> None:
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        return
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not read run manifest %s: %s", manifest_path, exc)
        return

    manifest["updated_at"] = datetime.now().isoformat()
    stats = manifest.setdefault("statistics", {})
    stats.update(
        {
            "cap_trip_replay_input_pmids": summary["input_pmids"],
            "cap_trip_replay_extracted_pmids": summary["papers_extracted"],
            "cap_trip_replay_failures": summary["failures"],
            "cap_trip_replay_last_completed_at": summary["completed_at"],
        }
    )
    for key in (
        "dense_table_overflow_records",
        "scanner_cap_trips",
        "table_merge_cap_trips",
        "table_candidates_omitted_after_dedupe",
        "missing_supplement_ref_pmids",
        "migration_successful",
        "migration_failed",
    ):
        if key in summary:
            stats[key] = summary[key]

    locations = manifest.setdefault("output_locations", {})
    locations["cap_trip_replay_summary"] = summary["summary_path"]
    for key in (
        "cap_trip_replay_pmid_list",
        "cap_trip_replay_fetch_input",
        "dense_table_overflow",
        "dense_table_overflow_tsv",
        "extraction_retry_summary",
        "sqlite_database",
    ):
        if summary.get(key):
            locations[key] = summary[key]

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--run-dir", required=True, type=Path, help="Existing run dir")
    parser.add_argument(
        "--pmid-file",
        type=Path,
        default=None,
        help="Explicit PMID list. Defaults to PMIDs from --cap-qc.",
    )
    parser.add_argument(
        "--cap-qc",
        type=Path,
        default=None,
        help="dense_table_overflow.json. Defaults to <run-dir>/pmid_status/.",
    )
    parser.add_argument(
        "--mode",
        choices=["table", "scanner", "any"],
        default="table",
        help="Which cap-trip records to replay when --pmid-file is omitted.",
    )
    parser.add_argument(
        "--fetch-supplements",
        action="store_true",
        help="Run fetch_paywalled.py first for selected PMIDs with missing supplement refs.",
    )
    parser.add_argument(
        "--fetch-all-input",
        action="store_true",
        help="With --fetch-supplements, fetch every input PMID, not just missing-supplement rows.",
    )
    parser.add_argument("--fetch-timeout-s", type=int, default=120)
    parser.add_argument("--fetch-cookie-timeout-s", type=float, default=8.0)
    parser.add_argument("--fetch-allow", action="append", default=None)
    parser.add_argument("--fetch-no-cookies", action="store_true")
    parser.add_argument("--fetch-headful", action="store_true")
    parser.add_argument("--fetch-no-chrome-channel", action="store_true")
    parser.add_argument("--browser-profile-dir", type=Path, default=None)
    parser.add_argument(
        "--no-fold-supplements",
        action="store_true",
        help="Skip folding on-disk supplements into FULL_CONTEXT before replay.",
    )
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--retry-attempts", type=int, default=1)
    parser.add_argument("--retry-workers", type=int, default=1)
    parser.add_argument("--retry-backoff-seconds", type=float, default=30.0)
    parser.add_argument("--no-migrate", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    initialize_runtime()

    gene = args.gene.upper()
    args.run_dir = args.run_dir.expanduser().resolve()
    run_dir = args.run_dir
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    pmid_status_dir = run_dir / "pmid_status"
    pmid_status_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(
        level=logging.INFO,
        log_file=run_dir / f"{gene}_cap_trip_replay.log",
    )

    if not run_dir.is_dir():
        logger.error("Run directory not found: %s", run_dir)
        return 2
    if not harvest_dir.is_dir():
        logger.error("Full-text directory not found: %s", harvest_dir)
        return 2

    cap_qc_path = args.cap_qc or (pmid_status_dir / "dense_table_overflow.json")
    cap_qc_summary, cap_records = _read_cap_qc(cap_qc_path)
    input_pmids = (
        _read_pmids(args.pmid_file.expanduser())
        if args.pmid_file
        else _select_cap_pmids(cap_records, args.mode)
    )
    if not input_pmids:
        logger.info("No cap-trip PMIDs selected for %s", gene)
        return 0

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pmid_list_path = pmid_status_dir / f"cap_trip_replay_pmids_{timestamp}.txt"
    pmid_list_path.write_text("\n".join(input_pmids) + "\n", encoding="utf-8")

    fetch_returncode: int | None = None
    fetch_input_path: Path | None = None
    fetch_pmids: list[str] = []
    if args.fetch_supplements:
        fetch_pmids = (
            input_pmids
            if args.fetch_all_input
            else _pmids_with_missing_supplements(input_pmids, cap_records, harvest_dir)
        )
        if fetch_pmids:
            fetch_input_path = pmid_status_dir / f"cap_trip_fetch_input_{timestamp}.csv"
            _write_fetch_input(fetch_input_path, fetch_pmids)
            fetch_returncode = _run_supplement_fetch(args, fetch_input_path)
        else:
            logger.info("No selected cap-trip PMIDs have missing supplement refs")

    folded_pmids: list[str] = []
    if not args.no_fold_supplements:
        folded_pmids = _fold_target_supplements(harvest_dir, input_pmids)

    disease, _disease_terms = _load_disease_context(run_dir)
    records = _abstract_records(run_dir)
    abstract_only_pmids = [pmid for pmid in input_pmids if pmid in records]

    logger.info(
        "Force-replaying %d %s cap-trip PMIDs for %s in %s",
        len(input_pmids),
        args.mode,
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
        force_pmids=input_pmids,
        max_workers=args.max_workers,
        retry_failed_extractions=True,
        extraction_retry_attempts=args.retry_attempts,
        extraction_retry_max_workers=args.retry_workers,
        extraction_retry_backoff_seconds=args.retry_backoff_seconds,
    )
    if not extract_result.success:
        logger.error("Cap-trip replay failed: %s", extract_result.error)
        return 3

    dense_report_path, dense_tsv_path, dense_summary = _refresh_run_level_cap_qc(
        run_dir=run_dir,
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene=gene,
    )

    summary: dict[str, Any] = {
        "gene": gene,
        "run_dir": str(run_dir),
        "completed_at": datetime.now().isoformat(),
        "mode": args.mode,
        "input_pmids": len(input_pmids),
        "input_pmid_list": input_pmids,
        "cap_qc_input": str(cap_qc_path),
        "cap_qc_input_summary": cap_qc_summary,
        "cap_trip_replay_pmid_list": str(pmid_list_path),
        "fetch_requested": args.fetch_supplements,
        "fetch_pmids": fetch_pmids,
        "fetch_returncode": fetch_returncode,
        "folded_supplement_pmids": folded_pmids,
        "papers_extracted": extract_result.stats.get("papers_extracted", 0),
        "failures": extract_result.stats.get("failures", 0),
        "initial_failures": extract_result.stats.get("initial_failures", 0),
        "total_variants": extract_result.stats.get("total_variants", 0),
        "extract_stats": extract_result.stats,
        "summary_path": "",
    }
    if fetch_input_path:
        summary["cap_trip_replay_fetch_input"] = str(fetch_input_path)
    if extract_result.data.get("retry_report_path"):
        summary["extraction_retry_summary"] = extract_result.data["retry_report_path"]
    if dense_report_path:
        summary["dense_table_overflow"] = str(dense_report_path)
    if dense_tsv_path:
        summary["dense_table_overflow_tsv"] = str(dense_tsv_path)
    for key in (
        "dense_table_overflow_records",
        "scanner_cap_trips",
        "table_merge_cap_trips",
        "table_candidates_omitted_after_dedupe",
        "missing_supplement_ref_pmids",
    ):
        summary[key] = dense_summary.get(
            {
                "dense_table_overflow_records": "records",
                "scanner_cap_trips": "scanner_cap_trips",
                "table_merge_cap_trips": "table_merge_cap_trips",
                "table_candidates_omitted_after_dedupe": "table_candidates_omitted_after_dedupe",
                "missing_supplement_ref_pmids": "missing_supplement_ref_pmids",
            }[key],
            0,
        )

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

    summary_path = pmid_status_dir / f"cap_trip_replay_summary_{timestamp}.json"
    summary["summary_path"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _update_manifest(run_dir, summary)

    logger.info("Cap-trip replay summary: %s", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
