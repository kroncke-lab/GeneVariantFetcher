#!/usr/bin/env python3
"""Run priority-gated extraction from an existing GVF run directory."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.steps import (
    aggregate_data,
    extract_variants,
    migrate_to_sqlite,
    preprocess_papers,
)
from utils.bootstrap import initialize_runtime
from utils.logging_utils import setup_logging

logger = logging.getLogger("gvf.priority_extraction")


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
    except Exception:
        return None, []
    disease = context.get("prompt_disease") or context.get("requested_disease")
    terms = context.get("disease_terms") or []
    return disease, [str(term) for term in terms if str(term).strip()]


def _update_manifest(run_dir: Path, **updates: Any) -> None:
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        return
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        logger.warning("Could not read run manifest %s: %s", manifest_path, exc)
        return

    manifest["updated_at"] = datetime.now().isoformat()
    if status := updates.get("status"):
        manifest["status"] = status
    config = manifest.setdefault("config", {})
    config.update(updates.get("config", {}))
    stats = manifest.setdefault("statistics", {})
    stats.update(updates.get("statistics", {}))
    locations = manifest.setdefault("output_locations", {})
    locations.update(updates.get("output_locations", {}))
    warnings = manifest.setdefault("warnings", [])
    for warning in updates.get("warnings", []):
        warnings.append({"timestamp": datetime.now().isoformat(), "message": warning})

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Priority-gated extraction from an existing GVF run directory."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--run-dir", required=True, type=Path, help="Existing run dir")
    parser.add_argument("--disease", default=None, help="Disease prompt override")
    parser.add_argument(
        "--pmid-file",
        type=Path,
        default=None,
        help=(
            "Optional PMID list to use instead of the run's filtered_pmids.txt. "
            "Useful for targeted retries of failed/missing extraction outputs."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        required=True,
        help="Number of ranked candidates to submit to LLM extraction",
    )
    parser.add_argument(
        "--priority-offset",
        type=int,
        default=0,
        help="Skip this many ranked candidates before selecting --top-n.",
    )
    parser.add_argument(
        "--tier-threshold",
        type=int,
        default=1,
        help="Model cascade threshold passed to extraction",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Optional extraction worker count; default is provider-aware",
    )
    parser.add_argument(
        "--skip-preprocess",
        action="store_true",
        help="Skip deterministic preprocessing before extraction",
    )
    parser.add_argument(
        "--no-migrate",
        action="store_true",
        help="Skip aggregation and SQLite migration after extraction",
    )
    parser.add_argument(
        "--triage-mode",
        choices=["deterministic", "hybrid", "llm"],
        default=None,
        help="Optional cheap triage pass after priority ranking before extraction.",
    )
    parser.add_argument(
        "--triage-model",
        default=None,
        help="Model for LLM triage; defaults to the Tier-2 model.",
    )
    parser.add_argument(
        "--triage-include-defer",
        action="store_true",
        help="Also extract triage 'defer' papers.",
    )
    parser.add_argument(
        "--triage-max-llm",
        type=int,
        default=None,
        help="Optional cap on LLM triage calls.",
    )
    parser.add_argument(
        "--no-retry-failed-extractions",
        action="store_true",
        help="Disable the default retry pass for failed extraction calls.",
    )
    parser.add_argument(
        "--retry-attempts",
        type=int,
        default=1,
        help="Number of retry passes for failed extraction calls.",
    )
    parser.add_argument(
        "--retry-workers",
        type=int,
        default=1,
        help="Worker count for extraction retry passes.",
    )
    parser.add_argument(
        "--retry-backoff-seconds",
        type=float,
        default=30.0,
        help="Cooldown before each failed-extraction retry pass.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    initialize_runtime()

    gene = args.gene.upper()
    run_dir = args.run_dir.expanduser().resolve()
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    priority_dir = run_dir / "extraction_priority"
    triage_dir = run_dir / "extraction_triage"
    log_file = run_dir / f"{gene}_priority_extraction.log"

    setup_logging(level=logging.INFO, log_file=log_file)

    if not run_dir.is_dir():
        logger.error("Run directory not found: %s", run_dir)
        return 2
    if not harvest_dir.is_dir():
        logger.error("Full-text directory not found: %s", harvest_dir)
        return 2
    if args.top_n <= 0:
        logger.error("--top-n must be positive")
        return 2

    context_disease, disease_terms = _load_disease_context(run_dir)
    disease = args.disease or context_disease
    if disease and disease not in disease_terms:
        disease_terms = [disease, *disease_terms]

    records = _abstract_records(run_dir)
    pmid_source = args.pmid_file or (run_dir / "pmid_status" / "filtered_pmids.txt")
    filtered_pmids = _read_pmids(pmid_source)
    if not filtered_pmids:
        filtered_pmids = sorted(records)
    abstract_only_pmids = [pmid for pmid in filtered_pmids if pmid in records]

    initial_locations = {"extraction_priority_dir": str(priority_dir)}
    if args.triage_mode:
        initial_locations["extraction_triage_dir"] = str(triage_dir)

    _update_manifest(
        run_dir,
        status="priority_extracting_variants",
        config={
            "extraction_top_n": args.top_n,
            "extraction_priority_offset": args.priority_offset,
            "priority_extraction": True,
            "extraction_triage_mode": args.triage_mode,
            "extraction_triage_model": args.triage_model,
            "extraction_triage_include_defer": args.triage_include_defer,
            "extraction_triage_max_llm": args.triage_max_llm,
            "extraction_retry_failed_extractions": not args.no_retry_failed_extractions,
            "extraction_retry_attempts": args.retry_attempts,
            "extraction_retry_workers": args.retry_workers,
            "extraction_retry_backoff_seconds": args.retry_backoff_seconds,
            "priority_pmid_file": str(args.pmid_file) if args.pmid_file else None,
        },
        output_locations=initial_locations,
    )

    logger.info("Priority extraction for %s in %s", gene, run_dir)
    logger.info("Candidate cap: top %s", args.top_n)
    logger.info("Priority offset: %s", args.priority_offset)
    if args.triage_mode:
        logger.info(
            "Cheap triage enabled: mode=%s model=%s include_defer=%s max_llm=%s",
            args.triage_mode,
            args.triage_model or "(tier-2 default)",
            args.triage_include_defer,
            args.triage_max_llm,
        )
    logger.info(
        "Disease terms: %s", "; ".join(disease_terms) if disease_terms else "(none)"
    )
    logger.info("Abstract records available: %s", len(records))
    logger.info("PMID surface: %s from %s", len(filtered_pmids), pmid_source)

    if not args.skip_preprocess:
        preprocess = preprocess_papers(harvest_dir=harvest_dir, gene_symbol=gene)
        logger.info("Preprocess stats: %s", preprocess.stats)

    extract_result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol=gene,
        disease=disease,
        abstract_records=records,
        abstract_only_pmids=abstract_only_pmids,
        candidate_pmids=filtered_pmids,
        tier_threshold=args.tier_threshold,
        max_workers=args.max_workers,
        priority_top_n=args.top_n,
        priority_offset=args.priority_offset,
        priority_report_dir=priority_dir,
        priority_disease_terms=disease_terms,
        triage_mode=args.triage_mode,
        triage_model=args.triage_model,
        triage_report_dir=triage_dir,
        triage_include_defer=args.triage_include_defer,
        triage_max_llm_candidates=args.triage_max_llm,
        retry_failed_extractions=not args.no_retry_failed_extractions,
        extraction_retry_attempts=args.retry_attempts,
        extraction_retry_max_workers=args.retry_workers,
        extraction_retry_backoff_seconds=args.retry_backoff_seconds,
    )
    if not extract_result.success:
        _update_manifest(
            run_dir,
            status="priority_extraction_failed",
            warnings=[f"Priority extraction failed: {extract_result.error}"],
        )
        logger.error("Extraction failed: %s", extract_result.error)
        return 3

    output_locations = {
        "extractions_dir": str(extraction_dir),
        "extraction_priority_dir": str(priority_dir),
    }
    if extract_result.data.get("triage_report_dir"):
        output_locations["extraction_triage_dir"] = extract_result.data[
            "triage_report_dir"
        ]
    if extract_result.data.get("retry_report_path"):
        output_locations["extraction_retry_summary"] = extract_result.data[
            "retry_report_path"
        ]
    if extract_result.data.get("dense_table_overflow_report"):
        output_locations["dense_table_overflow"] = extract_result.data[
            "dense_table_overflow_report"
        ]
    stats = {
        "papers_extracted": extract_result.stats.get("papers_extracted", 0),
        "papers_extraction_failed": extract_result.stats.get("failures", 0),
        "papers_extraction_initial_failures": extract_result.stats.get(
            "initial_failures", 0
        ),
        "extraction_retry_attempted": extract_result.stats.get(
            "extraction_retry_attempted", 0
        ),
        "extraction_retry_succeeded": extract_result.stats.get(
            "extraction_retry_succeeded", 0
        ),
        "extraction_retry_failed": extract_result.stats.get(
            "extraction_retry_failed", 0
        ),
        "extraction_retry_skipped": extract_result.stats.get(
            "extraction_retry_skipped", 0
        ),
        "dense_table_overflow_records": extract_result.stats.get(
            "dense_table_overflow_records", 0
        ),
        "scanner_cap_trips": extract_result.stats.get("scanner_cap_trips", 0),
        "table_merge_cap_trips": extract_result.stats.get("table_merge_cap_trips", 0),
        "table_candidates_omitted_after_dedupe": extract_result.stats.get(
            "table_candidates_omitted_after_dedupe", 0
        ),
        "missing_supplement_ref_pmids": extract_result.stats.get(
            "missing_supplement_ref_pmids", 0
        ),
        "total_variants_found": extract_result.stats.get("total_variants", 0),
        "candidate_pmid_count": extract_result.stats.get("candidate_pmid_count"),
        "extraction_priority_candidates": extract_result.stats.get(
            "priority_candidates", 0
        ),
        "extraction_priority_selected": extract_result.stats.get(
            "priority_selected", 0
        ),
        "extraction_priority_offset": extract_result.stats.get("priority_offset", 0),
        "extraction_priority_fulltext_selected": extract_result.stats.get(
            "priority_fulltext_selected", 0
        ),
        "extraction_priority_abstract_selected": extract_result.stats.get(
            "priority_abstract_selected", 0
        ),
        "extraction_triage_total": extract_result.stats.get("triage_total", 0),
        "extraction_triage_extract_now": extract_result.stats.get(
            "triage_extract_now", 0
        ),
        "extraction_triage_defer": extract_result.stats.get("triage_defer", 0),
        "extraction_triage_skip": extract_result.stats.get("triage_skip", 0),
        "extraction_triage_selected": extract_result.stats.get(
            "triage_selected_for_extraction", 0
        ),
    }

    if not args.no_migrate:
        aggregate_result = aggregate_data(
            extraction_dir=extraction_dir,
            gene_symbol=gene,
            output_path=run_dir,
        )
        stats["variants_with_penetrance_data"] = aggregate_result.stats.get(
            "variants_aggregated", 0
        )
        output_locations["penetrance_summary"] = str(
            run_dir / f"{gene}_penetrance_summary.json"
        )

        db_path = run_dir / f"{gene}.db"
        migrate_result = migrate_to_sqlite(
            extraction_dir=extraction_dir, db_path=db_path
        )
        stats["migration_successful"] = migrate_result.stats.get("successful", 0)
        stats["migration_failed"] = migrate_result.stats.get("failed", 0)
        output_locations["sqlite_database"] = str(db_path)

    summary_path = priority_dir / "priority_extraction_summary.json"
    priority_dir.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(
            {
                "gene": gene,
                "run_dir": str(run_dir),
                "top_n": args.top_n,
                "priority_offset": args.priority_offset,
                "triage_mode": args.triage_mode,
                "triage_model": args.triage_model,
                "triage_include_defer": args.triage_include_defer,
                "pmid_file": str(args.pmid_file) if args.pmid_file else None,
                "disease": disease,
                "disease_terms": disease_terms,
                "stats": stats,
                "completed_at": datetime.now().isoformat(),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    output_locations["priority_extraction_summary"] = str(summary_path)

    _update_manifest(
        run_dir,
        status="priority_extraction_complete",
        statistics=stats,
        output_locations=output_locations,
    )
    logger.info("Priority extraction complete: %s", stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
