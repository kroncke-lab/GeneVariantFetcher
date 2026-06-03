#!/usr/bin/env python3
"""Summarize acquisition worklist selection and actual full-text outcomes.

This bridges two artifacts:

* acquisition worklist output: what we intended to fetch/replay.
* ``fetch_paywalled.py`` ``summary.json`` or output directory: what actually
  landed as usable text.

When a gold/report denominator is supplied, the output reports PMID recall for
both the selected acquisition queue and the successfully downloaded/extracted
full-text PMIDs. Without gold, it reports the same PMIDs as worklist coverage so
new-gene runs can still distinguish selected vs actually usable acquisition
outcomes. It can also emit a ``source_override`` CSV for ``refresh_run_db.py`` so
successful fetches can be replayed from a staged extraction directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from pipeline.source_quality import is_usable_fulltext_source  # noqa: E402
    from scripts.recall_audit.common import (  # noqa: E402
        normalize_pmid,
        parse_int,
        read_csv_rows,
        repo_path,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from pipeline.source_quality import is_usable_fulltext_source  # type: ignore # noqa: E402
    from common import (  # type: ignore  # noqa: E402
        normalize_pmid,
        parse_int,
        read_csv_rows,
        repo_path,
        write_csv_rows,
    )


FETCH_SUCCESS_OUTCOMES = {
    "success",
    "success_from_output_dir",
    "success_via_pmc",
    "success_via_scholar_pdf",
    "success_via_elsevier_api",
    "success_via_publisher_api",
    "success_via_springer_api",
    "success_via_wiley_api",
    "success_supplement_only",
}

SOURCE_OVERRIDE_FIELDS = [
    "gene",
    "pmid",
    "action",
    "route",
    "available_context_path",
    "available_context_bytes",
    "fetch_outcome",
    "fetch_reason",
    "missing_rows",
    "missing_distinct_variants",
    "notes",
]


def _load_fetch_rows(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        rows = payload.get("rows") or payload.get("results") or []
        if isinstance(rows, list):
            return [row for row in rows if isinstance(row, dict)]
    raise ValueError(f"Unsupported fetch summary shape: {path}")


def _load_fetch_output_dir(path: Path | None) -> list[dict[str, Any]]:
    """Read partial/interrupted fetch_paywalled.py output directories."""
    if path is None:
        return []
    output_dir = repo_path(path).expanduser()
    rows_by_pmid: dict[str, dict[str, Any]] = {}
    for result_path in sorted(output_dir.glob("*/result.json")):
        try:
            payload = json.loads(result_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(payload, dict):
            continue
        pmid = normalize_pmid(payload.get("pmid") or result_path.parent.name)
        if not pmid:
            continue
        full_context = repo_path(
            payload.get("canonical_full_context_path")
            or output_dir / f"{pmid}_FULL_CONTEXT.md"
        ).expanduser()
        notes = payload.get("notes") or []
        reason = payload.get("error") or "; ".join(str(item) for item in notes)
        if full_context.exists() and is_usable_fulltext_source(full_context):
            outcome = "success_from_output_dir"
        elif payload.get("error"):
            outcome = "error"
        else:
            outcome = "empty"
        rows_by_pmid[pmid] = {
            "pmid": pmid,
            "outcome": outcome,
            "canonical_full_context_path": str(full_context),
            "reason": reason,
            "publisher": payload.get("publisher") or "",
            "final_url": payload.get("final_url") or "",
        }

    for context_path in sorted(output_dir.glob("*_FULL_CONTEXT.md")):
        pmid = normalize_pmid(context_path.name.removesuffix("_FULL_CONTEXT.md"))
        if not pmid or pmid in rows_by_pmid:
            continue
        if is_usable_fulltext_source(context_path):
            rows_by_pmid[pmid] = {
                "pmid": pmid,
                "outcome": "success_from_output_dir",
                "canonical_full_context_path": str(context_path),
                "reason": "usable flat FULL_CONTEXT recovered before summary write",
            }
    return list(rows_by_pmid.values())


def _pmids_from_refresh_summary(path: Path | None) -> tuple[set[str], set[str]]:
    """Return (attempted PMIDs, successful PMIDs) from refresh_run_db summary."""
    if path is None:
        return set(), set()
    summary_path = repo_path(path).expanduser()
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return set(), set()
    replay = payload.get("replay") or {}
    if not isinstance(replay, dict):
        return set(), set()

    attempted = {
        normalize_pmid(pmid)
        for pmid in replay.get("attempted_pmids", [])
        if normalize_pmid(pmid)
    }
    candidates_csv = str(payload.get("candidates_csv") or "").strip()
    if not attempted and candidates_csv:
        candidate_path = repo_path(candidates_csv).expanduser()
        if candidate_path.exists():
            attempted = {
                normalize_pmid(row.get("pmid"))
                for row in read_csv_rows(candidate_path)
                if normalize_pmid(row.get("pmid"))
            }

    successful = {
        normalize_pmid(pmid)
        for pmid in replay.get("successful_pmids", [])
        if normalize_pmid(pmid)
    }
    if successful:
        return attempted or successful, successful

    failed = {
        normalize_pmid(row.get("pmid"))
        for row in replay.get("errors", [])
        if isinstance(row, dict) and normalize_pmid(row.get("pmid"))
    }
    for key in ("failed_pmids", "gated_pmids"):
        failed.update(
            normalize_pmid(pmid) for pmid in replay.get(key, []) if normalize_pmid(pmid)
        )
    for key in ("gated_regressions", "gated_explosions"):
        failed.update(
            normalize_pmid(row.get("pmid"))
            for row in replay.get(key, [])
            if isinstance(row, dict) and normalize_pmid(row.get("pmid"))
        )
    return attempted, attempted - failed


def _gold_pmids_from_report(report: Path, gene: str) -> set[str]:
    return {
        normalize_pmid(row.get("pmid"))
        for row in read_csv_rows(report)
        if str(row.get("gene") or "").upper() == gene.upper()
        and normalize_pmid(row.get("pmid"))
    }


def _gold_pmids_from_gold_csv(path: Path) -> set[str]:
    return {
        normalize_pmid(row.get("pmid"))
        for row in read_csv_rows(path)
        if normalize_pmid(row.get("pmid"))
    }


def _recall_item(
    pmids: set[str],
    *,
    gold_total: int | None,
    missing_rows: int = 0,
    missing_distinct_variants: int = 0,
) -> dict[str, Any]:
    return {
        "pmids": len(pmids),
        "gold_pmids": gold_total,
        "recall": round(len(pmids) / gold_total, 6) if gold_total else None,
        "missing_rows": missing_rows,
        "missing_distinct_variants": missing_distinct_variants,
    }


def _coverage_item(
    pmids: set[str],
    *,
    total_pmids: int,
    missing_rows: int = 0,
    missing_distinct_variants: int = 0,
) -> dict[str, Any]:
    return {
        "pmids": len(pmids),
        "total_pmids": total_pmids,
        "coverage": round(len(pmids) / total_pmids, 6) if total_pmids else None,
        "missing_rows": missing_rows,
        "missing_distinct_variants": missing_distinct_variants,
    }


def _row_counts(
    rows_by_pmid: dict[str, dict[str, str]], pmids: set[str]
) -> tuple[int, int]:
    missing_rows = 0
    missing_distinct = 0
    for pmid in pmids:
        row = rows_by_pmid.get(pmid) or {}
        missing_rows += parse_int(row.get("missing_rows")) or 0
        missing_distinct += parse_int(row.get("missing_distinct_variants")) or 0
    return missing_rows, missing_distinct


def _fetch_source_path(row: dict[str, Any]) -> Path | None:
    for key in ("canonical_path", "path", "canonical_full_context_path"):
        value = str(row.get(key) or "").strip()
        if value:
            return repo_path(value).expanduser()
    return None


def successful_fetch_sources(
    fetch_rows: list[dict[str, Any]],
) -> dict[str, tuple[Path, dict[str, Any]]]:
    """Return PMID -> (usable full-text path, fetch row)."""
    successes: dict[str, tuple[Path, dict[str, Any]]] = {}
    for row in fetch_rows:
        pmid = normalize_pmid(row.get("pmid") or row.get("PMID"))
        if not pmid or str(row.get("outcome") or "") not in FETCH_SUCCESS_OUTCOMES:
            continue
        path = _fetch_source_path(row)
        if path is None or not path.exists() or not is_usable_fulltext_source(path):
            continue
        current = successes.get(pmid)
        current_size = current[0].stat().st_size if current else -1
        if path.stat().st_size > current_size:
            successes[pmid] = (path, row)
    return successes


def build_summary(
    *,
    gene: str,
    worklist_rows: list[dict[str, str]],
    fetch_rows: list[dict[str, Any]],
    refresh_attempted_pmids: set[str] | None = None,
    refresh_success_pmids: set[str] | None = None,
    gold_pmids: set[str] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows_by_pmid = {
        normalize_pmid(row.get("pmid")): row
        for row in worklist_rows
        if normalize_pmid(row.get("pmid"))
    }
    selected_fetch = {
        pmid for pmid, row in rows_by_pmid.items() if row.get("action") == "fetch"
    }
    selected_source = {
        pmid
        for pmid, row in rows_by_pmid.items()
        if row.get("action") in {"fetch", "refresh_replay", "manual_or_blocked"}
    }
    attempted = {
        normalize_pmid(row.get("pmid") or row.get("PMID"))
        for row in fetch_rows
        if normalize_pmid(row.get("pmid") or row.get("PMID"))
    }
    successes = successful_fetch_sources(fetch_rows)
    successful_pmids = set(successes)

    selected_fetch_rows, selected_fetch_distinct = _row_counts(
        rows_by_pmid, selected_fetch
    )
    selected_source_rows, selected_source_distinct = _row_counts(
        rows_by_pmid, selected_source
    )
    attempted_rows, attempted_distinct = _row_counts(rows_by_pmid, attempted)
    success_rows, success_distinct = _row_counts(rows_by_pmid, successful_pmids)
    refresh_attempted_rows, refresh_attempted_distinct = _row_counts(
        rows_by_pmid, refresh_attempted_pmids or set()
    )
    refresh_success_rows, refresh_success_distinct = _row_counts(
        rows_by_pmid, refresh_success_pmids or set()
    )

    gold_total = len(gold_pmids) if gold_pmids is not None else None
    worklist_total = len(rows_by_pmid)
    by_outcome = Counter(str(row.get("outcome") or "unknown") for row in fetch_rows)
    by_success_outcome = Counter(
        str(row.get("outcome") or "unknown")
        for _pmid, (_path, row) in successes.items()
    )
    denominator_text = "Gold PMIDs" if gold_pmids is not None else "Worklist PMIDs"

    source_override_rows: list[dict[str, Any]] = []
    for pmid in sorted(successful_pmids):
        path, fetch_row = successes[pmid]
        worklist_row = rows_by_pmid.get(pmid) or {}
        source_override_rows.append(
            {
                "gene": gene.upper(),
                "pmid": pmid,
                "action": "refresh_replay",
                "route": "fetched_fulltext",
                "available_context_path": str(path),
                "available_context_bytes": path.stat().st_size,
                "fetch_outcome": fetch_row.get("outcome") or "",
                "fetch_reason": fetch_row.get("reason") or "",
                "missing_rows": worklist_row.get("missing_rows") or 0,
                "missing_distinct_variants": worklist_row.get(
                    "missing_distinct_variants"
                )
                or 0,
                "notes": "usable full text downloaded/extracted by fetch_paywalled",
            }
        )

    summary = {
        "gene": gene.upper(),
        "gold_pmids": gold_total,
        "worklist_pmids": len(rows_by_pmid),
        "fetch_summary_pmids": len(attempted),
        "refresh_attempted_pmids": (
            len(refresh_attempted_pmids)
            if refresh_attempted_pmids is not None
            else None
        ),
        "refresh_success_pmids": (
            len(refresh_success_pmids) if refresh_success_pmids is not None else None
        ),
        "denominator": "gold_pmids" if gold_pmids is not None else "worklist_pmids",
        "pmid_recall": {
            "selected_for_fetch_download": _recall_item(
                selected_fetch,
                gold_total=gold_total,
                missing_rows=selected_fetch_rows,
                missing_distinct_variants=selected_fetch_distinct,
            ),
            "selected_for_source_acquisition_or_binding": _recall_item(
                selected_source,
                gold_total=gold_total,
                missing_rows=selected_source_rows,
                missing_distinct_variants=selected_source_distinct,
            ),
            "fetch_attempted": _recall_item(
                attempted,
                gold_total=gold_total,
                missing_rows=attempted_rows,
                missing_distinct_variants=attempted_distinct,
            ),
            "usable_fulltext_downloaded": _recall_item(
                successful_pmids,
                gold_total=gold_total,
                missing_rows=success_rows,
                missing_distinct_variants=success_distinct,
            ),
        },
        "pmid_coverage": {
            "selected_for_fetch_download": _coverage_item(
                selected_fetch,
                total_pmids=worklist_total,
                missing_rows=selected_fetch_rows,
                missing_distinct_variants=selected_fetch_distinct,
            ),
            "selected_for_source_acquisition_or_binding": _coverage_item(
                selected_source,
                total_pmids=worklist_total,
                missing_rows=selected_source_rows,
                missing_distinct_variants=selected_source_distinct,
            ),
            "fetch_attempted": _coverage_item(
                attempted,
                total_pmids=worklist_total,
                missing_rows=attempted_rows,
                missing_distinct_variants=attempted_distinct,
            ),
            "usable_fulltext_downloaded": _coverage_item(
                successful_pmids,
                total_pmids=worklist_total,
                missing_rows=success_rows,
                missing_distinct_variants=success_distinct,
            ),
        },
        "fetch_outcomes": dict(sorted(by_outcome.items())),
        "usable_fulltext_outcomes": dict(sorted(by_success_outcome.items())),
        "notes": {
            "selected_for_fetch_download": (
                f"{denominator_text} queued for download/refetch."
            ),
            "selected_for_source_acquisition_or_binding": (
                f"{denominator_text} queued for download/refetch, blocked manual "
                "acquisition, or replay from already available full text."
            ),
            "usable_fulltext_downloaded": (
                "Fetch summary PMIDs whose outcome was successful and whose "
                "FULL_CONTEXT file exists and passes the usable-fulltext gate."
            ),
        },
    }
    if refresh_attempted_pmids is not None:
        summary["pmid_recall"]["source_refresh_attempted"] = _recall_item(
            refresh_attempted_pmids,
            gold_total=gold_total,
            missing_rows=refresh_attempted_rows,
            missing_distinct_variants=refresh_attempted_distinct,
        )
        summary["pmid_coverage"]["source_refresh_attempted"] = _coverage_item(
            refresh_attempted_pmids,
            total_pmids=worklist_total,
            missing_rows=refresh_attempted_rows,
            missing_distinct_variants=refresh_attempted_distinct,
        )
        summary["notes"]["source_refresh_attempted"] = (
            "PMIDs handed to refresh_run_db from acquisition or source-binding output."
        )
    if refresh_success_pmids is not None:
        summary["pmid_recall"]["source_refresh_successful"] = _recall_item(
            refresh_success_pmids,
            gold_total=gold_total,
            missing_rows=refresh_success_rows,
            missing_distinct_variants=refresh_success_distinct,
        )
        summary["pmid_coverage"]["source_refresh_successful"] = _coverage_item(
            refresh_success_pmids,
            total_pmids=worklist_total,
            missing_rows=refresh_success_rows,
            missing_distinct_variants=refresh_success_distinct,
        )
        summary["notes"]["source_refresh_successful"] = (
            "Refresh PMIDs that produced accepted extraction JSON; this is the "
            "strict actual full-text downloaded-and-extracted count."
        )
    return summary, source_override_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--worklist", required=True, type=Path, help="Worklist CSV")
    parser.add_argument(
        "--fetch-summary",
        type=Path,
        default=None,
        help="fetch_paywalled.py summary.json. Optional before fetch has run.",
    )
    parser.add_argument(
        "--fetch-output-dir",
        type=Path,
        default=None,
        help=(
            "fetch_paywalled.py output directory. Useful when a long fetch was "
            "interrupted before summary.json was written."
        ),
    )
    parser.add_argument(
        "--refresh-summary",
        type=Path,
        default=None,
        help=(
            "Optional refresh_run_db refresh_summary.json used to report PMIDs "
            "that were actually accepted into refreshed extraction JSON."
        ),
    )
    denominator = parser.add_mutually_exclusive_group(required=False)
    denominator.add_argument(
        "--report",
        type=Path,
        help="paper_disagreement_report.csv used to count gold PMIDs for recall.",
    )
    denominator.add_argument(
        "--gold",
        type=Path,
        help="Normalized gold recall input CSV used to count gold PMIDs for recall.",
    )
    parser.add_argument("--out", required=True, type=Path, help="Summary JSON path")
    parser.add_argument(
        "--source-override-out",
        type=Path,
        default=None,
        help="Optional CSV of successful fetches for refresh_run_db.py.",
    )
    args = parser.parse_args()

    gene = args.gene.upper()
    worklist = repo_path(args.worklist).expanduser()
    fetch_summary = (
        repo_path(args.fetch_summary).expanduser() if args.fetch_summary else None
    )
    fetch_output_dir = (
        repo_path(args.fetch_output_dir).expanduser() if args.fetch_output_dir else None
    )
    refresh_summary = (
        repo_path(args.refresh_summary).expanduser() if args.refresh_summary else None
    )
    gold_pmids: set[str] | None
    if args.report:
        gold_pmids = _gold_pmids_from_report(repo_path(args.report).expanduser(), gene)
    elif args.gold:
        gold_pmids = _gold_pmids_from_gold_csv(repo_path(args.gold).expanduser())
    else:
        gold_pmids = None

    refresh_attempted_pmids, refresh_success_pmids = _pmids_from_refresh_summary(
        refresh_summary
    )
    summary, source_override_rows = build_summary(
        gene=gene,
        worklist_rows=read_csv_rows(worklist),
        fetch_rows=[
            *_load_fetch_rows(fetch_summary),
            *_load_fetch_output_dir(fetch_output_dir),
        ],
        refresh_attempted_pmids=(refresh_attempted_pmids if refresh_summary else None),
        refresh_success_pmids=(refresh_success_pmids if refresh_summary else None),
        gold_pmids=gold_pmids,
    )

    out = repo_path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    if args.source_override_out:
        write_csv_rows(
            source_override_rows,
            SOURCE_OVERRIDE_FIELDS,
            repo_path(args.source_override_out).expanduser(),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
