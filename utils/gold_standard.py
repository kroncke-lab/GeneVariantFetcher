"""Shared fail-closed semantics for adjudicated count columns in gold CSVs."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, TypeVar


GOLD_COUNT_FIELDS = ("carriers", "affected", "unaffected")

# Status values are data-schema vocabulary, not free-form notes. Adding a new
# status requires code review so a typo or work-in-progress marker cannot
# silently remove a gold assertion.
GOLD_V2_STATUSES = frozenset(
    {
        "adjudicated_current_study_cohort",
        "adjudicated_current_treatment_cohort",
        "adjudicated_genotyped_carrier_cohort",
        "adjudicated_null_unaffected",
        # The paper never reports a per-variant phenotype: carriers stand,
        # affected and unaffected are explicit nulls (SCN5A 30059973, 2026-09-08).
        "adjudicated_phenotype_not_reported",
        # The paper prints a per-variant phenotype split; unaffected is its
        # negative-phenotype column, affected every phenotype-positive carrier.
        "adjudicated_source_phenotype_partition",
        "adjudicated_variant_carrier_count",
        "confirmed_original",
        "confirmed_original_pedigree_derived",
        "confirmed_original_table_derived",
        "excluded_duplicate_current_cohort",
        "schema_policy_paper_defined_affected",
    }
)
EXCLUDED_GOLD_V2_STATUSES = frozenset({"excluded_duplicate_current_cohort"})

T = TypeVar("T")

# Append-only record of adjudications that changed a gold file's bytes. A
# registry pins the digest of the gold it was built from; after an approved
# adjudication the live file legitimately differs, and the chain here is what
# lets a provenance check tell an approved revision from silent drift.
GOLD_REVISIONS_LOG = (
    Path(__file__).resolve().parents[1]
    / "gene_variant_fetcher_gold_standard"
    / "gold_revisions.jsonl"
)


def gold_digest_lineage(path: Path, *, log_path: Path | None = None) -> list[str]:
    """Return the file's current sha256 followed by every recorded prior digest.

    Entries in ``gold_revisions.jsonl`` chain ``previous_sha256`` to ``sha256``
    for one relative path; the chain must end at the current file, otherwise
    the log does not describe the file on disk and only the current digest is
    returned (so an unrecorded edit still fails a pinned-digest check).
    """
    current = hashlib.sha256(path.read_bytes()).hexdigest()
    lineage = [current]
    log = GOLD_REVISIONS_LOG if log_path is None else log_path
    if not log.is_file():
        return lineage
    repo_root = log.resolve().parents[1]
    try:
        relative = path.resolve().relative_to(repo_root).as_posix()
    except ValueError:
        relative = path.as_posix()
    entries = []
    for line in log.read_text().splitlines():
        if line.strip():
            entry = json.loads(line)
            if entry.get("path") == relative:
                entries.append(entry)
    by_new = {str(e.get("sha256")): e for e in entries}
    cursor = current
    while cursor in by_new:
        previous = str(by_new[cursor].get("previous_sha256") or "")
        if not previous or previous in lineage:
            break
        lineage.append(previous)
        cursor = previous
    return lineage


def _clean_status(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if text in {"<NA>", "NaN", "nan"}:
        return ""
    return text


def gold_v2_status(row: Mapping[str, Any]) -> str:
    """Return a validated v2 status, raising on open-vocabulary values."""

    status = _clean_status(row.get("gold_v2_status"))
    if status and status not in GOLD_V2_STATUSES:
        allowed = ", ".join(sorted(GOLD_V2_STATUSES))
        raise ValueError(
            f"Unknown gold_v2_status {status!r}; expected one of: {allowed}"
        )
    return status


def gold_row_excluded(row: Mapping[str, Any]) -> bool:
    """Return whether a curator status removes the row from every score."""

    return gold_v2_status(row) in EXCLUDED_GOLD_V2_STATUSES


def authoritative_gold_count(
    row: Mapping[str, Any],
    field: str,
    *,
    parser: Callable[[Any], T],
) -> T:
    """Resolve one count without ever backfilling an explicit v2 null."""

    if field not in GOLD_COUNT_FIELDS:
        raise ValueError(f"Unsupported gold count field: {field!r}")
    source = f"gold_v2_{field}" if gold_v2_status(row) else field
    return parser(row.get(source))
