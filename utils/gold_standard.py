"""Shared fail-closed semantics for adjudicated count columns in gold CSVs."""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
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
