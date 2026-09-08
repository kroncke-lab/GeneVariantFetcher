"""Bounded, read-only coverage audit for deterministic table extraction.

This records what the extractor supplied and flags printed clinical table
headings for review. It does not establish table relevance, eligibility, or
field completeness, and never maps symptoms/ECGs/events onto a disease endpoint.
"""

from __future__ import annotations

import hashlib
import math
import re
from typing import Any

_TABLE_HEADING = re.compile(
    r"^(?:#{1,6}\s*)?(?:\*\*)?\s*"
    r"(?:(?:online|supplementary|supplemental|supp\.?)\s+|e)?"
    r"table\s+[a-z]?\d+[a-z]?(?:\s*[:.]\s*|\s*$)",
    re.IGNORECASE,
)
_CLINICAL_SIGNAL = re.compile(
    r"\b(?:clinical|phenotyp\w*|symptom\w*|asymptomatic|affected|unaffected|"
    r"ecg|qtc|diagnos\w*|cardiac\s+events?|follow[- ]?up)\b",
    re.IGNORECASE,
)
_FIELDS = ("total_carriers_observed", "affected_count", "unaffected_count")
MAX_CANDIDATES = 20


def audit_table_phenotype_coverage(
    extracted_data: dict[str, Any], source_text: str
) -> dict[str, Any]:
    """Describe pre-guard output and labelled clinical-table candidates.

    Candidate discovery is deliberately a heading/header heuristic, not proof
    that a source table holds a recoverable variant-specific count. Missing or
    unlabelled tables cannot be ruled out even when the candidate list is empty.
    No source excerpt, prompt, or API request is constructed by this audit.
    """
    raw_rows = extracted_data.get("variants") or []
    rows = [v for v in raw_rows if isinstance(v, dict)]
    supplied = {field: 0 for field in _FIELDS}
    for row in rows:
        counts = row.get("penetrance_data") or {}
        if not isinstance(counts, dict):
            counts = {}
        for field in _FIELDS:
            value = counts.get(field)
            if field == "total_carriers_observed" and value is None:
                patients = row.get("patients") or {}
                value = patients.get("count") if isinstance(patients, dict) else None
            # Count explicit zero as supplied; null and malformed values are
            # not observations. This is extraction coverage, not trust status.
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                if value >= 0 and (
                    isinstance(value, int)
                    or (math.isfinite(value) and value.is_integer())
                ):
                    supplied[field] += 1

    # PDF text contains form feeds. They are page separators, not editor line
    # breaks; splitlines() would silently shift every later source locator.
    lines = source_text.split("\n")
    headings = [i for i, line in enumerate(lines) if _TABLE_HEADING.match(line.strip())]
    candidates = []
    candidate_count = 0
    for position, start in enumerate(headings):
        stop = min(
            start + 12,
            headings[position + 1] if position + 1 < len(headings) else len(lines),
        )
        signals = sorted(
            {
                m.group(0).lower()
                for line in lines[start:stop]
                for m in _CLINICAL_SIGNAL.finditer(line)
            }
        )
        if not signals:
            continue
        candidate_count += 1
        if len(candidates) >= MAX_CANDIDATES:
            continue
        candidates.append(
            {
                "source_line_1based": start + 1,
                "caption": " ".join(lines[start].split())[:240],
                "clinical_signals": signals[:12],
                "signals_truncated": len(signals) > 12,
            }
        )
    missing = [field for field in _FIELDS[1:] if supplied[field] < len(rows)]
    return {
        "protocol_version": "table_phenotype_coverage_v1",
        "stage": "after_derivation_before_phenotype_guard",
        "source_sha256": hashlib.sha256(source_text.encode("utf-8")).hexdigest(),
        "source_representation": "exact_input_text",
        "variant_rows": len(rows),
        "supplied_rows": supplied,
        "missing_phenotype_fields": missing,
        "clinical_table_candidate_count": candidate_count,
        "clinical_table_candidates": candidates,
        "candidates_truncated": candidate_count > MAX_CANDIDATES,
        "status": "review_candidates"
        if candidates and missing
        else "completeness_not_established",
        "interpretation": (
            "Heading/header candidates require source review for variant ownership, "
            "target disease, endpoint and timepoint. Supplied rows are pre-guard, "
            "not validated coverage. No candidate list establishes completeness."
        ),
    }
