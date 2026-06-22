"""Refuse per-variant count assignments whose provenance is not per-variant.

Sibling of :mod:`pipeline.count_outlier_guard`. Where the outlier guard catches
suspicious values via a statistical (per-paper-median) rule, this classifier
catches them via the LLM's own self-reported attribution: each variant's
``count_provenance`` block (populated by the prompts in
:mod:`pipeline.prompts`) records which raw column the count came from and
classifies it as one of a closed vocabulary of count types.

Only ``per_variant_carrier`` should populate ``patients.count`` /
``penetrance_data.*``. Any other count type (cohort total, screened N,
proband-only, family count, case/control arm, etc.) is study-level, not
per-variant, and assigning it to a single row is the dominant source of the
MAE outliers documented in ``docs/RECALL_STATUS.md`` (KCNQ1 29622001 G589D
7 vs 453 is the canonical example — the 453 is the cohort total).

Pure internal-consistency check, gold-free. Works on new gene-diseases.

Public API mirrors the outlier guard for consistency:
    detect_misclassified_counts(variants, ...) -> list[ClassifierAnnotation]
    enforce_per_variant_policy(variants, annotations, policy) -> ClassifierResult
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from pipeline.count_outlier_guard import COUNT_FIELDS, _clear_field, _read_field

ACCEPTED_COUNT_TYPE = "per_variant_carrier"

# Closed vocabulary from the prompt schema; mirrors pipeline/prompts.py.
KNOWN_COUNT_TYPES = {
    "per_variant_carrier",
    "family_count",
    "proband_count",
    "cohort_total",
    "screened_N",
    "case",
    "control",
    "unaffected_control",
    "unknown",
}

# Map logical field name -> the keys the LLM uses inside count_provenance.
PROVENANCE_KEYS = {
    "carriers": ("carriers_column_label", "carriers_count_type"),
    "affected": ("affected_column_label", "affected_count_type"),
    "unaffected": ("unaffected_column_label", "unaffected_count_type"),
}

_SUSPICIOUS_LABEL_PATTERNS = (
    (
        "row identifier",
        re.compile(
            r"\b(?:adult|child|patient|proband|subject|case|individual|"
            r"participant|family|kindred|pedigree|sample)\s*(?:number|no\.?|id)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "population allele frequency/count",
        re.compile(
            r"\b(?:maf|minor allele frequency|allele frequency|eaf|frequency|"
            r"occurrences?|allele count|alleles|gnomad|exac)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "genotype observation count",
        re.compile(r"^\s*(?:het|hom|homo)\s*$|\bgenotype\b", re.IGNORECASE),
    ),
    (
        "cohort denominator",
        re.compile(
            r"\b(?:total\s+cases?|total\s+controls?|sample size|screened|"
            r"n\s+tested|cohort size)\b",
            re.IGNORECASE,
        ),
    ),
)


def _suspicious_column_label_reason(value: Any) -> Optional[str]:
    if not isinstance(value, str) or not value.strip():
        return None
    for reason, pattern in _SUSPICIOUS_LABEL_PATTERNS:
        if pattern.search(value):
            return reason
    return None


@dataclass
class ClassifierAnnotation:
    """One (variant_index, field) refused because count_type is not per-variant."""

    variant_index: int
    field: str
    value: int
    declared_count_type: Optional[str]
    column_label: Optional[str]
    reason: str


@dataclass
class ClassifierResult:
    flagged: int = 0
    cleared: int = 0
    annotations: list[ClassifierAnnotation] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "flagged": self.flagged,
            "cleared": self.cleared,
            "annotations": [
                {
                    "variant_index": a.variant_index,
                    "field": a.field,
                    "value": a.value,
                    "declared_count_type": a.declared_count_type,
                    "column_label": a.column_label,
                    "reason": a.reason,
                }
                for a in self.annotations
            ],
        }


def _normalize_count_type(value: Any) -> Optional[str]:
    if not isinstance(value, str):
        return None
    v = value.strip()
    if not v:
        return None
    # Allow the LLM to be slightly sloppy about capitalization of the underscore form.
    canonical = v.replace(" ", "_").replace("-", "_")
    # Preserve the LLM's value, but expose a canonical-cased version if it matches.
    for known in KNOWN_COUNT_TYPES:
        if canonical.lower() == known.lower():
            return known
    return canonical  # unknown value passes through so callers can log it


def detect_misclassified_counts(
    variants: list[dict[str, Any]],
    *,
    fields: Optional[Iterable[str]] = None,
) -> list[ClassifierAnnotation]:
    """Flag (variant_index, field) where the LLM-declared count_type is not
    ``per_variant_carrier`` AND the variant currently has a non-null count.

    Variants without a ``count_provenance`` block are silently skipped — the
    classifier is opt-in based on the LLM having populated the new schema.
    """
    field_keys = list(fields) if fields else list(PROVENANCE_KEYS.keys())
    annotations: list[ClassifierAnnotation] = []

    for idx, variant in enumerate(variants):
        if not isinstance(variant, dict):
            continue
        provenance = variant.get("count_provenance")
        if not isinstance(provenance, dict):
            continue
        for fkey in field_keys:
            label_key, type_key = PROVENANCE_KEYS[fkey]
            declared = _normalize_count_type(provenance.get(type_key))
            column_label = provenance.get(label_key) or None
            suspicious_label = _suspicious_column_label_reason(column_label)
            if (
                declared is None or declared == ACCEPTED_COUNT_TYPE
            ) and not suspicious_label:
                continue
            paths = COUNT_FIELDS.get(fkey)
            if not paths:
                continue
            current_value = _read_field(variant, paths)
            if current_value is None:
                # No value to refuse — provenance says non-per-variant and
                # the count is already null, so nothing to do.
                continue
            annotations.append(
                ClassifierAnnotation(
                    variant_index=idx,
                    field=fkey,
                    value=current_value,
                    declared_count_type=declared or "suspicious_column_label",
                    column_label=column_label,
                    reason=(
                        f"column label looks like {suspicious_label}; value "
                        f"{current_value} would be a non-per-variant assignment"
                        if suspicious_label
                        else (
                            f"declared count_type {declared!r} is not "
                            f"{ACCEPTED_COUNT_TYPE!r}; value {current_value} "
                            f"would be a non-per-variant assignment"
                        )
                    ),
                )
            )
    return annotations


def enforce_per_variant_policy(
    variants: list[dict[str, Any]],
    annotations: list[ClassifierAnnotation],
    *,
    policy: str = "flag",
) -> ClassifierResult:
    """Apply a policy to misclassified count assignments.

    Policies:
        "off"   - no-op, returns annotations only.
        "flag"  - annotate variants with `count_classifier_flags` metadata;
                  preserve raw values. Default and safest.
        "clear" - also null the offending count fields (raw preserved under
                  count_classifier_flags). Same shape as the outlier guard's
                  "clear" so the two policies compose cleanly.
    """
    if policy not in {"off", "flag", "clear"}:
        raise ValueError(f"unknown count-classifier policy: {policy!r}")
    result = ClassifierResult(annotations=list(annotations))
    if not annotations or policy == "off":
        return result

    for ann in annotations:
        if ann.variant_index >= len(variants):
            continue
        variant = variants[ann.variant_index]
        if not isinstance(variant, dict):
            continue
        flags_bucket = variant.setdefault("count_classifier_flags", {})
        flags_bucket[ann.field] = {
            "raw": ann.value,
            "declared_count_type": ann.declared_count_type,
            "column_label": ann.column_label,
            "reason": ann.reason,
            "policy": policy,
        }
        result.flagged += 1
        if policy == "clear":
            paths = COUNT_FIELDS[ann.field]
            _clear_field(variant, paths)
            result.cleared += 1
    return result
