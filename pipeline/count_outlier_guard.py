"""Per-PMID count-outlier detection and policy enforcement.

The 2026-05-27 MAE audit revealed that the largest rows-mode count errors
all follow one pattern: a study-wide N, cohort total, or screened-population
count gets assigned to a single per-variant row. Examples (matched-row, gold
vs extracted): KCNQ1 29622001 G589D 7 vs 453; KCNQ1 17192539 G589D 50 vs
243; SCN5A 16453024 S1103Y 1 vs 137; KCNH2 19160088 L552S 2 vs 92.

This module flags such rows by comparing each row's count against the
per-paper median (gold-free, internal-consistency only) and exposes a
policy step that can either flag-only or clear the flagged counts.

The guard is required to be gold-free by construction so it works on new
gene-diseases for which no curated answer exists.

Public API:
    detect_count_outliers(variants, ...) -> list[OutlierAnnotation]
    apply_outlier_policy(variants, annotations, policy) -> dict[str, int]

The variant dict shape matches `pipeline.extraction` output:
    variant["patients"]["count"]                 (carriers)
    variant["penetrance_data"]["affected_count"]
    variant["penetrance_data"]["unaffected_count"]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median
from typing import Any, Iterable, Optional

DEFAULT_MIN_VARIANTS = 4
DEFAULT_MULTIPLIER_THRESHOLD = 10.0
DEFAULT_ABSOLUTE_THRESHOLD = 50

COUNT_FIELDS = {
    "carriers": ("patients", "count"),
    "affected": ("penetrance_data", "affected_count"),
    "unaffected": ("penetrance_data", "unaffected_count"),
}


@dataclass
class OutlierAnnotation:
    """A single (variant_index, field) flagged as a likely study-wide-N reuse."""

    variant_index: int
    field: str
    value: int
    median: float
    multiplier: float
    reason: str


@dataclass
class GuardResult:
    """Returned from apply_outlier_policy for telemetry."""

    flagged: int = 0
    cleared: int = 0
    skipped_small_paper: bool = False
    annotations: list[OutlierAnnotation] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "flagged": self.flagged,
            "cleared": self.cleared,
            "skipped_small_paper": self.skipped_small_paper,
            "annotations": [
                {
                    "variant_index": a.variant_index,
                    "field": a.field,
                    "value": a.value,
                    "median": a.median,
                    "multiplier": a.multiplier,
                    "reason": a.reason,
                }
                for a in self.annotations
            ],
        }


def _read_nested(variant: dict[str, Any], path: tuple[str, str]) -> Optional[int]:
    container = variant.get(path[0])
    if not isinstance(container, dict):
        return None
    raw = container.get(path[1])
    try:
        return int(raw) if raw is not None else None
    except (TypeError, ValueError):
        return None


def _write_nested(variant: dict[str, Any], path: tuple[str, str], value: Any) -> None:
    container = variant.setdefault(path[0], {})
    if not isinstance(container, dict):
        return
    container[path[1]] = value


def detect_count_outliers(
    variants: list[dict[str, Any]],
    *,
    min_variants: int = DEFAULT_MIN_VARIANTS,
    multiplier_threshold: float = DEFAULT_MULTIPLIER_THRESHOLD,
    absolute_threshold: int = DEFAULT_ABSOLUTE_THRESHOLD,
    fields: Optional[Iterable[str]] = None,
) -> list[OutlierAnnotation]:
    """Flag per-variant count values that look like study-wide-N reuse.

    A value is flagged when, within a paper with at least ``min_variants``
    variants, the value is both (a) greater than ``absolute_threshold`` and
    (b) greater than ``multiplier_threshold`` * the median non-null,
    non-zero value across the other variants for the same field.

    Returns a list of OutlierAnnotation; empty when the paper is too small or
    when no value exceeds both thresholds. Pure read-only — does not mutate
    ``variants``.
    """
    if not isinstance(variants, list) or len(variants) < min_variants:
        return []

    field_keys = list(fields) if fields else list(COUNT_FIELDS.keys())
    annotations: list[OutlierAnnotation] = []

    for fkey in field_keys:
        path = COUNT_FIELDS.get(fkey)
        if path is None:
            continue
        values: list[tuple[int, int]] = []
        for idx, variant in enumerate(variants):
            if not isinstance(variant, dict):
                continue
            v = _read_nested(variant, path)
            if v is None or v <= 0:
                continue
            values.append((idx, v))
        if len(values) < min_variants:
            continue
        med = median(v for _, v in values)
        if med <= 0:
            continue
        for idx, v in values:
            if v <= absolute_threshold:
                continue
            mult = v / med if med else float("inf")
            if mult < multiplier_threshold:
                continue
            annotations.append(
                OutlierAnnotation(
                    variant_index=idx,
                    field=fkey,
                    value=v,
                    median=med,
                    multiplier=mult,
                    reason=(
                        f"value {v} > {absolute_threshold} and >{multiplier_threshold:.1f}x "
                        f"per-paper median {med:.1f}"
                    ),
                )
            )
    return annotations


def apply_outlier_policy(
    variants: list[dict[str, Any]],
    annotations: list[OutlierAnnotation],
    *,
    policy: str = "flag",
) -> GuardResult:
    """Apply a policy to flagged outliers.

    Policies:
        "off"   - no-op, returns annotations only.
        "flag"  - adds `count_outlier_flags` metadata on each affected variant;
                  preserves the raw value. Default and safest.
        "clear" - replaces the flagged value with None and records the prior
                  value under `count_outlier_flags[field]["raw"]`. Use only
                  when the goal is to drop suspected study-wide-N reuse from
                  matched-row count comparisons.
    """
    if policy not in {"off", "flag", "clear"}:
        raise ValueError(f"unknown count-outlier policy: {policy!r}")
    result = GuardResult(annotations=list(annotations))
    if not annotations or policy == "off":
        return result

    for ann in annotations:
        if ann.variant_index >= len(variants):
            continue
        variant = variants[ann.variant_index]
        if not isinstance(variant, dict):
            continue
        flags_bucket = variant.setdefault("count_outlier_flags", {})
        flags_bucket[ann.field] = {
            "raw": ann.value,
            "median": ann.median,
            "multiplier": ann.multiplier,
            "reason": ann.reason,
            "policy": policy,
        }
        result.flagged += 1
        if policy == "clear":
            path = COUNT_FIELDS[ann.field]
            _write_nested(variant, path, None)
            result.cleared += 1
    return result
