"""Evidence cards for high-uncertainty count assignments (#6).

An ``EvidenceCard`` packages everything a verifier needs to decide whether a
single (PMID × variant × count field) row should be trusted: the extracted
value, the LLM's self-reported provenance, the source text around the value,
optional gold count, and the list of triggers that pulled the row into
review.

Two modes are supported by the trigger detector:

* **gold-mode** (``GOLD_ERROR_THRESHOLD`` default 10): fires when a matched
  row's |gold - extracted| ≥ threshold.
* **no-gold mode** (works on unknown gene-diseases per
  ``feedback-no-gold-turnkey``): fires when the row is flagged by the
  count-outlier guard, the count classifier, or has a "large" value with
  ``unknown`` / missing ``count_provenance``.

Verification is intentionally decoupled: ``EvidenceCard.run_heuristic_verdict``
applies a simple rule-based pass that flags surrounding text patterns like
"Total cohort", "N studied", "screened" as REJECT signals; other verifiers
(LLM-based, claim-verifier-based) can be plugged in later by walking the
``cards`` produced here.

Design rule: gold-mode is for measurement (audit-time confidence checks on
the validation surface). No-gold mode is the primary product path. Every
card carries both kinds of trigger so the same artifact can be reviewed in
either mode.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from pipeline.count_outlier_guard import COUNT_FIELDS, _read_field

GOLD_ERROR_THRESHOLD_DEFAULT = 10
LARGE_VALUE_THRESHOLD_DEFAULT = 50
SOURCE_EXCERPT_LINES_BEFORE = 3
SOURCE_EXCERPT_LINES_AFTER = 3

VERDICT_CONFIRM = "confirm"
VERDICT_REJECT = "reject"
VERDICT_WITHHOLD = "withhold"

# Patterns that strongly suggest the surrounding context describes a
# study-wide / cohort / screened count rather than a per-variant carrier
# count. Used by run_heuristic_verdict to flag REJECT signals.
REJECT_TEXT_PATTERNS = (
    re.compile(r"\btotal\s+cohort\b", re.IGNORECASE),
    re.compile(r"\bcohort\s+(of|total|size)\b", re.IGNORECASE),
    re.compile(r"\bN\s*=\s*\d+", re.IGNORECASE),
    re.compile(
        r"\bnumber\s+(of\s+)?(patients|subjects|individuals)\s+studied\b", re.IGNORECASE
    ),
    re.compile(r"\bscreened\b", re.IGNORECASE),
    re.compile(r"\benrolled\b", re.IGNORECASE),
    re.compile(r"\b(study|sample)\s+(size|population)\b", re.IGNORECASE),
    re.compile(r"\boverall\s+(N|n)\b"),
    re.compile(r"\b(in\s+)?total\b", re.IGNORECASE),
)


@dataclass
class EvidenceCard:
    """One row-under-review: PMID × variant × field."""

    pmid: str
    variant_index: int
    variant_identifier: str
    field: str  # "carriers" | "affected" | "unaffected"
    extracted_value: Optional[int]
    gold_value: Optional[int]
    error: Optional[int]  # |gold - extracted| when both present
    declared_count_type: Optional[str]
    column_label: Optional[str]
    source_excerpt: str
    triggers: list[str] = field(default_factory=list)
    verdict: Optional[str] = None  # set by a verifier
    verdict_reason: Optional[str] = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "pmid": self.pmid,
            "variant_index": self.variant_index,
            "variant_identifier": self.variant_identifier,
            "field": self.field,
            "extracted_value": self.extracted_value,
            "gold_value": self.gold_value,
            "error": self.error,
            "declared_count_type": self.declared_count_type,
            "column_label": self.column_label,
            "source_excerpt": self.source_excerpt,
            "triggers": list(self.triggers),
            "verdict": self.verdict,
            "verdict_reason": self.verdict_reason,
        }

    def run_heuristic_verdict(self) -> None:
        """Apply rule-based verifier; populate ``verdict`` and ``verdict_reason``.

        Rules (in order):
          1. If ``declared_count_type`` ∈ {cohort_total, screened_N} → REJECT.
          2. If the source excerpt contains any REJECT_TEXT_PATTERNS → REJECT.
          3. If the extracted value matches an explicit "N=NUMBER" in the
             surrounding text AND that number is also the extracted value
             → REJECT (LLM grabbed the study-wide N).
          4. Otherwise → WITHHOLD (let a downstream verifier decide; do NOT
             auto-confirm).

        Never returns CONFIRM — that requires positive per-variant evidence
        which the heuristic cannot guarantee.
        """
        ct = (self.declared_count_type or "").lower()
        if ct in {"cohort_total", "screened_n"}:
            self.verdict = VERDICT_REJECT
            self.verdict_reason = f"declared count_type={self.declared_count_type!r}"
            return
        for pat in REJECT_TEXT_PATTERNS:
            m = pat.search(self.source_excerpt or "")
            if m:
                self.verdict = VERDICT_REJECT
                self.verdict_reason = (
                    f"source excerpt matches study-wide pattern {pat.pattern!r}"
                )
                return
        # Rule 3: explicit N=<value> match
        if self.extracted_value is not None:
            for m in re.finditer(r"\bN\s*=\s*(\d+)", self.source_excerpt or ""):
                try:
                    if int(m.group(1)) == self.extracted_value:
                        self.verdict = VERDICT_REJECT
                        self.verdict_reason = (
                            f"extracted value {self.extracted_value} matches an "
                            f"explicit N=N in the surrounding source text"
                        )
                        return
                except ValueError:
                    continue
        self.verdict = VERDICT_WITHHOLD
        self.verdict_reason = "no positive per-variant evidence; defer to verifier"


def _excerpt_around(value: int, source_text: str) -> str:
    """Return a small excerpt around the first occurrence of ``value`` in the
    source markdown (best-effort; if not found, return the first 1000 chars
    so the verifier has SOMETHING to look at).
    """
    if not source_text:
        return ""
    needle = str(value)
    lines = source_text.splitlines()
    for i, line in enumerate(lines):
        if needle in line:
            start = max(0, i - SOURCE_EXCERPT_LINES_BEFORE)
            end = min(len(lines), i + SOURCE_EXCERPT_LINES_AFTER + 1)
            return "\n".join(lines[start:end])
    return source_text[:1000]


def _classifier_or_outlier_flagged(
    variant: dict[str, Any], field_name: str
) -> list[str]:
    """Return existing per-variant flags for ``field_name`` set by the
    outlier guard (#3) or the classifier (#5).
    """
    triggers: list[str] = []
    outlier_flags = variant.get("count_outlier_flags") or {}
    if isinstance(outlier_flags, dict) and field_name in outlier_flags:
        triggers.append("outlier_guard")
    classifier_flags = variant.get("count_classifier_flags") or {}
    if isinstance(classifier_flags, dict) and field_name in classifier_flags:
        triggers.append("classifier")
    return triggers


def _provenance_count_type(variant: dict[str, Any], field_name: str) -> Optional[str]:
    provenance = variant.get("count_provenance")
    if not isinstance(provenance, dict):
        return None
    key = {
        "carriers": "carriers_count_type",
        "affected": "affected_count_type",
        "unaffected": "unaffected_count_type",
    }.get(field_name)
    if not key:
        return None
    v = provenance.get(key)
    return v if isinstance(v, str) and v.strip() else None


def _provenance_column(variant: dict[str, Any], field_name: str) -> Optional[str]:
    provenance = variant.get("count_provenance")
    if not isinstance(provenance, dict):
        return None
    key = {
        "carriers": "carriers_column_label",
        "affected": "affected_column_label",
        "unaffected": "unaffected_column_label",
    }.get(field_name)
    if not key:
        return None
    v = provenance.get(key)
    return v if isinstance(v, str) and v.strip() else None


def _variant_identifier(variant: dict[str, Any]) -> str:
    for k in ("protein_notation", "cdna_notation", "genomic_position"):
        v = variant.get(k)
        if isinstance(v, str) and v.strip():
            return v
    return "<no-identifier>"


def build_evidence_cards(
    *,
    pmid: str,
    variants: list[dict[str, Any]],
    source_text: str = "",
    gold_counts_by_variant: Optional[dict[str, dict[str, int]]] = None,
    gold_error_threshold: int = GOLD_ERROR_THRESHOLD_DEFAULT,
    large_value_threshold: int = LARGE_VALUE_THRESHOLD_DEFAULT,
    fields: Optional[Iterable[str]] = None,
) -> list[EvidenceCard]:
    """Build evidence cards for any (variant × field) row that warrants review.

    Triggers (dual-mode; a single card may carry multiple):

    * ``gold_error_ge_<N>`` — gold-mode: |gold - extracted| ≥ threshold.
    * ``outlier_guard``    — variant has ``count_outlier_flags[field]``.
    * ``classifier``       — variant has ``count_classifier_flags[field]``.
    * ``unknown_provenance_large`` — declared count_type is ``unknown`` (or
      provenance absent) AND extracted value ≥ ``large_value_threshold``.

    Rows without any trigger do not produce a card. The caller passes
    ``gold_counts_by_variant`` only in gold-mode; in no-gold mode it is
    ``None`` and the gold-error trigger never fires.

    Returns: list of EvidenceCard. Verdicts are NOT populated — call
    ``card.run_heuristic_verdict()`` or hand the card to a downstream
    LLM-based verifier.
    """
    field_keys = list(fields) if fields else list(COUNT_FIELDS.keys())
    cards: list[EvidenceCard] = []
    has_gold = gold_counts_by_variant is not None

    for idx, variant in enumerate(variants):
        if not isinstance(variant, dict):
            continue
        var_id = _variant_identifier(variant)
        for fkey in field_keys:
            paths = COUNT_FIELDS.get(fkey)
            if not paths:
                continue
            extracted = _read_field(variant, paths)

            # No-gold triggers (always evaluated)
            triggers = _classifier_or_outlier_flagged(variant, fkey)
            ct = _provenance_count_type(variant, fkey)
            if (
                extracted is not None
                and extracted >= large_value_threshold
                and (ct is None or ct.lower() == "unknown")
            ):
                triggers.append("unknown_provenance_large")

            # Gold trigger
            gold_value: Optional[int] = None
            err: Optional[int] = None
            if has_gold:
                gold_for_variant = gold_counts_by_variant.get(var_id, {})
                gv_raw = gold_for_variant.get(fkey)
                try:
                    gold_value = int(gv_raw) if gv_raw is not None else None
                except (TypeError, ValueError):
                    gold_value = None
                if gold_value is not None and extracted is not None:
                    err = abs(gold_value - extracted)
                    if err >= gold_error_threshold:
                        triggers.append(f"gold_error_ge_{gold_error_threshold}")

            if not triggers:
                continue

            excerpt = (
                _excerpt_around(extracted, source_text) if extracted is not None else ""
            )
            cards.append(
                EvidenceCard(
                    pmid=pmid,
                    variant_index=idx,
                    variant_identifier=var_id,
                    field=fkey,
                    extracted_value=extracted,
                    gold_value=gold_value,
                    error=err,
                    declared_count_type=ct,
                    column_label=_provenance_column(variant, fkey),
                    source_excerpt=excerpt,
                    triggers=triggers,
                )
            )
    return cards
