"""Refuse phenotype counts that were copied from the carrier total.

The extraction contract says a wrong affected/unaffected number is worse than
a missing one, and forbids two fabrications that keep showing up as gold-120
MAE:

* ``affected = carriers`` when the paper never counted a phenotype split
* ``unaffected = 0`` derived from that copy ("everyone is affected")
* an ``affected`` integer that cannot be true, because the emitted triple does
  not close (``affected + unaffected != carriers``)
* a bare ``affected = 0`` with no sourced phenotype column behind it

Those integers alone do not name a paper or a variant. A real split
(``affected != carriers``, or a positive ``unaffected``) is left alone. A
sourced unaffected column that actually reports 0 is left alone. This module
is gold-free so it transfers to genes with no curated answer.

Public API::

    sanitize_copied_phenotype(variant) -> list[PhenotypeClear]
    apply_phenotype_count_guard(variants) -> GuardSummary

Both accept the nested extraction shape (``patients`` / ``penetrance_data``)
and the flat prediction / figure-reader shape (``carriers`` / ``affected`` /
``unaffected``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from pipeline.count_outlier_guard import COUNT_FIELDS, _clear_field, _read_field
from utils.source_layers import source_layer_tokens

# Copying a lone carrier onto affected is often a true one-proband case
# report. Copying N>=2 is the "default the whole family to affected" pattern
# that produced the largest gold-120 phenotype errors (figure pedigrees).
COPIED_AFFECTED_MIN_CARRIERS = 2

# Labels that actually count a phenotype class, not the carrier total.
_AFFECTED_LABEL_MARKERS = (
    "affected",
    "symptomatic",
    "phenotype",
    "disease",
    "lqts",
    "lqt",
    "brs",
    "cpvt",
    "case",
    "cases",
    "proband",
    "probands",
)
_UNAFFECTED_LABEL_MARKERS = (
    "unaffected",
    "asymptomatic",
    "control",
    "controls",
    "healthy",
    "normal",
    "silent",
)

# count_type values that do not source an integer; "per_variant_carrier" is
# the LLM default and is not evidence of a phenotype column.
_UNSOURCED_TYPES = frozenset({"", "unknown", "per_variant_carrier"})


@dataclass(frozen=True)
class PhenotypeClear:
    """One field nulled because it was a copied/derived phenotype count."""

    field: str
    value: int
    reason: str


@dataclass
class GuardSummary:
    flagged: int = 0
    cleared: int = 0
    annotations: list[dict[str, Any]] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "flagged": self.flagged,
            "cleared": self.cleared,
            "annotations": list(self.annotations),
        }


def _is_count(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _coerce_count(value: Any) -> Optional[int]:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 0 else None
    if isinstance(value, float) and value.is_integer() and value >= 0:
        return int(value)
    if isinstance(value, str):
        cleaned = value.strip().replace(",", "")
        if cleaned.isdigit():
            return int(cleaned)
    return None


def _provenance(variant: dict[str, Any]) -> dict[str, Any]:
    raw = variant.get("count_provenance")
    return raw if isinstance(raw, dict) else {}


def _label_looks_like(label: Any, markers: Iterable[str]) -> bool:
    text = str(label or "").strip().lower()
    if not text:
        return False
    return any(marker in text for marker in markers)


def _count_type(provenance: dict[str, Any], field: str) -> str:
    key = {
        "affected": "affected_count_type",
        "unaffected": "unaffected_count_type",
        "carriers": "carriers_count_type",
    }.get(field, "")
    return str(provenance.get(key) or "").strip().lower()


def _column_label(provenance: dict[str, Any], field: str) -> str:
    key = {
        "affected": "affected_column_label",
        "unaffected": "unaffected_column_label",
        "carriers": "carriers_column_label",
    }.get(field, "")
    return str(provenance.get(key) or "").strip()


def _phenotype_is_sourced(variant: dict[str, Any], field: str) -> bool:
    """True only when a *phenotype* column, not the carrier total, was named."""
    provenance = _provenance(variant)
    label = _column_label(provenance, field)
    carrier_label = _column_label(provenance, "carriers")
    markers = (
        _AFFECTED_LABEL_MARKERS if field == "affected" else _UNAFFECTED_LABEL_MARKERS
    )
    if label and label != carrier_label and _label_looks_like(label, markers):
        return True
    declared = _count_type(provenance, field)
    if declared and declared not in _UNSOURCED_TYPES:
        if field == "affected" and declared in {"case", "proband_count"}:
            return True
        if field == "unaffected" and declared in {
            "control",
            "unaffected_control",
        }:
            return True
    return False


def read_phenotype_counts(variant: dict[str, Any]) -> dict[str, Optional[int]]:
    """Read carriers/affected/unaffected from nested or flat variant dicts."""
    if not isinstance(variant, dict):
        return {"carriers": None, "affected": None, "unaffected": None}

    nested = {
        field: _read_field(variant, paths) for field, paths in COUNT_FIELDS.items()
    }
    if any(value is not None for value in nested.values()):
        return nested

    return {
        "carriers": _coerce_count(variant.get("carriers")),
        "affected": _coerce_count(variant.get("affected")),
        "unaffected": _coerce_count(variant.get("unaffected")),
    }


def _write_flat(variant: dict[str, Any], field: str) -> None:
    if field in variant:
        variant[field] = None


def _layer_tokens(variant: dict[str, Any]) -> set[str]:
    return source_layer_tokens(variant.get("source_layer"))


def _is_copied_phenotype(
    carriers: Optional[int], affected: Optional[int], unaffected: Optional[int]
) -> bool:
    return (
        _is_count(carriers)
        and carriers > 0
        and _is_count(affected)
        and affected == carriers
        and (unaffected is None or unaffected == 0)
    )


def phenotype_fields_to_clear(variant: dict[str, Any]) -> list[PhenotypeClear]:
    """Return phenotype fields that look copied or derived, not counted."""
    if not isinstance(variant, dict):
        return []

    counts = read_phenotype_counts(variant)
    carriers = counts["carriers"]
    affected = counts["affected"]
    unaffected = counts["unaffected"]
    copied = _is_copied_phenotype(carriers, affected, unaffected)
    figure = "figure" in _layer_tokens(variant)
    aff_sourced = _phenotype_is_sourced(variant, "affected")
    una_sourced = _phenotype_is_sourced(variant, "unaffected")

    clears: list[PhenotypeClear] = []
    seen: set[str] = set()

    def add(field_name: str, value: Optional[int], reason: str) -> None:
        if field_name in seen or not _is_count(value):
            return
        seen.add(field_name)
        clears.append(PhenotypeClear(field=field_name, value=value, reason=reason))

    # Family/cohort copy: N>=2 carriers dumped onto affected with no
    # split. The implied una=0 rides along; a lone-proband 1/1/0 case
    # report is not this pattern and is left alone.
    if (
        copied
        and carriers is not None
        and carriers >= COPIED_AFFECTED_MIN_CARRIERS
        and not aff_sourced
    ):
        add("affected", affected, "copied_carriers_onto_affected")
        if unaffected == 0 and not una_sourced:
            add("unaffected", unaffected, "implied_unaffected_zero")

    # Figure pedigrees: vision copies every symbol onto affected. Fail
    # closed on any N — a one-symbol pedigree still has no counted split.
    if copied and figure and not aff_sourced:
        add("affected", affected, "figure_copied_phenotype")
        if unaffected == 0 and not una_sourced:
            add("unaffected", unaffected, "figure_copied_phenotype")

    # Incomplete figure phenotype: an affected integer with no carrier
    # total cannot be a counted split. Keep carriers if present.
    if (
        figure
        and _is_count(affected)
        and carriers is None
        and (unaffected is None or unaffected == 0)
        and not aff_sourced
    ):
        add("affected", affected, "figure_incomplete_phenotype")

    # Partition that does not close. When all three integers are emitted for
    # one cohort, ``affected + unaffected`` must equal ``carriers``; the
    # extraction contract already asks the model to run this self-check. A
    # triple that fails it contains at least one wrong number, and on the
    # gold-120 lock the wrong one is ``affected`` every time (3/3), while the
    # companion ``unaffected`` is exact on two of those three rows. So refuse
    # only ``affected`` -- nulling the pair destroys counted values that the
    # paper really did report.
    #
    # Underfill is not the same defect as overflow, but both are refused here:
    # ``carriers = affected + unaffected + unassessed`` has no ``unassessed``
    # slot in this schema, so a short partition is indistinguishable from a
    # miscount and must not be published as a phenotype split.
    if (
        _is_count(carriers)
        and _is_count(affected)
        and _is_count(unaffected)
        and affected + unaffected != carriers
        and not aff_sourced
    ):
        add("affected", affected, "partition_does_not_close")

    # A counted zero is a positive clinical claim ("this family is entirely
    # non-penetrant"), not an abstention, so it needs a real phenotype column
    # behind it. ``unaffected = 0`` is deliberately NOT refused here: it is the
    # ordinary single-proband case-report shape and clearing it was measured
    # net-negative.
    if _is_count(affected) and affected == 0 and not aff_sourced:
        add("affected", affected, "unsourced_zero_affected")

    return clears


def _apply_clears(variant: dict[str, Any], clears: list[PhenotypeClear]) -> None:
    for item in clears:
        paths = COUNT_FIELDS.get(item.field)
        if paths:
            _clear_field(variant, paths)
        _write_flat(variant, item.field)


def sanitize_copied_phenotype(variant: dict[str, Any]) -> list[PhenotypeClear]:
    """Null copied phenotype fields on one variant. Returns what was cleared."""
    clears = phenotype_fields_to_clear(variant)
    if not clears:
        return []
    _apply_clears(variant, clears)
    flags = variant.setdefault("phenotype_count_flags", {})
    for item in clears:
        flags[item.field] = {
            "raw": item.value,
            "reason": item.reason,
            "policy": "clear",
        }
    return clears


def apply_phenotype_count_guard(
    variants: list[dict[str, Any]],
) -> GuardSummary:
    """Sanitize a paper's variants in place. Gold-free; always-on contract."""
    summary = GuardSummary()
    if not isinstance(variants, list):
        return summary
    for idx, variant in enumerate(variants):
        if not isinstance(variant, dict):
            continue
        clears = sanitize_copied_phenotype(variant)
        if not clears:
            continue
        summary.flagged += len(clears)
        summary.cleared += len(clears)
        for item in clears:
            summary.annotations.append(
                {
                    "variant_index": idx,
                    "field": item.field,
                    "value": item.value,
                    "reason": item.reason,
                }
            )
    return summary
