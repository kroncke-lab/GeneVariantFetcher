"""Somatic-vs-germline QC for variant records (gold-free, cancer-gene-safe).

For cardiac/Mendelian genes a variant listed in a clinical table is almost
always a constitutional (germline) carrier. Cancer genes (BRCA2, etc.) are
different: their literature mixes germline hereditary cohorts with SOMATIC
tumor-sequencing — tumor biopsies, ctDNA/cfDNA, cell lines, xenografts, loss of
heterozygosity. A somatic tumor variant is NOT a heritable carrier, so counting
it as one corrupts carrier/penetrance numbers.

This module flags each variant record as ``germline`` / ``somatic`` /
``ambiguous`` / ``unknown`` from its own evidence text — no gold standard, no
network — so a cold-start cancer-gene run can MEASURE (and later gate) somatic
contamination. It is **flag-only by default** (it never silently drops counts),
consistent with the project rule that new validators annotate before they
mutate; an explicit ``drop_somatic`` policy is available for a gated pass.

Public API:
    classify_text(text) -> ClassifierHit
    classify_variant(variant, paper_context="") -> ClassifierHit
    annotate_variants(variants, paper_context="", policy="flag") -> dict  # summary

The variant dict shape matches ``pipeline.extraction`` output: each variant may
carry ``source_location`` (str), ``key_quotes`` (list[str]), and
``individual_records`` (list of objects with ``evidence_sentence`` /
``phenotype_details``). All are optional; missing text yields ``unknown``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

# High-precision SOMATIC (tumor-derived, non-heritable) indicators. Chosen to be
# safe as case-insensitive substrings (stems match plurals, e.g. "tumor" ->
# "tumors"/"tumoral"); short ambiguous acronyms are intentionally omitted.
SOMATIC_TERMS: tuple[str, ...] = (
    "somatic",
    "tumor",
    "tumour",
    "ctdna",
    "cfdna",
    "circulating tumor",
    "cell line",
    "cell-line",
    "xenograft",
    "loss of heterozygosity",
    "biopsy",
    "biopsies",
    "tcga",
    "tumor tissue",
    "tumour tissue",
    "tumor specimen",
    "tumor sample",
    "metastas",
    "metastat",
    "ffpe",
    "formalin-fixed",
    "tumor-normal",
    "tumor-only",
    "neoplas",
    "somatic mutation",
    "tumor dna",
)

# GERMLINE / constitutional / heritable-carrier indicators.
GERMLINE_TERMS: tuple[str, ...] = (
    "germline",
    "constitutional",
    "inherited",
    "hereditary",
    "carrier",
    "proband",
    "pedigree",
    "familial",
    "family history",
    "family member",
    "segregat",
    "peripheral blood",
    "blood sample",
    "blood-derived",
    "leukocyte",
    "lymphocyte",
    "saliva",
    "heterozygous carrier",
    "de novo",
    "unaffected relative",
    "autosomal dominant",
    "autosomal recessive",
)

GERMLINE = "germline"
SOMATIC = "somatic"
AMBIGUOUS = "ambiguous"
UNKNOWN = "unknown"

# Labels a downstream gate should treat as "possibly not a heritable carrier".
CONTAMINATION_LABELS = (SOMATIC, AMBIGUOUS)


@dataclass
class ClassifierHit:
    """The somatic/germline verdict for one record plus its supporting terms."""

    label: str
    somatic_terms: list[str]
    germline_terms: list[str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "somatic_terms": self.somatic_terms,
            "germline_terms": self.germline_terms,
        }


def _matches(text_lower: str, terms: Iterable[str]) -> list[str]:
    return sorted({t for t in terms if t in text_lower})


def classify_text(text: str | None) -> ClassifierHit:
    """Classify a blob of evidence text as somatic / germline / ambiguous / unknown.

    Both signals present -> ``ambiguous`` (a paper discussing tumor AND germline
    cohorts); exactly one -> that label; neither -> ``unknown``. The matched
    terms are returned so a human or a downstream gate can adjudicate.
    """
    text_lower = (text or "").lower()
    somatic = _matches(text_lower, SOMATIC_TERMS)
    germline = _matches(text_lower, GERMLINE_TERMS)
    if somatic and germline:
        label = AMBIGUOUS
    elif somatic:
        label = SOMATIC
    elif germline:
        label = GERMLINE
    else:
        label = UNKNOWN
    return ClassifierHit(label=label, somatic_terms=somatic, germline_terms=germline)


def _variant_text(variant: dict[str, Any], paper_context: str = "") -> str:
    """Collect every text field on a variant record that carries cohort context."""
    parts: list[str] = [paper_context or ""]
    sl = variant.get("source_location")
    if isinstance(sl, str):
        parts.append(sl)
    quotes = variant.get("key_quotes")
    if isinstance(quotes, list):
        parts.extend(str(q) for q in quotes)
    records = variant.get("individual_records")
    if isinstance(records, list):
        for rec in records:
            if not isinstance(rec, dict):
                continue
            for key in ("evidence_sentence", "phenotype_details"):
                val = rec.get(key)
                if isinstance(val, str):
                    parts.append(val)
    return "\n".join(p for p in parts if p)


def classify_variant(variant: dict[str, Any], paper_context: str = "") -> ClassifierHit:
    """Classify one variant record using all of its evidence text + paper context."""
    return classify_text(_variant_text(variant, paper_context))


def annotate_variants(
    variants: list[dict[str, Any]],
    paper_context: str = "",
    policy: str = "flag",
) -> dict[str, Any]:
    """Annotate (and optionally filter) variants by somatic/germline status.

    ``policy="flag"`` (default) writes ``variant["somatic_germline_flag"]`` on
    every variant and changes nothing else. ``policy="drop_somatic"`` ALSO
    removes records labelled ``somatic`` (ambiguous/unknown are kept — never
    dropped on a guess). Returns a summary with per-label counts and the
    ``somatic_fraction`` contamination estimate (somatic+ambiguous / total).
    """
    if policy not in ("flag", "drop_somatic"):
        raise ValueError(f"unknown policy: {policy!r}")

    counts = {GERMLINE: 0, SOMATIC: 0, AMBIGUOUS: 0, UNKNOWN: 0}
    kept: list[dict[str, Any]] = []
    dropped = 0
    for variant in variants:
        hit = classify_variant(variant, paper_context)
        variant["somatic_germline_flag"] = hit.to_dict()
        counts[hit.label] += 1
        if policy == "drop_somatic" and hit.label == SOMATIC:
            dropped += 1
            continue
        kept.append(variant)

    if policy == "drop_somatic":
        variants[:] = kept

    total = sum(counts.values())
    contaminated = counts[SOMATIC] + counts[AMBIGUOUS]
    return {
        "policy": policy,
        "total": total,
        "counts": counts,
        "dropped_somatic": dropped,
        "somatic_fraction": (contaminated / total) if total else 0.0,
    }
