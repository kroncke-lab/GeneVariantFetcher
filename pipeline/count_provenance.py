"""Code-owned provenance markers for deterministic phenotype count lanes.

These marker values are deliberately absent from the extraction schema. Model
output is scrubbed before verification and again at the common extraction
boundary, so downstream trust decisions can use the markers as evidence that a
Python audit actually ran rather than as an LLM-authored assertion.
"""

from __future__ import annotations

from typing import Any


PATIENT_ROW_PHENOTYPE_SOURCE = "patient_row_phenotype_v2"
SOURCE_BOUND_PHENOTYPE_SOURCE = "source_bound_phenotype_v1"

CODE_OWNED_PHENOTYPE_SOURCES = frozenset(
    {PATIENT_ROW_PHENOTYPE_SOURCE, SOURCE_BOUND_PHENOTYPE_SOURCE}
)


def strip_model_code_owned_phenotype_sources(
    extracted_data: dict[str, Any],
) -> int:
    """Remove code-owned phenotype source stamps from untrusted model output."""

    removed = 0
    variants = extracted_data.get("variants")
    if not isinstance(variants, list):
        return removed
    for variant in variants:
        if not isinstance(variant, dict):
            continue
        provenance = variant.get("count_provenance")
        if not isinstance(provenance, dict):
            continue
        for key in ("affected_source", "unaffected_source"):
            source = str(provenance.get(key) or "").strip().lower()
            if source in CODE_OWNED_PHENOTYPE_SOURCES:
                provenance.pop(key, None)
                removed += 1
    return removed
