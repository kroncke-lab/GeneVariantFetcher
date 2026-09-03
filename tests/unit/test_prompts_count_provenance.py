"""Ensure both extraction prompts request count_provenance.

These checks guard the schema contract that downstream classification
(#5) and the evidence-card validator (#6) depend on. If a future edit
removes the count_provenance section, these tests fail loudly rather
than silently regressing the LLM's output shape.
"""

import re

from pipeline.prompts import (
    COMPACT_EXTRACTION_SYSTEM_PROMPT,
    COMPACT_EXTRACTION_USER_PROMPT,
    EXTRACTION_SYSTEM_PROMPT,
    EXTRACTION_USER_PROMPT,
)

# The contract is over the full instruction payload the model sees, so each
# mode is checked as system + user combined; the split itself is pinned by
# test_system_prompts_are_byte_stable below.
PROMPTS = {
    "COMPACT_EXTRACTION_PROMPT": (
        COMPACT_EXTRACTION_SYSTEM_PROMPT + COMPACT_EXTRACTION_USER_PROMPT
    ),
    "EXTRACTION_PROMPT": EXTRACTION_SYSTEM_PROMPT + EXTRACTION_USER_PROMPT,
}


def test_system_prompts_are_byte_stable():
    """Provider prompt caching requires the system block to be identical
    across papers: no per-paper str.format field may appear in it, and the
    user template must still carry all of them."""
    for name, system in (
        ("COMPACT_EXTRACTION_SYSTEM_PROMPT", COMPACT_EXTRACTION_SYSTEM_PROMPT),
        ("EXTRACTION_SYSTEM_PROMPT", EXTRACTION_SYSTEM_PROMPT),
    ):
        for field in ("{gene_symbol}", "{title}", "{full_text}", "{pmid}"):
            assert field not in system, f"{name} must not vary per paper ({field})"
    for name, user in (
        ("COMPACT_EXTRACTION_USER_PROMPT", COMPACT_EXTRACTION_USER_PROMPT),
        ("EXTRACTION_USER_PROMPT", EXTRACTION_USER_PROMPT),
    ):
        for field in ("{gene_symbol}", "{title}", "{full_text}", "{pmid}"):
            assert field in user, f"{name} lost required field {field}"


REQUIRED_COUNT_TYPES = {
    "per_variant_carrier",
    "family_count",
    "proband_count",
    "cohort_total",
    "screened_N",
    "case",
    "control",
    "unaffected_control",
    "derived_from_patient_rows",
    "closed_variant_partition",
    "unknown",
}


def test_both_prompts_declare_count_provenance_schema():
    for name, prompt in PROMPTS.items():
        assert "count_provenance" in prompt, (
            f"{name} is missing the count_provenance block — #5/#6 depend on it"
        )
        # Field names per logical count
        for label in (
            "carriers_column_label",
            "carriers_count_type",
            "affected_column_label",
            "affected_count_type",
            "unaffected_column_label",
            "unaffected_count_type",
        ):
            assert label in prompt, f"{name} missing required key '{label}'"


def test_both_prompts_enumerate_count_types():
    for name, prompt in PROMPTS.items():
        # Some prompts use "(same enum)" shorthand for the second/third occurrences;
        # the full enum must appear at least once.
        for ct in REQUIRED_COUNT_TYPES:
            assert ct in prompt, (
                f"{name} missing count_type value '{ct}' — schema must list "
                f"all canonical types so the LLM has a closed vocabulary"
            )


def test_both_prompts_have_provenance_null_rule():
    """The rule that says 'if count_type is cohort_total/screened_N/unknown
    AND value is large, leave count NULL' is the protective heuristic
    that makes #4 useful even before #5/#6 ship."""
    pattern = re.compile(
        r"cohort_total.+screened_N.+unknown.+large.+NULL",
        re.DOTALL,
    )
    for name, prompt in PROMPTS.items():
        assert pattern.search(prompt), (
            f"{name} missing the 'leave count NULL when count_type is "
            f"cohort_total/screened_N/unknown and value is large' rule"
        )


def test_both_prompts_declare_fact_level_provenance_schema():
    for name, prompt in PROMPTS.items():
        assert "fact_provenance" in prompt, (
            f"{name} is missing the fact_provenance block required for "
            "row/paragraph-level evidence"
        )
        for label in (
            "fact_type",
            "fact_value",
            "source_table",
            "source_row",
            "source_column",
            "source_section",
            "source_paragraph",
            "evidence_quote",
        ):
            assert label in prompt, f"{name} missing required key '{label}'"


def test_both_prompts_gate_patient_row_phenotype_derivation():
    for name, prompt in PROMPTS.items():
        for phrase in (
            "COMPLETE, variant-specific table",
            "MUST populate affected, unaffected",
            "require ALL selected cells to be explicitly negative",
            "phenotype_derivation",
            "derived_from_patient_rows",
            "uncertain_count",
            "unaffected = total - affected",
            "predicate_tallies",
            "patient_row_phenotype_v2",
        ):
            assert phrase in prompt, f"{name} missing derived-count guard {phrase!r}"


def test_both_prompts_separate_diagnosis_from_asymptomatic_status():
    for name, prompt in PROMPTS.items():
        assert "can be affected while asymptomatic" in prompt, name
        assert "closed_variant_partition" in prompt, name
