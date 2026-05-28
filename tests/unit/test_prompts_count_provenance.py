"""Ensure both extraction prompts request count_provenance.

These checks guard the schema contract that downstream classification
(#5) and the evidence-card validator (#6) depend on. If a future edit
removes the count_provenance section, these tests fail loudly rather
than silently regressing the LLM's output shape.
"""

import re

from pipeline.prompts import COMPACT_EXTRACTION_PROMPT, EXTRACTION_PROMPT

PROMPTS = {
    "COMPACT_EXTRACTION_PROMPT": COMPACT_EXTRACTION_PROMPT,
    "EXTRACTION_PROMPT": EXTRACTION_PROMPT,
}

REQUIRED_COUNT_TYPES = {
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


def test_both_prompts_declare_count_provenance_schema():
    for name, prompt in PROMPTS.items():
        assert (
            "count_provenance" in prompt
        ), f"{name} is missing the count_provenance block — #5/#6 depend on it"
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
