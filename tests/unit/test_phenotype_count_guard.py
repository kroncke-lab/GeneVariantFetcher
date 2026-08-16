"""Gold-free tests for copied affected/unaffected refusal.

Examples use made-up family sizes and generic notations. They pin the
contract (do not copy carriers onto affected; do not derive una=0), not
any gold-120 PMID.
"""

from __future__ import annotations

from pipeline.phenotype_count_guard import (
    apply_phenotype_count_guard,
    phenotype_fields_to_clear,
    sanitize_copied_phenotype,
)
from pipeline.steps import _apply_phenotype_count_guard


def test_real_phenotype_split_is_left_alone():
    variant = {
        "patients": {"count": 13},
        "penetrance_data": {"affected_count": 6, "unaffected_count": 7},
        "source_layer": "llm_text",
    }
    assert phenotype_fields_to_clear(variant) == []
    sanitize_copied_phenotype(variant)
    assert variant["penetrance_data"]["affected_count"] == 6
    assert variant["penetrance_data"]["unaffected_count"] == 7


def test_family_copy_clears_both_phenotype_fields():
    """N>=2 carriers dumped onto affected with implied una=0 is not a split."""
    variant = {
        "patients": {"count": 7},
        "penetrance_data": {
            "total_carriers_observed": 7,
            "affected_count": 7,
            "unaffected_count": 0,
        },
        "source_layer": "llm_text",
    }
    reasons = {item.field: item.reason for item in sanitize_copied_phenotype(variant)}
    assert reasons["affected"] == "copied_carriers_onto_affected"
    assert reasons["unaffected"] == "implied_unaffected_zero"
    assert variant["patients"]["count"] == 7
    assert variant["penetrance_data"]["total_carriers_observed"] == 7
    assert variant["penetrance_data"]["affected_count"] is None
    assert variant["penetrance_data"]["unaffected_count"] is None


def test_figure_copy_clears_even_for_a_single_carrier():
    variant = {
        "carriers": 1,
        "affected": 1,
        "unaffected": 0,
        "source_layer": "figure",
    }
    reasons = {item.field: item.reason for item in sanitize_copied_phenotype(variant)}
    assert reasons["affected"] == "figure_copied_phenotype"
    assert variant["carriers"] == 1
    assert variant["affected"] is None
    assert variant["unaffected"] is None


def test_single_text_carrier_is_left_alone():
    """A one-proband 1/1/0 case report is not the family-copy pattern."""
    variant = {
        "patients": {"count": 1},
        "penetrance_data": {"affected_count": 1, "unaffected_count": 0},
        "source_layer": "llm_text",
    }
    assert phenotype_fields_to_clear(variant) == []


def test_figure_affected_without_carriers_is_incomplete():
    variant = {
        "carriers": None,
        "affected": 7,
        "unaffected": None,
        "source_layer": "figure",
    }
    clears = sanitize_copied_phenotype(variant)
    assert [item.field for item in clears] == ["affected"]
    assert variant["affected"] is None


def test_sourced_unaffected_zero_is_kept():
    variant = {
        "patients": {"count": 4},
        "penetrance_data": {"affected_count": 4, "unaffected_count": 0},
        "count_provenance": {
            "carriers_column_label": "N carriers",
            "affected_column_label": "LQTS",
            "unaffected_column_label": "Unaffected",
            "affected_count_type": "case",
            "unaffected_count_type": "unaffected_control",
        },
        "source_layer": "llm_table",
    }
    assert phenotype_fields_to_clear(variant) == []


def test_carrier_column_relabelled_as_affected_is_not_sourced():
    variant = {
        "patients": {"count": 5},
        "penetrance_data": {"affected_count": 5, "unaffected_count": 0},
        "count_provenance": {
            "carriers_column_label": "N",
            "affected_column_label": "N",
            "affected_count_type": "per_variant_carrier",
            "unaffected_count_type": "per_variant_carrier",
        },
        "source_layer": "llm_table",
    }
    fields = {item.field for item in phenotype_fields_to_clear(variant)}
    assert fields == {"affected", "unaffected"}


def test_apply_guard_is_idempotent_and_records_flags():
    variants = [
        {
            "patients": {"count": 8},
            "penetrance_data": {"affected_count": 8, "unaffected_count": 0},
            "source_layer": "llm_table",
        }
    ]
    first = apply_phenotype_count_guard(variants)
    second = apply_phenotype_count_guard(variants)
    assert first.cleared == 2
    assert second.cleared == 0
    assert "affected" in variants[0]["phenotype_count_flags"]
    assert variants[0]["phenotype_count_flags"]["affected"]["raw"] == 8


def test_steps_hook_runs_when_hygiene_policies_are_off():
    data = {
        "variants": [
            {
                "patients": {"count": 6},
                "penetrance_data": {"affected_count": 6, "unaffected_count": 0},
                "source_layer": "figure",
            }
        ],
        "extraction_metadata": {},
    }
    _apply_phenotype_count_guard(data)
    variant = data["variants"][0]
    assert variant["penetrance_data"]["affected_count"] is None
    assert variant["penetrance_data"]["unaffected_count"] is None
    assert variant["patients"]["count"] == 6
    assert data["extraction_metadata"]["phenotype_count_guard"]["cleared"] == 2


def test_partition_that_does_not_close_refuses_affected_only():
    """affected + unaffected must equal carriers when all three are emitted.

    The companion ``unaffected`` is often the correct half of a bad triple, so
    only ``affected`` is refused; nulling the pair destroys counted values the
    paper really did report.
    """
    variant = {
        "carriers": 31,
        "affected": 11,
        "unaffected": 18,
        "source_layer": "llm_text",
    }
    clears = phenotype_fields_to_clear(variant)
    assert [(item.field, item.reason) for item in clears] == [
        ("affected", "partition_does_not_close")
    ]


def test_partition_underfill_is_also_refused():
    """A short partition is indistinguishable from a miscount.

    There is no ``unassessed`` slot in this schema, so 2 + 1 != 4 cannot be
    published as a phenotype split.
    """
    variant = {"carriers": 4, "affected": 2, "unaffected": 1}
    assert [item.reason for item in phenotype_fields_to_clear(variant)] == [
        "partition_does_not_close"
    ]


def test_closed_partition_is_left_alone():
    variant = {"carriers": 4, "affected": 2, "unaffected": 2}
    assert phenotype_fields_to_clear(variant) == []


def test_partition_check_skips_rows_that_abstain_on_a_field():
    """Leaving people unassessed (unaffected null) is not a failed partition."""
    variant = {"carriers": 3, "affected": 2, "unaffected": None}
    assert phenotype_fields_to_clear(variant) == []


def test_sourced_affected_column_survives_a_failed_partition():
    variant = {
        "carriers": 31,
        "affected": 11,
        "unaffected": 18,
        "count_provenance": {
            "carriers_column_label": "Carriers",
            "affected_column_label": "Affected (LQTS)",
        },
    }
    assert phenotype_fields_to_clear(variant) == []


def test_unsourced_zero_affected_is_refused():
    """A counted zero is a positive claim and needs a phenotype column."""
    variant = {"carriers": 1, "affected": 0, "unaffected": 1}
    assert [
        (item.field, item.reason) for item in phenotype_fields_to_clear(variant)
    ] == [("affected", "unsourced_zero_affected")]


def test_sourced_zero_affected_is_kept():
    variant = {
        "carriers": 6,
        "affected": 0,
        "unaffected": 6,
        "count_provenance": {
            "carriers_column_label": "Carriers",
            "affected_column_label": "Symptomatic",
        },
    }
    assert phenotype_fields_to_clear(variant) == []


def test_zero_unaffected_is_not_refused_by_the_zero_rule():
    """``unaffected=0`` is the ordinary single-proband case-report shape.

    Clearing it was measured net-negative, so the zero rule is deliberately
    asymmetric and applies to ``affected`` only.
    """
    variant = {"carriers": 1, "affected": 1, "unaffected": 0}
    assert phenotype_fields_to_clear(variant) == []
