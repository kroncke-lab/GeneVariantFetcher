"""Gold-free tests for copied affected/unaffected refusal.

Examples use made-up family sizes and generic notations. They pin the
contract (do not copy carriers onto affected; do not derive una=0), not
any gold-120 PMID.
"""

from __future__ import annotations

from pipeline.phenotype_count_guard import (
    apply_phenotype_count_guard,
    phenotype_fields_to_clear,
    read_phenotype_counts,
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


def test_explicit_carrier_denominator_beats_legacy_patient_subset_mirror():
    variant = {
        "patients": {"count": 2},
        "penetrance_data": {
            "total_carriers_observed": 4,
            "affected_count": 2,
            "unaffected_count": 2,
            "uncertain_count": 0,
        },
    }

    assert read_phenotype_counts(variant) == {
        "carriers": 4,
        "affected": 2,
        "unaffected": 2,
        "uncertain": 0,
    }
    assert phenotype_fields_to_clear(variant) == []


def test_closed_partition_is_left_alone():
    variant = {"carriers": 4, "affected": 2, "unaffected": 2}
    assert phenotype_fields_to_clear(variant) == []


def test_closed_three_way_derived_partition_is_left_alone():
    variant = {
        "carriers": 185,
        "affected": 97,
        "unaffected": 62,
        "uncertain": 26,
        "count_provenance": {
            "affected_count_type": "derived_from_patient_rows",
            "unaffected_count_type": "derived_from_patient_rows",
        },
    }
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


def test_unsourced_zero_unaffected_with_null_total_is_refused():
    variant = {"carriers": None, "affected": None, "unaffected": 0}

    assert [
        (item.field, item.reason) for item in phenotype_fields_to_clear(variant)
    ] == [("unaffected", "unsourced_zero_without_closed_cohort")]


def test_unsourced_zero_unaffected_with_missing_affected_is_refused():
    variant = {"carriers": 4, "affected": None, "unaffected": 0}

    assert [
        (item.field, item.reason) for item in phenotype_fields_to_clear(variant)
    ] == [("unaffected", "unsourced_zero_without_closed_cohort")]


def test_llm_verdict_alone_does_not_source_an_open_set_zero():
    variant = {
        "carriers": None,
        "affected": None,
        "unaffected": 0,
        "claim_verification": {"field_verdicts": {"unaffected": "directly_supported"}},
    }

    assert [
        (item.field, item.reason) for item in phenotype_fields_to_clear(variant)
    ] == [("unaffected", "unsourced_zero_without_closed_cohort")]


def test_guard_keeps_equal_affected_count_with_verified_same_population_source():
    source = (
        "The current study population consisted of 30 carriers of the D1790G "
        "SCN5A mutation.\n"
        "The study population consisted of 30 LQT3 patients."
    )
    variant = {
        "source_notation": "D1790G",
        "carriers": 30,
        "affected": 30,
        "unaffected": None,
        "fact_provenance": [
            {
                "fact_type": "total_carriers_observed",
                "fact_value": "30",
                "evidence_quote": (
                    "The current study population consisted of 30 carriers of "
                    "the D1790G SCN5A mutation"
                ),
            },
            {
                "fact_type": "affected_count",
                "fact_value": "30",
                "evidence_quote": "The study population consisted of 30 LQT3 patients",
            },
        ],
    }

    result = apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert result.cleared == 0
    assert variant["affected"] == 30
    assert variant["source_verified_claims"]["affected"]["value"] == 30


def test_guard_promotes_missing_affected_for_same_closed_carrier_population():
    source = (
        "The current study population consisted of 30 carriers of the D1790G "
        "SCN5A mutation who were members of four families.\n"
        "The study population consisted of 30 LQT3 patients. Prior to enrolment, "
        "27 patients were asymptomatic and 3 patients were symptomatic.\n"
        "The study population comprised 30 D1790G carriers who were treated "
        "with flecainide and followed longitudinally."
    )
    variant = {
        "protein_notation": "p.Asp1790Gly",
        "source_notation": "D1790G",
        "patients": {"count": 30},
        "penetrance_data": {
            "total_carriers_observed": 30,
            "affected_count": None,
            "unaffected_count": None,
            "uncertain_count": 0,
        },
    }

    result = apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert result.cleared == 0
    assert variant["penetrance_data"]["affected_count"] == 30
    assert variant["penetrance_data"]["unaffected_count"] is None
    assert variant["source_bound_phenotype_promotion"]["fields_filled"] == ["affected"]
    assert variant["source_verified_claims"]["affected"]["value"] == 30


def test_guard_keeps_existing_equal_n_when_source_detector_validates_it():
    """A detected closed cohort verifies a present value, not only a null.

    This covers extraction records whose fact-provenance quote is not an exact
    source substring even though the complete locked source independently
    binds the same N carriers and N disease patients.
    """
    source = (
        "The current study population consisted of 30 carriers of the D1790G "
        "SCN5A mutation who were members of four families.\n"
        "The study population consisted of 30 LQT3 patients. Prior to enrolment, "
        "27 patients were asymptomatic and 3 patients were symptomatic."
    )
    variant = {
        "protein_notation": "p.Asp1790Gly",
        "source_notation": "D1790G",
        "patients": {"count": 30},
        "penetrance_data": {
            "total_carriers_observed": 30,
            "affected_count": 30,
            "unaffected_count": 0,
            "uncertain_count": 0,
        },
        "count_provenance": {
            "carriers_count_type": "per_variant_carrier",
            "affected_count_type": "per_variant_carrier",
            "unaffected_count_type": None,
        },
    }

    result = apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert result.cleared == 1
    assert variant["penetrance_data"]["affected_count"] == 30
    # The source does not explicitly report a zero; it remains refused.
    assert variant["penetrance_data"]["unaffected_count"] is None
    assert variant["source_verified_claims"]["affected"] == {
        "value": 30,
        "method": "identity_disease_coreference_equal_n",
        "quote": (
            "The current study population consisted of 30 carriers of the "
            "D1790G SCN5A mutation who were members of four families. || The "
            "study population consisted of 30 LQT3 patients."
        ),
        "promoted": False,
    }
    assert "source_bound_phenotype_promotion" not in variant


def test_guard_replaces_closed_symptom_split_for_verified_disease_cohort():
    source = (
        "The current study population consisted of 30 carriers of the D1790G "
        "SCN5A mutation who were members of four families.\n"
        "The study population consisted of 30 LQT3 patients. Prior to enrolment, "
        "27 patients were asymptomatic and 3 patients were symptomatic."
    )
    variant = {
        "protein_notation": "p.Asp1790Gly",
        "source_notation": "D1790G",
        "patients": {"count": 30},
        "penetrance_data": {
            "total_carriers_observed": 30,
            "affected_count": 3,
            "unaffected_count": 27,
            "uncertain_count": 0,
        },
    }

    result = apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert result.cleared == 0
    assert variant["penetrance_data"]["affected_count"] == 30
    assert variant["penetrance_data"]["unaffected_count"] is None
    assert variant["source_verified_claims"]["affected"] == {
        "value": 30,
        "method": "identity_disease_coreference_equal_n",
        "quote": (
            "The current study population consisted of 30 carriers of the "
            "D1790G SCN5A mutation who were members of four families. || The "
            "study population consisted of 30 LQT3 patients."
        ),
        "promoted": False,
        "corrected": True,
    }
    assert variant["source_bound_phenotype_correction"]["fields_corrected"] == {
        "affected": {"old": 3, "new": 30},
        "unaffected": {"old": 27, "new": None},
    }


def test_guard_promotes_lqt3_when_disease_context_is_acronym_only():
    source = (
        "The current study population consisted of 30 carriers of the D1790G "
        "SCN5A mutation.\n"
        "The study population consisted of 30 LQT3 patients."
    )
    variant = {
        "protein_notation": "p.Asp1790Gly",
        "penetrance_data": {
            "total_carriers_observed": 30,
            "affected_count": None,
            "unaffected_count": None,
            "uncertain_count": 0,
        },
    }

    apply_phenotype_count_guard([variant], source_text=source, disease="LQTS")

    assert variant["penetrance_data"]["affected_count"] == 30
    assert variant["source_bound_phenotype_promotion"]["fields_filled"] == ["affected"]


def test_guard_promotes_nonzero_remainder_from_explicit_closed_partition():
    source = (
        "Overall, 24 out of 29 KCNQ1 dup12 mutation carriers were affected "
        "by an LQT syndrome, suggesting an incomplete penetrance of about 80%."
    )
    variant = {
        "protein_notation": "p.R360_Q361dupQKQR",
        "source_notation": "KCNQ1 dup12; p.R360_Q361dupQKQR",
        "patients": {"count": 29},
        "penetrance_data": {
            "total_carriers_observed": 29,
            "affected_count": 24,
            "unaffected_count": None,
            "uncertain_count": 0,
        },
    }

    result = apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert result.cleared == 0
    assert variant["penetrance_data"] == {
        "total_carriers_observed": 29,
        "affected_count": 24,
        "unaffected_count": 5,
        "uncertain_count": 0,
    }
    assert variant["source_bound_phenotype_promotion"]["fields_filled"] == [
        "unaffected"
    ]
    assert variant["count_provenance"]["unaffected_count_type"] == (
        "closed_variant_partition"
    )
    assert variant["count_provenance"]["unaffected_source"] == (
        "source_bound_phenotype_v1"
    )


def test_guard_recovers_explicit_carrier_of_whom_partition():
    source = (
        "We found 8 carriers of the new RYR2 C2277R variant, 7 of whom "
        "exhibited the CPVT phenotype according to EST results."
    )
    variant = {
        "protein_notation": "p.Cys2277Arg",
        "source_notation": "C2277R",
        "penetrance_data": {
            "total_carriers_observed": 8,
            "affected_count": None,
            "unaffected_count": None,
            "uncertain_count": 0,
        },
    }

    result = apply_phenotype_count_guard([variant], source_text=source, disease="CPVT")

    assert result.cleared == 0
    assert variant["penetrance_data"] == {
        "total_carriers_observed": 8,
        "affected_count": 7,
        "unaffected_count": 1,
        "uncertain_count": 0,
    }
    provenance = variant["count_provenance"]
    assert provenance["affected_source"] == "source_bound_phenotype_v1"
    assert provenance["unaffected_source"] == "source_bound_phenotype_v1"


def test_guard_recovers_all_carriers_with_disease_defining_manifestation():
    source = (
        "The 6 Y652X carriers all manifested prolonged QTc intervals. "
        "The penetrance in this family was therefore assumed to be 100%."
    )
    variant = {
        "protein_notation": "p.Tyr652Ter",
        "source_notation": "Y652X",
        "penetrance_data": {
            "total_carriers_observed": 6,
            "affected_count": None,
            "unaffected_count": None,
            "uncertain_count": 0,
        },
    }

    result = apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert result.cleared == 0
    assert variant["penetrance_data"]["affected_count"] == 6
    assert variant["penetrance_data"]["unaffected_count"] is None
    provenance = variant["count_provenance"]
    assert provenance["affected_count_type"] == "closed_variant_partition"
    assert provenance["affected_source"] == "source_bound_phenotype_v1"


def test_model_declared_closed_partition_does_not_source_zero():
    variant = {
        "carriers": 6,
        "affected": 0,
        "unaffected": 6,
        "count_provenance": {"affected_count_type": "closed_variant_partition"},
    }

    assert [
        (item.field, item.reason) for item in phenotype_fields_to_clear(variant)
    ] == [("affected", "unsourced_zero_affected")]


def test_closed_partition_can_fill_all_missing_nonzero_counts():
    source = "Seven out of 10 p.Arg12Trp carriers were affected by long QT syndrome."
    variant = {
        "protein_notation": "p.Arg12Trp",
        "penetrance_data": {
            "total_carriers_observed": None,
            "affected_count": None,
            "unaffected_count": None,
            "uncertain_count": None,
        },
    }

    # Word numerals are deliberately not parsed: exact integer provenance is
    # required for a deterministic promotion.
    apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )
    assert variant["penetrance_data"]["total_carriers_observed"] is None

    source = "7 out of 10 p.Arg12Trp carriers were affected by long QT syndrome."
    apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )
    assert variant["penetrance_data"]["total_carriers_observed"] == 10
    assert variant["penetrance_data"]["affected_count"] == 7
    assert variant["penetrance_data"]["unaffected_count"] == 3


def test_closed_partition_does_not_turn_symptom_subset_into_disease_split():
    source = (
        "6 out of 20 p.Arg12Trp carriers had cardiac events after stopping "
        "therapy; all were long QT syndrome patients."
    )
    variant = {
        "protein_notation": "p.Arg12Trp",
        "carriers": 20,
        "affected": None,
        "unaffected": None,
    }

    apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert variant["affected"] is None
    assert variant["unaffected"] is None
    assert "source_bound_phenotype_promotion" not in variant


def test_closed_partition_refuses_assessed_subset_and_conflicting_integer():
    source = (
        "7 out of 10 evaluated p.Arg12Trp carriers were affected by long QT syndrome."
    )
    subset = {
        "protein_notation": "p.Arg12Trp",
        "carriers": None,
        "affected": None,
        "unaffected": None,
    }
    apply_phenotype_count_guard(
        [subset], source_text=source, disease="Long QT syndrome"
    )
    assert subset["carriers"] is None

    source = "7 out of 10 p.Arg12Trp carriers were affected by long QT syndrome."
    conflict = {
        "protein_notation": "p.Arg12Trp",
        "carriers": 11,
        "affected": None,
        "unaffected": None,
    }
    apply_phenotype_count_guard(
        [conflict], source_text=source, disease="Long QT syndrome"
    )
    assert conflict == {
        "protein_notation": "p.Arg12Trp",
        "carriers": 11,
        "affected": None,
        "unaffected": None,
    }


def test_closed_partition_never_manufactures_zero_remainder():
    source = "10 out of 10 p.Arg12Trp carriers were affected by long QT syndrome."
    variant = {
        "protein_notation": "p.Arg12Trp",
        "carriers": 10,
        "affected": None,
        "unaffected": None,
    }

    apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert variant["affected"] == 10
    assert variant["unaffected"] is None


def test_equal_n_across_unlinked_sentences_is_not_same_population():
    source = (
        "We enrolled 30 LQT3 patients from the referral clinic. "
        "Separately, 30 carriers of p.Asp1790Gly were present in a registry."
    )
    variant = {
        "protein_notation": "p.Asp1790Gly",
        "carriers": 30,
        "affected": None,
        "unaffected": None,
    }

    apply_phenotype_count_guard(
        [variant], source_text=source, disease="Long QT syndrome"
    )

    assert variant["affected"] is None
