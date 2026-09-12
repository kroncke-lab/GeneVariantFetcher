"""Focused count-ownership and canonical-identity checks for the new union."""

import numpy as np
import pandas as pd
import pytest

from rebuild_union import (
    build_union,
    fit_empirical,
    normalized_protein,
    posterior_table,
)


TX = "ENST00000403799.8"


def population_row(identity="7-100-A-T", **changes):
    row = {
        "gene": "GCK",
        "variant_id": identity,
        "population_eligible": True,
        "canonical_transcript_id": TX,
        "transcript_id": TX,
        "protein_key": "A10V",
        "canonical_wt_status": "match",
        "hgvsc": "c.29C>T",
        "vclass": "missense",
        "aa_pos": 10,
        "aa_ref": "A",
        "aa_alt": "V",
        "coding_or_splice": True,
        "gnomad_carriers": 7,
        "joint_ac": 8,
        "joint_hom": 1,
    }
    row.update(changes)
    return row


def clinical_row(key="A10V", **changes):
    row = {
        "gene": "GCK",
        "key": key,
        "vclass": "missense",
        "aa_pos": 10,
        "aa_ref": "A",
        "aa_alt": "V",
        "canonical_wt_status": "match",
        "canonical_point_key": "A10V",
        "protein_point_join_allowed": True,
        "identity_quarantine_recommended": False,
        "identity_quarantine_reason": "",
        "cdna": "",
        "feature_transcript": TX,
        "affected_literature": 4,
        "unaffected_literature": 1,
    }
    row.update(changes)
    return row


def population(*rows):
    return pd.DataFrame(rows, columns=list(population_row()))


def clinical(*rows):
    return pd.DataFrame(rows, columns=list(clinical_row()))


def test_population_only_variants_are_included_without_any_clinical_rows():
    p = population(population_row(), population_row("7-110-C-G", protein_key="A20G"))
    units, ledger, membership = build_union(p, clinical())
    assert len(units) == len(membership) == 2
    assert ledger.empty
    assert units.origin.eq("population_only").all()
    assert units.affected.eq(0).all()
    assert units.gnomad_carriers.sum() == 14


def test_literature_only_counts_survive_without_population_matches():
    units, ledger, membership = build_union(population(), clinical(clinical_row()))
    assert membership.empty
    assert ledger.prior_literature_eligible.all()
    assert units.iloc[0].origin == "literature_only"
    assert units.iloc[0].affected == 4
    assert units.iloc[0].unaffected == 1


def test_one_clinical_protein_aggregate_owns_multiple_dna_alleles_once():
    p = population(
        population_row(),
        population_row(
            "7-101-G-A", hgvsc="c.30G>A", gnomad_carriers=3, joint_ac=3, joint_hom=0
        ),
        population_row("7-150-C-T", protein_key="R20X", hgvsc="c.58C>T", aa_pos=20),
    )
    units, ledger, membership = build_union(p, clinical(clinical_row()))
    aggregate = units.loc[units.literature_key.eq("A10V")].iloc[0]
    assert aggregate.unit_grain == "clinical_protein_aggregate"
    assert aggregate.n_member_alleles == 2
    assert aggregate.affected == 4
    assert aggregate.unaffected_literature == 1
    assert aggregate.gnomad_carriers == 10
    assert aggregate.unaffected == 11
    assert aggregate.n == 15
    assert units.affected.sum() == 4
    assert len(membership) == 3 and membership.variant_id.is_unique
    assert ledger.prior_literature_eligible.all()


@pytest.mark.parametrize(
    "source,expected",
    [
        ("p.R186*", "R186X"),
        ("R186X", "R186X"),
        ("G258G", "G258="),
        ("p.G258=", "G258="),
    ],
)
def test_stop_and_synonymous_point_notation_equivalence(source, expected):
    assert normalized_protein(source) == expected


@pytest.mark.parametrize(
    "pop_key,lit_key,vclass,ref,alt",
    [
        ("G258G", "G258=", "synonymous", "G", "G"),
        ("R186*", "R186X", "nonsense", "R", "*"),
    ],
)
def test_equivalent_point_notations_join_in_actual_union(
    pop_key, lit_key, vclass, ref, alt
):
    p = population(
        population_row(protein_key=pop_key, vclass=vclass, aa_ref=ref, aa_alt=alt)
    )
    lit = clinical(
        clinical_row(
            lit_key, canonical_point_key=lit_key, vclass=vclass, aa_ref=ref, aa_alt=alt
        )
    )
    units, ledger, _ = build_union(p, lit)
    assert len(units) == 1
    assert ledger.iloc[0].join_method == "canonical_protein_aggregate"
    assert units.iloc[0].affected == 4


def test_distinct_clinical_keys_claiming_same_allele_are_quarantined_not_summed():
    units, ledger, membership = build_union(
        population(population_row()),
        clinical(clinical_row(), clinical_row("alias_A10V", affected_literature=9)),
    )
    assert len(units) == len(membership) == 1
    assert units.iloc[0].origin == "population_only"
    assert units.iloc[0].affected == 0
    assert not ledger.prior_literature_eligible.any()
    assert set(ledger.join_exclusion_reason) == {
        "multiple_literature_keys_match_same_population_allele"
    }


def test_canonical_cdna_match_requires_correct_population_transcript():
    units, ledger, _ = build_union(
        population(population_row(transcript_id="ENST_OTHER.1")),
        clinical(
            clinical_row(
                "c.29C>T",
                canonical_point_key="",
                protein_point_join_allowed=False,
                cdna="c.29C>T",
            )
        ),
    )
    assert len(units) == 2
    assert ledger.iloc[0].join_method == "literature_only_unmatched"
    pop_unit = units.loc[units.origin.eq("population_only")].iloc[0]
    assert pop_unit.protein_key == ""
    assert pop_unit.canonical_wt_status == "unavailable"


@pytest.mark.parametrize(
    "change",
    [
        {"feature_transcript": "ENST_OTHER.1"},
        {"cdna": "ENST_OTHER.1:c.29C>T"},
        {"cdna": "ENST_OTHER.1%3Ac.29C%3ET"},
    ],
)
def test_incompatible_clinical_transcript_is_not_erased_before_cdna_join(change):
    units, ledger, _ = build_union(
        population(population_row()), clinical(clinical_row(**change))
    )
    assert len(units) == 1 and units.iloc[0].origin == "population_only"
    assert not ledger.iloc[0].prior_literature_eligible
    assert "transcript" in ledger.iloc[0].join_exclusion_reason


def test_matching_explicit_transcript_cdna_is_accepted_without_protein_key():
    units, ledger, _ = build_union(
        population(population_row(hgvsc=f"{TX}:c.29C>T")),
        clinical(
            clinical_row(
                "c.29C>T",
                cdna=f"{TX}%3Ac.29C%3ET",
                canonical_point_key="",
                protein_point_join_allowed=False,
            )
        ),
    )
    assert len(units) == 1
    assert ledger.iloc[0].join_method == "exact_canonical_cdna"


def test_cdna_prefix_on_population_annotation_cannot_be_silently_discarded():
    units, ledger, _ = build_union(
        population(population_row(hgvsc="ENST_OTHER.1:c.29C>T")),
        clinical(
            clinical_row(
                "c.29C>T",
                cdna="c.29C>T",
                canonical_point_key="",
                protein_point_join_allowed=False,
            )
        ),
    )
    assert len(units) == 2
    assert ledger.iloc[0].join_method == "literature_only_unmatched"


def test_multiple_genomic_matches_for_cdna_only_identity_are_quarantined():
    units, ledger, membership = build_union(
        population(population_row(), population_row("7-101-G-A", protein_key="G11D")),
        clinical(
            clinical_row(
                "c.29C>T",
                cdna="c.29C>T",
                canonical_point_key="",
                protein_point_join_allowed=False,
            )
        ),
    )
    assert len(units) == len(membership) == 2
    assert units.affected.sum() == 0
    assert (
        ledger.iloc[0].join_exclusion_reason
        == "canonical_cdna_matches_multiple_genomic_alleles"
    )


def test_conflicting_cdna_cannot_override_canonical_protein_identity():
    units, ledger, _ = build_union(
        population(population_row(protein_key="G11D")),
        clinical(clinical_row(cdna="c.29C>T")),
    )
    assert len(units) == 1 and units.iloc[0].affected == 0
    assert ledger.iloc[0].join_exclusion_reason == "canonical_cdna_and_protein_conflict"


def test_existing_identity_quarantine_does_not_consume_population_allele():
    units, ledger, membership = build_union(
        population(population_row()),
        clinical(
            clinical_row(
                identity_quarantine_recommended=True,
                identity_quarantine_reason="canonical_protein_identity_invalid",
            )
        ),
    )
    assert len(membership) == 1
    assert units.iloc[0].origin == "population_only"
    assert not ledger.iloc[0].prior_literature_eligible


def test_unknown_legacy_cdna_stays_unmatched_rather_than_becoming_a_transcript_error():
    _, ledger, _ = build_union(
        population(population_row()),
        clinical(
            clinical_row(
                "IVS3+1G>A",
                cdna="IVS3+1G>A",
                canonical_point_key="",
                protein_point_join_allowed=False,
            )
        ),
    )
    assert ledger.iloc[0].prior_literature_eligible
    assert ledger.iloc[0].join_method == "literature_only_unmatched"


def test_population_qc_exclusions_and_duplicate_alleles_do_not_add_counts():
    units, _, membership = build_union(
        population(
            population_row(), population_row("7-200-C-T", population_eligible=False)
        ),
        clinical(clinical_row()),
    )
    assert len(membership) == 1 and units.gnomad_carriers.sum() == 7
    with pytest.raises(ValueError, match="Duplicate population allele"):
        build_union(
            population(population_row(), population_row()), clinical(clinical_row())
        )


def test_population_counts_require_ac_minus_homozygotes_and_input_counts_nonnegative():
    with pytest.raises(ValueError, match="AC minus homozygotes"):
        build_union(
            population(population_row(gnomad_carriers=8)), clinical(clinical_row())
        )
    with pytest.raises(ValueError, match="clinical unaffected_literature"):
        build_union(
            population(population_row()),
            clinical(clinical_row(unaffected_literature=-1)),
        )


def test_posterior_alpha_gets_affected_beta_gets_both_unaffected_sources():
    units, _, _ = build_union(population(population_row()), clinical(clinical_row()))
    posterior = posterior_table(units, {"alpha_empirical": 2, "beta_empirical": 3})
    row = posterior.iloc[0]
    assert row.posterior_alpha == 2 + 4
    assert row.posterior_beta == 3 + 1 + 7
    assert row.posterior_mean == pytest.approx(6 / 17)


def test_empirical_prior_uses_saturating_variant_weights_not_pooled_counts():
    units = pd.DataFrame({"affected": [1, 6, 0], "n": [1, 10, 100]})
    parameters = fit_empirical(units)
    expected = np.average(units.affected / units.n, weights=1 - 1 / (units.n + 0.01))
    assert parameters["mean"] == pytest.approx(expected)
    assert parameters["mean"] != pytest.approx(units.affected.sum() / units.n.sum())
    assert parameters["alpha_empirical"] > 0 and parameters["beta_empirical"] > 0
    with pytest.raises(ValueError, match="Invalid empirical Beta moments"):
        fit_empirical(pd.DataFrame({"affected": [0, 0], "n": [10, 20]}))
