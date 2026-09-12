"""Offline validation of fast global donor exclusion and fixed-target LOO."""

import numpy as np
import pandas as pd
from types import SimpleNamespace

from run_structure import resolve_am

from structure_statistics import (
    density_values,
    exclude_donor_values,
    variant_loo_comparison,
)


def test_fast_global_exclusion_matches_explicit_matrix_and_ignores_heldout_value():
    rng = np.random.default_rng(2)
    weights = rng.uniform(size=(18, 18))
    np.fill_diagonal(weights, 0)
    weights /= weights.sum(axis=1, keepdims=True)
    values = rng.uniform(size=18)
    for heldout in [0, 8, 17]:
        explicit = weights.copy()
        explicit[:, heldout] = 0
        explicit /= explicit.sum(axis=1, keepdims=True)
        expected = explicit @ values
        actual = exclude_donor_values(
            weights, values, heldout, density_values(weights, values)
        )
        np.testing.assert_allclose(actual, expected, atol=2e-15)
        changed = values.copy()
        changed[heldout] = 0.99
        np.testing.assert_allclose(
            exclude_donor_values(weights, changed, heldout), expected, atol=2e-15
        )


def test_near_dominant_and_unsupported_rows_keep_correct_remaining_support():
    weights = np.array([[0, 1, 0], [1 - 1e-15, 0, 1e-15], [1, 0, 0]])
    result = exclude_donor_values(weights, np.array([0.8, 0.2, 0.3]), 0)
    np.testing.assert_allclose(result[:2], [0.2, 0.3])
    assert np.isnan(result[2])


def test_outer_loo_prediction_ignores_own_outcome_and_preserves_same_residue_variants():
    n = 14
    positions = np.arange(n) // 2
    weights = np.exp(-abs(positions[:, None] - positions[None, :]))
    np.fill_diagonal(weights, 0)
    weights /= weights.sum(axis=1, keepdims=True)
    assert weights[0, 1] > 0
    variants = pd.DataFrame(
        {
            "variant_id": [f"v{i}" for i in range(n)],
            "unit_id": [f"v{i}" for i in range(n)],
            "literature_key": [f"A{i}V" for i in range(n)],
            "protein_key": [f"A{i}V" for i in range(n)],
            "origin": "population_only",
            "posterior_mean": np.linspace(0.02, 0.8, n),
            "am": np.linspace(0.1, 0.9, n),
            "affected": np.arange(n),
            "n": 20,
            "unaffected": 20 - np.arange(n),
        }
    )
    targets = variants.variant_id.iloc[:12].tolist()
    before = variant_loo_comparison(
        variants, weights, weights, prior_mean=0.2, target_ids=targets
    )
    changed = variants.copy()
    changed.loc[0, ["posterior_mean", "affected", "unaffected"]] = [0.99, 19, 1]
    after = variant_loo_comparison(
        changed, weights, weights, prior_mean=0.2, target_ids=targets
    )
    np.testing.assert_allclose(
        before.query("variant_id=='v0'").prediction,
        after.query("variant_id=='v0'").prediction,
        atol=1e-12,
    )
    assert set(before.variant_id) == set(targets)
    assert before.training_target_count.eq(11).all()
    assert before.eligible_donor_pool.eq(14).all()


def test_am_conflicts_do_not_trigger_archived_fallback_and_genomic_scores_take_priority():
    clinical = SimpleNamespace(member_alleles="allele1", literature_key="A1V")
    scores = {
        "allele1": {
            "alphamissense": 0.2,
            "alphamissense_version": "v1",
            "alphamissense_status": "available",
        }
    }
    assert resolve_am(clinical, scores, {"A1V": 0.9}) == (0.2, "exact_genomic_members")
    scores["allele1"]["alphamissense_status"] = "version_or_value_conflict"
    scores["allele1"]["alphamissense"] = np.nan
    value, source = resolve_am(clinical, scores, {"A1V": 0.9})
    assert np.isnan(value) and source == "member_version_or_value_conflict"
    assert resolve_am(clinical, {}, {"A1V": 0.9}) == (
        0.9,
        "archived_clinical_key_fallback",
    )


def test_distinct_member_allele_scores_are_not_arbitrarily_collapsed():
    clinical = SimpleNamespace(member_alleles="allele1;allele2", literature_key="A1V")
    scores = {
        "allele1": {
            "alphamissense": 0.2,
            "alphamissense_version": "v1",
            "alphamissense_status": "available",
        },
        "allele2": {
            "alphamissense": 0.8,
            "alphamissense_version": "v1",
            "alphamissense_status": "available",
        },
    }
    value, source = resolve_am(clinical, scores, {"A1V": 0.9})
    assert np.isnan(value) and source == "different_member_scores_or_versions"
