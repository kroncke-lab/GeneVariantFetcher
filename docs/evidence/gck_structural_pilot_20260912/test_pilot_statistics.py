"""Offline checks of the pilot's outer variant exclusion and comparison math."""

import numpy as np
import pandas as pd

from pilot_statistics import (
    comparison_metrics,
    density_values,
    fractional_logistic,
    remove_donor,
    variant_loo_comparison,
)


def test_global_exclusion_removes_heldout_from_every_training_row():
    weights = np.array([[0, 0.6, 0.4], [0.9, 0, 0.1], [0.3, 0.7, 0]])
    removed = remove_donor(weights, 0)
    np.testing.assert_allclose(removed, [[0, 0.6, 0.4], [0, 0, 1], [0, 1, 0]])
    # Both training features ignore the target's value, even when it dominated
    # one of the uncorrected local neighborhoods.
    before = density_values(removed, np.array([0.1, 0.7, 0.8]))
    after = density_values(removed, np.array([0.99, 0.7, 0.8]))
    np.testing.assert_array_equal(before, after)
    np.testing.assert_array_equal(weights[:, 0], [0, 0.9, 0.3])


def test_unsupported_density_stays_missing_after_exclusion():
    weights = np.array([[0, 1], [1, 0]], dtype=float)
    density = density_values(remove_donor(weights, 0), np.array([0.2, 0.8]))
    assert density[0] == 0.8
    assert np.isnan(density[1])


def test_intercept_fit_recovers_equal_variant_posterior_mean():
    y = np.array([0.1, 0.2, 0.9])
    result = fractional_logistic(np.empty((3, 0)), y, np.empty((1, 0)))
    np.testing.assert_allclose(result, [y.mean()], atol=1e-8)


def test_heldout_prediction_does_not_use_own_posterior_or_counts():
    n = 12
    positions = np.arange(n)
    weights = np.exp(-abs(positions[:, None] - positions[None, :]))
    np.fill_diagonal(weights, 0)
    weights /= weights.sum(axis=1, keepdims=True)
    variants = pd.DataFrame(
        {
            "key": [f"A{i}V" for i in range(n)],
            "posterior_empirical_mean": np.linspace(0.2, 0.8, n),
            "am": np.linspace(0.1, 0.9, n),
            "affected": np.arange(n),
            "n": 20,
            "unaffected_total": 20 - np.arange(n),
        }
    )
    before = variant_loo_comparison(variants, weights, weights, prior_mean=0.5)
    changed = variants.copy()
    changed.loc[0, ["posterior_empirical_mean", "affected", "unaffected_total"]] = [
        0.99,
        19,
        1,
    ]
    after = variant_loo_comparison(changed, weights, weights, prior_mean=0.5)
    first_before = before.loc[before.key.eq("A0V")].set_index("model").prediction
    first_after = after.loc[after.key.eq("A0V")].set_index("model").prediction
    np.testing.assert_allclose(first_before, first_after, atol=1e-12)
    before["support_stratum"] = "all"
    metrics = comparison_metrics(before, 3.8)
    assert (metrics.n_variants == n).all()
    assert np.isfinite(metrics.mean_beta_binomial_nll).all()
