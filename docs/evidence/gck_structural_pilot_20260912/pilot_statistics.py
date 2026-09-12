"""Fixed-setting, variant-excluded diagnostics for the GCK structural pilot.

These are internal comparisons in a literature-selected count dataset. Shared
empirical hyperparameters are deliberately held fixed, as requested; this is
not independent validation of a population penetrance model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import betabinom, spearmanr


def remove_donor(weights: np.ndarray, donor_index: int) -> np.ndarray:
    """Remove a held-out donor globally, then renormalize each supported row.

    Rows already exclude their own variant. This function is for the monomer
    pilot; multicontext assemblies must re-normalize within each context first.
    """
    result = np.array(weights, dtype=float, copy=True)
    result[:, donor_index] = 0.0
    sums = result.sum(axis=1)
    np.divide(result, sums[:, None], out=result, where=sums[:, None] > 0)
    return result


def density_values(weights: np.ndarray, values: np.ndarray) -> np.ndarray:
    result = weights @ values
    result[weights.sum(axis=1) == 0] = np.nan
    return result


def fractional_logistic(
    train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, ridge: float = 1.0
) -> np.ndarray:
    """Equal-variant fractional logistic fit; fixed ridge, train-only scaling.

    Labels are empirical posterior means, not binomial carrier counts. Penalize
    slopes only. This small pre-specified model has no bandwidth/ridge tuning.
    """
    mean = train_x.mean(axis=0)
    scale = train_x.std(axis=0)
    scale = np.where(scale > 1e-12, scale, 1.0)
    x = np.column_stack([np.ones(len(train_x)), (train_x - mean) / scale])
    xt = np.column_stack([np.ones(len(test_x)), (test_x - mean) / scale])
    initial = np.zeros(x.shape[1])
    y_mean = float(train_y.mean())
    initial[0] = np.log(y_mean / (1 - y_mean))

    def objective(coef):
        eta = x @ coef
        loss = np.sum(np.logaddexp(0, eta) - train_y * eta)
        loss += 0.5 * ridge * np.sum(coef[1:] ** 2)
        gradient = x.T @ (expit(eta) - train_y)
        gradient[1:] += ridge * coef[1:]
        return loss, gradient

    fit = minimize(objective, initial, jac=True, method="L-BFGS-B")
    if not fit.success:
        raise RuntimeError(f"Fractional logistic fit failed: {fit.message}")
    return expit(xt @ fit.x)


def variant_loo_comparison(
    variants: pd.DataFrame,
    weights: np.ndarray,
    sequence_weights: np.ndarray,
    *,
    prior_mean: float,
) -> pd.DataFrame:
    """Same-target outer LOO, removing target i from ALL training densities.

    Both square matrices must have rows and donor columns in variants order.
    The full mapped donor pool remains usable even if an AM value is missing.
    Evaluation/training target rows use the intersection of available features.
    """
    if weights.shape != (len(variants), len(variants)):
        raise ValueError("Expected a square, variant-aligned monomer weight matrix")
    if sequence_weights.shape != weights.shape:
        raise ValueError("Sequence matrix shape differs")
    if not np.all(np.diag(weights) == 0) or not np.all(np.diag(sequence_weights) == 0):
        raise ValueError("LOO matrices must exclude each target's own identity")
    y = variants.posterior_empirical_mean.to_numpy()
    am = variants.am.to_numpy(dtype=float)
    density = density_values(weights, y)
    sequence = density_values(sequence_weights, y)
    common = np.flatnonzero(np.isfinite(am + density + sequence))
    if len(common) < 10:
        raise ValueError("Too few common variants for the pilot comparison")
    rows = []
    for heldout in common:
        training = common[common != heldout]
        d_without_i = density_values(remove_donor(weights, heldout), y)
        s_without_i = density_values(remove_donor(sequence_weights, heldout), y)
        if not np.isfinite(d_without_i[training] + s_without_i[training]).all():
            raise ValueError(
                "A training density lost all support after outer exclusion"
            )
        # The held-out row is already variant-excluded in the original matrices.
        features = {
            "intercept_only": (np.empty((len(training), 0)), np.empty((1, 0))),
            "density_fit": (d_without_i[training, None], density[[heldout], None]),
            "am_fit": (am[training, None], am[[heldout], None]),
            "am_plus_density": (
                np.column_stack([am[training], d_without_i[training]]),
                np.array([[am[heldout], density[heldout]]]),
            ),
            "sequence_fit": (s_without_i[training, None], sequence[[heldout], None]),
            "am_plus_sequence": (
                np.column_stack([am[training], s_without_i[training]]),
                np.array([[am[heldout], sequence[heldout]]]),
            ),
        }
        predictions = {
            "empirical_prior": prior_mean,
            "density_unfitted": density[heldout],
            "sequence_unfitted": sequence[heldout],
        }
        for model, (train_x, test_x) in features.items():
            predictions[model] = float(
                fractional_logistic(train_x, y[training], test_x)[0]
            )
        for model, prediction in predictions.items():
            row = variants.iloc[heldout]
            rows.append(
                {
                    "key": row.key,
                    "model": model,
                    "prediction": prediction,
                    "empirical_posterior": y[heldout],
                    "observed_fraction": row.affected / row.n,
                    "affected": row.affected,
                    "unaffected": row.unaffected_total,
                    "n": row.n,
                    "density": density[heldout],
                    "am": am[heldout],
                }
            )
    return pd.DataFrame(rows)


def comparison_metrics(predictions: pd.DataFrame, strength: float) -> pd.DataFrame:
    """Variant-equal errors and Beta-binomial NLL at the fixed prior strength.

    The count score describes this ascertained count experiment under the
    adopted unaffected assumption; it is not a population-outcome likelihood.
    """
    rows = []
    for (model, stratum), frame in (
        (key, frame) for key, frame in predictions.groupby(["model", "support_stratum"])
    ):
        p = frame.prediction.to_numpy()
        if not ((p > 0) & (p < 1)).all():
            raise ValueError(
                "Predictive means must be interior for Beta-binomial score"
            )
        observed = frame.observed_fraction.to_numpy()
        empirical = frame.empirical_posterior.to_numpy()
        alpha, beta = p * strength, (1 - p) * strength
        nll = -betabinom.logpmf(frame.affected, frame.n, alpha, beta)
        if not np.isfinite(nll).all():
            raise ValueError("Nonfinite count score; verify integer counts")
        rho = (
            float(spearmanr(p, observed).statistic)
            if np.ptp(p) > 1e-12 and np.ptp(observed) > 0
            else np.nan
        )
        rows.append(
            {
                "model": model,
                "support_stratum": stratum,
                "n_variants": len(frame),
                "mae_empirical_posterior": np.mean(abs(p - empirical)),
                "mse_empirical_posterior": np.mean((p - empirical) ** 2),
                "mae_observed_fraction": np.mean(abs(p - observed)),
                "mse_observed_fraction": np.mean((p - observed) ** 2),
                "mean_beta_binomial_nll": np.mean(nll),
                "spearman_observed_fraction": rho,
                "mean_prediction": p.mean(),
                "mean_observed_fraction": observed.mean(),
            }
        )
    return pd.DataFrame(rows)
