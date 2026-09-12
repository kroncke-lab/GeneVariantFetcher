"""Variant-only outer LOO for the population-inclusive monomer experiment."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


OLD = Path(__file__).resolve().parent.parent / "gck_structural_pilot_20260912"
spec = importlib.util.spec_from_file_location(
    "frozen_gck_statistics", OLD / "pilot_statistics.py"
)
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
fractional_logistic = legacy.fractional_logistic
comparison_metrics = legacy.comparison_metrics
density_values = legacy.density_values


def exclude_donor_values(weights, values, heldout, baseline=None):
    """Globally remove one donor without copying/re-multiplying a square matrix.

    This is a monomer-only shortcut. The row's own variant is already excluded.
    Near-total donor dominance uses direct remaining-weight sums to avoid
    subtractive cancellation. Truly unsupported rows remain missing.
    """
    sums = weights.sum(axis=1)
    numerator = weights @ values if baseline is None else np.nan_to_num(baseline) * sums
    heldout_weight = weights[:, heldout]
    denominator = sums - heldout_weight
    remainder = numerator - heldout_weight * values[heldout]
    result = np.full(len(values), np.nan)
    stable = denominator > np.maximum(1e-8 * sums, 0)
    np.divide(remainder, denominator, out=result, where=stable)
    for row in np.flatnonzero(~stable):
        remaining = weights[row].copy()
        remaining[heldout] = 0
        total = remaining.sum()
        if total > 0:
            result[row] = remaining @ values / total
    return result


def variant_loo_comparison(
    variants,
    weights,
    sequence_weights,
    *,
    prior_mean,
    target_ids=None,
    cohort="all_common",
):
    """Fixed settings, equal-variant regression labels, globally excluded donors.

    Optional target_ids restrict both fitting targets and evaluation targets;
    the expanded eligible donor pool remains available in every feature. Common
    empirical hyperparameters are fixed from the full locus dataset by design.
    """
    n = len(variants)
    if weights.shape != (n, n) or sequence_weights.shape != (n, n):
        raise ValueError("Expected two square matrices aligned to the variant table")
    if np.any(np.diag(weights) != 0) or np.any(np.diag(sequence_weights) != 0):
        raise ValueError("Each row must exclude its own variant identity")
    y = variants.posterior_mean.to_numpy(dtype=float)
    am = variants.am.to_numpy(dtype=float)
    density = density_values(weights, y)
    sequence = density_values(sequence_weights, y)
    mask = np.isfinite(am + density + sequence)
    if target_ids is not None:
        unknown = set(target_ids) - set(variants.variant_id)
        if unknown:
            raise ValueError(
                f"Unknown requested comparison identities: {sorted(unknown)}"
            )
        mask &= variants.variant_id.isin(target_ids).to_numpy()
    common = np.flatnonzero(mask)
    if len(common) < 10:
        raise ValueError("Too few common targets for the prespecified comparison")
    records = []
    for heldout in common:
        training = common[common != heldout]
        d = exclude_donor_values(weights, y, heldout, density)
        s = exclude_donor_values(sequence_weights, y, heldout, sequence)
        if not np.isfinite(d[training] + s[training]).all():
            raise ValueError("Outer exclusion removed all support from a training row")
        features = {
            "intercept_only": (np.empty((len(training), 0)), np.empty((1, 0))),
            "density_fit": (d[training, None], density[[heldout], None]),
            "am_fit": (am[training, None], am[[heldout], None]),
            "am_plus_density": (
                np.column_stack([am[training], d[training]]),
                np.array([[am[heldout], density[heldout]]]),
            ),
            "sequence_fit": (s[training, None], sequence[[heldout], None]),
            "am_plus_sequence": (
                np.column_stack([am[training], s[training]]),
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
        row = variants.iloc[heldout]
        for model, value in predictions.items():
            records.append(
                {
                    "key": row.unit_id,
                    "variant_id": row.variant_id,
                    "literature_key": row.literature_key,
                    "protein_key": row.protein_key,
                    "origin": row.origin,
                    "cohort": cohort,
                    "model": model,
                    "prediction": value,
                    "empirical_posterior": y[heldout],
                    "observed_fraction": row.affected / row.n,
                    "affected": row.affected,
                    "unaffected": row.unaffected,
                    "n": row.n,
                    "density": density[heldout],
                    "am": am[heldout],
                    "training_target_count": len(training),
                    "eligible_donor_pool": n,
                }
            )
    return pd.DataFrame(records)
