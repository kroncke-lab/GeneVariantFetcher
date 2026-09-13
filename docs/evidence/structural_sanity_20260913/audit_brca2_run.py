"""Independent saved-array, global exclusion, and heldout prediction audit.

Uses fresh PPA reference exclusions for selected targets, a separate Newton
solver for selected regressions, and the log-Beta formula for count scores.
The existing full gene/missense prior remains fixed by design.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from scipy.special import betaln, expit, gammaln


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE / "analysis/BRCA2"
RAW = REPO / "results/structural_sanity_20260913/BRCA2"
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density, sigmoid_kernel

spec = importlib.util.spec_from_file_location("audited_runner", HERE / "run_brca2.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(4 * 1024 * 1024):
            value.update(block)
    return value.hexdigest()


def read_shards(stem):
    files = sorted(OUT.glob(f"{stem}.part*.csv.gz"))
    assert files, f"No {stem} shards available yet"
    return pd.concat(
        [pd.read_csv(p, low_memory=False) for p in files], ignore_index=True
    )


def selected_targets(primary, geometry):
    experiment = set(
        geometry.loc[
            geometry.geometry_state.eq("structured")
            & geometry.frame_id.isin(["7LDG", "8PBC"]),
            "canonical_pos",
        ]
    )
    categories = np.where(
        primary.canonical_pos.isin(experiment),
        "experimental",
        np.where(
            primary.density_source.eq("structured"), "alphafold", primary.density_source
        ),
    )
    rows = []
    for kind in ["experimental", "alphafold", "polymer"]:
        subset = primary.loc[(categories == kind) & primary.density.notna()]
        assert len(subset) >= 2
        for offset in [len(subset) // 4, 3 * len(subset) // 4]:
            row = subset.iloc[offset]
            rows.append({"variant_id": row.variant_id, "geometry_kind": kind})
    return rows


def newton_prediction(train_x, labels, test_x):
    """Independent Newton solution of the fixed ridge=1 fractional likelihood."""
    mean, scale = train_x.mean(axis=0), train_x.std(axis=0)
    scale[scale <= 1e-12] = 1.0
    x = np.column_stack([np.ones(len(labels)), (train_x - mean) / scale])
    test = np.r_[1.0, (test_x - mean) / scale]
    coef = np.zeros(x.shape[1])
    coef[0] = np.log(labels.mean() / (1 - labels.mean()))
    penalty = np.diag(np.r_[0.0, np.ones(x.shape[1] - 1)])

    def objective(b):
        eta = x @ b
        return np.sum(np.logaddexp(0, eta) - labels * eta) + 0.5 * b @ penalty @ b

    for _ in range(100):
        probability = expit(x @ coef)
        gradient = x.T @ (probability - labels) + penalty @ coef
        if abs(gradient).max() < 1e-9:
            break
        hessian = (x.T * (probability * (1 - probability))) @ x + penalty
        step = np.linalg.solve(hessian, gradient)
        multiplier, loss = 1.0, objective(coef)
        while objective(coef - multiplier * step) > loss + 1e-12:
            multiplier /= 2
            assert multiplier > 1e-10
        coef -= multiplier * step
    assert abs(x.T @ (expit(x @ coef) - labels) + penalty @ coef).max() < 1e-6
    return float(expit(test @ coef))


def feature_audit(variants, primary, geometry, exclusion, cache):
    ids = variants.variant_id.tolist()
    y = variants.posterior_mean.to_numpy()
    assert primary.variant_id.tolist() == ids
    assert (
        hashlib.sha256("\n".join(ids).encode()).hexdigest()
        == cache["variant_ids_sha256"]
    )
    assert sha(RAW / "all_global_exclusions.float64") == cache["exclusion_cache_sha256"]
    np.testing.assert_allclose(
        np.diag(exclusion), primary.density, atol=1e-12, equal_nan=True
    )
    unsupported = primary.density.isna().to_numpy()
    assert np.isnan(exclusion[unsupported]).all()
    max_density_error, max_variance_error, nonzero_pairs = 0.0, 0.0, 0
    expected_start = 0
    for name, expected_hash in sorted(cache["raw_weight_shards"].items()):
        assert sha(RAW / name) == expected_hash
        start, end = [int(value) for value in Path(name).stem.split("_")[1:]]
        assert start == expected_start
        expected_start = end
        weights = load_npz(RAW / name)
        assert weights.shape == (end - start, len(ids))
        assert (weights.data > 0).all() and np.isfinite(weights.data).all()
        assert not np.asarray(
            weights[np.arange(end - start), np.arange(start, end)]
        ).any()
        sums = np.asarray(weights.sum(axis=1)).ravel()
        np.testing.assert_allclose(sums, ~unsupported[start:end], atol=1e-12)
        density = weights @ y
        density[sums == 0] = np.nan
        np.testing.assert_allclose(
            density, primary.density.iloc[start:end], atol=1e-12, equal_nan=True
        )
        variance = weights.multiply(weights) @ variants.posterior_variance.to_numpy()
        variance[sums == 0] = np.nan
        np.testing.assert_allclose(
            variance,
            primary.density_conditional_variance.iloc[start:end],
            atol=1e-12,
            equal_nan=True,
        )
        max_density_error = max(
            max_density_error,
            float(np.nanmax(abs(density - primary.density.iloc[start:end]))),
        )
        max_variance_error = max(
            max_variance_error,
            float(
                np.nanmax(
                    abs(variance - primary.density_conditional_variance.iloc[start:end])
                )
            ),
        )
        nonzero_pairs += weights.nnz
    assert expected_start == len(ids)
    assert nonzero_pairs == cache["positive_primary_donor_pairs"]
    selected = selected_targets(primary, geometry)
    targets = [row["variant_id"] for row in selected]
    target_indices = pd.Index(ids).get_indexer(targets)
    # Missing/ambiguous rows cannot make any edge; omitting them only removes
    # empty contexts and makes the independent enumeration substantially faster.
    usable = geometry.loc[geometry.geometry_state.isin(["structured", "idr"])]
    checks = []
    for heldout in [None, *targets[::2]]:
        fresh = empirical_variant_density(
            variants,
            usable,
            target_ids=targets,
            include_context_weights=False,
            excluded_variant_ids=[] if heldout is None else [heldout],
            backend="reference",
        )
        stored = (
            primary.density.iloc[target_indices].to_numpy()
            if heldout is None
            else exclusion[target_indices, ids.index(heldout)]
        )
        np.testing.assert_allclose(
            fresh.summary.density, stored, atol=1e-11, equal_nan=True
        )
        checks.append(
            {
                "heldout": heldout,
                "targets": len(targets),
                "max_abs_error": float(
                    np.nanmax(abs(fresh.summary.density.to_numpy() - stored))
                ),
            }
        )
    scores = pd.read_csv(
        HERE.parent
        / "population_inclusive_penetrance_20260912/predictors/population_predictors.csv.gz"
    )
    scores = scores.loc[scores.gene.eq("BRCA2")].set_index("variant_id")
    conflict_units = []
    for row in variants.itertuples():
        members = (
            [] if pd.isna(row.member_alleles) else str(row.member_alleles).split(";")
        )
        conflicts = [
            member
            for member in members
            if member in scores.index
            and "conflict" in str(scores.loc[member, "alphamissense_status"])
        ]
        if conflicts:
            assert (
                pd.isna(row.am) and row.am_source == "member_version_or_value_conflict"
            )
            conflict_units.append(row.variant_id)
    assert not (
        variants.origin.eq("population_only")
        & variants.am_source.eq("archived_clinical_key_fallback")
    ).any()
    return {
        "variants": len(ids),
        "unsupported": int(unsupported.sum()),
        "all_cache_diagonals_match_primary": True,
        "all_unsupported_cache_rows_nan": True,
        "all_saved_weight_self_entries_zero": True,
        "checked_weight_shards": len(cache["raw_weight_shards"]),
        "positive_pairs": nonzero_pairs,
        "max_weight_reconstruction_error": max_density_error,
        "max_variance_reconstruction_error": max_variance_error,
        "fresh_reference_exclusions": checks,
        "selected_targets": selected,
        "am_conflict_units_kept_missing": len(conflict_units),
        "population_only_archived_am_fallbacks": 0,
        "am_sources": variants.am_source.value_counts().to_dict(),
        "exclusion_cache_sha256": cache["exclusion_cache_sha256"],
    }


def prediction_audit(variants, primary, exclusion, prior, selected):
    predictions = read_shards("loo_predictions")
    assert not predictions.duplicated(["variant_id", "cohort", "model"]).any()
    ids = pd.Index(variants.variant_id)
    y, am, positions = (
        variants.posterior_mean.to_numpy(),
        variants.am.to_numpy(),
        variants.canonical_pos.to_numpy(),
    )
    d0 = primary.density.to_numpy()
    sequence_sums, sequence_numerators = np.zeros(len(ids)), np.zeros(len(ids))
    for start in range(0, len(ids), 128):
        end = min(start + 128, len(ids))
        weights = sigmoid_kernel(
            3.8 * np.sqrt(abs(positions[start:end, None] - positions[None, :]))
        )
        weights[np.arange(end - start), np.arange(start, end)] = 0
        sequence_sums[start:end], sequence_numerators[start:end] = (
            weights.sum(axis=1),
            weights @ y,
        )
    sequence0 = sequence_numerators / sequence_sums
    refits = []
    for cohort, rows in predictions.groupby("cohort"):
        common = np.flatnonzero(
            np.isfinite(d0) & (np.isfinite(am) if cohort == "am_common" else True)
        )
        sets = [set(group.variant_id) for _, group in rows.groupby("model")]
        assert all(values == set(ids[common]) for values in sets)
        for item in selected:
            heldout = ids.get_loc(item["variant_id"])
            if heldout not in common:
                continue
            k = sigmoid_kernel(3.8 * np.sqrt(abs(positions - positions[heldout])))
            k[heldout] = 0
            sequence = (sequence_numerators - k * y[heldout]) / (sequence_sums - k)
            density = np.asarray(exclusion[:, heldout])
            training = common[common != heldout]
            valid = np.isfinite(density[training] + sequence[training])
            lost = int((~valid).sum())
            training = training[valid]
            saved = rows.loc[rows.variant_id.eq(item["variant_id"])].set_index("model")
            assert saved.training_targets.eq(len(training)).all()
            assert saved.training_targets_lost_support.eq(lost).all()
            calculated = {
                "empirical_prior": prior["mean"],
                "intercept_only": y[training].mean(),
                "density_unfitted": d0[heldout],
                "sequence_unfitted": sequence0[heldout],
            }
            features = {
                "density_fit": (density[:, None], np.array([d0[heldout]])),
                "sequence_fit": (sequence[:, None], np.array([sequence0[heldout]])),
            }
            if cohort == "am_common":
                features.update(
                    {
                        "am_fit": (am[:, None], np.array([am[heldout]])),
                        "am_plus_density": (
                            np.column_stack([am, density]),
                            np.array([am[heldout], d0[heldout]]),
                        ),
                        "am_plus_sequence": (
                            np.column_stack([am, sequence]),
                            np.array([am[heldout], sequence0[heldout]]),
                        ),
                    }
                )
            for name, (train_x, test_x) in features.items():
                calculated[name] = newton_prediction(
                    train_x[training], y[training], test_x
                )
            for name, expected in calculated.items():
                actual = float(saved.loc[name, "prediction"])
                assert abs(expected - actual) < 5e-5, (item, name, expected, actual)
                refits.append(
                    {
                        **item,
                        "cohort": cohort,
                        "model": name,
                        "saved_prediction": actual,
                        "independent_prediction": expected,
                        "absolute_difference": abs(expected - actual),
                    }
                )
    recomputed = []
    metrics = pd.read_csv(OUT / "loo_metrics.csv")
    for _, metric in metrics.iterrows():
        rows = predictions.loc[
            predictions.cohort.eq(metric.cohort) & predictions.model.eq(metric.model)
        ]
        if metric.support_stratum != "all":
            rows = rows.loc[rows.density_source.eq(metric.support_stratum)]
        p, affected, n = (
            rows.prediction.to_numpy(),
            rows.affected.to_numpy(),
            rows.n.to_numpy(),
        )
        np.testing.assert_allclose(affected + rows.unaffected, n)
        np.testing.assert_allclose(rows.observed_fraction, affected / n)
        expected_alpha, expected_beta = (
            p * prior["strength"],
            (1 - p) * prior["strength"],
        )
        logpmf = gammaln(n + 1) - gammaln(affected + 1) - gammaln(n - affected + 1)
        logpmf += betaln(
            affected + expected_alpha, n - affected + expected_beta
        ) - betaln(expected_alpha, expected_beta)
        values = {
            "mae_empirical_posterior": np.mean(abs(p - rows.empirical_posterior)),
            "mse_empirical_posterior": np.mean((p - rows.empirical_posterior) ** 2),
            "mae_observed_fraction": np.mean(abs(p - rows.observed_fraction)),
            "mse_observed_fraction": np.mean((p - rows.observed_fraction) ** 2),
            "mean_beta_binomial_nll": -logpmf.mean(),
        }
        assert len(rows) == metric.n_variants
        for name, expected in values.items():
            np.testing.assert_allclose(expected, metric[name], rtol=1e-8, atol=1e-9)
        recomputed.append(
            {
                "cohort": metric.cohort,
                "model": metric.model,
                "support_stratum": metric.support_stratum,
                **values,
            }
        )
    pd.DataFrame(refits).to_csv(OUT / "audit_selected_predictions.csv", index=False)
    pd.DataFrame(recomputed).to_csv(OUT / "audit_recomputed_metrics.csv", index=False)
    return {
        "prediction_rows": len(predictions),
        "cohorts": predictions.cohort.value_counts().to_dict(),
        "same_target_sets_across_models": True,
        "independent_refit_rows": len(refits),
        "maximum_independent_prediction_error": max(
            row["absolute_difference"] for row in refits
        ),
        "count_score_metric_rows_recomputed": len(recomputed),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features-only", action="store_true")
    parser.add_argument("--predictions-only", action="store_true")
    args = parser.parse_args()
    assert not (args.features_only and args.predictions_only)
    variants, prior, _ = runner.frozen.load_variants("BRCA2")
    primary = read_shards("primary_density")
    assert primary.variant_id.tolist() == variants.variant_id.tolist()
    geometry = pd.read_csv(
        HERE / "geometry/BRCA2/primary_geometry.csv.gz", low_memory=False
    ).fillna({"idr_segment": ""})
    cache = json.loads((OUT / "cache_manifest.json").read_text())
    exclusion = np.memmap(
        RAW / "all_global_exclusions.float64",
        dtype="float64",
        mode="r",
        shape=(len(variants), len(variants)),
        order="F",
    )
    receipt_path = OUT / "independent_audit.json"
    if args.predictions_only:
        receipt = json.loads(receipt_path.read_text())
        assert receipt["features"]["exclusion_cache_sha256"] == sha(
            RAW / "all_global_exclusions.float64"
        )
    else:
        receipt = {
            "features": feature_audit(variants, primary, geometry, exclusion, cache)
        }
    if not args.features_only:
        receipt["predictions"] = prediction_audit(
            variants, primary, exclusion, prior, receipt["features"]["selected_targets"]
        )
    receipt["full_class_prior_fixed_by_design"] = True
    receipt["script_sha256"], receipt["runner_sha256"] = (
        sha(__file__),
        sha(HERE / "run_brca2.py"),
    )
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
