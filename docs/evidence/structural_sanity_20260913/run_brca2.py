"""Rebuild BRCA2 neighborhoods with validated local 3D and canonical polymer.

Full weight arrays and the outer-LOO cache stay in ignored results/. Compact
summaries, complete predictions, sample donor contexts and provenance are saved.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PREVIOUS = HERE.parent / "missense_structural_extension_20260912"
spec = importlib.util.spec_from_file_location(
    "frozen_extension", PREVIOUS / "run_structure.py"
)
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)
save = frozen.save
shards = frozen.shards
GEOMETRY = HERE / "geometry/BRCA2"
OUT = HERE / "analysis/BRCA2"
RAW = REPO / "results/structural_sanity_20260913/BRCA2"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def loo_predictions(
    variants, primary, exclusion, sequence_sums, sequence_density, prior, kernel
):
    y = variants.posterior_mean.to_numpy()
    am = variants.am.to_numpy()
    positions = variants.canonical_pos.to_numpy()
    d0 = primary.density.to_numpy()
    cohorts = {"all_supported": np.flatnonzero(np.isfinite(d0))}
    common = np.flatnonzero(np.isfinite(d0 + am))
    if len(common) >= 10:
        cohorts["am_common"] = common
    records = []
    for cohort, common in cohorts.items():
        for number, heldout in enumerate(common):
            training = common[common != heldout]
            density = np.asarray(exclusion[:, heldout])
            w = (
                kernel(3.8 * np.sqrt(abs(positions - positions[heldout])))
                / sequence_sums
            )
            w[heldout] = 0
            assert (1 - w).min() > 1e-8
            sequence = (sequence_density - w * y[heldout]) / (1 - w)
            valid = np.isfinite(density[training] + sequence[training])
            lost = int((~valid).sum())
            training = training[valid]
            features = {
                "density_fit": (density[training, None], d0[[heldout], None]),
                "sequence_fit": (
                    sequence[training, None],
                    sequence_density[[heldout], None],
                ),
            }
            if cohort == "am_common":
                features.update(
                    {
                        "am_fit": (am[training, None], am[[heldout], None]),
                        "am_plus_density": (
                            np.column_stack([am[training], density[training]]),
                            np.array([[am[heldout], d0[heldout]]]),
                        ),
                        "am_plus_sequence": (
                            np.column_stack([am[training], sequence[training]]),
                            np.array([[am[heldout], sequence_density[heldout]]]),
                        ),
                    }
                )
            predictions = {
                "empirical_prior": prior["mean"],
                "intercept_only": y[training].mean(),
                "density_unfitted": d0[heldout],
                "sequence_unfitted": sequence_density[heldout],
            }
            for model, (tx, vx) in features.items():
                predictions[model] = float(
                    frozen.statistics.fractional_logistic(tx, y[training], vx)[0]
                )
            row = variants.iloc[heldout]
            for model, prediction in predictions.items():
                records.append(
                    {
                        "variant_id": row.variant_id,
                        "protein_key": row.protein_key,
                        "origin": row.origin,
                        "cohort": cohort,
                        "model": model,
                        "prediction": prediction,
                        "empirical_posterior": row.posterior_mean,
                        "observed_fraction": row.affected / row.n,
                        "affected": row.affected,
                        "unaffected": row.unaffected,
                        "n": row.n,
                        "density_source": primary.iloc[heldout].density_source,
                        "training_targets": len(training),
                        "training_targets_lost_support": lost,
                    }
                )
            if number % 250 == 0:
                print(
                    f"{cohort}: exact outer LOO {number + 1}/{len(common)}", flush=True
                )
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--draws", type=int, default=8192)
    parser.add_argument("--loo-only", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
    from alphafold_rin.empirical_density import (
        empirical_variant_density,
        sigmoid_kernel,
    )
    from outer_loo import all_excluded_density

    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    variants, prior, sources = frozen.load_variants("BRCA2")
    ids = variants.variant_id.tolist()
    n = len(ids)
    y = variants.posterior_mean.to_numpy()
    positions = variants.canonical_pos.to_numpy()
    source_hashes = {str(p.relative_to(HERE.parent)): digest(p) for p in sources}
    geometry_files = [
        GEOMETRY / f"{name}_geometry.csv.gz"
        for name in ["primary", "experimental_idr", "all_unresolved_polymer"]
    ]
    source_hashes.update(
        {str(p.relative_to(HERE.parent)): digest(p) for p in geometry_files}
    )
    sequence_sums = np.zeros(n)
    sequence_density = np.zeros(n)
    for start in range(0, n, args.batch_size):
        end = min(start + args.batch_size, n)
        w = sigmoid_kernel(
            3.8 * np.sqrt(abs(positions[start:end, None] - positions[None, :]))
        )
        w[np.arange(end - start), np.arange(start, end)] = 0
        sequence_sums[start:end] = w.sum(axis=1)
        sequence_density[start:end] = w @ y / sequence_sums[start:end]
    exclusion_path = RAW / "all_global_exclusions.float64"
    if args.loo_only:
        exclusion = np.memmap(
            exclusion_path, dtype="float64", mode="r", shape=(n, n), order="F"
        )
        cache = json.loads((OUT / "cache_manifest.json").read_text())
        assert (
            cache["variant_ids_sha256"]
            == hashlib.sha256("\n".join(ids).encode()).hexdigest()
        )
        assert cache["source_hashes"] == source_hashes
        assert cache["exclusion_cache_sha256"] == digest(exclusion_path)
        primary = pd.concat(
            [pd.read_csv(p) for p in sorted(OUT.glob("primary_density.part*.csv.gz"))],
            ignore_index=True,
        )
    else:
        exclusion = np.memmap(
            exclusion_path, dtype="float64", mode="w+", shape=(n, n), order="F"
        )
        exclusion[:] = np.nan
        rng = np.random.default_rng(20260913)
        draws = rng.beta(
            variants.posterior_alpha.to_numpy()[:, None],
            variants.posterior_beta.to_numpy()[:, None],
            size=(n, args.draws),
        )
        prior_retention = prior["strength"] / (
            prior["strength"] + variants.n.to_numpy()
        )
        fraction = variants.affected.to_numpy() / variants.n.to_numpy()
        donor_variance = variants.posterior_variance.to_numpy()
        identity_cols = [
            "variant_id",
            "protein_key",
            "origin",
            "canonical_pos",
            "affected",
            "unaffected",
            "n",
            "posterior_mean",
            "am",
        ]
        tables, raw_hashes = [], {}
        scenarios = [
            ("primary", "com", 3),
            ("primary", "com", 2),
            ("primary", "com", 5),
            ("primary", "ca", 3),
            ("experimental_idr", "com", 3),
            ("all_unresolved_polymer", "com", 3),
        ]
        for policy, metric, half in scenarios:
            geometry = pd.read_csv(
                GEOMETRY / f"{policy}_geometry.csv.gz", low_memory=False
            ).fillna({"idr_segment": ""})
            if metric == "ca":
                geometry["geometry_state"] = geometry.ca_geometry_state
            scenario = f"{policy}_{metric}_h{half}"
            is_primary = scenario == "primary_com_h3"
            batches = []
            positive_pairs = 0
            for start in range(0, n, args.batch_size):
                end = min(start + args.batch_size, n)
                batch_ids = ids[start:end]
                result = empirical_variant_density(
                    variants,
                    geometry,
                    target_ids=batch_ids,
                    metric=metric,
                    half_distance=half,
                    include_context_weights=False,
                    include_context_model=is_primary,
                )
                summary = result.summary.drop(
                    columns=["canonical_pos"], errors="ignore"
                )
                weights = result.donor_weights.reindex(
                    index=batch_ids, columns=ids
                ).to_numpy()
                assert not weights[np.arange(end - start), np.arange(start, end)].any()
                predicted = weights @ y
                predicted[weights.sum(axis=1) == 0] = np.nan
                np.testing.assert_allclose(
                    predicted, summary.density, atol=1e-12, equal_nan=True
                )
                if is_primary:
                    sparse = csr_matrix(weights)
                    positive_pairs += sparse.nnz
                    weight_path = RAW / f"weights_{start:05d}_{end:05d}.npz"
                    save_npz(weight_path, sparse)
                    raw_hashes[weight_path.name] = digest(weight_path)
                    summary["neighborhood_prior_retention"] = weights @ prior_retention
                    summary["prior_component"] = (
                        summary.neighborhood_prior_retention * prior["mean"]
                    )
                    summary["counts_component"] = weights @ (
                        (1 - prior_retention) * fraction
                    )
                    summary["raw_count_neighborhood"] = weights @ fraction
                    summary["density_minus_prior"] = summary.density - prior["mean"]
                    summary["density_conditional_variance"] = (
                        sparse.multiply(sparse) @ donor_variance
                    )
                    simulated = sparse @ draws
                    lo, hi = np.quantile(simulated, [0.025, 0.975], axis=1)
                    summary["density_lower_95"], summary["density_upper_95"] = lo, hi
                    summary.loc[
                        summary.density.isna(),
                        [
                            "prior_component",
                            "counts_component",
                            "neighborhood_prior_retention",
                            "raw_count_neighborhood",
                            "density_conditional_variance",
                            "density_lower_95",
                            "density_upper_95",
                        ],
                    ] = np.nan
                    np.testing.assert_allclose(
                        summary.prior_component + summary.counts_component,
                        summary.density,
                        atol=1e-12,
                        equal_nan=True,
                    )
                    matrix = (
                        all_excluded_density(result.context_model)
                        .reindex(index=batch_ids, columns=ids)
                        .to_numpy()
                    )
                    exclusion[start:end, :] = matrix
                    for heldout in [start, min(end - 1, n - 1)]:
                        exact = result.context_model.density(
                            excluded_variant_ids=[ids[heldout]]
                        ).reindex(batch_ids)
                        np.testing.assert_allclose(
                            matrix[:, heldout], exact, atol=1e-11, equal_nan=True
                        )
                batches.append(summary)
                if start % (args.batch_size * 10) == 0:
                    print(f"{scenario}: {end}/{n} targets", flush=True)
            complete = variants[identity_cols].merge(
                pd.concat(batches, ignore_index=True),
                on="variant_id",
                validate="one_to_one",
            )
            complete["scenario"] = scenario
            tables.append(complete)
            if is_primary:
                primary = complete
                shards(primary, OUT, "primary_density", rows=1500)
                exclusion.flush()
                cache = {
                    "variant_ids_sha256": hashlib.sha256(
                        "\n".join(ids).encode()
                    ).hexdigest(),
                    "source_hashes": source_hashes,
                    "exclusion_cache_sha256": digest(exclusion_path),
                    "raw_weight_shards": raw_hashes,
                    "positive_primary_donor_pairs": positive_pairs,
                    "draws": args.draws,
                    "seed": 20260913,
                }
                (OUT / "cache_manifest.json").write_text(
                    json.dumps(cache, indent=2) + "\n"
                )
            print(
                f"{scenario}: {complete.density.notna().sum()}/{n} supported",
                flush=True,
            )
        scenarios = pd.concat(tables, ignore_index=True)
        shards(scenarios, OUT, "density_scenarios", rows=2500)
        summaries = []
        for scenario, group in scenarios.groupby("scenario"):
            paired = (
                primary[["variant_id", "density"]]
                .merge(
                    group[["variant_id", "density"]],
                    on="variant_id",
                    suffixes=("_primary", "_alternative"),
                )
                .dropna()
            )
            delta = abs(paired.density_primary - paired.density_alternative)
            summaries.append(
                {
                    "scenario": scenario,
                    "supported": int(group.density.notna().sum()),
                    "polymer": int(group.density_source.eq("polymer").sum()),
                    "structured": int(group.density_source.eq("structured").sum()),
                    "median_density": group.density.median(),
                    "shared": len(paired),
                    "mean_absolute_change": delta.mean(),
                    "max_absolute_change": delta.max(),
                }
            )
        save(pd.DataFrame(summaries), OUT / "sensitivity.csv")
        old = pd.concat(
            [
                pd.read_csv(p)
                for p in sorted(
                    (PREVIOUS / "analysis/BRCA2").glob("primary_density.part*.csv.gz")
                )
            ]
        )
        paired = old[["variant_id", "density"]].merge(
            primary[["variant_id", "density"]],
            on="variant_id",
            suffixes=("_old", "_new"),
        )
        np.testing.assert_allclose(
            paired.loc[paired.density_old.notna(), "density_old"],
            paired.loc[paired.density_old.notna(), "density_new"],
            atol=1e-12,
        )
        save(paired, OUT / "old_vs_new.csv.gz")
        selected = []
        for source, group in primary.loc[primary.density.notna()].groupby(
            "density_source"
        ):
            selected.extend(group.variant_id.iloc[[0, len(group) // 2, -1]].tolist())
        geometry = pd.read_csv(
            GEOMETRY / "primary_geometry.csv.gz", low_memory=False
        ).fillna({"idr_segment": ""})
        sample = empirical_variant_density(
            variants, geometry, target_ids=selected, include_context_weights=True
        )
        shards(sample.context_weights, OUT, "sample_donor_contexts", rows=3000)
        del draws
    predictions = loo_predictions(
        variants,
        primary,
        exclusion,
        sequence_sums,
        sequence_density,
        prior,
        sigmoid_kernel,
    )
    shards(predictions, OUT, "loo_predictions", rows=3000)
    metrics = []
    for cohort, group in predictions.groupby("cohort"):
        for source, subset in [("all", group), *group.groupby("density_source")]:
            table = frozen.statistics.comparison_metrics(
                subset.assign(support_stratum=source), prior["strength"]
            )
            table.loc[
                table.model.eq("intercept_only"), "spearman_observed_fraction"
            ] = np.nan
            table["cohort"] = cohort
            metrics.append(table)
    save(pd.concat(metrics, ignore_index=True), OUT / "loo_metrics.csv")
    checks = {
        "prior": prior,
        "variants": n,
        "supported": int(primary.density.notna().sum()),
        "density_sources": primary.density_source.value_counts().to_dict(),
        "alpha_affected_beta_all_unaffected": True,
        "same_residue_variants_retained": True,
        "exact_global_variant_loo": True,
        "full_class_prior_fixed": True,
        "empirical_neighborhood_score_is_not_variant_disease_probability": True,
        "source_hashes": source_hashes,
        "runner_sha256": digest(__file__),
        "outer_loo_sha256": digest(HERE / "outer_loo.py"),
    }
    for relative, expected in source_hashes.items():
        assert digest(HERE.parent / relative) == expected
    (OUT / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print("BRCA2 corrected run complete", flush=True)


if __name__ == "__main__":
    main()
