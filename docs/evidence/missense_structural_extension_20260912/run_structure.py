"""Extend the frozen class-specific missense experiment to verified assemblies.

Inputs are immutable count/posterior snapshots. Coordinate frames remain separate;
this explicitly partial-structure experiment is not the strict full-unit pipeline.
"""

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
REPO = HERE.parents[2]
CLASS = EVIDENCE / "class_matched_penetrance_20260912"
POP = EVIDENCE / "population_inclusive_penetrance_20260912"
PILOT = EVIDENCE / "gck_structural_pilot_20260912"
spec = importlib.util.spec_from_file_location(
    "frozen_statistics", PILOT / "pilot_statistics.py"
)
statistics = importlib.util.module_from_spec(spec)
spec.loader.exec_module(statistics)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = frame.to_csv(index=False, lineterminator="\n").encode()
    if path.suffix == ".gz":
        blob = gzip.compress(blob, mtime=0)
    if len(blob) > 1_150_000:
        raise ValueError(f"Shard too large: {path} ({len(blob)})")
    path.write_bytes(blob)


def shards(frame, folder, stem, rows=4000):
    for number, start in enumerate(range(0, len(frame), rows), 1):
        save(
            frame.iloc[start : start + rows], folder / f"{stem}.part{number:03d}.csv.gz"
        )


def load_variants(gene):
    paths = sorted(
        (CLASS / "analysis/empirical_posteriors").glob(f"{gene}.part*.csv.gz")
    )
    data = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    data = data.loc[data.variant_type.eq("missense")].reset_index(drop=True)
    prior_path = CLASS / "analysis/empirical_prior_comparison.csv"
    priors = pd.read_csv(prior_path)
    prior = priors.loc[priors.gene.eq(gene) & priors.scope.eq("canonical_missense")]
    assert len(prior) == 1
    prior = prior.iloc[0].to_dict()
    assert len(data) == prior["variants"] and not data.unit_id.duplicated().any()
    assert data.canonical_wt_status.eq("match").all()
    assert (data.n > 0).all() and not ((data.aa_pos == 1) & data.aa_ref.eq("M")).any()
    np.testing.assert_allclose(
        data.unaffected, data.unaffected_literature + data.gnomad_carriers
    )
    np.testing.assert_allclose(
        data.posterior_alpha, prior["alpha_empirical"] + data.affected
    )
    np.testing.assert_allclose(
        data.posterior_beta, prior["beta_empirical"] + data.unaffected
    )
    np.testing.assert_allclose(
        data.posterior_mean,
        data.posterior_alpha / (data.posterior_alpha + data.posterior_beta),
    )
    data["variant_id"] = data.unit_id
    data["canonical_variant_id"] = data.unit_id
    data["canonical_pos"] = data.aa_pos.astype(int)
    data["donor_eligible"] = True
    predictor_path = POP / "predictors/population_predictors.csv.gz"
    archive_path = EVIDENCE / "prior_count_predictors_20260912/variant_counts.csv.gz"
    scores = (
        pd.read_csv(predictor_path)
        .query("gene == @gene")
        .set_index("variant_id")
        .to_dict("index")
    )
    archived = (
        pd.read_csv(archive_path)
        .query("gene == @gene")
        .set_index("key")
        .alphamissense.to_dict()
    )
    resolved = []
    for row in data.itertuples():
        members = (
            [] if pd.isna(row.member_alleles) else str(row.member_alleles).split(";")
        )
        records = [scores[member] for member in members if member in scores]
        finite = [r for r in records if np.isfinite(r["alphamissense"])]
        pairs = {(r["alphamissense_version"], r["alphamissense"]) for r in finite}
        if any("conflict" in str(r["alphamissense_status"]) for r in records):
            resolved.append((np.nan, "member_version_or_value_conflict"))
        elif len(pairs) > 1:
            resolved.append((np.nan, "different_member_scores_or_versions"))
        elif finite and len(finite) == len(members):
            resolved.append(
                (float(finite[0]["alphamissense"]), "exact_genomic_members")
            )
        elif finite:
            resolved.append((np.nan, "partial_member_coverage"))
        elif np.isfinite(archived.get(row.literature_key, np.nan)):
            resolved.append(
                (float(archived[row.literature_key]), "archived_clinical_key_fallback")
            )
        else:
            resolved.append((np.nan, "missing"))
    data["am"] = [x[0] for x in resolved]
    data["am_source"] = [x[1] for x in resolved]
    return data, prior, paths + [prior_path, predictor_path, archive_path]


def comparison(variants, primary, context_models, sequence_weights, prior):
    """Exact variant-only outer LOO; all context normalizers recomputed.

    The class prior stays fixed. Held-out identity is removed from every training
    neighborhood, retaining other alleles at its residue. Counts do not weight
    regression labels. No feature parameter is selected from evaluation scores.
    """
    ids = variants.variant_id.tolist()
    id_to_i = {key: i for i, key in enumerate(ids)}
    y = variants.posterior_mean.to_numpy()
    am = variants.am.to_numpy()
    target_ids = primary.variant_id.tolist()
    target_index = np.array([id_to_i[key] for key in target_ids])
    d0 = primary.density.to_numpy()
    s0 = sequence_weights @ y
    cohorts = {"all_supported": np.flatnonzero(np.isfinite(d0))}
    am_common = np.flatnonzero(np.isfinite(d0 + am[target_index]))
    if len(am_common) >= 10:
        cohorts["am_common"] = am_common
    records = []
    for cohort, common in cohorts.items():
        if len(common) < 10:
            continue
        for progress, heldout in enumerate(common):
            donor_i = target_index[heldout]
            train = common[common != heldout]
            excluded = pd.concat(
                [
                    model.density(excluded_variant_ids=[ids[donor_i]])
                    for model in context_models
                ]
            )
            d = excluded.reindex(target_ids).to_numpy()
            w_i = sequence_weights[:, donor_i]
            denominator = 1 - w_i
            s = (s0 - w_i * y[donor_i]) / denominator
            # A tiny isolated IDR may lose its final donor. Drop that training
            # target from every compared model in this fold, never impute zero.
            valid_training = np.isfinite(d[train] + s[train])
            lost_support = int((~valid_training).sum())
            train = train[valid_training]
            if len(train) < 5:
                raise ValueError("Too few supported training targets after exclusion")
            features = {
                "intercept_only": (np.empty((len(train), 0)), np.empty((1, 0))),
                "density_fit": (d[train, None], d0[[heldout], None]),
                "sequence_fit": (s[train, None], s0[[heldout], None]),
            }
            if cohort == "am_common":
                local_am = am[target_index]
                features.update(
                    {
                        "am_fit": (local_am[train, None], local_am[[heldout], None]),
                        "am_plus_density": (
                            np.column_stack([local_am[train], d[train]]),
                            np.array([[local_am[heldout], d0[heldout]]]),
                        ),
                        "am_plus_sequence": (
                            np.column_stack([local_am[train], s[train]]),
                            np.array([[local_am[heldout], s0[heldout]]]),
                        ),
                    }
                )
            predictions = {
                "empirical_prior": prior["mean"],
                "density_unfitted": d0[heldout],
                "sequence_unfitted": s0[heldout],
            }
            for model, (tx, vx) in features.items():
                predictions[model] = float(
                    statistics.fractional_logistic(tx, y[target_index[train]], vx)[0]
                )
            row = variants.iloc[donor_i]
            for model, value in predictions.items():
                records.append(
                    {
                        "gene": row.gene,
                        "variant_id": row.variant_id,
                        "protein_key": row.protein_key,
                        "origin": row.origin,
                        "cohort": cohort,
                        "model": model,
                        "prediction": value,
                        "empirical_posterior": row.posterior_mean,
                        "observed_fraction": row.affected / row.n,
                        "affected": row.affected,
                        "unaffected": row.unaffected,
                        "n": row.n,
                        "training_targets": len(train),
                        "training_targets_lost_support": lost_support,
                        "eligible_donor_pool": len(variants),
                    }
                )
            if progress % 100 == 0:
                print(
                    f"{variants.gene.iloc[0]} {cohort}: LOO {progress + 1}/{len(common)}",
                    flush=True,
                )
    return pd.DataFrame(records)


def run_gene(gene, config, args, density_function, kernel):
    out = args.output_dir / gene
    out.mkdir(parents=True, exist_ok=True)
    variants, prior, sources = load_variants(gene)
    shards(variants, out, "input_variants", rows=1500)
    ids = variants.variant_id.tolist()
    id_i = {key: i for i, key in enumerate(ids)}
    y = variants.posterior_mean.to_numpy()
    rng = np.random.default_rng(20260912)
    beta_draws = rng.beta(
        variants.posterior_alpha.to_numpy()[:, None],
        variants.posterior_beta.to_numpy()[:, None],
        size=(len(ids), args.draws),
    )
    tables, models, weight_rows = [], [], []
    primary_geometry = None
    exclusion_checks = []
    for frame in config["frames"]:
        path = HERE / frame["path"]
        geometry = pd.read_csv(path).fillna({"idr_segment": ""})
        sources.append(path)
        metric_path = frame.get("ca_path")
        if metric_path:
            sources.append(HERE / metric_path)
        for metric in ["com", "ca"]:
            g = (
                pd.read_csv(HERE / metric_path).fillna({"idr_segment": ""})
                if metric == "ca" and metric_path
                else geometry.copy()
            )
            if metric == "ca" and "ca_geometry_state" in g:
                g["geometry_state"] = g.ca_geometry_state
            for half in [2.0, 3.0, 5.0]:
                scenario = f"{frame['name']}_{metric}_h{half:g}"
                is_primary = scenario == config["primary"]
                if is_primary:
                    primary_geometry = g
                eligible_pos = set(
                    g.loc[g.geometry_state.isin(["structured", "idr"]), "canonical_pos"]
                )
                active_ids = variants.loc[
                    variants.canonical_pos.isin(eligible_pos), "variant_id"
                ].tolist()
                summaries = []
                for start in range(0, len(active_ids), args.batch_size):
                    batch = active_ids[start : start + args.batch_size]
                    result = density_function(
                        variants,
                        g,
                        half_distance=half,
                        metric=metric,
                        target_ids=batch,
                        include_context_weights=False,
                        include_context_model=is_primary,
                    )
                    summary = result.summary.drop(
                        columns=["canonical_pos"], errors="ignore"
                    ).copy()
                    weights = result.donor_weights.reindex(
                        index=batch, columns=ids, fill_value=0
                    ).to_numpy()
                    for local_i, key in enumerate(batch):
                        assert weights[local_i, id_i[key]] == 0
                    np.testing.assert_allclose(
                        weights.sum(axis=1),
                        summary.density.notna().astype(int),
                        atol=1e-12,
                    )
                    density = weights @ y
                    density[weights.sum(axis=1) == 0] = np.nan
                    np.testing.assert_allclose(
                        summary.density, density, atol=1e-12, equal_nan=True
                    )
                    if is_primary:
                        models.append(result.context_model)
                        sampled = weights @ beta_draws
                        bounds = np.quantile(sampled, [0.025, 0.975], axis=1)
                        summary["density_lower_95"], summary["density_upper_95"] = (
                            bounds
                        )
                        summary["density_conditional_variance"] = (
                            weights**2
                        ) @ variants.posterior_variance.to_numpy()
                        summary.loc[
                            summary.density.isna(),
                            [
                                "density_lower_95",
                                "density_upper_95",
                                "density_conditional_variance",
                            ],
                        ] = np.nan
                        row_idx, col_idx = np.nonzero(weights)
                        weight_rows.append(
                            pd.DataFrame(
                                {
                                    "target_id": np.array(batch)[row_idx],
                                    "donor_id": np.array(ids)[col_idx],
                                    "normalized_weight": weights[row_idx, col_idx],
                                }
                            )
                        )
                    summaries.append(summary)
                complete = variants[
                    [
                        "gene",
                        "variant_id",
                        "protein_key",
                        "origin",
                        "canonical_pos",
                        "posterior_mean",
                        "am",
                    ]
                ].merge(
                    pd.concat(summaries),
                    on="variant_id",
                    how="left",
                    validate="one_to_one",
                )
                (
                    complete["scenario"],
                    complete["frame"],
                    complete["metric"],
                    complete["half_distance"],
                ) = scenario, frame["name"], metric, half
                complete.loc[complete.density.isna(), "density_source"] = "unavailable"
                tables.append(complete)
                print(
                    f"{gene} {scenario}: {complete.density.notna().sum()}/{len(variants)} supported",
                    flush=True,
                )
    scenarios = pd.concat(tables, ignore_index=True)
    primary = scenarios.loc[scenarios.scenario.eq(config["primary"])].reset_index(
        drop=True
    )
    shards(scenarios, out, "density_scenarios")
    shards(primary, out, "primary_density", rows=2000)
    weights_long = pd.concat(weight_rows, ignore_index=True)
    shards(weights_long, out, "primary_weights", rows=35000)
    del weight_rows
    active_primary = primary.loc[
        primary.variant_id.isin(pd.concat([model.density() for model in models]).index)
    ].reset_index(drop=True)
    for heldout in active_primary.variant_id.iloc[
        [0, len(active_primary) // 2]
    ].tolist():
        check_ids = active_primary.variant_id.iloc[
            :: max(1, len(active_primary) // 5)
        ].tolist()[:5]
        actual = density_function(
            variants,
            primary_geometry,
            target_ids=check_ids,
            excluded_variant_ids=[heldout],
            include_context_weights=False,
            backend="reference",
        )
        predicted = pd.concat(
            [model.density(excluded_variant_ids=[heldout]) for model in models]
        ).reindex(check_ids)
        np.testing.assert_allclose(
            actual.summary.density, predicted, atol=1e-11, equal_nan=True
        )
        exclusion_checks.append(heldout)
    positions = variants.canonical_pos.to_numpy()
    target_positions = active_primary.canonical_pos.to_numpy()
    sequence = kernel(
        3.8 * np.sqrt(abs(target_positions[:, None] - positions[None, :]))
    )
    for row_i, key in enumerate(active_primary.variant_id):
        sequence[row_i, id_i[key]] = 0
    sequence /= sequence.sum(axis=1, keepdims=True)
    predictions = comparison(variants, active_primary, models, sequence, prior)
    shards(predictions, out, "loo_predictions")
    metric_tables = []
    for cohort, rows in predictions.groupby("cohort"):
        metrics = statistics.comparison_metrics(
            rows.assign(support_stratum="all"), prior["strength"]
        )
        metrics.loc[
            metrics.model.eq("intercept_only"), "spearman_observed_fraction"
        ] = np.nan
        metrics["gene"], metrics["cohort"] = gene, cohort
        metric_tables.append(metrics)
    metrics = pd.concat(metric_tables, ignore_index=True)
    save(metrics, out / "loo_metrics.csv")
    sensitivities = []
    for scenario, rows in scenarios.groupby("scenario"):
        pair = (
            primary[["variant_id", "density"]]
            .merge(
                rows[["variant_id", "density"]],
                on="variant_id",
                suffixes=("_primary", "_alternative"),
            )
            .dropna()
        )
        delta = abs(pair.density_primary - pair.density_alternative)
        sensitivities.append(
            {
                "gene": gene,
                "scenario": scenario,
                "supported": int(rows.density.notna().sum()),
                "shared": len(pair),
                "median_density": rows.density.median(),
                "mean_absolute_change": delta.mean(),
                "max_absolute_change": delta.max(),
            }
        )
    save(pd.DataFrame(sensitivities), out / "sensitivity.csv")
    checks = {
        "gene": gene,
        "prior": prior,
        "variants": len(variants),
        "population_only": int(variants.origin.eq("population_only").sum()),
        "primary": config["primary"],
        "supported": int(primary.density.notna().sum()),
        "median_density": primary.density.median(),
        "median_kish_donor_n": primary.kish_donor_n.median(),
        "normalized_nonzero_donor_pairs": len(weights_long),
        "global_exclusion_reference_checks": exclusion_checks,
        "am_sources": variants.am_source.value_counts().to_dict(),
        "variant_loo": True,
        "gene_missense_prior_fixed": True,
        "gnomad_assumed_unaffected": True,
        "count_multiplier": False,
        "draws": args.draws,
        "seed": 20260912,
        "scope": config["scope"],
        "inputs": {str(p.relative_to(EVIDENCE)): sha(p) for p in sources},
    }
    (out / "checks.json").write_text(
        json.dumps(checks, indent=2, allow_nan=False) + "\n"
    )
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genes", nargs="+", default=["HNF1A", "LDLR", "KCNQ1", "BRCA2"]
    )
    parser.add_argument("--config", type=Path, default=HERE / "geometry_config.json")
    parser.add_argument("--output-dir", type=Path, default=HERE / "analysis")
    parser.add_argument(
        "--ppa-src", type=Path, default=REPO.parent / "ProteinProximityAnalysis/src"
    )
    parser.add_argument("--draws", type=int, default=8192)
    parser.add_argument("--batch-size", type=int, default=128)
    args = parser.parse_args()
    if args.draws < 1000 or args.batch_size < 1:
        raise ValueError("Invalid simulation or batch setting")
    os.environ.setdefault("MPLCONFIGDIR", str(REPO / "tmp/matplotlib"))
    sys.path.insert(0, str(args.ppa_src))
    from alphafold_rin.empirical_density import (
        empirical_variant_density,
        sigmoid_kernel,
    )

    config = json.loads(args.config.read_text())
    for gene in args.genes:
        run_gene(gene, config[gene], args, empirical_variant_density, sigmoid_kernel)


if __name__ == "__main__":
    main()
