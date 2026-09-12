"""Run the population-inclusive GCK structural analysis from frozen snapshots."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from structure_statistics import (
    comparison_metrics,
    density_values,
    exclude_donor_values,
    variant_loo_comparison,
)


HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "gck_structural_pilot_20260912"
REPO = HERE.parents[2]
PRIMARY = "1V4S_com_h3"
MAX_BYTES = 1_200_000


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_csv(frame, path):
    text = frame.to_csv(index=False, lineterminator="\n").encode()
    data = gzip.compress(text, mtime=0) if path.suffix == ".gz" else text
    if len(data) > MAX_BYTES:
        raise ValueError(
            f"Artifact exceeds size limit: {path.name} ({len(data)} bytes)"
        )
    path.write_bytes(data)


def save_shards(frame, folder, stem, rows=4000):
    paths = []
    for start in range(0, len(frame), rows):
        part = frame.iloc[start : start + rows]
        path = folder / f"{stem}.part{len(paths) + 1:03d}.csv.gz"
        save_csv(part, path)
        paths.append(path)
    return paths


def resolve_am(row, population_scores, archived):
    members = [] if pd.isna(row.member_alleles) else str(row.member_alleles).split(";")
    records = [
        population_scores[member] for member in members if member in population_scores
    ]
    if any("conflict" in str(record["alphamissense_status"]) for record in records):
        return np.nan, "member_version_or_value_conflict"
    available = [record for record in records if np.isfinite(record["alphamissense"])]
    pairs = {
        (record["alphamissense_version"], record["alphamissense"])
        for record in available
    }
    if len(pairs) > 1:
        return np.nan, "different_member_scores_or_versions"
    if available and len(available) == len(members):
        return float(available[0]["alphamissense"]), "exact_genomic_members"
    if available:
        return np.nan, "partial_member_coverage"
    if pd.notna(row.literature_key) and row.literature_key in archived:
        value = archived[row.literature_key]
        if np.isfinite(value):
            return float(value), "archived_clinical_key_fallback"
    return np.nan, "missing"


def load_inputs():
    shards = sorted((HERE / "analysis/empirical_posteriors").glob("GCK.part*.csv.gz"))
    if not shards:
        raise ValueError("Final full-locus GCK posterior shards are required")
    full = pd.concat([pd.read_csv(path) for path in shards], ignore_index=True)
    priors_path = HERE / "analysis/empirical_prior_comparison.csv"
    priors = pd.read_csv(priors_path)
    prior = priors.loc[
        priors.gene.eq("GCK") & priors.scope.eq("all_observed_full_locus")
    ]
    if len(prior) != 1:
        raise ValueError(
            "Require one full-locus GCK empirical prior; no footprint fallback"
        )
    moments = prior.iloc[0].to_dict()
    if full.unit_id.duplicated().any() or not full.gene.eq("GCK").all():
        raise ValueError("Nonunique or incorrect GCK identities")
    for column in ["affected", "unaffected", "n"]:
        values = full[column].to_numpy()
        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or np.any(values != np.floor(values))
        ):
            raise ValueError(f"Invalid count column: {column}")
    if not (full.n > 0).all():
        raise ValueError("Unobserved alleles cannot become empirical negatives")
    np.testing.assert_allclose(
        full.unaffected, full.unaffected_literature + full.gnomad_carriers
    )
    np.testing.assert_allclose(full.n, full.affected + full.unaffected)
    np.testing.assert_allclose(full.alpha_empirical, moments["alpha_empirical"])
    np.testing.assert_allclose(full.beta_empirical, moments["beta_empirical"])
    np.testing.assert_allclose(
        full.posterior_alpha, full.alpha_empirical + full.affected
    )
    np.testing.assert_allclose(
        full.posterior_beta, full.beta_empirical + full.unaffected
    )
    np.testing.assert_allclose(
        full.posterior_mean,
        full.posterior_alpha / (full.posterior_alpha + full.posterior_beta),
    )
    eligible = full.vclass.eq("missense") & full.canonical_wt_status.eq("match")
    eligible &= ~((full.aa_pos == 1) & full.aa_ref.eq("M"))
    variants = full.loc[eligible].copy().reset_index(drop=True)
    variants["variant_id"] = variants.unit_id
    variants["canonical_variant_id"] = variants.unit_id
    variants["canonical_pos"] = variants.aa_pos.astype(int)
    variants["donor_eligible"] = True
    predictor_path = HERE / "predictors/population_predictors.csv.gz"
    scores = pd.read_csv(predictor_path).query("gene == 'GCK'")
    scores = scores.set_index("variant_id").to_dict("index")
    archive_path = OLD / "eligibility/archived_alphamissense.csv"
    archived = pd.read_csv(archive_path).set_index("key").alphamissense.to_dict()
    resolved = [resolve_am(row, scores, archived) for row in variants.itertuples()]
    variants["am"] = [pair[0] for pair in resolved]
    variants["am_source"] = [pair[1] for pair in resolved]
    variants["am_old_archived"] = variants.literature_key.map(archived)
    old_predictions_path = OLD / "analysis/GCK_variant_loo_predictions.csv"
    old_keys = set(pd.read_csv(old_predictions_path).key)
    controlled = variants.loc[variants.literature_key.isin(old_keys)]
    if set(controlled.literature_key) != old_keys or len(controlled) != 242:
        raise ValueError(
            "The requested original 242 clinical comparison targets are not preserved"
        )
    sources = [
        *shards,
        priors_path,
        predictor_path,
        archive_path,
        old_predictions_path,
        OLD / "pilot_statistics.py",
        HERE / "structure_statistics.py",
        Path(__file__),
        *sorted((OLD / "geometry").glob("*_canonical_geometry.csv")),
        OLD / "geometry/1V4T_missing_loop_control.csv",
    ]
    hashes = {str(path.relative_to(HERE.parent)): digest(path) for path in sources}
    return full, variants, moments, controlled.variant_id.tolist(), hashes


def make_plots(variants, primary, scenarios, metrics, out, moments):
    os.environ.setdefault("MPLCONFIGDIR", str(REPO / "tmp/matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    colors = {
        "population_only": "#2c7f93",
        "literature_only": "#b65a37",
        "literature_and_population": "#7564a6",
    }
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), constrained_layout=True)
    fig.suptitle("GCK structural density with population-only variants", fontsize=17)
    ax = axes[0, 0]
    bins = np.linspace(0, 1, 21)
    ax.hist(
        primary.posterior_mean,
        bins=bins,
        color="#9ab4c9",
        alpha=0.75,
        label="Own empirical posterior",
    )
    ax.hist(
        primary.density.dropna(),
        bins=bins,
        histtype="step",
        linewidth=2,
        color="#a94b2d",
        label="Variant-excluded density",
    )
    ax.axvline(
        moments["mean"],
        color="#263c50",
        linestyle="--",
        label="Full-locus empirical prior",
    )
    ax.set(
        xlabel="Probability / density",
        ylabel="Eligible missense variants",
        title="Same variant pool; different evidence",
    )
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for origin, group in primary.groupby("origin"):
        ax.scatter(
            group.posterior_mean,
            group.density,
            s=15,
            alpha=0.7,
            color=colors[origin],
            label=origin.replace("_", " "),
        )
    ax.plot([0, 1], [0, 1], ":", color="gray")
    ax.set(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="Own empirical posterior mean",
        ylabel="Variant-excluded density",
        title="Population and clinical sources",
    )
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    valid = primary.dropna(subset=["density"])
    ax.vlines(
        valid.aa_pos,
        valid.density_lower_95,
        valid.density_upper_95,
        color="#537c91",
        alpha=0.15,
        linewidth=0.7,
    )
    ax.scatter(valid.aa_pos, valid.density, s=10, color="#315d72")
    ax.set(
        xlim=(1, 465),
        ylim=(0, 1),
        xlabel="Canonical GCK residue",
        ylabel="Density and conditional 95% interval",
        title="1V4S biological monomer · COM · h = 3 Å",
    )
    ax = axes[1, 1]
    other = scenarios.query("scenario == '1V4T_com_h3'")[["variant_id", "density"]]
    pair = (
        primary[["variant_id", "density"]]
        .merge(other, on="variant_id", suffixes=("_closed", "_open"))
        .dropna()
    )
    ax.scatter(
        pair.density_closed, pair.density_open, s=15, color="#98633e", alpha=0.65
    )
    ax.plot([0, 1], [0, 1], ":", color="gray")
    ax.set(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="1V4S density",
        ylabel="1V4T density",
        title=f"Conformation sensitivity · {len(pair)} shared variants",
    )
    fig.savefig(out / "GCK_POPULATION_STRUCTURAL_DENSITY.png", dpi=160)
    plt.close(fig)
    labels = {
        "empirical_prior": "Empirical prior",
        "intercept_only": "Fitted intercept",
        "density_unfitted": "Density (unfitted)",
        "density_fit": "Density fit",
        "am_fit": "AlphaMissense",
        "am_plus_density": "AM + density",
        "sequence_fit": "Sequence control",
        "am_plus_sequence": "AM + sequence",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    for row, cohort in enumerate(["all_common", "original_242"]):
        table = (
            metrics.loc[metrics.cohort.eq(cohort) & metrics.support_stratum.eq("all")]
            .set_index("model")
            .loc[list(labels)]
        )
        for ax, metric, title in zip(
            axes[row],
            [
                "mae_empirical_posterior",
                "mse_observed_fraction",
                "mean_beta_binomial_nll",
            ],
            [
                "MAE vs empirical posterior",
                "MSE vs observed fraction",
                "Mean count negative log score",
            ],
        ):
            ax.barh(
                range(len(table)),
                table[metric],
                color=[
                    "#a94b2d" if model == "am_plus_density" else "#63869a"
                    for model in table.index
                ],
            )
            ax.set_yticks(
                range(len(table)), [labels[model] for model in table.index], fontsize=9
            )
            ax.invert_yaxis()
            ax.set(
                title=title,
                xlabel=f"{cohort.replace('_', ' ')} · {int(table.n_variants.iloc[0])} targets · lower is better",
            )
    fig.suptitle(
        "Variant-only outer LOO · fixed full-locus empirical hyperparameters",
        fontsize=15,
    )
    fig.savefig(out / "GCK_STRUCTURAL_LOO_COMPARISON.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "structural")
    parser.add_argument(
        "--ppa-src", type=Path, default=REPO.parent / "ProteinProximityAnalysis/src"
    )
    parser.add_argument("--draws", type=int, default=8192)
    args = parser.parse_args()
    if args.draws < 1000:
        raise ValueError("At least 1000 conditional draws required")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.ppa_src.resolve()))
    from alphafold_rin.empirical_density import (
        empirical_variant_density,
        sigmoid_kernel,
    )

    full, variants, moments, old_ids, hashes = load_inputs()
    ids = variants.variant_id.tolist()
    y = variants.posterior_mean.to_numpy()
    save_csv(variants, out / "GCK_structural_input_variants.csv.gz")
    rng = np.random.default_rng(20260912)
    draws = rng.beta(
        variants.posterior_alpha.to_numpy()[:, None],
        variants.posterior_beta.to_numpy()[:, None],
        size=(len(variants), args.draws),
    )
    tables, geometries = [], {}
    primary_result = primary_weights = None
    for structure in ["1V4S", "1V4T", "AF_P35557"]:
        geometry = pd.read_csv(OLD / f"geometry/{structure}_canonical_geometry.csv")
        if len(geometry) != 465 or geometry.canonical_pos.duplicated().any():
            raise ValueError(
                "Expected the frozen 465-residue canonical monomer manifest"
            )
        geometries[structure] = geometry
        for metric in ["com", "ca"]:
            for half in [2.0, 3.0, 5.0]:
                scenario = f"{structure}_{metric}_h{half:g}"
                result = empirical_variant_density(
                    variants,
                    geometry,
                    half_distance=half,
                    metric=metric,
                    include_context_weights=scenario == PRIMARY,
                )
                weights = result.donor_weights.reindex(
                    index=ids, columns=ids, fill_value=0
                ).to_numpy()
                np.testing.assert_array_equal(np.diag(weights), 0)
                summary = variants[
                    [
                        "variant_id",
                        "unit_id",
                        "literature_key",
                        "protein_key",
                        "origin",
                        "aa_pos",
                        "posterior_mean",
                    ]
                ].merge(result.summary, on="variant_id", validate="one_to_one")
                summary["scenario"], summary["structure"], summary["metric"] = (
                    scenario,
                    structure,
                    metric,
                )
                np.testing.assert_allclose(
                    summary.density, density_values(weights, y), equal_nan=True
                )
                if metric == "com" and half == 3:
                    sampled = weights @ draws
                    bounds = np.quantile(sampled, [0.025, 0.975], axis=1)
                    summary["density_lower_95"], summary["density_upper_95"] = bounds
                    summary["density_conditional_variance"] = (
                        weights**2
                    ) @ variants.posterior_variance.to_numpy()
                    summary["mc_mean_error"] = sampled.mean(axis=1) - summary.density
                    summary.loc[
                        summary.density.isna(),
                        [
                            "density_lower_95",
                            "density_upper_95",
                            "density_conditional_variance",
                        ],
                    ] = np.nan
                    if scenario == PRIMARY:
                        primary_result, primary_weights = result, weights
                tables.append(summary)
                print(
                    f"{scenario}: {summary.density.notna().sum()}/{len(variants)} supported",
                    flush=True,
                )
    for scenario, geometry, scale, exponent in [
        (
            "1V4T_loop_missing_control",
            pd.read_csv(OLD / "geometry/1V4T_missing_loop_control.csv"),
            3.8,
            0.5,
        ),
        ("1V4T_polymer_legacy_parameters", geometries["1V4T"], 5.5, 0.55),
    ]:
        result = empirical_variant_density(
            variants,
            geometry,
            polymer_scale=scale,
            polymer_exponent=exponent,
            include_context_weights=False,
        )
        summary = variants[
            [
                "variant_id",
                "unit_id",
                "literature_key",
                "protein_key",
                "origin",
                "aa_pos",
                "posterior_mean",
            ]
        ].merge(result.summary, on="variant_id", validate="one_to_one")
        summary["scenario"], summary["structure"], summary["metric"] = (
            scenario,
            "1V4T",
            "com",
        )
        tables.append(summary)
    scenarios = pd.concat(tables, ignore_index=True)
    save_shards(scenarios, out, "GCK_density_scenarios")
    primary = scenarios.loc[scenarios.scenario.eq(PRIMARY)].reset_index(drop=True)
    save_csv(primary, out / "GCK_primary_variant_density.csv")
    for start in range(0, len(ids), 40):
        context = primary_result.context_weights.loc[
            primary_result.context_weights.variant_id.isin(ids[start : start + 40])
        ]
        save_csv(
            context, out / f"GCK_primary_donor_contexts_{start // 40 + 1:03d}.csv.gz"
        )
    save_shards(
        primary_result.donor_weights.reset_index(),
        out,
        "GCK_primary_normalized_weights",
        rows=40,
    )
    far = primary_result.context_weights.query("distance > 20")
    assert len(far) and (far.kernel > 0).all()
    pos = variants.canonical_pos.to_numpy()
    sequence_weights = sigmoid_kernel(3.8 * np.sqrt(abs(pos[:, None] - pos[None, :])))
    np.fill_diagonal(sequence_weights, 0)
    sequence_weights /= sequence_weights.sum(axis=1, keepdims=True)
    predictions = [
        variant_loo_comparison(
            variants, primary_weights, sequence_weights, prior_mean=moments["mean"]
        )
    ]
    controlled_variants = variants.copy()
    controlled_variants.loc[controlled_variants.variant_id.isin(old_ids), "am"] = (
        controlled_variants.loc[
            controlled_variants.variant_id.isin(old_ids), "am_old_archived"
        ]
    )
    predictions.append(
        variant_loo_comparison(
            controlled_variants,
            primary_weights,
            sequence_weights,
            prior_mean=moments["mean"],
            target_ids=old_ids,
            cohort="original_242",
        )
    )
    predictions = pd.concat(predictions, ignore_index=True)
    support = primary.set_index("variant_id").kish_donor_n
    predictions["kish_donor_n"] = predictions.variant_id.map(support)
    metric_tables = []
    for cohort, group in predictions.groupby("cohort"):
        all_rows = group.assign(support_stratum="all")
        kish = group.assign(
            support_stratum=np.where(group.kish_donor_n < 5, "kish_lt_5", "kish_ge_5")
        )
        origins = group.assign(support_stratum="origin:" + group.origin)
        result = comparison_metrics(
            pd.concat([all_rows, kish, origins]), moments["strength"]
        )
        result["cohort"] = cohort
        metric_tables.append(result)
    metrics = pd.concat(metric_tables, ignore_index=True)
    save_csv(predictions, out / "GCK_variant_loo_predictions.csv.gz")
    save_csv(metrics, out / "GCK_variant_loo_metrics.csv")
    old_metrics = pd.read_csv(OLD / "analysis/GCK_variant_loo_metrics.csv")
    controlled = metrics.query(
        "cohort == 'original_242' and support_stratum == 'all'"
    ).merge(
        old_metrics.query("support_stratum == 'all'"),
        on="model",
        suffixes=("_population", "_old"),
        validate="one_to_one",
    )
    save_csv(controlled, out / "GCK_original_242_metric_comparison.csv")
    sensitivity = []
    for scenario, frame in scenarios.groupby("scenario"):
        paired = (
            primary[["variant_id", "density"]]
            .merge(
                frame[["variant_id", "density"]],
                on="variant_id",
                suffixes=("_primary", "_scenario"),
            )
            .dropna()
        )
        delta = abs(paired.density_primary - paired.density_scenario)
        sensitivity.append(
            {
                "scenario": scenario,
                "supported_variants": int(frame.density.notna().sum()),
                "shared_with_primary": len(paired),
                "median_density": frame.density.median(),
                "mean_abs_change_from_primary": delta.mean(),
                "max_abs_change_from_primary": delta.max(),
            }
        )
    save_csv(pd.DataFrame(sensitivity), out / "GCK_sensitivity_summary.csv")
    exclusion_checks = []
    base = density_values(primary_weights, y)
    for i in [0, len(ids) // 2, len(ids) - 1]:
        actual = empirical_variant_density(
            variants,
            geometries["1V4S"],
            excluded_variant_ids=[ids[i]],
            include_context_weights=False,
        )
        np.testing.assert_allclose(
            actual.summary.density,
            exclude_donor_values(primary_weights, y, i, base),
            atol=1e-12,
            equal_nan=True,
        )
        exclusion_checks.append(ids[i])
    make_plots(variants, primary, scenarios, metrics, out, moments)
    for relative, expected in hashes.items():
        if digest(HERE.parent / relative) != expected:
            raise ValueError(f"Input changed during structural run: {relative}")
    checks = {
        "prior": moments,
        "full_locus_count_units": len(full),
        "eligible_missense": len(variants),
        "origin_counts": variants.origin.value_counts().to_dict(),
        "primary_supported": int(primary.density.notna().sum()),
        "comparison_targets": predictions.groupby("cohort")
        .variant_id.nunique()
        .to_dict(),
        "am_sources": variants.am_source.value_counts().to_dict(),
        "scenario_count": len(tables),
        "draws": args.draws,
        "seed": 20260912,
        "alpha_adds_affected_beta_adds_all_unaffected": True,
        "gnomad_all_assumed_unaffected": True,
        "full_locus_hyperparameters_fixed": True,
        "zero_self_weight_all_scenarios": True,
        "positive_donor_weights_beyond_20": len(far),
        "global_outer_exclusion_engine_checks": exclusion_checks,
        "spatial_weighting": "equal_variant_kernel_weights_no_carrier_multiplier",
        "identity": "union unit_id; distinct population DNA alleles retained; clinical aggregate members consolidated upstream",
        "original_242_control": "Same original clinical fitting/evaluation targets and archived AM scores; expanded donors, refreshed counts and full-locus empirical prior. This comparison changes data as explicitly stated.",
        "endpoint": "pooled_GCK_clinical_evidence_not_disease_specific",
        "interval_scope": "independent donor Beta draws conditional on fixed counts, hyperparameters, identity, and geometry",
        "ppa_module_sha256": digest(
            args.ppa_src / "alphafold_rin/empirical_density.py"
        ),
        "input_hashes": hashes,
        "output_hashes": {
            path.name: digest(path)
            for path in sorted(out.iterdir())
            if path.is_file() and path.name != "run_checks.json"
        },
    }
    (out / "run_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    if any(path.stat().st_size > MAX_BYTES for path in out.iterdir() if path.is_file()):
        raise ValueError("Oversized structural artifact")
    print(metrics.query("support_stratum == 'all'").to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
