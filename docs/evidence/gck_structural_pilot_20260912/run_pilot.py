"""Run the frozen GCK empirical-posterior structural pilot offline.

Geometry acquisition is a separate source-backed step (geometry/build_geometry.py).
This runner consumes committed count/identity/coordinate snapshots and the new
explicit PPA empirical-density API, preserving historical PPA behavior.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution

from pilot_statistics import (
    comparison_metrics,
    density_values,
    remove_donor,
    variant_loo_comparison,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PREFLIGHT = HERE.parent / "structural_density_plan_20260912"
PRIMARY = "1V4S_com_h3"


def save_csv(frame, path):
    if str(path).endswith(".gz"):
        frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
    else:
        frame.to_csv(path, index=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs():
    """Recompute the full 398-row prior before any geometry/consequence filter."""
    spec = importlib.util.spec_from_file_location(
        "empirical_preflight", PREFLIGHT / "full_empirical_universe.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frame = pd.read_csv(PREFLIGHT / "full_universe_source_counts.csv.gz")
    frame = frame.loc[frame.gene.eq("GCK")].copy().reset_index(drop=True)
    if frame.key.duplicated().any() or len(frame) != 398:
        raise ValueError("Frozen GCK universe changed")
    frame["unaffected_total"] = frame.unaffected + frame.gnomad_added
    frame["n"] = frame.affected + frame.unaffected_total
    if not np.isfinite(frame.n).all() or (frame.n <= 0).any():
        raise ValueError("Invalid carrier denominator")
    if not np.allclose(frame.affected + frame.unaffected, frame.literature_n):
        raise ValueError("Literature outcome partition differs")
    for column in ["affected", "unaffected_total", "n"]:
        if not np.equal(frame[column], np.floor(frame[column])).all():
            raise ValueError(f"Count column {column} is not integer-valued")
    moments, posterior = module.fit_moments(
        frame.affected.to_numpy(), frame.n.to_numpy(), module.PRIMARY
    )
    alpha, beta = moments["alpha_empirical"], moments["beta_empirical"]
    frame["alpha_empirical"] = alpha
    frame["beta_empirical"] = beta
    frame["alpha_posterior_empirical"] = alpha + frame.affected
    frame["beta_posterior_empirical"] = beta + frame.unaffected_total
    frame["posterior_empirical_mean"] = posterior
    frame["posterior_empirical_variance"] = (
        posterior * (1 - posterior) / (moments["strength"] + frame.n + 1)
    )
    for label, quantile in [("lower", 0.025), ("upper", 0.975)]:
        frame[f"posterior_empirical_{label}_95"] = beta_distribution.ppf(
            quantile, frame.alpha_posterior_empirical, frame.beta_posterior_empirical
        )
    old = pd.read_csv(PREFLIGHT / "full_universe_empirical_posteriors.csv.gz")
    old = old.loc[old.gene.eq("GCK")].set_index("key").loc[frame.key]
    for column in [
        "alpha_empirical",
        "beta_empirical",
        "alpha_posterior_empirical",
        "beta_posterior_empirical",
        "posterior_empirical_mean",
    ]:
        np.testing.assert_allclose(frame[column], old[column], rtol=1e-12)
    audit = pd.read_csv(HERE / "eligibility/variant_eligibility.csv")
    predictions = pd.read_csv(HERE / "eligibility/archived_alphamissense.csv")
    extra = [c for c in audit if c not in frame or c == "key"]
    frame = frame.merge(audit[extra], on="key", validate="one_to_one")
    frame = frame.merge(
        predictions[["key", "alphamissense"]].rename(columns={"alphamissense": "am"}),
        on="key",
        validate="one_to_one",
    )
    if len(frame) != 398:
        raise ValueError("Identity/predictor join lost rows")
    if not frame.geometry_eligible_pre_structure.isin([True, False]).all():
        raise ValueError("Eligibility must be explicit Boolean")
    selected = (
        frame.loc[frame.geometry_eligible_pre_structure].copy().reset_index(drop=True)
    )
    selected["variant_id"] = "GCK:" + selected.key
    selected["canonical_variant_id"] = selected.variant_id
    selected["canonical_pos"] = selected.aa_pos.astype(int)
    selected["posterior_alpha"] = selected.alpha_posterior_empirical
    selected["posterior_beta"] = selected.beta_posterior_empirical
    selected["donor_eligible"] = True
    return frame, selected, moments


def make_plots(full, primary, scenarios, metrics, out, moments):
    os.environ.setdefault("MPLCONFIGDIR", str(REPO / "tmp/matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.8), constrained_layout=True)
    fig.suptitle(
        "GCK structural pilot · pooled clinical evidence",
        fontsize=17,
        fontweight="bold",
    )
    valid = primary.loc[primary.density.notna()]
    ax = axes[0, 0]
    bins = np.linspace(0, 1, 21)
    ax.hist(
        full.posterior_empirical_mean,
        bins=bins,
        color="#9ab4c9",
        alpha=0.75,
        label=f"Empirical posterior · all {len(full)} variants",
    )
    ax.hist(
        valid.density,
        bins=bins,
        histtype="step",
        linewidth=2.2,
        color="#aa4c27",
        label=f"Variant-excluded density · {len(valid)} mapped missense",
    )
    ax.axvline(
        moments["mean"],
        color="#293d50",
        linestyle="--",
        label="Shared empirical prior mean",
    )
    ax.set(
        xlabel="Probability / density",
        ylabel="Variants",
        title="Prior, own counts and neighbor evidence",
    )
    ax.legend(fontsize=8, loc="upper left")
    ax = axes[0, 1]
    dots = ax.scatter(
        valid.posterior_empirical_mean,
        valid.density,
        c=valid.kish_donor_n,
        cmap="viridis",
        s=24,
        alpha=0.8,
    )
    ax.plot([0, 1], [0, 1], ":", color="grey")
    ax.set(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="Own empirical posterior mean",
        ylabel="Variant-excluded structural density",
        title="Nearby variants need not agree",
    )
    fig.colorbar(dots, ax=ax, label="Kish effective donor variants")
    ax = axes[1, 0]
    ordered = valid.sort_values("aa_pos")
    ax.errorbar(
        ordered.aa_pos,
        ordered.density,
        yerr=np.vstack(
            [
                ordered.density - ordered.density_lower_95,
                ordered.density_upper_95 - ordered.density,
            ]
        ),
        fmt="none",
        alpha=0.22,
        color="#3b667f",
        linewidth=0.8,
    )
    ax.scatter(ordered.aa_pos, ordered.density, s=13, color="#3b667f")
    ax.set(
        xlim=(1, 465),
        ylim=(0, 1),
        xlabel="Canonical GCK residue (P35557-1)",
        ylabel="Density with conditional 95% interval",
        title="Biological monomer 1V4S · side-chain COM · h = 3 Å",
    )
    ax = axes[1, 1]
    primary_values = scenarios.loc[scenarios.scenario.eq(PRIMARY), ["key", "density"]]
    other = scenarios.loc[scenarios.scenario.eq("1V4T_com_h3"), ["key", "density"]]
    pair = primary_values.merge(other, on="key", suffixes=("_closed", "_open")).dropna()
    ax.scatter(pair.density_closed, pair.density_open, s=21, color="#aa4c27", alpha=0.7)
    ax.plot([0, 1], [0, 1], ":", color="grey")
    ax.set(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="1V4S density",
        ylabel="1V4T density",
        title=f"Conformation sensitivity · {len(pair)} shared variants",
    )
    fig.savefig(out / "GCK_STRUCTURAL_PILOT.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    fig.suptitle(
        "Variant-only LOO diagnostics · same targets · full-dataset hyperparameters fixed",
        fontsize=13,
    )
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
    table = (
        metrics.loc[metrics.support_stratum.eq("all")]
        .set_index("model")
        .loc[list(labels)]
    )
    for ax, column, title in zip(
        axes,
        ["mae_empirical_posterior", "mse_observed_fraction", "mean_beta_binomial_nll"],
        [
            "MAE vs empirical posterior",
            "MSE vs observed fraction",
            "Mean count negative log score",
        ],
    ):
        ax.barh(
            range(len(table)),
            table[column],
            color=[
                "#aa4c27" if x == "am_plus_density" else "#688da5" for x in table.index
            ],
        )
        ax.set_yticks(range(len(table)), [labels[x] for x in table.index], fontsize=9)
        ax.invert_yaxis()
        ax.set(title=title, xlabel="Lower is better · exploratory")
    fig.savefig(out / "GCK_INTERNAL_COMPARISON.png", dpi=180)
    plt.close(fig)

    fig = plt.figure(figsize=(11.5, 5.8), constrained_layout=True)
    fig.suptitle(
        "GCK biological monomers · spatial view of variant-excluded density",
        fontsize=15,
    )
    for panel, structure in enumerate(["1V4S", "1V4T"], start=1):
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        geometry = pd.read_csv(HERE / f"geometry/{structure}_canonical_geometry.csv")
        xyz = geometry[["ca_x", "ca_y", "ca_z"]].to_numpy()
        # NaN canonical rows break the backbone at missing/disordered segments.
        ax.plot(
            xyz[:, 0], xyz[:, 1], xyz[:, 2], color="#b6bdc2", linewidth=0.7, alpha=0.7
        )
        values = (
            scenarios.loc[scenarios.scenario.eq(f"{structure}_com_h3")]
            .groupby("aa_pos")
            .density.mean()
        )
        coords = geometry.merge(
            values.rename("density"), left_on="canonical_pos", right_index=True
        ).dropna(subset=["density", "com_x"])
        dots = ax.scatter(
            coords.com_x,
            coords.com_y,
            coords.com_z,
            c=coords.density,
            cmap="viridis",
            vmin=0,
            vmax=1,
            s=24,
            depthshade=False,
        )
        ax.set_title(
            f"{structure} · mean across assayed variants at each residue", fontsize=10
        )
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
    fig.colorbar(
        dots, ax=fig.axes, shrink=0.7, label="Mean variant-excluded density (h = 3 Å)"
    )
    fig.text(
        0.08,
        0.025,
        "Measured coordinates only; the disordered loop has no invented 3D position.",
        fontsize=10,
    )
    fig.savefig(out / "GCK_SPATIAL_DENSITY.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "analysis")
    parser.add_argument(
        "--ppa-src", type=Path, default=REPO.parent / "ProteinProximityAnalysis/src"
    )
    parser.add_argument("--draws", type=int, default=8192)
    args = parser.parse_args()
    if args.draws < 1000:
        raise ValueError("Use at least 1000 draws for the pilot conditional intervals")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.ppa_src.resolve()))
    from alphafold_rin.empirical_density import empirical_variant_density

    full, variants, moments = load_inputs()
    save_csv(full, out / "GCK_empirical_posteriors_and_eligibility.csv")
    ids = variants.variant_id.tolist()
    y = variants.posterior_empirical_mean.to_numpy()
    rng = np.random.default_rng(20260912)
    draws = rng.beta(
        variants.posterior_alpha.to_numpy()[:, None],
        variants.posterior_beta.to_numpy()[:, None],
        size=(len(variants), args.draws),
    )
    scenario_tables = []
    primary_result = None
    primary_weights = None
    geometries = {}
    for structure in ["1V4S", "1V4T", "AF_P35557"]:
        geometry = pd.read_csv(HERE / f"geometry/{structure}_canonical_geometry.csv")
        geometries[structure] = geometry
        if len(geometry) != 465 or geometry.canonical_pos.duplicated().any():
            raise ValueError("Expected one 465-residue canonical monomer manifest")
        for metric in ["com", "ca"]:
            for h in [2.0, 3.0, 5.0]:
                scenario = f"{structure}_{metric}_h{h:g}"
                result = empirical_variant_density(
                    variants, geometry, half_distance=h, metric=metric
                )
                weights = result.donor_weights.reindex(
                    index=ids, columns=ids, fill_value=0
                ).to_numpy()
                np.testing.assert_allclose(np.diag(weights), 0, atol=0)
                summary = result.summary.copy()
                summary = variants[
                    ["variant_id", "key", "aa_pos", "posterior_empirical_mean"]
                ].merge(summary, on="variant_id", validate="one_to_one")
                summary["scenario"] = scenario
                summary["structure"] = structure
                summary["metric"] = metric
                summary["half_distance"] = h
                np.testing.assert_allclose(
                    summary.density, density_values(weights, y), equal_nan=True
                )
                if metric == "com" and h == 3:
                    sampled = weights @ draws
                    bounds = np.quantile(sampled, [0.025, 0.975], axis=1)
                    summary["density_lower_95"] = bounds[0]
                    summary["density_upper_95"] = bounds[1]
                    summary.loc[
                        summary.density.isna(), ["density_lower_95", "density_upper_95"]
                    ] = np.nan
                    analytic_variance = (
                        weights**2
                    ) @ variants.posterior_empirical_variance.to_numpy()
                    summary["density_conditional_variance"] = analytic_variance
                    summary.loc[
                        summary.density.isna(), "density_conditional_variance"
                    ] = np.nan
                    mc_mean = sampled.mean(axis=1)
                    summary["mc_mean_error"] = mc_mean - summary.density
                    if scenario == PRIMARY:
                        primary_result, primary_weights = result, weights
                scenario_tables.append(summary)
                print(
                    f"{scenario}: {summary.density.notna().sum()}/{len(summary)} supported",
                    flush=True,
                )
    for label, geometry, scale, exponent in [
        (
            "1V4T_loop_missing_control",
            pd.read_csv(HERE / "geometry/1V4T_missing_loop_control.csv"),
            3.8,
            0.5,
        ),
        ("1V4T_polymer_legacy_parameters", geometries["1V4T"], 5.5, 0.55),
    ]:
        result = empirical_variant_density(
            variants, geometry, polymer_scale=scale, polymer_exponent=exponent
        )
        summary = variants[
            ["variant_id", "key", "aa_pos", "posterior_empirical_mean"]
        ].merge(result.summary, on="variant_id", validate="one_to_one")
        summary["scenario"] = label
        summary["structure"] = "1V4T"
        summary["metric"] = "com"
        summary["half_distance"] = 3.0
        summary["polymer_scale"] = scale
        summary["polymer_exponent"] = exponent
        scenario_tables.append(summary)
    scenarios = pd.concat(scenario_tables, ignore_index=True)
    save_csv(scenarios, out / "GCK_density_scenarios.csv.gz")
    primary = (
        scenarios.loc[scenarios.scenario.eq(PRIMARY)].copy().reset_index(drop=True)
    )
    save_csv(primary, out / "GCK_primary_variant_density.csv")
    # The complete donor table makes the long-distance tail directly inspectable.
    # Bounded chunks keep every target-donor pair inspectable without a single
    # oversized Git artifact. No donor/weight threshold is applied.
    for start in range(0, len(ids), 80):
        chunk = primary_result.context_weights.loc[
            primary_result.context_weights.variant_id.isin(ids[start : start + 80])
        ]
        save_csv(
            chunk, out / f"GCK_primary_donor_contexts_{start // 80 + 1:02d}.csv.gz"
        )
    save_csv(
        primary_result.donor_weights.rename_axis("variant_id").reset_index(),
        out / "GCK_primary_normalized_donor_weights.csv.gz",
    )
    # Simple sequence-neighborhood comparator across the canonical protein;
    # this is a control, not the IDR-specific geometry fallback.
    pos = variants.canonical_pos.to_numpy()
    distances = 3.8 * np.sqrt(np.abs(pos[:, None] - pos[None, :]))
    exp_negative = np.exp(-np.log(3) * distances / 3)
    seq_weights = 2 * exp_negative / (1 + exp_negative)
    np.fill_diagonal(seq_weights, 0)
    seq_weights /= seq_weights.sum(axis=1, keepdims=True)
    predictions = variant_loo_comparison(
        variants, primary_weights, seq_weights, prior_mean=moments["mean"]
    )
    supports = primary.set_index("key").kish_donor_n
    predictions["kish_donor_n"] = predictions.key.map(supports)
    predictions["support_stratum"] = "all"
    stratified = predictions.copy()
    stratified["support_stratum"] = np.where(
        stratified.kish_donor_n < 5, "kish_lt_5", "kish_ge_5"
    )
    metrics = comparison_metrics(
        pd.concat([predictions, stratified]), moments["strength"]
    )
    save_csv(predictions, out / "GCK_variant_loo_predictions.csv")
    save_csv(metrics, out / "GCK_variant_loo_metrics.csv")
    calibration = predictions.copy()
    calibration["prediction_bin"] = pd.cut(
        calibration.prediction, bins=np.linspace(0, 1, 6), include_lowest=True
    )
    calibration = (
        calibration.groupby(["model", "prediction_bin"], observed=True)
        .agg(
            n_variants=("key", "size"),
            mean_prediction=("prediction", "mean"),
            mean_empirical_posterior=("empirical_posterior", "mean"),
            mean_observed_fraction=("observed_fraction", "mean"),
        )
        .reset_index()
    )
    save_csv(calibration, out / "GCK_internal_calibration.csv")
    sensitivity = []
    base = primary.set_index("key").density
    for scenario, table in scenarios.groupby("scenario"):
        pair = pd.concat(
            [base.rename("primary"), table.set_index("key").density], axis=1
        ).dropna()
        sensitivity.append(
            {
                "scenario": scenario,
                "supported_variants": int(table.density.notna().sum()),
                "shared_with_primary": len(pair),
                "median_density": table.density.median(),
                "mean_abs_change_from_primary": np.mean(
                    abs(pair.primary - pair.density)
                ),
                "max_abs_change_from_primary": np.max(abs(pair.primary - pair.density)),
            }
        )
    save_csv(pd.DataFrame(sensitivity), out / "GCK_sensitivity_summary.csv")
    # Recompute selected outer exclusions through the engine, independent of
    # the matrix shortcut used for each regression fold.
    heldout_checks = []
    for i in [0, len(variants) // 2, len(variants) - 1]:
        excluded = empirical_variant_density(
            variants, geometries["1V4S"], excluded_variant_ids=[ids[i]]
        )
        actual = excluded.donor_weights.reindex(
            index=ids, columns=ids, fill_value=0
        ).to_numpy()
        np.testing.assert_allclose(actual, remove_donor(primary_weights, i), atol=1e-12)
        heldout_checks.append(ids[i])
    checks = {
        "prior": moments,
        "full_count_variants": len(full),
        "eligible_missense": len(variants),
        "primary_supported": int(primary.density.notna().sum()),
        "same_comparison_variants": int(predictions.key.nunique()),
        "scenario_count": len(scenario_tables),
        "draws": args.draws,
        "seed": 20260912,
        "full_prior_reproduces_frozen_preflight": True,
        "zero_self_weight_all_scenarios": True,
        "outer_exclusion_matches_engine_targets": heldout_checks,
        "gnomad_assumed_unaffected": True,
        "full_dataset_hyperparameters_fixed": True,
        "endpoint": "pooled_GCK_clinical_evidence_not_disease_specific",
        "identity_grain": "archived_protein_key_aggregate",
        "interval_scope": "independent donor Beta draws conditional on fixed counts, hyperparameters and geometry",
        "ppa_module_sha256": digest(
            args.ppa_src / "alphafold_rin/empirical_density.py"
        ),
        "input_hashes": {
            str(path.relative_to(HERE.parent)): digest(path)
            for path in [
                PREFLIGHT / "full_universe_source_counts.csv.gz",
                HERE / "eligibility/variant_eligibility.csv",
                HERE / "eligibility/archived_alphamissense.csv",
                *sorted((HERE / "geometry").glob("*_canonical_geometry.csv")),
            ]
        },
    }
    (out / "run_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    make_plots(full, primary, scenarios, metrics, out, moments)
    print(
        metrics.loc[metrics.support_stratum.eq("all")].to_string(index=False),
        flush=True,
    )


if __name__ == "__main__":
    main()
