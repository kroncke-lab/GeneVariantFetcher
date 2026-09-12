"""Replay the frozen GCK structural workflow with a missense-only prior.

The earlier runner/geometry/identities/AM scores remain unchanged on disk. This
adapter supplies class-specific posteriors and labels, then uses its exact
20-scenario, variant-only outer-LOO implementation in a new output directory.
"""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "population_inclusive_penetrance_20260912"
GEOMETRY_SOURCE = HERE.parent / "gck_structural_pilot_20260912"


def import_previous(name, path):
    sys.path.insert(0, str(PREVIOUS))
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


replay = import_previous("frozen_population_structure", PREVIOUS / "run_structure.py")
digest = replay.digest
MAX_BYTES = replay.MAX_BYTES


def load_inputs():
    paths = sorted((HERE / "analysis/empirical_posteriors").glob("GCK.part*.csv.gz"))
    if not paths:
        raise ValueError("Class-matched GCK empirical posterior shards are required")
    raw = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    if "variant_type" not in raw:
        raise ValueError(
            "An explicit variant_type is required to separate missense and nonsense"
        )
    full = raw.loc[raw.variant_type.eq("missense")].copy()
    if (
        not full.vclass.eq("missense").all()
        or not full.canonical_wt_status.eq("match").all()
    ):
        raise ValueError(
            "Only canonical-WT-matching missense can enter this structural run"
        )
    if ((full.aa_pos == 1) & full.aa_ref.eq("M")).any():
        raise ValueError("Start-loss cannot enter the missense structural donor pool")
    prior_path = HERE / "analysis/empirical_prior_comparison.csv"
    table = pd.read_csv(prior_path)
    prior = table.loc[table.gene.eq("GCK") & table.scope.eq("canonical_missense")]
    if len(prior) != 1:
        raise ValueError(
            "Require one canonical_missense GCK prior; no broader-class fallback"
        )
    moments = prior.iloc[0].to_dict()
    template_path = PREVIOUS / "structural/GCK_structural_input_variants.csv.gz"
    template = pd.read_csv(template_path)
    if (
        full.unit_id.duplicated().any()
        or set(full.unit_id) != set(template.unit_id)
        or len(full) != 634
    ):
        raise ValueError("The frozen 634 missense identities must remain unchanged")
    full = (
        full.set_index("unit_id", drop=False)
        .loc[template.unit_id]
        .reset_index(drop=True)
    )
    for column in [
        "gene",
        "unit_id",
        "literature_key",
        "protein_key",
        "origin",
        "member_alleles",
        "aa_pos",
        "aa_ref",
        "aa_alt",
    ]:
        pd.testing.assert_series_equal(
            full[column].fillna("").reset_index(drop=True),
            template[column].fillna("").reset_index(drop=True),
            check_names=False,
        )
    for column in [
        "affected",
        "unaffected_literature",
        "gnomad_carriers",
        "gnomad_ac",
        "unaffected",
        "n",
    ]:
        np.testing.assert_array_equal(full[column], template[column])
    np.testing.assert_allclose(full.alpha_empirical, moments["alpha_empirical"])
    np.testing.assert_allclose(full.beta_empirical, moments["beta_empirical"])
    np.testing.assert_allclose(
        full.posterior_alpha, full.alpha_empirical + full.affected
    )
    np.testing.assert_allclose(
        full.posterior_beta,
        full.beta_empirical + full.unaffected_literature + full.gnomad_carriers,
    )
    np.testing.assert_allclose(
        full.posterior_mean,
        full.posterior_alpha / (full.posterior_alpha + full.posterior_beta),
    )
    variants = full.copy()
    for column in [
        "variant_id",
        "canonical_variant_id",
        "canonical_pos",
        "donor_eligible",
        "am",
        "am_source",
        "am_old_archived",
    ]:
        variants[column] = template[column].to_numpy()
    old_predictions_path = GEOMETRY_SOURCE / "analysis/GCK_variant_loo_predictions.csv"
    old_keys = set(pd.read_csv(old_predictions_path).key)
    controlled = variants.loc[variants.literature_key.isin(old_keys)]
    if len(controlled) != 242 or set(controlled.literature_key) != old_keys:
        raise ValueError(
            "The original 242 clinical comparison identities must remain unchanged"
        )
    sources = [
        *paths,
        prior_path,
        template_path,
        PREVIOUS / "run_structure.py",
        PREVIOUS / "structure_statistics.py",
        GEOMETRY_SOURCE / "pilot_statistics.py",
        old_predictions_path,
        Path(__file__),
        *sorted((GEOMETRY_SOURCE / "geometry").glob("*_canonical_geometry.csv")),
        GEOMETRY_SOURCE / "geometry/1V4T_missing_loop_control.csv",
    ]
    hashes = {str(path.relative_to(HERE.parent)): digest(path) for path in sources}
    return full, variants, moments, controlled.variant_id.tolist(), hashes


def make_plots(variants, primary, scenarios, metrics, out, moments):
    os.environ.setdefault("MPLCONFIGDIR", str(HERE.parents[2] / "tmp/matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), constrained_layout=True)
    fig.suptitle("GCK structural density · missense-only empirical prior", fontsize=17)
    ax = axes[0, 0]
    bins = np.linspace(0, 1, 21)
    ax.hist(
        primary.posterior_mean,
        bins=bins,
        alpha=0.65,
        color="#86a7c1",
        label="Own empirical posterior",
    )
    ax.hist(
        primary.density.dropna(),
        bins=bins,
        histtype="step",
        linewidth=2,
        color="#a64e2c",
        label="Variant-excluded density",
    )
    ax.axvline(
        moments["mean"],
        color="#263c50",
        linestyle="--",
        label="Missense-only prior mean",
    )
    ax.set(
        xlabel="Probability / density",
        ylabel="Missense variants",
        title="One variant class throughout",
    )
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    colors = {
        "population_only": "#2c7f93",
        "literature_only": "#b65a37",
        "literature_and_population": "#7564a6",
    }
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
        alpha=0.15,
        linewidth=0.7,
        color="#527991",
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
        title=f"State sensitivity · {len(pair)} shared variants",
    )
    fig.savefig(out / "GCK_CLASS_MATCHED_STRUCTURAL_DENSITY.png", dpi=160)
    plt.close(fig)
    labels = {
        "empirical_prior": "Missense prior",
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
                    "#a64e2c" if model == "am_plus_density" else "#63869a"
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
        "Variant-only outer LOO · fixed missense-only empirical hyperparameters",
        fontsize=15,
    )
    fig.savefig(out / "GCK_STRUCTURAL_LOO_COMPARISON.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--output-dir", type=Path, default=HERE / "structural")
    options, _ = parser.parse_known_args()
    replay.HERE = HERE
    replay.load_inputs = load_inputs
    replay.make_plots = make_plots
    replay.main()
    out = options.output_dir.resolve()
    check_path = out / "run_checks.json"
    checks = json.loads(check_path.read_text())
    checks["missense_prior_count_units"] = checks.pop("full_locus_count_units")
    checks["missense_hyperparameters_fixed"] = checks.pop(
        "full_locus_hyperparameters_fixed"
    )
    checks["prior_scope"] = "canonical_missense"
    checks["nonsense_structural_donors"] = 0
    checks["same_counts_identity_AM_geometry_as_population_run"] = True
    checks["original_242_control"] = (
        "Same original clinical fitting/evaluation targets and archived AM scores; same expanded donors and counts; shared prior now calibrated only on canonical missense units."
    )
    checks["reused_runner_sha256"] = digest(PREVIOUS / "run_structure.py")
    checks["adapter_sha256"] = digest(__file__)
    check_path.write_text(json.dumps(checks, indent=2) + "\n")


if __name__ == "__main__":
    main()
