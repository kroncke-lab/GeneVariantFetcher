#!/usr/bin/env python3
"""Independent, fixed-input GCK prior/smoothing/endpoint diagnostics.

Run with numpy, pandas, scipy and matplotlib. No network or accepted-fit edits.
All numerical controls preserve the canonical 634-unit missense inventory.
Endpoint quarantine is an explicitly separate donor-only sensitivity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
BASE = REPO / "docs/evidence/class_matched_penetrance_20260912/structural"
SEED = 20260913


def save(frame, name):
    kwargs = {"index": False, "lineterminator": "\n", "float_format": "%.12g"}
    if name.endswith(".gz"):
        kwargs["compression"] = {"method": "gzip", "mtime": 0}
    frame.to_csv(HERE / name, **kwargs)


def describe(values):
    a = np.asarray(values)
    a = a[np.isfinite(a)]
    return dict(
        count=len(a),
        mean=float(a.mean()),
        sd=float(a.std(ddof=1)),
        **dict(
            zip(
                ["min", "p05", "p25", "median", "p75", "p95", "max"],
                np.quantile(a, [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]).tolist(),
            )
        ),
    )


def moments(a, n, mode):
    f = a / n
    w = np.ones(len(n)) if mode == "equal_variant" else 1 - 1 / (n + 0.01)
    mean = np.average(f, weights=w)
    denom = len(n) if mode == "historical" else w.sum()
    variance = np.sum(w * (f - mean) ** 2) / denom
    strength = mean * (1 - mean) / variance - 1
    assert strength > 0
    return mean * strength, (1 - mean) * strength, mean, variance, strength


def control_statistics(labels, weights, sequence, mask, baseline):
    """Each column is one joint donor+target label assignment; W diagonal=0."""
    if labels.ndim == 1:
        labels = labels[:, None]
    y = labels[mask]
    density = weights[mask] @ labels
    seq = sequence[mask] @ labels
    baseline_mse = np.mean((y - baseline) ** 2, axis=0)
    mse = np.mean((y - density) ** 2, axis=0)
    seq_mse = np.mean((y - seq) ** 2, axis=0)
    return {
        "mse_density": mse,
        "gain_vs_fixed_prior": baseline_mse - mse,
        "gain_vs_sequence": seq_mse - mse,
        "density_sd": np.std(density, axis=0, ddof=1),
    }


def permutations(d, w, s, mask, mu, repeats):
    # Move complete (A,U,n,posterior) records, including their target labels.
    # Shuffling donor labels alone against unchanged targets would leak labels.
    rng = np.random.default_rng(SEED)
    p, f = d.posterior_mean.to_numpy(), d.observed_fraction.to_numpy()
    eligible = np.flatnonzero(mask)
    strata = {
        "unconditional": pd.Series("all", index=d.index),
        "within_origin": d.origin,
        "within_origin_n_bin": d.origin + ":" + d.n_bin.astype(str),
    }
    observed = {
        "posterior": control_statistics(p, w, s, mask, mu),
        "count_only": control_statistics(f, w, s, mask, mu),
    }
    rows, summary = [], []
    for name, groups in strata.items():
        members = [
            eligible[groups.iloc[eligible].to_numpy() == group]
            for group in sorted(set(groups.iloc[eligible]))
        ]
        scheme = []
        for start in range(0, repeats, 100):
            size = min(100, repeats - start)
            idx = np.tile(np.arange(len(d))[:, None], (1, size))
            for col in range(size):
                for group in members:
                    idx[group, col] = rng.permutation(group)
            for outcome, labels in [("posterior", p[idx]), ("count_only", f[idx])]:
                stat = control_statistics(labels, w, s, mask, mu)
                block = pd.DataFrame(stat)
                block["scheme"], block["outcome"] = name, outcome
                block["replicate"] = np.arange(start + 1, start + size + 1)
                scheme.append(block)
        scheme = pd.concat(scheme, ignore_index=True)
        rows.append(scheme)
        for outcome in observed:
            for metric, value in observed[outcome].items():
                value = float(value[0])
                null = scheme.loc[scheme.outcome.eq(outcome), metric].to_numpy()
                summary.append(
                    dict(
                        scheme=name,
                        outcome=outcome,
                        metric=metric,
                        observed=value,
                        null_mean=float(null.mean()),
                        null_p025=float(np.quantile(null, 0.025)),
                        null_p975=float(np.quantile(null, 0.975)),
                        p_null_at_least_observed=float(
                            (1 + (null >= value).sum()) / (repeats + 1)
                        ),
                        repeats=repeats,
                        strata=len(members),
                    )
                )
    save(pd.concat(rows, ignore_index=True), "permutation_replicates.csv.gz")
    save(pd.DataFrame(summary), "permutation_summary.csv")
    return pd.DataFrame(summary)


def figure(d, mask, mu, summary):
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 160,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    v = d.loc[mask]
    ax.scatter(
        d.canonical_pos,
        d.posterior_mean,
        s=11,
        alpha=0.45,
        color="#767d87",
        label="Own count-updated posterior (634)",
    )
    ax.scatter(
        v.canonical_pos,
        v.density_reconstructed,
        s=13,
        alpha=0.72,
        color="#1865a2",
        label="Neighborhood score; own variant excluded (614)",
    )
    ax.axhline(
        mu,
        color="#b84a32",
        ls="--",
        lw=1.5,
        label=f"Shared missense prior mean: {mu:.1%}",
    )
    ax.set(
        xlabel="Canonical GCK residue",
        ylabel="Posterior or neighborhood score",
        ylim=(0, 1),
        title="A  The neighborhood score is not variant penetrance",
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax = axes[0, 1]
    keys = ["D217N", "E279Q", "W99R", "V389L"]
    examples = d.set_index("protein_key").loc[keys]
    x = np.arange(len(keys))
    ax.bar(
        x - 0.18,
        examples.posterior_mean,
        width=0.36,
        color="#767d87",
        label="Own posterior",
    )
    ax.bar(
        x + 0.18,
        examples.density_reconstructed,
        width=0.36,
        color="#1865a2",
        label="Neighborhood score",
    )
    for i, (_, row) in enumerate(examples.iterrows()):
        for dx, col in [(-0.18, "posterior_mean"), (0.18, "density_reconstructed")]:
            ax.text(i + dx, row[col] + 0.02, f"{row[col]:.1%}", ha="center", fontsize=8)
    ax.set(
        xticks=x,
        xticklabels=[
            "D217N\nWT-like",
            "E279Q\nWT-like",
            "W99R\nActivating / HH",
            "V389L\nActivating / HH",
        ],
        ylim=(0, 1),
        ylabel="Frozen model value",
        title="B  Functional controls expose the interpretation error",
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax = axes[1, 0]
    curves = [
        ("posterior_mean", "Own posterior", "#767d87"),
        ("density_reconstructed", "Posterior neighborhood", "#1865a2"),
        ("count_only_density", "Count-only neighborhood (diagnostic)", "#d28227"),
    ]
    bins = np.linspace(0, 1, 26)
    for col, label, color in curves:
        ax.hist(
            v[col],
            bins=bins,
            weights=np.ones(len(v)) / len(v),
            histtype="step",
            lw=2,
            label=label,
            color=color,
        )
    ax.axvline(mu, color="#b84a32", ls="--", lw=1)
    ax.set(
        xlabel="Posterior, score, or count fraction",
        ylabel="Fraction of supported variants",
        title="C  Prior shrinkage compresses neighborhood differences",
    )
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(fontsize=8, frameon=False)
    ax = axes[1, 1]
    z = summary.query("outcome == 'posterior' and metric == 'gain_vs_fixed_prior'")
    null = pd.read_csv(HERE / "permutation_replicates.csv.gz")
    for name, color, label in [
        ("unconditional", "#767d87", "Shuffle all supported labels"),
        ("within_origin", "#d28227", "Shuffle within source origin"),
        ("within_origin_n_bin", "#3d9477", "Shuffle within origin + count bin"),
    ]:
        vals = null.loc[
            null.scheme.eq(name) & null.outcome.eq("posterior"), "gain_vs_fixed_prior"
        ]
        ax.hist(vals, bins=40, histtype="step", lw=1.7, color=color, label=label)
    ax.axvline(z.observed.iloc[0], color="#1865a2", lw=2, label="Observed GCK")
    ax.set(
        xlabel="Posterior MSE improvement over constant prior",
        ylabel="Permutation replicates",
        title="D  Does geography add information beyond source selection?",
    )
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        "GCK sanity check: shared prior + mixed clinical counts + local smoothing",
        fontsize=15,
        y=0.99,
    )
    fig.text(
        0.02,
        0.015,
        "Frozen 634 canonical missense units; gnomAD carriers assumed unaffected. Primary: 1V4S side-chain COM, sigmoid half-weight 3 Å, positive tail.\n"
        "HH = hyperinsulinemic hypoglycemia. These pooled model values are not calibrated MODY risks. WT-like activity: Gersing et al. 2023.",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.075, 1, 0.96], h_pad=2.2)
    fig.savefig(HERE / "GCK_SANITY.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    ax = axes[0]
    subsets = [
        d.loc[mask],
        d.loc[mask & d.origin.eq("population_only")],
        d.loc[mask & ~d.origin.eq("population_only")],
    ]
    labels = [
        f"All supported\n(n={len(subsets[0])})",
        f"Population-only targets\n(n={len(subsets[1])})",
        f"Literature-linked targets\n(n={len(subsets[2])})",
    ]
    prior = np.array([x.prior_additive_component.mean() for x in subsets])
    counts = np.array([x.affected_additive_component.mean() for x in subsets])
    ax.bar(labels, prior, color="#b9cde0", label="Explicit prior term")
    ax.bar(labels, counts, bottom=prior, color="#1865a2", label="Affected-count term")
    for i, (pr, co) in enumerate(zip(prior, counts)):
        ax.text(i, pr / 2, f"{pr:.1%}", ha="center", va="center")
        ax.text(i, pr + co / 2, f"{co:.1%}", ha="center", va="center", color="white")
        ax.text(i, pr + co + 0.015, f"{pr + co:.1%}", ha="center")
    ax.set(
        ylim=(0, 0.65),
        ylabel="Mean neighborhood score",
        title="A  Exact additive decomposition of the frozen score",
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    ax = axes[1]
    examples = d.set_index("protein_key").loc[
        ["W99C", "W99L", "V389D", "V389F", "D217N", "E279Q"]
    ]
    x = np.arange(len(examples))
    ax.bar(
        x - 0.18,
        examples.density_reconstructed,
        width=0.36,
        color="#1865a2",
        label="Frozen primary",
    )
    ax.bar(
        x + 0.18,
        examples.density_HH_units_quarantined_fixed_prior,
        width=0.36,
        color="#d28227",
        label="Quarantine W99R / V389L donors",
    )
    for i, (_, row) in enumerate(examples.iterrows()):
        change = 100 * row.HH_quarantine_density_change
        change_label = f"{change:.2f} pp" if abs(change) >= 0.005 else "<0.01 pp"
        ax.text(
            i,
            row.density_reconstructed + 0.015,
            change_label,
            ha="center",
            fontsize=8,
        )
    ax.set(
        xticks=x,
        xticklabels=examples.index,
        ylim=(0, 0.65),
        ylabel="Neighborhood score",
        title="B  Endpoint errors matter locally, especially at the same residue",
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    fig.suptitle(
        "GCK: prior contribution and a conservative endpoint sensitivity",
        fontsize=15,
        y=0.98,
    )
    fig.text(
        0.02,
        0.025,
        "A: D = Σ w [μκ/(κ+n) + A/(κ+n)]. Terms add exactly; this is not a causal or variance attribution.\n"
        "B: Remove two HH-linked donor units; keep the frozen prior and all other counts. Their archived A=14 includes 13 confirmed HH observations + 1 unassessed carrier.\n"
        "Remaining observations are still a pooled, incompletely adjudicated endpoint set. No observations are converted to unaffected; same-residue variants remain eligible.",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.16, 1, 0.94], w_pad=2)
    fig.savefig(HERE / "GCK_PRIOR_ENDPOINT_SENSITIVITY.png")
    plt.close(fig)


def main(repeats):
    assert repeats > 0
    endpoint_path = HERE / "endpoint_review.csv"
    frozen_endpoint_path = HERE / "frozen_endpoint_observations.csv"
    endpoint = pd.read_csv(endpoint_path).merge(
        pd.read_csv(frozen_endpoint_path),
        left_on=["protein_key", "pmid"],
        right_on=["key", "pmid"],
        validate="one_to_one",
    )
    assert len(endpoint) == 6
    np.testing.assert_array_equal(endpoint.frozen_affected, endpoint.affected)
    np.testing.assert_array_equal(endpoint.frozen_unaffected, endpoint.unaffected)
    hh_endpoint = endpoint.loc[endpoint.protein_key.isin(["W99R", "V389L"])]
    np.testing.assert_array_equal(
        hh_endpoint.confirmed_HH_observations
        + hh_endpoint.phenotype_unknown_observations,
        hh_endpoint.frozen_affected,
    )
    assert hh_endpoint.confirmed_HH_observations.sum() == 13
    assert hh_endpoint.phenotype_unknown_observations.sum() == 1
    input_path = BASE / "GCK_structural_input_variants.csv.gz"
    density_path = BASE / "GCK_primary_variant_density.csv"
    d = pd.read_csv(input_path).set_index("variant_id", drop=False)
    published = pd.read_csv(density_path).set_index("variant_id").loc[d.index]
    paths = sorted(BASE.glob("GCK_primary_normalized_weights.part*.csv.gz"))
    weights = pd.concat([pd.read_csv(p, index_col=0) for p in paths]).loc[
        d.index, d.index
    ]
    w = weights.to_numpy()
    mask = published.density.notna().to_numpy()
    assert len(d) == 634 and mask.sum() == 614 and d.index.is_unique
    assert np.all(w >= 0) and np.max(np.abs(np.diag(w))) == 0
    np.testing.assert_allclose(w.sum(axis=1), mask.astype(float), atol=1e-12)
    assert not w[:, ~mask].any()
    a, n = d.affected.to_numpy(), d.n.to_numpy()
    np.testing.assert_allclose(n, a + d.unaffected)
    alpha, beta, mu, var, strength = moments(a, n, "historical")
    np.testing.assert_allclose(d.alpha_empirical, alpha, atol=1e-12)
    np.testing.assert_allclose(d.beta_empirical, beta, atol=1e-12)
    p = (alpha + a) / (strength + n)
    np.testing.assert_allclose(p, d.posterior_mean, atol=1e-12)
    reconstructed = w @ p
    np.testing.assert_allclose(reconstructed[mask], published.density[mask], atol=1e-12)
    pos = d.canonical_pos.to_numpy()
    dist = 3.8 * np.sqrt(abs(pos[:, None] - pos[None, :]))
    s = 2 / (1 + np.exp(np.log(3) * dist / 3))
    np.fill_diagonal(s, 0)
    s /= s.sum(axis=1, keepdims=True)
    d["observed_fraction"] = a / n
    d["own_prior_retention"] = strength / (strength + n)
    d["historical_fit_weight"] = 1 - 1 / (n + 0.01)
    d["n_bin"] = pd.cut(
        n, bins=[0, 1, 2, 5, 10, np.inf], labels=["1", "2", "3-5", "6-10", ">10"]
    )
    diagnostic = {
        "density_reconstructed": reconstructed,
        "count_only_density": w @ (a / n),
        "prior_mixture_weight": w @ (strength / (strength + n)),
        "prior_additive_component": mu * (w @ (strength / (strength + n))),
        "affected_additive_component": w @ (a / (strength + n)),
        "n_weighted_posterior_density": np.divide(w @ (n * p), w @ n, where=mask),
        "n_weighted_count_density": np.divide(w @ a, w @ n, where=mask),
        "population_only_donor_share": w @ d.origin.eq("population_only").astype(float),
        "singleton_donor_share": w @ (n == 1).astype(float),
    }
    for col, vals in diagnostic.items():
        d[col] = np.where(mask, vals, np.nan)
    d["prior_additive_share_of_score"] = (
        d.prior_additive_component / d.density_reconstructed
    )
    d["kish_donor_n"] = published.kish_donor_n
    np.testing.assert_allclose(
        (d.prior_additive_component + d.affected_additive_component)[mask],
        reconstructed[mask],
    )
    assert np.max(abs((w @ np.full(len(d), mu))[mask] - mu)) < 1e-12

    # All five retained W99R/V389L source rows concern HH; one V389L carrier's
    # biochemical status is unknown. Remove whole endpoint-incompatible units,
    # not unknown->unaffected, and do not relabel the remaining pool MODY-only.
    quarantine = d.protein_key.isin(["W99R", "V389L"]).to_numpy()
    assert quarantine.sum() == 2 and a[quarantine].sum() == 14
    qw = w.copy()
    qw[:, quarantine] = 0
    rowsum = qw.sum(axis=1)
    qw[mask] /= rowsum[mask, None]
    d["density_HH_units_quarantined_fixed_prior"] = np.where(mask, qw @ p, np.nan)
    d["quarantined_HH_donor_weight"] = np.where(
        mask, w[:, quarantine].sum(axis=1), np.nan
    )
    d["HH_quarantine_density_change"] = (
        d.density_HH_units_quarantined_fixed_prior - d.density_reconstructed
    )

    save(d.reset_index(drop=True), "variant_diagnostics.csv.gz")
    cols = [
        "variant_id",
        "protein_key",
        "origin",
        "affected",
        "unaffected_literature",
        "gnomad_carriers",
        "n",
        "posterior_mean",
        "density_reconstructed",
        "prior_mixture_weight",
        "density_HH_units_quarantined_fixed_prior",
    ]
    examples = d.loc[
        d.protein_key.isin(
            [
                "D217N",
                "E279Q",
                "W99R",
                "V389L",
                "G68D",
                "V455E",
                "G261R",
                "R191W",
                "F150S",
            ]
        ),
        cols,
    ]
    save(examples.sort_values("protein_key"), "variant_examples.csv")
    groups = []
    for field in ["origin", "n_bin"]:
        for group, x in d.groupby(field, observed=True):
            groups.append(
                dict(
                    grouping=field,
                    group=str(group),
                    units=len(x),
                    affected=x.affected.sum(),
                    unaffected_literature=x.unaffected_literature.sum(),
                    gnomad_carriers=x.gnomad_carriers.sum(),
                    n=x.n.sum(),
                    singleton_units=x.n.eq(1).sum(),
                    fit_weight=x.historical_fit_weight.sum(),
                    fit_weight_share=x.historical_fit_weight.sum()
                    / d.historical_fit_weight.sum(),
                    raw_fraction_mean=x.observed_fraction.mean(),
                    posterior_mean=x.posterior_mean.mean(),
                )
            )
    save(pd.DataFrame(groups), "count_and_fit_weight_groups.csv")
    distributions = [
        dict(quantity=c, **describe(d.loc[mask, c]))
        for c in [
            "posterior_mean",
            "observed_fraction",
            *diagnostic,
            "prior_additive_share_of_score",
            "kish_donor_n",
            "density_HH_units_quarantined_fixed_prior",
            "HH_quarantine_density_change",
        ]
    ]
    save(pd.DataFrame(distributions), "supported_distributions.csv")
    sensitivity = []
    for name in ["historical", "normalized_weighted_mse", "equal_variant"]:
        al, be, m, v, k = moments(a, n, name)
        post = (al + a) / (k + n)
        density = (w @ post)[mask]
        sensitivity.append(
            dict(
                diagnostic=name,
                alpha=al,
                beta=be,
                mean=m,
                variance=v,
                strength=k,
                singleton_unaffected=al / (k + 1),
                singleton_prior_retention=k / (k + 1),
                posterior_mean=post.mean(),
                density_mean=density.mean(),
                density_p05=np.quantile(density, 0.05),
                density_p95=np.quantile(density, 0.95),
            )
        )
    for k in [0.1, 1.0, 10.0]:
        al, be = mu * k, (1 - mu) * k
        post = (al + a) / (k + n)
        density = (w @ post)[mask]
        sensitivity.append(
            dict(
                diagnostic=f"fixed_mean_strength_{k:g}",
                alpha=al,
                beta=be,
                mean=mu,
                variance=mu * (1 - mu) / (k + 1),
                strength=k,
                singleton_unaffected=al / (k + 1),
                singleton_prior_retention=k / (k + 1),
                posterior_mean=post.mean(),
                density_mean=density.mean(),
                density_p05=np.quantile(density, 0.05),
                density_p95=np.quantile(density, 0.95),
            )
        )
    save(pd.DataFrame(sensitivity), "prior_sensitivities_DIAGNOSTIC.csv")
    metrics_path = BASE / "GCK_variant_loo_metrics.csv"
    metrics = pd.read_csv(metrics_path)
    save(
        metrics.loc[
            metrics.cohort.eq("all_common") & metrics.support_stratum.eq("all")
        ],
        "frozen_same_row_loo_metrics.csv",
    )
    summary = permutations(d, w, s, mask, mu, repeats)
    checks = dict(
        input_variants=len(d),
        supported=int(mask.sum()),
        unsupported=int((~mask).sum()),
        prior_alpha=alpha,
        prior_beta=beta,
        prior_mean=mu,
        prior_strength=strength,
        A_total=float(a.sum()),
        U_total=float(d.unaffected.sum()),
        n_total=float(n.sum()),
        pooled_A_over_N=float(a.sum() / n.sum()),
        equal_variant_raw_fraction=float((a / n).mean()),
        equal_variant_posterior=float(p.mean()),
        singleton_units=int((n == 1).sum()),
        population_only_singletons=int(
            (d.origin.eq("population_only") & d.n.eq(1)).sum()
        ),
        population_only_units=int(d.origin.eq("population_only").sum()),
        maximum_density_reconstruction_error=float(
            max(abs(reconstructed[mask] - published.density[mask]))
        ),
        maximum_self_weight=float(np.diag(w).max()),
        permutation_seed=SEED,
        permutation_repeats_per_scheme=repeats,
        fixed_prior_null_density=mu,
        quarantine_units=int(quarantine.sum()),
        quarantine_frozen_A=int(a[quarantine].sum()),
        quarantine_source_confirmed_HH_observations=13,
        quarantine_unassessed_carriers=1,
        quarantine_max_absolute_change=float(abs(d.HH_quarantine_density_change).max()),
        source_hashes={
            str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                input_path,
                density_path,
                metrics_path,
                *paths,
                frozen_endpoint_path,
                endpoint_path,
            ]
        },
    )
    (HERE / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    figure(d, mask, mu, summary)
    print(
        json.dumps({k: v for k, v in checks.items() if k != "source_hashes"}, indent=2)
    )
    print(
        summary.query(
            "outcome == 'posterior' and metric == 'gain_vs_fixed_prior'"
        ).to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=2000)
    args = parser.parse_args()
    main(args.permutations)
