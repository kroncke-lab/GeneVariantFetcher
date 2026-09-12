"""Fit separate gene-by-type empirical priors from the frozen observed union.

Missense and nonsense are the requested types. Nonsense includes equivalent
nonsense/stop_gained labels only. Other consequences never influence these fits.
The historical mean/MSE method, all observed counts, and identity policy stay
unchanged, making this an isolated correction of the prior's variant class.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "population_inclusive_penetrance_20260912"
spec = importlib.util.spec_from_file_location(
    "frozen_union_methods", PREVIOUS / "rebuild_union.py"
)
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
GENES = legacy.GENES
TYPES = {"missense": {"missense"}, "nonsense": {"nonsense", "stop_gained"}}


def select_type(units, variant_type):
    if variant_type not in TYPES:
        raise ValueError(f"Unsupported prior type: {variant_type}")
    selected = units.loc[
        units.vclass.isin(TYPES[variant_type]) & units.canonical_wt_status.eq("match")
    ].copy()
    selected["variant_type"] = variant_type
    if selected.unit_id.duplicated().any():
        raise ValueError("Nonunique type-matched identities")
    return selected


def equal_variant_moments(frame):
    """Explicit diagnostic of equal-weight raw fractions, not the primary fit."""
    y = frame.affected.to_numpy() / frame.n.to_numpy()
    mean, variance = float(y.mean()), float(y.var())
    if not 0 < variance < mean * (1 - mean):
        raise ValueError("Raw-fraction moments do not identify a proper Beta")
    strength = mean * (1 - mean) / variance - 1
    return {
        "mean": mean,
        "variance": variance,
        "alpha_empirical": mean * strength,
        "beta_empirical": (1 - mean) * strength,
        "strength": strength,
    }


def fit_types(units):
    moments, posteriors, counts, examples = [], [], [], []
    for gene in GENES:
        for kind in TYPES:
            frame = select_type(units.loc[units.gene.eq(gene)], kind)
            if len(frame) < 2:
                raise ValueError(f"Insufficient observed units for {gene} {kind}")
            primary = legacy.fit_empirical(frame)
            weights = 1 - 1 / (frame.n + 0.01)
            pop = frame.origin.eq("population_only")
            record = {
                "gene": gene,
                "variant_type": kind,
                "variants": len(frame),
                "population_only": int(pop.sum()),
                "affected": int(frame.affected.sum()),
                "unaffected_literature": int(frame.unaffected_literature.sum()),
                "unaffected_gnomad": int(frame.gnomad_carriers.sum()),
            }
            for scope, params in [
                (f"canonical_{kind}", primary),
                (f"canonical_{kind}_equal_raw_fraction", equal_variant_moments(frame)),
                (
                    f"canonical_{kind}_normalized_weighted_mse",
                    legacy.fit_empirical(frame, "normalized"),
                ),
            ]:
                moments.append(
                    {
                        **record,
                        "scope": scope,
                        **params,
                        "unaffected_singleton_posterior": params["alpha_empirical"]
                        / (params["strength"] + 1),
                        "singleton_prior_retention": params["strength"]
                        / (params["strength"] + 1),
                    }
                )
            post = legacy.posterior_table(frame, primary)
            post["class_scope"] = f"canonical_{kind}"
            posteriors.append(post)
            counts.append(
                {
                    **record,
                    "with_population_members": int(frame.gnomad_carriers.gt(0).sum()),
                    "population_only_pct": 100 * pop.mean(),
                    "n1_units": int(frame.n.eq(1).sum()),
                    "n1_pct": 100 * frame.n.eq(1).mean(),
                    "population_only_n1_units": int((pop & frame.n.eq(1)).sum()),
                    "population_only_n1_pct": 100 * frame.loc[pop, "n"].eq(1).mean(),
                    "affected_singletons": int(
                        (frame.n.eq(1) & frame.affected.eq(1)).sum()
                    ),
                    "unaffected_singletons": int(
                        (frame.n.eq(1) & frame.affected.eq(0)).sum()
                    ),
                    "n1_weight_share_pct": 100
                    * weights.loc[frame.n.eq(1)].sum()
                    / weights.sum(),
                    "total_prior_fit_weight": weights.sum(),
                    "posterior_median": post.posterior_mean.median(),
                    "posterior_below_0_10_pct": 100
                    * post.posterior_mean.lt(0.1).mean(),
                }
            )
            for affected, unaffected in [
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 5),
                (0, 10),
                (0, 100),
                (1, 0),
                (1, 1),
            ]:
                examples.append(
                    {
                        "gene": gene,
                        "variant_type": kind,
                        "affected": affected,
                        "unaffected": unaffected,
                        "prior_mean": primary["mean"],
                        "prior_strength": primary["strength"],
                        "posterior_mean": (primary["alpha_empirical"] + affected)
                        / (primary["strength"] + affected + unaffected),
                        "example_only_not_added_to_fit": True,
                    }
                )
    return (
        pd.DataFrame(moments),
        pd.concat(posteriors, ignore_index=True),
        pd.DataFrame(counts),
        pd.DataFrame(examples),
    )


def make_plots(priors, posteriors, examples, out):
    os.environ.setdefault("MPLCONFIGDIR", str(HERE.parents[2] / "tmp/matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(2, 5, figsize=(16, 7.4), constrained_layout=True)
    for row, kind in enumerate(TYPES):
        for col, gene in enumerate(GENES):
            ax = axes[row, col]
            d = posteriors.loc[
                posteriors.gene.eq(gene) & posteriors.variant_type.eq(kind)
            ]
            p = priors.loc[
                priors.gene.eq(gene) & priors.scope.eq("canonical_" + kind)
            ].iloc[0]
            ax.hist(
                [
                    d.loc[d.origin.eq("population_only"), "posterior_mean"],
                    d.loc[d.origin.ne("population_only"), "posterior_mean"],
                ],
                bins=np.linspace(0, 1, 21),
                stacked=True,
                color=["#78a6b8", "#bd724b"],
                label=["Population only", "Contains literature"],
            )
            ax.axvline(
                p["mean"], color="#152f4a", ls="--", lw=1.6, label="Type-specific prior"
            )
            ax.set(
                title=f"{gene} · {kind}\n{len(d):,} units; prior {p['mean']:.1%}",
                xlabel="Empirical posterior mean",
                yscale="log",
                ylim=(0.8, None),
                xlim=(0, 1),
            )
            if col == 0:
                ax.set_ylabel("Variant units (log scale)")
    axes[0, 4].legend(fontsize=8, loc="upper right")
    fig.suptitle(
        "Separate empirical priors for missense and nonsense · gnomAD carriers assumed unaffected",
        fontsize=14,
    )
    fig.savefig(out / "TYPE_MATCHED_POSTERIORS.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for ax, kind in zip(axes, TYPES):
        for gene in GENES:
            d = examples.loc[
                examples.gene.eq(gene)
                & examples.variant_type.eq(kind)
                & examples.affected.eq(0)
                & examples.unaffected.le(10)
            ]
            ax.plot(d.unaffected, d.posterior_mean, marker="o", ms=4, label=gene)
        ax.set(
            title=kind.capitalize(),
            xlabel="Unaffected carrier observations (A = 0)",
            ylabel="Posterior mean",
            ylim=(0, 1),
            xticks=[0, 1, 2, 5, 10],
        )
    axes[1].legend(fontsize=8)
    fig.suptitle(
        "One unaffected carrier updates the matching prior; it does not set risk to zero",
        fontsize=14,
    )
    fig.savefig(out / "UNAFFECTED_SINGLETON_UPDATE.png", dpi=160)
    plt.close(fig)


def main():
    paths = sorted((PREVIOUS / "analysis/union_counts").glob("*.csv.gz"))
    previous_manifest = json.loads((PREVIOUS / "artifact_manifest.json").read_text())
    for path in paths:
        expected = previous_manifest["files"][str(path.relative_to(PREVIOUS))]["sha256"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Frozen union changed: {path}")
    units = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    priors, posterior, counts, examples = fit_types(units)
    out = HERE / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    legacy.save(priors, out / "empirical_prior_comparison.csv")
    legacy.save_shards(posterior, out / "empirical_posteriors", "")
    legacy.save(counts, out / "type_count_summary.csv")
    legacy.save(examples, out / "singleton_update_examples.csv")
    catalog = (
        units.groupby(["gene", "vclass", "canonical_wt_status"], dropna=False)
        .agg(
            units=("unit_id", "size"),
            affected=("affected", "sum"),
            unaffected=("unaffected", "sum"),
        )
        .reset_index()
    )
    legacy.save(catalog, out / "source_type_inventory.csv")
    old = pd.concat(
        [
            pd.read_csv(p)
            for p in (PREVIOUS / "analysis/empirical_posteriors").glob("*.csv.gz")
        ],
        ignore_index=True,
    )
    comparison = posterior[["gene", "variant_type", "unit_id", "posterior_mean"]].merge(
        old[["unit_id", "posterior_mean", "alpha_empirical", "beta_empirical"]],
        on="unit_id",
        suffixes=("_class", "_full_locus"),
        validate="one_to_one",
    )
    comparison["full_locus_prior_mean"] = comparison.alpha_empirical / (
        comparison.alpha_empirical + comparison.beta_empirical
    )
    legacy.save_shards(comparison, out / "same_variant_comparison", "")
    # Class selection cannot change any observed count or allele membership.
    assert posterior.unit_id.is_unique
    aligned = units.set_index("unit_id").loc[posterior.unit_id]
    np.testing.assert_array_equal(
        aligned.member_alleles.fillna("").to_numpy(),
        posterior.member_alleles.fillna("").to_numpy(),
    )
    for column in [
        "affected",
        "unaffected_literature",
        "gnomad_carriers",
        "unaffected",
        "n",
    ]:
        np.testing.assert_array_equal(
            aligned[column].to_numpy(), posterior[column].to_numpy()
        )
    all_sources = [
        *paths,
        PREVIOUS / "rebuild_union.py",
        PREVIOUS / "artifact_manifest.json",
        *sorted((PREVIOUS / "analysis/empirical_posteriors").glob("*.csv.gz")),
    ]
    checks = {
        "source_union_units": len(units),
        "included_type_units": len(posterior),
        "by_type": posterior.variant_type.value_counts().to_dict(),
        "prior_count": 10,
        "separate_prior_per_gene_and_type": True,
        "nonsense_labels": ["nonsense", "stop_gained"],
        "excluded_from_both_priors": [
            "synonymous",
            "frameshift",
            "splice",
            "noncoding",
            "stop_lost",
            "start_lost",
            "unresolved_canonical_annotation",
        ],
        "canonical_WT_match_required": True,
        "observed_counts_and_membership_unchanged": True,
        "gnomad_all_assumed_unaffected": True,
        "historical_weight_formula_unchanged": "1 - 1/(n + 0.01)",
        "historical_variance_denominator": "number of variant units within gene and type",
        "singleton_weight": 1 - 1 / 1.01,
        "no_hypothetical_A0_U0_units_added": True,
        "input_sha256": {
            str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in all_sources
        },
    }
    (out / "fit_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    make_plots(priors, posterior, examples, out)
    print(
        priors.loc[
            priors.scope.isin(["canonical_missense", "canonical_nonsense"]),
            [
                "gene",
                "variant_type",
                "variants",
                "mean",
                "alpha_empirical",
                "beta_empirical",
                "strength",
                "unaffected_singleton_posterior",
            ],
        ].to_string(index=False),
        flush=True,
    )


if __name__ == "__main__":
    main()
