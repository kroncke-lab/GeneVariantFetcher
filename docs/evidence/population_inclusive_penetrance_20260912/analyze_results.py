"""Summarize the frozen population-inclusive prior/count/predictor comparison."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, spearmanr

from rebuild_union import GENES, save, save_shards

HERE = Path(__file__).resolve().parent
OUT = HERE / "analysis"
METRICS = {
    "alphamissense": "AlphaMissense ↑ pathogenicity",
    "gpn_star_m447_llr_calibrated": "GPN-Star M447 signed LLR",
    "alphagenome_avi": "AlphaGenome AVI ↑ impact",
}


def resolve_unit_scores(units, snapshot):
    """Only compatible exact genomic members; a conflict cannot become a score."""
    by_id = {(r.gene, r.variant_id): r for r in snapshot.itertuples()}
    rows = []
    for row in units.itertuples():
        members = [] if pd.isna(row.member_alleles) else row.member_alleles.split(";")
        source = [
            by_id[(row.gene, allele)]
            for allele in members
            if (row.gene, allele) in by_id
        ]
        output = {"gene": row.gene, "unit_id": row.unit_id}
        for metric in METRICS:
            conflicts = any(
                "conflict" in str(getattr(r, metric + "_status")) for r in source
            )
            values = {
                (str(getattr(r, metric + "_version")), float(getattr(r, metric)))
                for r in source
                if pd.notna(getattr(r, metric))
            }
            output[metric] = (
                next(iter(values))[1] if len(values) == 1 and not conflicts else np.nan
            )
            output[metric + "_status"] = (
                "member_conflict"
                if conflicts or len(values) > 1
                else "available"
                if len(values) == 1
                else "no_resolved_member_score"
            )
            output[metric + "_scored_members"] = sum(
                pd.notna(getattr(r, metric)) for r in source
            )
        rows.append(output)
    return units.merge(
        pd.DataFrame(rows), on=["gene", "unit_id"], validate="one_to_one"
    )


def main():
    full = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted((OUT / "empirical_posteriors").glob("*.csv.gz"))
        ],
        ignore_index=True,
    )
    prior = pd.read_csv(OUT / "empirical_prior_comparison.csv")
    current = prior.loc[prior.scope.eq("all_observed_full_locus")].set_index("gene")
    old = pd.read_csv(
        HERE.parent
        / "structural_density_plan_20260912/full_universe_empirical_posteriors.csv.gz"
    )
    old_moments = old.groupby("gene")[["alpha_empirical", "beta_empirical"]].first()
    old_moments["mean"] = old_moments.alpha_empirical / old_moments.sum(axis=1)
    summaries, distributions = [], []
    for gene in GENES:
        f = full.loc[full.gene.eq(gene)]
        old_gene = old.loc[old.gene.eq(gene)]
        pair = f.loc[f.literature_key.notna()].merge(
            old_gene[["key", "posterior_empirical_mean"]],
            left_on="literature_key",
            right_on="key",
            validate="one_to_one",
        )
        summaries.append(
            {
                "gene": gene,
                "old_variants": len(old_gene),
                "new_units": len(f),
                "old_prior_alpha": old_moments.loc[gene, "alpha_empirical"],
                "old_prior_beta": old_moments.loc[gene, "beta_empirical"],
                "old_prior_mean": old_moments.loc[gene, "mean"],
                "new_prior_alpha": current.loc[gene, "alpha_empirical"],
                "new_prior_beta": current.loc[gene, "beta_empirical"],
                "new_prior_mean": current.loc[gene, "mean"],
                "new_median_posterior": f.posterior_mean.median(),
                "new_posterior_below_0_01_pct": 100 * f.posterior_mean.lt(0.01).mean(),
                "new_posterior_below_0_10_pct": 100 * f.posterior_mean.lt(0.1).mean(),
                "new_posterior_above_0_50_pct": 100 * f.posterior_mean.ge(0.5).mean(),
                "same_clinical_units": len(pair),
                "same_clinical_old_median": pair.posterior_empirical_mean.median(),
                "same_clinical_new_median": pair.posterior_mean.median(),
                "same_clinical_old_below_0_10_pct": 100
                * pair.posterior_empirical_mean.lt(0.1).mean(),
                "same_clinical_new_below_0_10_pct": 100
                * pair.posterior_mean.lt(0.1).mean(),
            }
        )
        for label, group in [("all", f), *list(f.groupby("origin"))]:
            for field in [
                "affected",
                "unaffected_literature",
                "gnomad_carriers",
                "unaffected",
                "posterior_mean",
            ]:
                x = group[field]
                distributions.append(
                    {
                        "gene": gene,
                        "origin": label,
                        "quantity": field,
                        "units": len(x),
                        "sum": x.sum(),
                        "zero_pct": 100 * x.eq(0).mean(),
                        "below_0_01_pct": 100 * x.lt(0.01).mean(),
                        "below_0_10_pct": 100 * x.lt(0.1).mean(),
                        "median": x.median(),
                        "p90": x.quantile(0.9),
                        "p99": x.quantile(0.99),
                        "max": x.max(),
                    }
                )
    save(pd.DataFrame(summaries), OUT / "before_after_summary.csv")
    save(pd.DataFrame(distributions), OUT / "count_posterior_distributions.csv")
    snapshot = pd.read_csv(HERE / "predictors/population_predictors.csv.gz")
    scored = resolve_unit_scores(full, snapshot)
    save_shards(
        scored.loc[scored[list(METRICS)].notna().any(axis=1)],
        OUT / "resolved_predictor_units",
        "",
    )
    associations = []
    for gene in GENES:
        f = scored.loc[scored.gene.eq(gene)]
        for metric in METRICS:
            for label, group in [
                ("all_available", f),
                ("three_predictor_common_rows", f.dropna(subset=list(METRICS))),
            ]:
                selected = group.dropna(subset=[metric])
                associations.append(
                    {
                        "gene": gene,
                        "predictor": metric,
                        "scope": label,
                        "all_gene_units": len(f),
                        "available_units": len(selected),
                        "population_only_units": int(
                            selected.origin.eq("population_only").sum()
                        ),
                        "posterior_spearman": float(
                            spearmanr(
                                selected[metric], selected.posterior_mean
                            ).statistic
                        )
                        if len(selected) > 2 and selected[metric].nunique() > 1
                        else np.nan,
                        "prior_spearman": np.nan,
                        "prior_reason": "shared empirical prior is constant within gene; no predictor used in fit",
                    }
                )
    save(pd.DataFrame(associations), OUT / "predictor_associations.csv")
    os.environ.setdefault("MPLCONFIGDIR", str(HERE.parents[2] / "tmp/matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    bins = np.linspace(0, 1, 26)
    fig, axes = plt.subplots(1, 5, figsize=(16, 4.2), constrained_layout=True)
    for ax, gene in zip(axes, GENES):
        f = full.loc[full.gene.eq(gene)]
        ax.hist(
            [
                f.loc[f.origin.eq("population_only"), "posterior_mean"],
                f.loc[f.origin.ne("population_only"), "posterior_mean"],
            ],
            bins=bins,
            stacked=True,
            color=["#76a6bc", "#ce784b"],
            label=["Population only", "Contains literature"],
        )
        ax.axvline(
            current.loc[gene, "mean"],
            color="#172f4a",
            ls="--",
            lw=1.5,
            label="New prior mean",
        )
        ax.set(
            title=f"{gene} · {len(f):,} units",
            xlabel="Empirical posterior mean",
            yscale="log",
            ylim=(0.8, None),
            xlim=(0, 1),
        )
    axes[0].set_ylabel("Variant units (log scale)")
    axes[-1].legend(fontsize=8, loc="upper right")
    fig.suptitle(
        "Population-inclusive posteriors · every included gnomAD carrier assumed unaffected",
        fontsize=14,
    )
    fig.savefig(OUT / "PENETRANCE_HISTOGRAMS.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(5, 3, figsize=(13, 14), constrained_layout=True)
    for i, gene in enumerate(GENES):
        f = full.loc[full.gene.eq(gene)]
        o, n = old_moments.loc[gene], current.loc[gene]
        ax = axes[i, 0]
        x = (bins[:-1] + bins[1:]) / 2
        ax.plot(
            x,
            np.diff(beta.cdf(bins, o.alpha_empirical, o.beta_empirical)) * 100,
            color="#b9764b",
            label=f"Old mean {o['mean']:.3f}",
        )
        ax.plot(
            x,
            np.diff(beta.cdf(bins, n.alpha_empirical, n.beta_empirical)) * 100,
            color="#246982",
            label=f"New mean {n['mean']:.3f}",
        )
        ax.set(
            title=f"{gene} · prior Beta probability mass",
            xlabel="Probability",
            ylabel="Mass per 0.04 bin (%)",
            xlim=(0, 1),
        )
        ax.legend(fontsize=8)
        ax = axes[i, 1]
        cbins = np.arange(-0.05, 7, 0.25)
        ax.hist(
            [
                np.log10(1 + f.affected),
                np.log10(1 + f.unaffected_literature),
                np.log10(1 + f.gnomad_carriers),
            ],
            bins=cbins,
            histtype="step",
            linewidth=1.6,
            label=["Affected", "Literature unaffected", "gnomAD unaffected"],
            color=["#b84b37", "#886ba0", "#267d8f"],
        )
        ax.set(
            xlabel="log10(1 + carrier observations)",
            ylabel="Variant units",
            yscale="log",
            title=f"{gene} · observed counts",
        )
        if i == 0:
            ax.legend(fontsize=8)
        ax = axes[i, 2]
        pair = f.loc[f.literature_key.notna()].merge(
            old.loc[old.gene.eq(gene), ["key", "posterior_empirical_mean"]],
            left_on="literature_key",
            right_on="key",
            validate="one_to_one",
        )
        ax.scatter(
            pair.posterior_empirical_mean,
            pair.posterior_mean,
            s=8,
            alpha=0.35,
            color="#276d86",
            rasterized=True,
        )
        ax.plot([0, 1], [0, 1], ":", color="grey")
        ax.set(
            xlabel="Old empirical posterior",
            ylabel="New empirical posterior",
            title=f"{gene} · same {len(pair):,} clinical keys",
            xlim=(0, 1),
            ylim=(0, 1),
        )
    fig.suptitle(
        "What changed: prior distribution, observed counts, and the same clinical variants",
        fontsize=15,
    )
    fig.savefig(OUT / "PRIORS_AND_COUNTS.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(5, 3, figsize=(13, 14), constrained_layout=True)
    for i, gene in enumerate(GENES):
        f = scored.loc[scored.gene.eq(gene)]
        for j, (metric, label) in enumerate(METRICS.items()):
            ax = axes[i, j]
            d = f.dropna(subset=[metric])
            for origin, color in [
                ("population_only", "#76a6bc"),
                ("clinical", "#ce784b"),
            ]:
                group = d.loc[
                    d.origin.eq("population_only")
                    if origin == "population_only"
                    else d.origin.ne("population_only")
                ]
                ax.scatter(
                    group[metric],
                    group.posterior_mean,
                    s=8,
                    alpha=0.3,
                    color=color,
                    rasterized=True,
                    label="Population only"
                    if origin == "population_only"
                    else "Contains literature",
                )
            ax.axhline(
                current.loc[gene, "mean"],
                color="#172f4a",
                ls="--",
                lw=1.2,
                label="Shared prior",
            )
            ax.set(
                xlabel=label,
                ylabel="Empirical posterior mean",
                ylim=(0, 1),
                title=f"{gene} · {len(d):,} resolved units",
            )
            if not len(d):
                if metric == "alphamissense":
                    ax.set_xlim(0, 1)
                ax.text(
                    0.5,
                    0.5,
                    "No resolved scores\n(version conflicts retained as missing)",
                    ha="center",
                    transform=ax.transAxes,
                    fontsize=9,
                )
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Predictors versus count-updated posteriors · none enters the empirical prior",
        fontsize=15,
    )
    fig.savefig(OUT / "POSTERIOR_VS_PREDICTORS.png", dpi=150)
    plt.close(fig)
    (OUT / "analysis_checks.json").write_text(
        json.dumps(
            {
                "union_units": len(full),
                "unique_units": bool(full.unit_id.is_unique),
                "genes": GENES,
                "predictor_join": "exact genomic members only; disagreements and absent scores stay missing; no archive fallback in this five-gene comparison",
                "predictor_inventory_scope": "original canonical gene footprint; region-only extra alleles unscored",
                "prior_has_no_predictor_input": True,
                "figure_count": 3,
            },
            indent=2,
        )
        + "\n"
    )
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
