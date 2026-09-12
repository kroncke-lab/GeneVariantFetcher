"""Compile five-gene scientific figures and source-stratified LOO readout."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from validate_outputs import read_parts


HERE = Path(__file__).resolve().parent
CLASS = HERE.parent / "class_matched_penetrance_20260912"
os.environ.setdefault("MPLCONFIGDIR", str(HERE.parents[2] / "tmp/matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

spec = importlib.util.spec_from_file_location(
    "frozen_metrics", HERE.parent / "gck_structural_pilot_20260912/pilot_statistics.py"
)
metrics_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(metrics_module)
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
LENGTHS = {"HNF1A": 631, "GCK": 465, "LDLR": 860, "BRCA2": 3418, "KCNQ1": 676}
SOURCE_COLORS = {"structured": "#176c91", "polymer": "#c67b20", "mixed": "#816495"}


def main():
    out = HERE / "analysis"
    priors = (
        pd.read_csv(CLASS / "analysis/empirical_prior_comparison.csv")
        .query("scope == 'canonical_missense'")
        .set_index("gene")
    )
    rows, all_primary, all_metrics, stratified = [], {}, [], []
    for gene in GENES:
        if gene == "GCK":
            primary = pd.read_csv(CLASS / "structural/GCK_primary_variant_density.csv")
            primary["canonical_pos"] = primary.aa_pos
            predictions = (
                pd.read_csv(CLASS / "structural/GCK_variant_loo_predictions.csv.gz")
                .query("cohort == 'all_common'")
                .copy()
            )
            predictions["cohort"] = "am_common"
            total = 634
            frame = "1V4S"
        else:
            primary = read_parts(out / gene, "primary_density")
            predictions = read_parts(out / gene, "loo_predictions")
            total = len(primary)
            frame = str(primary.frame.iloc[0])
        supported = primary.loc[primary.density.notna()]
        all_primary[gene] = primary
        rows.append(
            {
                "gene": gene,
                "missense_variants": total,
                "supported": len(supported),
                "structured_only": int(supported.density_source.eq("structured").sum()),
                "polymer_only": int(supported.density_source.eq("polymer").sum()),
                "mixed": int(supported.density_source.eq("mixed").sum()),
                "unavailable": total - len(supported),
                "prior_mean": priors.loc[gene, "mean"],
                "median_own_posterior_supported": supported.posterior_mean.median(),
                "median_density": supported.density.median(),
                "median_kish_donor_n": supported.kish_donor_n.median(),
                "median_weight_share_beyond_20": supported.weight_share_beyond_20.median(),
                "primary_frame": frame,
            }
        )
        predictions = predictions.merge(
            primary[["variant_id", "density_source"]],
            on="variant_id",
            validate="many_to_one",
        )
        for cohort, group in predictions.groupby("cohort"):
            table = metrics_module.comparison_metrics(
                group.assign(support_stratum="all"), priors.loc[gene, "strength"]
            )
            table.loc[
                table.model.eq("intercept_only"), "spearman_observed_fraction"
            ] = np.nan
            table["gene"], table["cohort"] = gene, cohort
            all_metrics.append(table)
            for source, segment in group.groupby("density_source"):
                table = metrics_module.comparison_metrics(
                    segment.assign(support_stratum=source), priors.loc[gene, "strength"]
                )
                table.loc[
                    table.model.eq("intercept_only"), "spearman_observed_fraction"
                ] = np.nan
                table["gene"], table["cohort"] = gene, cohort
                stratified.append(table)
    overview = pd.DataFrame(rows)
    overview.to_csv(out / "five_gene_summary.csv", index=False, lineterminator="\n")
    pd.concat(all_metrics).to_csv(
        out / "five_gene_loo_metrics.csv", index=False, lineterminator="\n"
    )
    pd.concat(stratified).to_csv(
        out / "geometry_stratified_loo_metrics.csv", index=False, lineterminator="\n"
    )
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(1, 5, figsize=(17, 4.5), constrained_layout=True)
    for ax, gene in zip(axes, GENES):
        data = all_primary[gene]
        supported = data.loc[data.density.notna()]
        bins = np.linspace(0, 1, 21)
        ax.hist(
            supported.posterior_mean,
            bins=bins,
            color="#91b1c3",
            alpha=0.65,
            label="Own posterior",
        )
        ax.hist(
            supported.density,
            bins=bins,
            histtype="step",
            linewidth=2,
            color="#ad5426",
            label="Variant-excluded density",
        )
        ax.axvline(
            priors.loc[gene, "mean"],
            color="#263b49",
            linestyle="--",
            linewidth=1.3,
            label="Missense prior",
        )
        ax.set(
            xlim=(0, 1),
            title=f"{gene}\n{len(supported):,} / {len(data):,} supported",
            xlabel="Probability / density",
        )
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.set_xticks([0, 0.5, 1])
    axes[0].set_ylabel("Variants (separate y scales)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle(
        "Missense priors and 3D / polymer neighborhoods · same supported variants in each panel",
        fontsize=15,
    )
    fig.savefig(out / "MISSENSE_STRUCTURE_DISTRIBUTIONS.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(5, 1, figsize=(14, 12), constrained_layout=True)
    for ax, gene in zip(axes, GENES):
        data = all_primary[gene]
        supported = data.loc[data.density.notna()]
        for source, group in supported.groupby("density_source"):
            color = SOURCE_COLORS.get(source, "#666666")
            ax.vlines(
                group.canonical_pos,
                group.density_lower_95,
                group.density_upper_95,
                color=color,
                alpha=0.13,
                linewidth=0.6,
            )
            ax.scatter(
                group.canonical_pos,
                group.density,
                color=color,
                s=8,
                alpha=0.6,
                label=source,
            )
        unavailable = data.loc[data.density.isna()]
        ax.scatter(
            unavailable.canonical_pos,
            np.full(len(unavailable), -0.06),
            color="#b3b3b3",
            s=7,
            marker="|",
            label="Unavailable",
        )
        ax.axhline(
            priors.loc[gene, "mean"], color="#333333", linewidth=1, linestyle="--"
        )
        ax.set(
            xlim=(1, LENGTHS[gene]),
            ylim=(-0.1, 1),
            ylabel=gene,
            xlabel="Canonical residue",
        )
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.set_yticks([0, 0.5, 1])
        ax.grid(axis="y", alpha=0.12)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle(
        "Variant-excluded density along each protein\nThin bars: conditional 95% intervals; dashed line: gene-specific missense prior",
        fontsize=15,
    )
    fig.savefig(out / "MISSENSE_STRUCTURE_RESIDUE_MAP.png", dpi=150)
    plt.close(fig)
    receipt = {
        "schemaVersion": 1,
        "items": [
            {
                "id": "missense-structure",
                "title": "Missense structural analysis",
                "queries": [
                    {
                        "id": "five-gene-summary",
                        "source": {
                            "label": "Frozen missense counts and verified structural geometry",
                            "files": [
                                {"label": "five_gene_summary.csv"},
                                {"label": "five_gene_loo_metrics.csv"},
                            ],
                            "metricDefinitions": [
                                {
                                    "label": "Density",
                                    "definition": "Equal-variant sigmoid-weighted empirical posterior mean, excluding the target identity across all assembly copies.",
                                }
                            ],
                            "filters": [
                                "Gene-specific missense prior",
                                "gnomAD carriers assumed unaffected",
                                "Other variants at target residue retained",
                            ],
                            "caveats": [
                                "BRCA2 covers only experimentally resolved domains",
                                "HNF1A has a partial DNA-bound dimer; KCNQ1 is an engineered partial channel tetramer",
                                "LDLR uses a predicted monomer with uncertain interdomain packing",
                                "Polymer-only estimates are distinct from resolved 3D contacts",
                                "Priors fixed in internal variant-only LOO; not independent clinical validation",
                            ],
                        },
                        "columns": [
                            "Gene",
                            "Missense variants",
                            "Supported",
                            "Structured only",
                            "Polymer only",
                            "Mixed",
                            "Median density",
                        ],
                        "rows": [
                            {
                                "Gene": r["gene"],
                                "Missense variants": r["missense_variants"],
                                "Supported": r["supported"],
                                "Structured only": r["structured_only"],
                                "Polymer only": r["polymer_only"],
                                "Mixed": r["mixed"],
                                "Median density": f"{100 * r['median_density']:.2f}%",
                            }
                            for r in rows
                        ],
                        "preview": {
                            "kind": "aggregate",
                            "note": "All five genes; GCK reuses the frozen class-specific rerun",
                            "totalRows": 5,
                        },
                        "methods": [
                            {
                                "language": "calculation",
                                "code": "posterior = Beta(alpha_missense + affected, beta_missense + literature_unaffected + gnomAD_carriers); kernel(d) = 2 / (1 + exp(log(3) * d / 3)); density_i = weighted mean of posterior_j for j != i",
                            }
                        ],
                    }
                ],
            }
        ],
    }
    common_metrics = pd.concat(all_metrics).query(
        "cohort == 'am_common' and model in ['am_fit', 'am_plus_density', 'am_plus_sequence']"
    )
    receipt["items"][0]["queries"].append(
        {
            "id": "paired-loo-comparison",
            "source": {
                "label": "Paired variant-only outer-LOO predictions",
                "files": [
                    {"label": "five_gene_loo_metrics.csv"},
                    {"label": "run_structure.py"},
                ],
                "filters": [
                    "Identical available target variants within each gene",
                    "Gene-specific missense priors held fixed",
                    "Target identity removed globally from every training neighborhood",
                ],
                "caveats": [
                    "BRCA2 AM comparison contains only 44 literature-only archive fallbacks",
                    "Shared empirical priors and selected literature counts prevent claims of independent clinical validation",
                ],
            },
            "columns": [
                "gene",
                "model",
                "n_variants",
                "mae_empirical_posterior",
                "mse_observed_fraction",
            ],
            "rows": common_metrics[
                [
                    "gene",
                    "model",
                    "n_variants",
                    "mae_empirical_posterior",
                    "mse_observed_fraction",
                ]
            ].to_dict("records"),
            "preview": {
                "kind": "aggregate",
                "note": "Same-target model comparisons; lower error is better",
                "totalRows": len(common_metrics),
            },
            "methods": [
                {
                    "language": "calculation",
                    "code": "MAE = mean(abs(LOO_prediction - empirical_posterior_mean)); each variant has equal weight.",
                }
            ],
        }
    )
    (out / "sources_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(overview.to_string(index=False))


if __name__ == "__main__":
    main()
