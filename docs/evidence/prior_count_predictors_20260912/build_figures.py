#!/usr/bin/env python3
"""Plot frozen fitted priors against counts and genomic scores; gnomAD = unaffected."""

from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / "grant_e2e_20260909" / "analysis"
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
BLUE, GOLD, GRAY = "#2875a5", "#c58220", "#8b959e"


def rho(x, y):
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan
    return float(spearmanr(x, y).statistic)


def setup(ncols, title, subtitle):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5})
    fig, axes = plt.subplots(5, ncols, figsize=(13.2, 14.8), squeeze=False)
    fig.subplots_adjust(
        left=0.08, right=0.975, top=0.885, bottom=0.085, hspace=0.47, wspace=0.26
    )
    fig.suptitle(title, x=0.08, y=0.978, ha="left", fontsize=20)
    fig.text(0.08, 0.946, subtitle, fontsize=11.5)
    for i, gene in enumerate(GENES):
        for ax in axes[i]:
            ax.set_axisbelow(True)
            ax.grid(color="#e4e7ea", linewidth=0.55)
            ax.spines[["top", "right"]].set_visible(False)
        axes[i, 0].text(
            -0.16,
            1.12,
            gene,
            transform=axes[i, 0].transAxes,
            fontsize=13.5,
            weight="bold",
        )
    return fig, axes


def load():
    frames = {}
    for gene in GENES:
        d = pd.read_csv(SOURCE / f"{gene}_protocol" / "posterior_variants.csv")
        assert d.key.is_unique
        d["prior"] = d.insilico_prior_mean
        d["A"] = d.affected
        d["U_literature"] = d.literature_n - d.affected
        d["U"] = d.U_literature + d.gnomad_added
        d["fraction"] = d.A / (d.A + d.U)
        assert np.allclose(d.n, d.A + d.U)
        assert (d.U >= 0).all()
        frames[gene] = d
    return frames


def distributions(frames):
    fig, axes = setup(
        2,
        "Prior distributions and affected / unaffected counts",
        "Five archived protocol fits • each variant has equal weight • gnomAD counts assumed unaffected",
    )
    fig.subplots_adjust(top=0.865)
    prior_bins = np.linspace(0, 1, 26)
    edges = np.array([-0.5, 0.5, 1.5, 4.5, 9.5, 49.5, 99.5, 999.5, np.inf])
    labels = ["0", "1", "2–4", "5–9", "10–49", "50–99", "100–999", "1,000+"]
    for i, gene in enumerate(GENES):
        d = frames[gene]
        ax, bx = axes[i]
        weights = np.full(len(d), 100 / len(d))
        ax.hist(
            d.prior,
            bins=prior_bins,
            weights=weights,
            histtype="step",
            color=GOLD,
            linewidth=1.8,
            label="Fitted prior",
        )
        ax.hist(
            d.insilico_post_mean,
            bins=prior_bins,
            weights=weights,
            histtype="step",
            color="#555d65",
            linewidth=1.5,
            linestyle=":",
            label="Final posterior",
        )
        ax.hist(
            d.fraction,
            bins=prior_bins,
            weights=weights,
            histtype="step",
            color=BLUE,
            linewidth=1.5,
            linestyle="--",
            label="Count-only fraction A / (A + U)",
        )
        ax.set(
            xlim=(0, 1),
            ylim=(0, 100),
            xlabel="Probability / affected fraction",
            ylabel="Percent of variants",
        )
        ax.text(
            0.04,
            0.88,
            f"n = {len(d):,}; prior median = {d.prior.median():.3f}",
            transform=ax.transAxes,
            fontsize=9.5,
        )
        x = np.arange(len(labels))
        for shift, col, color, label, hatch in [
            (-0.25, "A", GOLD, "Affected A", None),
            (0, "U_literature", GRAY, "Literature unaffected", "//"),
            (0.25, "U", BLUE, "Total unaffected U", None),
        ]:
            counts, _ = np.histogram(d[col], edges)
            assert counts.sum() == len(d)
            bx.bar(
                x + shift,
                100 * counts / len(d),
                width=0.25,
                color=color,
                label=label,
                hatch=hatch,
                linewidth=0.4,
                edgecolor="white",
            )
        bx.set(
            xticks=x,
            xticklabels=labels,
            ylim=(0, 100),
            xlabel="Number of affected / unaffected observations",
            ylabel="Percent of variants",
        )
        bx.tick_params(axis="x", labelsize=8.5)
    axes[0, 0].set_title(
        "Distribution of prior vs count-only fraction", pad=14, fontsize=12
    )
    axes[0, 1].set_title(
        "Count distributions (literature U is a subset of total U)",
        pad=14,
        fontsize=11.5,
    )
    handles, labels0 = axes[0, 0].get_legend_handles_labels()
    handles1, labels1 = axes[0, 1].get_legend_handles_labels()
    fig.legend(
        handles + handles1,
        labels0 + labels1,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.928),
        ncol=3,
        frameon=False,
        fontsize=10,
    )
    fig.text(
        0.08,
        0.045,
        "A = literature affected; U = literature unaffected + gnomAD. A 1/1 case contributes count-only fraction 1, even with little evidence.",
        fontsize=10,
    )
    fig.text(
        0.08,
        0.025,
        "Prior means are the archived feature-conditioned values (S = 10). These are distributions across variant keys, not across people.",
        fontsize=10,
    )
    fig.savefig(OUT / "PRIOR_COUNT_DISTRIBUTIONS.png", dpi=170, facecolor="white")
    plt.close(fig)


def scatter_groups(ax, d, x, y="prior"):
    groups = [
        (d.vclass.eq("missense"), BLUE, "o", "Missense"),
        (d.is_truncating.astype(bool), GOLD, "^", "Truncating"),
        (
            ~d.vclass.eq("missense") & ~d.is_truncating.astype(bool),
            GRAY,
            "s",
            "Other consequences",
        ),
    ]
    for mask, color, marker, label in groups:
        ax.scatter(
            d.loc[mask, x],
            d.loc[mask, y],
            s=12,
            color=color,
            marker=marker,
            alpha=0.45,
            linewidths=0,
            label=label,
            rasterized=True,
        )
    ax.set(ylim=(-0.025, 1.025), yticks=[0, 0.5, 1])


def vs_counts(frames):
    fig, axes = setup(
        3,
        "Fitted priors versus affected and unaffected evidence",
        "All retained variants • gnomAD counts assumed unaffected • y-axis = prior mean in every panel",
    )
    specs = [
        ("A", "Affected count A"),
        ("U", "Unaffected count U"),
        ("fraction", "Count-only affected fraction A / (A + U)"),
    ]
    maxima = {col: max(d[col].max() for d in frames.values()) for col, _ in specs}
    for i, gene in enumerate(GENES):
        d = frames[gene]
        for j, (col, label) in enumerate(specs):
            ax = axes[i, j]
            scatter_groups(ax, d, col)
            ax.set_xlabel(label, fontsize=9.5)
            if j == 0:
                ax.set_ylabel("Prior mean")
            if col == "fraction":
                ax.plot([0, 1], [0, 1], color="#454b50", linewidth=0.7, linestyle="--")
                ax.set_xlim(-0.025, 1.025)
            else:
                ax.set_xscale("symlog", linthresh=1)
                ax.set_xlim(-0.12, maxima[col] * 1.3)
                ticks = [
                    v
                    for v in [0, 1, 10, 100, 1000, 10000, 100000, 1000000]
                    if v <= maxima[col] * 1.3
                ]
                ax.set_xticks(ticks)
                ax.xaxis.set_major_formatter(
                    FuncFormatter(
                        lambda x, _: (
                            f"{x / 1e6:g}m"
                            if x >= 1e6
                            else f"{x / 1000:g}k"
                            if x >= 1000
                            else f"{x:g}"
                        )
                    )
                )
            r = rho(d[col], d.prior)
            ax.text(
                0.035,
                0.925,
                f"n = {len(d):,}; Spearman ρ = {r:+.2f}",
                transform=ax.transAxes,
                fontsize=9,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.82,
                    "edgecolor": "none",
                    "pad": 1,
                },
            )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.928),
        ncol=3,
        frameon=False,
    )
    fig.text(
        0.08,
        0.045,
        "Count axes retain zero and become logarithmic above 1; each column uses the same scale across genes. Points are not jittered.",
        fontsize=10,
    )
    fig.text(
        0.08,
        0.025,
        "Shared horizontal bands reflect common feature priors; correlation describes the frozen fit and is not held-out prediction performance.",
        fontsize=10,
    )
    fig.savefig(OUT / "PRIOR_VS_COUNTS.png", dpi=170, facecolor="white")
    plt.close(fig)


def vs_predictors(frames):
    scores = pd.read_csv(OUT / "predictor_scores.csv.gz")
    assert not scores.duplicated(["gene", "key"]).any()
    specs = [
        ("alphamissense_frozen", "AlphaMissense (higher = more impact)"),
        ("gpn_star_m447_llr_calibrated", "GPN-Star-M signed LLR (lower = more impact)"),
        ("alphagenome_avi", "AlphaGenome AVI (higher = more impact)"),
    ]
    combined = {}
    for gene, d in frames.items():
        d = d.rename(columns={"alphamissense": "alphamissense_frozen"})
        sc = scores.loc[
            scores.gene.eq(gene),
            [
                "key",
                "mapped_allele_count",
                "gpn_star_m447_llr_calibrated",
                "alphagenome_avi",
            ],
        ]
        m = d.merge(sc, on="key", how="left", validate="one_to_one", indicator=True)
        assert m._merge.eq("both").all()
        combined[gene] = m
    bounds = {}
    for col, _ in specs:
        vals = pd.concat([d[col] for d in combined.values()]).dropna()
        if vals.empty:
            bounds[col] = (0, 1)
        elif col == "alphamissense_frozen":
            bounds[col] = (-0.02, 1.02)
        else:
            lo, hi = vals.min(), vals.max()
            pad = max((hi - lo) * 0.04, 0.01)
            bounds[col] = (lo - pad, hi + pad)
    fig, axes = setup(
        3,
        "Fitted priors versus AlphaMissense, GPN-Star and AlphaGenome",
        "Frozen in-silico prior • genomic scores joined to the existing variant set • all score ranges shown without clipping",
    )
    metrics = []
    for i, gene in enumerate(GENES):
        d = combined[gene]
        for j, (col, label) in enumerate(specs):
            ax = axes[i, j]
            valid = d[col].notna() & np.isfinite(d[col])
            sub = d.loc[valid]
            scatter_groups(ax, sub, col)
            ax.set(xlim=bounds[col], xlabel=label)
            ax.xaxis.label.set_fontsize(9)
            if j == 0:
                ax.set_ylabel("Prior mean")
            r = rho(sub[col], sub.prior)
            ax.text(
                0.035,
                0.925,
                f"n = {len(sub):,} / {len(d):,}; ρ = {r:+.2f}",
                transform=ax.transAxes,
                fontsize=9,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.82,
                    "edgecolor": "none",
                    "pad": 1,
                },
            )
            metrics.append(
                {
                    "gene": gene,
                    "comparison": "all_available",
                    "predictor": col,
                    "n": len(sub),
                    "total_variants": len(d),
                    "spearman_rho": r,
                }
            )
        common = d.vclass.eq("missense") & d.mapped_allele_count.eq(1)
        for col, _ in specs:
            common &= d[col].notna() & np.isfinite(d[col])
        for col, _ in specs:
            sub = d.loc[common]
            metrics.append(
                {
                    "gene": gene,
                    "comparison": "same_single_allele_missense_variants",
                    "predictor": col,
                    "n": len(sub),
                    "total_variants": len(d),
                    "spearman_rho": rho(sub[col], sub.prior),
                }
            )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.928),
        ncol=3,
        frameon=False,
    )
    fig.text(
        0.08,
        0.045,
        "GPN-Star / AVI scalars require unambiguous allele scores; missing scores are omitted, not zero. AlphaMissense is the original prior input.",
        fontsize=9.5,
    )
    fig.text(
        0.08,
        0.025,
        "The prior uses AlphaMissense, consequence class and missingness. AVI itself combines AlphaGenome and AlphaMissense; agreement is descriptive.",
        fontsize=9.5,
    )
    fig.savefig(OUT / "PRIOR_VS_PREDICTORS.png", dpi=170, facecolor="white")
    plt.close(fig)
    pd.DataFrame(metrics).to_csv(OUT / "predictor_correlations.csv", index=False)
    return metrics


def main():
    if not (OUT / "predictor_scores.csv.gz").exists():
        raise FileNotFoundError(
            "The committed predictor_scores.csv.gz snapshot is required."
        )
    frames = load()
    distributions(frames)
    vs_counts(frames)
    metrics = vs_predictors(frames)
    (OUT / "figure_check.json").write_text(
        json.dumps(
            {
                "gnomad_assumed_unaffected": True,
                "row_counts": {g: len(d) for g, d in frames.items()},
                "predictor_metrics": metrics,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
