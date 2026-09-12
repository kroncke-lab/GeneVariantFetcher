#!/usr/bin/env python3
"""Replot archived protocol scores by their actual evidence, without refitting."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / "grant_e2e_20260909" / "analysis"
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]


def main():
    frames = {
        gene: pd.read_csv(SOURCE / f"{gene}_protocol" / "posterior_variants.csv")
        for gene in GENES
    }
    d = frames["BRCA2"]
    singleton = d["n"].eq(1) & d["affected"].eq(1) & d["gnomad_added"].eq(0)
    masks = [singleton, d["gnomad_added"].eq(0) & ~singleton, d["gnomad_added"].gt(0)]
    assert sum(int(m.sum()) for m in masks) == len(d)
    labels = [
        "One affected observation; no added gnomAD count",
        "Other observations; no added gnomAD count",
        "Positive gnomAD count added as unaffected",
    ]
    colors = ["#d48c32", "#aab2b9", "#306fa2"]
    bins = np.linspace(0, 1, 51)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(13.8, 5.8), gridspec_kw={"width_ratios": [1.5, 1]}
    )
    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.24, top=0.78, wspace=0.27)
    fig.suptitle(
        "Sparse observations and shared priors shape the protocol histogram",
        x=0.06,
        ha="left",
        y=0.98,
        fontsize=17,
    )
    fig.text(
        0.06,
        0.915,
        "Archived 2026-09-09 run • literature-selected variants • diagnostic scores, not calibrated disease risks",
        fontsize=11,
    )
    ax.hist(
        [d.loc[m, "insilico_post_mean"] for m in masks],
        bins=bins,
        stacked=True,
        color=colors,
        label=labels,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.set(
        title="BRCA2: 5,962 variants, separated by available counts",
        xlabel="Archived protocol posterior mean",
        ylabel="Number of variants",
        xlim=(0, 1),
    )
    ax.legend(loc="upper right", fontsize=8.3, frameon=False)
    shares = np.array([100 * frames[g]["n"].lt(10).mean() for g in GENES])
    y = np.arange(len(GENES))
    bx.barh(y, shares, color="#d48c32", height=0.57)
    bx.set(
        yticks=y,
        yticklabels=GENES,
        xlim=(0, 100),
        xlabel="Percent of plotted variants",
        title="The prior outweighs the counts for most variants",
    )
    bx.invert_yaxis()
    for i, share in enumerate(shares):
        count = int(frames[GENES[i]]["n"].lt(10).sum())
        bx.text(share + 1.3, i, f"{share:.0f}%", va="center", fontsize=11)
        bx.text(
            2,
            i,
            f"{count:,} / {len(frames[GENES[i]]):,}",
            va="center",
            color="#15191c",
            fontsize=9,
        )
    for axis in (ax, bx):
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_axisbelow(True)
        axis.grid(axis="y" if axis is ax else "x", color="#e5e7e9", linewidth=0.7)
    fig.text(
        0.06,
        0.14,
        "Prior weight = 10 / (10 + n). Bars at right count n < 10; n includes gnomAD allele counts treated as unaffected by the archived model.",
        fontsize=10,
    )
    fig.text(
        0.06,
        0.095,
        "A single affected observation with prior 0.003 gives (1 + 10 × 0.003) / 11 ≈ 0.094; with prior 0.255 it gives ≈ 0.323.",
        fontsize=10,
    )
    fig.text(
        0.06,
        0.05,
        "gnomAD phenotypes are unknown. Zero added count can mean an absent annotation or failed join; it does not prove population absence.",
        fontsize=10,
    )
    fig.savefig(OUT / "PROTOCOL_DIAGNOSTIC.png", dpi=180, facecolor="white")
    plt.close(fig)
    rows = []
    for label, mask in zip(labels, masks):
        counts, _ = np.histogram(d.loc[mask, "insilico_post_mean"], bins=bins)
        rows.extend(
            {
                "stratum": label,
                "bin_left": left,
                "bin_right": right,
                "count": int(count),
            }
            for left, right, count in zip(bins[:-1], bins[1:], counts)
        )
    pd.DataFrame(rows).to_csv(OUT / "diagnostic_histogram_bins.csv", index=False)
    assert sum(r["count"] for r in rows) == len(d)


if __name__ == "__main__":
    main()
