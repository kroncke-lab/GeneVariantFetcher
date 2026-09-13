"""Show corrected BRCA2 coverage, local score and its prior contribution."""

import os
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE.parents[2] / "tmp/matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    out = HERE / "analysis"
    data = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted((out / "BRCA2").glob("primary_density.part*.csv.gz"))
        ],
        ignore_index=True,
    )
    geometry = pd.read_csv(
        HERE / "geometry/BRCA2/primary_geometry.csv.gz", low_memory=False
    )
    structured = (
        geometry.loc[geometry.geometry_state.eq("structured")]
        .groupby("canonical_pos")
        .geometry_source.first()
    )
    data["source"] = np.where(
        data.density_source.eq("polymer"), "Polymer", "Unavailable"
    )
    data.loc[data.density_source.eq("structured"), "source"] = data.loc[
        data.density_source.eq("structured"), "canonical_pos"
    ].map(structured)
    # Source names come from the validated geometry ledger, not coordinate magnitude.
    for source in data.source.unique():
        if "experiment" in str(source).lower():
            data.loc[data.source.eq(source), "source"] = "Experimental 3D"
        elif "alpha" in str(source).lower() or "af_" in str(source).lower():
            data.loc[data.source.eq(source), "source"] = "AF local 3D"
    supported = data.loc[data.density.notna()]
    prior = (
        pd.read_csv(
            HERE.parent
            / "class_matched_penetrance_20260912/analysis/empirical_prior_comparison.csv"
        )
        .query("gene == 'BRCA2' and scope == 'canonical_missense'")
        .iloc[0]["mean"]
    )
    old = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted(
                (
                    HERE.parent
                    / "missense_structural_extension_20260912/analysis/BRCA2"
                ).glob("primary_density.part*.csv.gz")
            )
        ]
    )
    palette = {
        "Polymer": "#bf791e",
        "Experimental 3D": "#2c5578",
        "AF local 3D": "#278f8e",
        "Unavailable": "#aaaaaa",
    }
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(
        3, 1, figsize=(14, 9), sharex=True, constrained_layout=True
    )
    old = old.loc[old.density.notna()]
    axes[0].scatter(
        old.canonical_pos, old.density, s=10, color=palette["Experimental 3D"]
    )
    axes[0].set(
        title=f"Previous: {len(old):,} / {len(data):,} variants · two small experimental regions",
        ylabel="Neighborhood score",
        ylim=(0, 0.75),
    )
    for source, group in supported.groupby("source"):
        color = palette.get(source, "#75528b")
        axes[1].vlines(
            group.canonical_pos,
            group.density_lower_95,
            group.density_upper_95,
            color=color,
            alpha=0.08,
            linewidth=0.5,
        )
        axes[1].scatter(
            group.canonical_pos,
            group.density,
            color=color,
            s=7,
            alpha=0.6,
            label=f"{source}: {len(group):,}",
        )
        axes[2].scatter(
            group.canonical_pos, group.density_minus_prior, color=color, s=7, alpha=0.6
        )
    unavailable = data.loc[data.density.isna()]
    axes[1].scatter(
        unavailable.canonical_pos,
        np.full(len(unavailable), -0.025),
        color=palette["Unavailable"],
        marker="|",
        s=12,
        label=f"Unavailable: {len(unavailable):,}",
    )
    axes[1].set(
        title=f"Corrected: {len(supported):,} / {len(data):,} variants · separate local frames and one canonical polymer layer",
        ylabel="Neighborhood score",
        ylim=(-0.04, 0.75),
    )
    for ax in axes[:2]:
        ax.axhline(prior, color="#263747", linestyle="--", linewidth=1)
        ax.grid(axis="y", alpha=0.15)
    axes[2].axhline(0, color="#263747", linestyle="--", linewidth=1)
    axes[2].set(
        title="Local departure from the shared missense prior",
        xlabel="Canonical BRCA2 residue",
        ylabel="Score − prior",
        xlim=(1, 3418),
    )
    axes[2].grid(axis="y", alpha=0.15)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle(
        "BRCA2 coverage correction · empirical-posterior neighborhood scores\nDashed line: fixed prior; thin bars: conditional 95% intervals. Scores are not individual cancer-risk estimates.",
        fontsize=14,
    )
    fig.savefig(out / "BRCA2_POLYMER_CORRECTION.png", dpi=150)
    plt.close(fig)
    pd.DataFrame(data.groupby("source").size(), columns=["variants"]).to_csv(
        out / "BRCA2_figure_source_counts.csv", lineterminator="\n"
    )


if __name__ == "__main__":
    main()
