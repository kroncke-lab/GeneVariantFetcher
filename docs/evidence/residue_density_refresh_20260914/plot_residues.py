"""Render inspectable residue plots; distinguish posterior features from counts."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
REPO = HERE.parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(REPO / "tmp/matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator, PercentFormatter


PAIRS = {
    "HNF1A": (
        "MODY3",
        "missense_structural_extension_20260912/geometry/HNF1A/P20823-1.fasta",
    ),
    "GCK": (
        "GCK-MODY / mild hyperglycemia",
        "gck_structural_pilot_20260912/geometry/P35557-1.fasta",
    ),
    "LDLR": (
        "Familial hypercholesterolemia",
        "missense_structural_extension_20260912/geometry/LDLR/P01130-1.fasta",
    ),
    "BRCA2": (
        "Hereditary breast/ovarian cancer susceptibility",
        "structural_sanity_20260913/geometry/BRCA2/P51587-1.fasta",
    ),
    "KCNQ1": (
        "Long QT syndrome type 1",
        "missense_structural_extension_20260912/geometry/KCNQ1/P51787-1.fasta",
    ),
}
BLUE, GOLD, PURPLE = "#175F9B", "#AD6919", "#80528E"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "pdf.fonttype": 42,
        "axes.titleweight": "bold",
    }
)
HASHES = {}


def read(path):
    HASHES[str(path.relative_to(EVIDENCE))] = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()
    return pd.read_csv(path, low_memory=False)


def read_primary(gene, suffix=""):
    paths = sorted(
        (HERE / "analysis" / (gene + suffix)).glob("density_h3.part*.csv.gz")
    )
    assert paths, gene + suffix
    frame = pd.concat([read(path) for path in paths], ignore_index=True)
    assert frame.variant_id.is_unique and frame.gene.eq(gene).all()
    return frame


def previous_density(gene):
    if gene == "BRCA2":
        paths = sorted(
            (EVIDENCE / "brca2_distance_prior_audit_20260913/corrected").glob(
                "primary_density.part*.csv.gz"
            )
        )
    elif gene == "GCK":
        paths = [
            EVIDENCE
            / "class_matched_penetrance_20260912/structural/GCK_primary_variant_density.csv"
        ]
    else:
        paths = sorted(
            (EVIDENCE / "missense_structural_extension_20260912/analysis" / gene).glob(
                "primary_density.part*.csv.gz"
            )
        )
    return pd.concat([read(path) for path in paths], ignore_index=True)


def residue_table(gene, frame, fasta, prior):
    HASHES[str(fasta.relative_to(EVIDENCE))] = hashlib.sha256(
        fasta.read_bytes()
    ).hexdigest()
    sequence = "".join(
        line.strip()
        for line in fasta.read_text().splitlines()
        if not line.startswith(">")
    )
    assert frame.canonical_pos.between(1, len(sequence)).all()
    result = []
    for position in range(1, len(sequence) + 1):
        group = frame.loc[frame.canonical_pos.eq(position)]
        supported = group.loc[group.density.notna()]
        sources = supported.density_source.unique()
        result.append(
            dict(
                gene=gene,
                intended_disease=PAIRS[gene][0],
                canonical_pos=position,
                canonical_aa=sequence[position - 1],
                variants=len(group),
                supported_variants=len(supported),
                status="supported"
                if len(supported)
                else "observed_without_density"
                if len(group)
                else "no_observed_variant",
                density_source=sources[0]
                if len(sources) == 1
                else "mixed"
                if len(sources)
                else "missing",
                density=supported.density.mean(),
                density_min=supported.density.min(),
                density_max=supported.density.max(),
                raw_variant_fraction_density=supported.raw_variant_fraction_density.mean(),
                raw_kernel_pooled_fraction=supported.raw_kernel_pooled_fraction.mean(),
                raw_kernel_one_prior_posterior=supported.raw_kernel_one_prior_posterior.mean(),
                own_posterior_mean=group.posterior_mean.mean(),
                own_observed_fraction_mean=(group.affected / group.n).mean(),
                affected=group.affected.sum(),
                unaffected_literature=group.unaffected_literature.sum(),
                gnomad_carriers=group.gnomad_carriers.sum(),
                unaffected=group.unaffected.sum(),
                mean_prior_retention=supported.neighborhood_prior_retention.mean(),
                mean_prior_component=supported.prior_component.mean(),
                mean_counts_component=supported.counts_component.mean(),
                mean_kish_donor_n=supported.kish_donor_n.mean(),
                mean_raw_kernel_mass=supported.raw_kernel_total_weight.mean(),
                mean_weight_share_beyond_20=supported.weight_share_beyond_20.mean(),
                no_affected_donor_variants=int(
                    supported.affected_donor_count.eq(0).sum()
                ),
                prior_mean=prior["mean"],
                alpha_empirical=prior.alpha_empirical,
                beta_empirical=prior.beta_empirical,
            )
        )
    data = pd.DataFrame(result)
    assert data.variants.sum() == len(frame)
    assert data.supported_variants.sum() == frame.density.notna().sum()
    return data


def axes_style(ax, length, ylabel):
    ax.set_xlim(1, length)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(7, integer=True))
    ticks = [int(x) for x in ax.get_xticks() if length * 0.05 < x < length * 0.95]
    ax.set_xticks([1, *ticks, length])
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.16, linewidth=0.65)
    ax.set_axisbelow(True)


def primary_panel(ax, gene, data, detailed=False):
    x = data.canonical_pos.to_numpy()
    y = data.density.to_numpy()
    ax.fill_between(
        x, data.density_min, data.density_max, color=BLUE, alpha=0.11, linewidth=0
    )
    ax.plot(
        x,
        y,
        color=BLUE,
        marker=".",
        markersize=2.8 if gene == "BRCA2" else 3.4,
        linewidth=0.85,
        label="Posterior neighborhood · primary",
    )
    ax.plot(
        x,
        data.raw_variant_fraction_density,
        color=GOLD,
        linestyle="--",
        marker=".",
        markersize=1.7,
        linewidth=0.6,
        alpha=0.7,
        label="Mean observed variant fraction · diagnostic",
    )
    ax.axhline(
        data.prior_mean.iloc[0], color="#434950", linestyle=(0, (4, 3)), linewidth=0.8
    )
    ax.text(
        0.99,
        0.96,
        f"Prior {data.prior_mean.iloc[0]:.1%} · {data.density.notna().sum():,}/{len(data):,} residues",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
    )
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    axes_style(ax, len(data), "Neighborhood score")
    ax.set_title(f"{gene} · {PAIRS[gene][0]}", loc="left", fontsize=12, pad=8)
    # Support track is outside the score scale; missing positions are not zero.
    for source, color, marker in [
        ("structured", BLUE, "|"),
        ("polymer", GOLD, "|"),
        ("mixed", PURPLE, "|"),
        ("missing", "#B1B4B8", "."),
    ]:
        positions = data.loc[data.density_source.eq(source), "canonical_pos"]
        ax.scatter(
            positions,
            np.full(len(positions), -0.055),
            transform=ax.get_xaxis_transform(),
            color=color,
            marker=marker,
            s=10,
            linewidths=0.65,
            clip_on=False,
        )
    if detailed:
        ax.set_xlabel(
            "Canonical residue · track: blue 3D, gold polymer, purple mixed, grey unestimated",
            fontsize=9,
        )


def plot_individual(gene, data, frame):
    fig, axes = plt.subplots(
        4, 1, figsize=(12, 11.6), gridspec_kw={"height_ratios": [1.05, 1, 0.7, 0.7]}
    )
    primary_panel(axes[0], gene, data, detailed=True)
    axes[0].legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.19),
        ncol=2,
        fontsize=8.4,
        frameon=False,
    )
    log_columns = [
        "density",
        "raw_variant_fraction_density",
        "raw_kernel_one_prior_posterior",
    ]
    positive_values = data[log_columns].to_numpy()
    positive_values = positive_values[
        np.isfinite(positive_values) & (positive_values > 0)
    ]
    zero_display = min(0.00001, positive_values.min() / 3)
    for column, label, color, marker in [
        ("density", "Posterior neighborhood", BLUE, "."),
        ("raw_variant_fraction_density", "Mean observed variant fraction", GOLD, "."),
        (
            "raw_kernel_one_prior_posterior",
            "Pooled counts + one prior / context",
            PURPLE,
            ".",
        ),
    ]:
        value = data[column]
        positive = value.gt(0)
        axes[1].plot(
            data.canonical_pos[positive],
            value[positive],
            linestyle="none",
            marker=marker,
            ms=3,
            color=color,
            label=label,
        )
        zeros = value.eq(0)
        axes[1].scatter(
            data.canonical_pos[zeros],
            np.full(zeros.sum(), zero_display),
            marker="v",
            color=color,
            s=14,
        )
    axes[1].set_yscale("log")
    axes[1].set_ylim(zero_display * 0.65, 1)
    axes[1].yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value * 100:g}%")
    )
    axes[1].axhline(0.001, linestyle="--", color="#50565D", linewidth=0.8)
    axes[1].text(
        1.005,
        0.001,
        "0.1%",
        transform=axes[1].get_yaxis_transform(),
        va="center",
        fontsize=8,
    )
    axes[1].legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=3,
        fontsize=8,
        frameon=False,
    )
    axes_style(axes[1], len(data), "Neighborhood score · log scale")
    axes[2].scatter(
        frame.canonical_pos,
        frame.posterior_mean,
        facecolors="none",
        edgecolors=BLUE,
        linewidths=0.55,
        s=10,
        alpha=0.65,
        label="Each variant's own count posterior",
    )
    axes[2].plot(
        data.canonical_pos,
        data.own_posterior_mean,
        color="#454C53",
        linewidth=0.6,
        marker=".",
        markersize=2.5,
        label="Own posterior mean at residue",
    )
    axes[2].set_ylim(0, 1)
    axes[2].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    axes[2].legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=2,
        fontsize=8,
        frameon=False,
    )
    axes_style(axes[2], len(data), "Own posterior")
    count_ax = axes[3]
    for column, color, label, marker in [
        ("affected", "#343A40", "Affected observations", "o"),
        (
            "unaffected",
            BLUE,
            "Unaffected observations (gnomAD assumed unaffected)",
            ".",
        ),
    ]:
        positive = data[column].gt(0)
        count_ax.vlines(
            data.canonical_pos[positive],
            1,
            data.loc[positive, column],
            color=color,
            alpha=0.25,
            linewidth=0.6,
        )
        count_ax.scatter(
            data.canonical_pos[positive],
            data.loc[positive, column],
            color=color,
            s=8,
            marker=marker,
            label=label,
        )
    count_ax.set_yscale("log")
    count_ax.set_ylim(0.7, None)
    count_ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=2,
        fontsize=8,
        frameon=False,
    )
    axes_style(count_ax, len(data), "Observed counts / residue")
    count_ax.set_xlabel("Canonical residue")
    fig.suptitle(
        f"{gene} · missense penetrance-density profile",
        x=0.08,
        ha="left",
        fontsize=16,
        weight="bold",
        y=0.987,
    )
    fig.text(
        0.08,
        0.018,
        "Variant-only leave-one-out; same-residue alternatives retained. Half weight at 3 Å; positive tail; same-IDR polymer distances.\nLine = mean of supported variant scores; pale band = between-variant range, not uncertainty. Gaps are unestimated.\n"
        f"Log panel: ▼ = exact zero, displayed at {zero_display * 100:.3g}%, below every positive value.\n"
        "These are descriptive features, not calibrated disease probabilities. Source exclusions are documented in the accompanying report.",
        fontsize=8.5,
    )
    fig.subplots_adjust(left=0.12, right=0.95, top=0.865, bottom=0.14, hspace=0.51)
    for ext in ["png", "pdf"]:
        fig.savefig(HERE / "plots" / f"{gene}_RESIDUE_DENSITY.{ext}", dpi=150)
    plt.close(fig)


def main():
    (HERE / "plots").mkdir(exist_ok=True)
    (HERE / "tables").mkdir(exist_ok=True)
    tables, metrics, prior_rows, original_frames = {}, [], [], {}
    for gene, (_, fasta) in PAIRS.items():
        frame = read_primary(gene)
        original_frames[gene] = frame
        before = previous_density(gene)
        prior = read(HERE / "analysis" / gene / "prior_comparison.csv")
        prior_rows.append(prior)
        current = prior.loc[
            prior.scenario.eq("refreshed") & prior.variant_type.eq("missense")
        ].iloc[0]
        data = residue_table(gene, frame, EVIDENCE / fasta, current)
        tables[gene] = data
        data.to_csv(
            HERE / "tables" / f"{gene}_residue_density.csv",
            index=False,
            lineterminator="\n",
        )
        plot_individual(gene, data, frame)
        for label, source in [("previous", before), ("refreshed", frame)]:
            columns = (
                ["density"]
                if label == "previous"
                else [
                    "density",
                    "raw_variant_fraction_density",
                    "raw_kernel_pooled_fraction",
                    "raw_kernel_one_prior_posterior",
                ]
            )
            for column in columns:
                value = source[column].dropna()
                metrics.append(
                    dict(
                        gene=gene,
                        scenario=label,
                        score=column,
                        variants=len(source),
                        supported=len(value),
                        supported_residues=source.loc[
                            source[column].notna(), "canonical_pos"
                        ].nunique(),
                        mean=value.mean(),
                        median=value.median(),
                        minimum=value.min(),
                        maximum=value.max(),
                        at_or_below_0p1_percent=int(value.le(0.001).sum()),
                        exact_zero=int(value.eq(0).sum()),
                    )
                )
    fig, axes = plt.subplots(5, 1, figsize=(12, 14))
    for axis, (gene, data) in zip(axes, tables.items()):
        primary_panel(axis, gene, data)
        axis.set_xlabel("Canonical residue", fontsize=9, labelpad=9)
    handles = [
        Line2D(
            [],
            [],
            color=BLUE,
            marker=".",
            linewidth=1,
            label="Posterior neighborhood · primary",
        ),
        Line2D(
            [],
            [],
            color=GOLD,
            linestyle="--",
            linewidth=1,
            label="Observed-fraction neighborhood · diagnostic",
        ),
        Line2D(
            [],
            [],
            color="#434950",
            linestyle="--",
            linewidth=1,
            label="Gene/missense prior",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.953),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "Missense penetrance density by canonical residue",
        x=0.08,
        ha="left",
        y=0.983,
        fontsize=17,
        weight="bold",
    )
    fig.text(
        0.08,
        0.015,
        "Variant-only leave-one-out · h=3 Å · gnomAD assumed unaffected. Track: blue 3D, gold polymer, purple mixed, grey unestimated.\nMeans of supported variant scores; pale range is between variants, not a confidence interval. Lines break at unestimated residues.\nSource and endpoint corrections applied. These neighborhood features are not calibrated individual disease probabilities.",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.08, right=0.97, top=0.905, bottom=0.085, hspace=0.58)
    for ext in ["png", "pdf"]:
        fig.savefig(HERE / "plots" / f"ALL_GENES_RESIDUE_DENSITY.{ext}", dpi=140)
    plt.close(fig)
    pd.DataFrame(metrics).to_csv(
        HERE / "tables/score_summary.csv", index=False, lineterminator="\n"
    )
    pd.concat(prior_rows, ignore_index=True).to_csv(
        HERE / "tables/prior_summary.csv", index=False, lineterminator="\n"
    )
    for gene, endpoint in [("GCK", "diabetes_proxy"), ("KCNQ1", "cardiac_events")]:
        sensitivity = read_primary(gene, "_" + endpoint)
        merged = original_frames[gene][["variant_id", "density"]].merge(
            sensitivity[["variant_id", "density"]],
            on="variant_id",
            suffixes=("_primary", "_sensitivity"),
            validate="one_to_one",
        )
        delta = abs(merged.density_primary - merged.density_sensitivity)
        result = dict(
            endpoint=endpoint,
            common_variants=len(merged),
            common_supported=int(delta.notna().sum()),
            primary_median=original_frames[gene].density.median(),
            sensitivity_median=sensitivity.density.median(),
            mean_abs_delta_common=delta.mean(),
            max_abs_delta_common=delta.max(),
        )
        (HERE / "tables" / f"{gene}_{endpoint}_sensitivity.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
    HASHES[str(Path(__file__).relative_to(EVIDENCE))] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    (HERE / "plot_input_hashes.json").write_text(json.dumps(HASHES, indent=2) + "\n")


if __name__ == "__main__":
    main()
