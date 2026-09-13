"""Plot residue summaries of the frozen variant-only LOO penetrance densities."""

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
from matplotlib.ticker import MaxNLocator, PercentFormatter


SOURCES = EVIDENCE / "structural_sanity_20260913/analysis/variant_diagnostics"
GENES = {
    "HNF1A": (
        "MODY3",
        "missense_structural_extension_20260912/geometry/HNF1A/P20823-1.fasta",
    ),
    "GCK": (
        "GCK-MODY (MODY2)",
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
STYLE = {
    "structured": ("#2c5578", "o", "3D: experimental or local AF"),
    "polymer": ("#bf791e", "^", "Polymer"),
    "mixed": ("#825487", "s", "Mixed contexts"),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_csv(data, path):
    data.to_csv(path, index=False, lineterminator="\n")
    assert path.stat().st_size < 1_200_000


def load_residues(gene, disease, fasta, hashes):
    hashes[str(fasta.relative_to(EVIDENCE))] = sha(fasta)
    sequence = "".join(
        line.strip()
        for line in fasta.read_text().splitlines()
        if not line.startswith(">")
    )
    paths = sorted(SOURCES.glob(f"{gene}.part*.csv.gz"))
    hashes.update({str(path.relative_to(EVIDENCE)): sha(path) for path in paths})
    variants = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    assert variants.gene.eq(gene).all() and variants.unit_id.is_unique
    assert variants.aa_pos.notna().all()
    assert np.equal(variants.aa_pos, variants.aa_pos.astype(int)).all()
    variants["canonical_pos"] = variants.aa_pos.astype(int)
    assert variants.canonical_pos.between(1, len(sequence)).all()
    assert variants.prior_mean.nunique() == 1
    assert variants.density.dropna().between(0, 1).all()
    prior = float(variants.prior_mean.iloc[0])
    rows = []
    for pos in range(1, len(sequence) + 1):
        group = variants.loc[variants.canonical_pos.eq(pos)]
        supported = group.loc[group.density.notna()]
        categories = supported.density_source.unique()
        assert len(categories) <= 1, (gene, pos, categories)
        status = (
            "supported"
            if len(supported)
            else "observed_without_density"
            if len(group)
            else "no_observed_variant"
        )
        rows.append(
            {
                "gene": gene,
                "intended_disease": disease,
                "canonical_pos": pos,
                "n_variants": len(group),
                "n_supported_variants": len(supported),
                "density_mean": supported.density.mean(),
                "density_min": supported.density.min(),
                "density_max": supported.density.max(),
                "density_source": categories[0] if len(categories) else "missing",
                "prior_mean": prior,
                "own_posterior_mean": group.posterior_mean.mean(),
                "affected": group.affected.sum(),
                "unaffected": group.unaffected.sum(),
                "gnomad_carriers": group.gnomad_carriers.sum(),
                "status": status,
            }
        )
    residues = pd.DataFrame(rows)
    assert residues.n_variants.sum() == len(variants)
    assert residues.n_supported_variants.sum() == variants.density.notna().sum()
    return residues


def draw_panel(ax, gene, disease, data):
    for source, (color, marker, _) in STYLE.items():
        group = data.loc[data.density_source.eq(source)]
        ax.vlines(
            group.canonical_pos,
            group.density_min,
            group.density_max,
            color=color,
            linewidth=0.65,
            alpha=0.35,
        )
        ax.scatter(
            group.canonical_pos,
            group.density_mean,
            color=color,
            marker=marker,
            s=10 if gene == "BRCA2" else 15,
            linewidths=0,
            alpha=0.82,
            zorder=3,
        )
    missing = data.loc[data.density_mean.isna()]
    # Missingness is drawn outside the score scale, never as a zero-valued score.
    ax.scatter(
        missing.canonical_pos,
        np.full(len(missing), -0.035),
        marker="|",
        color="#a0a5ab",
        s=17,
        linewidths=0.7,
    )
    prior = data.prior_mean.iloc[0]
    ax.axhline(prior, color="#303b46", linestyle="--", linewidth=0.9, zorder=2)
    ax.text(
        0.992,
        prior + 0.012,
        f"Prior {prior:.1%}",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        color="#303b46",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1},
    )
    supported = int(data.density_mean.notna().sum())
    units = int(data.n_supported_variants.sum())
    ax.set_title(
        f"{gene} · {disease}", loc="left", fontsize=13, fontweight="bold", pad=10
    )
    ax.text(
        1,
        1.035,
        f"{supported:,}/{len(data):,} residues · {units:,} variants",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="#4f5963",
    )
    ax.set(
        xlim=(1, len(data)),
        ylim=(-0.065, 1),
        ylabel="Penetrance density",
        xlabel="Canonical residue",
    )
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.xaxis.set_major_locator(MaxNLocator(8, integer=True, prune="both"))
    ticks = ax.get_xticks()
    ticks = [int(x) for x in ticks if 1 < x < len(data)]
    # Keep sequence endpoints legible without closely adjacent ticks.
    ticks = [x for x in ticks if x > len(data) * 0.04 and x < len(data) * 0.96]
    ax.set_xticks([1, *ticks, len(data)])
    ax.grid(axis="y", color="#e5e7ea", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#8c939b")


def legend(fig, y):
    handles = [
        Line2D(
            [],
            [],
            color=color,
            marker=marker,
            linestyle="none",
            markersize=5,
            label=label,
        )
        for color, marker, label in STYLE.values()
    ]
    handles.extend(
        [
            Line2D(
                [],
                [],
                color="#303b46",
                linestyle="--",
                linewidth=1,
                label="Missense prior",
            ),
            Line2D(
                [],
                [],
                color="#a0a5ab",
                marker="|",
                linestyle="none",
                markersize=7,
                label="No supported estimate",
            ),
        ]
    )
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        ncol=5,
        frameon=False,
        fontsize=9,
        handletextpad=0.5,
        columnspacing=1.4,
    )


def export(fig, stem):
    for ext in ["png", "pdf"]:
        path = HERE / f"{stem}.{ext}"
        kwargs = (
            {"dpi": 170}
            if ext == "png"
            else {"metadata": {"CreationDate": None, "ModDate": None}}
        )
        fig.savefig(path, **kwargs)
        assert path.stat().st_size < 1_200_000, (path, path.stat().st_size)
    plt.close(fig)


def main():
    (HERE / "tables").mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    hashes, tables, coverage = {}, {}, []
    for gene, (disease, fasta_name) in GENES.items():
        data = load_residues(gene, disease, EVIDENCE / fasta_name, hashes)
        tables[gene] = data
        save_csv(data, HERE / "tables" / f"{gene}_residue_density.csv")
        coverage.append(
            {
                "gene": gene,
                "intended_disease": disease,
                "canonical_residues": len(data),
                "variant_units": int(data.n_variants.sum()),
                "supported_variants": int(data.n_supported_variants.sum()),
                "supported_residues": int(data.density_mean.notna().sum()),
                "structured_residues": int(data.density_source.eq("structured").sum()),
                "polymer_residues": int(data.density_source.eq("polymer").sum()),
                "mixed_residues": int(data.density_source.eq("mixed").sum()),
                "observed_without_density": int(
                    data.status.eq("observed_without_density").sum()
                ),
                "no_observed_variant": int(data.status.eq("no_observed_variant").sum()),
                "prior_mean": float(data.prior_mean.iloc[0]),
            }
        )
        fig, ax = plt.subplots(figsize=(13.5, 4.6))
        fig.subplots_adjust(left=0.07, right=0.985, bottom=0.22, top=0.78)
        draw_panel(ax, gene, disease, data)
        legend(fig, 0.94)
        fig.text(
            0.07,
            0.105,
            "Dot = mean variant-only LOO density at that residue; span = variant range, not a confidence interval.",
            fontsize=10,
        )
        fig.text(
            0.07,
            0.052,
            "Intended disease label; counts remain pooled, especially GCK. gnomAD assumed unaffected. Neighborhood scores are not individual disease risks.",
            fontsize=9,
            color="#4f5963",
        )
        export(fig, f"{gene}_PENETRANCE_DENSITY_BY_RESIDUE")
    fig, axes = plt.subplots(5, 1, figsize=(13.5, 16))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.085, top=0.915, hspace=0.62)
    fig.suptitle(
        "Penetrance density by residue", fontsize=22, y=0.984, fontweight="bold"
    )
    fig.text(
        0.5,
        0.96,
        "Missense variants · corrected BRCA2 polymer coverage · fixed variant-only leave-one-out scores",
        ha="center",
        fontsize=11,
    )
    legend(fig, 0.948)
    for ax, (gene, (disease, _)) in zip(axes, GENES.items(), strict=True):
        draw_panel(ax, gene, disease, tables[gene])
    fig.text(
        0.075,
        0.047,
        "Dot = equal mean over supported variants at that residue. Thin span = between-variant range, not a confidence interval.",
        fontsize=11,
    )
    fig.text(
        0.075,
        0.03,
        "Disease names identify intended pairs; current clinical counts are pooled, particularly for GCK. gnomAD carriers are assumed unaffected.",
        fontsize=10,
        color="#4f5963",
    )
    fig.text(
        0.075,
        0.014,
        "These neighborhood features are not individual disease-risk estimates. Missing residues remain unestimated; no interpolation or extra smoothing.",
        fontsize=10,
        color="#4f5963",
    )
    export(fig, "ALL_GENES_PENETRANCE_DENSITY_BY_RESIDUE")
    save_csv(pd.DataFrame(coverage), HERE / "tables/coverage.csv")
    (HERE / "inputs.json").write_text(
        json.dumps(
            {
                "input_hashes_relative_to_evidence": hashes,
                "plot_script_sha256": sha(Path(__file__)),
                "aggregation": "Equal mean of supported variant-only LOO densities within each canonical residue; min/max are between-variant range, not uncertainty bounds. No new fitting, donor changes or count updates.",
                "y_scale": [0, 1],
            },
            indent=2,
        )
        + "\n"
    )
    print(pd.DataFrame(coverage).to_string(index=False))


if __name__ == "__main__":
    main()
