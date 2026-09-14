"""Recompute every residue-table field and final plot summaries from variant rows."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_refresh import GEOMETRY, REFRESH, EVIDENCE, read, close, sha

HERE = Path(__file__).resolve().parent
FASTAS = {
    "GCK": "gck_structural_pilot_20260912/geometry/P35557-1.fasta",
    "HNF1A": "missense_structural_extension_20260912/geometry/HNF1A/P20823-1.fasta",
    "LDLR": "missense_structural_extension_20260912/geometry/LDLR/P01130-1.fasta",
    "KCNQ1": "missense_structural_extension_20260912/geometry/KCNQ1/P51787-1.fasta",
    "BRCA2": "structural_sanity_20260913/geometry/BRCA2/P51587-1.fasta",
}


def main():
    receipts = []
    metrics = pd.read_csv(REFRESH / "tables/score_summary.csv")
    for gene in GEOMETRY:
        data = read(REFRESH / "analysis" / gene, "density_h3")
        residue = pd.read_csv(REFRESH / "tables" / f"{gene}_residue_density.csv")
        sequence = "".join(
            line.strip()
            for line in (EVIDENCE / FASTAS[gene]).read_text().splitlines()
            if not line.startswith(">")
        )
        assert residue.canonical_pos.tolist() == list(range(1, len(sequence) + 1))
        assert "".join(residue.canonical_aa) == sequence
        assert residue.gene.eq(gene).all()
        counts = (
            data.groupby("canonical_pos")
            .size()
            .reindex(residue.canonical_pos, fill_value=0)
        )
        support = (
            data.loc[data.density.notna()]
            .groupby("canonical_pos")
            .size()
            .reindex(residue.canonical_pos, fill_value=0)
        )
        close(counts, residue.variants)
        close(support, residue.supported_variants)
        assert np.array_equal(
            residue.status,
            np.where(
                support > 0,
                "supported",
                np.where(counts > 0, "observed_without_density", "no_observed_variant"),
            ),
        )
        assert residue.density.isna().eq(support.eq(0).to_numpy()).all()
        data["own_observed_fraction"] = data.affected / data.n
        means = {
            "density": "density",
            "raw_variant_fraction_density": "raw_variant_fraction_density",
            "raw_kernel_pooled_fraction": "raw_kernel_pooled_fraction",
            "raw_kernel_one_prior_posterior": "raw_kernel_one_prior_posterior",
            "own_posterior_mean": "posterior_mean",
            "own_observed_fraction_mean": "own_observed_fraction",
            "mean_prior_retention": "neighborhood_prior_retention",
            "mean_prior_component": "prior_component",
            "mean_counts_component": "counts_component",
            "mean_kish_donor_n": "kish_donor_n",
            "mean_raw_kernel_mass": "raw_kernel_total_weight",
            "mean_weight_share_beyond_20": "weight_share_beyond_20",
        }
        grouped = data.groupby("canonical_pos")
        for displayed, original in means.items():
            scoped = (
                grouped
                if displayed.startswith("own_")
                else data.loc[data.density.notna()].groupby("canonical_pos")
            )
            close(
                scoped[original].mean().reindex(residue.canonical_pos),
                residue[displayed],
            )
        close(grouped.density.min().reindex(residue.canonical_pos), residue.density_min)
        close(grouped.density.max().reindex(residue.canonical_pos), residue.density_max)
        for field in [
            "affected",
            "unaffected_literature",
            "gnomad_carriers",
            "unaffected",
        ]:
            close(
                grouped[field].sum().reindex(residue.canonical_pos, fill_value=0),
                residue[field],
            )
        zeros = (
            data.loc[data.affected_donor_count.eq(0) & data.density.notna()]
            .groupby("canonical_pos")
            .size()
        )
        close(
            zeros.reindex(residue.canonical_pos, fill_value=0),
            residue.no_affected_donor_variants,
        )
        for row in residue.itertuples():
            sources = data.loc[
                data.canonical_pos.eq(row.canonical_pos) & data.density.notna(),
                "density_source",
            ].unique()
            expected = (
                sources[0]
                if len(sources) == 1
                else "mixed"
                if len(sources)
                else "missing"
            )
            assert row.density_source == expected
        prior = pd.read_csv(REFRESH / "analysis" / gene / "prior_comparison.csv")
        prior = prior.loc[
            prior.scenario.eq("refreshed") & prior.variant_type.eq("missense")
        ].iloc[0]
        close(residue.prior_mean, prior["mean"])
        close(residue.alpha_empirical, prior.alpha_empirical)
        close(residue.beta_empirical, prior.beta_empirical)
        for row in metrics.loc[
            metrics.gene.eq(gene) & metrics.scenario.eq("refreshed")
        ].itertuples():
            values = data[row.score].dropna()
            close(
                [row.mean, row.median, row.minimum, row.maximum],
                [values.mean(), values.median(), values.min(), values.max()],
            )
            assert [
                row.variants,
                row.supported,
                row.supported_residues,
                row.at_or_below_0p1_percent,
                row.exact_zero,
            ] == [
                len(data),
                len(values),
                data.loc[data[row.score].notna(), "canonical_pos"].nunique(),
                values.le(0.001).sum(),
                values.eq(0).sum(),
            ]
        receipts.append(
            dict(
                gene=gene,
                residue_rows=len(residue),
                supported_residues=int((support > 0).sum()),
                observed_without_density=int(((counts > 0) & support.eq(0)).sum()),
                no_observed_variant=int(counts.eq(0).sum()),
                variant_multiplicity_preserved=True,
                all_numeric_fields_match=True,
                zero_missing_and_unobserved_distinguished=True,
            )
        )
    for gene, endpoint in [("GCK", "diabetes_proxy"), ("KCNQ1", "cardiac_events")]:
        primary = read(REFRESH / "analysis" / gene, "density_h3")
        proxy = read(REFRESH / "analysis" / (gene + "_" + endpoint), "density_h3")
        common = primary[["variant_id", "density"]].merge(
            proxy[["variant_id", "density"]],
            on="variant_id",
            suffixes=("_primary", "_proxy"),
            validate="one_to_one",
        )
        delta = (common.density_primary - common.density_proxy).abs()
        sensitivity = json.loads(
            (REFRESH / "tables" / f"{gene}_{endpoint}_sensitivity.json").read_text()
        )
        assert sensitivity["common_variants"] == len(common)
        assert sensitivity["common_supported"] == delta.notna().sum()
        close(
            [
                sensitivity[k]
                for k in [
                    "primary_median",
                    "sensitivity_median",
                    "mean_abs_delta_common",
                    "max_abs_delta_common",
                ]
            ],
            [
                primary.density.median(),
                proxy.density.median(),
                delta.mean(),
                delta.max(),
            ],
        )
    hashes = json.loads((REFRESH / "plot_input_hashes.json").read_text())
    for relative, digest in hashes.items():
        assert sha(EVIDENCE / relative) == digest, relative
    result = dict(
        genes=receipts,
        plot_input_hashes_checked=len(hashes),
        gck_and_kcnq1_sensitivity_comparisons_verified=True,
        audit_script_sha256=sha(__file__),
        plot_script_sha256=sha(REFRESH / "plot_residues.py"),
    )
    (HERE / "plot_table_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
