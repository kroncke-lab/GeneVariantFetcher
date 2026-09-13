#!/usr/bin/env python3
"""Audit residue plot data against frozen variants, without importing plot code."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EVIDENCE = REPO / "docs/evidence"
CLASS = EVIDENCE / "class_matched_penetrance_20260912"
EXTENSION = EVIDENCE / "missense_structural_extension_20260912"
SANITY = EVIDENCE / "structural_sanity_20260913"
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
LENGTHS = dict(HNF1A=631, GCK=465, LDLR=860, BRCA2=3418, KCNQ1=676)


def read_shards(paths):
    paths = sorted(paths)
    assert paths, "No source shards found"
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True), paths


def audit_gene(gene):
    fasta = (
        EVIDENCE
        / "population_inclusive_penetrance_20260912/population"
        / f"{gene}_canonical.fasta"
    )
    sequence = "".join(
        line.strip()
        for line in fasta.read_text().splitlines()
        if not line.startswith(">")
    )
    assert len(sequence) == LENGTHS[gene]
    if gene == "GCK":
        paths = [CLASS / "structural/GCK_structural_input_variants.csv.gz"]
        inputs, paths = read_shards(paths)
        density_paths = [CLASS / "structural/GCK_primary_variant_density.csv"]
    else:
        inputs, paths = read_shards(
            (EXTENSION / "analysis" / gene).glob("input_variants.part*.csv.gz")
        )
        source_folder = SANITY if gene == "BRCA2" else EXTENSION
        density_paths = list(
            (source_folder / "analysis" / gene).glob("primary_density.part*.csv.gz")
        )
    density, density_paths = read_shards(density_paths)
    assert inputs.variant_id.is_unique and density.variant_id.is_unique
    assert set(inputs.variant_id) == set(density.variant_id)
    d = inputs.set_index("variant_id")
    den = density.set_index("variant_id").loc[d.index]
    np.testing.assert_array_equal(d.canonical_pos, den.canonical_pos)
    d["density"] = den.density
    d["density_source"] = den.density_source
    d["canonical_pos"] = d.canonical_pos.astype(int)
    assert d.canonical_pos.between(1, len(sequence)).all()
    assert [sequence[p - 1] for p in d.canonical_pos] == d.aa_ref.to_list()

    # Check the compact plotting input independently against older source counts
    # and latest density files, so a shared summary error cannot silently pass.
    compact, compact_paths = read_shards(
        (SANITY / "analysis/variant_diagnostics").glob(f"{gene}.part*.csv.gz")
    )
    assert compact.unit_id.is_unique and set(compact.unit_id) == set(d.index)
    compact = compact.set_index("unit_id").loc[d.index]
    for col in [
        "affected",
        "unaffected",
        "gnomad_carriers",
        "posterior_mean",
        "density",
    ]:
        np.testing.assert_allclose(
            compact[col], d[col], rtol=1e-10, atol=1e-11, equal_nan=True
        )
    np.testing.assert_array_equal(compact.aa_pos, d.canonical_pos)
    prior_path = CLASS / "analysis/empirical_prior_comparison.csv"
    priors = pd.read_csv(prior_path)
    prior = priors.loc[
        priors.gene.eq(gene) & priors.scope.eq("canonical_missense"), "mean"
    ].item()
    np.testing.assert_allclose(compact.prior_mean, prior, atol=1e-11)

    table_path = HERE / "tables" / f"{gene}_residue_density.csv"
    table = pd.read_csv(table_path)
    assert len(table) == len(sequence)
    np.testing.assert_array_equal(table.canonical_pos, np.arange(1, len(sequence) + 1))
    assert table.gene.eq(gene).all() and table.intended_disease.nunique() == 1
    assert table.intended_disease.notna().all()
    expected = []
    for pos in range(1, len(sequence) + 1):
        x = d.loc[d.canonical_pos.eq(pos)]
        supported = x.loc[x.density.notna()]
        sources = set(supported.density_source)
        assert len(sources) <= 1, (gene, pos, sources)
        source = next(iter(sources)) if sources else "missing"
        expected.append(
            dict(
                canonical_pos=pos,
                n_variants=len(x),
                n_supported_variants=len(supported),
                density_mean=supported.density.mean(),
                density_min=supported.density.min(),
                density_max=supported.density.max(),
                density_source=source,
                prior_mean=prior,
                own_posterior_mean=x.posterior_mean.mean(),
                affected=x.affected.sum(),
                unaffected=x.unaffected.sum(),
                gnomad_carriers=x.gnomad_carriers.sum(),
                status="supported"
                if len(supported)
                else "observed_without_density"
                if len(x)
                else "no_observed_variant",
            )
        )
    expected = pd.DataFrame(expected)
    for col in ["n_variants", "n_supported_variants", "density_source", "status"]:
        np.testing.assert_array_equal(
            table[col], expected[col], err_msg=f"{gene} {col}"
        )
    errors = {}
    for col in [
        "density_mean",
        "density_min",
        "density_max",
        "prior_mean",
        "own_posterior_mean",
        "affected",
        "unaffected",
        "gnomad_carriers",
    ]:
        np.testing.assert_allclose(
            table[col],
            expected[col],
            atol=1e-10,
            rtol=1e-10,
            equal_nan=True,
            err_msg=f"{gene} {col}",
        )
        errors[col] = float(np.nanmax(abs(table[col] - expected[col])))
    assert (
        not table.loc[
            table.status.ne("supported"), ["density_mean", "density_min", "density_max"]
        ]
        .notna()
        .any()
        .any()
    )
    assert (
        not table.loc[table.status.eq("no_observed_variant"), "own_posterior_mean"]
        .notna()
        .any()
    )
    assert table.n_variants.sum() == len(d)
    assert table.n_supported_variants.sum() == d.density.notna().sum()
    for col in ["affected", "unaffected", "gnomad_carriers"]:
        np.testing.assert_allclose(table[col].sum(), d[col].sum())
    if gene == "HNF1A":
        assert table.loc[
            table.density_source.eq("mixed"), "canonical_pos"
        ].to_list() == [181, 200]
        assert (
            table.loc[table.density_source.eq("mixed"), "n_supported_variants"].sum()
            == 7
        )
    assert b"\r\n" not in table_path.read_bytes()
    assert table_path.stat().st_size < 1_200_000
    files = [fasta, *paths, *density_paths, *compact_paths, prior_path, table_path]
    return dict(
        gene=gene,
        intended_disease=table.intended_disease.iloc[0],
        canonical_length=len(sequence),
        variants=len(d),
        supported_variants=int(d.density.notna().sum()),
        observed_positions=int(d.canonical_pos.nunique()),
        residue_status_counts={
            k: int(v) for k, v in table.status.value_counts().items()
        },
        supported_residue_sources={
            k: int(v)
            for k, v in table.loc[table.status.eq("supported"), "density_source"]
            .value_counts()
            .items()
        },
        residue_mean_of_density=float(table.density_mean.mean()),
        maximum_absolute_errors=errors,
        sha256={
            str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
    )


if __name__ == "__main__":
    results = [audit_gene(gene) for gene in GENES]
    report = dict(
        passed=True,
        genes=results,
        method="Independent complete-residue reconstruction from frozen variant inputs and latest primary density files; plot builder is not imported.",
        aggregation="Equal mean of supported variant-unit LOO density values; min/max is between-variant spread, not a confidence interval. No new fit or exclusions.",
        labels="Intended gene-disease pair labels do not imply the pooled counts have been endpoint-adjudicated or that these scores are calibrated disease risks.",
        hnf1a_mixed="Residues181 and200 retain the internally mixed 3D/polymer classification of all seven frozen variant scores; no source category is averaged twice or reassigned.",
        missingness="Every canonical residue present; observed-but-unsupported and no-observed-variant are distinct. Neither has a fabricated zero density.",
    )
    (HERE / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            [
                {
                    k: v
                    for k, v in r.items()
                    if k not in ["sha256", "maximum_absolute_errors"]
                }
                for r in results
            ],
            indent=2,
        )
    )
