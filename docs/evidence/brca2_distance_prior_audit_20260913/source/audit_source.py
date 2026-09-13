#!/usr/bin/env python3
"""Reconstruct BRCA2 source corrections and count diagnostics from frozen evidence.

No network, source DB mutation, model refit, or population-union mutation.
Committed DB/clinical-table snapshots are tied to original hashes in provenance.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
EVIDENCE = REPO / "docs/evidence"
CLASS = EVIDENCE / "class_matched_penetrance_20260912"
POP = EVIDENCE / "population_inclusive_penetrance_20260912"


def save(table, name):
    kwargs = dict(index=False, lineterminator="\n", float_format="%.12g")
    if name.endswith(".gz"):
        kwargs["compression"] = dict(method="gzip", mtime=0)
    table.to_csv(HERE / name, **kwargs)
    assert (HERE / name).stat().st_size < 1_200_000


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    observations = pd.read_csv(HERE / "original_observations.csv.gz")
    provenance = pd.read_csv(HERE / "paper_source_db_records.csv.gz")
    provenance = provenance.set_index("db_variant_id")
    assert provenance.index.is_unique
    catalogue = provenance.loc[
        provenance.key_quotes.fillna("").str.contains("VariationID")
        & provenance.count_provenance.fillna("").str.contains(
            "implicit one carrier per clinical row"
        )
    ]
    selected = observations.loc[
        observations.pmid.eq(40664060) & observations.variant_id.isin(catalogue.index)
    ].copy()
    assert (
        len(selected) == 4070
        and selected.affected.eq(1).all()
        and selected.unaffected.eq(0).all()
    )
    proof = catalogue.loc[selected.variant_id]
    assert (
        proof.affected_count.eq(1).all() and proof.total_carriers_observed.eq(1).all()
    )
    assert proof.unaffected_count.isna().all()
    assert proof.affected_fact_kind.eq("synthesized").all()
    quarantine = pd.read_csv(HERE / "BRCA2_catalog_quarantine_observations.csv")
    pd.testing.assert_frame_equal(
        selected[["key", "pmid", "variant_id"]].reset_index(drop=True),
        quarantine[["original_key", "pmid", "original_observation_variant_id"]].rename(
            columns={
                "original_key": "key",
                "original_observation_variant_id": "variant_id",
            }
        ),
        check_dtype=False,
    )
    np.testing.assert_array_equal(quarantine.source_row, proof.affected_fact_source_row)
    np.testing.assert_array_equal(quarantine.A_to_remove, selected.affected)
    np.testing.assert_array_equal(quarantine.U_to_remove, selected.unaffected)

    # Match the complete source observation identity, never variant ID alone:
    # one DB variant can appear in many PMIDs and those other counts survive.
    excluded = observations.pmid.eq(40664060) & observations.variant_id.isin(
        selected.variant_id
    )
    remaining = observations.loc[~excluded].copy()
    assert len(remaining) == len(observations) - 4070
    paper = pd.read_csv(HERE / "paper_primary_table2_BRCA2.csv")
    assert len(paper) == 15
    clinical = paper.loc[paper.effect.eq("Missense")].copy()
    assert len(clinical) == 5 and clinical.endpoint.eq("BC").all()
    clinical["key"] = clinical.protein.str.removeprefix("p.")
    fasta = POP / "population/BRCA2_canonical.fasta"
    sequence = "".join(fasta.read_text().splitlines()[1:])
    for key in clinical.key:
        ref, pos, _ = re.fullmatch(r"([A-Z])(\d+)([A-Z])", key).groups()
        assert sequence[int(pos) - 1] == ref
    restoration = pd.read_csv(HERE / "BRCA2_clinical_restoration_observations.csv")
    assert set(clinical.key) == set(restoration.protein_key)
    assert restoration.A_to_restore.eq(1).all() and restoration.U_to_restore.eq(0).all()
    clinical_counts = clinical.assign(affected=1.0, unaffected=0.0, pmid=40664060)
    expected_counts = pd.concat(
        [
            remaining[["key", "affected", "unaffected", "pmid"]],
            clinical_counts[["key", "affected", "unaffected", "pmid"]],
        ]
    )
    expected = expected_counts.groupby("key")[["affected", "unaffected"]].sum()
    corrected = pd.read_csv(
        HERE / "BRCA2_literature_identity_corrected.csv.gz"
    ).set_index("key")
    old_flags = pd.read_csv(POP / "audit/literature_identity_flags.csv.gz")
    old_flags = old_flags.loc[old_flags.gene.eq("BRCA2")].set_index("key")
    assert set(corrected.index) == set(old_flags.index) | set(clinical.key)
    np.testing.assert_allclose(
        corrected.affected_literature,
        expected.affected.reindex(corrected.index).fillna(0),
    )
    np.testing.assert_allclose(
        corrected.unaffected_literature,
        expected.unaffected.reindex(corrected.index).fillna(0),
    )
    np.testing.assert_allclose(
        corrected.n_literature,
        corrected.affected_literature + corrected.unaffected_literature,
    )
    np.testing.assert_array_equal(
        corrected.drop_clinical_key_zero_clinical_evidence, corrected.n_literature.eq(0)
    )
    assert len(corrected) == 6110 and corrected.n_literature.eq(0).sum() == 3644
    actual_pmids = expected_counts.groupby("key").pmid.agg(
        lambda s: ";".join(map(str, sorted(set(s))))
    )
    assert (
        corrected.pmids.fillna("")
        .eq(actual_pmids.reindex(corrected.index).fillna(""))
        .all()
    )

    paths = sorted((CLASS / "analysis/empirical_posteriors").glob("BRCA2.part*.csv.gz"))
    full = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    assert full.unit_id.is_unique and full.gene.eq("BRCA2").all()
    assert full.canonical_wt_status.eq("match").all()
    assert set(full.variant_type) == {"missense", "nonsense"}
    assert full.n.gt(0).all()
    for actual, target in [
        (full.unaffected, full.unaffected_literature + full.gnomad_carriers),
        (full.posterior_alpha, full.alpha_empirical + full.affected),
        (full.posterior_beta, full.beta_empirical + full.unaffected),
        (
            full.posterior_mean,
            full.posterior_alpha / (full.posterior_alpha + full.posterior_beta),
        ),
    ]:
        np.testing.assert_allclose(actual, target, atol=2e-12)
    d = full.loc[full.variant_type.eq("missense")].copy()
    assert len(d) == 6656
    d["fit_weight"] = 1 - 1 / (d.n + 0.01)
    d["prior_numerator"] = d.fit_weight * d.affected / d.n
    d["catalog_A_removed"] = d.literature_key.map(
        selected.groupby("key").affected.sum()
    ).fillna(0)
    d["count_pattern"] = np.select(
        [
            d.affected.eq(0) & d.unaffected.eq(1),
            d.affected.eq(0) & d.unaffected.gt(1),
            d.affected.eq(1) & d.unaffected.eq(0),
            d.affected.eq(1) & d.unaffected.gt(0),
        ],
        ["A0_U1", "A0_Ugt1", "A1_U0", "A1_Upositive"],
        default="Agt1",
    )
    groups = []
    for label in ["origin", "count_pattern"]:
        for value, x in d.groupby(label):
            groups.append(
                dict(
                    grouping=label,
                    group=value,
                    units=len(x),
                    affected=x.affected.sum(),
                    unaffected_literature=x.unaffected_literature.sum(),
                    gnomad_carriers=x.gnomad_carriers.sum(),
                    historical_fit_weight=x.fit_weight.sum(),
                    fit_weight_share=x.fit_weight.sum() / d.fit_weight.sum(),
                    prior_mean_additive_component=x.prior_numerator.sum()
                    / d.fit_weight.sum(),
                    catalog_affected_removed=x.catalog_A_removed.sum(),
                )
            )
    save(pd.DataFrame(groups), "missense_count_weight_groups.csv")
    save(
        d.nlargest(20, "gnomad_carriers")[
            [
                "unit_id",
                "protein_key",
                "affected",
                "unaffected_literature",
                "gnomad_carriers",
                "n",
                "fit_weight",
                "posterior_mean",
            ]
        ],
        "largest_population_denominators.csv",
    )
    keys = set(d.literature_key.dropna())
    obs = observations.loc[observations.key.isin(keys)].copy()
    obs["A1U0"] = obs.affected.eq(1) & obs.unaffected.eq(0)
    by_pmid = (
        obs.groupby("pmid")
        .agg(
            observations=("key", "size"),
            unique_keys=("key", "nunique"),
            affected=("affected", "sum"),
            unaffected=("unaffected", "sum"),
            A1U0_rows=("A1U0", "sum"),
        )
        .reset_index()
        .sort_values("observations", ascending=False)
    )
    save(by_pmid, "missense_observations_by_pmid.csv")
    queue = obs.loc[obs.pmid.isin([36385461, 33054725])].copy()
    queue["review_status"] = "pending_individual_source_adjudication_no_count_change"
    queue["review_reason"] = queue.pmid.map(
        {
            36385461: "Variant inventory plus implicit carrier inference; inspect each selected source row and clinical evidence elsewhere; possible cohort overlap.",
            33054725: "Actual tumor sample IDs but germline/somatic status, mixed cancer endpoint and cohort overlap require review.",
        }
    )
    queue["primary_url"] = queue.pmid.map(
        {
            36385461: "https://pmc.ncbi.nlm.nih.gov/articles/PMC10098510/",
            33054725: "https://pmc.ncbi.nlm.nih.gov/articles/PMC7556962/",
        }
    )
    assert len(queue) == 227 and queue.affected.sum() == 227
    save(queue, "remaining_BRCA2_source_queue.csv.gz")
    save(
        d.loc[
            d.affected.eq(0) & d.gnomad_carriers.gt(0),
            [
                "unit_id",
                "protein_key",
                "origin",
                "affected",
                "unaffected_literature",
                "gnomad_carriers",
                "n",
            ],
        ],
        "zero_affected_population_units.csv.gz",
    )
    comparison = pd.read_csv(CLASS / "analysis/empirical_prior_comparison.csv")
    prior = comparison.loc[
        comparison.gene.eq("BRCA2") & comparison.scope.eq("canonical_missense")
    ].iloc[0]
    mean = d.prior_numerator.sum() / d.fit_weight.sum()
    mse = np.sum(d.fit_weight * (d.affected / d.n - mean) ** 2) / len(d)
    strength = mean * (1 - mean) / mse - 1
    np.testing.assert_allclose(
        [mean, strength], [prior["mean"], prior.strength], atol=2e-12
    )
    summaries = dict(
        passed=True,
        original_missense_units=len(d),
        original_A=float(d.affected.sum()),
        original_U_literature=float(d.unaffected_literature.sum()),
        original_U_gnomad=int(d.gnomad_carriers.sum()),
        A1U0_units=int((d.affected.eq(1) & d.unaffected.eq(0)).sum()),
        A1U0_units_from_catalog_only=int(
            (d.affected.eq(1) & d.unaffected.eq(0) & d.catalog_A_removed.eq(1)).sum()
        ),
        A0U1_units=int((d.affected.eq(0) & d.unaffected.eq(1)).sum()),
        A0_gnomad_observed_units=int(
            (d.affected.eq(0) & d.gnomad_carriers.gt(0)).sum()
        ),
        A0_gnomad_carriers=int(d.loc[d.affected.eq(0), "gnomad_carriers"].sum()),
        top10_gnomad_count_share=float(
            d.nlargest(10, "gnomad_carriers").gnomad_carriers.sum()
            / d.gnomad_carriers.sum()
        ),
        top10_gnomad_prior_fit_weight_share=float(
            d.nlargest(10, "gnomad_carriers").fit_weight.sum() / d.fit_weight.sum()
        ),
        frozen_prior_mean=float(mean),
        frozen_prior_strength=float(strength),
        all_type_catalog_observations_quarantined=len(selected),
        canonical_missense_catalog_A_removed=float(d.catalog_A_removed.sum()),
        genuine_missense_clinical_rows_restored=5,
        corrected_clinical_keys=len(corrected),
        corrected_zero_clinical_evidence_keys=int(corrected.n_literature.eq(0).sum()),
        pending_truncated_header_observations=1,
        retained_missense_observations_in_other_source_queues=227,
        union_rule="Drop clinical keys with corrected A+Ulit=0 before rebuild_union. Retain their observed genomic population alleles and all gnomAD counts. Rejoin five genuine case observations; never retain orphaned clinical protein aggregates as extra population units.",
        orientation="No A/U swap or lost gnomAD denominator; alpha adds A, beta adds Ulit+Ugn.",
        source_error="Implicit clinical-row carrier inference was applied to a ClinVar variant catalog; this is not a validated human carrier count or assay subject.",
        scope_limit="Only the proven PMID40664060 catalog block is adjudicated. Other papers and the pending out-of-scope frameshift may require further curation; corrected counts are not fully endpoint/population calibrated.",
        hashes={
            str(p.relative_to(REPO)): sha(p)
            for p in [
                *paths,
                HERE / "original_observations.csv.gz",
                HERE / "paper_source_db_records.csv.gz",
                HERE / "paper_primary_table2_BRCA2.csv",
                HERE / "BRCA2_literature_identity_corrected.csv.gz",
            ]
        },
    )
    (HERE / "checks.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summaries.items() if k != "hashes"}, indent=2))


if __name__ == "__main__":
    main()
