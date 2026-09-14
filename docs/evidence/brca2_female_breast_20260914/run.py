"""Reuse the validated union, prior, density and plotting implementations."""

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from fetch_population import HERE, REPO, module

BASE = HERE.parent / "residue_density_refresh_20260914"
refresh = module("frozen_refresh", BASE / "run_refresh.py")
plot = module("frozen_plot", BASE / "plot_residues.py")
legacy = refresh.legacy
legacy.GENES = ["BRCA2"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    population, sources = legacy.load_population()
    annotated = population.loc[
        population.population_eligible & population.canonical_wt_status.eq("match")
    ].copy()
    annotated["dna"] = annotated.hgvsc.map(legacy.normalized_cdna)
    dna_types = (
        annotated.loc[annotated.dna.ne("")].groupby("dna").vclass.agg(set).to_dict()
    )
    receipt = json.loads((HERE / "population_sex_counts_receipt.json").read_text())
    assert receipt["complete_inventory"] and receipt["alleles"] == 4419
    sexes = pd.read_csv(HERE / "population_sex_counts.csv.gz").set_index(
        "variant_id", verify_integrity=True
    )
    rows, memberships, exclusions = [], [], []
    selected_missense = None
    out = HERE / "analysis/BRCA2_female_breast"
    out.mkdir(parents=True, exist_ok=True)
    sources += [
        HERE / "population_sex_counts.csv.gz",
        HERE / "population_sex_counts_receipt.json",
    ]
    for scope, file in [
        ("source_corrected_all_sex", "clinical_source_corrected_all_sex.csv.gz"),
        ("female_breast", "clinical_female_breast.csv.gz"),
    ]:
        path = HERE / file
        clinical = pd.read_csv(path, low_memory=False)
        sources.append(path)
        for kind, accepted in [
            ("missense", {"missense"}),
            ("nonsense", {"nonsense", "stop_gained"}),
        ]:
            pop = population.loc[
                population.population_eligible
                & population.vclass.isin(accepted)
                & population.canonical_wt_status.eq("match")
            ].copy()
            clin = clinical.loc[clinical.vclass.isin(accepted)].copy()
            # Equal protein stop labels can arise from frameshifts and nonsense
            # SNVs. Restrict population matches to the requested type and reject
            # clinical DNA annotations that explicitly match another consequence.
            for i, row in clin.iterrows():
                cdna = legacy.normalized_cdna(row.cdna)
                dna = dna_types.get(cdna, set())
                # A protein stop/point shorthand does not make a DNA indel a
                # nonsense/missense substitution, even when gnomAD lacks it.
                # Keep unambiguous equal-length delins within one codon (MNVs);
                # other indels require their own consequence adjudication.
                mnv = re.fullmatch(r"c\.(\d+)(?:_(\d+))?del[ACGT]*ins([ACGT]+)", cdna)
                point_mnv = bool(mnv) and (
                    int(mnv[2] or mnv[1]) - int(mnv[1]) + 1 == len(mnv[3])
                    and (int(mnv[1]) - 1) // 3 == (int(mnv[2] or mnv[1]) - 1) // 3
                )
                indel = bool(re.search("del|dup|ins", cdna)) and not point_mnv
                reason = (
                    "clinical_indel_requires_separate_consequence_curation"
                    if indel
                    else "explicit_clinical_DNA_conflicts_with_variant_type"
                    if dna and not dna <= accepted
                    else ""
                )
                if reason:
                    exclusions.append(
                        dict(
                            scope=scope,
                            key=row.key,
                            cdna=row.cdna,
                            requested_type=kind,
                            observed_DNA_types=";".join(sorted(dna)),
                            reason=reason,
                        )
                    )
                    clin.loc[i, "identity_quarantine_recommended"] = True
                    clin.loc[i, "identity_quarantine_reason"] = reason
            # Resolve identity with the full observed DNA inventory of this
            # type, then substitute XX counts before any prior or density fit.
            union, ledger, members = legacy.build_union(pop, clin)
            if scope == "female_breast":
                exact = sexes.reindex(pop.variant_id)
                assert exact.all_carriers.notna().all()
                np.testing.assert_array_equal(exact.all_carriers, pop.gnomad_carriers)
                counts = (
                    members.merge(
                        sexes.reset_index(), on="variant_id", validate="one_to_one"
                    )
                    .groupby("unit_id")[["XX_ac", "XX_homozygote_count", "XX_carriers"]]
                    .sum()
                )
                for field, col in [
                    ("gnomad_ac", "XX_ac"),
                    ("gnomad_homozygotes", "XX_homozygote_count"),
                    ("gnomad_carriers", "XX_carriers"),
                ]:
                    union[field] = union.unit_id.map(counts[col]).fillna(0)
                union["unaffected"] = (
                    union.unaffected_literature + union.gnomad_carriers
                )
                union["n"] = union.affected + union.unaffected
            units = refresh.class_fit.select_type(union, kind)
            zero_n = units.n.eq(0)
            for row in units.loc[zero_n].itertuples():
                exclusions.append(
                    dict(
                        scope=scope,
                        key=row.key,
                        requested_type=kind,
                        reason="no_observed_female_carriers",
                        observed_DNA_types=kind,
                        cdna="",
                    )
                )
            units = units.loc[~zero_n].copy()
            assert units.n.gt(0).all()
            units["population_sex_scope"] = "XX" if scope == "female_breast" else "all"
            units.loc[units.gnomad_carriers.eq(0), "origin"] = "literature_only"
            assert (
                units.gnomad_carriers == units.gnomad_ac - units.gnomad_homozygotes
            ).all()
            prior = legacy.fit_empirical(units)
            rows.append(
                dict(
                    gene="BRCA2",
                    scenario="refreshed",
                    input_scope=scope,
                    variant_type=kind,
                    variants=len(units),
                    affected=units.affected.sum(),
                    unaffected_literature=units.unaffected_literature.sum(),
                    gnomad_carriers=units.gnomad_carriers.sum(),
                    **prior,
                )
            )
            if scope == "female_breast":
                post = legacy.posterior_table(units, prior)
                refresh.shards(post, out, kind + "_posteriors")
                refresh.save(ledger, out / f"{kind}_clinical_join_ledger.csv.gz")
                members = members.merge(
                    sexes.reset_index()[
                        ["variant_id", "all_carriers", "XX_carriers", "XY_carriers"]
                    ],
                    on="variant_id",
                    validate="one_to_one",
                )
                members["variant_type"] = kind
                members["observed_female_unit"] = members.unit_id.isin(units.unit_id)
                assert len(members) == len(pop)
                memberships.append(members)
                if kind == "missense":
                    selected_missense = post
    priors = pd.DataFrame(rows)
    refresh.save(priors, HERE / "prior_comparison.csv")
    members = pd.concat(memberships, ignore_index=True)
    assert members.variant_id.is_unique and len(members) == len(sexes)
    assert set(members.variant_id) == set(sexes.index)
    refresh.save(members, HERE / "population_membership.csv.gz")
    refresh.save(pd.DataFrame(exclusions), HERE / "excluded_units.csv.gz")
    print(
        priors[
            [
                "input_scope",
                "variant_type",
                "variants",
                "affected",
                "unaffected_literature",
                "gnomad_carriers",
                "mean",
                "alpha_empirical",
                "beta_empirical",
            ]
        ].to_string(index=False),
        flush=True,
    )
    if not args.counts_only:
        female_prior = priors.loc[priors.input_scope.eq("female_breast")]
        frame, geometry_sources = refresh.run_density(
            "BRCA2", selected_missense, female_prior, out, [3, 2, 5]
        )
        sources += geometry_sources
        plot.HERE = HERE
        label, fasta = plot.PAIRS["BRCA2"]
        plot.PAIRS["BRCA2"] = ("Female breast cancer — verified clinical subset", fasta)
        (HERE / "plots").mkdir(exist_ok=True)
        current = female_prior.loc[female_prior.variant_type.eq("missense")].iloc[0]
        table = plot.residue_table("BRCA2", frame, HERE.parent / fasta, current)
        refresh.save(table, HERE / "BRCA2_residue_density.csv.gz")
        plot.plot_individual("BRCA2", table, frame)
        sources.append(HERE.parent / fasta)
    sources += [
        Path(__file__).resolve(),
        BASE / "run_refresh.py",
        BASE / "plot_residues.py",
        HERE.parent / "population_inclusive_penetrance_20260912/rebuild_union.py",
        HERE.parent / "class_matched_penetrance_20260912/fit_class_priors.py",
        HERE.parent / "brca2_distance_prior_audit_20260913/recompute_brca2.py",
        REPO.parent / "ProteinProximityAnalysis/src/alphafold_rin/empirical_density.py",
        REPO.parent / "ProteinProximityAnalysis/src/alphafold_rin/empirical_context.py",
    ]
    hashes = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources
    }
    (HERE / "run_input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main()
