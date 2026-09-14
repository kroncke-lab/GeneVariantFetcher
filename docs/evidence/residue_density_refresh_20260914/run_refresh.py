"""Rebuild five gene/type priors and residue neighborhoods from reviewed counts.

Canonical observed missense donors, variant-only LOO, fixed h=3 positive sigmoid,
independent geometry contexts and same-IDR polymer rules are unchanged. Explicit
count diagnostics are companion estimands, not replacement clinical risks.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
REPO = HERE.parents[2]
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
POP = EVIDENCE / "population_inclusive_penetrance_20260912"
EXTENSION = EVIDENCE / "missense_structural_extension_20260912"
PREVIOUS = EVIDENCE / "brca2_distance_prior_audit_20260913"
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


previous = module("frozen_brca_refresh", PREVIOUS / "recompute_brca2.py")
legacy = previous.legacy
class_fit = previous.class_fit


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = frame.to_csv(index=False, lineterminator="\n").encode()
    if path.suffix == ".gz":
        blob = gzip.compress(blob, mtime=0)
    assert len(blob) < 1_150_000, (path, len(blob))
    path.write_bytes(blob)


def shards(frame, folder, stem, rows=1600):
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.glob(f"{stem}.part*.csv.gz"):
        old.unlink()
    for number, start in enumerate(range(0, len(frame), rows), 1):
        save(
            frame.iloc[start : start + rows], folder / f"{stem}.part{number:03d}.csv.gz"
        )


def geometry_for(gene):
    if gene == "GCK":
        return (
            EVIDENCE
            / "gck_structural_pilot_20260912/geometry/1V4S_canonical_geometry.csv"
        )
    if gene == "BRCA2":
        return (
            EVIDENCE
            / "structural_sanity_20260913/geometry/BRCA2/primary_geometry.csv.gz"
        )
    config = json.loads((EXTENSION / "geometry_config.json").read_text())[gene]
    selected = next(
        x for x in config["frames"] if config["primary"].startswith(x["name"] + "_")
    )
    return EXTENSION / selected["path"]


def original_clinical(gene):
    path = (
        PREVIOUS / "source/BRCA2_literature_identity_corrected.csv.gz"
        if gene == "BRCA2"
        else POP / "audit/literature_identity_flags.csv.gz"
    )
    data = pd.read_csv(path, low_memory=False)
    return data.loc[data.gene.eq(gene)].copy(), path


def build_inputs(gene, population, suffix=""):
    out = HERE / "analysis" / (gene + suffix)
    original, original_path = original_clinical(gene)
    corrected_path = HERE / "source" / gene / f"clinical_input{suffix}.csv.gz"
    if gene == "BRCA2":
        corrected_path = (
            HERE / "source/BRCA2/BRCA2_literature_identity_corrected.csv.gz"
        )
    if gene != "LDLR" and not corrected_path.exists():
        raise FileNotFoundError(
            f"Source-adjudicated clinical input required: {corrected_path}"
        )
    corrected = (
        pd.read_csv(corrected_path, low_memory=False)
        if corrected_path.exists()
        else original.copy()
    )
    assert corrected.gene.eq(gene).all() and corrected.key.is_unique
    legacy.GENES = [gene]
    pop = population.loc[population.gene.eq(gene)].copy()
    old_clinical = original.loc[
        (original.affected_literature + original.unaffected_literature).gt(0)
    ].copy()
    clinical = corrected.loc[
        (corrected.affected_literature + corrected.unaffected_literature).gt(0)
    ].copy()
    before, old_ledger, old_members = legacy.build_union(pop, old_clinical)
    after, new_ledger, new_members = legacy.build_union(pop, clinical)
    assert before.gnomad_carriers.sum() == after.gnomad_carriers.sum()
    crosswalk = old_members.merge(
        new_members,
        on=["gene", "variant_id"],
        suffixes=("_before", "_after"),
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    assert crosswalk._merge.eq("both").all()
    shards(crosswalk, out, "population_membership", rows=7000)
    shards(new_ledger, out, "clinical_join_ledger")
    shards(after, out, "complete_union", rows=7000)
    rows, missense = [], None
    for scenario, universe in [("previous", before), ("refreshed", after)]:
        for kind in ["missense", "nonsense"]:
            units = class_fit.select_type(universe, kind)
            p = legacy.fit_empirical(units)
            row = dict(
                gene=gene,
                input_scope=suffix.removeprefix("_") if suffix else "primary",
                scenario=scenario,
                variant_type=kind,
                variants=len(units),
                affected=units.affected.sum(),
                unaffected_literature=units.unaffected_literature.sum(),
                gnomad_carriers=units.gnomad_carriers.sum(),
                population_only_units=int(units.origin.eq("population_only").sum()),
                zero_affected_units=int(units.affected.eq(0).sum()),
                zero_affected_singletons=int(
                    (units.affected.eq(0) & units.n.eq(1)).sum()
                ),
                **p,
            )
            rows.append(row)
            if scenario == "refreshed":
                post = legacy.posterior_table(units, p)
                shards(post, out, kind + "_posteriors")
                if kind == "missense":
                    missense = post
    priors = pd.DataFrame(rows)
    save(priors, out / "prior_comparison.csv")
    facts = dict(
        gene=gene,
        population_alleles=len(crosswalk),
        population_alleles_preserved_once=True,
        population_alleles_reassigned=int(
            crosswalk.unit_id_before.ne(crosswalk.unit_id_after).sum()
        ),
        gnomad_carriers_all_consequences=int(after.gnomad_carriers.sum()),
        zero_clinical_evidence_keys_excluded=len(corrected) - len(clinical),
        old_retained_clinical_keys=len(old_clinical),
        new_retained_clinical_keys=len(clinical),
        no_zero_observation_donors=bool(after.n.gt(0).all()),
    )
    (out / "union_checks.json").write_text(json.dumps(facts, indent=2) + "\n")
    sources = [original_path] + ([corrected_path] if corrected_path.exists() else [])
    return missense, priors, out, sources


def run_density(gene, variants, priors, out, halves):
    variants = variants.copy().reset_index(drop=True)
    variants["variant_id"] = variants.unit_id
    variants["canonical_variant_id"] = variants.unit_id
    variants["canonical_pos"] = variants.aa_pos.astype(int)
    variants["donor_eligible"] = True
    assert (
        variants.canonical_wt_status.eq("match").all()
        and variants.vclass.eq("missense").all()
    )
    ids = variants.variant_id.tolist()
    geometry_path = geometry_for(gene)
    geometry = pd.read_csv(geometry_path, low_memory=False).fillna({"idr_segment": ""})
    p = priors.loc[
        priors.scenario.eq("refreshed") & priors.variant_type.eq("missense")
    ].iloc[0]
    a, u, n = [variants[c].to_numpy() for c in ["affected", "unaffected", "n"]]
    y = variants.posterior_mean.to_numpy()
    retention = p.strength / (p.strength + n)
    raw_dir = REPO / "results/residue_density_refresh_20260914" / out.name
    raw_dir.mkdir(parents=True, exist_ok=True)
    tables, raw_hashes = [], {}
    for half in halves:
        batches = []
        for start in range(0, len(ids), 128):
            targets = ids[start : start + 128]
            r = empirical_variant_density(
                variants,
                geometry,
                target_ids=targets,
                metric="com",
                half_distance=half,
                include_context_weights=False,
                include_context_model=half == 3,
            )
            weights = r.donor_weights.reindex(index=targets, columns=ids).to_numpy()
            assert not weights[
                np.arange(len(targets)), np.arange(start, start + len(targets))
            ].any()
            supported = weights.sum(axis=1) > 0
            np.testing.assert_allclose(weights.sum(axis=1)[supported], 1, atol=1e-12)
            s = r.summary.copy()
            estimate = weights @ y
            estimate[~supported] = np.nan
            np.testing.assert_allclose(estimate, s.density, equal_nan=True, atol=1e-12)
            s["half_distance"] = half
            s["raw_variant_fraction_density"] = np.where(
                supported, weights @ (a / n), np.nan
            )
            if half == 3:
                assert list(r.context_model.donor_ids) == ids
                diagnostics = dict(
                    neighborhood_prior_retention=weights @ retention,
                    prior_component=weights @ retention * p["mean"],
                    counts_component=weights @ ((1 - retention) * a / n),
                    affected_donor_weight_share=weights @ (a > 0),
                    affected_donor_count=(weights > 0) @ (a > 0).astype(int),
                    conditional_donor_variance=weights**2
                    @ variants.posterior_variance.to_numpy(),
                )
                diagnostics.update(
                    previous.context_count_diagnostics(r.context_model, a, u, p)
                )
                for name, value in diagnostics.items():
                    s[name] = np.where(supported, value, np.nan)
                np.testing.assert_allclose(
                    s.prior_component + s.counts_component, s.density, equal_nan=True
                )
                path = raw_dir / f"weights_{start:05d}_{start + len(targets):05d}.npz"
                save_npz(path, csr_matrix(weights))
                raw_hashes[str(path.relative_to(REPO))] = sha(path)
            identity = variants.iloc[start : start + len(targets)][
                [
                    "gene",
                    "variant_id",
                    "protein_key",
                    "origin",
                    "affected",
                    "unaffected_literature",
                    "gnomad_carriers",
                    "unaffected",
                    "n",
                    "posterior_mean",
                    "posterior_lower_95",
                    "posterior_upper_95",
                ]
            ]
            batches.append(identity.merge(s, on="variant_id", validate="one_to_one"))
            print(
                f"{out.name} h={half:g}: {start + len(targets)}/{len(ids)}", flush=True
            )
        table = pd.concat(batches, ignore_index=True)
        shards(table, out, f"density_h{half:g}")
        tables.append(table)
    all_scenarios = pd.concat(tables, ignore_index=True)
    primary = all_scenarios.loc[all_scenarios.half_distance.eq(3)].copy()
    summary = []
    for half in halves:
        frame = all_scenarios.loc[all_scenarios.half_distance.eq(half)]
        merged = primary[["variant_id", "density"]].merge(
            frame[["variant_id", "density"]],
            on="variant_id",
            suffixes=("_primary", "_alternative"),
            validate="one_to_one",
        )
        delta = (merged.density_alternative - merged.density_primary).abs()
        summary.append(
            dict(
                gene=gene,
                half_distance=half,
                supported=frame.density.notna().sum(),
                median=frame.density.median(),
                minimum=frame.density.min(),
                maximum=frame.density.max(),
                mean_abs_delta_h3=delta.mean(),
                max_abs_delta_h3=delta.max(),
                below_or_equal_0p1_percent=int(frame.density.le(0.001).sum()),
            )
        )
    save(pd.DataFrame(summary), out / "bandwidth_sensitivity.csv")
    sources = [geometry_path]
    # Unchanged genes must exactly reproduce the already-validated primary.
    if (
        gene in {"HNF1A", "LDLR", "KCNQ1"}
        and not (HERE / "source" / gene / "clinical_input.csv.gz").exists()
    ):
        old_paths = sorted(
            (EXTENSION / "analysis" / gene).glob("primary_density.part*.csv.gz")
        )
        old = pd.concat([pd.read_csv(path) for path in old_paths], ignore_index=True)
        common = primary.merge(
            old[["variant_id", "density"]],
            on="variant_id",
            suffixes=("", "_previous"),
            validate="one_to_one",
        )
        assert len(common) == len(primary) == len(old)
        np.testing.assert_allclose(
            common.density, common.density_previous, atol=1e-12, equal_nan=True
        )
        sources += old_paths
    receipt = dict(
        raw_weight_shards=raw_hashes,
        variant_ids_sha256=hashlib.sha256("\n".join(ids).encode()).hexdigest(),
        own_variant_weights_zero=True,
        normalized_weights_sum_to_one=True,
        alpha_adds_affected=True,
        beta_adds_literature_and_gnomad_unaffected=True,
        posterior_decomposition_reproduces=True,
    )
    (out / "density_checks.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return primary, sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", nargs="+", choices=GENES, default=GENES)
    parser.add_argument("--halves", nargs="+", type=float, default=[3, 2, 5])
    parser.add_argument("--gck-proxy", action="store_true")
    parser.add_argument("--sensitivity", choices=["diabetes_proxy", "cardiac_events"])
    args = parser.parse_args()
    assert 3 in args.halves
    legacy.GENES = args.genes
    population, population_sources = legacy.load_population()
    for gene in args.genes:
        sensitivity = "diabetes_proxy" if args.gck_proxy else args.sensitivity
        suffix = "_" + sensitivity if sensitivity else ""
        assert (
            not sensitivity
            or (gene == "GCK" and sensitivity == "diabetes_proxy")
            or (gene == "KCNQ1" and sensitivity == "cardiac_events")
        )
        variants, priors, out, sources = build_inputs(gene, population, suffix)
        _, density_sources = run_density(gene, variants, priors, out, args.halves)
        sources += (
            population_sources
            + density_sources
            + [
                Path(__file__),
                PREVIOUS / "recompute_brca2.py",
                POP / "rebuild_union.py",
                EVIDENCE / "class_matched_penetrance_20260912/fit_class_priors.py",
                REPO.parent
                / "ProteinProximityAnalysis/src/alphafold_rin/empirical_density.py",
                REPO.parent
                / "ProteinProximityAnalysis/src/alphafold_rin/empirical_context.py",
            ]
        )
        hashes = {os.path.relpath(path, EVIDENCE): sha(path) for path in sources}
        (out / "input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main()
