"""Rebuild BRCA2 from adjudicated clinical keys; keep geometry and h=3 fixed.

The posterior-neighborhood remains the primary historical estimand. Raw variant
fractions and absolute-kernel count pooling are explicitly separate diagnostics.
No disease-risk calibration or additional outcome-model fitting is performed.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
REPO = HERE.parents[2]
POP = EVIDENCE / "population_inclusive_penetrance_20260912"
OLD = EVIDENCE / "structural_sanity_20260913"
OUT = HERE / "corrected"
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy = load_module("frozen_union", POP / "rebuild_union.py")
legacy.GENES = ["BRCA2"]
class_fit = load_module(
    "frozen_class", EVIDENCE / "class_matched_penetrance_20260912/fit_class_priors.py"
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, name):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = frame.to_csv(index=False, lineterminator="\n").encode()
    if path.suffix == ".gz":
        blob = gzip.compress(blob, mtime=0)
    assert len(blob) < 1_150_000, (name, len(blob))
    path.write_bytes(blob)


def shards(frame, stem, size=2000):
    for number, start in enumerate(range(0, len(frame), size), 1):
        save(frame.iloc[start : start + size], f"{stem}.part{number:03d}.csv.gz")


def rebuild():
    population, sources = legacy.load_population()
    literature_path = POP / "audit/literature_identity_flags.csv.gz"
    original = pd.read_csv(literature_path, low_memory=False)
    original = original.loc[original.gene.eq("BRCA2")].copy()
    before, old_ledger, old_members = legacy.build_union(population, original)
    frozen_paths = sorted(
        (
            EVIDENCE / "class_matched_penetrance_20260912/analysis/empirical_posteriors"
        ).glob("BRCA2.part*.csv.gz")
    )
    frozen = pd.concat([pd.read_csv(p) for p in frozen_paths], ignore_index=True)
    compare_cols = ["affected", "unaffected_literature", "gnomad_carriers", "n"]
    for kind in ["missense", "nonsense"]:
        actual = class_fit.select_type(before, kind).set_index("unit_id").sort_index()
        expected = frozen.loc[frozen.variant_type.eq(kind)].set_index("unit_id")
        expected = expected.sort_index()
        assert actual.index.equals(expected.index)
        np.testing.assert_allclose(actual[compare_cols], expected[compare_cols])

    removal_path = HERE / "source/BRCA2_catalog_quarantine_observations.csv"
    restoration_path = HERE / "source/BRCA2_clinical_restoration_observations.csv"
    removal = pd.read_csv(removal_path)
    restoration = pd.read_csv(restoration_path)
    deltas = removal.groupby("original_key")[["A_to_remove", "U_to_remove"]].sum()
    missing = deltas.index.difference(original.key)
    # Retained observation evidence and clinical eligibility are different stages.
    save(
        deltas.loc[missing].reset_index(), "quarantine_keys_outside_clinical_input.csv"
    )
    clinical = original.set_index("key", drop=False).copy()
    common = clinical.index.intersection(deltas.index)
    clinical.loc[common, "affected_literature"] -= deltas.loc[common, "A_to_remove"]
    clinical.loc[common, "unaffected_literature"] -= deltas.loc[common, "U_to_remove"]
    assert (clinical[["affected_literature", "unaffected_literature"]] >= 0).all().all()
    template = original.loc[original.key.eq("V2076I")].iloc[0]
    for row in restoration.itertuples(index=False):
        key = row.protein_key
        if key not in clinical.index:
            record = template.to_dict()
            record.update(
                key=key,
                aa_pos=int(key[1:-1]),
                aa_ref=key[0],
                aa_alt=key[-1],
                cdna=row.source_cdna,
                canonical_wt=key[0],
                canonical_point_key=key,
                archived_vf_variant_ids="",
                archived_n_alleles=0,
                archived_allele_match="source_adjudicated_clinical_variant",
                affected_literature=0,
                unaffected_literature=0,
                pmids=str(row.pmid),
            )
            clinical.loc[key] = record
        clinical.loc[key, "affected_literature"] += row.A_to_restore
        clinical.loc[key, "unaffected_literature"] += row.U_to_restore
    clinical["n_literature"] = (
        clinical.affected_literature + clinical.unaffected_literature
    )
    # Independently calculated deltas must match the source auditor's complete
    # key ledger; then use that ledger's corrected PMID/identity provenance.
    corrected_path = HERE / "source/BRCA2_literature_identity_corrected.csv.gz"
    adjudicated = pd.read_csv(corrected_path, low_memory=False).set_index(
        "key", drop=False
    )
    assert clinical.index.sort_values().equals(adjudicated.index.sort_values())
    for field in ["affected_literature", "unaffected_literature", "n_literature"]:
        np.testing.assert_allclose(
            clinical.sort_index()[field], adjudicated.sort_index()[field]
        )
    clinical = adjudicated.copy()
    clinical["count_source_correction"] = (
        "PMID40664060 ClinVar catalog quarantined; Table 2 patient rows restored"
    )
    save(
        clinical.reset_index(drop=True), "clinical_key_counts_before_zero_filter.csv.gz"
    )
    removed = clinical.loc[clinical.n_literature.eq(0)].reset_index(drop=True)
    save(removed, "zero_clinical_evidence_keys.csv.gz")
    clinical = clinical.loc[clinical.n_literature.gt(0)].reset_index(drop=True)
    after, new_ledger, new_members = legacy.build_union(population, clinical)
    membership = old_members.merge(
        new_members,
        on=["gene", "variant_id"],
        suffixes=("_before", "_after"),
        validate="one_to_one",
        how="outer",
        indicator=True,
    )
    assert membership._merge.eq("both").all()
    assert before.gnomad_carriers.sum() == after.gnomad_carriers.sum()
    assert (after.n > 0).all()
    shards(membership, "population_membership_crosswalk", size=8000)
    shards(new_ledger, "clinical_join_ledger")
    shards(after, "complete_union", size=8000)
    records, posterior = [], {}
    for scenario, universe in [("before", before), ("corrected", after)]:
        for kind in ["missense", "nonsense"]:
            units = class_fit.select_type(universe, kind)
            parameters = legacy.fit_empirical(units)
            records.append(
                dict(
                    scenario=scenario,
                    variant_type=kind,
                    variants=len(units),
                    affected=float(units.affected.sum()),
                    unaffected_literature=float(units.unaffected_literature.sum()),
                    gnomad_carriers=float(units.gnomad_carriers.sum()),
                    population_only_units=int(units.origin.eq("population_only").sum()),
                    no_affected_units=int(units.affected.eq(0).sum()),
                    unaffected_singletons=int(
                        (units.affected.eq(0) & units.n.eq(1)).sum()
                    ),
                    **parameters,
                )
            )
            if scenario == "corrected":
                posterior[kind] = legacy.posterior_table(units, parameters)
                shards(posterior[kind], f"{kind}_posteriors")
    priors = pd.DataFrame(records)
    save(priors, "prior_comparison.csv")
    facts = {
        "uncorrected_union_matches_frozen_classes": True,
        "catalog_observations_quarantined": len(removal),
        "catalog_keys_absent_from_prior_clinical_input": len(missing),
        "clinical_observations_restored": len(restoration),
        "zero_evidence_clinical_keys_removed": len(removed),
        "population_alleles_preserved": len(membership),
        "population_alleles_reassigned": int(
            membership.unit_id_before.ne(membership.unit_id_after).sum()
        ),
        "gnomad_carriers_all_consequences_before_and_after": int(
            after.gnomad_carriers.sum()
        ),
        "old_clinical_keys": len(original),
        "corrected_clinical_keys": len(clinical),
        "zero_observation_donors": 0,
    }
    (OUT / "union_checks.json").write_text(json.dumps(facts, indent=2) + "\n")
    sources += [
        literature_path,
        removal_path,
        restoration_path,
        corrected_path,
        *frozen_paths,
    ]
    return posterior["missense"], priors, sources


def context_count_diagnostics(model, affected, unaffected, prior):
    valid = np.isfinite(model.context_log_sums) & np.isfinite(model.context_densities)
    rows = model.context_target_rows[valid]
    geometry_rows = model.context_geometry_rows[valid]
    raw = np.exp(
        model.position_log_kernels[
            geometry_rows[:, None], model.donor_position_indices[None, :]
        ]
    )
    ids = np.asarray(model.donor_ids)
    raw[model.target_canonical_ids[rows, None] == ids[None, :]] = 0
    a, u = raw @ affected, raw @ unaffected
    kappa, alpha = prior.strength, prior.alpha_empirical
    values = {
        "raw_kernel_affected": a,
        "raw_kernel_unaffected": u,
        "raw_kernel_pooled_fraction": a / (a + u),
        "raw_kernel_one_prior_posterior": (alpha + a) / (kappa + a + u),
        "raw_kernel_total_weight": raw.sum(axis=1),
    }
    counts = np.bincount(rows, minlength=len(model.target_ids))
    result = {}
    for key, value in values.items():
        sums = np.bincount(rows, weights=value, minlength=len(model.target_ids))
        result[key] = np.divide(
            sums, counts, out=np.full(len(counts), np.nan), where=counts > 0
        )
    return result


def density(variants, priors):
    variants = variants.copy().reset_index(drop=True)
    variants["variant_id"] = variants.unit_id
    variants["canonical_variant_id"] = variants.unit_id
    variants["canonical_pos"] = variants.aa_pos.astype(int)
    variants["donor_eligible"] = True
    ids = variants.variant_id.tolist()
    p = priors.loc[
        priors.scenario.eq("corrected") & priors.variant_type.eq("missense")
    ].iloc[0]
    old_prior = priors.loc[
        priors.scenario.eq("before") & priors.variant_type.eq("missense")
    ].iloc[0]
    geometry_path = OLD / "geometry/BRCA2/primary_geometry.csv.gz"
    geometry = pd.read_csv(geometry_path, low_memory=False).fillna({"idr_segment": ""})
    a, u, n = (variants[c].to_numpy() for c in ["affected", "unaffected", "n"])
    y = variants.posterior_mean.to_numpy()
    fraction = a / n
    retention = p.strength / (p.strength + n)
    old_prior_new_counts = (old_prior.alpha_empirical + a) / (old_prior.strength + n)
    batches = []
    for start in range(0, len(ids), 128):
        targets = ids[start : start + 128]
        result = empirical_variant_density(
            variants,
            geometry,
            target_ids=targets,
            metric="com",
            half_distance=3,
            include_context_weights=False,
            include_context_model=True,
        )
        weights = result.donor_weights.reindex(index=targets, columns=ids).to_numpy()
        assert list(result.context_model.donor_ids) == ids
        assert not weights[
            np.arange(len(targets)), np.arange(start, start + len(targets))
        ].any()
        supported = weights.sum(axis=1) > 0
        np.testing.assert_allclose(weights.sum(axis=1)[supported], 1, atol=1e-12)
        s = result.summary.copy()
        predicted = weights @ y
        predicted[~supported] = np.nan
        np.testing.assert_allclose(predicted, s.density, equal_nan=True, atol=1e-12)
        diagnostics = dict(
            raw_variant_fraction_density=weights @ fraction,
            neighborhood_prior_retention=weights @ retention,
            prior_component=weights @ retention * p["mean"],
            counts_component=weights @ ((1 - retention) * fraction),
            old_prior_corrected_counts_density=weights @ old_prior_new_counts,
            affected_donor_weight_share=weights @ (a > 0),
            affected_donor_count=(weights > 0) @ (a > 0).astype(int),
        )
        diagnostics.update(context_count_diagnostics(result.context_model, a, u, p))
        for key, value in diagnostics.items():
            s[key] = np.where(supported, value, np.nan)
        np.testing.assert_allclose(
            s.prior_component + s.counts_component, s.density, equal_nan=True
        )
        identity = variants.iloc[start : start + len(targets)][
            [
                "variant_id",
                "protein_key",
                "origin",
                "affected",
                "unaffected_literature",
                "gnomad_carriers",
                "unaffected",
                "n",
                "posterior_mean",
            ]
        ]
        batches.append(identity.merge(s, on="variant_id", validate="one_to_one"))
        print(f"corrected h3 density: {start + len(targets)}/{len(ids)}", flush=True)
    result = pd.concat(batches, ignore_index=True)
    shards(result, "primary_density", size=1600)
    segments = geometry.loc[
        geometry.geometry_state.eq("idr"), ["canonical_pos", "idr_segment"]
    ].drop_duplicates()
    assert segments.canonical_pos.is_unique
    segment_units = result.merge(
        segments, on="canonical_pos", how="inner", validate="many_to_one"
    )
    segment_counts = (
        segment_units.groupby("idr_segment")
        .agg(
            variants=("variant_id", "size"),
            affected=("affected", "sum"),
            unaffected=("unaffected", "sum"),
            supported=("density", "count"),
            mean_density=("density", "mean"),
            min_density=("density", "min"),
            max_density=("density", "max"),
            mean_raw_variant_fraction_density=("raw_variant_fraction_density", "mean"),
            mean_raw_kernel_one_prior_posterior=(
                "raw_kernel_one_prior_posterior",
                "mean",
            ),
        )
        .reset_index()
    )
    save(segment_counts, "idr_segment_counts.csv")
    no_case = result.loc[result.affected_donor_count.eq(0)]
    save(no_case, "targets_without_affected_donors.csv.gz")
    old_paths = sorted((OLD / "analysis/BRCA2").glob("primary_density.part*.csv.gz"))
    old = pd.concat([pd.read_csv(path) for path in old_paths], ignore_index=True)
    metrics = []
    for name, data in [("before", old), ("corrected", result)]:
        columns = ["density"]
        if name == "corrected":
            columns += [
                "old_prior_corrected_counts_density",
                "raw_variant_fraction_density",
                "raw_kernel_pooled_fraction",
                "raw_kernel_one_prior_posterior",
            ]
        for column in columns:
            value = data[column].dropna()
            metrics.append(
                dict(
                    scenario=name,
                    score=column,
                    variants=len(data),
                    supported=len(value),
                    mean=value.mean(),
                    median=value.median(),
                    minimum=value.min(),
                    maximum=value.max(),
                    at_or_below_0p1_percent=int(value.le(0.001).sum()),
                    at_or_below_1_percent=int(value.le(0.01).sum()),
                    exact_zero=int(value.eq(0).sum()),
                )
            )
    save(pd.DataFrame(metrics), "density_comparison.csv")
    return result, old, geometry_path, old_paths


def plots(corrected, old, priors):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42}
    )
    c = corrected.groupby("canonical_pos").agg(
        density=("density", "mean"),
        raw_variant_fraction_density=("raw_variant_fraction_density", "mean"),
        raw_kernel_pooled_fraction=("raw_kernel_pooled_fraction", "mean"),
        raw_kernel_one_prior_posterior=("raw_kernel_one_prior_posterior", "mean"),
        old_prior_corrected_counts_density=(
            "old_prior_corrected_counts_density",
            "mean",
        ),
        affected=("affected", "sum"),
        unaffected=("unaffected", "sum"),
        supported_variants=("density", "count"),
        variants=("variant_id", "size"),
    )
    before = old.groupby("canonical_pos").density.mean()
    c = c.reindex(range(1, 3419))
    c["before_density"] = before.reindex(c.index)
    save(c.rename_axis("canonical_pos").reset_index(), "BRCA2_residue_density.csv")
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(13, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [1.25, 1.1, 0.8]},
    )
    colors = {
        "old": "#92979E",
        "corrected": "#135DAD",
        "raw": "#BB571C",
        "pooled": "#7E3C96",
    }
    axes[0].plot(
        c.index,
        100 * c.before_density,
        ".",
        ms=2.4,
        color=colors["old"],
        alpha=0.6,
        label="Previous posterior neighborhood",
    )
    axes[0].plot(
        c.index,
        100 * c.density,
        ".",
        ms=2.8,
        color=colors["corrected"],
        label="Corrected posterior neighborhood",
    )
    axes[0].plot(
        c.index,
        100 * c.raw_variant_fraction_density,
        ".",
        ms=2.2,
        color=colors["raw"],
        alpha=0.7,
        label="Corrected mean of observed variant fractions",
    )
    axes[0].set_ylabel("Neighborhood score (%)")
    axes[0].set_ylim(bottom=0)
    axes[0].legend(loc="upper left", fontsize=9, frameon=False, ncol=1)
    floor = 0.001  # Percent: 0.00001 fraction. Zero markers are explicitly labeled.
    for column, label, color in [
        ("density", "Posterior neighborhood", colors["corrected"]),
        (
            "raw_variant_fraction_density",
            "Mean observed variant fraction",
            colors["raw"],
        ),
        (
            "raw_kernel_one_prior_posterior",
            "Raw-kernel pooled counts + one prior / context",
            colors["pooled"],
        ),
    ]:
        value = c[column] * 100
        positive = value.gt(0)
        axes[1].plot(
            c.index[positive], value[positive], ".", ms=2.7, color=color, label=label
        )
        zero = value.eq(0)
        axes[1].scatter(
            c.index[zero], np.full(zero.sum(), floor), marker="v", s=16, color=color
        )
    axes[1].set_yscale("log")
    axes[1].set_ylim(floor * 0.75, 100)
    axes[1].axhline(0.1, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].text(3430, 0.1, "0.1%", va="center", fontsize=9)
    axes[1].set_ylabel("Score (%) — logarithmic scale")
    axes[1].legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=3,
        fontsize=8,
        frameon=False,
    )
    axes[1].text(
        0.005,
        0.015,
        "▼ at 0.001% means exactly zero, not a fitted floor",
        transform=axes[1].transAxes,
        fontsize=8,
    )
    for column, color, label in [
        ("affected", "#A62930", "Affected observations"),
        (
            "unaffected",
            "#387D48",
            "Unaffected observations (gnomAD assumed unaffected)",
        ),
    ]:
        value = c[column]
        valid = value.gt(0)
        axes[2].vlines(
            c.index[valid], 1, value[valid], color=color, alpha=0.35, linewidth=0.6
        )
        axes[2].scatter(c.index[valid], value[valid], s=4, color=color, label=label)
    axes[2].set_yscale("log")
    axes[2].set_ylim(0.7, None)
    axes[2].set_ylabel("Observed counts / residue")
    axes[2].legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        fontsize=8,
        frameon=False,
    )
    axes[2].set_xlabel("BRCA2 canonical residue (P51587; 1–3418)")
    for axis in axes:
        axis.set_xlim(1, 3418)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.15)
    fig.suptitle(
        "BRCA2 • corrected source counts, unchanged distance weighting",
        x=0.075,
        ha="left",
        fontsize=16,
        weight="bold",
    )
    fig.text(
        0.075,
        0.935,
        "Missense only • variant-only leave-one-out • h = 3 Å • same-segment polymer for candidate IDRs",
        fontsize=10,
    )
    fig.text(
        0.075,
        0.025,
        "One point = mean of supported variant scores at a residue. Gaps are unestimated; points are not connected across gaps.\nThese are descriptive neighborhood features, not calibrated individual cancer risks. Count-pooling diagnostic changes the estimand.",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.105, hspace=0.23)
    fig.savefig(OUT / "BRCA2_CORRECTED_RESIDUE_DENSITY.png", dpi=160)
    fig.savefig(OUT / "BRCA2_CORRECTED_RESIDUE_DENSITY.pdf")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    variants, priors, sources = rebuild()
    result, old, geometry_path, old_paths = density(variants, priors)
    plots(result, old, priors)
    sources += [
        geometry_path,
        *old_paths,
        Path(__file__),
        POP / "rebuild_union.py",
        EVIDENCE / "class_matched_penetrance_20260912/fit_class_priors.py",
    ]
    receipt = {str(path.relative_to(EVIDENCE)): sha(path) for path in sources}
    (OUT / "input_hashes.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(priors.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
