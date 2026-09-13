"""Independent corrected-union, empirical moment, and context-count audit."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
EVIDENCE = AUDIT.parent
REPO = EVIDENCE.parents[1]
OUT = AUDIT / "corrected"
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density

spec = importlib.util.spec_from_file_location(
    "corrected_inputs", AUDIT / "recompute_brca2.py"
)
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def read(stem):
    return pd.concat(
        [
            pd.read_csv(p, low_memory=False)
            for p in sorted(OUT.glob(f"{stem}.part*.csv.gz"))
        ],
        ignore_index=True,
    )


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    population, _ = source.legacy.load_population()
    pop = population.loc[population.population_eligible].set_index("variant_id")
    union, crosswalk = read("complete_union"), read("population_membership_crosswalk")
    priors, primary = pd.read_csv(OUT / "prior_comparison.csv"), read("primary_density")
    frozen = pd.concat(
        [
            pd.read_csv(p, low_memory=False)
            for p in sorted(
                (
                    EVIDENCE
                    / "class_matched_penetrance_20260912/analysis/empirical_posteriors"
                ).glob("BRCA2.part*.csv.gz")
            )
        ],
        ignore_index=True,
    )
    assert union.unit_id.is_unique
    assert crosswalk.variant_id.is_unique and crosswalk._merge.eq("both").all()
    assert set(crosswalk.variant_id) == set(pop.index)
    assert (union.n > 0).all()
    np.testing.assert_allclose(
        union.unaffected, union.unaffected_literature + union.gnomad_carriers
    )
    np.testing.assert_allclose(union.n, union.affected + union.unaffected)
    exploded = union.loc[
        union.member_alleles.notna(), ["unit_id", "member_alleles"]
    ].copy()
    exploded["variant_id"] = exploded.member_alleles.str.split(";")
    exploded = exploded.explode("variant_id").drop(columns="member_alleles")
    assert exploded.variant_id.is_unique
    assert set(exploded.variant_id) == set(pop.index)
    members = exploded.set_index("variant_id").reindex(crosswalk.variant_id)
    assert members.unit_id.tolist() == crosswalk.unit_id_after.tolist()
    counts = exploded.merge(
        pop[["gnomad_carriers"]],
        left_on="variant_id",
        right_index=True,
        validate="one_to_one",
    )
    expected = (
        counts.groupby("unit_id")
        .gnomad_carriers.sum()
        .reindex(union.unit_id, fill_value=0)
    )
    np.testing.assert_array_equal(expected, union.gnomad_carriers)
    assert union.gnomad_carriers.sum() == pop.gnomad_carriers.sum()

    original = pd.read_csv(
        EVIDENCE
        / "population_inclusive_penetrance_20260912/audit/literature_identity_flags.csv.gz",
        low_memory=False,
    )
    original = original.loc[original.gene.eq("BRCA2")].set_index("key")
    removal = pd.read_csv(AUDIT / "source/BRCA2_catalog_quarantine_observations.csv")
    restoration = pd.read_csv(
        AUDIT / "source/BRCA2_clinical_restoration_observations.csv"
    )
    adjudicated = pd.read_csv(
        AUDIT / "source/BRCA2_literature_identity_corrected.csv.gz", low_memory=False
    ).set_index("key")
    corrected_keys = pd.read_csv(
        OUT / "clinical_key_counts_before_zero_filter.csv.gz", low_memory=False
    ).set_index("key")
    assert corrected_keys.index.equals(adjudicated.index)
    expected_keys = original[["affected_literature", "unaffected_literature"]].copy()
    for key, group in removal.groupby("original_key"):
        expected_keys.loc[key] -= group[["A_to_remove", "U_to_remove"]].sum().to_numpy()
    for row in restoration.itertuples():
        if row.protein_key not in expected_keys.index:
            expected_keys.loc[row.protein_key] = [0, 0]
        expected_keys.loc[row.protein_key] += [row.A_to_restore, row.U_to_restore]
    assert set(expected_keys.index) == set(adjudicated.index)
    for field in ["affected_literature", "unaffected_literature"]:
        np.testing.assert_array_equal(
            expected_keys[field].reindex(adjudicated.index), adjudicated[field]
        )
        np.testing.assert_array_equal(adjudicated[field], corrected_keys[field])
    np.testing.assert_allclose(
        corrected_keys.n_literature,
        corrected_keys.affected_literature + corrected_keys.unaffected_literature,
    )
    zero_keys = set(corrected_keys.index[corrected_keys.n_literature.eq(0)])
    assert not set(union.literature_key.dropna()).intersection(zero_keys)
    changed = crosswalk.loc[crosswalk.unit_id_before.ne(crosswalk.unit_id_after)]
    released = changed.loc[changed.unit_id_after.str.startswith("BRCA2|g:")]
    assert released.unit_id_before.isin(["BRCA2|lit:" + key for key in zero_keys]).all()
    assert released.unit_id_after.eq("BRCA2|g:" + released.variant_id).all()
    captured = changed.loc[~changed.index.isin(released.index)]
    assert captured.unit_id_after.isin("BRCA2|lit:" + restoration.protein_key).all()

    old_nonsense = frozen.loc[frozen.variant_type.eq("nonsense")].set_index("unit_id")
    new_nonsense = read("nonsense_posteriors")
    class_changes = crosswalk.loc[
        crosswalk.unit_id_before.isin(old_nonsense.index)
        & ~crosswalk.unit_id_after.isin(new_nonsense.unit_id)
    ].copy()
    after = union.set_index("unit_id")
    class_rows = []
    for row in class_changes.itertuples():
        old, new, allele = (
            old_nonsense.loc[row.unit_id_before],
            after.loc[row.unit_id_after],
            pop.loc[row.variant_id],
        )
        delta = len(allele.alt) - len(allele.ref)
        assert new.vclass == allele.vclass == "frameshift" and delta % 3 != 0
        class_rows.append(
            {
                "variant_id": row.variant_id,
                "unit_id_before": row.unit_id_before,
                "unit_id_after": row.unit_id_after,
                "old_class": old.vclass,
                "new_class": new.vclass,
                "protein_key": new.protein_key,
                "ref": allele.ref,
                "alt": allele.alt,
                "net_length_change": delta,
                "source_hgvsc": allele.hgvsc,
                "source_hgvsp": allele.hgvsp,
                "source_transcript": allele.transcript_id,
                "canonical_transcript": allele.canonical_transcript_id,
                "gnomad_carriers": int(new.gnomad_carriers),
            }
        )
    classification = pd.DataFrame(class_rows)
    assert len(classification) == 4 and classification.gnomad_carriers.sum() == 511
    classification.to_csv(HERE / "released_frameshift_alleles.csv", index=False)

    orientation = []
    for kind in ["missense", "nonsense"]:
        post = read(f"{kind}_posteriors")
        p = priors.loc[
            priors.scenario.eq("corrected") & priors.variant_type.eq(kind)
        ].iloc[0]
        fraction = post.affected.to_numpy() / post.n.to_numpy()
        weights = 1 - 1 / (post.n.to_numpy() + 0.01)
        mean = np.dot(weights, fraction) / weights.sum()
        variance = np.dot(weights, (fraction - mean) ** 2) / len(post)
        strength = mean * (1 - mean) / variance - 1
        alpha, beta = mean * strength, (1 - mean) * strength
        np.testing.assert_allclose(
            [mean, variance, strength, alpha, beta],
            p[
                ["mean", "variance", "strength", "alpha_empirical", "beta_empirical"]
            ].to_numpy(dtype=float),
            atol=1e-13,
        )
        np.testing.assert_allclose(
            post.posterior_alpha, alpha + post.affected, atol=1e-13
        )
        np.testing.assert_allclose(
            post.posterior_beta, beta + post.unaffected, atol=1e-13
        )
        np.testing.assert_allclose(
            post.unaffected, post.unaffected_literature + post.gnomad_carriers
        )
        np.testing.assert_allclose(
            post.posterior_mean,
            (alpha + post.affected) / (strength + post.n),
            atol=1e-13,
        )
        assert len(post) == p.variants
        orientation.append(
            {
                "variant_type": kind,
                "variants": len(post),
                "mean": mean,
                "alpha_affected": alpha,
                "beta_unaffected": beta,
                "strength": strength,
                "affected": int(post.affected.sum()),
                "unaffected_literature": int(post.unaffected_literature.sum()),
                "gnomad_carriers": int(post.gnomad_carriers.sum()),
            }
        )

    variants = read("missense_posteriors")
    assert primary.variant_id.tolist() == variants.unit_id.tolist()
    variants["variant_id"] = variants.unit_id
    variants["canonical_variant_id"] = variants.unit_id
    variants["canonical_pos"] = variants.aa_pos.astype(int)
    variants["donor_eligible"] = True
    geometry = pd.read_csv(
        EVIDENCE / "structural_sanity_20260913/geometry/BRCA2/primary_geometry.csv.gz",
        low_memory=False,
    ).fillna({"idr_segment": ""})
    experimental = set(
        geometry.loc[
            geometry.frame_id.isin(["7LDG", "8PBC"])
            & geometry.geometry_state.eq("structured"),
            "canonical_pos",
        ]
    )
    groups = {
        "experimental": primary.canonical_pos.isin(experimental),
        "alphafold": primary.density_source.eq("structured")
        & ~primary.canonical_pos.isin(experimental),
        "polymer": primary.density_source.eq("polymer"),
        "zero_affected_donors": primary.affected_donor_count.eq(0),
    }
    selected = []
    for mask in groups.values():
        eligible = primary.loc[mask & primary.density.notna()]
        selected += eligible.variant_id.iloc[[0, len(eligible) // 2]].tolist()
    selected = list(dict.fromkeys(selected))
    reference = empirical_variant_density(
        variants,
        geometry.loc[geometry.geometry_state.isin(["structured", "idr"])],
        target_ids=selected,
        backend="reference",
        include_context_weights=True,
    )
    v = variants.set_index("variant_id")
    p = priors.loc[
        priors.scenario.eq("corrected") & priors.variant_type.eq("missense")
    ].iloc[0]
    oldp = priors.loc[
        priors.scenario.eq("before") & priors.variant_type.eq("missense")
    ].iloc[0]
    contexts = []
    for (target, context), rows in reference.context_weights.groupby(
        ["variant_id", "context_id"]
    ):
        assert target not in set(rows.donor_id) and rows.donor_id.is_unique
        donors = v.loc[rows.donor_id]
        raw = 2 / (1 + np.exp(np.log(3) * rows.distance.to_numpy() / 3))
        q = raw / raw.sum()
        np.testing.assert_allclose(q, rows.context_normalized_weight, atol=1e-12)
        affected, unaffected, n = [
            donors[field].to_numpy() for field in ["affected", "unaffected", "n"]
        ]
        a, u = raw @ affected, raw @ unaffected
        retention = p.strength / (p.strength + n)
        contexts.append(
            {
                "variant_id": target,
                "context_id": context,
                "density": q @ donors.posterior_mean,
                "raw_variant_fraction_density": q @ (affected / n),
                "prior_component": q @ retention * p["mean"],
                "counts_component": q @ ((1 - retention) * affected / n),
                "old_prior_corrected_counts_density": q
                @ ((oldp.alpha_empirical + affected) / (oldp.strength + n)),
                "raw_kernel_affected": a,
                "raw_kernel_unaffected": u,
                "raw_kernel_pooled_fraction": a / (a + u),
                "raw_kernel_one_prior_posterior": (p.alpha_empirical + a)
                / (p.strength + a + u),
                "raw_kernel_total_weight": raw.sum(),
            }
        )
    context_table = pd.DataFrame(contexts)
    computed = context_table.groupby("variant_id").mean(numeric_only=True)
    saved = primary.set_index("variant_id").loc[computed.index, computed.columns]
    np.testing.assert_allclose(computed, saved, atol=2e-11, rtol=1e-12)
    context_table.to_csv(HERE / "selected_context_formula_checks.csv", index=False)
    direct = computed.reset_index()
    direct.to_csv(HERE / "selected_target_formula_checks.csv", index=False)
    density_comparison = pd.read_csv(OUT / "density_comparison.csv")
    for row in density_comparison.loc[
        density_comparison.scenario.eq("corrected")
    ].itertuples():
        values = primary[row.score].dropna()
        np.testing.assert_allclose(
            [values.mean(), values.median(), values.min(), values.max()],
            [row.mean, row.median, row.minimum, row.maximum],
            atol=1e-13,
        )
        assert [
            len(values),
            values.le(0.001).sum(),
            values.le(0.01).sum(),
            values.eq(0).sum(),
        ] == [
            row.supported,
            row.at_or_below_0p1_percent,
            row.at_or_below_1_percent,
            row.exact_zero,
        ]
    residue = pd.read_csv(OUT / "BRCA2_residue_density.csv").set_index("canonical_pos")
    assert residue.index.tolist() == list(range(1, 3419))
    grouped = primary.groupby("canonical_pos")
    mean_fields = [
        "density",
        "raw_variant_fraction_density",
        "raw_kernel_pooled_fraction",
        "raw_kernel_one_prior_posterior",
        "old_prior_corrected_counts_density",
    ]
    for field in mean_fields:
        expected_residue = grouped[field].mean().reindex(residue.index)
        np.testing.assert_allclose(expected_residue, residue[field], atol=1e-13)
    for field in ["affected", "unaffected"]:
        np.testing.assert_allclose(
            grouped[field].sum().reindex(residue.index), residue[field]
        )
    np.testing.assert_allclose(
        grouped.density.count().reindex(residue.index), residue.supported_variants
    )
    np.testing.assert_allclose(grouped.size().reindex(residue.index), residue.variants)
    old_primary = pd.concat(
        [
            pd.read_csv(path, low_memory=False)
            for path in sorted(
                (EVIDENCE / "structural_sanity_20260913/analysis/BRCA2").glob(
                    "primary_density.part*.csv.gz"
                )
            )
        ],
        ignore_index=True,
    )
    np.testing.assert_allclose(
        old_primary.groupby("canonical_pos").density.mean().reindex(residue.index),
        residue.before_density,
        atol=1e-13,
    )
    no_affected = primary.loc[
        primary.density.notna() & primary.affected_donor_count.eq(0)
    ]
    assert no_affected.raw_variant_fraction_density.eq(0).all()
    segment = primary.loc[primary.canonical_pos.isin([865, 866])]
    report_claims = {
        "supported_targets": int(primary.density.notna().sum()),
        "supported_residues": int(residue.density.notna().sum()),
        "supported_geometry": primary.loc[primary.density.notna()]
        .density_source.value_counts()
        .to_dict(),
        "no_affected_donor_targets": len(no_affected),
        "no_affected_donor_residues": no_affected.canonical_pos.nunique(),
        "no_affected_donor_density_minimum": no_affected.density.min(),
        "no_affected_donor_density_maximum": no_affected.density.max(),
        "residues_865_866": {
            "variants": len(segment),
            "affected": float(segment.affected.sum()),
            "unaffected": float(segment.unaffected.sum()),
            "mean_density": segment.density.mean(),
            "before_mean_density": old_primary.loc[
                old_primary.canonical_pos.isin([865, 866])
            ].density.mean(),
            "mean_raw_variant_fraction_density": segment.raw_variant_fraction_density.mean(),
        },
        "zero_affected_singleton_posterior": p.alpha_empirical / (p.strength + 1),
        "zero_affected_unaffected_count_needed_for_mean_at_most_0p1_percent": int(
            np.ceil(p.alpha_empirical / 0.001 - p.strength)
        ),
    }
    receipt = {
        "all_population_alleles_preserved_once": len(pop),
        "all_consequence_gnomad_carriers_preserved": int(pop.gnomad_carriers.sum()),
        "clinical_count_deltas_match_source_ledger": True,
        "zero_clinical_evidence_keys_removed": len(zero_keys),
        "population_alleles_released_to_genomic_units": len(released),
        "population_alleles_captured_by_restored_patient_keys": len(captured),
        "released_frameshift_alleles_correctly_excluded_from_nonsense": len(
            classification
        ),
        "reclassified_frameshift_carriers_preserved_in_full_union": int(
            classification.gnomad_carriers.sum()
        ),
        "independent_moments_and_orientation": orientation,
        "reference_targets": len(selected),
        "reference_contexts": len(context_table),
        "maximum_selected_formula_error": float(abs(computed - saved).max().max()),
        "absolute_kernel_count_pooling_verified": True,
        "one_prior_per_supported_context_then_equal_context_mean_verified": True,
        "corrected_metric_rows_recomputed": int(
            density_comparison.scenario.eq("corrected").sum()
        ),
        "all_residue_table_rows_checked": len(residue),
        "report_numerical_claims": report_claims,
        "script_sha256": sha(__file__),
        "runner_sha256": sha(AUDIT / "recompute_brca2.py"),
    }
    if (OUT / "input_hashes.json").exists():
        inputs = json.loads((OUT / "input_hashes.json").read_text())
        for relative, expected_hash in inputs.items():
            assert sha(EVIDENCE / relative) == expected_hash, relative
        receipt["final_input_hashes_verified"] = len(inputs)
    else:
        raise RuntimeError(
            "Wait for corrected/input_hashes.json before finalizing audit"
        )
    (HERE / "audit.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
