"""Independent count, geometry, weight and residue checks of the frozen refresh."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.sparse import load_npz

HERE = Path(__file__).resolve().parent
REFRESH = HERE.parent
EVIDENCE = REFRESH.parent
REPO = EVIDENCE.parents[1]
PPA = REPO.parent / "ProteinProximityAnalysis/src/alphafold_rin"
sys.path.insert(0, str(PPA.parent))
from alphafold_rin.empirical_density import empirical_variant_density

spec = importlib.util.spec_from_file_location(
    "refresh_inputs", REFRESH / "run_refresh.py"
)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

GEOMETRY = {
    "GCK": "gck_structural_pilot_20260912/geometry/1V4S_canonical_geometry.csv",
    "BRCA2": "structural_sanity_20260913/geometry/BRCA2/primary_geometry.csv.gz",
    "HNF1A": "missense_structural_extension_20260912/geometry/HNF1A/8PI8_canonical_geometry.csv",
    "LDLR": "missense_structural_extension_20260912/geometry/LDLR/AF_P01130_canonical_geometry.csv",
    "KCNQ1": "missense_structural_extension_20260912/geometry/KCNQ1/9U7F_canonical_geometry.csv",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(folder, stem):
    paths = sorted(folder.glob(f"{stem}.part*.csv.gz"))
    assert paths, (folder, stem)
    return pd.concat(
        [pd.read_csv(p, low_memory=False) for p in paths], ignore_index=True
    )


def close(a, b):
    np.testing.assert_allclose(a, b, rtol=1e-11, atol=2e-12, equal_nan=True)


def clinical_check(gene):
    original_path = (
        EVIDENCE
        / "brca2_distance_prior_audit_20260913/source/BRCA2_literature_identity_corrected.csv.gz"
        if gene == "BRCA2"
        else EVIDENCE
        / "population_inclusive_penetrance_20260912/audit/literature_identity_flags.csv.gz"
    )
    original = pd.read_csv(original_path, low_memory=False)
    original = original.loc[original.gene.eq(gene)].set_index("key")
    expected = original[["affected_literature", "unaffected_literature"]].copy()
    out = REFRESH / "source" / gene
    if gene == "BRCA2":
        decisions = pd.read_csv(out / "BRCA2_observation_decisions.csv.gz")
        for key, group in decisions.groupby("original_key"):
            expected.loc[key] += [
                group.A_to_restore.sum() - group.A_to_remove.sum(),
                group.U_to_restore.sum() - group.U_to_remove.sum(),
            ]
        corrected = pd.read_csv(
            out / "BRCA2_literature_identity_corrected.csv.gz", low_memory=False
        )
    elif gene == "GCK":
        decisions = pd.read_csv(out / "GCK_observation_corrections.csv")
        for key, group in decisions.groupby("key"):
            expected.loc[key] += [
                group.affected_restored.sum() - group.affected_removed.sum(),
                group.unaffected_restored.sum() - group.unaffected_removed.sum(),
            ]
        proxy = pd.read_csv(
            out / "clinical_input_diabetes_proxy.csv.gz", low_memory=False
        ).set_index("key")
        close(expected.reindex(proxy.index), proxy[expected.columns])
        pending = pd.read_csv(out / "GCK_diabetes_endpoint_pending.csv")
        for key, group in pending.groupby("key"):
            expected.loc[key] -= [group.affected.sum(), group.unaffected.sum()]
        corrected = pd.read_csv(out / "clinical_input.csv.gz", low_memory=False)
    elif (
        gene in {"HNF1A", "KCNQ1", "LDLR"} and (out / "clinical_input.csv.gz").exists()
    ):
        decision_paths = [
            out / name
            for name in [
                "observation_decisions.csv",
                "observation_decisions.csv.gz",
                "observation_corrections.csv",
            ]
            if (out / name).exists()
        ]
        assert len(decision_paths) == 1, decision_paths
        decisions = pd.read_csv(decision_paths[0]).rename(
            columns={
                "affected_removed": "A_to_remove",
                "unaffected_removed": "U_to_remove",
                "affected_restored": "A_to_restore",
                "unaffected_restored": "U_to_restore",
            }
        )
        for key, group in decisions.groupby("key"):
            expected.loc[key] += [
                group.A_to_restore.sum() - group.A_to_remove.sum(),
                group.U_to_restore.sum() - group.U_to_remove.sum(),
            ]
        corrected = pd.read_csv(out / "clinical_input.csv.gz", low_memory=False)
    else:
        corrected = original.reset_index()
    corrected = corrected.set_index("key")
    assert set(expected.index) == set(corrected.index)
    close(expected.reindex(corrected.index), corrected[expected.columns])
    assert corrected[expected.columns].ge(0).all().all()
    return corrected, {
        "source_key_delta_arithmetic_verified": True,
        "clinical_affected_delta": float(
            expected.affected_literature.sum() - original.affected_literature.sum()
        ),
        "clinical_unaffected_delta": float(
            expected.unaffected_literature.sum() - original.unaffected_literature.sum()
        ),
    }


def verify(gene, population):
    out = REFRESH / "analysis" / gene
    pop = population.loc[
        population.gene.eq(gene) & population.population_eligible
    ].set_index("variant_id")
    union = read(out, "complete_union")
    membership = read(out, "population_membership")
    assert union.unit_id.is_unique and membership.variant_id.is_unique
    assert membership._merge.eq("both").all() and set(membership.variant_id) == set(
        pop.index
    )
    assert union.n.ge(1).all() and union[["affected", "unaffected"]].ge(0).all().all()
    close(union.unaffected, union.unaffected_literature + union.gnomad_carriers)
    close(union.n, union.affected + union.unaffected)
    members = union.loc[
        union.member_alleles.notna(), ["unit_id", "member_alleles"]
    ].copy()
    members["variant_id"] = members.member_alleles.str.split(";")
    members = members.explode("variant_id")
    assert members.variant_id.is_unique and set(members.variant_id) == set(pop.index)
    assert (
        members.set_index("variant_id").unit_id.reindex(membership.variant_id).tolist()
        == membership.unit_id_after.tolist()
    )
    counts = members.merge(
        pop[["gnomad_carriers"]],
        left_on="variant_id",
        right_index=True,
        validate="one_to_one",
    )
    close(
        counts.groupby("unit_id")
        .gnomad_carriers.sum()
        .reindex(union.unit_id, fill_value=0),
        union.gnomad_carriers,
    )
    assert union.gnomad_carriers.sum() == pop.gnomad_carriers.sum()
    corrected, receipt = clinical_check(gene)
    retained = corrected.affected_literature.add(corrected.unaffected_literature).gt(0)
    ledger = read(out, "clinical_join_ledger")
    assert set(ledger.key) == set(corrected.index[retained])
    excluded_keys = ledger.loc[~ledger.prior_literature_eligible]
    assert excluded_keys.join_exclusion_reason.notna().all()
    assert set(union.literature_key.dropna()) == set(
        ledger.loc[ledger.prior_literature_eligible, "key"]
    )
    receipt["clinical_identity_keys_excluded"] = len(excluded_keys)
    receipt["clinical_identity_exclusion_reasons"] = (
        excluded_keys.join_exclusion_reason.value_counts().to_dict()
    )
    clinical = union.loc[union.literature_key.notna()].set_index("literature_key")
    close(clinical.affected, corrected.affected_literature.reindex(clinical.index))
    close(
        clinical.unaffected_literature,
        corrected.unaffected_literature.reindex(clinical.index),
    )

    if gene == "BRCA2":
        old_nonsense = read(
            EVIDENCE / "brca2_distance_prior_audit_20260913/corrected",
            "nonsense_posteriors",
        )
        new_nonsense = read(out, "nonsense_posteriors")
        moved = (
            membership.loc[
                membership.unit_id_before.isin(old_nonsense.unit_id)
                & ~membership.unit_id_after.isin(new_nonsense.unit_id)
            ]
            .merge(
                union[["unit_id", "vclass", "gnomad_carriers", "protein_key"]],
                left_on="unit_id_after",
                right_on="unit_id",
                validate="many_to_one",
            )
            .merge(
                pop[
                    [
                        "ref",
                        "alt",
                        "hgvsc",
                        "hgvsp",
                        "transcript_id",
                        "canonical_transcript_id",
                    ]
                ],
                left_on="variant_id",
                right_index=True,
                validate="one_to_one",
            )
        )
        moved["net_length_change"] = moved.alt.str.len() - moved.ref.str.len()
        assert moved.vclass.eq("frameshift").all()
        assert moved.net_length_change.mod(3).ne(0).all()
        assert len(moved) == 3 and moved.gnomad_carriers.sum() == 4
        moved.to_csv(HERE / "BRCA2_released_frameshift_alleles.csv", index=False)
        receipt["released_frameshift_units_outside_nonsense"] = 3
        receipt["released_frameshift_carriers_preserved"] = 4

    priors = pd.read_csv(out / "prior_comparison.csv")
    prior_rows = []
    for kind in ["missense", "nonsense"]:
        post = read(out, kind + "_posteriors")
        expected = union.loc[
            union.vclass.isin([kind, "stop_gained"] if kind == "nonsense" else [kind])
            & union.canonical_wt_status.eq("match")
        ]
        assert set(post.unit_id) == set(expected.unit_id)
        p = priors.loc[
            priors.scenario.eq("refreshed") & priors.variant_type.eq(kind)
        ].iloc[0]
        y, w = post.affected / post.n, 1 - 1 / (post.n + 0.01)
        mu = np.dot(w, y) / w.sum()
        variance = np.dot(w, (y - mu) ** 2) / len(post)
        k = mu * (1 - mu) / variance - 1
        alpha, beta = mu * k, (1 - mu) * k
        assert alpha > 0 and beta > 0 and np.isfinite([alpha, beta]).all()
        close(
            [mu, variance, k, alpha, beta],
            p[
                ["mean", "variance", "strength", "alpha_empirical", "beta_empirical"]
            ].to_numpy(float),
        )
        close(post.posterior_alpha, alpha + post.affected)
        close(
            post.posterior_beta,
            beta + post.unaffected_literature + post.gnomad_carriers,
        )
        close(post.posterior_mean, (alpha + post.affected) / (k + post.n))
        # Nondegenerate orientation and posterior monotonicity, with fixed prior.
        swapped_mean = np.dot(w, 1 - y) / w.sum()
        swapped_v = np.dot(w, (1 - y - swapped_mean) ** 2) / len(post)
        swapped_k = swapped_mean * (1 - swapped_mean) / swapped_v - 1
        close([swapped_mean * swapped_k, (1 - swapped_mean) * swapped_k], [beta, alpha])
        assert np.all(
            (alpha + post.affected + 1) / (k + post.n + 1) > post.posterior_mean
        )
        assert np.all((alpha + post.affected) / (k + post.n + 1) < post.posterior_mean)
        prior_rows.append(
            dict(
                gene=gene,
                variant_type=kind,
                variants=len(post),
                affected=post.affected.sum(),
                unaffected_literature=post.unaffected_literature.sum(),
                gnomad_carriers=post.gnomad_carriers.sum(),
                mean=mu,
                alpha=alpha,
                beta=beta,
                strength=k,
                zero_case_singleton_posterior=alpha / (k + 1),
            )
        )

    variants, density = read(out, "missense_posteriors"), read(out, "density_h3")
    ids = variants.unit_id.tolist()
    assert density.variant_id.tolist() == ids
    variants["variant_id"] = variants.unit_id
    variants["canonical_variant_id"] = variants.unit_id
    variants["canonical_pos"] = variants.aa_pos.astype(int)
    variants["donor_eligible"] = True
    p = priors.loc[
        priors.scenario.eq("refreshed") & priors.variant_type.eq("missense")
    ].iloc[0]
    a, u, n, y = [
        variants[c].to_numpy(float)
        for c in ["affected", "unaffected", "n", "posterior_mean"]
    ]
    geometry_path = EVIDENCE / GEOMETRY[gene]
    assert geometry_path.resolve() == runner.geometry_for(gene).resolve()
    geometry = pd.read_csv(geometry_path, low_memory=False).fillna({"idr_segment": ""})
    assert not geometry.duplicated(["frame_id", "chain", "canonical_pos"]).any()
    canonical = geometry.groupby("canonical_pos").aa_ref.agg(lambda s: set(s))
    for v in variants.itertuples():
        assert canonical.loc[v.canonical_pos] == {v.aa_ref}
    assert (
        density.half_distance.eq(3).all() and density.coordinate_metric.eq("com").all()
    )
    assert (
        density.polymer_scale.eq(3.8).all() and density.polymer_exponent.eq(0.5).all()
    )
    weight_receipt = json.loads((out / "density_checks.json").read_text())
    close(density.posterior_mean, variants.posterior_mean)
    same_position_rows, max_weight_error = 0, 0.0
    positions = variants.canonical_pos.to_numpy()
    cache = sorted(weight_receipt["raw_weight_shards"])
    cursor = 0
    for rel in cache:
        path = REPO / rel
        assert sha(path) == weight_receipt["raw_weight_shards"][rel]
        w = load_npz(path).toarray()
        batch = density.iloc[cursor : cursor + len(w)]
        assert w.shape == (len(batch), len(ids))
        assert np.isfinite(w).all() and (w >= 0).all()
        assert (w[np.arange(len(w)), np.arange(cursor, cursor + len(w))] == 0).all()
        supported = w.sum(axis=1) > 0
        assert np.array_equal(supported, batch.density.notna())
        close(w.sum(axis=1)[supported], 1)
        estimate = w @ y
        max_weight_error = max(
            max_weight_error,
            float(
                np.max(
                    np.abs(estimate[supported] - batch.density.to_numpy()[supported]),
                    initial=0,
                )
            ),
        )
        close(estimate[supported], batch.density.to_numpy()[supported])
        for name, vector in [
            ("raw_variant_fraction_density", a / n),
            ("neighborhood_prior_retention", p.strength / (p.strength + n)),
            ("prior_component", p.alpha_empirical / (p.strength + n)),
            ("counts_component", a / (p.strength + n)),
            ("affected_donor_weight_share", (a > 0).astype(float)),
        ]:
            close((w @ vector)[supported], batch[name].to_numpy()[supported])
            assert batch.loc[~supported, name].isna().all()
        close(
            ((w > 0) @ (a > 0).astype(int))[supported],
            batch.affected_donor_count.to_numpy()[supported],
        )
        close(
            (w**2 @ variants.posterior_variance.to_numpy())[supported],
            batch.conditional_donor_variance.to_numpy()[supported],
        )
        close(
            1 / np.sum(w[supported] ** 2, axis=1),
            batch.kish_donor_n.to_numpy()[supported],
        )
        same = (positions[None, :] == positions[cursor : cursor + len(w), None]) & (
            w > 0
        )
        close(same.sum(axis=1), batch.same_residue_donors)
        same_position_rows += int((same.sum(axis=1) > 0).sum())
        cursor += len(w)
    assert cursor == len(ids)
    close(density.prior_component + density.counts_component, density.density)
    assert (
        density.loc[density.density.isna(), "raw_kernel_one_prior_posterior"]
        .isna()
        .all()
    )

    selected = []
    for source_name in ["structured", "polymer", "mixed"]:
        candidates = density.loc[
            density.density_source.eq(source_name) & density.density.notna()
        ]
        if len(candidates):
            selected.append(candidates.iloc[len(candidates) // 2].variant_id)
    for mask in [
        density.same_residue_donors.gt(0),
        density.affected_donor_count.eq(0),
        density.density.isna(),
    ]:
        candidates = density.loc[mask]
        if len(candidates):
            selected.append(candidates.iloc[0].variant_id)
    selected = list(dict.fromkeys(selected))
    if gene in {"GCK", "LDLR"}:
        supported_pool = density.loc[density.density.notna()]
        for index in np.linspace(0, len(supported_pool) - 1, 8).astype(int):
            if len(set(selected) & set(supported_pool.variant_id)) >= 8:
                break
            selected.append(supported_pool.iloc[index].variant_id)
        selected = list(dict.fromkeys(selected))
    reference = empirical_variant_density(
        variants,
        geometry.loc[geometry.geometry_state.isin(["structured", "idr"])],
        target_ids=selected,
        backend="reference",
        include_context_weights=True,
    )
    current = density.set_index("variant_id").loc[selected]
    close(reference.summary.density, current.density)
    lookup = geometry.set_index(["frame_id", "chain", "canonical_pos"])
    v = variants.set_index("variant_id")
    context_rows = []
    for (target, context), pairs in reference.context_weights.groupby(
        ["variant_id", "context_id"]
    ):
        assert target not in set(pairs.donor_id) and pairs.donor_id.is_unique
        for pair in pairs.itertuples():
            tc = lookup.loc[(pair.frame_id, pair.chain, pair.canonical_pos)]
            dc = lookup.loc[
                (pair.donor_frame_id, pair.donor_chain, pair.donor_canonical_pos)
            ]
            assert pair.frame_id == pair.donor_frame_id
            if pair.source == "polymer":
                assert tc.geometry_state == dc.geometry_state == "idr"
                assert (
                    pair.chain == pair.donor_chain and tc.idr_segment == dc.idr_segment
                )
                distance = 3.8 * np.sqrt(
                    abs(pair.canonical_pos - pair.donor_canonical_pos)
                )
            else:
                assert tc.geometry_state == dc.geometry_state == "structured"
                distance = np.linalg.norm(
                    tc[["com_x", "com_y", "com_z"]].to_numpy(float)
                    - dc[["com_x", "com_y", "com_z"]].to_numpy(float)
                )
            close(distance, pair.distance)
        raw = 2 / (1 + np.exp(np.log(3) * pairs.distance.to_numpy() / 3))
        q = raw / raw.sum()
        assert (raw > 0).all()
        close(q, pairs.context_normalized_weight)
        donors = v.loc[pairs.donor_id]
        ak, uk = raw @ donors.affected, raw @ donors.unaffected
        context_rows.append(
            dict(
                gene=gene,
                variant_id=target,
                context_id=context,
                density=q @ donors.posterior_mean,
                raw_variant_fraction_density=q @ (donors.affected / donors.n),
                raw_kernel_affected=ak,
                raw_kernel_unaffected=uk,
                raw_kernel_total_weight=raw.sum(),
                raw_kernel_pooled_fraction=ak / (ak + uk),
                raw_kernel_one_prior_posterior=(p.alpha_empirical + ak)
                / (p.strength + ak + uk),
            )
        )
    context_table = pd.DataFrame(context_rows)
    computed = context_table.groupby("variant_id").mean(numeric_only=True)
    close(
        computed, density.set_index("variant_id").loc[computed.index, computed.columns]
    )
    # Distinct labels expose self inclusion; fixed-prior target changes do not alter its LOO score.
    supported_ids = current.index[current.density.notna()]
    target = supported_ids[0]
    changed = variants.copy()
    changed.loc[changed.variant_id.eq(target), "posterior_alpha"] += 123
    changed = changed.drop(columns=["posterior_mean"], errors="ignore")
    perturbed = empirical_variant_density(
        changed, geometry, target_ids=[target], include_context_weights=False
    )
    close(perturbed.summary.density.iloc[0], current.loc[target, "density"])
    close(perturbed.donor_weights, reference.donor_weights.loc[[target]])
    np.testing.assert_allclose(
        reference.donor_weights.loc[supported_ids].to_numpy() @ np.full(len(ids), 0.37),
        0.37,
    )

    old_error = None
    if (
        not (REFRESH / "source" / gene / "clinical_input.csv.gz").exists()
        and gene != "BRCA2"
    ):
        old = read(
            EVIDENCE / "missense_structural_extension_20260912/analysis" / gene,
            "primary_density",
        ).set_index("variant_id")
        assert set(old.index) == set(ids)
        close(density.density, old.density.reindex(ids))
        old_error = float(
            (density.set_index("variant_id").density - old.density).abs().max()
        )
        for kind in ["missense", "nonsense"]:
            pair = priors.loc[priors.variant_type.eq(kind)].set_index("scenario")
            cols = [
                "variants",
                "affected",
                "unaffected_literature",
                "gnomad_carriers",
                "mean",
                "variance",
                "strength",
                "alpha_empirical",
                "beta_empirical",
            ]
            close(
                pair.loc["previous", cols].to_numpy(float),
                pair.loc["refreshed", cols].to_numpy(float),
            )

    no_a = density.loc[density.density.notna() & density.affected_donor_count.eq(0)]
    assert no_a.raw_variant_fraction_density.eq(0).all()
    regions = []
    for (frame, chain, segment), rows in geometry.loc[
        geometry.geometry_state.eq("idr")
    ].groupby(["frame_id", "chain", "idr_segment"]):
        local = density.loc[density.canonical_pos.isin(rows.canonical_pos)]
        if len(local) and local.affected.sum() == 0:
            regions.append(
                dict(
                    gene=gene,
                    frame_id=frame,
                    chain=chain,
                    idr_segment=segment,
                    variants=len(local),
                    affected=0,
                    unaffected=local.unaffected.sum(),
                    supported=local.density.notna().sum(),
                    density_mean=local.density.mean(),
                )
            )
    receipt.update(
        dict(
            gene=gene,
            population_alleles=len(pop),
            population_carriers=int(pop.gnomad_carriers.sum()),
            population_reassignments=int(
                membership.unit_id_before.ne(membership.unit_id_after).sum()
            ),
            missense_donors=len(variants),
            supported_targets=int(density.density.notna().sum()),
            unsupported_targets=int(density.density.isna().sum()),
            supported_residues=int(
                density.loc[density.density.notna(), "canonical_pos"].nunique()
            ),
            geometry_support=density.density_source.value_counts().to_dict(),
            targets_with_same_residue_donors=same_position_rows,
            no_affected_donor_targets=len(no_a),
            no_affected_donor_residues=int(no_a.canonical_pos.nunique()),
            raw_fraction_zeros=int(density.raw_variant_fraction_density.eq(0).sum()),
            primary_zeros=int(density.density.eq(0).sum()),
            primary_median=float(density.density.median()),
            no_affected_density_min=float(no_a.density.min()) if len(no_a) else None,
            no_affected_density_max=float(no_a.density.max()) if len(no_a) else None,
            prior_fraction_of_score_median=float(
                (density.prior_component / density.density).median()
            ),
            all_weight_shards_checked=len(cache),
            max_weight_reproduction_error=max_weight_error,
            reference_targets=len(selected),
            reference_contexts=len(context_table),
            fixed_prior_target_count_perturbation_invariant=True,
            unchanged_gene_max_density_error=old_error,
        )
    )
    inputs = json.loads((out / "input_hashes.json").read_text())
    for rel, digest in inputs.items():
        assert sha(EVIDENCE / rel) == digest, rel
    receipt["verified_input_hashes"] = len(inputs)
    receipt["run_sha256"] = sha(REFRESH / "run_refresh.py")
    (HERE / f"{gene}_audit.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    )
    pd.DataFrame(prior_rows).to_csv(HERE / f"{gene}_prior_checks.csv", index=False)
    context_table.to_csv(HERE / f"{gene}_context_checks.csv", index=False)
    pd.DataFrame(
        regions,
        columns=[
            "gene",
            "frame_id",
            "chain",
            "idr_segment",
            "variants",
            "affected",
            "unaffected",
            "supported",
            "density_mean",
        ],
    ).to_csv(HERE / f"{gene}_zero_case_idr_segments.csv", index=False)
    print(json.dumps(receipt, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genes", nargs="+", default=list(GEOMETRY), choices=list(GEOMETRY)
    )
    args = parser.parse_args()
    runner.legacy.GENES = args.genes
    population, _ = runner.legacy.load_population()
    for gene in args.genes:
        verify(gene, population)
    (HERE / "audit_source_hashes.json").write_text(
        json.dumps(
            {
                str(path): sha(path)
                for path in [
                    Path(__file__),
                    PPA / "empirical_density.py",
                    PPA / "empirical_context.py",
                    REFRESH / "run_refresh.py",
                ]
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
