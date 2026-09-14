"""Check saved sex/count provenance and independently reconstruct plot features."""

import hashlib
import json

import numpy as np
import pandas as pd
from scipy.sparse import load_npz, vstack

from fetch_population import HERE, REPO


def shards(stem):
    paths = sorted((HERE / "analysis/BRCA2_female_breast").glob(stem + ".part*.csv.gz"))
    assert paths, stem
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)


def main():
    sex = pd.read_csv(HERE / "population_sex_counts.csv.gz")
    members = pd.read_csv(HERE / "population_membership.csv.gz")
    obs = pd.read_csv(HERE / "observation_eligibility.csv.gz")
    priors = pd.read_csv(HERE / "prior_comparison.csv")
    assert sex.variant_id.is_unique and len(sex) == 4419
    assert members.variant_id.is_unique and set(members.variant_id) == set(
        sex.variant_id
    )
    probes = json.loads(
        (HERE.parent / "brca2_sex_endpoint_audit_20260914/sex_probe.json").read_text()
    )
    for probe in probes:
        row = sex.set_index("variant_id").loc[probe["variant_id"]]
        for counts in probe["assays"]["joint"]["sex"]:
            for field in ["ac", "an", "homozygote_count"]:
                assert row[counts["id"] + "_" + field] == counts[field]
    for col in ["ac", "an", "homozygote_count", "carriers"]:
        np.testing.assert_array_equal(
            sex["XX_" + col] + sex["XY_" + col], sex["all_" + col]
        )
    for scope in ["all", "XX", "XY"]:
        np.testing.assert_array_equal(
            sex[scope + "_ac"] - sex[scope + "_homozygote_count"],
            sex[scope + "_carriers"],
        )
    assert obs.source_row_id.is_unique
    retained = obs.loc[(obs.female_affected + obs.female_unaffected).gt(0)]
    assert set(retained.pmid) == {
        28283652,
        28222693,
        29755871,
        38863777,
        33278427,
        37719058,
        37816281,
        17100994,
        30287823,
    }
    for key, pmid in [
        ("E1581D", 17100994),
        ("V109G", 16672066),
        ("D2723H", 20927582),
        ("K2729N", 28222693),
    ]:
        assert not (retained.key.eq(key) & retained.pmid.eq(pmid)).any()
    han = retained.loc[retained.pmid.eq(17100994)]
    assert (
        len(han) == 1
        and han.female_affected.sum() == 0
        and han.female_unaffected.sum() == 5
    )
    japan = retained.loc[retained.pmid.eq(30287823)]
    assert japan.source_row_id.str.startswith("female_SD1:").all()
    assert (
        not retained.loc[retained.source_row_id.str.startswith("original:")]
        .pmid.eq(30287823)
        .any()
    )
    summary = []
    for kind in ["missense", "nonsense"]:
        post = shards(kind + "_posteriors")
        assert post.unit_id.is_unique and post.n.gt(0).all()
        assert post.canonical_wt_status.eq("match").all()
        assert set(post.vclass) <= (
            {"missense"} if kind == "missense" else {"nonsense", "stop_gained"}
        )
        m = members.loc[members.variant_type.eq(kind)]
        assert len(m) == (4231 if kind == "missense" else 188)
        expected = m.groupby("unit_id").XX_carriers.sum()
        np.testing.assert_array_equal(
            post.unit_id.map(expected).fillna(0), post.gnomad_carriers
        )
        assert post.gnomad_carriers.sum() == m.XX_carriers.sum()
        ledger = pd.read_csv(
            HERE / f"analysis/BRCA2_female_breast/{kind}_clinical_join_ledger.csv.gz"
        )
        eligible = ledger.loc[ledger.prior_literature_eligible]
        if kind == "nonsense":
            assert not eligible.key.eq("N1098X").any()
            excluded = ledger.loc[ledger.key.eq("N1098X")]
            assert len(excluded) == 1 and excluded.affected_literature.sum() == 3
        assert post.affected.sum() == eligible.affected_literature.sum()
        assert post.unaffected_literature.sum() == eligible.unaffected_literature.sum()
        np.testing.assert_array_equal(
            post.unaffected_literature + post.gnomad_carriers, post.unaffected
        )
        np.testing.assert_array_equal(post.affected + post.unaffected, post.n)
        p = priors.loc[
            priors.input_scope.eq("female_breast") & priors.variant_type.eq(kind)
        ].iloc[0]
        y = post.affected / post.n
        w = 1 - 1 / (post.n + 0.01)
        mu = (w * y).sum() / w.sum()
        variance = (w * (y - mu) ** 2).sum() / len(post)
        strength = mu * (1 - mu) / variance - 1
        np.testing.assert_allclose(
            [mu, variance, mu * strength, (1 - mu) * strength],
            [p["mean"], p.variance, p.alpha_empirical, p.beta_empirical],
            rtol=1e-12,
        )
        np.testing.assert_allclose(
            post.posterior_alpha, p.alpha_empirical + post.affected
        )
        np.testing.assert_allclose(
            post.posterior_beta, p.beta_empirical + post.unaffected
        )
        np.testing.assert_allclose(
            post.posterior_mean,
            post.posterior_alpha / (post.posterior_alpha + post.posterior_beta),
        )
        summary.append(
            dict(
                variant_type=kind,
                variants=len(post),
                affected=int(post.affected.sum()),
                clinical_unaffected=int(post.unaffected_literature.sum()),
                XX_unaffected=int(post.gnomad_carriers.sum()),
                clinical_units=int(
                    post.affected.add(post.unaffected_literature).gt(0).sum()
                ),
                clinical_only_A_positive_U_zero=int(
                    (post.affected.gt(0) & post.unaffected.eq(0)).sum()
                ),
                posterior_max=float(post.posterior_mean.max()),
            )
        )
    post = shards("missense_posteriors")
    density = shards("density_h3")
    assert list(density.variant_id) == list(post.unit_id)
    weights_dir = REPO / "results/residue_density_refresh_20260914/BRCA2_female_breast"
    W = vstack([load_npz(p) for p in sorted(weights_dir.glob("weights_*.npz"))]).tocsr()
    assert W.shape == (len(post), len(post)) and not W.diagonal().any()
    mass = np.asarray(W.sum(axis=1)).ravel()
    supported = mass > 0
    np.testing.assert_allclose(mass[supported], 1, atol=1e-12)
    assert np.array_equal(supported, density.density.notna())
    for column, values in [
        ("density", post.posterior_mean),
        ("raw_variant_fraction_density", post.affected / post.n),
    ]:
        rebuilt = W @ values.to_numpy()
        rebuilt[~supported] = np.nan
        np.testing.assert_allclose(rebuilt, density[column], atol=1e-12, equal_nan=True)
    np.testing.assert_allclose(
        density.prior_component + density.counts_component,
        density.density,
        atol=1e-12,
        equal_nan=True,
    )
    row, col = W.nonzero()
    same_residue = post.aa_pos.to_numpy()[row] == post.aa_pos.to_numpy()[col]
    assert same_residue.any()
    table = pd.read_csv(HERE / "BRCA2_residue_density.csv.gz")
    assert len(table) == 3418 and table.canonical_pos.is_unique
    assert table.variants.sum() == len(post)
    grouped = (
        density.groupby("canonical_pos").density.mean().reindex(table.canonical_pos)
    )
    np.testing.assert_allclose(grouped, table.density, atol=1e-12, equal_nan=True)
    assert table.loc[table.supported_variants.eq(0), "density"].isna().all()
    peak = int(table.loc[table.raw_variant_fraction_density.idxmax(), "canonical_pos"])
    peak_weights = np.asarray(W[post.index[post.aa_pos.eq(peak)]].mean(axis=0)).ravel()
    contributions = peak_weights * post.affected.to_numpy() / post.n.to_numpy()
    donor = int(np.argmax(contributions))
    peak_donor = dict(
        residue=peak,
        donor=str(post.iloc[donor].protein_key),
        affected=int(post.iloc[donor].affected),
        unaffected=int(post.iloc[donor].unaffected),
        normalized_weight=float(peak_weights[donor]),
        gold_contribution=float(contributions[donor]),
    )
    # Source hashes are receipts, not live lookup pointers. Verify local inputs
    # and raw replies have not drifted between extraction, fitting and review.
    checked = 0
    for name in [
        "clinical_input_hashes.json",
        "run_input_hashes.json",
        "population_sex_counts_receipt.json",
    ]:
        hashes = json.loads((HERE / name).read_text())
        if "input_hashes" in hashes:
            hashes = hashes["input_hashes"]
        for name, expected in hashes.items():
            assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected, (
                name
            )
            checked += 1
    result = dict(
        passed=True,
        previous_six_sex_probes_reproduced=True,
        input_hashes_checked=checked,
        types=summary,
        supported_variants=int(supported.sum()),
        supported_residues=int(table.supported_variants.gt(0).sum()),
        missing_residues=int(table.supported_variants.eq(0).sum()),
        same_residue_nonself_weight_pairs=int(same_residue.sum()),
        median_residue_blue=float(table.density.median()),
        median_residue_gold=float(table.raw_variant_fraction_density.median()),
        gold_zero_residues=int(table.raw_variant_fraction_density.eq(0).sum()),
        gold_at_or_below_0p1_percent_residues=int(
            table.raw_variant_fraction_density.le(0.001).sum()
        ),
        largest_gold_residue_main_donor=peak_donor,
        max_residue_blue=float(table.density.max()),
        max_residue_gold=float(table.raw_variant_fraction_density.max()),
        median_weight_share_beyond_20=float(table.mean_weight_share_beyond_20.median()),
    )
    (HERE / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
