"""Verify pilot identities, uncertainty, donor coverage and optional replay."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reproduction", type=Path)
    args = parser.parse_args()
    out = HERE / "analysis"
    full = pd.read_csv(out / "GCK_empirical_posteriors_and_eligibility.csv")
    primary = pd.read_csv(out / "GCK_primary_variant_density.csv")
    scenarios = pd.read_csv(out / "GCK_density_scenarios.csv.gz")
    weights = pd.read_csv(
        out / "GCK_primary_normalized_donor_weights.csv.gz"
    ).set_index("variant_id")
    pairs = pd.concat(
        [
            pd.read_csv(path)
            for path in sorted(out.glob("GCK_primary_donor_contexts_*.csv.gz"))
        ],
        ignore_index=True,
    )
    assert len(full) == 398 and not full.key.duplicated().any()
    assert len(primary) == 249 and not primary.key.duplicated().any()
    assert scenarios.scenario.nunique() == 20
    assert not scenarios.duplicated(["scenario", "variant_id"]).any()
    assert len(pairs) == 245 * 244
    assert not pairs.duplicated(["variant_id", "donor_id"]).any()
    assert (pairs.variant_id != pairs.donor_id).all()
    assert (pairs.kernel > 0).all() and (pairs.distance > 20).any()
    assert (pairs.loc[pairs.distance > 20, "target_normalized_weight"] > 0).all()
    assert (pairs.donor_chain == "A").all() and (pairs.donor_frame_id == "1V4S").all()
    assert (pairs.groupby("variant_id").size() == 244).all()
    reconstructed = pairs.pivot(
        index="variant_id", columns="donor_id", values="target_normalized_weight"
    )
    reconstructed = reconstructed.reindex(
        index=weights.index, columns=weights.columns
    ).fillna(0)
    np.testing.assert_allclose(reconstructed, weights, atol=1e-14)
    assert np.all(np.diag(weights.to_numpy()) == 0)
    assert (primary.same_residue_donors > 0).sum() == 121
    assert set(primary.loc[primary.density.isna(), "key"]) == {
        "A11T",
        "D4N",
        "M462I",
        "M462V",
    }
    for scenario in ["1V4S_com_h3", "1V4T_com_h3", "AF_P35557_com_h3"]:
        subset = scenarios.loc[scenarios.scenario.eq(scenario)]
        unavailable = subset.density.isna()
        for field in [
            "density_lower_95",
            "density_upper_95",
            "density_conditional_variance",
        ]:
            assert subset.loc[unavailable, field].isna().all()
            assert np.isfinite(subset.loc[~unavailable, field]).all()
        assert (subset.loc[~unavailable, "density_conditional_variance"] > 0).all()
    open_primary = scenarios.loc[scenarios.scenario.eq("1V4T_com_h3")]
    assert open_primary.density_source.eq("polymer").sum() == 16
    assert (
        open_primary.loc[open_primary.density_source.eq("polymer"), "aa_pos"]
        .between(157, 179)
        .all()
    )
    predictions = pd.read_csv(out / "GCK_variant_loo_predictions.csv")
    keysets = predictions.groupby("model").key.apply(set)
    assert len(keysets) == 9 and all(keys == keysets.iloc[0] for keys in keysets)
    assert len(keysets.iloc[0]) == 242
    checks = json.loads((out / "run_checks.json").read_text())
    repo = HERE.parents[2]
    module = (
        repo.parent / "ProteinProximityAnalysis/src/alphafold_rin/empirical_density.py"
    )
    assert digest(module) == checks["ppa_module_sha256"]
    for name, expected in checks["input_hashes"].items():
        assert digest(HERE.parent / name) == expected, name
    geometry_checks = json.loads((HERE / "geometry/geometry_checks.json").read_text())
    for name, expected in geometry_checks["outputs_sha256"].items():
        assert digest(HERE / "geometry" / name) == expected, name
    replay = {}
    if args.reproduction:
        originals = {path.name: path for path in out.iterdir() if path.is_file()}
        reproduced = {
            path.name: path for path in args.reproduction.iterdir() if path.is_file()
        }
        assert set(originals) == set(reproduced)
        for name, original in originals.items():
            replay[name] = digest(original) == digest(reproduced[name])
        assert all(replay.values()), [
            name for name, equal in replay.items() if not equal
        ]
    result = {
        "full_variant_count": len(full),
        "eligible_variant_count": len(primary),
        "primary_supported": 245,
        "primary_target_donor_pairs": len(pairs),
        "all_primary_geometric_donors_retained": True,
        "far_donors_have_positive_kernel_and_normalized_weight": True,
        "donor_chunks_reconstruct_normalized_matrix": True,
        "zero_self_weight": True,
        "other_same_residue_variants_retained": True,
        "missing_density_has_missing_interval_and_variance": True,
        "same_242_comparison_variants_all_nine_models": True,
        "ppa_source_and_input_and_geometry_hashes_match": True,
        "byte_identical_reproduction": replay,
    }
    (HERE / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
