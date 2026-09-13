"""Offline checks against the frozen PPA model and independent reference path."""

import copy
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2].parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density
from outer_loo import all_excluded_density


def fixture(copies=4):
    rng = np.random.default_rng(73)
    positions = [1, 1, 2, 3, 4, 5, 6, 7, 8, 999]
    data = pd.DataFrame(
        {
            "variant_id": [f"v{i}" for i in range(len(positions))],
            "canonical_pos": positions,
            "aa_ref": "A",
            "donor_eligible": True,
            "posterior_alpha": rng.uniform(0.1, 3, len(positions)),
            "posterior_beta": rng.uniform(0.1, 3, len(positions)),
        }
    )
    data["canonical_variant_id"] = data.variant_id
    alias = data.iloc[[0]].copy()
    alias["variant_id"] = "alias"
    data = pd.concat([data, alias], ignore_index=True)
    rows = []
    for frame in ["first", "other"]:
        for chain in range(copies):
            xyz = rng.normal(size=(8, 3)) * 20
            g = pd.DataFrame(
                {
                    "frame_id": frame,
                    "chain": str(chain),
                    "canonical_pos": np.arange(1, 9),
                    "aa_ref": "A",
                    "geometry_state": "structured",
                    "idr_segment": "",
                    "com_x": xyz[:, 0],
                    "com_y": xyz[:, 1],
                    "com_z": xyz[:, 2],
                }
            )
            g.loc[g.canonical_pos.isin([5, 6]), "geometry_state"] = "idr"
            g.loc[g.geometry_state.eq("idr"), "idr_segment"] = "tail"
            g.loc[g.canonical_pos.eq(7), "geometry_state"] = "ambiguous"
            g.loc[g.canonical_pos.eq(8), "geometry_state"] = "missing"
            if chain == 1:
                g.loc[g.canonical_pos.eq(2), "geometry_state"] = "missing"
            rows.append(g)
    return data, pd.concat(rows, ignore_index=True).sample(frac=1, random_state=4)


@pytest.mark.parametrize("copies", [1, 2, 4])
@pytest.mark.parametrize("block_size", [1, 3, 512])
def test_all_exclusions_match_model_and_fresh_reference(copies, block_size):
    data, geometry = fixture(copies)
    targets = data.variant_id.tolist()[::-1]
    cached = empirical_variant_density(
        data,
        geometry,
        target_ids=targets,
        include_context_model=True,
        include_context_weights=False,
    )
    result = all_excluded_density(cached.context_model, donor_batch_size=block_size)
    assert result.index.tolist() == targets
    assert result.columns.tolist() == list(cached.context_model.donor_ids)
    assert "alias" not in result.columns
    np.testing.assert_allclose(result.loc["v0"], result.loc["alias"], equal_nan=True)
    for donor in result.columns:
        np.testing.assert_allclose(
            result[donor],
            cached.context_model.density([donor]),
            atol=3e-14,
            rtol=2e-13,
            equal_nan=True,
        )
    for donor in ["v0", "v2", "v5"]:
        reference = empirical_variant_density(
            data,
            geometry,
            target_ids=targets,
            excluded_variant_ids=[donor],
            include_context_weights=False,
            backend="reference",
        )
        np.testing.assert_allclose(
            result[donor],
            reference.summary.density,
            atol=3e-14,
            rtol=2e-13,
            equal_nan=True,
        )
    assert result.loc["v9"].isna().all()


def test_two_variant_polymer_loses_support_without_zero_imputation():
    data, geometry = fixture(1)
    data = data.loc[data.variant_id.isin(["v5", "v6"])].copy()
    geometry = geometry.loc[geometry.canonical_pos.isin([5, 6])].copy()
    model = empirical_variant_density(
        data,
        geometry,
        include_context_model=True,
        include_context_weights=False,
    ).context_model
    result = all_excluded_density(model)
    assert result.loc["v5", "v6"] != result.loc["v5", "v6"]
    assert result.loc["v6", "v5"] != result.loc["v6", "v5"]
    assert result.loc["v5", "v5"] == pytest.approx(model.density()["v5"])
    assert result.loc["v6", "v6"] == pytest.approx(model.density()["v6"])


def test_no_mutation_and_target_chunks_match():
    data, geometry = fixture()
    all_targets = empirical_variant_density(
        data,
        geometry,
        include_context_model=True,
        include_context_weights=False,
    ).context_model
    before = copy.deepcopy(all_targets)
    result = all_excluded_density(all_targets)
    for name, value in vars(before).items():
        actual = getattr(all_targets, name)
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(actual, value)
        else:
            assert actual == value
    requested = ["alias", "v1", "v3", "v9"]
    partial = empirical_variant_density(
        data,
        geometry,
        include_context_model=True,
        include_context_weights=False,
        target_ids=requested,
    ).context_model
    pd.testing.assert_frame_equal(all_excluded_density(partial), result.loc[requested])


def test_fully_unmapped_context_batch():
    data, geometry = fixture()
    model = empirical_variant_density(
        data,
        geometry,
        target_ids=["v9"],
        include_context_model=True,
        include_context_weights=False,
    ).context_model
    assert all_excluded_density(model).isna().all().all()


def test_dominant_fallback_preserves_remote_surviving_donor():
    data, geometry = fixture(1)
    data = data.loc[data.variant_id.isin(["v0", "v1", "v2"])].copy()
    geometry = geometry.loc[
        geometry.frame_id.eq("first") & geometry.canonical_pos.isin([1, 2])
    ].copy()
    geometry[["com_y", "com_z"]] = 0.0
    geometry["com_x"] = (geometry.canonical_pos - 1) * 3000.0
    cached = empirical_variant_density(
        data,
        geometry,
        include_context_model=True,
        include_context_weights=False,
    )
    result = all_excluded_density(cached.context_model)
    assert result.loc["v0", "v1"] == pytest.approx(
        cached.donors.loc["v2", "posterior_mean"]
    )
    for donor in result.columns:
        np.testing.assert_allclose(
            result[donor], cached.context_model.density([donor]), atol=3e-14
        )


@pytest.mark.parametrize("batch_size", [0, -1, 1.5, True])
def test_invalid_block_size(batch_size):
    with pytest.raises(ValueError, match="positive integer"):
        all_excluded_density(None, donor_batch_size=batch_size)
