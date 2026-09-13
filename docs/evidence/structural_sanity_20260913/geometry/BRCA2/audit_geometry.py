"""Independent parsed-geometry and selected PPA reference checks; no outcomes fit."""

import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
OLD = HERE.parents[2] / "missense_structural_extension_20260912"
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density  # noqa: E402


def main():
    ledger = pd.read_csv(HERE / "canonical_confidence_map.csv.gz")
    fragments = pd.read_csv(HERE / "fragment_confidence.csv.gz")
    assert ledger.canonical_pos.tolist() == list(range(1, 3419))
    aggregate = fragments.groupby("canonical_pos").plddt.agg(
        ["min", "max", "median", "count"]
    )
    for source, column in [
        ("min", "min_plddt"),
        ("max", "max_plddt"),
        ("median", "median_plddt"),
        ("count", "overlapping_AF_fragments"),
    ]:
        np.testing.assert_allclose(aggregate[source], ledger[column], atol=1e-12)
    primary = pd.read_csv(HERE / "primary_geometry.csv.gz", low_memory=False).fillna(
        {"idr_segment": ""}
    )
    previous = pd.read_csv(
        OLD / "geometry/BRCA2/BRCA2_canonical_geometry.csv.gz", low_memory=False
    )
    experiments = primary.loc[
        primary.geometry_source.eq("experimental_biological_assembly")
    ]
    key = ["frame_id", "chain", "canonical_pos"]
    comparison_columns = ["aa_ref", "geometry_state", "ca_geometry_state"] + [
        f"{metric}_{axis}" for metric in ["com", "ca"] for axis in "xyz"
    ]
    pd.testing.assert_frame_equal(
        experiments.set_index(key)[comparison_columns].sort_index(),
        previous.set_index(key)[comparison_columns].sort_index(),
        check_exact=True,
    )
    source_groups = (
        primary.loc[primary.geometry_state.isin(["structured", "idr"])]
        .groupby("canonical_pos")
        .geometry_source.nunique()
    )
    assert source_groups.eq(1).all()
    row_counts = primary.groupby(["frame_id", "chain"]).size()
    assert row_counts.eq(3418).all() and not primary.duplicated(key).any()
    coords = [f"com_{axis}" for axis in "xyz"]
    assert (
        primary.loc[primary.geometry_state.eq("structured"), coords].notna().all().all()
    )
    assert (
        primary.loc[~primary.geometry_state.eq("structured"), coords].isna().all().all()
    )
    af = primary.loc[
        primary.geometry_source.eq("AF_fragment_local")
        & primary.geometry_state.eq("structured")
    ]
    assert af.plddt.ge(70).all()
    assert not set(af.canonical_pos) & set(
        experiments.loc[experiments.ca_geometry_state.eq("structured"), "canonical_pos"]
    )
    polymer = primary.loc[primary.geometry_state.eq("idr")]
    assert not polymer.canonical_pos.duplicated().any()
    assert (
        polymer.frame_id.eq("BRCA2_canonical_polymer").all()
        and polymer.chain.eq("A").all()
    )
    assert polymer.plddt.lt(50).all()
    for _, block in polymer.groupby("idr_segment"):
        assert block.canonical_pos.tolist() == list(
            range(block.canonical_pos.min(), block.canonical_pos.max() + 1)
        )
    # Read exactly the previously frozen count/posterior universe.
    variants = pd.concat(
        [
            pd.read_csv(path)
            for path in sorted(
                (OLD / "analysis/BRCA2").glob("input_variants.part*.csv.gz")
            )
        ],
        ignore_index=True,
    )
    routing = ledger.set_index("canonical_pos").exclusive_primary_route
    target_ids = []
    for route in ["experimental", "AF_ordered", "polymer"]:
        positions = set(routing.index[routing.eq(route)])
        candidates = variants.loc[
            variants.canonical_pos.isin(positions) & variants.canonical_pos.ne(2322)
        ]
        target_ids.extend(
            candidates.iloc[[0, len(candidates) // 2]].variant_id.tolist()
        )
    actual = empirical_variant_density(variants, primary, target_ids=target_ids)
    reference = empirical_variant_density(
        variants, primary, target_ids=target_ids, backend="reference"
    )
    np.testing.assert_allclose(
        actual.summary.density, reference.summary.density, atol=1e-11, equal_nan=True
    )
    np.testing.assert_allclose(
        actual.donor_weights, reference.donor_weights, atol=1e-11
    )
    contexts = reference.context_weights
    assert contexts.frame_id.eq(contexts.donor_frame_id).all()
    assert not contexts.variant_id.eq(contexts.donor_id).any()
    assert contexts.log_kernel.notna().all()
    poly_contexts = contexts.loc[contexts.source.eq("polymer")]
    assert poly_contexts.chain.eq(poly_contexts.donor_chain).all()
    np.testing.assert_allclose(
        poly_contexts.distance,
        3.8
        * np.sqrt(abs(poly_contexts.canonical_pos - poly_contexts.donor_canonical_pos)),
        atol=1e-12,
    )
    segment_map = polymer.set_index("canonical_pos").idr_segment
    assert all(
        segment_map[a] == segment_map[b]
        for a, b in zip(
            poly_contexts.canonical_pos, poly_contexts.donor_canonical_pos, strict=True
        )
    )
    # Structured domains retain their old donor contexts exactly.
    old_density = (
        pd.concat(
            [
                pd.read_csv(path)
                for path in sorted(
                    (OLD / "analysis/BRCA2").glob("primary_density.part*.csv.gz")
                )
            ]
        )
        .set_index("variant_id")
        .density
    )
    new_density = actual.summary.set_index("variant_id").density
    np.testing.assert_allclose(
        new_density.loc[target_ids[:2]], old_density.loc[target_ids[:2]], atol=1e-12
    )
    file_checks = []
    for path in sorted(HERE.glob("*.csv*")):
        payload = (
            gzip.decompress(path.read_bytes())
            if path.suffix == ".gz"
            else path.read_bytes()
        )
        assert b"\r" not in payload and path.stat().st_size < 1_150_000
        file_checks.append(
            {
                "file": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "LF": True,
            }
        )
    report = {
        "canonical_positions": len(ledger),
        "fragment_rows": len(fragments),
        "confidence_aggregates_reconstructed": True,
        "experimental_coordinates_parsed_exactly_unchanged": True,
        "full_canonical_rows_per_frame_chain": True,
        "exclusive_source_assignments": True,
        "AF_individual_fragment_threshold70": True,
        "polymer_consensus_max_below50": True,
        "one_canonical_polymer_context_per_position": True,
        "no_cross_frame_donor_distances": True,
        "same_segment_same_chain_polymer": True,
        "variant_only_self_exclusion": True,
        "selected_density_and_weights_reference_agreement": True,
        "selected_experimental_densities_unchanged": True,
        "reference_checked_target_ids": target_ids,
        "reference_density_rows": actual.summary[
            ["variant_id", "density", "density_source", "donor_count"]
        ].to_dict("records"),
        "files": file_checks,
    }
    (HERE / "geometry_audit.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(report["reference_density_rows"], indent=2))


if __name__ == "__main__":
    main()
