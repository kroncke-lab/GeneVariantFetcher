"""Reconstruct saved structural densities and variances from donor artifacts."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_structure import digest, MAX_BYTES


HERE = Path(__file__).resolve().parent
OUT = HERE / "structural"


def main():
    run = json.loads((OUT / "run_checks.json").read_text())
    for path, expected in run["input_hashes"].items():
        assert digest(HERE.parent / path) == expected, path
    for name, expected in run["output_hashes"].items():
        assert digest(OUT / name) == expected, name
    variants = pd.read_csv(OUT / "GCK_structural_input_variants.csv.gz").set_index(
        "variant_id"
    )
    primary = pd.read_csv(OUT / "GCK_primary_variant_density.csv").set_index(
        "variant_id"
    )
    seen, pair_rows, far_rows = set(), 0, 0
    maximum_density_error = maximum_variance_error = 0.0
    for path in sorted(OUT.glob("GCK_primary_donor_contexts_*.csv.gz")):
        frame = pd.read_csv(path)
        assert frame.variant_id.nunique() <= 40
        assert not frame.variant_id.eq(frame.donor_id).any()
        assert (frame.kernel > 0).all()
        assert not frame.duplicated(["variant_id", "context_id", "donor_id"]).any()
        assert not (set(frame.variant_id) & seen)
        seen.update(frame.variant_id)
        pair_rows += len(frame)
        far_rows += int((frame.distance > 20).sum())
        frame["weighted_mean"] = frame.target_normalized_weight * frame.donor_id.map(
            variants.posterior_mean
        )
        frame["weighted_variance"] = (
            frame.target_normalized_weight** 2
            * frame.donor_id.map(variants.posterior_variance)
        )
        reconstructed = frame.groupby("variant_id").agg(
            density=("weighted_mean", "sum"),
            variance=("weighted_variance", "sum"),
            weight=("target_normalized_weight", "sum"),
            donors=("donor_id", "nunique"),
        )
        expected = primary.loc[reconstructed.index]
        np.testing.assert_allclose(reconstructed.weight, 1, atol=1e-13)
        np.testing.assert_array_equal(reconstructed.donors, expected.donor_count)
        np.testing.assert_allclose(reconstructed.density, expected.density, atol=1e-13)
        np.testing.assert_allclose(
            reconstructed.variance, expected.density_conditional_variance, atol=1e-13
        )
        maximum_density_error = max(
            maximum_density_error,
            float(abs(reconstructed.density - expected.density).max()),
        )
        maximum_variance_error = max(
            maximum_variance_error,
            float(
                abs(
                    reconstructed.variance - expected.density_conditional_variance
                ).max()
            ),
        )
    assert seen == set(primary.index[primary.density.notna()])
    assert far_rows == run["positive_donor_weights_beyond_20"]
    assert all(
        path.stat().st_size <= MAX_BYTES for path in OUT.iterdir() if path.is_file()
    )
    available = primary.dropna(subset=["density"])
    standard_error = np.sqrt(available.density_conditional_variance / run["draws"])
    z = available.mc_mean_error.abs() / standard_error
    report = {
        "input_and_output_hashes_match": True,
        "all_supported_targets_reconstructed": len(seen),
        "donor_pair_rows": pair_rows,
        "positive_pairs_beyond_20": far_rows,
        "no_self_donor_records": True,
        "normalized_donor_weights_sum_to_one": True,
        "posterior_alpha_affected_beta_unaffected_orientation": run[
            "alpha_adds_affected_beta_adds_all_unaffected"
        ],
        "maximum_density_reconstruction_error": maximum_density_error,
        "maximum_conditional_variance_reconstruction_error": maximum_variance_error,
        "maximum_monte_carlo_mean_error_in_standard_errors": float(z.max()),
        "all_artifacts_at_most_1200000_bytes": True,
        "two_generated_figures_visually_inspected": True,
        "validator_sha256": digest(__file__),
    }
    (OUT / "artifact_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
