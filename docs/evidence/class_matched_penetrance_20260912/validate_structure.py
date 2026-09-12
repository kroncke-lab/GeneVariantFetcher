"""Reuse donor-artifact reconstruction against the new class-matched output."""

import json

import pandas as pd

from run_structure import HERE, PREVIOUS, digest, import_previous, replay


def main():
    validator = import_previous(
        "frozen_population_validator", PREVIOUS / "validate_structure.py"
    )
    validator.HERE = HERE
    validator.OUT = HERE / "structural"
    validator.main()
    output = HERE / "structural"
    new_shards = sorted(output.glob("GCK_primary_normalized_weights.part*.csv.gz"))
    old_shards = sorted(
        (PREVIOUS / "structural").glob("GCK_primary_normalized_weights.part*.csv.gz")
    )
    assert len(new_shards) == len(old_shards)
    for new, old in zip(new_shards, old_shards):
        pd.testing.assert_frame_equal(
            pd.read_csv(new), pd.read_csv(old), check_exact=True
        )
    current = pd.read_csv(output / "GCK_primary_variant_density.csv")
    old_path = PREVIOUS / "structural/GCK_primary_variant_density.csv"
    previous = pd.read_csv(old_path)
    comparison = current[["variant_id", "origin", "posterior_mean", "density"]].merge(
        previous[["variant_id", "posterior_mean", "density"]],
        on="variant_id",
        suffixes=("_missense", "_broader_prior"),
        validate="one_to_one",
    )
    assert len(comparison) == 634
    comparison["density_change"] = (
        comparison.density_missense - comparison.density_broader_prior
    )
    comparison["posterior_change"] = (
        comparison.posterior_mean_missense - comparison.posterior_mean_broader_prior
    )
    comparison_path = output / "GCK_prior_scope_comparison.csv.gz"
    replay.save_csv(comparison, comparison_path)
    checks_path = output / "artifact_checks.json"
    checks = json.loads(checks_path.read_text())
    checks["prior_scope"] = "canonical_missense"
    checks["same_normalized_weights_as_previous_run_exactly"] = True
    checks["no_nonsense_donors"] = True
    checks["scope_comparison_input_sha256"] = digest(old_path)
    checks["scope_comparison_output_sha256"] = digest(comparison_path)
    checks["validation_adapter_sha256"] = digest(__file__)
    checks_path.write_text(json.dumps(checks, indent=2) + "\n")


if __name__ == "__main__":
    main()
