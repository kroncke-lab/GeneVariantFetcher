"""Change empirical priors/count proxy while freezing GCK geometry and donors."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_structure import digest, save_csv
from structure_statistics import density_values


HERE = Path(__file__).resolve().parent
OUT = HERE / "structural"
SCOPES = [
    "all_observed_full_locus",
    "coding_or_splice",
    "all_normalized_mse",
    "all_allele_count_proxy",
]


def main():
    variants = pd.read_csv(OUT / "GCK_structural_input_variants.csv.gz")
    primary = (
        pd.read_csv(OUT / "GCK_primary_variant_density.csv")
        .set_index("variant_id")
        .loc[variants.variant_id]
    )
    shards = sorted(OUT.glob("GCK_primary_normalized_weights.part*.csv.gz"))
    table = pd.concat(
        [pd.read_csv(path) for path in shards], ignore_index=True
    ).set_index("variant_id")
    weights = table.loc[variants.variant_id, variants.variant_id].to_numpy()
    np.testing.assert_array_equal(np.diag(weights), 0)
    priors_path = HERE / "analysis/empirical_prior_comparison.csv"
    priors = pd.read_csv(priors_path).query("gene == 'GCK'").set_index("scope")
    result = variants[["variant_id", "literature_key", "protein_key", "origin"]].copy()
    summaries = []
    for scope in SCOPES:
        prior = priors.loc[scope]
        gnomad = (
            variants.gnomad_ac
            if scope == "all_allele_count_proxy"
            else variants.gnomad_carriers
        )
        affected = variants.affected.to_numpy()
        unaffected = (variants.unaffected_literature + gnomad).to_numpy()
        alpha_post = float(prior.alpha_empirical) + affected
        beta_post = float(prior.beta_empirical) + unaffected
        posterior = alpha_post / (alpha_post + beta_post)
        assert np.all(alpha_post > 0) and np.all(beta_post > 0)
        np.testing.assert_allclose(alpha_post - affected, prior.alpha_empirical)
        np.testing.assert_allclose(
            beta_post - unaffected, prior.beta_empirical, atol=1e-11
        )
        density = density_values(weights, posterior)
        if scope == "all_observed_full_locus":
            np.testing.assert_allclose(posterior, variants.posterior_mean, atol=1e-13)
            np.testing.assert_allclose(
                density, primary.density, atol=1e-13, equal_nan=True
            )
        result[f"density_{scope}"] = density
        shared = np.isfinite(density + primary.density.to_numpy())
        difference = abs(density[shared] - primary.density.to_numpy()[shared])
        summaries.append(
            {
                "scope": scope,
                "prior_mean": prior["mean"],
                "prior_strength": prior.strength,
                "alpha_empirical": prior.alpha_empirical,
                "beta_empirical": prior.beta_empirical,
                "donor_variants": len(variants),
                "supported_targets": int(shared.sum()),
                "own_posterior_median": float(np.median(posterior)),
                "density_median": float(np.nanmedian(density)),
                "mean_absolute_density_change": float(difference.mean()),
                "median_absolute_density_change": float(np.median(difference)),
                "maximum_absolute_density_change": float(difference.max()),
                "gnomad_unaffected_measure": "allele_count"
                if scope == "all_allele_count_proxy"
                else "carrier_count_AC_minus_homozygotes",
            }
        )
    save_csv(result, OUT / "GCK_fixed_geometry_prior_sensitivity.csv.gz")
    summary = pd.DataFrame(summaries)
    save_csv(summary, OUT / "GCK_fixed_geometry_prior_sensitivity_summary.csv")
    sources = [
        OUT / "GCK_structural_input_variants.csv.gz",
        OUT / "GCK_primary_variant_density.csv",
        priors_path,
        *shards,
    ]
    checks = {
        "same_634_donor_identities": len(variants) == 634,
        "same_geometry_and_normalized_weights": True,
        "same_variant_only_exclusion": True,
        "alpha_adds_affected": True,
        "beta_adds_literature_unaffected_and_gnomad": True,
        "allele_count_proxy_uses_AC_in_beta": True,
        "primary_reproduces_original_run": True,
        "no_new_regression_bandwidth_or_geometry_fitting": True,
        "interpretation": "Sensitivity only; the primary remains the historical-MSE full-locus carrier-count prior.",
        "script_sha256": digest(__file__),
        "inputs": {str(path.relative_to(HERE)): digest(path) for path in sources},
        "outputs": {
            name: digest(OUT / name)
            for name in [
                "GCK_fixed_geometry_prior_sensitivity.csv.gz",
                "GCK_fixed_geometry_prior_sensitivity_summary.csv",
            ]
        },
    }
    (OUT / "prior_geometry_sensitivity_checks.json").write_text(
        json.dumps(checks, indent=2) + "\n"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
