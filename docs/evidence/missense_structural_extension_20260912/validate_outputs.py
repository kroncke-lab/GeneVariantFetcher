"""Independently reconstruct saved structural means, variances and scores."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import betaln, gammaln


HERE = Path(__file__).resolve().parent


def read_parts(folder, stem):
    return pd.concat(
        [pd.read_csv(p) for p in sorted(folder.glob(f"{stem}.part*.csv.gz"))],
        ignore_index=True,
    )


def validate(gene):
    folder = HERE / "analysis" / gene
    checks = json.loads((folder / "checks.json").read_text())
    for relative, expected in checks["inputs"].items():
        assert (
            hashlib.sha256((HERE.parent / relative).read_bytes()).hexdigest()
            == expected
        ), relative
    data = read_parts(folder, "input_variants").set_index("variant_id")
    primary = read_parts(folder, "primary_density").set_index("variant_id")
    weights = read_parts(folder, "primary_weights")
    assert not weights.duplicated(["target_id", "donor_id"]).any()
    assert weights.target_id.ne(weights.donor_id).all()
    assert weights.normalized_weight.gt(0).all()
    assert set(data.index) == set(primary.index)
    assert data.variant_type.eq("missense").all()
    assert data.canonical_wt_status.eq("match").all()
    np.testing.assert_allclose(
        data.posterior_alpha, data.alpha_empirical + data.affected
    )
    np.testing.assert_allclose(
        data.posterior_beta,
        data.beta_empirical + data.unaffected_literature + data.gnomad_carriers,
    )
    a, b = data.posterior_alpha, data.posterior_beta
    mean = a / (a + b)
    variance = a * b / ((a + b) ** 2 * (a + b + 1))
    weights["weighted_mean"] = weights.normalized_weight * weights.donor_id.map(mean)
    weights["weighted_variance"] = weights.normalized_weight**2 * weights.donor_id.map(
        variance
    )
    reconstructed = weights.groupby("target_id")[
        ["normalized_weight", "weighted_mean", "weighted_variance"]
    ].sum()
    supported = primary.loc[primary.density.notna()]
    assert set(supported.index) == set(reconstructed.index)
    np.testing.assert_allclose(reconstructed.normalized_weight, 1, atol=1e-12)
    np.testing.assert_allclose(
        reconstructed.weighted_mean,
        primary.loc[reconstructed.index, "density"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        reconstructed.weighted_variance,
        primary.loc[reconstructed.index, "density_conditional_variance"],
        atol=1e-12,
    )
    assert primary.loc[primary.density.isna(), "density_lower_95"].isna().all()
    weights["target_pos"] = weights.target_id.map(data.canonical_pos)
    weights["donor_pos"] = weights.donor_id.map(data.canonical_pos)
    same_residue = int(weights.target_pos.eq(weights.donor_pos).sum())
    predictions = read_parts(folder, "loo_predictions")
    metrics = pd.read_csv(folder / "loo_metrics.csv")
    for (cohort, model), rows in predictions.groupby(["cohort", "model"]):
        expected = metrics.loc[
            metrics.cohort.eq(cohort) & metrics.model.eq(model)
        ].iloc[0]
        p = rows.prediction.to_numpy()
        own = rows.variant_id.map(mean).to_numpy()
        observed = rows.affected.to_numpy() / rows.n.to_numpy()
        np.testing.assert_allclose(
            np.mean(abs(p - own)), expected.mae_empirical_posterior, atol=1e-12
        )
        np.testing.assert_allclose(
            np.mean((p - observed) ** 2), expected.mse_observed_fraction, atol=1e-12
        )
        strength = checks["prior"]["strength"]
        alpha, beta = p * strength, (1 - p) * strength
        k, n = rows.affected.to_numpy(), rows.n.to_numpy()
        ll = (
            gammaln(n + 1)
            - gammaln(k + 1)
            - gammaln(n - k + 1)
            + betaln(k + alpha, n - k + beta)
            - betaln(alpha, beta)
        )
        np.testing.assert_allclose(
            -ll.mean(), expected.mean_beta_binomial_nll, atol=1e-10
        )
    for cohort, rows in predictions.groupby("cohort"):
        cohorts = [set(group.variant_id) for _, group in rows.groupby("model")]
        assert all(group == cohorts[0] for group in cohorts)
    for path in folder.iterdir():
        assert path.stat().st_size < 1_200_000, path.name
        if path.suffix == ".gz":
            assert b"\r\n" not in gzip.decompress(path.read_bytes())
    return {
        "gene": gene,
        "variants": len(data),
        "supported": len(supported),
        "positive_donor_pairs": len(weights),
        "same_residue_distinct_variant_pairs": same_residue,
        "means_and_variances_reconstructed": True,
        "metrics_independently_recomputed": True,
        "input_hashes_verified": len(checks["inputs"]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genes", nargs="+", default=["HNF1A", "LDLR", "KCNQ1", "BRCA2"]
    )
    args = parser.parse_args()
    report = [validate(gene) for gene in args.genes]
    destination = (
        HERE
        / "analysis"
        / (
            "validation.json"
            if len(args.genes) == 4
            else f"validation_{'_'.join(args.genes)}.json"
        )
    )
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
