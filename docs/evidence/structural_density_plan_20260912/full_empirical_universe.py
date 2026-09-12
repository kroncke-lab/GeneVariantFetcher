"""Reproduce empirical-prior preflight from the adjacent committed count snapshot.

The primary calculation uses the historical weight and MSE formula on exact A/n
fractions. It fits one shared prior per gene on every count-bearing variant,
then computes empirical Beta posteriors. No structural fit is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
PRIMARY = "historical_formula_raw_fraction"
METHODS = [
    PRIMARY,
    "historical_formula_endpoint_clipped",
    "historical_weight_normalized_mse",
    "carrier_weighted",
    "variant_equal_weight",
]


def fit_moments(affected, total, method):
    """Return shared empirical Beta parameters and diagnostic moments."""
    fraction = affected / total
    if method.startswith("historical"):
        weights = 1 - 1 / (total + 0.01)
    elif method == "carrier_weighted":
        weights = total
    else:
        weights = np.ones(len(total))
    values = (
        np.clip(fraction, 0.0005, 0.9995)
        if method == "historical_formula_endpoint_clipped"
        else fraction
    )
    mean = np.average(values, weights=weights)
    residual = values - mean
    variance = (
        (weights * residual**2).mean()
        if method.startswith("historical_formula")
        else np.average(residual**2, weights=weights)
    )
    mae = np.average(abs(residual), weights=weights)
    if not (0 < mean < 1 and 0 < variance < mean * (1 - mean)):
        raise ValueError(f"Invalid Beta moments: mean={mean}, variance={variance}")
    strength = mean * (1 - mean) / variance - 1
    alpha, beta = mean * strength, (1 - mean) * strength
    posterior = (alpha + affected) / (strength + total)
    if not np.all((posterior > 0) & (posterior < 1)):
        raise ValueError("Empirical posterior lies outside (0, 1)")
    return {
        "mean": mean,
        "variance": variance,
        "normalized_mae": mae,
        "alpha_empirical": alpha,
        "beta_empirical": beta,
        "strength": strength,
        "posterior_median": np.median(posterior),
        "posterior_lt_01_pct": 100 * np.mean(posterior < 0.1),
        "mean_weight": weights.mean(),
        "MAE_as_variance_valid": bool(mae < mean * (1 - mean)),
    }, posterior


def verify_moments():
    """Verify algebra and endpoint observations without adding pseudocounts."""
    alpha, beta = 2.0, 3.0
    strength = alpha + beta
    mean = alpha / strength
    variance = alpha * beta / (strength**2 * (strength + 1))
    recovered_strength = mean * (1 - mean) / variance - 1
    assert np.allclose(
        [mean * recovered_strength, (1 - mean) * recovered_strength],
        [alpha, beta],
    )
    # Adding 1 affected to Beta(2,3) gives Beta(3,3); n=0 leaves it unchanged.
    assert np.isclose((alpha + 1) / (strength + 1), 0.5)
    assert np.isclose((alpha + 0) / (strength + 0), mean)
    affected = np.array([0.0, 1.0, 4.0])
    total = np.array([1.0, 1.0, 5.0])
    result, _ = fit_moments(affected, total, PRIMARY)
    w = 1 - 1 / (total + 0.01)
    expected = np.average(affected / total, weights=w)
    assert np.isclose(result["mean"], expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((HERE / "source_count_provenance.json").read_text())
    source = HERE / provenance["snapshot_file"]
    snapshot_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    if snapshot_hash != provenance["snapshot_sha256"]:
        raise ValueError("Source-count snapshot SHA-256 differs from provenance")
    full = pd.read_csv(source)
    assert len(full) == provenance["snapshot_rows"]
    assert not full.duplicated(["gene", "key"]).any()
    assert full.key.notna().all()
    count_fields = ["affected", "literature_n", "unaffected", "gnomad_added"]
    assert np.isfinite(full[count_fields].to_numpy()).all()
    assert (full[count_fields] >= 0).all().all()
    assert np.allclose(full.affected + full.unaffected, full.literature_n)
    assert full.is_histogram.isin([True, False]).all()
    assert full.gene.unique().tolist() == GENES
    full["unaffected_total"] = full.unaffected + full.gnomad_added
    full["n"] = full.literature_n + full.gnomad_added
    assert np.allclose(full.affected + full.unaffected_total, full.n)
    assert np.all(full.n > 0)
    rows, totals, posteriors = [], [], []
    verify_moments()
    for gene in GENES:
        frame = full[full.gene == gene]
        histogram = frame[frame.is_histogram]
        restored = frame[~frame.is_histogram]
        totals.append(
            {
                "gene": gene,
                "full_variants": len(frame),
                "histogram_variants": len(histogram),
                "restored_rows": len(restored),
                "restored_classes": ";".join(sorted(restored.vclass.unique())),
                "full_affected": frame.affected.sum(),
                "full_literature_unaffected": frame.unaffected.sum(),
                "full_gnomad_unaffected": frame.gnomad_added.sum(),
                "full_unaffected": frame.unaffected_total.sum(),
                "restored_affected": restored.affected.sum(),
                "restored_literature_unaffected": restored.unaffected.sum(),
                "restored_gnomad_unaffected": restored.gnomad_added.sum(),
                "full_missense": (frame.vclass == "missense").sum(),
                "full_all_have_literature_carriers": bool(
                    (frame.literature_n > 0).all()
                ),
                "full_pooled_carrier_fraction": frame.affected.sum() / frame.n.sum(),
            }
        )
        subsets = [
            ("full_count_dataset", frame),
            ("histogram_subset", histogram),
            ("missense_only", frame[frame.vclass == "missense"]),
        ]
        for subset_name, subset in subsets:
            affected = subset.affected.to_numpy()
            total = subset.n.to_numpy()
            for method in METHODS:
                moments, posterior = fit_moments(affected, total, method)
                rows.append(
                    {
                        "gene": gene,
                        "universe": subset_name,
                        "method": method,
                        "n_variants": len(subset),
                        **moments,
                    }
                )
                if subset_name == "full_count_dataset" and method == PRIMARY:
                    output = subset[
                        [
                            "gene",
                            "key",
                            "vclass",
                            "aa_pos",
                            "aa_ref",
                            "aa_alt",
                            "affected",
                            "literature_n",
                            "unaffected",
                            "gnomad_added",
                            "unaffected_total",
                            "n",
                        ]
                    ].copy()
                    alpha = moments["alpha_empirical"]
                    beta = moments["beta_empirical"]
                    strength = moments["strength"]
                    output["alpha_empirical"] = alpha
                    output["beta_empirical"] = beta
                    output["alpha_posterior_empirical"] = alpha + affected
                    output["beta_posterior_empirical"] = beta + total - affected
                    output["posterior_empirical_mean"] = posterior
                    output["posterior_empirical_variance"] = (
                        posterior * (1 - posterior) / (strength + total + 1)
                    )
                    posteriors.append(output)
    pd.DataFrame(totals).to_csv(out / "full_universe_totals.csv", index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "full_universe_moments.csv", index=False)
    pd.concat(posteriors, ignore_index=True).to_csv(
        out / "full_universe_empirical_posteriors.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    checks = {
        "source_hashes": provenance["external_source_hashes"],
        "snapshot_sha256": snapshot_hash,
        "snapshot_rows": len(full),
        "count_partitions": True,
        "one_row_per_gene_key": True,
        "all_full_rows_positive_n": True,
        "all_candidate_moments_valid": True,
        "beta_moment_roundtrip": True,
        "exact_zero_and_one_fractions_preserved": True,
        "only_committed_snapshot_and_provenance_read": True,
    }
    (out / "full_universe_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(
        summary[
            (summary.universe == "full_count_dataset") & (summary.method == PRIMARY)
        ][["gene", "mean", "alpha_empirical", "beta_empirical", "strength"]].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
