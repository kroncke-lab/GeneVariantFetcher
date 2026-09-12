"""Independent count/weight audit of frozen missense and nonsense prior units.

This recipe reads the preceding population-inclusive union without changing it.
It independently implements the requested empirical-moment convention and two
diagnostic alternatives. It does not select or tune a replacement prior.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution

HERE = Path(__file__).resolve().parent
SOURCE = (
    HERE.parents[1] / "population_inclusive_penetrance_20260912/analysis/union_counts"
)
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
CLASSES = {"missense": {"missense"}, "nonsense": {"nonsense", "stop_gained"}}
METHODS = [
    "historical_weighted_mean_mse",
    "normalized_weighted_mse",
    "equal_variant_mean_mse",
]
BINS = [-0.5, 0.5, 1.5, 2.5, 4.5, 9.5, 19.5, 99.5, np.inf]
LABELS = ["0", "1", "2", "3-4", "5-9", "10-19", "20-99", "100+"]


def moments(affected, n, method):
    """Calculate from observed fractions; do not clip endpoints or prior shape."""
    y = affected / n
    w = 1 - 1 / (n + 0.01)
    if method == "equal_variant_mean_mse":
        mean = y.mean()
        variance = np.mean((y - mean) ** 2)
    else:
        mean = np.average(y, weights=w)
        weighted_ss = np.sum(w * (y - mean) ** 2)
        divisor = len(y) if method == "historical_weighted_mean_mse" else w.sum()
        variance = weighted_ss / divisor
    result = {
        "method": method,
        "mean": mean,
        "variance": variance,
        "mean_historical_weight": w.mean(),
        "n_units": len(y),
    }
    if not (0 < mean < 1 and 0 < variance < mean * (1 - mean)):
        return {**result, "fit_status": "invalid_beta_moments"}
    strength = mean * (1 - mean) / variance - 1
    alpha, beta = mean * strength, (1 - mean) * strength
    result.update(
        {
            "fit_status": "valid",
            "strength": strength,
            "alpha_empirical": alpha,
            "beta_empirical": beta,
            "A0_U1_posterior_mean": alpha / (strength + 1),
            "A1_U0_posterior_mean": (alpha + 1) / (strength + 1),
            "A0_U1_absolute_drop_from_prior": mean / (strength + 1),
            "A0_U1_relative_drop_from_prior": 1 / (strength + 1),
            "A0_U10_posterior_mean": alpha / (strength + 10),
            "prior_weight_after_one_observation": strength / (strength + 1),
        }
    )
    lower, upper = beta_distribution.ppf([0.025, 0.975], alpha, beta + 1)
    result["A0_U1_lower_95"], result["A0_U1_upper_95"] = lower, upper
    return result


def main():
    paths = sorted(SOURCE.glob("*.csv.gz"))
    if not paths:
        raise ValueError("Frozen union shards are missing")
    frame = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    assert len(frame) == 180613 and frame.unit_id.is_unique
    assert np.equal(frame.affected + frame.unaffected, frame.n).all()
    assert np.equal(
        frame.unaffected_literature + frame.gnomad_carriers, frame.unaffected
    ).all()
    assert (frame.n > 0).all()
    assert (frame.loc[frame.origin.eq("population_only"), "affected"] == 0).all()
    inventory = (
        frame.groupby(["gene", "vclass", "canonical_wt_status"], dropna=False)
        .agg(n_units=("unit_id", "size"), n_carrier_observations=("n", "sum"))
        .reset_index()
    )
    inventory["selected_type"] = np.select(
        [
            inventory.vclass.eq("missense") & inventory.canonical_wt_status.eq("match"),
            inventory.vclass.isin(CLASSES["nonsense"])
            & inventory.canonical_wt_status.eq("match"),
        ],
        ["missense", "nonsense"],
        default="excluded_other_class_or_unavailable_WT",
    )
    counts, fits, distributions, influence = [], [], [], []
    selected_ids = set()
    for gene in GENES:
        for variant_type, included in CLASSES.items():
            group = frame.loc[
                frame.gene.eq(gene)
                & frame.vclass.isin(included)
                & frame.canonical_wt_status.eq("match")
            ].copy()
            assert not selected_ids.intersection(group.unit_id)
            selected_ids.update(group.unit_id)
            n, a = group.n.to_numpy(float), group.affected.to_numpy(float)
            y, w = a / n, 1 - 1 / (n + 0.01)
            population = group.origin.eq("population_only").to_numpy()
            has_population = group.n_member_alleles.gt(0).to_numpy()
            singleton = n == 1
            base = {"gene": gene, "variant_type": variant_type}
            counts.append(
                {
                    **base,
                    "n_units": len(group),
                    "population_only_units": int(population.sum()),
                    "population_only_share": population.mean(),
                    "n1_units": int(singleton.sum()),
                    "n1_share": singleton.mean(),
                    "population_only_n1_units": int((population & singleton).sum()),
                    "population_only_n1_share": singleton[population].mean(),
                    "has_population_units": int(has_population.sum()),
                    "one_gnomad_carrier_units": int(
                        ((group.gnomad_carriers == 1) & has_population).sum()
                    ),
                    "one_gnomad_carrier_share_among_has_population": group.loc[
                        has_population, "gnomad_carriers"
                    ]
                    .eq(1)
                    .mean(),
                    "clinical_n1_units": int((~population & singleton).sum()),
                    "A0_U1_units": int(((a == 0) & singleton).sum()),
                    "A1_U0_units": int(((a == 1) & singleton).sum()),
                    "median_n": np.median(n),
                    "q25_n": np.quantile(n, 0.25),
                    "q75_n": np.quantile(n, 0.75),
                    "affected_sum": a.sum(),
                    "unaffected_sum": group.unaffected.sum(),
                    "gnomad_carriers_sum": group.gnomad_carriers.sum(),
                    "n1_historical_weight_share": w[singleton].sum() / w.sum(),
                    "n1_historical_affected_signal_share": (
                        w[singleton] * y[singleton]
                    ).sum()
                    / (w * y).sum(),
                    "pooled_affected_fraction": a.sum() / n.sum(),
                    "equal_variant_fraction_mean": y.mean(),
                }
            )
            current_fits = [moments(a, n, method) for method in METHODS]
            for fit in current_fits:
                fits.append({**base, **fit})
            historical, normalized, _ = current_fits
            if historical["fit_status"] == normalized["fit_status"] == "valid":
                assert np.isclose(
                    normalized["strength"] + 1, (historical["strength"] + 1) * w.mean()
                )
            for origin, mask in [
                ("all", np.ones(len(group), dtype=bool)),
                ("population_only", population),
                ("has_literature", ~population),
            ]:
                subset = group.loc[mask]
                for field in [
                    "n",
                    "affected",
                    "unaffected",
                    "gnomad_carriers",
                    "unaffected_literature",
                ]:
                    buckets = pd.cut(subset[field], BINS, labels=LABELS)
                    frequency = buckets.value_counts(sort=False)
                    for bucket, count in frequency.items():
                        distributions.append(
                            {
                                **base,
                                "origin_scope": origin,
                                "count_field": field,
                                "count_bucket": bucket,
                                "n_units": int(count),
                                "denominator_units": len(subset),
                                "share_units": count / len(subset)
                                if len(subset)
                                else np.nan,
                            }
                        )
            group["n_bucket"] = pd.cut(group.n, BINS, labels=LABELS)
            group["historical_weight"] = w
            group["fraction"] = y
            for bucket, subset in group.groupby("n_bucket", observed=True):
                influence.append(
                    {
                        **base,
                        "n_bucket": bucket,
                        "n_units": len(subset),
                        "unit_share": len(subset) / len(group),
                        "population_only_units": int(
                            subset.origin.eq("population_only").sum()
                        ),
                        "equal_variant_fraction_mean": subset.fraction.mean(),
                        "historical_weight_sum": subset.historical_weight.sum(),
                        "historical_weight_share": subset.historical_weight.sum()
                        / w.sum(),
                        "weighted_affected_fraction_sum": (
                            subset.historical_weight * subset.fraction
                        ).sum(),
                    }
                )
    assert len(selected_ids) == 11300
    for name, rows in [
        ("class_count_summary", counts),
        ("class_moment_comparison", fits),
        ("count_distributions", distributions),
        ("historical_weight_influence", influence),
    ]:
        pd.DataFrame(rows).to_csv(HERE / f"{name}.csv", index=False)
    inventory.to_csv(HERE / "consequence_inventory.csv", index=False)
    pd.DataFrame(
        [
            {
                "source_file": str(path.relative_to(HERE.parents[1])),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in paths
        ]
    ).to_csv(HERE / "frozen_source_hashes.csv", index=False)
    print(
        pd.DataFrame(fits)
        .query('method == "historical_weighted_mean_mse"')[
            [
                "gene",
                "variant_type",
                "n_units",
                "mean",
                "strength",
                "A0_U1_posterior_mean",
                "A0_U1_relative_drop_from_prior",
            ]
        ]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
