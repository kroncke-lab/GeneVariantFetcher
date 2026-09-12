#!/usr/bin/env python3
"""Describe frozen protocol priors and counts under gnomAD-as-unaffected.

No refit or source-data changes. Run with the repository .venv/bin/python.
One row is one retained literature variant key, not one individual.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

GENES = ("HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distribution(x: pd.Series) -> dict:
    x = x.dropna()
    if not len(x):
        return {"n": 0}
    return {
        "n": len(x),
        "min": float(x.min()),
        "q25": float(x.quantile(0.25)),
        "median": float(x.median()),
        "q75": float(x.quantile(0.75)),
        "max": float(x.max()),
        "mean": float(x.mean()),
    }


def rho(x: pd.Series, y: pd.Series) -> dict:
    valid = x.notna() & y.notna()
    a, b = x[valid], y[valid]
    value = None
    if len(a) >= 3 and a.nunique() > 1 and b.nunique() > 1:
        value = float(a.rank(method="average").corr(b.rank(method="average")))
    return {"n": int(valid.sum()), "spearman_rho": value}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "grant_e2e_20260909/analysis",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "date": "2026-09-12",
        "assumption": "Every gnomad_added observation is treated as unaffected, as requested.",
        "universe": "Frozen insilico protocol rows: retained literature variants with usable phenotype counts; population-only variants absent.",
        "unit": "Per-variant rows and variant-carrier observations; counts summed across variants are not unique individuals.",
        "prior": "insilico_prior_mean; feature-conditioned mean fitted to these counts using truncation, AlphaMissense, and AlphaMissense missingness.",
        "pseudo_counts": "Prior equivalent A=S*p0 and U=S*(1-p0), S=10; model weights, not observed people.",
        "mean_formula": "p_mean=(affected+S*p0)/(affected+unaffected_total+S). Saved posterior means use Monte Carlo draws, so differ slightly.",
        "correlation_note": "Descriptive Spearman rank correlation; no independent validation or causal interpretation. AlphaMissense is an input to this prior.",
        "genes": {},
    }
    (
        variants,
        gene_rows,
        counts_rows,
        strength_rows,
        am_rows,
        probability_rows,
        strata_rows,
    ) = ([] for _ in range(7))
    for gene in GENES:
        source = args.source_root / f"{gene}_protocol/posterior_variants.csv"
        fit_path = source.parent / "fit_summary.json"
        fit = json.loads(fit_path.read_text())
        d = pd.read_csv(source)
        strength = float(fit["prior_strength"])
        assert d.key.is_unique and strength == 10.0, gene
        assert np.allclose(d.n, d.literature_n + d.gnomad_added), gene
        assert (d.affected >= 0).all() and (d.affected <= d.literature_n).all(), gene
        assert (d.gnomad_added >= 0).all() and (d.n > 0).all(), gene
        for column in ("affected", "literature_n", "gnomad_added", "n"):
            assert np.allclose(d[column], d[column].round()), (gene, column)
        d["gene"] = gene
        d["unaffected_literature"] = d.literature_n - d.affected
        d["unaffected_total"] = d.unaffected_literature + d.gnomad_added
        d["prior_mean"] = d.insilico_prior_mean
        d["prior_affected"] = strength * d.prior_mean
        d["prior_unaffected"] = strength * (1 - d.prior_mean)
        d["prior_strength"] = strength
        d["prior_weight"] = strength / (strength + d.n)
        d["count_fraction"] = d.affected / d.n
        d["posterior_analytic"] = (d.affected + d.prior_affected) / (d.n + strength)
        d["posterior_mean"] = d.insilico_post_mean
        d["posterior_minus_prior"] = d.posterior_mean - d.prior_mean
        d["feature_stratum"] = np.select(
            [d.is_truncating, ~d.is_truncating & d.alphamissense.notna()],
            ["truncating", "nontruncating_AM_present"],
            default="nontruncating_AM_missing",
        )
        assert np.allclose(d.count_fraction, d.observed), gene
        assert np.allclose(d.affected + d.unaffected_total, d.n), gene
        numeric = [
            "prior_mean",
            "prior_affected",
            "prior_unaffected",
            "prior_weight",
            "affected",
            "unaffected_literature",
            "gnomad_added",
            "unaffected_total",
            "n",
            "count_fraction",
            "posterior_mean",
            "posterior_minus_prior",
        ]
        g = {
            "source_sha256": digest(source),
            "fit_summary_sha256": digest(fit_path),
            "variants": len(d),
            "prior_strength": strength,
            "distributions": {col: distribution(d[col]) for col in numeric},
            "sums": {
                col: float(d[col].sum())
                for col in [
                    "affected",
                    "unaffected_literature",
                    "gnomad_added",
                    "unaffected_total",
                    "n",
                ]
            },
            "pooled_count_fraction": float(d.affected.sum() / d.n.sum()),
            "pooled_literature_fraction": float(
                d.affected.sum() / d.literature_n.sum()
            ),
            "gnomad_positive_variants": int((d.gnomad_added > 0).sum()),
            "unaffected_literature_positive_variants": int(
                (d.unaffected_literature > 0).sum()
            ),
            "zero_affected_variants": int((d.affected == 0).sum()),
            "zero_unaffected_variants": int((d.unaffected_total == 0).sum()),
            "prior_majority_variants": int((d.prior_weight > 0.5).sum()),
            "prior_90pct_variants": int((d.prior_weight >= 0.9).sum()),
            "count_majority_variants": int((d.n > strength).sum()),
            "am_present": int(d.alphamissense.notna().sum()),
            "truncating": int(d.is_truncating.sum()),
            "gnomad_concentration": {
                f"top_{k}_variant_share": float(
                    d.gnomad_added.nlargest(k).sum() / d.gnomad_added.sum()
                )
                for k in (1, 5, 10)
            },
            "largest_gnomad_rows": json.loads(
                d.nlargest(5, "gnomad_added")[
                    [
                        "key",
                        "affected",
                        "unaffected_literature",
                        "gnomad_added",
                        "prior_mean",
                        "posterior_mean",
                    ]
                ].to_json(orient="records")
            ),
            "prior_count_correlations": {
                col: rho(d.prior_mean, d[col])
                for col in [
                    "affected",
                    "unaffected_literature",
                    "gnomad_added",
                    "unaffected_total",
                    "n",
                    "count_fraction",
                ]
            },
            "analytic_saved_abs_diff_max": float(
                (d.posterior_analytic - d.posterior_mean).abs().max()
            ),
        }
        result["genes"][gene] = g
        gene_rows.append(
            {
                "gene": gene,
                "variants": len(d),
                **{col + "_sum": value for col, value in g["sums"].items()},
                **{
                    col + "_" + stat: g["distributions"][col][stat]
                    for col in [
                        "prior_mean",
                        "count_fraction",
                        "posterior_mean",
                        "affected",
                        "unaffected_total",
                        "n",
                        "prior_weight",
                    ]
                    for stat in ["q25", "median", "q75"]
                },
                "pooled_count_fraction": g["pooled_count_fraction"],
                "prior_majority_pct": 100 * g["prior_majority_variants"] / len(d),
                "count_majority_pct": 100 * g["count_majority_variants"] / len(d),
                "gnomad_positive_pct": 100 * g["gnomad_positive_variants"] / len(d),
            }
        )
        for metric in [
            "affected",
            "unaffected_literature",
            "gnomad_added",
            "unaffected_total",
            "n",
        ]:
            categories = pd.cut(
                d[metric],
                [-0.1, 0, 1, 4, 9, np.inf],
                labels=["0", "1", "2-4", "5-9", "10+"],
            )
            for band, n in categories.value_counts(sort=False).items():
                counts_rows.append(
                    {
                        "gene": gene,
                        "metric": metric,
                        "band": band,
                        "variants": int(n),
                        "pct": 100 * n / len(d),
                    }
                )
        weights = {
            "prior >=90% (n<=1)": d.prior_weight >= 0.9,
            "prior >50% to <90% (n=2-9)": (d.prior_weight > 0.5)
            & (d.prior_weight < 0.9),
            "equal weight (n=10)": d.prior_weight == 0.5,
            "counts >50% to <90% (n=11-89)": (d.prior_weight < 0.5)
            & (d.prior_weight > 0.1),
            "counts >=90% (n>=90)": d.prior_weight <= 0.1,
        }
        assert sum(int(mask.sum()) for mask in weights.values()) == len(d)
        for label, mask in weights.items():
            strength_rows.append(
                {
                    "gene": gene,
                    "band": label,
                    "variants": int(mask.sum()),
                    "pct": 100 * mask.mean(),
                }
            )
        am_masks = {
            "all_AM_present": d.alphamissense.notna(),
            "missense_AM_present": (d.vclass == "missense") & d.alphamissense.notna(),
            "nontruncating_AM_present": ~d.is_truncating & d.alphamissense.notna(),
            "truncating_AM_present": d.is_truncating & d.alphamissense.notna(),
        }
        for label, mask in am_masks.items():
            for target in ["prior_mean", "count_fraction", "posterior_mean"]:
                am_rows.append(
                    {
                        "gene": gene,
                        "stratum": label,
                        "target": target,
                        **rho(d.loc[mask, "alphamissense"], d.loc[mask, target]),
                    }
                )
        for label, group in d.groupby("feature_stratum"):
            strata_rows.append(
                {
                    "gene": gene,
                    "stratum": label,
                    "variants": len(group),
                    **{
                        col + "_" + stat: distribution(group[col])[stat]
                        for col in ["prior_mean", "posterior_mean", "count_fraction"]
                        for stat in ["min", "q25", "median", "q75", "max"]
                    },
                }
            )
        for metric in ["prior_mean", "count_fraction", "posterior_mean"]:
            n, bins = np.histogram(d[metric], bins=np.linspace(0, 1, 21))
            for index, count in enumerate(n):
                probability_rows.append(
                    {
                        "gene": gene,
                        "metric": metric,
                        "bin_left": bins[index],
                        "bin_right": bins[index + 1],
                        "variants": int(count),
                        "pct": 100 * count / len(d),
                    }
                )
        variants.append(
            d[
                [
                    "gene",
                    "key",
                    "vclass",
                    "is_truncating",
                    "clinvar_simple",
                    "alphamissense",
                    "feature_stratum",
                    "affected",
                    "unaffected_literature",
                    "gnomad_added",
                    "unaffected_total",
                    "n",
                    "prior_mean",
                    "prior_affected",
                    "prior_unaffected",
                    "prior_strength",
                    "prior_weight",
                    "count_fraction",
                    "posterior_mean",
                    "posterior_analytic",
                    "posterior_minus_prior",
                ]
            ]
        )
    (args.output_dir / "counts_summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    tables = {
        "gene_summary": gene_rows,
        "count_distributions": counts_rows,
        "evidence_strength": strength_rows,
        "alphamissense_comparison": am_rows,
        "probability_distributions": probability_rows,
        "feature_strata": strata_rows,
    }
    for name, rows in tables.items():
        pd.DataFrame(rows).to_csv(args.output_dir / f"{name}.csv", index=False)
    variant_csv = pd.concat(variants, ignore_index=True).to_csv(index=False)
    (args.output_dir / "variant_counts.csv.gz").write_bytes(
        gzip.compress(variant_csv.encode("utf-8"), mtime=0)
    )
    print(
        pd.DataFrame(gene_rows)[
            [
                "gene",
                "variants",
                "prior_mean_median",
                "count_fraction_median",
                "posterior_mean_median",
                "affected_sum",
                "unaffected_literature_sum",
                "gnomad_added_sum",
                "pooled_count_fraction",
                "prior_majority_pct",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
