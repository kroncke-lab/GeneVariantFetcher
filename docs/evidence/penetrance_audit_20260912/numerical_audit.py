#!/usr/bin/env python3
"""Reproduce the 2026-09-12 numerical audit of frozen protocol posteriors.

Run from the GVF checkout with .venv/bin/python and this script's path.
Only reads source data; writes numerical_summary.json and strata.csv beside itself.
Minimal feature-join snapshots beside the script make the audit self-contained.
If missing, these are created from fit_summary.json's recorded input path;
--feature-root can point to another directory containing GENE_protocol subfolders.
The no-gnomAD calculation holds fitted priors fixed: it is a diagnostic, not a refit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

GENES = ("HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1")
MODELS = ("base", "class", "insilico", "clinvar", "insilico_clinvar")
POSTERIOR = "insilico_post_mean"
PRIOR = "insilico_prior_mean"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def feature_snapshot(
    gene: str, source: Path, keys: pd.Series
) -> tuple[pd.DataFrame | None, dict]:
    """Read the frozen minimal snapshot, creating it only when not yet preserved."""
    folder = Path(__file__).resolve().parent / "feature_join"
    snapshot = folder / f"{gene}.csv"
    manifest = folder / f"{gene}.source.json"
    if not snapshot.is_file():
        if not source.is_file():
            return None, {"unavailable_source": str(source)}
        full = pd.read_csv(source)
        assert full.key.is_unique, (gene, "duplicate original feature keys")
        columns = [
            "key",
            "allele_match",
            "clinvar_simple",
            "n_alleles",
            "clinvar_alleles_with_record",
            "single_proband_rows",
            "single_proband_share",
            "gnomad_ac",
            "gnomad_an",
            "gnomad_dataset",
        ]
        minimal = full.loc[full.key.isin(keys), columns].copy()
        assert len(minimal) == len(keys), (gene, "incomplete feature snapshot")
        folder.mkdir(parents=True, exist_ok=True)
        minimal.to_csv(snapshot, index=False)
        provenance = {
            "original_source": str(source),
            "original_source_sha256": digest(source),
            "snapshot": str(snapshot.relative_to(Path(__file__).resolve().parent)),
            "snapshot_sha256": digest(snapshot),
            "rows": len(minimal),
            "columns": columns,
        }
        manifest.write_text(json.dumps(provenance, indent=2) + "\n")
    provenance = json.loads(manifest.read_text())
    assert digest(snapshot) == provenance["snapshot_sha256"], (
        gene,
        "snapshot hash mismatch",
    )
    return pd.read_csv(snapshot), provenance


def stats(d: pd.DataFrame) -> dict:
    x = d[POSTERIOR]
    return {
        "variants": len(d),
        "posterior_median": float(x.median()) if len(d) else None,
        "below_0_10": int((x < 0.10).sum()),
        "at_or_above_0_50": int((x >= 0.50).sum()),
        "at_or_above_0_90": int((x >= 0.90).sum()),
        "gnomad_zero": int((d.gnomad_added == 0).sum()),
        "gnomad_positive": int((d.gnomad_added > 0).sum()),
        "gnomad_count_sum": float(d.gnomad_added.sum()),
        "affected_sum": float(d.affected.sum()),
        "literature_carrier_sum": float(d.literature_n.sum()),
        "literature_n_median": float(d.literature_n.median()) if len(d) else None,
        "n_one": int((d.n == 1).sum()),
        "affected_one_n_one": int(((d.affected == 1) & (d.n == 1)).sum()),
        "n_at_most_5": int((d.n <= 5).sum()),
        "literature_n_at_most_5": int((d.literature_n <= 5).sum()),
        "all_literature_carriers_affected": int((d.affected == d.literature_n).sum()),
        "zero_affected": int((d.affected == 0).sum()),
        "truncating": int(d.is_truncating.sum()),
        "missing_alphamissense": int(d.alphamissense.isna().sum()),
        "clinvar_none": int((d.clinvar_simple == "none").sum()),
        "prior_median": float(d[PRIOR].median()) if len(d) else None,
        "prior_weight_median": float(d.prior_weight.median()) if len(d) else None,
        "posterior_interval_width_median": float(
            (d.insilico_post_q975 - d.insilico_post_q025).median()
        )
        if len(d)
        else None,
    }


def records(d: pd.DataFrame, limit: int = 5) -> list[dict]:
    cols = [
        "key",
        "vclass",
        "is_truncating",
        "clinvar_simple",
        "alphamissense",
        "affected",
        "literature_n",
        "gnomad_added",
        "n",
        "n_papers",
        PRIOR,
        POSTERIOR,
        "insilico_post_q025",
        "insilico_post_q975",
        "pmids",
    ]
    return json.loads(d.loc[:, cols].head(limit).to_json(orient="records"))


def auc_rank(scores: pd.Series, labels: pd.Series) -> float | None:
    """Mann-Whitney AUC, with average ranks for tied ClinVar scores."""
    n_positive = int(labels.sum())
    n_negative = len(labels) - n_positive
    if not n_positive or not n_negative:
        return None
    ranks = scores.rank(method="average")
    return float(
        (ranks[labels].sum() - n_positive * (n_positive + 1) / 2)
        / (n_positive * n_negative)
    )


def same_row_auc(d: pd.DataFrame, source_dir: Path) -> dict:
    path = source_dir / "oof_predictions.csv"
    oof = pd.read_csv(path)
    assert oof.key.is_unique, (path, "duplicate OOF keys")
    paired = d.merge(
        oof[
            ["key", "affected", "n", "observed", "clinvar_simple", "insilico_oof_mean"]
        ],
        on="key",
        how="outer",
        suffixes=("", "_oof"),
        validate="one_to_one",
        indicator=True,
    )
    assert (paired._merge == "both").all(), (path, "different variant universes")
    for col in ("affected", "n", "observed"):
        assert np.allclose(paired[col], paired[col + "_oof"]), (path, col)
    assert (paired.clinvar_simple == paired.clinvar_simple_oof).all(), path
    ordinal = paired.clinvar_simple.map(
        {
            "P/LP": 2.0,
            "Conflicting": 1.0,
            "VUS": 1.0,
            "none": 1.0,
            "Other": 1.0,
            "B/LB": 0.0,
        }
    )
    mask = (
        (paired.n >= 5)
        & (paired.clinvar_simple != "none")
        & np.isfinite(ordinal)
        & np.isfinite(paired.insilico_oof_mean)
        & np.isfinite(paired.observed)
    )
    labels = paired.loc[mask, "observed"] >= 0.5
    model_auc = auc_rank(paired.loc[mask, "insilico_oof_mean"], labels)
    clinvar_auc = auc_rank(ordinal[mask], labels)
    old_path = source_dir / "classification_metrics.csv"
    old = pd.read_csv(old_path)
    old = old[
        (old.high_threshold == 0.5)
        & (old.min_carriers == 5)
        & old.scorer.isin(["clinvar_ordinal", "insilico_oof_prior"])
    ]
    return {
        "note": "Point AUC on identical rows; no significance claim. The outcome is the synthetic observed affected/(literature+gnomAD) ratio, not validated population penetrance.",
        "score": "insilico_oof_mean (prior prediction, excluding target outcome update)",
        "label": "observed >= 0.5; n >= 5; finite scores; ClinVar class != none",
        "oof_source_sha256": digest(path),
        "original_metrics_sha256": digest(old_path),
        "variants": int(mask.sum()),
        "high": int(labels.sum()),
        "insilico_oof_prior_auc": model_auc,
        "clinvar_ordinal_auc": clinvar_auc,
        "difference_model_minus_clinvar": model_auc - clinvar_auc
        if model_auc is not None and clinvar_auc is not None
        else None,
        "matched_keys_sha256": hashlib.sha256(
            "\n".join(sorted(paired.loc[mask, "key"])).encode()
        ).hexdigest(),
        "original_unmatched_comparison": json.loads(
            old[["scorer", "n_variants", "n_high", "auc"]].to_json(orient="records")
        ),
    }


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
    parser.add_argument("--feature-root", type=Path)
    args = parser.parse_args()
    out = {
        "audit_date": "2026-09-12",
        "posterior_column_reproducing_user_table": POSTERIOR,
        "unit": "one retained variant key; gnomAD sum is allele/carrier observations, not unique people",
        "fraction_note": "All shares use retained variant rows, never summed carrier counts.",
        "clinvar_none_note": "Missing joined classification is not proof no ClinVar record exists.",
        "fixed_prior_diagnostic_note": (
            "(affected + strength*prior_mean)/(n+strength) is the analytic mean of the "
            "documented Beta update. Saved posteriors are sampled and differ slightly. "
            "Removing gnomAD below holds fitted priors fixed; it does not remove indirect "
            "gnomAD influence on the fitted prior and is not a population-risk estimate."
        ),
        "genes": {},
    }
    strata = []
    for gene in GENES:
        source_dir = args.source_root / f"{gene}_protocol"
        source = source_dir / "posterior_variants.csv"
        fit_path = source_dir / "fit_summary.json"
        fit = json.loads(fit_path.read_text())
        d = pd.read_csv(source)
        assert d.key.is_unique, (gene, "duplicate keys")
        assert np.allclose(d.n, d.literature_n + d.gnomad_added), gene
        assert (d.affected >= 0).all() and (d.affected <= d.literature_n).all(), gene
        strength = float(fit["prior_strength"])
        d["prior_weight"] = strength / (strength + d.n)
        analytic = (d.affected + strength * d[PRIOR]) / (d.n + strength)
        no_gnomad = (d.affected + strength * d[PRIOR]) / (d.literature_n + strength)
        summary = {
            "source_csv": str(source.resolve()),
            "source_sha256": digest(source),
            "fit_summary_sha256": digest(fit_path),
            "prior_strength": strength,
            "pooled_rate": fit["pooled_rate"],
            "all": stats(d),
            "histogram_table": {
                "variants": len(d),
                "median": float(d[POSTERIOR].median()),
                "below_0_10_pct": float((d[POSTERIOR] < 0.1).mean() * 100),
                "at_or_above_0_50_pct": float((d[POSTERIOR] >= 0.5).mean() * 100),
                "at_or_above_0_90_pct": float((d[POSTERIOR] >= 0.9).mean() * 100),
            },
            "models": {},
            "bands": {},
            "missing_feature_groups": [],
            "gnomad_concentration": {},
            "fixed_prior_diagnostic": {
                "analytic_vs_sampled_mean_abs_difference_max": float(
                    (analytic - d[POSTERIOR]).abs().max()
                ),
                "analytic_vs_sampled_mean_abs_difference_median": float(
                    (analytic - d[POSTERIOR]).abs().median()
                ),
                "saved_mean_below_0_10": int((d[POSTERIOR] < 0.1).sum()),
                "analytic_mean_below_0_10": int((analytic < 0.1).sum()),
                "without_direct_gnomad_mean_below_0_10": int((no_gnomad < 0.1).sum()),
                "analytic_low_crosses_above_0_10_without_gnomad": int(
                    ((analytic < 0.1) & (no_gnomad >= 0.1)).sum()
                ),
            },
            "examples": {},
            "same_row_auc": same_row_auc(d, source_dir),
        }
        synonymous_mismatch = d.loc[
            d.key.str.endswith("=") & (d.vclass != "synonymous")
        ]
        truncating_with_am = d.loc[d.is_truncating & d.alphamissense.notna()]
        summary["classification_consistency"] = {
            "note": "Flagged internal key/consequence or consequence/feature inconsistencies require allele-level adjudication; this check does not automatically relabel variants.",
            "synonymous_key_non_synonymous_vclass": len(synonymous_mismatch),
            "synonymous_mismatch_gnomad_count_sum": float(
                synonymous_mismatch.gnomad_added.sum()
            ),
            "synonymous_mismatch_rows": records(
                synonymous_mismatch, len(synonymous_mismatch)
            ),
            "truncating_with_alphamissense": len(truncating_with_am),
            "truncating_with_alphamissense_rows": records(
                truncating_with_am, len(truncating_with_am)
            ),
        }
        for model in MODELS:
            x = d[f"{model}_post_mean"]
            summary["models"][model] = {
                "median": float(x.median()),
                "below_0_10_pct": float((x < 0.1).mean() * 100),
                "at_or_above_0_50_pct": float((x >= 0.5).mean() * 100),
                "at_or_above_0_90_pct": float((x >= 0.9).mean() * 100),
            }
        masks = {
            "all": pd.Series(True, index=d.index),
            "below_0_10": d[POSTERIOR] < 0.1,
            "mode_0_09_to_0_11_inclusive": d[POSTERIOR].between(0.09, 0.11),
            "mode_0_30_to_0_36_inclusive": d[POSTERIOR].between(0.30, 0.36),
            "at_or_above_0_90": d[POSTERIOR] >= 0.9,
            "clinvar_none_at_or_above_0_90": (d.clinvar_simple == "none")
            & (d[POSTERIOR] >= 0.9),
            "below_0_10_gnomad_zero": (d[POSTERIOR] < 0.1) & (d.gnomad_added == 0),
            "mode_0_09_to_0_11_n_one": d[POSTERIOR].between(0.09, 0.11) & (d.n == 1),
            "mode_0_30_to_0_36_n_one": d[POSTERIOR].between(0.30, 0.36) & (d.n == 1),
            "gnomad_zero": d.gnomad_added == 0,
            "gnomad_positive": d.gnomad_added > 0,
            "n_one": d.n == 1,
            "literature_n_at_least_10": d.literature_n >= 10,
            "literature_has_unaffected": d.affected < d.literature_n,
        }
        for name, mask in masks.items():
            s = d.loc[mask]
            summary["bands"][name] = {
                **stats(s),
                "clinvar_counts": s.clinvar_simple.value_counts().to_dict(),
                "variant_class_counts": s.vclass.value_counts().to_dict(),
            }
            strata.append(
                {"gene": gene, "dimension": "band", "stratum": name, **stats(s)}
            )
        for column in ("clinvar_simple", "vclass", "is_truncating"):
            for value, s in d.groupby(column, dropna=False):
                strata.append(
                    {
                        "gene": gene,
                        "dimension": column,
                        "stratum": str(value),
                        **stats(s),
                    }
                )
        for (truncating, missing), s in d.groupby(
            ["is_truncating", d.alphamissense.isna()]
        ):
            value = f"truncating={truncating};alphamissense_missing={missing}"
            group = {
                "truncating": bool(truncating),
                "alphamissense_missing": bool(missing),
                **stats(s),
            }
            summary["missing_feature_groups"].append(group)
            strata.append(
                {
                    "gene": gene,
                    "dimension": "feature_status",
                    "stratum": value,
                    **stats(s),
                }
            )
        top = d.sort_values("gnomad_added", ascending=False)
        total = float(d.gnomad_added.sum())
        summary["gnomad_concentration"] = {
            "sum": total,
            "top_5_share": float(top.gnomad_added.head(5).sum() / total)
            if total
            else None,
            "top_10_share": float(top.gnomad_added.head(10).sum() / total)
            if total
            else None,
            "common_af_at_least_0_001_variants": int((d.gnomad_af_max >= 0.001).sum()),
            "top_5_examples": records(top),
        }
        for name in (
            "mode_0_09_to_0_11_inclusive",
            "mode_0_30_to_0_36_inclusive",
            "at_or_above_0_90",
        ):
            summary["examples"][name] = records(d.loc[masks[name] & (d.n == 1)])
        source_feature = (
            args.feature_root / f"{gene}_protocol/variants_features.csv"
            if args.feature_root
            else Path(fit["input"])
        )
        f, provenance = feature_snapshot(gene, source_feature, d.key)
        if f is not None:
            assert f.key.is_unique, (gene, "duplicate feature keys")
            fields = [
                "key",
                "allele_match",
                "clinvar_alleles_with_record",
                "single_proband_rows",
                "single_proband_share",
            ]
            j = d.merge(
                f[fields], on="key", how="left", validate="one_to_one", indicator=True
            )
            assert (j._merge == "both").all(), (gene, "unmatched source feature rows")
            none = j.loc[j.clinvar_simple == "none"]
            summary["optional_feature_join_audit"] = {
                **provenance,
                "clinvar_none_rows": len(none),
                "clinvar_none_with_no_warehouse_allele_match": int(
                    (none.allele_match == "none").sum()
                ),
                "clinvar_none_with_warehouse_allele_match": int(
                    (none.allele_match != "none").sum()
                ),
                "clinvar_none_allele_match_counts": none.allele_match.value_counts().to_dict(),
                "n_one_rows_with_single_proband_evidence": int(
                    ((j.n == 1) & (j.single_proband_rows > 0)).sum()
                ),
                "n_one_rows_without_single_proband_evidence": int(
                    ((j.n == 1) & ~(j.single_proband_rows > 0)).sum()
                ),
                "affected_one_n_one_rows_with_single_proband_evidence": int(
                    ((j.n == 1) & (j.affected == 1) & (j.single_proband_rows > 0)).sum()
                ),
                "affected_one_n_one_rows_without_single_proband_evidence": int(
                    (
                        (j.n == 1) & (j.affected == 1) & ~(j.single_proband_rows > 0)
                    ).sum()
                ),
                "all_rows_with_single_proband_share_one": int(
                    (j.single_proband_share == 1).sum()
                ),
                "all_rows_with_any_single_proband_evidence": int(
                    (j.single_proband_rows > 0).sum()
                ),
            }
            for name, mask in masks.items():
                s = j.loc[mask]
                summary["bands"][name]["single_proband_share_one"] = int(
                    (s.single_proband_share == 1).sum()
                )
        else:
            summary["optional_feature_join_audit"] = {
                "unavailable_source": str(source_feature)
            }
        out["genes"][gene] = summary
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "numerical_summary.json").write_text(
        json.dumps(out, indent=2, allow_nan=False) + "\n"
    )
    pd.DataFrame(strata).to_csv(args.output_dir / "strata.csv", index=False)
    print(
        pd.DataFrame(
            {g: s["histogram_table"] for g, s in out["genes"].items()}
        ).T.to_string()
    )


if __name__ == "__main__":
    main()
