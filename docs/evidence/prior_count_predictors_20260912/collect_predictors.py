#!/usr/bin/env python3
"""Read-only score extraction using the frozen protocol's existing allele IDs.

Run from GeneVariantFetcher with .venv/bin/python; no API jobs are submitted.
The base predictor columns are restricted to one mapped allele and an
unambiguous version/value. Means/ranges are descriptive sensitivity columns,
not an arbitrary selection of a nucleotide allele for a protein-level key.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

GENES = ("HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1")
MODELS = ("m447", "v100", "p243")
PATHOGENICITY = ["alphamissense", "alphagenome_avi"] + [
    f"gpn_star_{model}_{metric}"
    for model in MODELS
    for metric in ("llr_calibrated", "abs_llr_calibrated")
]
CONSERVATION = [f"gpn_star_{model}_entropy_calibrated" for model in MODELS]
METRICS = PATHOGENICITY + CONSERVATION
GPN_VERSION = "5c799b2ec6aa089f0caa8294ae72adb4510f81ae"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_ids(value):
    if pd.isna(value) or not str(value).strip():
        return []
    return sorted({int(x) for x in str(value).split(";")})


def finite(value):
    return value is not None and np.isfinite(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("/Users/kronckbm/GitRepos/variantFeatures/data/variants.db"),
    )
    parser.add_argument(
        "--protocol-root",
        type=Path,
        default=Path(
            "/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results"
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    frames = []
    inputs = {}
    for gene in GENES:
        path = args.protocol_root / f"{gene}_protocol/variants_features.csv"
        frame = pd.read_csv(path, dtype={"vf_variant_ids": str})
        assert frame.key.is_unique
        frame.insert(0, "gene", gene)
        frames.append(frame)
        inputs[gene] = {"path": str(path), "sha256": sha256(path), "rows": len(frame)}
    features = pd.concat(frames, ignore_index=True)
    features["mapped_ids"] = features.vf_variant_ids.map(parse_ids)
    ids = sorted({vid for group in features.mapped_ids for vid in group})
    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only=ON")
    conn.execute("BEGIN")  # One consistent read snapshot across all tables.
    query_time = datetime.now(timezone.utc).isoformat()
    scores = defaultdict(lambda: defaultdict(list))
    avi_contexts = defaultdict(list)
    variants = {}
    identity_genes = defaultdict(set)
    provenance = defaultdict(Counter)
    for start in range(0, len(ids), 400):
        chunk = ids[start : start + 400]
        marks = ",".join("?" for _ in chunk)
        predmarks = ",".join("?" for _ in PATHOGENICITY)
        query = f"""SELECT * FROM annotations_pathogenicity
                    WHERE variant_id IN ({marks}) AND predictor IN ({predmarks})"""
        for row in conn.execute(query, chunk + PATHOGENICITY):
            row = dict(row)
            scores[row["variant_id"]][row["predictor"]].append(row)
            provenance[row["predictor"]][
                (row["predictor_version"], row["source"], row["fetched_at"])
            ] += 1
        consmarks = ",".join("?" for _ in CONSERVATION)
        query = f"""SELECT * FROM annotations_conservation
                    WHERE variant_id IN ({marks}) AND metric IN ({consmarks})"""
        for row in conn.execute(query, chunk + CONSERVATION):
            row = dict(row)
            row["predictor_version"] = ""  # This table has no version column.
            scores[row["variant_id"]][row["metric"]].append(row)
            provenance[row["metric"]][("", row["source"], row["fetched_at"])] += 1
        query = f"""SELECT * FROM annotations_alphagenome_context
                    WHERE variant_id IN ({marks}) AND scorer='AVI_SCORE'"""
        for row in conn.execute(query, chunk):
            avi_contexts[row["variant_id"]].append(dict(row))
        for row in conn.execute(f"SELECT * FROM variants WHERE id IN ({marks})", chunk):
            variants[row["id"]] = dict(row)
        for row in conn.execute(
            f"SELECT DISTINCT variant_id,gene_symbol FROM variant_consequences WHERE variant_id IN ({marks})",
            chunk,
        ):
            identity_genes[row["variant_id"]].add(row["gene_symbol"])
    schema_version = conn.execute("PRAGMA user_version").fetchone()[0]
    conn.rollback()
    conn.close()

    # Resolve duplicate version/value records only if there is a single identity.
    # Validate AVI against the lossless scalar context, not molecular tracks.
    resolved = defaultdict(dict)
    allele_status = defaultdict(dict)
    avi_context_conflicts = []
    for vid in ids:
        for metric in METRICS:
            records = [r for r in scores[vid][metric] if finite(r["score"])]
            distinct = {(r["predictor_version"], r["score"]) for r in records}
            state = "missing"
            value = np.nan
            if len(distinct) == 1:
                state = "available"
                value = records[0]["score"]
            elif len(distinct) > 1:
                state = "version_or_value_conflict"
            if metric == "alphagenome_avi" and records:
                contexts = avi_contexts[vid]
                context_values = {
                    (r["dataset_version"], r["raw_score"])
                    for r in contexts
                    if finite(r["raw_score"])
                }
                if len(context_values) != 1 or distinct != context_values:
                    state = "context_projection_conflict"
                    value = np.nan
                    avi_context_conflicts.append(vid)
            resolved[vid][metric] = value
            allele_status[vid][metric] = state

    key_rows, allele_rows = [], []
    for record in features.to_dict("records"):
        row = {
            k: record[k]
            for k in ("gene", "key", "vclass", "allele_match", "vf_variant_ids")
        }
        group = record["mapped_ids"]
        row["mapped_allele_count"] = len(group)
        row["all_ids_found"] = all(vid in variants for vid in group)
        row["all_ids_match_gene"] = all(
            record["gene"] in identity_genes[vid] for vid in group
        )
        assert len(group) == record["n_alleles"], (record["gene"], record["key"])
        assert row["all_ids_found"] and row["all_ids_match_gene"], row
        for metric in METRICS:
            values = [
                resolved[vid][metric] for vid in group if finite(resolved[vid][metric])
            ]
            n = len(values)
            row[metric] = values[0] if len(group) == 1 and n == 1 else np.nan
            row[f"{metric}_mean"] = float(np.mean(values)) if n else np.nan
            row[f"{metric}_min"] = min(values) if n else np.nan
            row[f"{metric}_max"] = max(values) if n else np.nan
            row[f"{metric}_n_alleles_scored"] = n
            if not group:
                status = "no_allele_match"
            elif not n:
                conflicts = {allele_status[vid][metric] for vid in group} - {"missing"}
                status = (
                    ";".join(sorted(conflicts)) if conflicts else "missing_annotation"
                )
            elif n < len(group):
                status = "partial_allele_coverage"
            elif len(group) > 1:
                status = "multiple_alleles"
            else:
                status = "single_allele"
            row[f"{metric}_status"] = status
        key_rows.append(row)
        for vid in group:
            allele = {k: record[k] for k in ("gene", "key", "vclass", "allele_match")}
            allele["variant_id"] = vid
            allele.update(
                {
                    k: variants[vid][k]
                    for k in ("chromosome", "position", "ref", "alt", "variant_type")
                }
            )
            allele["key_mapped_allele_count"] = len(group)
            allele.update(resolved[vid])
            allele["alphagenome_context_count"] = len(avi_contexts[vid])
            allele_rows.append(allele)
    result = pd.DataFrame(key_rows)
    assert not result[["gene", "key"]].duplicated().any()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    compressed = {"method": "gzip", "mtime": 0}
    result.to_csv(
        args.output_dir / "predictor_scores.csv.gz", index=False, compression=compressed
    )
    pd.DataFrame(allele_rows).to_csv(
        args.output_dir / "predictor_alleles.csv.gz",
        index=False,
        compression=compressed,
    )
    conflicts = [
        dict(row, metric=metric, resolution_status=allele_status[vid][metric])
        for vid in ids
        for metric in METRICS
        if "conflict" in allele_status[vid][metric]
        for row in scores[vid][metric]
    ]
    pd.DataFrame(conflicts).to_csv(
        args.output_dir / "predictor_version_conflicts.csv.gz",
        index=False,
        compression=compressed,
    )
    coverage = {}
    for gene in GENES:
        subset = result[result.gene == gene]
        coverage[gene] = {
            "feature_keys": len(subset),
            "one_mapped_allele": int((subset.mapped_allele_count == 1).sum()),
            "multiple_mapped_alleles": int((subset.mapped_allele_count > 1).sum()),
            "no_mapped_allele": int((subset.mapped_allele_count == 0).sum()),
            "predictors": {
                metric: {
                    "single_allele_scored": int(subset[metric].notna().sum()),
                    "any_mapped_allele_scored": int(
                        (subset[f"{metric}_n_alleles_scored"] > 0).sum()
                    ),
                    "status": subset[f"{metric}_status"].value_counts().to_dict(),
                }
                for metric in METRICS
            },
        }
    metadata = {
        "queried_at_utc": query_time,
        "db_path": str(args.db),
        "db_schema_version": schema_version,
        "read_only": True,
        "unique_queried_variant_ids": len(ids),
        "output_files": [
            "predictor_scores.csv.gz",
            "predictor_alleles.csv.gz",
            "predictor_version_conflicts.csv.gz",
        ],
        "inputs": inputs,
        "selected_primary_predictors": [
            "alphamissense",
            "gpn_star_m447_llr_calibrated",
            "alphagenome_avi",
        ],
        "name_resolution": "User's gmt-star interpreted as GPN-Star based on installed project source; no GMT-Star source exists in local predictor contract.",
        "identity_policy": "Reuse frozen vf_variant_ids only; validate each ID's existence and gene membership. Do not add keys, remap variants, substitute proxy alleles, or silently choose one of multiple nucleotide alleles. Base scalar columns require exactly one mapped allele; descriptive means/min/max retain all observed allele scores for sensitivity only. Missing scores remain NaN, never zero.",
        "unit": "one row per frozen protocol feature key, including rows later excluded from the plotted protocol estimates; join on gene,key to select final analysis set",
        "version_policy": "Multiple distinct version/value records are conflicts. GPN-Star expected pinned revision is recorded below. Empty version means not supplied, not latest or known identical model version.",
        "gpn_star_artifact_revision": GPN_VERSION,
        "gpn_star_direction": "Signed calibrated LLR: more negative indicates stronger functional constraint/effect. -LLR is an optional display-only sign flip, not a probability. abs_llr_calibrated is independently supplied and is not abs(signed LLR). Entropy near 1 is neutral, lower is more constrained.",
        "alphamissense_direction": "Higher raw AlphaMissense scores indicate predicted pathogenicity; this is not penetrance.",
        "alphagenome_metric": "Atlas variant-level AVI_SCORE raw_score mirrored as alphagenome_avi and verified against annotations_alphagenome_context. Higher indicates larger predicted variant impact. No tissue, gene-expression, splice, or other molecular track is selected or maximized. AVI includes AlphaMissense information and therefore is not independent of AlphaMissense.",
        "avi_context_conflict_variant_ids": sorted(set(avi_context_conflicts)),
        "sources": {
            "gpn_star": "https://huggingface.co/datasets/songlab/gpn-star-scores",
            "alphagenome": "https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/",
            "local_gpn_contract": "/Users/kronckbm/GitRepos/variantFeatures/docs/GPN_STAR.md",
            "local_avi_contract": "/Users/kronckbm/GitRepos/variantFeatures/variantfeatures/handlers/alphagenome.py",
        },
        "observed_source_versions": {
            metric: [
                {"version": v, "source": s, "fetched_at": t, "annotation_rows": n}
                for (v, s, t), n in sorted(provenance[metric].items())
            ]
            for metric in METRICS
        },
        "coverage_feature_universe": coverage,
    }
    (args.output_dir / "predictor_sources.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "keys": len(result),
                "queried_ids": len(ids),
                "avi_context_conflicts": len(set(avi_context_conflicts)),
                "output_dir": str(args.output_dir),
            }
        )
    )
    for gene in GENES:
        c = coverage[gene]
        print(
            gene,
            c["feature_keys"],
            {
                p: c["predictors"][p]["single_allele_scored"]
                for p in metadata["selected_primary_predictors"]
            },
        )


if __name__ == "__main__":
    main()
