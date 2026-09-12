#!/usr/bin/env python3
"""Freeze existing warehouse predictors for exact GRCh38 population alleles.

The source SQLite connection is read-only and held in one read transaction.
No annotation jobs, imports, or downloads are started. Missing or conflicting
scores remain missing; no latest-version rule is used.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
GENES = ("GCK", "HNF1A", "LDLR", "BRCA2", "KCNQ1")
METRICS = ("alphamissense", "gpn_star_m447_llr_calibrated", "alphagenome_avi")
DEFAULT_DB = Path("/Users/kronckbm/GitRepos/variantFeatures/data/variants.db")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_csv(frame, path):
    text = frame.to_csv(index=False, lineterminator="\n")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(text.encode())


def storage_guard(db):
    if not db.is_file() or db.stat().st_size == 0:
        raise RuntimeError(
            "Existing, nonempty source database is required; no database will be created"
        )
    if db == DEFAULT_DB:
        root = db.parent.parent
        for relative in ("data", "annovar/humandb"):
            link = root / relative
            if not link.is_symlink() or not link.is_dir():
                raise RuntimeError(
                    f"Required mounted source symlink is unavailable: {link}"
                )


def finite(value):
    return value is not None and np.isfinite(value)


def resolve(records, contexts=None):
    """One version/value only; AVI must agree with authoritative scalar context."""
    observed = [record for record in records if finite(record["score"])]
    distinct = {(record["predictor_version"], record["score"]) for record in observed}
    value, version, status = np.nan, "", "missing_annotation"
    if len(distinct) == 1:
        version, value = next(iter(distinct))
        status = "available"
    elif len(distinct) > 1:
        status = "version_or_value_conflict"
    if contexts is not None and observed:
        context_values = {
            (record["dataset_version"], record["raw_score"])
            for record in contexts
            if finite(record["raw_score"])
        }
        if len(context_values) != 1 or distinct != context_values:
            value, version, status = np.nan, "", "context_projection_conflict"
    return value, version, status


def self_check():
    a = {"predictor_version": "v1", "score": 0.2}
    assert resolve([a]) == (0.2, "v1", "available")
    assert resolve([a, dict(a)])[2] == "available"
    assert (
        resolve([a, {"predictor_version": "v2", "score": 0.2}])[2]
        == "version_or_value_conflict"
    )
    assert (
        resolve([a, {"predictor_version": "v1", "score": 0.3}])[2]
        == "version_or_value_conflict"
    )
    assert resolve([a], [{"dataset_version": "v1", "raw_score": 0.2}])[2] == "available"
    assert (
        resolve([a], [{"dataset_version": "v1", "raw_score": 0.3}])[2]
        == "context_projection_conflict"
    )
    assert resolve([a], [])[2] == "context_projection_conflict"
    assert resolve([])[2] == "missing_annotation"


def normalize_inventory(path):
    frame = pd.read_csv(
        path,
        dtype={
            "gene": str,
            "variant_id": str,
            "chrom": str,
            "chromosome": str,
            "ref": str,
            "alt": str,
        },
    )
    aliases = {"chromosome": "chrom", "position": "pos"}
    frame = frame.rename(
        columns={
            old: new
            for old, new in aliases.items()
            if old in frame and new not in frame
        }
    )
    required = ["gene", "variant_id", "chrom", "pos", "ref", "alt"]
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f"Population inventory missing columns: {sorted(missing)}")
    if frame[required].isna().any().any():
        raise ValueError("Population allele identity cannot be missing")
    if not set(frame.gene).issubset(GENES):
        raise ValueError("Unexpected gene outside the five-gene analysis scope")
    if frame.duplicated(["gene", "variant_id"]).any():
        raise ValueError("Population gene/variant ID must be unique")
    if (frame.pos <= 0).any() or (frame.pos % 1 != 0).any():
        raise ValueError("GRCh38 positions must be positive integers")
    frame["pos"] = frame.pos.astype(int)
    frame["chrom"] = frame.chrom.str.removeprefix("chr")
    if not (
        frame.ref.str.fullmatch("[ACGT]+") & frame.alt.str.fullmatch("[ACGT]+")
    ).all():
        raise ValueError(
            "Population ref/alt must be explicit uppercase nucleotide sequences"
        )
    if frame.duplicated(["gene", "chrom", "pos", "ref", "alt"]).any():
        raise ValueError("Distinct inventory IDs refer to the same gene/genomic allele")
    return frame.loc[:, required].copy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument(
        "--population",
        type=Path,
        default=HERE.parent / "population/population_variants.csv.gz",
    )
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    self_check()
    storage_guard(args.db)
    population = normalize_inventory(args.population)
    input_hash = sha256(args.population)
    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only=ON")
    conn.execute("BEGIN")
    # This first read establishes the logical snapshot under any concurrent WAL writer.
    schema_version = conn.execute("PRAGMA user_version").fetchone()[0]
    query_time = datetime.now(timezone.utc).isoformat()
    source_stat = args.db.stat()
    identities = {}
    scanned = 0
    intervals = []
    for (gene, chrom), group in population.groupby(["gene", "chrom"], sort=True):
        lower, upper = int(group.pos.min()), int(group.pos.max())
        intervals.append({"gene": gene, "chrom": chrom, "start": lower, "stop": upper})
        for record in conn.execute(
            "SELECT id,chromosome,position,ref,alt FROM variants "
            "WHERE chromosome=? AND position BETWEEN ? AND ?",
            (chrom, lower, upper),
        ):
            scanned += 1
            key = (
                record["chromosome"],
                record["position"],
                record["ref"],
                record["alt"],
            )
            if key in identities and identities[key] != record["id"]:
                raise ValueError(f"Warehouse has nonunique canonical allele: {key}")
            identities[key] = record["id"]
    allele_keys = list(
        zip(population.chrom, population.pos, population.ref, population.alt)
    )
    population["vf_variant_id"] = pd.array(
        [identities.get(key) for key in allele_keys], dtype="Int64"
    )
    ids = sorted(population.vf_variant_id.dropna().astype(int).unique().tolist())
    scores = defaultdict(lambda: defaultdict(list))
    gene_membership = defaultdict(set)
    contexts = defaultdict(list)
    annotation_rows, context_rows = [], []
    provenance = defaultdict(Counter)
    for start in range(0, len(ids), 400):
        chunk = ids[start : start + 400]
        marks = ",".join("?" for _ in chunk)
        metrics = ",".join("?" for _ in METRICS)
        for record in conn.execute(
            "SELECT variant_id,predictor,predictor_version,score,source,fetched_at "
            f"FROM annotations_pathogenicity WHERE variant_id IN ({marks}) AND predictor IN ({metrics})",
            [*chunk, *METRICS],
        ):
            record = dict(record)
            scores[record["variant_id"]][record["predictor"]].append(record)
            annotation_rows.append(record)
            provenance[record["predictor"]][
                (record["predictor_version"], record["source"], record["fetched_at"])
            ] += 1
        for record in conn.execute(
            "SELECT variant_id,dataset_version,context_id,raw_score,source,fetched_at "
            f"FROM annotations_alphagenome_context WHERE variant_id IN ({marks}) AND scorer='AVI_SCORE'",
            chunk,
        ):
            record = dict(record)
            contexts[record["variant_id"]].append(record)
            context_rows.append(record)
        for record in conn.execute(
            f"SELECT DISTINCT variant_id,gene_symbol FROM variant_consequences WHERE variant_id IN ({marks})",
            chunk,
        ):
            gene_membership[record["variant_id"]].add(record["gene_symbol"])
        if start % 4000 == 0:
            print(
                f"Read predictor batches: {min(start + 400, len(ids))}/{len(ids)} IDs",
                flush=True,
            )
    conn.rollback()
    conn.close()
    result_rows = []
    for source in population.to_dict("records"):
        vid = source["vf_variant_id"]
        if pd.isna(vid):
            state = "no_exact_warehouse_allele"
            vid = None
        elif source["gene"] not in gene_membership[vid]:
            state = "warehouse_gene_membership_missing"
        else:
            state = "exact_allele_and_gene"
        row = {
            **source,
            "warehouse_status": state,
            "warehouse_genes": ";".join(sorted(g for g in gene_membership[vid] if g))
            if vid is not None
            else "",
        }
        for metric in METRICS:
            value, version, status = resolve(
                scores[vid][metric],
                contexts[vid] if metric == "alphagenome_avi" else None,
            )
            if state != "exact_allele_and_gene":
                value, version, status = np.nan, "", state
            row[metric] = value
            row[f"{metric}_version"] = version
            row[f"{metric}_status"] = status
        row["avi_scalar_context_rows"] = len(contexts[vid]) if vid is not None else 0
        result_rows.append(row)
    result = pd.DataFrame(result_rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    save_csv(result, args.output_dir / "population_predictors.csv.gz")
    annotations = pd.DataFrame(
        annotation_rows,
        columns=[
            "variant_id",
            "predictor",
            "predictor_version",
            "score",
            "source",
            "fetched_at",
        ],
    )
    annotations = annotations.sort_values(
        ["variant_id", "predictor", "predictor_version"]
    )
    save_csv(annotations, args.output_dir / "predictor_annotation_rows.csv.gz")
    context_table = pd.DataFrame(
        context_rows,
        columns=[
            "variant_id",
            "dataset_version",
            "context_id",
            "raw_score",
            "source",
            "fetched_at",
        ],
    )
    context_table = context_table.sort_values(
        ["variant_id", "dataset_version", "context_id"]
    )
    save_csv(context_table, args.output_dir / "avi_scalar_contexts.csv.gz")
    coverage = []
    for gene, group in result.groupby("gene", sort=True):
        for metric in METRICS:
            coverage.append(
                {
                    "gene": gene,
                    "predictor": metric,
                    "inventory_rows": len(group),
                    "exact_warehouse_alleles": int(group.vf_variant_id.notna().sum()),
                    "exact_allele_and_gene": int(
                        group.warehouse_status.eq("exact_allele_and_gene").sum()
                    ),
                    "available": int(group[metric].notna().sum()),
                    "version_or_value_conflict": int(
                        group[f"{metric}_status"].eq("version_or_value_conflict").sum()
                    ),
                    "context_projection_conflict": int(
                        group[f"{metric}_status"]
                        .eq("context_projection_conflict")
                        .sum()
                    ),
                    "missing_annotation": int(
                        group[f"{metric}_status"].eq("missing_annotation").sum()
                    ),
                }
            )
    coverage = pd.DataFrame(coverage)
    coverage.to_csv(args.output_dir / "coverage.csv", index=False)
    output_names = [
        "population_predictors.csv.gz",
        "predictor_annotation_rows.csv.gz",
        "avi_scalar_contexts.csv.gz",
        "coverage.csv",
    ]
    manifest = {
        "queried_at_utc": query_time,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "read_only": True,
        "sqlite_open_mode": "ro",
        "query_only": True,
        "one_read_transaction": True,
        "db_path": str(args.db),
        "db_schema_version": schema_version,
        "db_file_size_at_snapshot": source_stat.st_size,
        "db_file_mtime_ns_at_snapshot": source_stat.st_mtime_ns,
        "db_hash_policy": "No full database hash; a single logical read snapshot and exact queried annotation/context rows are preserved.",
        "population_input": {
            "path": str(args.population),
            "sha256": input_hash,
            "rows": len(population),
        },
        "collector_sha256": sha256(Path(__file__)),
        "self_checks_passed": True,
        "unique_exact_warehouse_variant_ids": len(ids),
        "warehouse_variant_rows_scanned": scanned,
        "coordinate_query_intervals": intervals,
        "predictors": list(METRICS),
        "predictor_table": "annotations_pathogenicity",
        "avi_validation_table": "annotations_alphagenome_context",
        "gpn_table_note": "Selected signed M447 LLR is stored in annotations_pathogenicity; conservation stores GPN entropy and is not queried here.",
        "identity_policy": "Exact GRCh38 chromosome/position/ref/alt and matching gene membership; only a chr prefix is removed. No allele normalization, protein-level collapse, proxy mapping, or transcript score substitution. Every inventory row is retained.",
        "score_policy": "One distinct finite (version,score) pair required. Multiple versions conflict even if numeric scores agree. No arbitrary latest version. Missing remains NaN. AVI must exactly agree with one distinct (dataset_version,raw_score) in the authoritative AVI_SCORE contexts.",
        "alphamissense_direction": "Higher means predicted pathogenicity, not penetrance.",
        "gpn_direction": "Signed calibrated M447 LLR: more negative implies greater predicted impact/constraint; not a probability.",
        "alphagenome_metric": "Scalar Atlas AVI_SCORE raw score, higher means greater predicted variant impact. AVI incorporates AlphaMissense; it is not an independent predictor of AlphaMissense.",
        "observed_source_versions": {
            metric: [
                {
                    "version": version,
                    "source": source,
                    "fetched_at": fetched,
                    "annotation_rows": count,
                }
                for (version, source, fetched), count in sorted(
                    provenance[metric].items(),
                    key=lambda item: tuple(str(x) for x in item[0]),
                )
            ]
            for metric in METRICS
        },
        "coverage": coverage.to_dict("records"),
        "output_sha256": {
            name: sha256(args.output_dir / name) for name in output_names
        },
    }
    (args.output_dir / "predictor_provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    if sha256(args.population) != input_hash:
        raise RuntimeError(
            "Population input changed during query; rerun against a frozen inventory"
        )
    print(coverage.to_string(index=False), flush=True)
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir),
                "rows": len(result),
                "warehouse_ids": len(ids),
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
