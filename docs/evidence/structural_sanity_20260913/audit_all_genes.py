"""Reconstruct fixed neighborhood densities and their prior/count components.

No model, count, prior or donor policy is changed. Every observed missense unit
is retained, including geometry-unsupported rows with missing density. Saved
weights are interpreted as variant weights, not independent carrier counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import load_npz


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
REPO = HERE.parents[2]
CLASS = EVIDENCE / "class_matched_penetrance_20260912"
OLD = EVIDENCE / "missense_structural_extension_20260912"
OUT = HERE / "analysis"
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
LIMIT = 1_200_000
SOURCES = {}
STATS = {}


def record(path, expected=None):
    path = Path(path).resolve()
    if path not in SOURCES:
        before = path.stat()
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Source changed while hashing: {path}")
        SOURCES[path] = sha
        STATS[path] = (after.st_size, after.st_mtime_ns)
    if expected is not None and SOURCES[path] != expected:
        raise ValueError(f"Source hash does not match recorded input: {path}")
    return path


def read_csv(path):
    return pd.read_csv(record(path))


def concat(paths):
    if not paths:
        raise ValueError("Required source shards are missing")
    return pd.concat([read_csv(p) for p in paths], ignore_index=True)


def write_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )
    if path.stat().st_size >= LIMIT:
        raise ValueError(f"Artifact too large: {path}")


def write_shards(frame, gene):
    folder = OUT / "variant_diagnostics"
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.glob(f"{gene}.part*.csv.gz"):
        old.unlink()
    for number, start in enumerate(range(0, len(frame), 1500), 1):
        write_csv(
            frame.iloc[start : start + 1500], folder / f"{gene}.part{number:03d}.csv.gz"
        )


def inputs(gene, require_corrected):
    units = concat(
        sorted((CLASS / "analysis/empirical_posteriors").glob(f"{gene}.part*.csv.gz"))
    )
    units = units.loc[units.variant_type.eq("missense")].copy()
    if units.unit_id.duplicated().any():
        raise ValueError(f"Duplicate canonical count unit: {gene}")
    if gene == "GCK":
        density = read_csv(CLASS / "structural/GCK_primary_variant_density.csv")
        mode = "frozen_class_specific_monomer"
    elif gene == "BRCA2" and (OUT / "BRCA2/cache_manifest.json").exists():
        density = concat(sorted((OUT / "BRCA2").glob("primary_density.part*.csv.gz")))
        mode = "corrected_local_3d_and_canonical_polymer"
    else:
        if gene == "BRCA2" and require_corrected:
            raise ValueError("Corrected BRCA2 primary and cache receipt are required")
        density = concat(
            sorted((OLD / f"analysis/{gene}").glob("primary_density.part*.csv.gz"))
        )
        mode = (
            "old_domain_only_baseline"
            if gene == "BRCA2"
            else "frozen_missense_extension"
        )
    if density.variant_id.duplicated().any() or set(density.variant_id) != set(
        units.unit_id
    ):
        raise ValueError(f"Incomplete or mismatched primary identity universe: {gene}")
    units = units.set_index("unit_id").loc[density.variant_id].reset_index()
    np.testing.assert_allclose(units.posterior_mean, density.posterior_mean, atol=1e-12)
    np.testing.assert_allclose(units.n, units.affected + units.unaffected)
    np.testing.assert_allclose(
        units.unaffected, units.unaffected_literature + units.gnomad_carriers
    )
    if (
        not (units.n > 0).all()
        or units.alpha_empirical.nunique() != 1
        or units.beta_empirical.nunique() != 1
    ):
        raise ValueError(f"Invalid observed counts or nonshared class prior: {gene}")
    return units, density, mode


def weighted_values(units, density, gene, mode):
    ids = pd.Index(units.unit_id)
    size = len(ids)
    strength = float(units.alpha_empirical.iloc[0] + units.beta_empirical.iloc[0])
    mean = float(units.alpha_empirical.iloc[0] / strength)
    n = units.n.to_numpy(float)
    affected = units.affected.to_numpy(float)
    retention = strength / (strength + n)
    count_component = affected / (strength + n)
    values = np.column_stack(
        [
            np.ones(size),
            retention,
            count_component,
            affected / n,
            units.posterior_mean,
            units.origin.eq("population_only"),
            n == 1,
        ]
    )
    result = np.zeros((size, values.shape[1]))
    pairs = 0
    if gene == "BRCA2" and mode == "corrected_local_3d_and_canonical_polymer":
        path = record(OUT / "BRCA2/cache_manifest.json")
        cache = json.loads(path.read_text())
        expected = hashlib.sha256("\n".join(ids).encode()).hexdigest()
        if cache["variant_ids_sha256"] != expected:
            raise ValueError(
                "Corrected BRCA2 weight order differs from target identities"
            )
        for relative, sha in cache["source_hashes"].items():
            record(EVIDENCE / relative, sha)
        covered = np.zeros(size, dtype=bool)
        for filename, sha in sorted(cache["raw_weight_shards"].items()):
            path = record(
                REPO / "results/structural_sanity_20260913/BRCA2" / filename, sha
            )
            _, start, stop = path.stem.split("_")
            start, stop = int(start), int(stop)
            weights = load_npz(path)
            if weights.shape != (stop - start, size) or covered[start:stop].any():
                raise ValueError("Incoherent or repeated corrected BRCA2 weight batch")
            if np.any(weights.data < 0) or not np.isfinite(weights.data).all():
                raise ValueError("Invalid corrected BRCA2 weights")
            if weights[np.arange(stop - start), np.arange(start, stop)].sum() != 0:
                raise ValueError("A BRCA2 target donated to itself")
            result[start:stop] = weights @ values
            pairs += weights.nnz
            covered[start:stop] = True
        if not covered.all() or pairs != cache["positive_primary_donor_pairs"]:
            raise ValueError("Corrected BRCA2 weights are incomplete")
    elif gene == "GCK":
        covered = np.zeros(size, dtype=bool)
        for path in sorted(
            (CLASS / "structural").glob("GCK_primary_normalized_weights.part*.csv.gz")
        ):
            table = read_csv(path).set_index("variant_id")
            if set(table.columns) != set(ids):
                raise ValueError("GCK donor identities differ")
            rows = ids.get_indexer(table.index)
            if (rows < 0).any() or covered[rows].any():
                raise ValueError("Invalid or repeated GCK target rows")
            weights = table.reindex(columns=ids).to_numpy(float)
            if (
                np.any(weights < 0)
                or not np.isfinite(weights).all()
                or np.any(weights[np.arange(len(rows)), rows] != 0)
            ):
                raise ValueError("Invalid GCK weights or target self-donation")
            result[rows] = weights @ values
            pairs += int((weights > 0).sum())
            covered[rows] = True
        if not covered.all():
            raise ValueError("GCK weights do not include every target")
    else:
        seen = np.zeros(size * size, dtype=bool)
        for path in sorted(
            (OLD / f"analysis/{gene}").glob("primary_weights.part*.csv.gz")
        ):
            table = read_csv(path)
            rows = ids.get_indexer(table.target_id)
            columns = ids.get_indexer(table.donor_id)
            weights = table.normalized_weight.to_numpy(float)
            if (rows < 0).any() or (columns < 0).any() or (rows == columns).any():
                raise ValueError(f"Unknown or self donor pair: {gene}")
            flat = rows * size + columns
            if seen[flat].any() or len(np.unique(flat)) != len(flat):
                raise ValueError(f"Duplicate donor pair: {gene}")
            seen[flat] = True
            if not np.isfinite(weights).all() or np.any(weights <= 0):
                raise ValueError(f"Invalid donor weights: {gene}")
            for column in range(values.shape[1]):
                result[:, column] += np.bincount(
                    rows, weights=weights * values[columns, column], minlength=size
                )
            pairs += len(table)
    supported = density.density.notna().to_numpy()
    np.testing.assert_allclose(result[supported, 0], 1, atol=1e-10)
    np.testing.assert_allclose(result[~supported, 0], 0, atol=1e-15)
    direct = mean * result[:, 1]
    reconstructed = direct + result[:, 2]
    np.testing.assert_allclose(
        reconstructed[supported], density.loc[supported, "density"], atol=1e-11
    )
    np.testing.assert_allclose(
        result[supported, 4], reconstructed[supported], atol=1e-11
    )
    if mode == "corrected_local_3d_and_canonical_polymer":
        for column, actual in [
            ("prior_component", direct),
            ("counts_component", result[:, 2]),
            ("neighborhood_prior_retention", result[:, 1]),
            ("raw_count_neighborhood", result[:, 3]),
        ]:
            np.testing.assert_allclose(
                density.loc[supported, column], actual[supported], atol=1e-11
            )
    result[~supported] = np.nan
    direct[~supported] = np.nan
    return result, direct, pairs


def audit_gene(gene, require_corrected):
    units, density, mode = inputs(gene, require_corrected)
    values, direct, pairs = weighted_values(units, density, gene, mode)
    strength = float(units.alpha_empirical.iloc[0] + units.beta_empirical.iloc[0])
    mean = float(units.alpha_empirical.iloc[0] / strength)
    keep = [
        "gene",
        "unit_id",
        "literature_key",
        "protein_key",
        "origin",
        "aa_pos",
        "affected",
        "unaffected_literature",
        "gnomad_carriers",
        "unaffected",
        "n",
        "posterior_mean",
        "posterior_lower_95",
        "posterior_upper_95",
    ]
    frame = units[keep].copy()
    frame["primary_source"] = mode
    frame["prior_mean"] = mean
    frame["prior_strength"] = strength
    frame["density"] = density.density.to_numpy()
    frame["density_source"] = density.density_source.fillna("unavailable").to_numpy()
    frame["neighborhood_prior_retention"] = values[:, 1]
    frame["prior_component"] = direct
    frame["counts_component"] = values[:, 2]
    frame["raw_count_neighborhood"] = values[:, 3]
    frame["prior_share_of_density"] = frame.prior_component / frame.density
    frame["density_minus_prior"] = frame.density - mean
    frame["density_minus_own_posterior"] = frame.density - frame.posterior_mean
    frame["population_only_donor_weight"] = values[:, 5]
    frame["singleton_donor_weight"] = values[:, 6]
    for column in [
        "kish_donor_n",
        "largest_donor_share",
        "same_residue_donors",
        "donor_count",
    ]:
        if column in density:
            frame[column] = density[column].to_numpy()
    frame["zero_affected_control"] = (
        frame.affected.eq(0)
        & frame.unaffected.ge(10)
        & frame.posterior_mean.le(0.01)
        & frame.density.ge(0.10)
    )
    frame["population_only_control"] = frame.zero_affected_control & frame.origin.eq(
        "population_only"
    )
    write_shards(frame, gene)
    supported = frame.loc[frame.density.notna()]
    counts = frame.density_source.value_counts()
    summary = {
        "gene": gene,
        "primary_source": mode,
        "variants": len(frame),
        "supported": len(supported),
        "structured": int(counts.get("structured", 0)),
        "polymer": int(counts.get("polymer", 0)),
        "mixed": int(counts.get("mixed", 0)),
        "unavailable": int(frame.density.isna().sum()),
        "prior_mean": mean,
        "prior_strength": strength,
        "median_own_posterior_supported": supported.posterior_mean.median(),
        "median_density": supported.density.median(),
        "median_raw_count_neighborhood": supported.raw_count_neighborhood.median(),
        "median_prior_retention": supported.neighborhood_prior_retention.median(),
        "median_prior_component": supported.prior_component.median(),
        "median_counts_component": supported.counts_component.median(),
        "median_prior_share_of_density": supported.prior_share_of_density.median(),
        "aggregate_prior_component_share": supported.prior_component.sum()
        / supported.density.sum(),
        "median_population_only_donor_weight": supported.population_only_donor_weight.median(),
        "median_singleton_donor_weight": supported.singleton_donor_weight.median(),
        "zero_affected_controls": int(frame.zero_affected_control.sum()),
        "population_only_controls": int(frame.population_only_control.sum()),
        "population_only_targets": int(frame.origin.eq("population_only").sum()),
        "positive_weight_pairs_checked": pairs,
        "max_decomposition_error": float(
            abs(
                supported.density
                - supported.prior_component
                - supported.counts_component
            ).max()
        ),
        "validation_scope": "Fixed-prior variant-only internal posterior-label evaluation; no independent disease-risk calibration",
    }
    examples = (
        frame.loc[frame.zero_affected_control]
        .sort_values(["density_minus_own_posterior", "unaffected"], ascending=False)
        .head(5)
        .copy()
    )
    examples.insert(
        1,
        "control_category",
        np.where(
            examples.population_only_control, "population_only", "other_zero_affected"
        ),
    )
    coverage = pd.DataFrame(
        [
            {
                "gene": gene,
                "primary_source": mode,
                "density_source": source,
                "variants": int(count),
            }
            for source, count in counts.items()
        ]
    )
    print(
        json.dumps(
            {
                k: summary[k]
                for k in [
                    "gene",
                    "primary_source",
                    "supported",
                    "median_density",
                    "median_prior_share_of_density",
                    "zero_affected_controls",
                ]
            }
        ),
        flush=True,
    )
    return summary, examples, coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-corrected-brca2", action="store_true")
    args = parser.parse_args()
    record(Path(__file__))
    summaries, controls, coverage = [], [], []
    for gene in GENES:
        summary, examples, sources = audit_gene(gene, args.require_corrected_brca2)
        summaries.append(summary)
        controls.append(examples)
        coverage.append(sources)
    write_csv(pd.DataFrame(summaries), OUT / "allgene_summary.csv")
    write_csv(
        pd.concat(controls, ignore_index=True), OUT / "strongest_control_examples.csv"
    )
    write_csv(
        pd.concat(coverage, ignore_index=True),
        OUT / "allgene_density_source_coverage.csv",
    )
    for path, signature in STATS.items():
        stat = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != signature:
            raise ValueError(f"Source changed during the audit: {path}")
    receipt = {
        "sources": {str(path): sha for path, sha in SOURCES.items()},
        "all_sources_unchanged_during_audit": True,
        "definitions": {
            "neighborhood_prior_retention": "sum_j W_ij*S/(S+n_j)",
            "prior_component": "mu*sum_j W_ij*S/(S+n_j)",
            "counts_component": "sum_j W_ij*A_j/(S+n_j)",
            "raw_count_neighborhood": "sum_j W_ij*A_j/n_j",
            "control": "A=0,U>=10,own posterior<=0.01,density>=0.10; diagnostic not a biological benign classification",
            "prior_share_of_density": "direct prior component / density; not total causal influence of prior choice",
        },
        "no_models_refit": True,
        "fixed_gene_by_missense_priors": True,
        "limitations": [
            "Neighborhood density is a posterior-derived feature, not the target variant's measured penetrance.",
            "Carriers from gnomAD retain the user-fixed unaffected assumption.",
            "Conditional intervals and internal LOO do not establish independent disease-outcome calibration.",
            "Medians of components need not sum to median density.",
        ],
    }
    (OUT / "allgene_audit_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
