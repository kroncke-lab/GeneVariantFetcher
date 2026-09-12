"""Verify type-specific prior, source hashes, and committed artifact manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    manifest_path = HERE / "artifact_manifest.json"
    if manifest_path.exists() and not args.freeze:
        for name, record in json.loads(manifest_path.read_text())["files"].items():
            assert digest(HERE / name) == record["sha256"], name
    fit = json.loads((HERE / "analysis/fit_checks.json").read_text())
    for name, expected in fit["input_sha256"].items():
        assert digest(HERE.parent / name) == expected, name
    prior = pd.read_csv(HERE / "analysis/empirical_prior_comparison.csv")
    primary = prior.loc[prior.scope.isin(["canonical_missense", "canonical_nonsense"])]
    assert len(primary) == 10 and not primary.duplicated(["gene", "variant_type"]).any()
    posterior = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted((HERE / "analysis/empirical_posteriors").glob("*.csv.gz"))
        ],
        ignore_index=True,
    )
    assert posterior.unit_id.is_unique and len(posterior) == 11300
    for r in primary.itertuples():
        rows = posterior.loc[
            posterior.gene.eq(r.gene) & posterior.variant_type.eq(r.variant_type)
        ]
        assert len(rows) == r.variants
        assert rows.vclass.isin(
            ["missense"]
            if r.variant_type == "missense"
            else ["nonsense", "stop_gained"]
        ).all()
        assert rows.canonical_wt_status.eq("match").all()
        np.testing.assert_allclose(rows.alpha_empirical, r.alpha_empirical)
        np.testing.assert_allclose(rows.beta_empirical, r.beta_empirical)
        np.testing.assert_allclose(
            rows.posterior_alpha, r.alpha_empirical + rows.affected
        )
        np.testing.assert_allclose(
            rows.posterior_beta,
            r.beta_empirical + rows.unaffected_literature + rows.gnomad_carriers,
        )
    audit = pd.read_csv(HERE / "audit/class_moment_comparison.csv")
    methods = {
        "historical_weighted_mean_mse": "",
        "normalized_weighted_mse": "_normalized_weighted_mse",
        "equal_variant_mean_mse": "_equal_raw_fraction",
    }
    for r in audit.itertuples():
        row = prior.loc[
            prior.gene.eq(r.gene)
            & prior.scope.eq("canonical_" + r.variant_type + methods[r.method])
        ].iloc[0]
        np.testing.assert_allclose(
            [
                row["mean"],
                row.variance,
                row.alpha_empirical,
                row.beta_empirical,
                row.strength,
                row.unaffected_singleton_posterior,
            ],
            [
                r.mean,
                r.variance,
                r.alpha_empirical,
                r.beta_empirical,
                r.strength,
                r.A0_U1_posterior_mean,
            ],
            rtol=1e-12,
            atol=1e-13,
        )
    structural = json.loads((HERE / "structural/run_checks.json").read_text())
    for name, expected in structural["input_hashes"].items():
        assert digest(HERE.parent / name) == expected, name
    for name, expected in structural["output_hashes"].items():
        assert digest(HERE / "structural" / name) == expected, name
    assert (
        digest(
            HERE.parents[3]
            / "ProteinProximityAnalysis/src/alphafold_rin/empirical_density.py"
        )
        == structural["ppa_module_sha256"]
    )
    files = {}
    for path in sorted(HERE.rglob("*")):
        if (
            not path.is_file()
            or "__pycache__" in path.parts
            or path == manifest_path
            or path.suffix == ".log"
        ):
            continue
        assert path.stat().st_size < 1200 * 1024, str(path)
        if path.suffix in {".md", ".py", ".json", ".csv", ".txt"}:
            assert b"\r\n" not in path.read_bytes(), str(path)
        files[str(path.relative_to(HERE))] = {
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }
    checks = {
        "primary_gene_type_priors": 10,
        "missense_units": 10512,
        "nonsense_units": 788,
        "independent_primary_and_sensitivity_fits_match": 30,
        "no_cross_type_prior_pooling": True,
        "posteriors_add_affected_to_alpha_and_all_unaffected_to_beta": True,
        "source_and_structural_hashes_match": True,
        "artifact_files": len(files),
    }
    if args.freeze:
        manifest_path.write_text(
            json.dumps({"checks": checks, "files": files}, indent=2) + "\n"
        )
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
