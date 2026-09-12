"""Verify committed artifacts and the population/count/posterior invariants."""

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
    parser.add_argument(
        "--freeze",
        action="store_true",
        help="Write artifact manifest after verification",
    )
    args = parser.parse_args()
    manifest_path = HERE / "artifact_manifest.json"
    if manifest_path.exists() and not args.freeze:
        manifest = json.loads(manifest_path.read_text())
        for name, expected in manifest["files"].items():
            assert digest(HERE / name) == expected["sha256"], name
    full = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted((HERE / "analysis/empirical_posteriors").glob("*.csv.gz"))
        ],
        ignore_index=True,
    )
    members = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted((HERE / "analysis/population_membership").glob("*.csv.gz"))
        ],
        ignore_index=True,
    )
    coverage = pd.read_csv(HERE / "analysis/union_coverage.csv")
    assert full.unit_id.is_unique and len(full) == coverage.union_units.sum()
    assert not members.duplicated(["gene", "variant_id"]).any()
    assert len(members) == coverage.eligible_population_alleles.sum()
    assert set(members.unit_id) <= set(full.unit_id)
    np.testing.assert_allclose(
        full.posterior_alpha, full.alpha_empirical + full.affected
    )
    np.testing.assert_allclose(
        full.posterior_beta,
        full.beta_empirical + full.unaffected_literature + full.gnomad_carriers,
    )
    np.testing.assert_allclose(
        full.posterior_mean,
        full.posterior_alpha / (full.posterior_alpha + full.posterior_beta),
    )
    assert (full.n > 0).all()
    population_only = full.loc[full.origin.eq("population_only")]
    assert (
        population_only.affected.eq(0).all()
        and population_only.gnomad_carriers.gt(0).all()
    )
    assert (
        population_only.posterior_mean
        < population_only.alpha_empirical
        / (population_only.alpha_empirical + population_only.beta_empirical)
    ).all()
    for r in coverage.itertuples():
        gene = full.loc[full.gene.eq(r.gene)]
        assert gene.affected.sum() == r.affected
        assert gene.unaffected_literature.sum() == r.unaffected_literature
        assert gene.gnomad_carriers.sum() == r.gnomad_unaffected_carrier_observations
    union_checks = json.loads((HERE / "analysis/union_checks.json").read_text())
    for name, expected in union_checks["input_sha256"].items():
        assert digest(HERE / name) == expected, name
    structural = json.loads((HERE / "structural/run_checks.json").read_text())
    for name, expected in structural["input_hashes"].items():
        assert digest(HERE.parent / name) == expected, name
    for name, expected in structural["output_hashes"].items():
        assert digest(HERE / "structural" / name) == expected, name
    module = (
        HERE.parents[3]
        / "ProteinProximityAnalysis/src/alphafold_rin/empirical_density.py"
    )
    assert digest(module) == structural["ppa_module_sha256"]
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
        if path.suffix in {".py", ".md", ".json", ".csv", ".txt"}:
            assert b"\r\n" not in path.read_bytes(), str(path)
        files[str(path.relative_to(HERE))] = {
            "sha256": digest(path),
            "bytes": path.stat().st_size,
        }
    checks = {
        "union_units": len(full),
        "population_alleles": len(members),
        "population_only_units": len(population_only),
        "unique_population_membership": True,
        "affected_and_unaffected_count_conservation": True,
        "posterior_alpha_adds_affected_beta_adds_all_unaffected": True,
        "all_population_only_rows_move_below_shared_prior": True,
        "input_and_structural_hashes_match": True,
        "source_module_matches_committed_PPA": True,
        "file_size_and_LF_checks": True,
        "artifact_files": len(files),
    }
    if args.freeze:
        manifest_path.write_text(
            json.dumps({"checks": checks, "files": files}, indent=2) + "\n"
        )
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
