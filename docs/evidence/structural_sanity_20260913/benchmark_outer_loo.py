"""Benchmark one BRCA2 target batch against PPA's single-exclusion method."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import empirical_variant_density
from outer_loo import all_excluded_density


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=int, default=128)
    parser.add_argument("--donor-batch-size", type=int, default=512)
    args = parser.parse_args()
    files = sorted(
        (
            HERE.parent
            / "class_matched_penetrance_20260912/analysis/empirical_posteriors"
        ).glob("BRCA2.part*.csv.gz")
    )
    data = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    data = data.loc[data.variant_type.eq("missense")].copy()
    data["variant_id"] = data.unit_id
    data["canonical_variant_id"] = data.unit_id
    data["canonical_pos"] = data.aa_pos.astype(int)
    data["donor_eligible"] = True
    geometry_file = HERE / "geometry/BRCA2/primary_geometry.csv.gz"
    geometry = pd.read_csv(geometry_file).fillna({"idr_segment": ""})
    targets = data.variant_id.iloc[: args.targets].tolist()
    started = perf_counter()
    model = empirical_variant_density(
        data,
        geometry,
        target_ids=targets,
        include_context_model=True,
        include_context_weights=False,
    ).context_model
    build_seconds = perf_counter() - started
    started = perf_counter()
    batched = all_excluded_density(model, donor_batch_size=args.donor_batch_size)
    batched_seconds = perf_counter() - started
    started = perf_counter()
    expected = np.column_stack(
        [model.density([donor]).to_numpy() for donor in model.donor_ids]
    )
    repeated_seconds = perf_counter() - started
    np.testing.assert_allclose(batched, expected, atol=1e-12, equal_nan=True)
    record = {
        "targets": len(targets),
        "donors": len(model.donor_ids),
        "contexts": len(model.context_densities),
        "donor_batch_size": args.donor_batch_size,
        "model_build_seconds": build_seconds,
        "batched_exclusion_seconds": batched_seconds,
        "repeated_exclusion_seconds": repeated_seconds,
        "speedup": repeated_seconds / batched_seconds,
        "maximum_absolute_difference": float(
            np.nanmax(abs(batched.to_numpy() - expected))
        ),
        "output_array_bytes": batched.to_numpy().nbytes,
        "helper_sha256": hashlib.sha256(
            (HERE / "outer_loo.py").read_bytes()
        ).hexdigest(),
        "geometry_sha256": hashlib.sha256(geometry_file.read_bytes()).hexdigest(),
    }
    (HERE / "outer_loo_benchmark.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
