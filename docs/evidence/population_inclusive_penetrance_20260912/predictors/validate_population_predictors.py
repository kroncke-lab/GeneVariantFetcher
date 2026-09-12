"""Validate population/predictor identity and resolve scores from saved rows."""

from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from collect_population_predictors import METRICS, normalize_inventory, resolve, sha256


HERE = Path(__file__).resolve().parent


def main():
    provenance = json.loads((HERE / "predictor_provenance.json").read_text())
    population_path = Path(provenance["population_input"]["path"])
    assert sha256(population_path) == provenance["population_input"]["sha256"]
    assert (
        sha256(HERE / "collect_population_predictors.py")
        == provenance["collector_sha256"]
    )
    for name, expected in provenance["output_sha256"].items():
        assert sha256(HERE / name) == expected, name
    population = normalize_inventory(population_path)
    result = pd.read_csv(
        HERE / "population_predictors.csv.gz",
        dtype={"chrom": str, **{f"{metric}_version": str for metric in METRICS}},
    )
    assert len(result) == len(population)
    assert not result.duplicated(["gene", "variant_id"]).any()
    identity = ["gene", "variant_id", "chrom", "pos", "ref", "alt"]
    pd.testing.assert_frame_equal(
        result[identity], population[identity], check_dtype=False
    )
    annotations = pd.read_csv(
        HERE / "predictor_annotation_rows.csv.gz", keep_default_na=False
    )
    contexts = pd.read_csv(HERE / "avi_scalar_contexts.csv.gz", keep_default_na=False)
    annotations["score"] = pd.to_numeric(annotations.score, errors="coerce")
    contexts["raw_score"] = pd.to_numeric(contexts.raw_score, errors="coerce")
    scores = defaultdict(lambda: defaultdict(list))
    avi = defaultdict(list)
    for row in annotations.to_dict("records"):
        scores[row["variant_id"]][row["predictor"]].append(row)
    for row in contexts.to_dict("records"):
        avi[row["variant_id"]].append(row)
    max_difference = {metric: 0.0 for metric in METRICS}
    for row in result.to_dict("records"):
        valid_identity = row["warehouse_status"] == "exact_allele_and_gene"
        if valid_identity:
            assert row["gene"] in row["warehouse_genes"].split(";")
            assert np.isfinite(row["vf_variant_id"])
        for metric in METRICS:
            if valid_identity:
                value, version, status = resolve(
                    scores[row["vf_variant_id"]][metric],
                    avi[row["vf_variant_id"]] if metric == "alphagenome_avi" else None,
                )
                assert row[f"{metric}_status"] == status
                if status == "available":
                    difference = abs(row[metric] - value)
                    max_difference[metric] = max(max_difference[metric], difference)
                    assert np.isclose(row[metric], value, rtol=1e-13, atol=1e-15)
                    saved_version = row[f"{metric}_version"]
                    saved_version = "" if pd.isna(saved_version) else saved_version
                    assert saved_version == version
                else:
                    assert pd.isna(row[metric])
            else:
                assert pd.isna(row[metric])
                assert row[f"{metric}_status"] == row["warehouse_status"]
    am = result.alphamissense.dropna()
    assert am.between(0, 1).all()
    report = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_and_output_hashes_match": True,
        "input_rows_preserved": len(result),
        "exact_identity_columns_unchanged": True,
        "score_resolution_reproduced_from_saved_rows": True,
        "max_absolute_score_roundtrip_difference": max_difference,
        "alphamissense_values_in_unit_interval": True,
        "gene_membership_validated_for_every_available_score": True,
        "missing_or_conflicting_scores_not_imputed": True,
        "validator_sha256": sha256(Path(__file__)),
    }
    (HERE / "predictor_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
