"""Read exact joint XX/XY counts for public BRCA2 reference variants."""

from datetime import UTC, datetime
import gzip
import hashlib
import json
from pathlib import Path
import time

import requests


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = REPO / "results/brca2_sex_endpoint_audit_20260914"
API = "https://gnomad.broadinstitute.org/api"
QUERY = """
query SexAudit($id: String!) {
  variant(variantId: $id, dataset: gnomad_r4) {
    variant_id
    exome { ac an homozygote_count populations { id ac an homozygote_count } }
    genome { ac an homozygote_count populations { id ac an homozygote_count } }
    joint { ac an homozygote_count populations { id ac an homozygote_count } }
  }
}
"""


def fetch_one(allele):
    payload = {"query": QUERY, "variables": {"id": allele}}
    path = RAW / (
        "sex_probe.json.gz" if allele == "13-32363389-G-T" else f"{allele}.json.gz"
    )
    if path.exists():
        result = json.loads(gzip.decompress(path.read_bytes()))
        assert result["request"] == payload
    else:
        response = requests.post(API, json=payload, timeout=45)
        response.raise_for_status()
        result = {
            "api": API,
            "fetched_at": datetime.now(UTC).isoformat(),
            "request": payload,
            "response": response.json(),
        }
        path.write_bytes(gzip.compress(json.dumps(result).encode(), mtime=0))
    body = result["response"]
    if body.get("errors"):
        raise ValueError(body["errors"])
    value = body["data"]["variant"]
    assert value["variant_id"] == allele
    summary = {}
    for assay in ["exome", "genome", "joint"]:
        record = value[assay]
        if record is None:
            summary[assay] = None
            continue
        all_sex = [p for p in record["populations"] if p["id"] in {"XX", "XY"}]
        unique = {}
        for row in all_sex:
            if row["id"] in unique:
                assert row == unique[row["id"]], (assay, "conflicting sex duplicates")
            unique[row["id"]] = row
        sex = list(unique.values())
        summary[assay] = {
            **{k: record[k] for k in ["ac", "an", "homozygote_count"]},
            "sex": sex,
            "identical_duplicate_sex_rows_removed": len(all_sex) - len(sex),
        }
        assert len(sex) == 2
        for field in ["ac", "an", "homozygote_count"]:
            assert sum(p[field] for p in sex) == record[field], (assay, field)
    receipt = {
        "variant_id": allele,
        "assays": summary,
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source": str(path.relative_to(REPO)),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fetched_at": result["fetched_at"],
        "interpretation": "XX/XY are inferred sex strata, not observed cancer status.",
    }
    return receipt


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    alleles = [
        "13-32363389-G-T",  # K2729N
        "13-32363369-G-C",  # D2723H
        "13-32380043-C-T",  # R3052W
        "13-32362595-G-C",  # W2626C
        "13-32332592-A-C",  # N372H
        "13-32363225-A-G",  # I2675V
    ]
    receipts = []
    for i, allele in enumerate(alleles):
        if i:
            time.sleep(7)
        record = fetch_one(allele)
        receipts.append(record)
        print(allele, json.dumps(record["assays"]["joint"]), flush=True)
    (HERE / "sex_probe.json").write_text(json.dumps(receipts, indent=2) + "\n")


if __name__ == "__main__":
    main()
