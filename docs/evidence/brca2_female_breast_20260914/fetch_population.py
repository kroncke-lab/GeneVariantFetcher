"""Fetch exact joint sex counts for the frozen canonical BRCA2 type inventory."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = REPO / "results/brca2_female_breast_20260914/population"
API = "https://gnomad.broadinstitute.org/api"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def inventory():
    path = HERE.parent / "population_inclusive_penetrance_20260912/rebuild_union.py"
    legacy = module("population_union", path)
    legacy.GENES = ["BRCA2"]
    population, sources = legacy.load_population()
    selected = (
        population.loc[
            population.population_eligible
            & population.vclass.isin(["missense", "nonsense", "stop_gained"])
            & population.canonical_wt_status.eq("match")
        ]
        .copy()
        .sort_values("variant_id")
    )
    assert selected.variant_id.is_unique
    return selected, legacy, sources + [path]


def fetch(ids):
    key = hashlib.sha256("\n".join(ids).encode()).hexdigest()[:20]
    path = RAW / f"{key}.json.gz"
    fields = [
        f'v{i}: variant(variantId: "{allele}", dataset: gnomad_r4) {{ variant_id joint {{ ac an homozygote_count populations {{ id ac an homozygote_count }} }} }}'
        for i, allele in enumerate(ids)
    ]
    request = {"query": "query FemaleCounts { " + " ".join(fields) + " }"}
    if path.exists():
        receipt = json.loads(gzip.decompress(path.read_bytes()))
        assert receipt["request"] == request
    else:
        for attempt in range(8):
            try:
                response = requests.post(API, json=request, timeout=90)
                response.raise_for_status()
                body = response.json()
                if body.get("errors"):
                    raise ValueError(body["errors"])
                receipt = dict(
                    api=API,
                    fetched_at=datetime.now(timezone.utc).isoformat(),
                    request=request,
                    response=body,
                )
                path.write_bytes(gzip.compress(json.dumps(receipt).encode(), mtime=0))
                break
            except (requests.RequestException, ValueError):
                if attempt == 7:
                    raise
                time.sleep(min(5 * 2**attempt, 60))
        time.sleep(5)
    body = receipt["response"]
    assert not body.get("errors")
    rows = []
    for i, allele in enumerate(ids):
        value = body["data"][f"v{i}"]
        assert value and value["variant_id"] == allele
        joint = value["joint"]
        assert joint is not None
        sex = {}
        repeats = 0
        for row in joint["populations"]:
            if row["id"] not in {"XX", "XY"}:
                continue
            if row["id"] in sex:
                assert row == sex[row["id"]], (allele, "conflicting sex rows")
                repeats += 1
            sex[row["id"]] = row
        assert set(sex) == {"XX", "XY"}
        result = dict(
            variant_id=allele,
            identical_sex_duplicates_removed=repeats,
            response_path=str(path.relative_to(REPO)),
        )
        for field in ["ac", "an", "homozygote_count"]:
            assert sum(s[field] for s in sex.values()) == joint[field], (allele, field)
            result["all_" + field] = joint[field]
            for label, counts in sex.items():
                result[label + "_" + field] = counts[field]
        for label in ["all", "XX", "XY"]:
            ac, hom = result[label + "_ac"], result[label + "_homozygote_count"]
            assert 0 <= 2 * hom <= ac <= result[label + "_an"]
            result[label + "_carriers"] = ac - hom
        rows.append(result)
    return rows, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--limit", type=int, help="Probe only; incomplete output is explicitly labeled."
    )
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    pop, _, sources = inventory()
    ids = pop.variant_id.tolist()
    if args.limit:
        ids = ids[: args.limit]
    batches = [ids[i : i + 20] for i in range(0, len(ids), 20)]
    rows, paths = [], []
    with ThreadPoolExecutor(max_workers=1) as pool:
        futures = [pool.submit(fetch, batch) for batch in batches]
        for future in as_completed(futures):
            data, path = future.result()
            rows.extend(data)
            paths.append(path)
            print(f"Exact sex counts: {len(rows)}/{len(ids)}", flush=True)
    counts = pd.DataFrame(rows).sort_values("variant_id")
    matched = counts.merge(pop, on="variant_id", validate="one_to_one")
    assert len(matched) == len(ids)
    for actual, expected in [
        ("all_ac", "joint_ac"),
        ("all_an", "joint_an"),
        ("all_homozygote_count", "joint_hom"),
        ("all_carriers", "gnomad_carriers"),
    ]:
        assert matched[actual].eq(matched[expected]).all(), (
            actual,
            "frozen count drift",
        )
    stem = "population_probe" if args.limit else "population_sex_counts"
    blob = counts.to_csv(index=False, lineterminator="\n").encode()
    (HERE / f"{stem}.csv.gz").write_bytes(gzip.compress(blob, mtime=0))
    hashes = {
        str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sources + sorted(paths) + [Path(__file__).resolve()]
    }
    (HERE / f"{stem}_receipt.json").write_text(
        json.dumps(
            dict(
                complete_inventory=len(ids) == len(pop),
                alleles=len(ids),
                expected_alleles=len(pop),
                all_sex_counts_match_frozen=True,
                xx_xy_reconcile=True,
                duplicate_rows_deduplicated=True,
                input_hashes=hashes,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
