"""Fetch all observed small variants within five exact genomic gene spans.

This is an additive full-locus extension, separate from the CDS+75bp snapshot.
Queries are bounded to 20 kb and adaptively split if an API limit is reached.
"""

from datetime import UTC, datetime
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import time

import requests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from build_population import check_counts, joint_filter

REPO = HERE.parents[4]
RAW = REPO / "results/population_inclusive_penetrance_20260912/raw/full_locus"
API = "https://gnomad.broadinstitute.org/api"
QUERY = """
query FullLocusPopulation($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    chrom start stop
    variants(dataset: gnomad_r4) {
      variant_id chrom pos ref alt
      exome { ac an homozygote_count filters }
      genome { ac an homozygote_count filters }
      joint { ac an homozygote_count filters }
    }
  }
}
"""
CHUNK_BP = 20000
_last_request_at = 0.0


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def fetch_region(gene, chrom, start, stop):
    global _last_request_at
    key = f"{gene}_{chrom}_{start}_{stop}"
    path = RAW / f"{key}.json.gz"
    payload = {
        "query": QUERY,
        "variables": {"chrom": chrom, "start": start, "stop": stop},
    }
    if path.exists():
        result = json.loads(gzip.decompress(path.read_bytes()))
        if result["request"] != payload:
            raise ValueError("Cached request scope differs")
        return [result]
    for attempt in range(5):
        delay = 7 - (time.monotonic() - _last_request_at)
        if delay > 0:
            time.sleep(delay)
        _last_request_at = time.monotonic()
        response = requests.post(API, json=payload, timeout=180)
        if response.status_code in (429, 502, 503, 504) and attempt < 4:
            time.sleep(min(15 * (attempt + 1), 45))
            continue
        response.raise_for_status()
        body = response.json()
        if body.get("errors"):
            errors = "; ".join(e.get("message", "") for e in body["errors"])
            save_json(RAW / f"{key}.error.json", {"request": payload, "response": body})
            if (
                any(
                    term in errors.lower()
                    for term in ("too many variants", "smaller region")
                )
                and stop > start
            ):
                midpoint = (start + stop) // 2
                return fetch_region(gene, chrom, start, midpoint) + fetch_region(
                    gene, chrom, midpoint + 1, stop
                )
            raise RuntimeError(errors)
        region = (body.get("data") or {}).get("region")
        if not region or (region["chrom"], region["start"], region["stop"]) != (
            chrom,
            start,
            stop,
        ):
            raise ValueError("Region response identity mismatch")
        result = {
            "gene": gene,
            "fetched_at": datetime.now(UTC).isoformat(),
            "api": API,
            "request": payload,
            "response": body,
        }
        path.write_bytes(gzip.compress(json.dumps(result).encode(), mtime=0))
        print(key, len(region["variants"]), "source rows", flush=True)
        return [result]
    raise RuntimeError("Request failed")


def row_from_variant(gene, variant, result, source_hashes):
    check_counts(variant)
    key = variant["variant_id"]
    chrom, pos, ref, alt = key.split("-")
    if (chrom, int(pos), ref, alt) != (
        variant["chrom"],
        variant["pos"],
        variant["ref"],
        variant["alt"],
    ):
        raise ValueError("Genomic allele identity disagrees with variant ID")
    region = result["request"]["variables"]
    if chrom != region["chrom"] or not region["start"] <= int(pos) <= region["stop"]:
        raise ValueError("Variant falls outside queried interval")
    row = {
        "gene": gene,
        "variant_id": key,
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "variant_type": "SNV"
        if len(ref) == len(alt) == 1
        else "MNV"
        if len(ref) == len(alt)
        else "indel",
        "gnomad_release": "4.1.1",
        "source_region": f"{chrom}:{region['start']}-{region['stop']}",
        "source_query_sha256": source_hashes[0],
        "source_response_sha256": source_hashes[1],
    }
    for assay in ("exome", "genome", "joint"):
        value = variant.get(assay)
        row[f"{assay}_present"] = value is not None
        for output, field in (("ac", "ac"), ("an", "an"), ("hom", "homozygote_count")):
            row[f"{assay}_{output}"] = value[field] if value is not None else ""
        row[f"{assay}_filters"] = (
            ";".join(value["filters"]) if value is not None else "ABSENT"
        )
    joint = variant.get("joint")
    row.update(
        joint_api_pass=joint is not None and not joint["filters"],
        reconstructed_joint_filter=joint_filter(variant),
        qc_pass=joint is not None and joint_filter(variant) == "PASS",
        gnomad_carriers=joint["ac"] - joint["homozygote_count"]
        if joint is not None
        else "",
        observed=joint is not None and joint["ac"] > 0,
    )
    row["population_eligible"] = row["observed"] and row["qc_pass"]
    return row


def encode_rows(rows):
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return gzip.compress(out.getvalue().encode(), mtime=0)


def encode_parts(rows):
    encoded = encode_rows(rows)
    if len(encoded) <= 1100000:
        return [(rows, encoded)]
    midpoint = len(rows) // 2
    if midpoint == 0:
        raise ValueError("One variant record exceeds artifact size limit")
    return encode_parts(rows[:midpoint]) + encode_parts(rows[midpoint:])


def build_gene(gene, chrom, start, stop, results, prior):
    cursor = start
    intervals = sorted(
        (r["request"]["variables"]["start"], r["request"]["variables"]["stop"])
        for r in results
    )
    for left, right in intervals:
        if left > cursor or left < start or right > stop:
            raise ValueError("Incomplete or out-of-scope interval coverage")
        cursor = max(cursor, right + 1)
    if cursor != stop + 1:
        raise ValueError("Final interval does not reach the gene boundary")
    by_id, duplicates, sources = {}, 0, []
    for result in results:
        params = result["request"]["variables"]
        raw_path = RAW / f"{gene}_{chrom}_{params['start']}_{params['stop']}.json.gz"
        sources.append(
            {
                "raw_file": str(raw_path.relative_to(REPO)),
                "sha256": sha(raw_path.read_bytes()),
                "fetched_at": result["fetched_at"],
                "region": params,
                "rows": len(result["response"]["data"]["region"]["variants"]),
            }
        )
        source_hashes = (
            sha(json.dumps(result["request"], sort_keys=True).encode()),
            sha(json.dumps(result["response"], sort_keys=True).encode()),
        )
        for variant in result["response"]["data"]["region"]["variants"]:
            row = row_from_variant(gene, variant, result, source_hashes)
            key = row["variant_id"]
            if key in by_id:
                relevant = [k for k in row if not k.startswith("source_")]
                if any(row[k] != by_id[key][k] for k in relevant):
                    raise ValueError(
                        f"Inconsistent duplicate allele across chunks: {key}"
                    )
                duplicates += 1
                continue
            by_id[key] = row
    overlap, missing, differences = 0, [], []
    for key, row in prior.items():
        pos = int(row["pos"])
        if row["gene"] != gene or not start <= pos <= stop:
            continue
        if key not in by_id:
            missing.append(key)
            continue
        overlap += 1
        current = by_id[key]
        for assay in ("exome", "genome", "joint"):
            for field in ("ac", "an", "hom", "filters"):
                col = f"{assay}_{field}"
                if str(current[col]) != row[col]:
                    differences.append(
                        {
                            "variant_id": key,
                            "column": col,
                            "prior": row[col],
                            "current": current[col],
                        }
                    )
    if missing or differences:
        save_json(
            HERE / f"{gene}_overlap_failure.json",
            {"missing": missing, "differences": differences},
        )
        raise ValueError(f"{gene} full-locus and coding-footprint sources disagree")
    rows = sorted(by_id.values(), key=lambda r: (r["pos"], r["ref"], r["alt"]))
    parts = encode_parts(rows)
    output_files = []
    for index, (part_rows, encoded) in enumerate(parts, 1):
        suffix = "" if len(parts) == 1 else f".part{index:03}"
        path = HERE / f"{gene}_population_variants{suffix}.csv.gz"
        path.write_bytes(encoded)
        output_files.append(
            {
                "file": path.name,
                "bytes": len(encoded),
                "sha256": sha(encoded),
                "rows": len(part_rows),
            }
        )
    eligible = [r for r in rows if r["population_eligible"]]
    report = {
        "gene": gene,
        "gnomad_release": "4.1.1",
        "dataset_selector": "gnomad_r4",
        "interval": {
            "chrom": chrom,
            "start": start,
            "stop": stop,
            "reference_genome": "GRCh38",
        },
        "source_rows": len(rows),
        "observed_joint_ac_positive": sum(r["observed"] for r in rows),
        "observed_qc_pass": len(eligible),
        "joint_carriers_qc_pass": sum(r["gnomad_carriers"] for r in eligible),
        "joint_ac_qc_pass": sum(r["joint_ac"] for r in eligible),
        "joint_hom_qc_pass": sum(r["joint_hom"] for r in eligible),
        "prior_footprint_overlap_rows": overlap,
        "all_overlap_alleles_counts_filters_exact_match": True,
        "deduplicated_identical_chunk_rows": duplicates,
        "successful_region_requests": len(results),
        "complete_interval_coverage_verified": True,
        "output_files": output_files,
        "output_bytes": sum(f["bytes"] for f in output_files),
        "source_scope": "All returned small variants whose POS is within the exact genomic gene span; no coding-consequence filter. Structural variants/CNVs require separate endpoints.",
        "annotation_policy": "No protein annotation inferred for region-only records. Join canonical annotations from the earlier snapshot by exact genomic allele ID; preserve unmapped state.",
        "carrier_policy": "joint AC minus joint homozygote count, once per allele; gnomAD carriers assumed unaffected",
        "qc_policy": "Every present assay must PASS; absent assay permitted. Preserve raw joint API filters separately.",
        "sources": sources,
    }
    save_json(HERE / f"{gene}_provenance.json", report)
    print(
        gene,
        "READY",
        len(rows),
        "source rows",
        len(eligible),
        "observed QC-pass",
        flush=True,
    )
    return report


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    with gzip.open(HERE.parent / "population_variants.csv.gz", "rt") as handle:
        prior_rows = list(csv.DictReader(handle))
    with (HERE.parent / "full_locus_availability.csv").open() as handle:
        interval_rows = {r["gene"]: r for r in csv.DictReader(handle)}
    save_json(
        HERE / "query_template.json",
        {"query": QUERY, "dataset_selector": "gnomad_r4", "chunk_bp": CHUNK_BP},
    )
    reports = []
    for gene in ("GCK", "HNF1A", "LDLR", "BRCA2", "KCNQ1"):
        interval = interval_rows[gene]
        chrom, start, stop = (
            interval["chrom"],
            int(interval["start_grch38"]),
            int(interval["stop_grch38"]),
        )
        results = []
        for left in range(start, stop + 1, CHUNK_BP):
            results.extend(
                fetch_region(gene, chrom, left, min(left + CHUNK_BP - 1, stop))
            )
        prior = {r["variant_id"]: r for r in prior_rows if r["gene"] == gene}
        reports.append(build_gene(gene, chrom, start, stop, results, prior))
        save_json(
            HERE / "run_progress.json",
            {"genes_completed": [r["gene"] for r in reports], "reports": reports},
        )
    columns = (
        "gene",
        "source_rows",
        "observed_joint_ac_positive",
        "observed_qc_pass",
        "joint_carriers_qc_pass",
        "joint_ac_qc_pass",
        "joint_hom_qc_pass",
        "prior_footprint_overlap_rows",
        "successful_region_requests",
        "all_overlap_alleles_counts_filters_exact_match",
        "complete_interval_coverage_verified",
        "deduplicated_identical_chunk_rows",
    )
    with (HERE / "full_locus_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: report[key] for key in columns} for report in reports)


if __name__ == "__main__":
    main()
