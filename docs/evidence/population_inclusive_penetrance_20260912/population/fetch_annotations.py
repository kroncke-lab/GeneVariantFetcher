"""Fetch canonical gnomAD annotation context for frozen, gene-wide populations."""

from datetime import UTC, datetime
import gzip
import json
from pathlib import Path
import shutil
import time

import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RAW = REPO / "results/population_inclusive_penetrance_20260912/raw"
PRIOR = REPO.parent / "VariantFeatures/data/gnomad/runs/r01_browser_20260911"
TRANSCRIPTS = {
    "GCK": "ENST00000403799",
    "HNF1A": "ENST00000257555",
    "LDLR": "ENST00000558518",
    "BRCA2": "ENST00000380152",
    "KCNQ1": "ENST00000155840",
}
ACCESSIONS = {
    "HNF1A": "P20823",
    "GCK": "P35557",
    "LDLR": "P01130",
    "BRCA2": "P51587",
    "KCNQ1": "P51787",
}
API = "https://gnomad.broadinstitute.org/api"
ANNOTATION = """
variant_id consequence transcript_id transcript_version hgvsc hgvsp
consequence_in_canonical_transcript
transcript_consequence {
  consequence_terms major_consequence gene_id gene_symbol
  transcript_id transcript_version is_canonical is_mane_select
  hgvsc hgvsp lof lof_filter lof_flags
}
"""
QUERY = """
query PopulationInclusiveAnnotation($symbol: String!, $transcript: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id gene_version symbol chrom start stop strand canonical_transcript_id
    mane_select_transcript { ensembl_id ensembl_version refseq_id refseq_version }
    exons { start stop }
    variants(dataset: gnomad_r4) {
      ANNOTATION
      exome { ac an homozygote_count filters }
      genome { ac an homozygote_count filters }
      joint { ac an homozygote_count filters }
    }
  }
  transcript(transcript_id: $transcript, reference_genome: GRCh38) {
    transcript_id transcript_version gene_id chrom start stop strand
    exons { start stop }
    variants(dataset: gnomad_r4) { ANNOTATION }
  }
}
""".replace("ANNOTATION", ANNOTATION)


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    for gene, accession in ACCESSIONS.items():
        path = RAW / f"{gene}_UniProt.json"
        if not path.exists():
            response = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=90
            )
            response.raise_for_status()
            path.write_bytes(response.content)
    for gene, transcript in TRANSCRIPTS.items():
        source = PRIOR / f"{gene}.json.gz"
        copied = RAW / f"{gene}_population_20260911.json.gz"
        if not copied.exists():
            shutil.copyfile(source, copied)
        destination = RAW / f"{gene}_annotations.json.gz"
        payload = {
            "query": QUERY,
            "variables": {"symbol": gene, "transcript": transcript},
        }
        (HERE / f"{gene}_annotation_request.json").write_text(
            json.dumps(payload, indent=2) + "\n"
        )
        if destination.exists():
            print(gene, "cached", flush=True)
            continue
        for attempt in range(5):
            response = requests.post(API, json=payload, timeout=180)
            if response.status_code in (429, 502, 503, 504) and attempt < 4:
                time.sleep(min(15 * (attempt + 1), 45))
                continue
            response.raise_for_status()
            body = response.json()
            if body.get("errors"):
                raise RuntimeError(body["errors"])
            break
        result = {
            "fetched_at": datetime.now(UTC).isoformat(),
            "api": API,
            "request": payload,
            "response": body,
        }
        destination.write_bytes(gzip.compress(json.dumps(result).encode(), mtime=0))
        print(
            gene,
            "gene_variants",
            len(body["data"]["gene"]["variants"]),
            "canonical_variants",
            len(body["data"]["transcript"]["variants"]),
            flush=True,
        )
        time.sleep(7)


if __name__ == "__main__":
    main()
