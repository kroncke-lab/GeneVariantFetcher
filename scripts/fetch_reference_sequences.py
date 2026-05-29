"""Populate data/reference_sequences/<GENE>.fasta from NCBI (Entrez).

Operational step for the reference-transcript validation gate
(pipeline/reference_validation.py). Reads data/reference_sequences/manifest.json,
fetches each gene's canonical RefSeq protein via Entrez efetch, validates the
fetched length against the manifest's expected_length (which matches
utils.variant_normalizer.PROTEIN_LENGTHS), and writes the FASTA only when the
length checks out — a wrong/stale accession is reported and skipped, never
trusted.

Requires network and an NCBI email (NCBI_EMAIL or ENTREZ_EMAIL); biopython is
already a project dependency. Not part of any automated flow — run once to seed
the cache, then the gate works offline.

Examples::

    python scripts/fetch_reference_sequences.py            # all manifest genes
    python scripts/fetch_reference_sequences.py --genes KCNH2 SCN5A
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

logger = logging.getLogger("fetch_reference_sequences")

_REF_DIR = Path(__file__).resolve().parents[1] / "data" / "reference_sequences"
_MANIFEST = _REF_DIR / "manifest.json"


def _load_manifest() -> dict:
    data = json.loads(_MANIFEST.read_text(encoding="utf-8"))
    return data.get("genes", {})


def _fetch_protein_fasta(accession: str, *, email: str, api_key: str | None) -> str:
    from Bio import Entrez  # imported lazily; only this script needs biopython

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    try:
        return handle.read()
    finally:
        handle.close()


def _sequence_from_fasta(text: str) -> str:
    return "".join(
        line.strip() for line in text.splitlines() if line and not line.startswith(">")
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--genes",
        nargs="*",
        default=None,
        help="Subset of manifest genes (default: all).",
    )
    p.add_argument(
        "--email", default=None, help="NCBI email (else NCBI_EMAIL / ENTREZ_EMAIL)."
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-fetch even if the FASTA already exists.",
    )
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    email = args.email or os.environ.get("NCBI_EMAIL") or os.environ.get("ENTREZ_EMAIL")
    if not email:
        build_parser().error(
            "set --email or NCBI_EMAIL/ENTREZ_EMAIL (NCBI requires it)"
        )
    api_key = os.environ.get("NCBI_API_KEY")

    manifest = _load_manifest()
    genes = args.genes or sorted(manifest)
    _REF_DIR.mkdir(parents=True, exist_ok=True)

    ok, skipped, failed = 0, 0, 0
    for gene in genes:
        gene = gene.upper()
        entry = manifest.get(gene)
        if not entry:
            logger.warning("%s: not in manifest; skipping", gene)
            skipped += 1
            continue
        out = _REF_DIR / f"{gene}.fasta"
        if out.exists() and not args.force:
            logger.info(
                "%s: %s already present (use --force to refetch)", gene, out.name
            )
            skipped += 1
            continue

        accession = entry["refseq_protein"]
        expected = entry.get("expected_length")
        try:
            fasta = _fetch_protein_fasta(accession, email=email, api_key=api_key)
        except Exception as exc:  # noqa: BLE001
            logger.error("%s: fetch failed for %s: %s", gene, accession, exc)
            failed += 1
            continue

        seq = _sequence_from_fasta(fasta)
        if expected and len(seq) != expected:
            logger.error(
                "%s: %s length %d != expected %d — wrong/stale accession, NOT writing",
                gene,
                accession,
                len(seq),
                expected,
            )
            failed += 1
            continue

        out.write_text(fasta, encoding="utf-8")
        logger.info("%s: wrote %s (%d aa) from %s", gene, out.name, len(seq), accession)
        ok += 1

    logger.info(
        "reference sequences: %d written, %d skipped, %d failed", ok, skipped, failed
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
