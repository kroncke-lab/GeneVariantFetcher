"""Populate gvf_data/reference_sequences/<GENE>.fasta from NCBI (Entrez).

Operational step for the reference-transcript validation gate
(pipeline/reference_validation.py). Reads
gvf_data/reference_sequences/manifest.json,
fetches each gene's canonical RefSeq protein via Entrez efetch, validates the
fetched length against the manifest's expected_length (which matches
utils.variant_normalizer.PROTEIN_LENGTHS), and writes the FASTA only when the
length checks out — a wrong/stale accession is reported and skipped, never
trusted.

Requires network and an NCBI email (NCBI_EMAIL or ENTREZ_EMAIL); biopython is
already a project dependency. Not part of any automated flow — run once to seed
the cache, then the gate works offline.

``--exons`` also fetches each gene's RefSeq mRNA feature table and writes
``exon_maps.json`` so structural matching can equate ``EXON N DELETION`` with
the in-frame protein range of that coding exon.

Examples::

    python scripts/fetch_reference_sequences.py            # all manifest genes
    python scripts/fetch_reference_sequences.py --genes KCNH2 SCN5A
    python scripts/fetch_reference_sequences.py --exons    # rebuild exon maps
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

logger = logging.getLogger("fetch_reference_sequences")

_REF_DIR = Path(__file__).resolve().parents[1] / "gvf_data" / "reference_sequences"
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


def _fetch_feature_table(accession: str, *, email: str, api_key: str | None) -> str:
    from Bio import Entrez

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="ft", retmode="text")
    try:
        return handle.read()
    finally:
        handle.close()


def _write_exon_maps(
    manifest: dict, genes: list[str], *, email: str, api_key: str | None
) -> int:
    from pipeline.reference_validation import load_reference_protein
    from utils.structural_alleles import build_exon_records, parse_ncbi_feature_table

    maps: dict = {}
    failed = 0
    for gene in genes:
        gene = gene.upper()
        entry = manifest.get(gene)
        if not entry or not entry.get("refseq_mrna"):
            logger.warning("%s: no refseq_mrna in manifest; skipping exon map", gene)
            failed += 1
            continue
        sequence = load_reference_protein(gene)
        if not sequence:
            logger.warning("%s: no cached protein FASTA; skipping exon map", gene)
            failed += 1
            continue
        try:
            table = _fetch_feature_table(
                entry["refseq_mrna"], email=email, api_key=api_key
            )
        except Exception as exc:  # noqa: BLE001
            logger.error("%s: feature-table fetch failed: %s", gene, exc)
            failed += 1
            continue
        exons, cds = parse_ncbi_feature_table(table)
        if not exons or cds is None:
            logger.error("%s: feature table missing exon/CDS spans", gene)
            failed += 1
            continue
        records = build_exon_records(exons, cds, sequence)
        maps[gene] = {
            "transcript": entry["refseq_mrna"],
            "protein": entry["refseq_protein"],
            "protein_length": len(sequence),
            "cds": {
                "transcript_start": cds[0],
                "transcript_end": cds[1],
                "coding_nt": sum(
                    max(0, min(end, cds[1]) - max(start, cds[0]) + 1)
                    for start, end in exons
                    if min(end, cds[1]) >= max(start, cds[0])
                ),
            },
            "exons": records,
        }
        in_frame = sum(1 for rec in records if rec.get("in_frame"))
        logger.info(
            "%s: %d exons (%d in-frame coding) from %s",
            gene,
            len(records),
            in_frame,
            entry["refseq_mrna"],
        )

    payload = {
        "_comment": (
            "Coding-exon protein spans for the canonical RefSeq transcripts in "
            "manifest.json. Built from NCBI feature tables (exon features in "
            "transcript order intersected with CDS). The terminal stop codon is "
            "excluded from protein_nt. in_frame exons may be equated to a "
            "protein-range deletion. Regenerate with "
            "scripts/fetch_reference_sequences.py --exons."
        ),
        "genes": maps,
    }
    out = _REF_DIR / "exon_maps.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    logger.info("wrote %s (%d genes)", out.name, len(maps))
    return failed


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
    p.add_argument(
        "--exons",
        action="store_true",
        help="Rebuild reference_sequences/exon_maps.json from NCBI feature tables.",
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
    if args.exons:
        failed += _write_exon_maps(
            manifest, [str(g).upper() for g in genes], email=email, api_key=api_key
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
