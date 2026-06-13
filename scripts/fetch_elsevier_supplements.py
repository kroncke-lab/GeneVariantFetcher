#!/usr/bin/env python3
"""Recover Elsevier mmc supplements that the full-text-API route never fetched.

The Elsevier full-text API returns body markdown only; a paper recovered via
that route lands without its ``mmc`` supplement files (where the mutation tables
usually live). This is the single largest supplement-recall bucket
(~31% of the gold-missing-variant gap; see docs/SUPPLEMENT_ACQUISITION_PLAN.md).

For each target paper this script: reads the DOI from the corpus
``{pmid}_artifacts.json``, fetches the authenticated full-text XML
(``harvesting/elsevier_api.ElsevierAPIClient``, API key + insttoken), extracts
the ``1-s2.0-<PII>-mmc<N>.<ext>`` references, and downloads them from the open
ScienceDirect CDN into ``corpus/<GENE>/<PMID>/<PMID>_supplements/``. The existing
re-fold (wired into ``scripts/refresh_run_db.py``) then makes them visible to
Tier-3 on the next refresh.

  # all Elsevier-DOI corpus papers for a gene that lack supplements on disk:
  python scripts/fetch_elsevier_supplements.py --gene SCN5A
  # a specific set:
  python scripts/fetch_elsevier_supplements.py --gene SCN5A --pmids 29325976,20129283
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from harvesting.elsevier_api import ElsevierAPIClient  # noqa: E402


def _doi_for(pdir: Path, pmid: str) -> str:
    art = pdir / f"{pmid}_artifacts.json"
    if not art.exists():
        return ""
    try:
        return (json.loads(art.read_text()).get("doi") or "").strip()
    except (json.JSONDecodeError, OSError):
        return ""


def _has_supplements_on_disk(pdir: Path, pmid: str) -> bool:
    sd = pdir / f"{pmid}_supplements"
    return sd.is_dir() and any(p.is_file() for p in sd.rglob("*"))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument("--pmids", default="", help="Comma-separated PMIDs (default: all).")
    ap.add_argument(
        "--include-with-supplements",
        action="store_true",
        help="Also (re)check papers that already have supplements on disk.",
    )
    ap.add_argument("--limit", type=int, default=0, help="Max papers to process.")
    args = ap.parse_args()

    client = ElsevierAPIClient(
        api_key=os.getenv("ELSEVIER_API_KEY"),
        insttoken=os.getenv("ELSEVIER_INSTTOKEN"),
    )
    if not client.is_available:
        raise SystemExit("ELSEVIER_API_KEY not set; cannot fetch supplements.")

    gene_dir = Path(args.corpus) / args.gene.upper()
    if not gene_dir.is_dir():
        raise SystemExit(f"No corpus dir for gene: {gene_dir}")

    want = {p.strip() for p in args.pmids.split(",") if p.strip()}
    targets: list[tuple[str, Path, str]] = []
    for pdir in sorted(gene_dir.iterdir()):
        if not pdir.is_dir():
            continue
        pmid = pdir.name
        if want and pmid not in want:
            continue
        doi = _doi_for(pdir, pmid)
        if not doi or not ElsevierAPIClient.is_elsevier_doi(doi):
            continue
        if not args.include_with_supplements and _has_supplements_on_disk(pdir, pmid):
            continue
        targets.append((pmid, pdir, doi))

    if args.limit:
        targets = targets[: args.limit]
    print(f"{args.gene}: {len(targets)} Elsevier-DOI paper(s) to check for supplements")

    papers_with_supp = 0
    total_files = 0
    for pmid, pdir, doi in targets:
        xml, err = client.get_fulltext_by_doi(doi)
        if not xml:
            print(f"  {pmid} ({doi}): no XML ({err})")
            continue
        saved = client.download_supplements(xml, pdir / f"{pmid}_supplements")
        if saved:
            papers_with_supp += 1
            total_files += len(saved)
            print(
                f"  {pmid} ({doi}): +{len(saved)} supplement(s) {[p.name for p in saved]}"
            )
        else:
            print(f"  {pmid} ({doi}): no mmc references in XML")

    print(
        f"\nDONE: recovered {total_files} supplement file(s) across "
        f"{papers_with_supp}/{len(targets)} papers. Re-fold + re-extract "
        f"(scripts/refresh_run_db.py folds automatically) to realize the variants."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
