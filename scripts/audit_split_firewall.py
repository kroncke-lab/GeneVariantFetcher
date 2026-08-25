#!/usr/bin/env python3
"""Audit a preregistered calibration/holdout packet tree for firewall breaches.

A multi-gene cohort shares papers: one PMID can be in the BRCA1 roster and the
BRCA2 roster at once. Split ranking is per gene, so an unpinned shared paper can
land in calibration for one gene and holdout for another. Tuning on the first
gene's calibration then burns the second gene's holdout before it is ever
scored, and the "never-tuned holdout" claim is false while every per-gene
manifest still looks correct on its own.

Two independent defects are reported:

``split`` breach
    The same PMID is calibration under one gene and holdout under another.

``source`` divergence
    The same PMID is bound to different source bytes under different genes, so
    the frozen SHA-256 guarantee is per gene-packet rather than per paper and
    curators are not reading the same document.

Exits 0 when clean, 1 when any breach is found, 2 on a usage error. It reads
manifests only and never writes, so it is safe to run against a sealed packet.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

SPLITS = ("calibration", "holdout")


def collect(root: Path) -> dict[str, dict[str, tuple[str, str, str]]]:
    """Map ``pmid -> gene -> (split, source_sha256, n_bytes)`` from manifests."""
    found: dict[str, dict[str, tuple[str, str, str]]] = defaultdict(dict)
    for manifest in sorted(root.glob("*/*/manifest.csv")):
        gene = manifest.parent.parent.name
        split = manifest.parent.name
        if split not in SPLITS:
            continue
        with manifest.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                pmid = (row.get("pmid") or "").strip()
                if not pmid:
                    continue
                found[pmid][gene] = (
                    (row.get("evaluation_split") or split).strip(),
                    (row.get("source_sha256") or "").strip(),
                    (row.get("n_bytes") or "").strip(),
                )
    return found


def audit(root: Path) -> tuple[list[str], list[str], int, int]:
    found = collect(root)
    shared = {pmid: genes for pmid, genes in found.items() if len(genes) > 1}
    split_breaches: list[str] = []
    source_breaches: list[str] = []
    for pmid, genes in sorted(shared.items()):
        splits = {gene: value[0] for gene, value in genes.items()}
        if len(set(splits.values())) > 1:
            detail = ", ".join(f"{gene}={splits[gene]}" for gene in sorted(splits))
            split_breaches.append(f"{pmid}: {detail}")
        shas = {value[1] for value in genes.values() if value[1]}
        if len(shas) > 1:
            detail = ", ".join(
                f"{gene}={genes[gene][2]}B/{genes[gene][1][:10]}"
                for gene in sorted(genes)
            )
            source_breaches.append(f"{pmid}: {detail}")
    return split_breaches, source_breaches, len(found), len(shared)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "root",
        type=Path,
        help="Packet tree root containing <GENE>/<split>/manifest.csv.",
    )
    ap.add_argument(
        "--allow-source-divergence",
        action="store_true",
        help="Report divergent per-gene source bytes without failing.",
    )
    args = ap.parse_args()
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        print(f"error: not a directory: {root}", file=sys.stderr)
        return 2

    split_breaches, source_breaches, total, shared = audit(root)
    if not total:
        print(f"error: no manifests under {root}", file=sys.stderr)
        return 2

    print(f"unique PMIDs: {total}    shared across genes: {shared}")
    if split_breaches:
        print(
            f"\nFIREWALL BREACH — {len(split_breaches)} PMID(s) in conflicting splits:"
        )
        for line in split_breaches:
            print(f"  {line}")
    else:
        print(
            "split firewall: OK (no PMID is calibration for one gene and holdout for another)"
        )

    if source_breaches:
        label = "WARNING" if args.allow_source_divergence else "SOURCE DIVERGENCE"
        print(
            f"\n{label} — {len(source_breaches)} shared PMID(s) bound to different bytes:"
        )
        for line in source_breaches:
            print(f"  {line}")
    else:
        print("source binding: OK (shared PMIDs agree on SHA-256)")

    failed = bool(split_breaches) or (
        bool(source_breaches) and not args.allow_source_divergence
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
