#!/usr/bin/env python3
"""Stage data-bearing corpus figures into the flat layout the figure-reader wants.

The corpus stores figures nested as ``corpus/<GENE>/<PMID>/<PMID>_figures/``,
but ``scripts/extract_figure_variants.py`` (and ``figure_variant_reader``) expect
a flat ``<pmc_dir>/<PMID>_figures/`` directory. This helper builds that flat view
with **symlinks** (no copying) and, crucially, links *only* the figures the
caption triage flagged as data-bearing (pedigrees / clinical panels) — so the
vision pass never spends a call on functional/structural images.

Reads ``corpus/figure_triage.csv`` (written by ``fetch_gold_figures.py``); a
figure is staged when its decision is ``fetch`` or ``have`` and the image is
actually on disk. No network, no LLM.

Usage::

    python scripts/stage_corpus_figures.py --gene KCNH2 --out /tmp/fig_stage/KCNH2
    python scripts/extract_figure_variants.py --gene KCNH2 \\
        --pmc-dir /tmp/fig_stage/KCNH2 --auto-pmids \\
        --out /tmp/fig_reads/KCNH2 --db /tmp/working/KCNH2.db
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def stage(gene: str, corpus: Path, triage_csv: Path, out: Path) -> tuple[int, int]:
    """Symlink data-bearing images for *gene* into ``out/<PMID>_figures/``.

    Returns ``(pmids_staged, images_staged)``.
    """
    by_pmid: dict[str, list[str]] = defaultdict(list)
    with triage_csv.open(newline="") as fh:
        for r in csv.DictReader(fh):
            if r["gene"].upper() != gene.upper():
                continue
            if r["decision"] not in ("fetch", "have"):
                continue
            src = corpus / r["gene"] / r["pmid"] / f"{r['pmid']}_figures" / r["image"]
            if src.exists():
                by_pmid[r["pmid"]].append(r["image"])

    out.mkdir(parents=True, exist_ok=True)
    n_img = 0
    for pmid, images in sorted(by_pmid.items()):
        fig_dir = out / f"{pmid}_figures"
        fig_dir.mkdir(exist_ok=True)
        for img in sorted(set(images)):
            src = corpus / gene / pmid / f"{pmid}_figures" / img
            link = fig_dir / img
            if link.exists() or link.is_symlink():
                continue
            link.symlink_to(src.resolve())
            n_img += 1
    return len(by_pmid), n_img


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument("--triage-csv", default=str(REPO / "corpus" / "figure_triage.csv"))
    ap.add_argument(
        "--out", required=True, help="Flat staging dir (gets <PMID>_figures/)."
    )
    args = ap.parse_args()

    triage = Path(args.triage_csv)
    if not triage.exists():
        raise SystemExit(
            f"triage CSV not found: {triage}\n"
            "Run scripts/fetch_gold_figures.py first to write it."
        )
    pmids, images = stage(args.gene.upper(), Path(args.corpus), triage, Path(args.out))
    print(
        f"{args.gene}: staged {images} data-bearing image(s) across {pmids} PMID(s) "
        f"-> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
