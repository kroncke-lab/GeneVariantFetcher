#!/usr/bin/env python3
"""Bridge the nested corpus to the flat harvest layout the refresh path expects.

``scripts/refresh_run_db.py`` discovers sources by globbing a *flat* harvest dir
(``<dir>/{PMID}_FULL_CONTEXT.md`` + ``<dir>/{PMID}_supplements/``), but the
consolidated source home is the *nested* corpus
(``corpus/<GENE>/<PMID>/{PMID}_FULL_CONTEXT.md`` + ``<PMID>_supplements/``); the
per-run ``pmc_fulltext/`` trees were removed. This builds the flat view so a
refresh/re-score can run against the corpus.

The FULL_CONTEXT markdown is **copied** (so the refresh-time supplement fold
writes to the copy, never mutating the corpus); the bulky ``_supplements/`` and
``_figures/`` directories are **symlinked** (read-only, no duplication). The
condensed ``_CLEANED.md`` / ``_DATA_ZONES.md`` are intentionally NOT copied so
the (fold-grown) FULL_CONTEXT is unambiguously the extraction source.

  python scripts/corpus_to_harvest.py --gene SCN5A --out /tmp/bridge/SCN5A
  python scripts/corpus_to_harvest.py --gene SCN5A --pmids 29325976,20129283 --out /tmp/bridge/SCN5A
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def build(gene: str, corpus: Path, out: Path, pmids: set[str]) -> tuple[int, int]:
    gene_dir = corpus / gene.upper()
    if not gene_dir.is_dir():
        raise SystemExit(f"No corpus dir for gene: {gene_dir}")
    out.mkdir(parents=True, exist_ok=True)
    n_papers = n_supp = 0
    for pdir in sorted(gene_dir.iterdir()):
        if not pdir.is_dir():
            continue
        pmid = pdir.name
        if pmids and pmid not in pmids:
            continue
        fc = pdir / f"{pmid}_FULL_CONTEXT.md"
        if not fc.is_file():
            continue
        shutil.copy2(fc, out / f"{pmid}_FULL_CONTEXT.md")
        n_papers += 1
        for sub in (f"{pmid}_supplements", f"{pmid}_figures"):
            src = pdir / sub
            link = out / sub
            if src.is_dir():
                if link.exists() or link.is_symlink():
                    if link.is_symlink() or link.is_file():
                        link.unlink()
                    else:
                        shutil.rmtree(link)
                link.symlink_to(src.resolve())
                if sub.endswith("_supplements"):
                    n_supp += 1
    return n_papers, n_supp


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument("--out", required=True, help="Flat harvest dir to build.")
    ap.add_argument(
        "--pmids",
        default="",
        help="Comma-separated PMIDs (default: all corpus papers).",
    )
    args = ap.parse_args()
    pmids = {p.strip() for p in args.pmids.split(",") if p.strip()}
    n_papers, n_supp = build(
        args.gene.upper(), Path(args.corpus), Path(args.out), pmids
    )
    print(
        f"{args.gene}: bridged {n_papers} paper(s) ({n_supp} with supplements) "
        f"-> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
