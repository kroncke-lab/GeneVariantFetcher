#!/usr/bin/env python3
"""QC: count papers whose on-disk supplements are NOT reflected in FULL_CONTEXT.

The re-extraction path discovers sources by globbing ``*_FULL_CONTEXT.md`` and
never re-reads ``{pmid}_supplements/``. A paper that has convertible supplement
files on disk but whose FULL_CONTEXT shows no sign of them is "genuinely
unparsed" — its supplement tables are invisible to Tier-3. This script counts
that gap per gene so the fold-wiring in ``scripts/refresh_run_db.py`` has a
denominator to be measured against.

Works on either layout:
  * the nested corpus: ``corpus/<GENE>/<PMID>/{<PMID>_FULL_CONTEXT.md, <PMID>_supplements/}``
  * a flat harvest dir: ``<dir>/{<PMID>_FULL_CONTEXT.md, <PMID>_supplements/}``

No network, no LLM. Reuses ``harvesting.supplement_fold`` to decide which files
are convertible (i.e. could carry variant tables).

  python scripts/recall_audit/supplement_fold_gap.py                 # corpus, all genes
  python scripts/recall_audit/supplement_fold_gap.py --harvest-dir <run>/pmc_fulltext
  python scripts/recall_audit/supplement_fold_gap.py --list          # print the unparsed PMIDs
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

from harvesting.supplement_fold import (  # noqa: E402
    FOLD_BEGIN,
    _CONVERTIBLE_SUFFIXES,
)

# FULL_CONTEXT markers that indicate supplement content is already present —
# either from our re-fold (FOLD_BEGIN) or the harvest-time fold ("# SUPPLEMENTAL
# FILE ...", emitted by both paths).
_PARSED_MARKERS = (FOLD_BEGIN, "# SUPPLEMENTAL FILE", "FOLDED SUPPLEMENTS")


def _has_convertible_supplements(supp_dir: Path) -> int:
    if not supp_dir.is_dir():
        return 0
    return sum(
        1
        for p in supp_dir.rglob("*")
        if p.is_file()
        and p.suffix.lower() in _CONVERTIBLE_SUFFIXES
        and "__MACOSX" not in p.parts
        and not p.name.startswith(".")
    )


def _full_context_reflects_supplements(fc: Path) -> bool:
    try:
        text = fc.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return any(m in text for m in _PARSED_MARKERS)


def _paper_dirs(corpus: Path, harvest_dir: Path | None, genes: list[str]):
    """Yield (gene, pmid, full_context_path, supplements_dir)."""
    if harvest_dir is not None:
        for fc in sorted(harvest_dir.glob("*_FULL_CONTEXT.md")):
            pmid = fc.name.replace("_FULL_CONTEXT.md", "")
            yield "(flat)", pmid, fc, harvest_dir / f"{pmid}_supplements"
        return
    for gene in genes:
        gdir = corpus / gene
        if not gdir.is_dir():
            continue
        for pdir in sorted(gdir.iterdir()):
            if not pdir.is_dir():
                continue
            pmid = pdir.name
            fc = pdir / f"{pmid}_FULL_CONTEXT.md"
            if fc.is_file():
                yield gene, pmid, fc, pdir / f"{pmid}_supplements"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument(
        "--harvest-dir",
        default=None,
        help="Flat run harvest dir instead of the corpus (overrides --corpus).",
    )
    ap.add_argument(
        "--genes",
        default="",
        help="Comma-separated genes (corpus mode; default: all corpus genes).",
    )
    ap.add_argument("--list", action="store_true", help="Print the unparsed PMIDs.")
    ap.add_argument("--out", default="", help="Optional CSV of unparsed papers.")
    args = ap.parse_args()

    corpus = Path(args.corpus)
    harvest_dir = Path(args.harvest_dir) if args.harvest_dir else None
    if harvest_dir is None:
        genes = (
            [g.strip().upper() for g in args.genes.split(",") if g.strip()]
            if args.genes
            else sorted(p.name for p in corpus.iterdir() if p.is_dir())
            if corpus.is_dir()
            else []
        )
    else:
        genes = []

    have_supp = Counter()
    unparsed = Counter()
    unparsed_rows: list[dict] = []
    for gene, pmid, fc, supp_dir in _paper_dirs(corpus, harvest_dir, genes):
        n_conv = _has_convertible_supplements(supp_dir)
        if n_conv == 0:
            continue
        have_supp[gene] += 1
        if not _full_context_reflects_supplements(fc):
            unparsed[gene] += 1
            unparsed_rows.append(
                {"gene": gene, "pmid": pmid, "convertible_supp_files": n_conv}
            )

    print("Supplement-fold gap (papers with convertible supplements on disk):")
    print(f"  {'gene':10s} {'have_supp':>10s} {'UNPARSED':>10s}")
    for g in sorted(set(have_supp) | set(unparsed)):
        print(f"  {g:10s} {have_supp[g]:>10d} {unparsed[g]:>10d}")
    print(
        f"  {'TOTAL':10s} {sum(have_supp.values()):>10d} {sum(unparsed.values()):>10d}"
        "   <- the no-network fold opportunity (these PMIDs lose supplement tables on replay)"
    )

    if args.list:
        for r in unparsed_rows:
            print(
                f"    UNPARSED {r['gene']}/{r['pmid']} ({r['convertible_supp_files']} files)"
            )
    if args.out and unparsed_rows:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(unparsed_rows[0].keys()))
            w.writeheader()
            w.writerows(unparsed_rows)
        print(f"  wrote {len(unparsed_rows)} unparsed papers -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
