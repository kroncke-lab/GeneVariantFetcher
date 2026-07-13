"""Fold on-disk ``{pmid}_supplements/`` files into their ``FULL_CONTEXT.md``.

Operational aid for the source re-parse path: supplements are converted at
harvest time but re-extraction/replay never re-reads ``{pmid}_supplements/``.
Run this before ``scripts/refresh_run_db.py`` to surface already-downloaded
supplement tables (no network). Non-destructive (backs up each FULL_CONTEXT once
to ``*.pre_fold_bak``) and idempotent. Gene scoping happens at parse time and the
refresh explosion gate guards the DB rebuild, so this only assembles text.

Examples::

    python scripts/fold_supplements.py --run-dir results/SCN5A/20260506_102238
    python scripts/fold_supplements.py --harvest-dir <dir> --pmids-file pmids.txt --dry-run
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from harvesting.supplement_fold import fold_supplements_into_full_context

logger = logging.getLogger("fold_supplements")

_SUPP_SUFFIX = "_supplements"


def discover_pmids(harvest_dir: Path) -> list[str]:
    """PMIDs that have BOTH a FULL_CONTEXT.md and a non-empty supplements dir."""
    pmids: set[str] = set()
    for d in harvest_dir.glob(f"*{_SUPP_SUFFIX}"):
        if d.is_dir():
            pmids.add(d.name[: -len(_SUPP_SUFFIX)])
    return sorted(p for p in pmids if (harvest_dir / f"{p}_FULL_CONTEXT.md").is_file())


def discover_corpus_papers(corpus: Path, genes: list[str]) -> list[tuple[str, Path]]:
    """Return ``(pmid, paper_dir)`` pairs from the nested corpus layout."""
    papers: list[tuple[str, Path]] = []
    selected_genes = genes or sorted(
        path.name for path in corpus.iterdir() if path.is_dir()
    )
    for gene in selected_genes:
        gene_dir = corpus / gene.upper()
        if not gene_dir.is_dir():
            continue
        for paper_dir in sorted(path for path in gene_dir.iterdir() if path.is_dir()):
            pmid = paper_dir.name
            if (paper_dir / f"{pmid}_FULL_CONTEXT.md").is_file() and (
                paper_dir / f"{pmid}{_SUPP_SUFFIX}"
            ).is_dir():
                papers.append((pmid, paper_dir))
    return papers


def _resolve_harvest_dir(args: argparse.Namespace) -> Path | None:
    if args.harvest_dir:
        return args.harvest_dir
    if args.run_dir:
        candidate = args.run_dir / "pmc_fulltext"
        return candidate if candidate.is_dir() else args.run_dir
    return None


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help="GVF run dir (uses its pmc_fulltext/).",
    )
    p.add_argument(
        "--harvest-dir",
        type=Path,
        default=None,
        help="Dir holding *_FULL_CONTEXT.md + *_supplements/.",
    )
    p.add_argument(
        "--corpus",
        type=Path,
        default=None,
        help="Nested corpus root containing <GENE>/<PMID>/ directories.",
    )
    p.add_argument(
        "--genes",
        default="",
        help="Comma-separated genes in --corpus mode (default: all).",
    )
    p.add_argument(
        "--pmids-file",
        type=Path,
        default=None,
        help="Optional file of PMIDs (one per line) to restrict folding to.",
    )
    p.add_argument(
        "--dry-run", action="store_true", help="Report counts; do not modify files."
    )
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    harvest_dir = _resolve_harvest_dir(args)
    if args.corpus is not None and harvest_dir is not None:
        parser.error("--corpus cannot be combined with --run-dir/--harvest-dir")
    if args.corpus is not None and not args.corpus.is_dir():
        parser.error(f"--corpus directory does not exist: {args.corpus}")
    if args.corpus is None and (harvest_dir is None or not harvest_dir.is_dir()):
        parser.error(
            "provide --corpus, --harvest-dir, or --run-dir containing pmc_fulltext/"
        )

    wanted: set[str] = set()
    if args.pmids_file:
        wanted = {
            line.strip()
            for line in args.pmids_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        }

    if args.corpus is not None:
        genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
        papers = discover_corpus_papers(args.corpus, genes)
        if wanted:
            papers = [(pmid, pdir) for pmid, pdir in papers if pmid in wanted]
    else:
        assert harvest_dir is not None
        pmids = discover_pmids(harvest_dir)
        if wanted:
            pmids = [pmid for pmid in pmids if pmid in wanted]
        papers = [(pmid, harvest_dir) for pmid in pmids]

    folded = 0
    for pmid, paper_dir in papers:
        if args.dry_run:
            supp = paper_dir / f"{pmid}{_SUPP_SUFFIX}"
            n = sum(1 for path in supp.rglob("*") if path.is_file())
            logger.info("[dry-run] %s: %d supplement file(s)", pmid, n)
            continue
        if fold_supplements_into_full_context(pmid, paper_dir) is not None:
            folded += 1

    if not args.dry_run:
        logger.info(
            "folded supplements into %d/%d FULL_CONTEXT file(s)",
            folded,
            len(papers),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
