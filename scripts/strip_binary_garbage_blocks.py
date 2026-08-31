#!/usr/bin/env python3
"""Strip provably-binary supplement blocks out of FULL_CONTEXT base text.

``fold_supplements_into_full_context`` regenerates only the sentinel-delimited
fold block, deliberately preserving harvest-time base content. But harvest-time
assembly itself appended ``# SUPPLEMENTAL FILE N:`` blocks directly, and before
the converter garbage gate landed, a failed ``.doc``/``.xml`` conversion could
put megabytes of raw container bytes there (LMNA 19220582 carried 57 MB of it;
six corpus papers were found in the 2026-08-31 audit). Those blocks sit outside
the fold sentinel, so a re-fold cannot remove them.

This script removes exactly the base supplement blocks whose text fails
:func:`harvesting.format_converters.looks_like_binary_garbage` — the same
detector the converters and the fold now use. Inside a failing block it first
tries ``Nested file:`` granularity so clean members of an expanded archive
survive. The sentinel fold block is never touched (a re-fold already
regenerates it cleanly), the original file is backed up once to
``*.pre_degarbage_bak``, and the on-disk supplement files are never deleted —
a future, better converter can re-fold them.

Usage::

    python scripts/strip_binary_garbage_blocks.py --corpus corpus --genes LMNA
    python scripts/strip_binary_garbage_blocks.py --harvest-dir <dir> --dry-run
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from harvesting.format_converters import looks_like_binary_garbage  # noqa: E402
from harvesting.supplement_fold import (  # noqa: E402
    FOLD_BEGIN,
    _NESTED_FILE_HEADING_RE,
    _SUPPLEMENT_HEADING_RE,
)
from scripts.fold_supplements import discover_corpus_papers, discover_pmids  # noqa: E402

logger = logging.getLogger("strip_binary_garbage_blocks")


def _strip_garbage_nested_files(block: str) -> tuple[str, list[str]]:
    """Remove only the failing ``Nested file:`` sub-blocks of one block."""
    matches = list(_NESTED_FILE_HEADING_RE.finditer(block))
    if not matches:
        return block, []
    removed: list[str] = []
    kept_parts: list[str] = [block[: matches[0].start()]]
    for idx, match in enumerate(matches):
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(block)
        sub = block[match.start() : end]
        label = match.group("label").strip()
        if looks_like_binary_garbage(sub):
            removed.append(label)
        else:
            kept_parts.append(sub)
    return "".join(kept_parts), removed


def strip_garbage_supplement_blocks(text: str) -> tuple[str, list[str]]:
    """Return ``(cleaned_text, removed_labels)`` for one FULL_CONTEXT body.

    Only base content is considered: the sentinel fold block (``FOLD_BEGIN``
    onward) is preserved verbatim, since a re-fold regenerates it with the
    converter garbage gate already applied. A base block spans from its
    ``# SUPPLEMENTAL FILE N:`` heading to the next heading, the fold sentinel,
    or end of text — the same span rule the fold module uses.
    """
    sentinel_start = text.find(FOLD_BEGIN)
    base = text if sentinel_start == -1 else text[:sentinel_start]
    tail = "" if sentinel_start == -1 else text[sentinel_start:]

    matches = list(_SUPPLEMENT_HEADING_RE.finditer(base))
    if not matches:
        return text, []

    removed: list[str] = []
    kept_parts: list[str] = [base[: matches[0].start()]]
    for idx, match in enumerate(matches):
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(base)
        block = base[match.start() : end]
        label = match.group("label").strip()
        if not looks_like_binary_garbage(block):
            kept_parts.append(block)
            continue
        # Prefer nested-file granularity: an expanded archive block may hold
        # one leaked member among clean ones.
        cleaned_block, nested_removed = _strip_garbage_nested_files(block)
        if nested_removed and not looks_like_binary_garbage(cleaned_block):
            kept_parts.append(cleaned_block)
            removed.extend(f"{label} :: {name}" for name in nested_removed)
        else:
            removed.append(label)

    if not removed:
        return text, []
    cleaned = "".join(kept_parts).rstrip()
    if tail:
        cleaned = cleaned + "\n\n" + tail.lstrip()
    return cleaned + "\n", removed


def process_paper(pmid: str, harvest_dir: Path, dry_run: bool) -> bool:
    """Degarbage one paper's FULL_CONTEXT. Returns True when it changed."""
    fc = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    if not fc.is_file():
        return False
    original = fc.read_text(encoding="utf-8", errors="replace")
    cleaned, removed = strip_garbage_supplement_blocks(original)
    if not removed:
        return False
    logger.info(
        "%s%s: removing %d binary-garbage block(s), %.2f MB -> %.2f MB: %s",
        "[dry-run] " if dry_run else "",
        pmid,
        len(removed),
        len(original) / 1e6,
        len(cleaned) / 1e6,
        "; ".join(removed),
    )
    if dry_run:
        return True
    backup = fc.parent / (fc.name + ".pre_degarbage_bak")
    if not backup.exists():
        backup.write_text(original, encoding="utf-8")
    fc.write_text(cleaned, encoding="utf-8")
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--harvest-dir", type=Path, help="Dir with *_FULL_CONTEXT.md.")
    parser.add_argument(
        "--corpus", type=Path, help="Nested corpus root of <GENE>/<PMID>/ dirs."
    )
    parser.add_argument(
        "--genes", default="", help="Comma-separated genes in --corpus mode."
    )
    parser.add_argument(
        "--pmids-file", type=Path, help="Optional file of PMIDs to restrict to."
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )

    only_pmids: set[str] | None = None
    if args.pmids_file:
        only_pmids = {
            line.strip()
            for line in args.pmids_file.read_text().splitlines()
            if line.strip()
        }

    papers: list[tuple[str, Path]] = []
    if args.corpus:
        if not args.corpus.is_dir():
            parser.error(f"--corpus not a directory: {args.corpus}")
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        papers = discover_corpus_papers(args.corpus, genes)
    elif args.harvest_dir:
        if not args.harvest_dir.is_dir():
            parser.error(f"--harvest-dir not a directory: {args.harvest_dir}")
        papers = [(p, args.harvest_dir) for p in discover_pmids(args.harvest_dir)]
    else:
        parser.error("one of --corpus or --harvest-dir is required")

    changed = 0
    for pmid, paper_dir in papers:
        if only_pmids is not None and pmid not in only_pmids:
            continue
        if process_paper(pmid, paper_dir, args.dry_run):
            changed += 1
    logger.info("degarbaged %d FULL_CONTEXT file(s)", changed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
