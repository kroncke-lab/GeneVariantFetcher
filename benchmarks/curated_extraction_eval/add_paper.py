#!/usr/bin/env python3
"""Add a paper to the curated benchmark in one step.

    python add_paper.py <GENE> <PMID> --strategy table --note "why it's here"

Appends a row to registry.tsv, reports whether the paper has a gold answer and
cached source, then rebuilds the fixture. Idempotent: re-adding an existing
(gene, pmid) is a no-op. Use this when you discover a problem paper you want the
benchmark to start guarding against.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import build_fixture as bf  # noqa: E402


def registry_has(gene: str, pmid: str) -> bool:
    for p in bf.load_registry():
        if p["gene"] == gene and p["pmid"] == pmid:
            return True
    return False


def gold_location(gene: str, pmid: str) -> str:
    if pmid in bf._read_gold_csv(bf.GOLD_SRC / f"{gene}_recall_input.csv"):
        return "repo gold standard"
    if pmid in bf._read_gold_csv(bf.OVERRIDES_DIR / f"{gene}_recall_input.csv"):
        return "gold_overrides"
    return ""


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("gene")
    ap.add_argument("pmid")
    ap.add_argument("--strategy", default="table", choices=sorted(bf.VALID_STRATEGIES))
    ap.add_argument(
        "--note", default="", help="One line on why this paper is in the set."
    )
    ap.add_argument(
        "--no-build", action="store_true", help="Append only; don't rebuild."
    )
    args = ap.parse_args()

    gene = args.gene.strip().upper()
    pmid = args.pmid.strip()
    if not pmid.isdigit():
        raise SystemExit(f"PMID should be digits, got {pmid!r}")

    if registry_has(gene, pmid):
        print(f"{gene}/{pmid} is already in registry.tsv — nothing to add.")
        return 0

    note = args.note or "added via add_paper.py (fill in why)"
    line = f"{gene}\t{pmid}\t{args.strategy}\t{note}\n"
    text = bf.REGISTRY.read_text()
    if text and not text.endswith("\n"):
        text += "\n"
    bf.REGISTRY.write_text(text + line)
    print(f"Appended to registry.tsv:\n  {line.strip()}")

    # Status checks
    gold = gold_location(gene, pmid)
    if gold:
        print(f"  gold: found in {gold} ✓")
    else:
        print(
            f"  gold: NONE — supply expected variants in "
            f"gold_overrides/{gene}_recall_input.csv (see gold_overrides/README.md),"
        )
        print(f"        otherwise this paper cannot be scored.")
    src = bf.CORPUS / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
    if src.exists():
        print(f"  source: cached at {src.relative_to(bf.REPO)} ✓")
    else:
        print(
            f"  source: NOT cached at corpus/{gene}/{pmid}/ — fetch it (e.g. a "
            f"gvf-run that includes this PMID) before relying on extract mode. "
            f"score mode can still score an existing DB without local source."
        )

    if args.no_build:
        print(
            "\n--no-build: skipping rebuild. Run `python build_fixture.py` when ready."
        )
        return 0

    print("\nRebuilding fixture...")
    rc = bf.main([])
    print("\nNext: `python run_benchmark.py` to score the updated set.")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
