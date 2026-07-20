#!/usr/bin/env python3
"""Offline shadow report for the final-check triage predicate (no LLM).

Runs the deterministic triage predicate over a migrated DB and prints the
full/cheap/skip distribution + the reasons, so the projected `@xhigh` cost
reduction and the risk of demoting completeness-critical papers can be judged
BEFORE spending on a live enforce run. This is the go/no-go artifact
(codex/grok recommendation): it never calls an LLM.

Usage:
    python scripts/final_check_triage_report.py --db GENE=path [--db ...] \
        [--corpus ./corpus] [--per-paper]
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _has_source_fn(gene: str, corpus: Path):
    def has_source(pmid: str) -> bool:
        d = corpus / gene / pmid
        if not d.is_dir():
            return False
        return any(d.glob("*FULL_CONTEXT*.md")) or any(d.glob("*DATA_ZONES*.md"))

    return has_source


def report_db(
    gene: str, db: Path, corpus: Path, per_paper: bool, zero_count_tier: str = "full"
) -> Counter:
    from pipeline.paper_final_check_triage import (
        collect_risk_views,
        decide_final_check_tier,
    )
    from pipeline.trust_gate import apply_trust_gate

    # Work on a copy so the source DB is never mutated; ensure trust columns.
    with tempfile.TemporaryDirectory() as td:
        work = Path(td) / "work.db"
        shutil.copy(db, work)
        try:
            apply_trust_gate(work)
        except Exception as exc:  # pragma: no cover - diagnostic path
            print(f"  (trust gate skipped: {exc})")
        views = collect_risk_views(work, has_source=_has_source_fn(gene, corpus))

    tiers: Counter = Counter()
    reasons: Counter = Counter()
    for pmid, view in sorted(views.items()):
        d = decide_final_check_tier(view, zero_count_source_tier=zero_count_tier)
        tiers[d.tier] += 1
        for r in d.reasons:
            reasons[r] += 1
        if per_paper:
            print(
                f"  {pmid:>10}  {d.tier:<5}  facts={view.n_count_facts:<3} "
                f"quar={int(view.any_quarantine)} src={int(view.has_usable_source)} "
                f"conf={view.extraction_confidence}  {','.join(d.reasons)}"
            )
    total = sum(tiers.values()) or 1
    full = tiers.get("full", 0)
    print(
        f"  {gene}: {total} papers -> full={full} ({full / total:.0%}) "
        f"cheap={tiers.get('cheap', 0)} skip={tiers.get('skip', 0)}"
    )
    print(f"    reasons: {dict(reasons)}")
    return tiers


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", action="append", required=True, help="GENE=path")
    ap.add_argument("--corpus", default="corpus")
    ap.add_argument("--per-paper", action="store_true")
    ap.add_argument(
        "--zero-count-tier",
        choices=("full", "cheap"),
        default="full",
        help="Lane for zero-count-with-source papers (the biggest cost lever).",
    )
    args = ap.parse_args()

    corpus = Path(args.corpus).expanduser()
    grand: Counter = Counter()
    for spec in args.db:
        gene, _, path = spec.partition("=")
        print(f"\n=== {gene} ({path}) [zero_count_tier={args.zero_count_tier}] ===")
        grand.update(
            report_db(
                gene.upper(),
                Path(path).expanduser(),
                corpus,
                args.per_paper,
                zero_count_tier=args.zero_count_tier,
            )
        )

    total = sum(grand.values()) or 1
    full = grand.get("full", 0)
    print(
        f"\nTOTAL: {total} papers -> full={full} ({full / total:.0%}, i.e. "
        f"~{1 - full / total:.0%} of @xhigh calls avoided before escalation) "
        f"cheap={grand.get('cheap', 0)} skip={grand.get('skip', 0)}"
    )


if __name__ == "__main__":
    main()
