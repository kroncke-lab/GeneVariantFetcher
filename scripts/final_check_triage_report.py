#!/usr/bin/env python3
"""Offline shadow report for the final-check triage predicate (no LLM).

Runs the deterministic triage predicate over a migrated DB and prints the
full/cheap/skip distribution + the reasons, so the projected `@xhigh` cost
reduction and the risk of demoting completeness-critical papers can be judged
BEFORE spending on a live enforce run. This is the go/no-go artifact
(codex/grok recommendation): it never calls an LLM.

Usage:
    # Arbitrary migrated DB(s):
    python scripts/final_check_triage_report.py --db GENE=path [--db ...] \
        [--corpus ./corpus] [--per-paper]

    # The 101-paper curated staging set (the calibration go/no-go artifact):
    python scripts/final_check_triage_report.py --staging [--per-paper] \
        [--zero-count-tier full|cheap]

In ``--staging`` mode the canonical per-gene DBs (the ones scored by
``benchmarks/curated_extraction_eval/run_benchmark.py``) are each restricted to
their manifest PMIDs, and ``has_source`` is taken from the manifest
(``source_status``) rather than probed from disk — every staging paper has
confirmed cached source, so this is the faithful signal the live pipeline sees.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def _has_source_fn(gene: str, corpus: Path):
    def has_source(pmid: str) -> bool:
        d = corpus / gene / pmid
        if not d.is_dir():
            return False
        return any(d.glob("*FULL_CONTEXT*.md")) or any(d.glob("*DATA_ZONES*.md"))

    return has_source


def report_db(
    gene: str,
    db: Path,
    corpus: Path,
    per_paper: bool,
    zero_count_tier: str = "full",
    *,
    pmid_filter: set[str] | None = None,
    source_pmids: set[str] | None = None,
) -> Counter:
    """Report the triage distribution for one gene's DB.

    ``pmid_filter`` restricts the report to a fixed PMID set (e.g. the staging
    manifest) so a large canonical DB is scored only over the papers of interest.
    ``source_pmids`` overrides the on-disk source probe with a known-good set (the
    manifest confirms every staging paper has cached source).
    """
    from pipeline.paper_final_check_triage import (
        collect_risk_views,
        decide_final_check_tier,
    )
    from pipeline.trust_gate import apply_trust_gate

    if source_pmids is not None:
        known = set(source_pmids)

        def has_source(pmid: str) -> bool:
            return pmid in known

    else:
        has_source = _has_source_fn(gene, corpus)

    # Work on a copy so the source DB is never mutated; ensure trust columns.
    with tempfile.TemporaryDirectory() as td:
        work = Path(td) / "work.db"
        shutil.copy(db, work)
        try:
            apply_trust_gate(work)
        except Exception as exc:  # pragma: no cover - diagnostic path
            print(f"  (trust gate skipped: {exc})")
        views = collect_risk_views(work, has_source=has_source)

    if pmid_filter is not None:
        keep = set(pmid_filter)
        missing = sorted(keep - set(views))
        views = {p: v for p, v in views.items() if p in keep}
        if missing:
            print(
                f"  ({gene}: {len(missing)} staging PMID(s) absent from the canonical "
                f"DB, excluded from routing: {','.join(missing)})"
            )

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


def _staging_targets() -> list[tuple[str, Path, set[str], set[str]]]:
    """Return (gene, db_path, staging_pmids, source_pmids) for the 101-paper set.

    DB paths + manifest are the single source of truth in the benchmark runner, so
    reuse them rather than re-declaring canonical paths here.
    """
    bench = REPO / "benchmarks" / "curated_extraction_eval"
    sys.path.insert(0, str(bench))
    import run_benchmark as rb  # noqa: E402  (path injected above)

    manifest = rb.load_manifest()  # keyed "GENE:PMID"
    by_gene: dict[str, dict[str, set[str]]] = {}
    for row in manifest.values():
        gene = row["gene"].strip().upper()
        pmid = row["pmid"].strip()
        g = by_gene.setdefault(gene, {"pmids": set(), "source": set()})
        g["pmids"].add(pmid)
        if (row.get("source_status") or "").strip().lower() == "ok":
            g["source"].add(pmid)

    dbs = rb.resolve_dbs(None)  # existence-gated canonical DBs
    targets: list[tuple[str, Path, set[str], set[str]]] = []
    for gene in sorted(by_gene):
        db = dbs.get(gene)
        if db is None:
            print(f"  (skip {gene}: no canonical DB on this checkout)")
            continue
        targets.append((gene, db, by_gene[gene]["pmids"], by_gene[gene]["source"]))
    return targets


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", action="append", help="GENE=path (arbitrary DB mode)")
    ap.add_argument(
        "--staging",
        action="store_true",
        help="Report over the 101-paper curated staging set (canonical DBs "
        "restricted to the benchmark manifest PMIDs).",
    )
    ap.add_argument("--corpus", default="corpus")
    ap.add_argument("--per-paper", action="store_true")
    ap.add_argument(
        "--zero-count-tier",
        choices=("full", "cheap"),
        default="full",
        help="Lane for zero-count-with-source papers (the biggest cost lever).",
    )
    args = ap.parse_args()

    if not args.staging and not args.db:
        ap.error("pass --staging or at least one --db GENE=path")

    corpus = Path(args.corpus).expanduser()
    grand: Counter = Counter()

    if args.staging:
        for gene, db, pmids, source in _staging_targets():
            print(
                f"\n=== {gene} ({len(pmids)} staging papers) "
                f"[zero_count_tier={args.zero_count_tier}] ==="
            )
            grand.update(
                report_db(
                    gene,
                    db,
                    corpus,
                    args.per_paper,
                    zero_count_tier=args.zero_count_tier,
                    pmid_filter=pmids,
                    source_pmids=source,
                )
            )
    else:
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
