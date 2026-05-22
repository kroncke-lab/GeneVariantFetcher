"""PMID-level recall test for the tier-1/tier-2 filter stack.

Reads a gene's filter_progress.jsonl (produced by ``cli.discover`` /
``pipeline.steps.filter_papers``), splits PMIDs into search-universe,
tier-1-pass, and tier-1+tier-2-pass sets, then compares each against the
gold-standard PMID list at ``gene_variant_fetcher_gold_standard/normalized/``.

Recall is bounded by search coverage: a PMID that the PubMed/PubMind/EuropePMC
search never returned cannot pass any downstream filter. The report breaks the
recall budget into:

  search_coverage = |universe ∩ gold| / |gold|     # upper bound for any filter
  tier1_recall    = |tier1_pass ∩ gold| / |gold|
  filter_recall   = |tier2_pass ∩ gold| / |gold|   # PMIDs queued for download
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLD_DIR = REPO_ROOT / "gene_variant_fetcher_gold_standard" / "normalized"
GENES = ("KCNH2", "KCNQ1", "RYR2", "SCN5A")


def load_gold_pmids(gene: str, gold_dir: Path) -> Set[str]:
    path = gold_dir / f"{gene}_recall_input.csv"
    if not path.exists():
        raise FileNotFoundError(f"Gold file not found: {path}")
    pmids: Set[str] = set()
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            v = (row.get("pmid") or "").strip()
            if v and v.isdigit():
                pmids.add(v)
    return pmids


def parse_filter_progress(
    path: Path,
) -> Tuple[Set[str], Set[str], Set[str], Dict[str, int]]:
    """Return (universe, tier1_pass, tier2_pass, stats).

    tier1_pass: PMIDs whose tier1 decision was 'pass' (or who had no abstract;
    KeywordFilter is fail-open in that case and the record reflects this).
    tier2_pass: PMIDs whose final_decision is PASS (passes both tiers).
    stats: misc bookkeeping counts.
    """
    universe: Set[str] = set()
    tier1_pass: Set[str] = set()
    tier2_pass: Set[str] = set()
    stats = defaultdict(int)

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                stats["malformed_lines"] += 1
                continue
            pmid = str(rec.get("pmid", "")).strip()
            if not pmid:
                continue
            universe.add(pmid)
            stats["records_total"] += 1

            final = (rec.get("final_decision") or "").upper()
            tier1 = rec.get("tier1") or {}
            t1_dec = (tier1.get("decision") or "").lower()

            # Tier-1 perspective: did this PMID pass Tier 1?
            # KeywordFilter is fail-open when abstract missing — final decision
            # then comes from Tier 2 alone. Treat "no abstract" pass-through as
            # tier1=pass.
            if t1_dec == "pass":
                tier1_pass.add(pmid)
            elif not tier1 and final == "PASS":
                # Some records lack tier1 if tier1 was skipped/disabled;
                # treat them as having passed (or not run).
                tier1_pass.add(pmid)
            # Tier-2 / final perspective:
            if final == "PASS":
                tier2_pass.add(pmid)
                stats["final_pass"] += 1
            elif final == "FAIL":
                stats["final_fail"] += 1
            else:
                stats["final_other"] += 1

    return universe, tier1_pass, tier2_pass, dict(stats)


def find_default_filter_progress(gene: str) -> Optional[Path]:
    """Best-effort: locate the most recent filter_progress.jsonl for a gene.

    Search order: recall_metrics/pmid_filter_recall_* (current test outputs),
    validation_runs/turnkey_e2e_*/, validation_runs/*/, results/.
    """
    candidates: List[Path] = []
    for root in (
        REPO_ROOT / "recall_metrics",
        REPO_ROOT / "validation_runs",
        REPO_ROOT / "results",
    ):
        if not root.exists():
            continue
        candidates.extend(root.glob(f"**/{gene}/**/pmid_status/filter_progress.jsonl"))
        candidates.extend(
            root.glob(f"**/results/{gene}/**/pmid_status/filter_progress.jsonl")
        )
    if not candidates:
        return None
    # Pick the one with the most records (the most complete run) and break
    # ties by mtime (newest wins).
    candidates = list({c.resolve() for c in candidates})

    def score(p: Path) -> Tuple[int, float]:
        try:
            n = sum(1 for _ in p.open())
        except OSError:
            n = 0
        return (n, p.stat().st_mtime)

    candidates.sort(key=score, reverse=True)
    return candidates[0]


def compute_per_gene(
    gene: str,
    progress_path: Path,
    gold_pmids: Set[str],
) -> Dict:
    universe, tier1_pass, tier2_pass, stats = parse_filter_progress(progress_path)
    gold_n = len(gold_pmids)

    universe_in_gold = universe & gold_pmids
    tier1_in_gold = tier1_pass & gold_pmids
    tier2_in_gold = tier2_pass & gold_pmids
    missing_from_search = gold_pmids - universe
    dropped_by_tier1 = (universe & gold_pmids) - tier1_pass
    dropped_by_tier2 = (tier1_pass & gold_pmids) - tier2_pass

    def pct(num: int) -> float:
        return round(100.0 * num / gold_n, 2) if gold_n else 0.0

    return {
        "gene": gene,
        "filter_progress_path": str(progress_path),
        "gold_pmids": gold_n,
        "search_universe": len(universe),
        "tier1_pass": len(tier1_pass),
        "tier2_pass_for_download": len(tier2_pass),
        "search_coverage": {
            "matched": len(universe_in_gold),
            "missing": len(missing_from_search),
            "recall_pct": pct(len(universe_in_gold)),
            "sample_missing": sorted(missing_from_search)[:20],
        },
        "tier1_recall": {
            "matched": len(tier1_in_gold),
            "dropped_by_tier1": len(dropped_by_tier1),
            "recall_pct": pct(len(tier1_in_gold)),
            "sample_dropped": sorted(dropped_by_tier1)[:20],
        },
        "filter_recall": {
            "matched": len(tier2_in_gold),
            "dropped_by_tier2": len(dropped_by_tier2),
            "recall_pct": pct(len(tier2_in_gold)),
            "sample_dropped": sorted(dropped_by_tier2)[:20],
        },
        "downstream_load": {
            "universe_size": len(universe),
            "queued_for_download": len(tier2_pass),
            "reduction_pct": (
                round(100.0 * (1 - len(tier2_pass) / len(universe)), 2)
                if universe
                else 0.0
            ),
        },
        "parser_stats": stats,
    }


def render_markdown(rows: List[Dict], aggregate: Dict, generated_at: str) -> str:
    lines = []
    lines.append("# PMID-Filter Recall — Tier-1+Tier-2 Screening")
    lines.append("")
    lines.append(f"_Generated {generated_at}_")
    lines.append("")
    lines.append(
        "Compares the set of PMIDs that the tier-1 keyword filter + tier-2 LLM "
        "filter would queue for full-text download against the gold-standard "
        "PMID list. Search coverage caps everything downstream — no filter can "
        "recover a PMID PubMed/PubMind never returned."
    )
    lines.append("")
    lines.append(
        "| Gene | Gold PMIDs | Search universe | Search coverage | After Tier-1 | After Tier-2 (queued) | Final recall |"
    )
    lines.append(
        "|------|-----------:|----------------:|----------------:|-------------:|----------------------:|-------------:|"
    )
    for r in rows:
        lines.append(
            f"| {r['gene']} | {r['gold_pmids']} | {r['search_universe']:,} | "
            f"{r['search_coverage']['recall_pct']}% ({r['search_coverage']['matched']}/{r['gold_pmids']}) | "
            f"{r['tier1_recall']['recall_pct']}% ({r['tier1_recall']['matched']}/{r['gold_pmids']}) | "
            f"{r['filter_recall']['matched']:,} | "
            f"**{r['filter_recall']['recall_pct']}%** ({r['filter_recall']['matched']}/{r['gold_pmids']}) |"
        )
    lines.append("")
    lines.append("## Aggregate (union of all gold sets)")
    lines.append("")
    a = aggregate
    lines.append(f"- Gold PMIDs (union, deduped): **{a['gold_pmids']}**")
    lines.append(
        f"- Reached by search: **{a['search_coverage']['matched']} ({a['search_coverage']['recall_pct']}%)**"
    )
    lines.append(
        f"- Survived Tier-1: **{a['tier1_recall']['matched']} ({a['tier1_recall']['recall_pct']}%)**"
    )
    lines.append(
        f"- Queued for download after Tier-2: **{a['filter_recall']['matched']} ({a['filter_recall']['recall_pct']}%)**"
    )
    lines.append("")
    lines.append("## Where the gold PMIDs are lost")
    lines.append("")
    lines.append(
        "| Gene | Missing from search | Dropped by Tier-1 | Dropped by Tier-2 |"
    )
    lines.append(
        "|------|--------------------:|------------------:|------------------:|"
    )
    for r in rows:
        lines.append(
            f"| {r['gene']} | {r['search_coverage']['missing']} | "
            f"{r['tier1_recall']['dropped_by_tier1']} | "
            f"{r['filter_recall']['dropped_by_tier2']} |"
        )
    lines.append("")
    lines.append("## Downstream load reduction")
    lines.append("")
    lines.append(
        "| Gene | PMIDs after search | PMIDs queued for download | Reduction |"
    )
    lines.append(
        "|------|-------------------:|--------------------------:|----------:|"
    )
    for r in rows:
        lines.append(
            f"| {r['gene']} | {r['search_universe']:,} | "
            f"{r['downstream_load']['queued_for_download']:,} | "
            f"{r['downstream_load']['reduction_pct']}% |"
        )
    lines.append("")
    lines.append("## Sources")
    lines.append("")
    for r in rows:
        lines.append(f"- **{r['gene']}**: `{r['filter_progress_path']}`")
    lines.append("")
    return "\n".join(lines)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--gene",
        action="append",
        metavar="GENE",
        help="Gene to evaluate (repeatable). Defaults to all four targets.",
    )
    p.add_argument(
        "--filter-progress",
        action="append",
        default=[],
        metavar="GENE=PATH",
        help=(
            "Explicit filter_progress.jsonl path per gene "
            "(e.g. KCNH2=validation_runs/turnkey.../filter_progress.jsonl). "
            "Repeatable. Genes without an override are auto-discovered."
        ),
    )
    p.add_argument(
        "--gold-dir",
        type=Path,
        default=GOLD_DIR,
        help=f"Directory of {{GENE}}_recall_input.csv (default: {GOLD_DIR})",
    )
    p.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output directory for per-gene JSON + summary report.",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    genes: Iterable[str] = args.gene if args.gene else GENES

    overrides: Dict[str, Path] = {}
    for spec in args.filter_progress:
        if "=" not in spec:
            print(
                f"warning: --filter-progress spec missing '=' : {spec!r}",
                file=sys.stderr,
            )
            continue
        g, _, p = spec.partition("=")
        overrides[g.strip().upper()] = Path(p.strip())

    args.out.mkdir(parents=True, exist_ok=True)
    per_gene: List[Dict] = []

    for gene in genes:
        gene = gene.upper()
        if gene in overrides:
            progress_path = overrides[gene]
        else:
            progress_path = find_default_filter_progress(gene)
        if not progress_path or not progress_path.exists():
            print(
                f"[{gene}] skipping — no filter_progress.jsonl found "
                "(pass --filter-progress {gene}=PATH or run `gvf discover` first)",
                file=sys.stderr,
            )
            continue
        gold = load_gold_pmids(gene, args.gold_dir)
        result = compute_per_gene(gene, progress_path, gold)
        per_gene.append(result)
        out_json = args.out / f"{gene}_pmid_recall.json"
        out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(
            f"[{gene}] gold={result['gold_pmids']} "
            f"universe={result['search_universe']} "
            f"coverage={result['search_coverage']['recall_pct']}% "
            f"tier1={result['tier1_recall']['recall_pct']}% "
            f"final={result['filter_recall']['recall_pct']}% "
            f"(queued={result['filter_recall']['matched']}/{result['gold_pmids']})"
        )

    # Aggregate over union of gold sets.
    union_gold: Set[str] = set()
    union_universe: Set[str] = set()
    union_tier1: Set[str] = set()
    union_tier2: Set[str] = set()
    for r in per_gene:
        gene = r["gene"]
        gold = load_gold_pmids(gene, args.gold_dir)
        union_gold |= gold
        u, t1, t2, _ = parse_filter_progress(Path(r["filter_progress_path"]))
        union_universe |= u
        union_tier1 |= t1
        union_tier2 |= t2

    def pct(num: int, denom: int) -> float:
        return round(100.0 * num / denom, 2) if denom else 0.0

    aggregate = {
        "gold_pmids": len(union_gold),
        "search_coverage": {
            "matched": len(union_universe & union_gold),
            "recall_pct": pct(len(union_universe & union_gold), len(union_gold)),
        },
        "tier1_recall": {
            "matched": len(union_tier1 & union_gold),
            "recall_pct": pct(len(union_tier1 & union_gold), len(union_gold)),
        },
        "filter_recall": {
            "matched": len(union_tier2 & union_gold),
            "recall_pct": pct(len(union_tier2 & union_gold), len(union_gold)),
        },
    }

    generated_at = datetime.now(tz=timezone.utc).isoformat(timespec="seconds")
    summary = {
        "generated_at": generated_at,
        "per_gene": per_gene,
        "aggregate": aggregate,
    }
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    md = render_markdown(per_gene, aggregate, generated_at)
    (args.out / "summary.md").write_text(md, encoding="utf-8")
    print(f"\nWrote {args.out / 'summary.md'} and summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
