#!/usr/bin/env python3
"""Score recall/MAE/RMSE/precision twice: on all gold, and on acquirable gold.

A single headline number conflates two very different failures. If a paper's
full text was never downloaded, no extractor could have found its variants —
that is an *acquisition* miss. If the text is on disk and complete and the
variants were still missed, that is an *extraction* miss. Only the second is a
prompt/model problem, and only the first is fixed by fetching more source.

So every metric is reported for two populations:

  ALL GOLD           every gold variant-paper row, nothing excluded. This is
                     the honest end-to-end number and the one to quote
                     externally.
  SOURCE-COMPLETE    the subset whose paper has real full text on disk (not the
                     abstract-only fallback) and no supplement that the paper's
                     own markup advertised but we never fetched. This is the
                     ceiling the current corpus can support.
  SOURCE-INCOMPLETE  the disjoint remainder. Reported so the two subsets add
                     back up to ALL GOLD and the contrast is visible rather
                     than inferred.

SOURCE-COMPLETE minus SOURCE-INCOMPLETE is the acquisition debt, quantified.

Reads the per-row ``discrepancies.csv`` that ``scripts/run_recall_suite.py``
already writes, so it adds no scoring logic of its own and cannot drift from
the canonical matcher.

  python scripts/recall_audit/source_stratified_metrics.py \
      --metrics-dir recall_metrics/current --corpus corpus
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from harvesting.supplement_link_resolver import (  # noqa: E402
    looks_like_supplement_file,
    parse_links_from_markdown,
)

# Written by orchestrator._build_abstract_only_fallback as the first line when
# every full-text route failed. Its presence means the body is a PubMed
# abstract, so any variant in a table or supplement was unreachable.
ABSTRACT_ONLY_MARKER = "# ABSTRACT-ONLY FALLBACK"

# A FULL_CONTEXT smaller than this is a stub, not an article body. Abstract-only
# fallbacks are ~2-4 KB; the corpus median is ~30 KB.
MIN_FULL_TEXT_CHARS = 5_000

COUNT_FIELDS = ("carriers", "affected", "unaffected")


@dataclass
class PaperSource:
    """What we actually hold on disk for one paper."""

    pmid: str
    gene: str
    exists: bool = False
    abstract_only: bool = False
    chars: int = 0
    linked_supplements: int = 0
    missing_supplements: int = 0

    @property
    def has_full_text(self) -> bool:
        return (
            self.exists and not self.abstract_only and self.chars >= MIN_FULL_TEXT_CHARS
        )

    @property
    def is_source_complete(self) -> bool:
        """Full text present AND nothing the paper advertised is still missing."""
        return self.has_full_text and self.missing_supplements == 0

    def reason(self) -> str:
        if not self.exists:
            return "no source on disk"
        if self.abstract_only:
            return "abstract-only fallback"
        if self.chars < MIN_FULL_TEXT_CHARS:
            return "source too short to be a body"
        if self.missing_supplements:
            return f"{self.missing_supplements} advertised supplement(s) not fetched"
        return "complete"


def classify_paper(corpus: Path, gene: str, pmid: str) -> PaperSource:
    out = PaperSource(pmid=pmid, gene=gene)
    paper_dir = corpus / gene / pmid
    full_context = paper_dir / f"{pmid}_FULL_CONTEXT.md"
    if not full_context.is_file():
        return out
    out.exists = True
    try:
        text = full_context.read_text(errors="replace")
    except OSError:
        return out
    out.chars = len(text)
    out.abstract_only = text.lstrip().startswith(ABSTRACT_ONLY_MARKER)

    links = [
        link
        for link in parse_links_from_markdown(text)
        if looks_like_supplement_file(link.href)
    ]
    out.linked_supplements = len(links)
    if links:
        on_disk = {
            p.name.lower()
            for d in paper_dir.glob(f"{pmid}_supplements")
            for p in d.rglob("*")
            if p.is_file()
        }
        out.missing_supplements = sum(
            1 for link in links if Path(link.href).name.lower() not in on_disk
        )
    return out


@dataclass
class Stratum:
    """Accumulated metrics for one population of gold rows."""

    label: str
    gold_rows: int = 0
    matched_rows: int = 0
    gold_pmids: set = field(default_factory=set)
    matched_pmids: set = field(default_factory=set)
    extras: int = 0  # DB-only rows on gold PMIDs (loose)
    counted_extras: int = 0  # ...of those, rows carrying at least one count
    abs_err: dict = field(default_factory=lambda: {f: 0.0 for f in COUNT_FIELDS})
    sq_err: dict = field(default_factory=lambda: {f: 0.0 for f in COUNT_FIELDS})
    n_err: dict = field(default_factory=lambda: {f: 0 for f in COUNT_FIELDS})

    def add_gold(self, pmid: str, matched: bool) -> None:
        self.gold_rows += 1
        self.gold_pmids.add(pmid)
        if matched:
            self.matched_rows += 1
            self.matched_pmids.add(pmid)

    def add_error(self, field_name: str, diff: float) -> None:
        self.abs_err[field_name] += abs(diff)
        self.sq_err[field_name] += diff * diff
        self.n_err[field_name] += 1

    def to_dict(self) -> dict:
        counts = {}
        for f in COUNT_FIELDS:
            n = self.n_err[f]
            counts[f] = {
                "n_matched": n,
                "mae": (self.abs_err[f] / n) if n else None,
                "rmse": math.sqrt(self.sq_err[f] / n) if n else None,
            }
        loose_denom = self.matched_rows + self.extras
        counted_denom = self.matched_rows + self.counted_extras
        return {
            "label": self.label,
            "gold_rows": self.gold_rows,
            "matched_rows": self.matched_rows,
            "row_recall": (self.matched_rows / self.gold_rows)
            if self.gold_rows
            else None,
            "gold_pmids": len(self.gold_pmids),
            "matched_pmids": len(self.matched_pmids),
            "pmid_recall": (len(self.matched_pmids) / len(self.gold_pmids))
            if self.gold_pmids
            else None,
            "extra_rows_on_gold_pmids": self.extras,
            "counted_extra_rows_on_gold_pmids": self.counted_extras,
            # Names and semantics mirror cli.compare_variants.compute_precision_summary:
            # a false-positive UPPER BOUND on curated papers, not clean precision.
            "loose_precision_vs_gold_pmids": (self.matched_rows / loose_denom)
            if loose_denom
            else None,
            "precision_vs_counted_gold_pmids": (self.matched_rows / counted_denom)
            if counted_denom
            else None,
            "counts": counts,
        }


def _num(value: str):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def score_gene(gene: str, rows_path: Path, corpus: Path) -> dict:
    """Split one gene's per-row comparison into the two strata.

    Two passes, because an extra DB row only counts against precision when it
    sits on a PMID the gold standard actually curated — and that set is not
    known until every gold row has been read.
    """
    strata = {
        "all_gold": Stratum("ALL GOLD"),
        "source_complete": Stratum("SOURCE-COMPLETE"),
        "source_incomplete": Stratum("SOURCE-INCOMPLETE"),
    }
    sources: dict[str, PaperSource] = {}

    def source_for(pmid: str) -> PaperSource:
        if pmid not in sources:
            sources[pmid] = classify_paper(corpus, gene, pmid)
        return sources[pmid]

    def is_true(row: dict, key: str) -> bool:
        return (row.get(key) or "").strip().lower() in ("true", "1", "yes")

    with rows_path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))

    gold_pmids = {
        (r.get("pmid") or "").strip()
        for r in rows
        if not is_true(r, "missing_in_excel") and (r.get("pmid") or "").strip()
    }

    for row in rows:
        pmid = (row.get("pmid") or "").strip()
        if not pmid:
            continue
        src = source_for(pmid)
        targets = [strata["all_gold"]]
        targets.append(
            strata["source_complete"]
            if src.is_source_complete
            else strata["source_incomplete"]
        )

        if is_true(row, "missing_in_excel"):
            # A DB row gold never curated cannot be adjudicated; only rows on
            # gold PMIDs are judgeable.
            if pmid not in gold_pmids:
                continue
            has_count = any(
                _num(row.get(k)) is not None
                for k in (
                    "sqlite_carriers_total",
                    "sqlite_affected",
                    "sqlite_unaffected",
                )
            )
            for s in targets:
                s.extras += 1
                if has_count:
                    s.counted_extras += 1
            continue

        matched = not is_true(row, "missing_in_sqlite") and (
            (row.get("match_type") or "").strip() != "none"
        )
        for s in targets:
            s.add_gold(pmid, matched)
        if not matched:
            continue

        for f in COUNT_FIELDS:
            gold = _num(
                row.get("excel_carriers_total" if f == "carriers" else f"excel_{f}")
            )
            got = _num(
                row.get("sqlite_carriers_total" if f == "carriers" else f"sqlite_{f}")
            )
            if gold is None or got is None:
                continue
            for s in targets:
                s.add_error(f, got - gold)

    acquisition = {
        "papers_seen": len(sources),
        "source_complete": sum(1 for s in sources.values() if s.is_source_complete),
        "no_source_on_disk": sum(1 for s in sources.values() if not s.exists),
        "abstract_only": sum(1 for s in sources.values() if s.abstract_only),
        "too_short": sum(
            1
            for s in sources.values()
            if s.exists and not s.abstract_only and s.chars < MIN_FULL_TEXT_CHARS
        ),
        "missing_supplements": sum(
            1 for s in sources.values() if s.has_full_text and s.missing_supplements
        ),
    }
    return {
        "gene": gene,
        "acquisition": acquisition,
        "strata": {k: v.to_dict() for k, v in strata.items()},
    }


def _pct(value) -> str:
    return "n/a" if value is None else f"{value * 100:.1f}%"


def _num_fmt(value) -> str:
    return "n/a" if value is None else f"{value:.3f}"


def render(results: list[dict]) -> str:
    lines: list[str] = []
    lines.append("# Source-stratified recall metrics\n")
    lines.append(
        "ALL GOLD is the honest end-to-end number. SOURCE-COMPLETE restricts to\n"
        "papers whose full text is on disk and whose advertised supplements were\n"
        "all fetched — the ceiling the current corpus can support. The gap\n"
        "between them is acquisition debt, not extraction error.\n"
    )

    lines.append("\n## Acquisition state of gold papers\n")
    lines.append(
        "| gene | gold papers | source-complete | no source | abstract-only "
        "| too short | missing suppl. |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for r in results:
        a = r["acquisition"]
        lines.append(
            f"| {r['gene']} | {a['papers_seen']} | {a['source_complete']} "
            f"| {a['no_source_on_disk']} | {a['abstract_only']} | {a['too_short']} "
            f"| {a['missing_supplements']} |"
        )

    for key, title in (
        ("all_gold", "ALL GOLD — every gold row"),
        ("source_complete", "SOURCE-COMPLETE — full text on disk, nothing missing"),
        ("source_incomplete", "SOURCE-INCOMPLETE — the acquisition-limited remainder"),
    ):
        lines.append(f"\n## {title}\n")
        lines.append(
            "| gene | rows | row recall | pmids | pmid recall | precision "
            "| carriers MAE / RMSE | affected MAE / RMSE | unaffected MAE / RMSE |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in results:
            s = r["strata"][key]
            c = s["counts"]
            lines.append(
                f"| {r['gene']} | {s['matched_rows']}/{s['gold_rows']} "
                f"| {_pct(s['row_recall'])} "
                f"| {s['matched_pmids']}/{s['gold_pmids']} "
                f"| {_pct(s['pmid_recall'])} "
                f"| {_pct(s['precision_vs_counted_gold_pmids'])} "
                f"| {_num_fmt(c['carriers']['mae'])} / {_num_fmt(c['carriers']['rmse'])} "
                f"| {_num_fmt(c['affected']['mae'])} / {_num_fmt(c['affected']['rmse'])} "
                f"| {_num_fmt(c['unaffected']['mae'])} / "
                f"{_num_fmt(c['unaffected']['rmse'])} |"
            )
    return "\n".join(lines) + "\n"


def aggregate(results: list[dict]) -> dict:
    """Pool every gene into one all-gene row per stratum."""
    pooled = {
        "all_gold": Stratum("ALL GOLD"),
        "source_complete": Stratum("SOURCE-COMPLETE"),
        "source_incomplete": Stratum("SOURCE-INCOMPLETE"),
    }
    for r in results:
        for key, s in r["strata"].items():
            p = pooled[key]
            p.gold_rows += s["gold_rows"]
            p.matched_rows += s["matched_rows"]
            p.extras += s["extra_rows_on_gold_pmids"]
            p.counted_extras += s["counted_extra_rows_on_gold_pmids"]
            # PMIDs are gene-scoped, so pooling by count is exact here.
            p.gold_pmids.update(f"{r['gene']}:{i}" for i in range(s["gold_pmids"]))
            p.matched_pmids.update(
                f"{r['gene']}:m{i}" for i in range(s["matched_pmids"])
            )
            for f in COUNT_FIELDS:
                c = s["counts"][f]
                n = c["n_matched"]
                if not n:
                    continue
                p.n_err[f] += n
                p.abs_err[f] += (c["mae"] or 0.0) * n
                p.sq_err[f] += ((c["rmse"] or 0.0) ** 2) * n
    return {k: v.to_dict() for k, v in pooled.items()}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--metrics-dir",
        type=Path,
        required=True,
        help="run_recall_suite.py output dir holding <GENE>/discrepancies.csv",
    )
    ap.add_argument("--corpus", type=Path, default=REPO / "corpus")
    ap.add_argument("--genes", default="", help="Comma-separated subset.")
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args()

    wanted = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
    results = []
    for gene_dir in sorted(args.metrics_dir.iterdir()):
        if not gene_dir.is_dir():
            continue
        gene = gene_dir.name.upper()
        if wanted and gene not in wanted:
            continue
        rows = gene_dir / "discrepancies.csv"
        if not rows.is_file():
            continue
        results.append(score_gene(gene, rows, args.corpus))

    if not results:
        print(f"No <GENE>/discrepancies.csv under {args.metrics_dir}")
        return 1

    payload = {"genes": results, "aggregate": aggregate(results)}
    md = render(results)

    agg = payload["aggregate"]
    md += "\n## All genes pooled\n\n"
    md += (
        "| stratum | rows | row recall | precision | carriers MAE / RMSE "
        "| affected MAE / RMSE | unaffected MAE / RMSE |\n"
    )
    md += "|---|---:|---:|---:|---:|---:|---:|\n"
    for key in ("all_gold", "source_complete", "source_incomplete"):
        s = agg[key]
        c = s["counts"]
        md += (
            f"| {s['label']} | {s['matched_rows']}/{s['gold_rows']} "
            f"| {_pct(s['row_recall'])} | {_pct(s['precision_vs_counted_gold_pmids'])} "
            f"| {_num_fmt(c['carriers']['mae'])} / {_num_fmt(c['carriers']['rmse'])} "
            f"| {_num_fmt(c['affected']['mae'])} / {_num_fmt(c['affected']['rmse'])} "
            f"| {_num_fmt(c['unaffected']['mae'])} / "
            f"{_num_fmt(c['unaffected']['rmse'])} |\n"
        )

    print(md)
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"Wrote {args.out_json}")
    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text(md, encoding="utf-8")
        print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
