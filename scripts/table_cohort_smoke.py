#!/usr/bin/env python3
"""Gold-free transfer smoke for the table-cohort phenotype projection.

Runs ``pipeline.table_cohort_phenotype`` over every archived per-PMID extraction
JSON it can find under the given roots (production ``gvf-run`` output trees),
using each record's own frozen source text, and lists every table the module
would classify as a case series or a control cohort -- caption, count column,
tier, rows -- plus a tally of refusal reasons. No LLM call, no database write.

The point is to audit the classifier on genes and paper shapes that have no
curated answer key (BRCA1/2, MYBPC3, BMPR2, APOE): every listed table must be
defensible from its caption and header alone.

Example::

    scripts/table_cohort_smoke.py --root results --genes BRCA1 BRCA2 MYBPC3 \
        --out docs/evidence/table_cohort_phenotype_20260907/transfer_smoke.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Optional

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
os.environ.setdefault("GVF_DISABLE_LOCAL_DATA", "1")

from pipeline.table_cohort_phenotype import (  # noqa: E402
    derive_table_cohort_phenotype_counts,
)


def source_text_for(extraction: dict[str, Any], path: Path, pmid: str) -> str:
    metadata = extraction.get("extraction_metadata") or {}
    candidates = [metadata.get("source_file")]
    run_dir = path.parent.parent
    for name in (f"{pmid}_FULL_CONTEXT.md", f"{pmid}_CLEANED.md"):
        candidates.append(run_dir / "pmc_fulltext" / name)
    for candidate in candidates:
        if not candidate:
            continue
        candidate = Path(candidate)
        if candidate.is_file():
            return candidate.read_text(errors="replace")
    return ""


def iter_extractions(roots: list[Path], genes: Optional[set[str]]):
    seen: set[str] = set()
    for root in roots:
        for path in sorted(root.rglob("*_PMID_*.json")):
            if "extractions" not in path.parts:
                continue
            gene, _, rest = path.name.partition("_PMID_")
            if genes and gene not in genes:
                continue
            pmid = rest[:-5]
            key = f"{gene}:{pmid}:{path.stat().st_size}"
            if key in seen:
                continue
            seen.add(key)
            yield gene, pmid, path


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, action="append", required=True)
    parser.add_argument("--genes", nargs="*")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--paper-tier", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args(argv)

    genes = set(args.genes) if args.genes else None
    tables: list[dict[str, Any]] = []
    reasons: Counter = Counter()
    papers = 0
    rows_applied = 0
    for gene, pmid, path in iter_extractions(args.root, genes):
        if args.limit and papers >= args.limit:
            break
        try:
            extraction = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if not isinstance(extraction, dict) or not extraction.get("variants"):
            continue
        papers += 1
        text = source_text_for(extraction, path, pmid)
        result = derive_table_cohort_phenotype_counts(
            extraction,
            text,
            gene_symbol=gene,
            disease=None,
            title=(extraction.get("paper_metadata") or {}).get("title"),
            enabled=True,
            allow_paper_tier=bool(args.paper_tier),
        )
        meta = result["extraction_metadata"].get(
            "table_cohort_phenotype_derivation", {}
        )
        for outcome in meta.get("outcomes", []):
            reasons[str(outcome.get("status"))] += 1
        for caption, info in (meta.get("tables") or {}).items():
            if info.get("role") is None:
                continue
            rows_applied += int(info.get("rows_applied") or 0) + int(
                info.get("rows_stamped") or 0
            )
            tables.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "run": str(path.parent.parent.relative_to(REPO))
                    if path.is_relative_to(REPO)
                    else str(path.parent.parent),
                    "role": info["role"],
                    "tier": info["tier"],
                    "caption": info["caption"],
                    "count_label": info["count_label"],
                    "count_unit": info["count_unit"],
                    "rows_seen": info["rows_seen"],
                    "rows_applied": info["rows_applied"],
                    "rows_stamped": info["rows_stamped"],
                    "quote": info["quote"],
                }
            )
    summary = {
        "roots": [str(r) for r in args.root],
        "genes": sorted(genes) if genes else None,
        "paper_tier": bool(args.paper_tier),
        "papers_scanned": papers,
        "tables_classified": len(tables),
        "rows_touched": rows_applied,
        "refusal_reasons": dict(reasons.most_common()),
        "tables": tables,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"papers={papers} tables_classified={len(tables)} rows_touched={rows_applied}"
    )
    for table in tables:
        print(
            f"  {table['gene']} {table['pmid']} [{table['role']}/{table['tier']}] "
            f"rows={table['rows_applied']}+{table['rows_stamped']} "
            f"col={table['count_label']!r} caption={table['caption'][:110]!r}"
        )
    print("refusals:", dict(reasons.most_common(12)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
