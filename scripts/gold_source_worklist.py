#!/usr/bin/env python3
"""Gold-PMID figure/supplement acquisition worklist + diagnosis.

For each gold-standard PMID, reports whether we have full text / figure images /
supplements, WHY any are missing, and the recommended fix — written to a CSV you
can act on. Reads only local artifacts (corpus + per-paper {PMID}_artifacts.json);
no network, no LLM.

Key finding it surfaces: papers fetched via the PMC **JATS-XML** route get figure
*captions* parsed but the figure *images are never downloaded* (images only come
from the PMC-HTML / browser routes). Those papers all have a PMCID, so the images
are fetchable from PMC — the single biggest, easiest figure-coverage win.

Usage:
  python scripts/gold_source_worklist.py                       # all gold genes
  python scripts/gold_source_worklist.py --genes KCNH2,SCN5A --out wl.csv
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIELDS = [
    "gene",
    "pmid",
    "pmcid",
    "doi",
    "source_route",
    "has_fulltext",
    "fig_captions",
    "fig_images",
    "figure_status",
    "supp_advertised",
    "supp_converted",
    "supplement_status",
    "recommended_action",
]


def gold_pmids(gene: str, gold_dir: Path) -> set[str]:
    p = gold_dir / "normalized" / f"{gene}_recall_input.csv"
    out: set[str] = set()
    if p.exists():
        for row in csv.DictReader(p.open(newline="")):
            if row.get("pmid"):
                out.add(row["pmid"].strip())
    return out


def classify(gene: str, pmid: str, corpus: Path) -> dict:
    art = corpus / gene / pmid / f"{pmid}_artifacts.json"
    ft = corpus / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
    has_ft = ft.exists() and ft.stat().st_size > 500
    d = {}
    if art.exists():
        try:
            d = json.loads(art.read_text())
        except (json.JSONDecodeError, OSError):
            d = {}
    mt, summ = d.get("main_text") or {}, d.get("summary") or {}
    src = mt.get("source") or ("(no artifacts)" if not d else "?")
    pmcid = d.get("pmcid") or ""
    fcap = int(mt.get("figure_captions_count") or 0)
    fimg = int(summ.get("figure_count") or 0)
    sdesc = int(mt.get("supplement_descriptions_count") or 0)
    sconv = int(summ.get("supplements_converted") or 0)

    if not has_ft and not d:
        fig_status, sup_status = "no_source", "no_source"
    else:
        if fimg > 0:
            fig_status = "downloaded"
        elif fcap > 0 and pmcid:
            fig_status = "in_pmc_not_downloaded"  # FIXABLE: fetch images from PMC
        elif fcap > 0:
            fig_status = "captioned_no_pmcid"
        else:
            fig_status = "none_advertised"
        if sconv > 0:
            sup_status = "converted"
        elif sdesc > 0:
            sup_status = "advertised_not_converted"  # FIXABLE: re-convert/refetch
        else:
            sup_status = "none_advertised"

    actions = []
    if fig_status == "in_pmc_not_downloaded":
        actions.append("fetch figure images from PMC (have PMCID)")
    elif fig_status == "captioned_no_pmcid":
        actions.append("fetch figures via publisher/browser route")
    if sup_status == "advertised_not_converted":
        actions.append("re-download/convert supplement (textutil/antiword fallback)")
    if fig_status == "no_source":
        actions.append("acquire full text first (not in corpus)")
    return {
        "gene": gene,
        "pmid": pmid,
        "pmcid": pmcid,
        "doi": d.get("doi") or "",
        "source_route": src,
        "has_fulltext": int(has_ft),
        "fig_captions": fcap,
        "fig_images": fimg,
        "figure_status": fig_status,
        "supp_advertised": sdesc,
        "supp_converted": sconv,
        "supplement_status": sup_status,
        "recommended_action": "; ".join(actions),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--genes",
        default="",
        help="Comma-separated genes (default: all with a gold set).",
    )
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument(
        "--gold-dir", default=str(REPO / "gene_variant_fetcher_gold_standard")
    )
    ap.add_argument("--out", default=str(REPO / "corpus" / "gold_source_worklist.csv"))
    args = ap.parse_args()
    corpus, gold_dir = Path(args.corpus), Path(args.gold_dir)
    genes = (
        [g.strip().upper() for g in args.genes.split(",") if g.strip()]
        if args.genes
        else sorted(
            p.name.split("_recall_input")[0]
            for p in (gold_dir / "normalized").glob("*_recall_input.csv")
        )
    )

    rows = [
        classify(g, pmid, corpus)
        for g in genes
        for pmid in sorted(gold_pmids(g, gold_dir))
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    fig = Counter(r["figure_status"] for r in rows)
    sup = Counter(r["supplement_status"] for r in rows)
    print(f"gold PMIDs: {len(rows)}  -> {out}")
    print("figure_status:", dict(fig.most_common()))
    print("supplement_status:", dict(sup.most_common()))
    print(
        f"\nTOP FIX: {fig.get('in_pmc_not_downloaded', 0)} gold papers have figures in PMC "
        "(captions parsed, PMCID present) that were never downloaded — fetchable from PMC."
    )
    print(
        f"         {sup.get('advertised_not_converted', 0)} gold papers advertise supplements "
        "that were not converted — re-run supplement conversion."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
