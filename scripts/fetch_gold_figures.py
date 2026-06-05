#!/usr/bin/env python3
"""Triage figure captions and fetch only the data-bearing figures, à la carte.

Most figures (current traces, confocal, topology schematics) carry no extractable
variant/patient data and aren't worth downloading. We already have the figure
*captions* for PMC JATS-XML papers (in FULL_CONTEXT.md), so we can score each
caption and pull ONLY the high-value figures — pedigrees (patient counts /
segregation) and clinical panels (ECG / phenotype) — from PMC.

Dry-run (default) writes a triage worklist CSV. `--fetch` downloads the selected
figures from PMC into corpus/<GENE>/<PMID>/<PMID>_figures/ (no LLM; one HTTP GET
per figure, rate-limited). Scope defaults to gold-standard PMIDs.

  python scripts/fetch_gold_figures.py                 # triage only
  python scripts/fetch_gold_figures.py --fetch --email you@example.com
"""

from __future__ import annotations

import argparse
import csv
import re
import time
from collections import Counter
from pathlib import Path
from urllib.request import Request, urlopen

REPO = Path(__file__).resolve().parents[1]

# Caption triage (lightweight "data scout"): pedigrees first (patient counts /
# segregation live there), then clinical panels; skip functional/structural.
PEDIGREE = ("pedigree", "family tree", "kindred", "proband", "segregat", " generation")
CLINICAL = (
    "ecg",
    "electrocardiogram",
    "qtc",
    "holter",
    "phenotyp",
    "carrier",
    "affected",
    "unaffected",
    "clinical features",
    "patient",
)
FUNCTIONAL = (
    "current trace",
    "tail current",
    "boltzmann",
    "voltage clamp",
    "voltage-clamp",
    "conductance",
    "kinetic",
    "confocal",
    "immunoblot",
    "western",
    "topology",
    "schematic",
    "alignment",
    "construct",
    "fluoresc",
    "crystal",
    "simulation",
    "representative current",
    "patch clamp",
    "patch-clamp",
)

FIG_RE = re.compile(r"^#{2,4}\s*(Figure|Fig\.?)\s*([0-9]+[A-Za-z]?)\b", re.I)
IMG_RE = re.compile(r"^_image_:\s*(\S+)")


def parse_figures(text: str) -> list[dict]:
    """Extract (label, caption, image_filename) from FULL_CONTEXT figure sections."""
    lines = text.replace("\r\n", "\n").split("\n")
    figs, cur = [], None
    for ln in lines:
        m = FIG_RE.match(ln.strip())
        if m:
            if cur:
                figs.append(cur)
            cur = {"label": f"Figure {m.group(2)}", "caption": "", "img": ""}
            continue
        if cur is None:
            continue
        mi = IMG_RE.match(ln.strip())
        if mi:
            # strip query string / fragment / path; cap length (some are full CDN URLs)
            base = (
                mi.group(1).split("?")[0].split("#")[0].rstrip("/").rsplit("/", 1)[-1]
            )
            cur["img"] = base[:120]
        elif ln.strip().startswith("#"):
            figs.append(cur)
            cur = None
        else:
            cur["caption"] += " " + ln.strip()
    if cur:
        figs.append(cur)
    for f in figs:
        f["caption"] = f["caption"].strip()
    return figs


def triage(caption: str) -> tuple[str, str]:
    low = caption.lower()
    if any(k in low for k in PEDIGREE):
        return "fetch", "pedigree (patient counts/segregation)"
    if any(k in low for k in CLINICAL):
        return "fetch", "clinical/phenotype panel"
    if any(k in low for k in FUNCTIONAL):
        return "skip", "functional/structural"
    return "skip", "other (not data-bearing)"


def gold_pmids(gene: str, gold_dir: Path) -> set[str]:
    p = gold_dir / "normalized" / f"{gene}_recall_input.csv"
    return (
        {r["pmid"].strip() for r in csv.DictReader(p.open(newline="")) if r.get("pmid")}
        if p.exists()
        else set()
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--genes", default="", help="Comma-separated genes (default: all gold genes)."
    )
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument(
        "--gold-dir", default=str(REPO / "gene_variant_fetcher_gold_standard")
    )
    ap.add_argument("--out", default=str(REPO / "corpus" / "figure_triage.csv"))
    ap.add_argument(
        "--fetch", action="store_true", help="Download the selected figures from PMC."
    )
    ap.add_argument(
        "--email", default="", help="Contact email for NCBI (recommended with --fetch)."
    )
    ap.add_argument(
        "--limit", type=int, default=0, help="Max figures to fetch (0 = all selected)."
    )
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

    rows, decisions = [], Counter()
    for gene in genes:
        for pmid in sorted(gold_pmids(gene, gold_dir)):
            d = corpus / gene / pmid
            ft = d / f"{pmid}_FULL_CONTEXT.md"
            figdir = d / f"{pmid}_figures"
            if not ft.exists():
                continue
            import json

            art = d / f"{pmid}_artifacts.json"
            pmcid = ""
            if art.exists():
                try:
                    pmcid = json.loads(art.read_text()).get("pmcid") or ""
                except (json.JSONDecodeError, OSError):
                    pass
            for f in parse_figures(ft.read_text(encoding="utf-8", errors="replace")):
                if not f["img"]:
                    continue
                already = (figdir / f["img"]).exists()
                decision, reason = triage(f["caption"])
                if already:
                    decision, reason = "have", "already downloaded"
                decisions[decision] += 1
                rows.append(
                    {
                        "gene": gene,
                        "pmid": pmid,
                        "pmcid": pmcid,
                        "figure": f["label"],
                        "image": f["img"],
                        "decision": decision,
                        "reason": reason,
                        "article_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
                        if pmcid
                        else "",
                        "caption": f["caption"][:160],
                    }
                )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["gene"])
        w.writeheader()
        w.writerows(rows)
    print(f"figures triaged: {len(rows)} -> {out}")
    print("decision:", dict(decisions.most_common()))
    to_fetch = [r for r in rows if r["decision"] == "fetch" and r["pmcid"]]
    print(
        f"selected for à-la-carte fetch: {len(to_fetch)} (vs {len(rows)} total figures)"
    )

    if args.fetch and to_fetch:
        ua = f"GeneVariantFetcher/figure-triage ({args.email or 'no-email'})"
        targets = to_fetch[: args.limit] if args.limit else to_fetch
        cdn_cache: dict[str, dict[str, str]] = {}
        ok = 0
        for r in targets:
            dst = corpus / r["gene"] / r["pmid"] / f"{r['pmid']}_figures" / r["image"]
            if dst.exists():
                continue
            # Resolve real CDN image URLs by scraping the PMC article HTML once
            # per paper (images live at cdn.ncbi.nlm.nih.gov/pmc/blobs/.../<file>
            # with a hashed path that can't be constructed blindly).
            if r["pmcid"] not in cdn_cache:
                cdn_cache[r["pmcid"]] = {}
                try:
                    req = Request(r["article_url"], headers={"User-Agent": ua})
                    page = urlopen(req, timeout=30).read().decode("utf-8", "replace")  # noqa: S310
                    for u in re.findall(
                        r"https://cdn\.ncbi\.nlm\.nih\.gov/pmc/blobs/[^\"' ]+", page
                    ):
                        cdn_cache[r["pmcid"]][u.rsplit("/", 1)[-1]] = u
                    time.sleep(0.4)
                except Exception as e:  # noqa: BLE001
                    print(f"  html fail {r['pmcid']}: {e}")
            url = cdn_cache[r["pmcid"]].get(r["image"])
            if not url:
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            try:
                data = urlopen(
                    Request(url, headers={"User-Agent": ua}), timeout=30
                ).read()  # noqa: S310
                if data[:3] in (b"\xff\xd8\xff", b"\x89PN", b"GIF"):
                    dst.write_bytes(data)
                    ok += 1
                time.sleep(0.4)
            except Exception as e:  # noqa: BLE001
                print(f"  img fail {r['pmid']} {r['image']}: {e}")
        print(f"fetched {ok}/{len(targets)} figure images into corpus")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
