#!/usr/bin/env python3
"""Rank papers for *manual* source acquisition by what a human download would buy.

Automated harvesting has a ceiling: some papers sit behind bot walls, some
publisher pages serve a landing stub, some PDFs render as glyph codes, and some
supplements were never linked. A person can fetch a few hundred of those by
hand, not thousands, so the question is *which* ones.

Yield is measured, not guessed, wherever gold exists: the source-presence
sweep (``scripts/recall_audit/gold_source_presence_sweep.py``) classifies every
gold row by whether its variant string is present in anything we hold on disk.
Rows in the four acquisition classes (``source_absent``,
``text_absent_stub_body``, ``text_absent_garbled_body``,
``text_absent_substitution``) are unreachable by any reading protocol until the
source is acquired; a paper's *hard yield* is the number of such rows, summed
over every gene that cites it, because one download serves all of them. The two
undecidable classes (figures present, non-searchable notation) form the wider
*possible yield* and are reported, never ranked on by default.

For a gene without gold (``--pmid-file`` mode) the same report is built from the
abstract-only acquisition expected-value predictor in ``scripts/acquisition_ev``.

Access is classified from what the corpus holds plus PubMed metadata and
Unpaywall open-access status, so the operator sees *why* a bot could not get the
paper and what to click:

* ``free_pdf_not_fetched``  -- Unpaywall lists an OA copy but nothing usable is on
  disk: a bot wall or a fetch failure; download the PDF at the listed URL.
* ``free_supplements``      -- OA body on disk but the gold variants are absent
  from every byte: fetch the supplements from the OA landing page or PMC.
* ``paywalled_pdf``         -- no OA copy and no usable body: institutional
  browser access for the PDF and its supplements.
* ``paywalled_supplements`` -- body on disk, variants absent, no OA copy.
* ``garbled_pdf``           -- the body is glyph codes; re-download a clean PDF.
* ``no_doi``                -- resolve the article by hand first.

Nothing is written to ``corpus/``; the output is a worklist for a person.

Example::

    scripts/recall_audit/rank_manual_acquisition.py \\
        --out-dir docs/evidence/manual_acquisition_20260908 --max-papers 200
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
os.environ.setdefault("GVF_DISABLE_LOCAL_DATA", "1")

HARD_CLASSES = (
    "source_absent",
    "text_absent_stub_body",
    "text_absent_garbled_body",
    "text_absent_substitution",
)
UNKNOWN_CLASSES = ("text_absent_figures_present", "text_absent_notation_inconclusive")
DEFAULT_GENES = ("KCNH2", "KCNQ1", "SCN5A", "RYR2", "BRCA2")
DEFAULT_SWEEP = REPO / "docs/evidence/gold_source_presence_sweep_20260903"
DEFAULT_EMAIL = "brett.kroncke@gmail.com"

WORKLIST_COLUMNS = [
    "rank",
    "pmid",
    "genes",
    "hard_yield",
    "positive_carrier_rows",
    "possible_yield",
    "gold_rows",
    "cumulative_hard_yield",
    "cumulative_hard_pct",
    "access_class",
    "manual_action",
    "corpus_state",
    "classes",
    "year",
    "journal",
    "publisher",
    "title",
    "doi",
    "doi_url",
    "pubmed_url",
    "pmc_id",
    "is_oa",
    "oa_status",
    "oa_pdf_url",
    "oa_host_type",
    "oa_landing_url",
    "inventory_status",
    "tranches",
    "body_chars",
    "figure_files",
    "supplement_files_on_disk",
    "unfetched_links",
    "predicted_ev",
]


def _pct(numerator: float, denominator: float) -> float:
    return round(100.0 * numerator / denominator, 1) if denominator else 0.0


# --------------------------------------------------------------------------- #
# Yield from gold (sweep) or from the abstract-only predictor
# --------------------------------------------------------------------------- #


def load_gold_yield(gold_rows: Path, genes: set[str]) -> dict[str, dict[str, Any]]:
    papers: dict[str, dict[str, Any]] = {}
    with gold_rows.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene"]
            if genes and gene not in genes:
                continue
            pmid = str(row["pmid"]).strip()
            paper = papers.setdefault(
                pmid,
                {
                    "pmid": pmid,
                    "genes": set(),
                    "gold_rows": 0,
                    "hard_yield": 0,
                    "positive_carrier_rows": 0,
                    "possible_yield": 0,
                    "classes": Counter(),
                    "inventory_status": set(),
                    "tranches": set(),
                    "body_chars": None,
                    "figure_files": None,
                },
            )
            paper["genes"].add(gene)
            paper["gold_rows"] += 1
            klass = row["class"]
            paper["classes"][klass] += 1
            paper["inventory_status"].add(row.get("inventory_status", ""))
            paper["tranches"].add(row.get("tranche", ""))
            if klass in HARD_CLASSES:
                paper["hard_yield"] += 1
                carriers = row.get("gold_carriers", "")
                if carriers not in ("", "0"):
                    paper["positive_carrier_rows"] += 1
            if klass in UNKNOWN_CLASSES:
                paper["possible_yield"] += 1
            for field in ("body_chars", "figure_files"):
                value = row.get(field)
                if value not in (None, ""):
                    try:
                        paper[field] = max(int(float(value)), paper[field] or 0)
                    except ValueError:
                        pass
    return papers


def load_predicted_yield(
    pmid_file: Path, gene: str, email: str
) -> dict[str, dict[str, Any]]:
    """Gold-free ranking: abstract-only expected value per paper."""
    from scripts.acquisition_ev import predict_yield

    pmids = [
        line.strip()
        for line in pmid_file.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    abstracts = predict_yield.fetch_abstracts(pmids, email)
    papers: dict[str, dict[str, Any]] = {}
    for pmid in pmids:
        record = abstracts.get(pmid, {})
        scored = predict_yield.score(predict_yield.compute_features(gene, record))
        papers[pmid] = {
            "pmid": pmid,
            "genes": {gene},
            "gold_rows": 0,
            "hard_yield": 0,
            "positive_carrier_rows": 0,
            "possible_yield": 0,
            "classes": Counter(),
            "inventory_status": set(),
            "tranches": set(),
            "body_chars": None,
            "figure_files": None,
            "predicted_ev": scored["ev_score"],
            "predicted_note": scored.get("note", ""),
            "title": record.get("title", ""),
            "journal": record.get("journal", ""),
            "year": record.get("year", ""),
        }
    return papers


# --------------------------------------------------------------------------- #
# Corpus state, PubMed metadata, Unpaywall
# --------------------------------------------------------------------------- #


def load_corpus_state(
    papers_json: Optional[Path], corpus_index: Optional[Path]
) -> tuple[dict[str, dict[str, Any]], dict[tuple[str, str], str]]:
    per_paper: dict[str, dict[str, Any]] = {}
    if papers_json and papers_json.is_file():
        per_paper = json.loads(papers_json.read_text())
    index_status: dict[tuple[str, str], str] = {}
    if corpus_index and corpus_index.is_file():
        with corpus_index.open(newline="") as handle:
            for row in csv.DictReader(handle):
                index_status[(row["gene"], str(row["pmid"]))] = row.get(
                    "full_text_status", ""
                )
    return per_paper, index_status


def corpus_summary(
    paper: dict[str, Any],
    per_paper: dict[str, dict[str, Any]],
    index_status: dict[tuple[str, str], str],
) -> dict[str, Any]:
    states: list[str] = []
    body_chars = 0
    figure_files = 0
    supplements = 0
    unfetched = 0
    garbled = False
    statuses: set[str] = set()
    for gene in sorted(paper["genes"]):
        entry = per_paper.get(f"{gene}:{paper['pmid']}") or {}
        if entry:
            if entry.get("body_state"):
                states.append(entry["body_state"])
            body_chars = max(body_chars, int(entry.get("body_chars") or 0))
            figure_files = max(figure_files, int(entry.get("figure_files") or 0))
            supplements = max(
                supplements,
                int(entry.get("supplement_converted_files") or 0)
                + int(entry.get("supplement_text_files") or 0)
                + int(entry.get("supplement_unsearchable_files") or 0)
                + int(entry.get("supplement_failed_files") or 0),
            )
            unfetched = max(unfetched, int(entry.get("unfetched_links") or 0))
            garbled = garbled or bool(entry.get("body_garbled"))
        status = index_status.get((gene, paper["pmid"]))
        if status:
            statuses.add(status)
    if not states:
        if body_chars:
            states.append("body on disk")
        elif paper.get("body_chars"):
            body_chars = int(paper["body_chars"])
            states.append("body on disk")
        else:
            states.append("no source on disk")
    return {
        "corpus_state": "; ".join(sorted(set(states)))
        + (f" [index: {', '.join(sorted(statuses))}]" if statuses else ""),
        "body_chars": body_chars,
        "figure_files": figure_files or int(paper.get("figure_files") or 0),
        "supplement_files_on_disk": supplements,
        "unfetched_links": unfetched,
        "garbled": garbled,
        "body_present": body_chars >= 6000,
    }


def load_cache(path: Path) -> dict[str, Any]:
    if path.is_file():
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def fetch_pubmed(pmids: list[str], cache: dict[str, Any], email: str) -> None:
    todo = [p for p in pmids if p not in cache.setdefault("pubmed", {})]
    if not todo:
        return
    from utils.pubmed_utils import batch_fetch_metadata

    records = batch_fetch_metadata(todo, email=email)
    for pmid in todo:
        record = records.get(pmid)
        if not record:
            cache["pubmed"][pmid] = {"missing": True}
            continue
        ids = record.get("ArticleIds") or {}
        doi = record.get("DOI") or ids.get("doi") or ""
        if isinstance(doi, list):
            doi = doi[0] if doi else ""
        pmc = ids.get("pmc") or ids.get("pmcid") or ""
        if isinstance(pmc, list):
            pmc = pmc[0] if pmc else ""
        pubdate = str(record.get("PubDate") or "")
        year_match = re.search(r"\b(19|20)\d{2}\b", pubdate)
        cache["pubmed"][pmid] = {
            "title": str(record.get("Title") or ""),
            "journal": str(record.get("Source") or ""),
            "full_journal": str(record.get("FullJournalName") or ""),
            "year": year_match.group(0) if year_match else pubdate,
            "doi": str(doi).strip(),
            "pmc_id": str(pmc).strip(),
        }


def fetch_unpaywall(dois: list[str], cache: dict[str, Any], email: str) -> None:
    todo = [d for d in dois if d and d not in cache.setdefault("unpaywall", {})]
    if not todo:
        return
    from harvesting.unpaywall_api import UnpaywallClient

    client = UnpaywallClient(email=email)
    for index, doi in enumerate(todo, 1):
        result, error = client.find_open_access(doi)
        if result is None:
            cache["unpaywall"][doi] = {"error": error or "lookup failed"}
        else:
            best = result.get("best_oa_location") or {}
            cache["unpaywall"][doi] = {
                "is_oa": bool(result.get("is_oa")),
                "oa_status": result.get("oa_status") or "",
                "publisher": result.get("publisher") or "",
                "journal_name": result.get("journal_name") or "",
                "pdf_url": result.get("pdf_url") or best.get("url_for_pdf") or "",
                "host_type": best.get("host_type") or "",
                "landing_url": best.get("url_for_landing_page")
                or best.get("url")
                or "",
            }
        if index % 25 == 0:
            print(f"  unpaywall {index}/{len(todo)}", file=sys.stderr)
        time.sleep(0.15)


# --------------------------------------------------------------------------- #
# Access class and the instruction to a human
# --------------------------------------------------------------------------- #


def classify_access(
    paper: dict[str, Any],
    corpus: dict[str, Any],
    meta: dict[str, Any],
    oa: dict[str, Any],
) -> tuple[str, str]:
    doi = meta.get("doi") or ""
    classes = paper["classes"]
    body_missing = (
        classes.get("source_absent", 0) + classes.get("text_absent_stub_body", 0) > 0
        or not corpus["body_present"]
    )
    variants_absent_from_body = classes.get("text_absent_substitution", 0) > 0
    if corpus.get("garbled") or classes.get("text_absent_garbled_body", 0) > 0:
        return (
            "garbled_pdf",
            "Re-download a clean PDF (current text is unmapped glyph codes), then re-run"
            " the corpus builder; supplements as well if the table is supplementary.",
        )
    if not doi:
        return (
            "no_doi",
            "Resolve the article by PMID/title first (PubMed record has no DOI), then"
            " fetch PDF and supplements by hand.",
        )
    is_oa = bool(oa.get("is_oa"))
    pdf = oa.get("pdf_url") or oa.get("landing_url") or f"https://doi.org/{doi}"
    if body_missing and is_oa:
        return (
            "free_pdf_not_fetched",
            f"Open-access copy exists ({oa.get('oa_status')}, {oa.get('host_type')})"
            f" but nothing usable is on disk: download the PDF at {pdf} and any"
            " supplements; a bot wall or fetch failure is the likely cause.",
        )
    if body_missing:
        return (
            "paywalled_pdf",
            f"No open-access copy: use institutional browser access at"
            f" https://doi.org/{doi} to download the PDF and its supplements.",
        )
    if variants_absent_from_body and is_oa:
        return (
            "free_supplements",
            f"Body is on disk but the gold variants are not in it: fetch the"
            f" supplementary tables from the open-access page ({pdf}) or PMC.",
        )
    if variants_absent_from_body:
        return (
            "paywalled_supplements",
            f"Body is on disk but the gold variants are not in it: download the"
            f" supplementary tables via institutional access at https://doi.org/{doi}.",
        )
    return (
        "review_needed",
        "Rows are behind the ceiling but the body is present and no supplement gap"
        " is evident; inspect the sweep classes before spending manual effort.",
    )


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def build_worklist(args: argparse.Namespace) -> dict[str, Any]:
    genes = set(args.genes or [])
    if args.pmid_file:
        if not args.gene:
            raise SystemExit("--pmid-file requires --gene")
        papers = load_predicted_yield(args.pmid_file, args.gene, args.email)
        rank_key = lambda p: (-float(p.get("predicted_ev") or 0.0), p["pmid"])  # noqa: E731
        yield_label = "predicted_ev"
    else:
        papers = load_gold_yield(args.gold_rows, genes)
        papers = {k: v for k, v in papers.items() if v["hard_yield"] > 0}
        rank_key = lambda p: (  # noqa: E731
            -p["hard_yield"],
            -p["possible_yield"],
            -p["gold_rows"],
            p["pmid"],
        )
        yield_label = "hard_yield"

    per_paper, index_status = load_corpus_state(args.papers_json, args.corpus_index)
    ranked = sorted(papers.values(), key=rank_key)
    selected = ranked[: args.max_papers] if args.max_papers else ranked

    cache_path = args.cache or (args.out_dir / "metadata_cache.json")
    cache = load_cache(cache_path)
    if not args.no_network:
        fetch_pubmed([p["pmid"] for p in selected], cache, args.email)
        dois = [(cache["pubmed"].get(p["pmid"]) or {}).get("doi", "") for p in selected]
        fetch_unpaywall([d for d in dois if d], cache, args.email)
        args.out_dir.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=1, sort_keys=True))

    total_hard = sum(p["hard_yield"] for p in ranked)
    rows: list[dict[str, Any]] = []
    cumulative = 0
    access_counter: Counter = Counter()
    for rank, paper in enumerate(selected, 1):
        meta = (cache.get("pubmed") or {}).get(paper["pmid"]) or {}
        oa = (cache.get("unpaywall") or {}).get(meta.get("doi") or "") or {}
        corpus = corpus_summary(paper, per_paper, index_status)
        access, action = classify_access(paper, corpus, meta, oa)
        access_counter[access] += 1
        cumulative += paper["hard_yield"]
        doi = meta.get("doi") or ""
        rows.append(
            {
                "rank": rank,
                "pmid": paper["pmid"],
                "genes": "+".join(sorted(paper["genes"])),
                "hard_yield": paper["hard_yield"],
                "positive_carrier_rows": paper["positive_carrier_rows"],
                "possible_yield": paper["possible_yield"],
                "gold_rows": paper["gold_rows"],
                "cumulative_hard_yield": cumulative,
                "cumulative_hard_pct": _pct(cumulative, total_hard),
                "access_class": access,
                "manual_action": action,
                "corpus_state": corpus["corpus_state"],
                "classes": "; ".join(
                    f"{k}={v}" for k, v in sorted(paper["classes"].items())
                ),
                "year": meta.get("year") or paper.get("year", ""),
                "journal": meta.get("journal") or paper.get("journal", ""),
                "publisher": oa.get("publisher", ""),
                "title": meta.get("title") or paper.get("title", ""),
                "doi": doi,
                "doi_url": f"https://doi.org/{doi}" if doi else "",
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/",
                "pmc_id": meta.get("pmc_id", ""),
                "is_oa": oa.get("is_oa", ""),
                "oa_status": oa.get("oa_status", ""),
                "oa_pdf_url": oa.get("pdf_url", ""),
                "oa_host_type": oa.get("host_type", ""),
                "oa_landing_url": oa.get("landing_url", ""),
                "inventory_status": "; ".join(
                    sorted(s for s in paper["inventory_status"] if s)
                ),
                "tranches": "; ".join(sorted(t for t in paper["tranches"] if t)),
                "body_chars": corpus["body_chars"],
                "figure_files": corpus["figure_files"],
                "supplement_files_on_disk": corpus["supplement_files_on_disk"],
                "unfetched_links": corpus["unfetched_links"],
                "predicted_ev": paper.get("predicted_ev", ""),
            }
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    worklist = args.out_dir / "manual_acquisition_worklist.csv"
    with worklist.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=WORKLIST_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    capture = {
        str(n): {
            "papers": min(n, len(ranked)),
            "hard_rows": sum(p["hard_yield"] for p in ranked[:n]),
            "pct_of_hard_rows": _pct(
                sum(p["hard_yield"] for p in ranked[:n]), total_hard
            ),
        }
        for n in (25, 50, 100, 150, 200, 250, 300)
    }
    summary = {
        "generated_from": {
            "gold_rows": str(args.gold_rows) if not args.pmid_file else None,
            "papers_json": str(args.papers_json) if args.papers_json else None,
            "corpus_index": str(args.corpus_index) if args.corpus_index else None,
            "pmid_file": str(args.pmid_file) if args.pmid_file else None,
            "genes": sorted(genes) if genes else None,
            "network": not args.no_network,
        },
        "yield_label": yield_label,
        "papers_with_hard_yield": len(ranked),
        "hard_rows_total": total_hard,
        "selected_papers": len(rows),
        "selected_hard_rows": cumulative,
        "selected_hard_pct": _pct(cumulative, total_hard),
        "capture_curve": capture,
        "access_classes_in_selection": dict(access_counter.most_common()),
        "hard_classes": HARD_CLASSES,
        "unknown_classes": UNKNOWN_CLASSES,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        f"# Manual acquisition worklist (top {len(rows)} of {len(ranked)} papers)",
        "",
        f"Ranked by gold rows behind the acquisition ceiling (`{yield_label}`), summed",
        "over every gene citing the paper. Cumulative % is of all",
        f"{total_hard} hard-ceiling rows across the ranked papers.",
        "",
        "| # | PMID | genes | hard | possible | cum % | access | year | journal | title |",
        "| ---: | --- | --- | ---: | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['rank']} | [{row['pmid']}]({row['pubmed_url']}) | {row['genes']} |"
            f" {row['hard_yield']} | {row['possible_yield']} | {row['cumulative_hard_pct']} |"
            f" {row['access_class']} | {row['year']} | {row['journal']} |"
            f" {str(row['title'])[:90]} |"
        )
    (args.out_dir / "manual_acquisition_worklist.md").write_text(
        "\n".join(lines) + "\n"
    )
    return summary


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--gold-rows", type=Path, default=DEFAULT_SWEEP / "gold_rows.tsv"
    )
    parser.add_argument(
        "--papers-json", type=Path, default=DEFAULT_SWEEP / "papers.json"
    )
    parser.add_argument(
        "--corpus-index", type=Path, default=REPO / "corpus" / "INDEX.csv"
    )
    parser.add_argument("--genes", nargs="*", default=list(DEFAULT_GENES))
    parser.add_argument(
        "--pmid-file",
        type=Path,
        help="gold-free mode: rank these PMIDs by predicted yield",
    )
    parser.add_argument("--gene", help="gene symbol for --pmid-file mode")
    parser.add_argument("--max-papers", type=int, default=200)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--email", default=DEFAULT_EMAIL)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--no-network", action="store_true")
    args = parser.parse_args(argv)
    summary = build_worklist(args)
    print(
        json.dumps({k: v for k, v in summary.items() if k != "capture_curve"}, indent=2)
    )
    print("capture:", json.dumps(summary["capture_curve"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
