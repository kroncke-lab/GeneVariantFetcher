#!/usr/bin/env python3
"""Build an adjudication sheet: the pipeline's extractions + provenance, for review.

Unlike the BLINDED curation packet (``build_curation_packet.py``), this sheet
SHOWS the pipeline's extracted values and the provenance behind each one so a
reviewer can verify *why* each extraction worked or didn't — and mark it
correct / wrong / not-in-paper, plus add any the pipeline MISSED.

Per extracted variant it surfaces: the cDNA/protein notation, carrier/affected/
unaffected counts, the count_provenance (which column/label/count-type each
number came from), the source_location, a verbatim key_quote, the model's own
notes, and the gold-free germline/somatic flag. Papers where the pipeline found
nothing get a "(none extracted)" row so the reviewer can confirm or correct.

Tradeoff: showing the pipeline's answers anchors the reviewer, so the recall it
yields is biased high — it measures precision + error analysis ("why"), not an
unbiased recall. For an unbiased recall number, use the blinded packet instead.

Emits a CSV (import into Google Sheets / Excel). Pair it with the curation
packet's papers/ folder so the reviewer can read each paper.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pipeline.somatic_germline_qc import classify_variant  # noqa: E402

PIPELINE_COLUMNS = [
    "pmid",
    "title",
    "pubmed_url",
    "fulltext_file",
    "pipeline_variant",
    "clinical_significance",
    "carriers",
    "affected",
    "unaffected",
    "germline_or_somatic",
    "source_location",
    "key_quote",
    "count_provenance",
    "pipeline_notes",
]
# Blank columns the reviewer fills in.
REVIEW_COLUMNS = [
    "verdict",  # correct | wrong_variant | wrong_count | not_in_paper | MISSED | unsure
    "correct_carriers",
    "correct_affected",
    "correct_unaffected",
    "reviewer_notes",
]


def _first(*vals):
    for v in vals:
        if v not in (None, ""):
            return v
    return ""


def _truncate(s, n=300):
    s = str(s or "").replace("\n", " ").strip()
    return s if len(s) <= n else s[: n - 1] + "…"


def _count_provenance(v: dict) -> str:
    cp = v.get("count_provenance") or {}
    bits = []
    for field in ("carriers", "affected", "unaffected"):
        label = cp.get(f"{field}_column_label")
        ctype = cp.get(f"{field}_count_type")
        if (label and label not in (None, "")) or (ctype and ctype != "unknown"):
            piece = f"{field}: {ctype or '?'}"
            if label:
                piece += f" ('{label}')"
            bits.append(piece)
    return "; ".join(bits)


def build_rows(records: list[dict]) -> list[dict]:
    """records: list of {pmid,title,url,fulltext_file,variants,abstract}. -> CSV rows."""
    rows: list[dict] = []
    for rec in records:
        base = {
            "pmid": rec["pmid"],
            "title": rec.get("title", ""),
            "pubmed_url": rec.get("url", ""),
            "fulltext_file": rec.get("fulltext_file", ""),
        }
        variants = rec.get("variants") or []
        if not variants:
            rows.append({**base, "pipeline_variant": "(none extracted)"})
            continue
        for v in variants:
            pen = v.get("penetrance_data") or {}
            pat = v.get("patients") or {}
            notation = " / ".join(
                x for x in (v.get("cdna_notation"), v.get("protein_notation")) if x
            )
            rows.append(
                {
                    **base,
                    "pipeline_variant": notation or "(no notation)",
                    "clinical_significance": v.get("clinical_significance") or "",
                    "carriers": _first(
                        pen.get("total_carriers_observed"), pat.get("count")
                    ),
                    "affected": _first(pen.get("affected_count")),
                    "unaffected": _first(pen.get("unaffected_count")),
                    "germline_or_somatic": classify_variant(
                        v, paper_context=rec.get("abstract", "")
                    ).label,
                    "source_location": _truncate(v.get("source_location"), 200),
                    "key_quote": _truncate((v.get("key_quotes") or [""])[0]),
                    "count_provenance": _count_provenance(v),
                    "pipeline_notes": _truncate(v.get("additional_notes")),
                }
            )
    return rows


def _load(run_dir: Path, manifest_csv: Path) -> list[dict]:
    ext_dir = run_dir / "extractions"
    abs_dir = run_dir / "abstract_json"
    records = []
    with manifest_csv.open(newline="", encoding="utf-8") as f:
        for m in csv.DictReader(f):
            pmid = m["pmid"]
            variants = []
            hits = list(ext_dir.glob(f"*{pmid}*.json"))
            if hits:
                try:
                    variants = (
                        json.loads(hits[0].read_text(encoding="utf-8")).get("variants")
                        or []
                    )
                except (OSError, json.JSONDecodeError):
                    pass
            abstract = ""
            aj = abs_dir / f"{pmid}.json"
            if aj.exists():
                try:
                    abstract = (
                        json.loads(aj.read_text(encoding="utf-8")).get("abstract") or ""
                    )
                except (OSError, json.JSONDecodeError):
                    pass
            records.append(
                {
                    "pmid": pmid,
                    "title": m.get("title", ""),
                    "url": m.get("pubmed_url", ""),
                    "fulltext_file": m.get("fulltext_file", ""),
                    "variants": variants,
                    "abstract": abstract,
                }
            )
    return records


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Packet manifest.csv (the sampled papers).",
    )
    ap.add_argument("--out", required=True, type=Path, help="Output adjudication CSV.")
    args = ap.parse_args()

    records = _load(args.run_dir.expanduser().resolve(), args.manifest.expanduser())
    rows = build_rows(records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=PIPELINE_COLUMNS + REVIEW_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    n_var = sum(1 for r in rows if r["pipeline_variant"] not in ("(none extracted)",))
    n_none = sum(1 for r in rows if r["pipeline_variant"] == "(none extracted)")
    print(f"Adjudication CSV: {args.out}")
    print(
        f"  {len(records)} papers, {n_var} extracted-variant rows, {n_none} no-variant papers"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
