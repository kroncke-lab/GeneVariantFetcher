#!/usr/bin/env python3
"""Score the pipeline against an assistant's completed curation packet.

Closes the cold-start measurement loop: takes the filled
``curation_template.csv`` from ``build_curation_packet.py`` and the gene's DB,
converts the human answers into a normalized gold recall input, and runs the
existing recall scorer (recall + precision + F2) restricted to the curated
PMIDs. The result is a real, defensible recall/F2 number for a gene that had no
gold standard.

Germline is the heritable-carrier target: ``somatic``-labelled rows are excluded
by default (they are not carriers); ``unknown``/blank are kept. ``NONE`` and
empty-variant rows are dropped from the gold variant list but their PMIDs still
count as curated (gold) PMIDs, so the pipeline IS penalized for false positives
there via the precision metric.

Usage:
    python scripts/score_curation_packet.py \
        --filled-csv curation_template_FILLED.csv \
        --db results/BRCA2/<...>/BRCA2.db --gene BRCA2 \
        --out recall_metrics/brca2_gold_50
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


def _norm(s: str | None) -> str:
    return (s or "").strip()


def convert(filled_csv: Path, include_somatic: bool) -> tuple[list[dict], dict]:
    """Filled curation CSV -> normalized recall rows + a coverage summary."""
    rows: list[dict] = []
    pmids: set[str] = set()
    dropped_somatic = none_rows = 0
    with filled_csv.open(newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            pmid = _norm(r.get("pmid"))
            if not pmid or pmid.upper() == "EXAMPLE":
                continue
            pmids.add(pmid)
            variant = _norm(r.get("variant"))
            if not variant or variant.upper() == "NONE":
                none_rows += 1
                continue
            status = _norm(r.get("germline_or_somatic")).lower()
            if status.startswith("som") and not include_somatic:
                dropped_somatic += 1
                continue
            rows.append(
                {
                    "variant": variant,
                    "pmid": pmid,
                    "carriers": _norm(r.get("carriers")),
                    "affected": _norm(r.get("affected")),
                    "unaffected": _norm(r.get("unaffected")),
                }
            )
    summary = {
        "curated_pmids": len(pmids),
        "gold_variant_rows": len(rows),
        "no_variant_papers": none_rows,
        "dropped_somatic": dropped_somatic,
        "include_somatic": include_somatic,
    }
    return rows, summary


def write_recall_input(rows: list[dict], gene: str, gold_dir: Path) -> Path:
    norm = gold_dir / "normalized"
    norm.mkdir(parents=True, exist_ok=True)
    path = norm / f"{gene}_recall_input.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f, fieldnames=["variant", "pmid", "carriers", "affected", "unaffected"]
        )
        w.writeheader()
        w.writerows(rows)
    return path


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--filled-csv", required=True, type=Path)
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument(
        "--include-somatic",
        action="store_true",
        help="Count somatic-labelled variants too (default: germline/unknown only).",
    )
    args = ap.parse_args()

    gene = args.gene.upper()
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    rows, summary = convert(args.filled_csv.expanduser(), args.include_somatic)
    if not rows:
        print(
            "error: no gold variant rows after conversion — is the CSV filled?",
            file=sys.stderr,
        )
        return 2
    gold_dir = out / "_gold"
    write_recall_input(rows, gene, gold_dir)
    (out / "curation_coverage.json").write_text(json.dumps(summary, indent=2))
    print("Curation coverage:", json.dumps(summary))

    cmd = [
        PY,
        str(REPO / "scripts" / "run_recall_suite.py"),
        "--score",
        "--genes",
        gene,
        "--db",
        f"{gene}={args.db.expanduser().resolve()}",
        "--gold-dir",
        str(gold_dir),
        "--outdir",
        str(out),
    ]
    print("Scoring:", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=str(REPO))
    if proc.returncode != 0:
        return proc.returncode

    summ_path = out / gene / "summary.json"
    if summ_path.exists():
        s = json.loads(summ_path.read_text())
        rec = s.get("recall", {})
        uv = rec.get("unique_variants", {})
        fb = (s.get("mae") and None) or s.get("aggregate_fbeta") or {}
        print("\n=== COLD-START RESULT (vs curated gold) ===")
        print(
            f"  unique-variant recall: {uv.get('matched')}/{uv.get('gold')} = {uv.get('recall')}"
        )
        prec = (s.get("precision") or {}).get("precision_vs_gold_pmids")
        print(f"  precision (vs gold PMIDs): {prec}")
        if fb:
            print(f"  F2 (conservative): {fb.get('fbeta_vs_gold_pmids')}")
        print(f"  full report: {summ_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
