#!/usr/bin/env python3
"""CLI: extract gene variants from figure images using a vision LLM.

Walks ``{pmc_fulltext_dir}/{pmid}_figures/`` for each requested PMID, asks a
multimodal model to read every image for variants of the given gene, and
writes one JSON report per PMID. Optionally appends the discovered variants
to a SQLite DB so they show up in recall comparisons.

Typical use::

    # Read figures for one PMID, just print JSON
    python scripts/extract_figure_variants.py \\
        --gene KCNH2 --pmid 24667783 \\
        --pmc-dir results/KCNH2/20260517_074737/pmc_fulltext \\
        --out /tmp/figure_reads

    # Process the top blockers and inject into the recall DB
    python scripts/extract_figure_variants.py \\
        --gene KCNH2 \\
        --pmid 29650123 --pmid 24667783 --pmid 19038855 \\
        --pmc-dir results/KCNH2/20260517_074737/pmc_fulltext \\
        --out recall_metrics/figure_reads_$(date +%Y%m%d) \\
        --db results/KCNH2/20260517_074737/KCNH2.db
"""

from __future__ import annotations

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path
from typing import Iterable, List

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from harvesting.figure_variant_reader import (  # noqa: E402
    PMIDFigureReport,
    find_pmid_figures,
    read_figures_for_pmid,
)


logger = logging.getLogger("extract_figure_variants")


def _load_pmid_list(args: argparse.Namespace) -> List[str]:
    pmids: List[str] = list(args.pmid or [])
    if args.pmid_file:
        for line in Path(args.pmid_file).read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                pmids.append(line.split()[0])
    if args.missing_csv:
        import csv

        with Path(args.missing_csv).open() as f:
            seen = set()
            for row in csv.DictReader(f):
                pid = (row.get("pmid") or "").strip()
                if pid and pid not in seen:
                    seen.add(pid)
                    pmids.append(pid)
    if not pmids:
        sys.exit("No PMIDs provided. Use --pmid, --pmid-file, or --missing-csv.")
    seen = set()
    out: List[str] = []
    for p in pmids:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _ensure_paper(con: sqlite3.Connection, pmid: str, gene: str) -> None:
    if con.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,)).fetchone():
        return
    con.execute(
        "INSERT OR IGNORE INTO papers (pmid, gene_symbol, extraction_summary) "
        "VALUES (?, ?, 'figure-reader stub')",
        (pmid, gene),
    )


def _ensure_variant(
    con: sqlite3.Connection, gene: str, cdna: str | None, protein: str | None
) -> int:
    row = con.execute(
        """SELECT variant_id FROM variants
           WHERE gene_symbol=? AND cdna_notation IS ?
             AND protein_notation IS ? AND genomic_position IS NULL""",
        (gene, cdna, protein),
    ).fetchone()
    if row:
        return row[0]
    cur = con.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES (?, ?, ?)",
        (gene, cdna, protein),
    )
    return int(cur.lastrowid)


def ingest_report(
    report: PMIDFigureReport, db_path: Path, source_tag: str = "figure-reader"
) -> int:
    """Write extracted variants to *db_path*. Returns number of new variant_papers rows."""
    if not report.distinct_variants:
        return 0
    con = sqlite3.connect(str(db_path))
    try:
        _ensure_paper(con, report.pmid, report.gene)
        added = 0
        for v in report.distinct_variants:
            cdna = (v.get("cdna") or "").strip() or None
            protein = (v.get("protein") or "").strip() or None
            if not (cdna or protein):
                continue
            vid = _ensure_variant(con, report.gene, cdna, protein)
            exists = con.execute(
                "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?",
                (vid, report.pmid),
            ).fetchone()
            if exists:
                continue
            note = json.dumps(
                {k: v.get(k) for k in ("carriers", "affected", "unaffected", "context")}
            )
            con.execute(
                """INSERT INTO variant_papers
                   (variant_id, pmid, source_location, additional_notes, key_quotes)
                   VALUES (?, ?, ?, ?, ?)""",
                (vid, report.pmid, source_tag, note, "[]"),
            )
            added += 1
        con.commit()
        return added
    finally:
        con.close()


def run(args: argparse.Namespace) -> int:
    pmids = _load_pmid_list(args)
    out_dir = Path(args.out).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    pmc_dir = Path(args.pmc_dir).expanduser()

    summary = {
        "gene": args.gene,
        "pmc_dir": str(pmc_dir),
        "model": args.model,
        "pmids": [],
    }

    total_variants = 0
    total_db_added = 0
    for pmid in pmids:
        images = find_pmid_figures(pmc_dir, pmid)
        if not images:
            logger.warning("PMID %s: no figures on disk under %s", pmid, pmc_dir)
            summary["pmids"].append(
                {"pmid": pmid, "figures": 0, "variants": 0, "db_added": 0}
            )
            continue

        if args.max_images is not None:
            images = images[: args.max_images]

        logger.info("PMID %s: reading %d figure(s)", pmid, len(images))
        report = read_figures_for_pmid(
            pmid=pmid,
            gene=args.gene,
            pmc_fulltext_dir=pmc_dir,
            model=args.model,
            max_images=args.max_images,
        )

        report_file = out_dir / f"{pmid}.json"
        report_file.write_text(json.dumps(report.to_dict(), indent=2))
        variants = report.distinct_variants
        total_variants += len(variants)
        db_added = 0
        if args.db and variants:
            db_added = ingest_report(report, Path(args.db))
            total_db_added += db_added
        logger.info(
            "PMID %s: %d distinct variant(s), %d added to DB",
            pmid,
            len(variants),
            db_added,
        )
        summary["pmids"].append(
            {
                "pmid": pmid,
                "figures": len(images),
                "variants": len(variants),
                "db_added": db_added,
            }
        )

    summary["total_variants"] = total_variants
    summary["total_db_added"] = total_db_added
    summary_file = out_dir / "summary.json"
    summary_file.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--gene", required=True, help="Gene symbol, e.g. KCNH2")
    p.add_argument(
        "--pmc-dir",
        required=True,
        help="Directory containing {PMID}_figures/ subdirectories",
    )
    p.add_argument(
        "--out", required=True, help="Output directory for per-PMID JSON reports"
    )
    p.add_argument(
        "--pmid",
        action="append",
        default=None,
        help="PMID to process (repeat for multiple)",
    )
    p.add_argument(
        "--pmid-file", help="File with one PMID per line ('#' comments allowed)"
    )
    p.add_argument(
        "--missing-csv",
        help="CSV with a 'pmid' column (e.g. recall missing_in_sqlite.csv)",
    )
    p.add_argument(
        "--db",
        help="Optional SQLite DB to insert extracted variants into",
    )
    p.add_argument(
        "--model",
        default=None,
        help="Vision model override (defaults to config.settings.get_vision_model)",
    )
    p.add_argument(
        "--max-images",
        type=int,
        default=None,
        help="Per-PMID image cap (default: all)",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
