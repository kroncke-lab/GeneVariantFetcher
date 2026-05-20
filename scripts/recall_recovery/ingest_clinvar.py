#!/usr/bin/env python3
"""Ingest ClinVar variants citing PMIDs already present in a recall DB.

For each PMID already present in the SQLite DB, query ClinVar via E-Utils for
variants in the target gene that cite that PMID. Insert anything new
into the DB as ``variant_papers`` rows tagged ``ClinVar (PMID citation)``.

Gene-agnostic: pass ``--gene`` and the matching DB. Gold-PMID mode still exists
as an explicit evaluation/debug option, but it is intentionally not the default
because using the recall target to choose enrichment PMIDs leaks gold-set
knowledge into the scored DB.

Usage::

    python scripts/recall_recovery/ingest_clinvar.py \\
        --gene KCNQ1 \\
        --db results/KCNQ1/<timestamp>/KCNQ1.db
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import sqlite3
import sys
import time
from pathlib import Path
from typing import Optional

import requests
from dotenv import load_dotenv

logger = logging.getLogger("ingest_clinvar")

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


# ---------------------------------------------------------------------------
# Variant-name parsing
# ---------------------------------------------------------------------------


def parse(name: str) -> tuple[Optional[str], Optional[str]]:
    """Extract cDNA + protein notation from a ClinVar variation title."""
    cdna_m = re.search(r"c\.([^\s)]+)", name)
    prot_m = re.search(r"p\.\(?([A-Za-z]+\d+[A-Za-z*?_=>\-]*)\)?", name)
    return (
        f"c.{cdna_m.group(1)}" if cdna_m else None,
        f"p.{prot_m.group(1)}" if prot_m else None,
    )


# ---------------------------------------------------------------------------
# NCBI client
# ---------------------------------------------------------------------------


def _eutils_params(email: str, api_key: str) -> dict:
    params = {"tool": "GVF", "email": email}
    if api_key:
        params["api_key"] = api_key
    return params


def _throttle(api_key: str) -> None:
    time.sleep(0.12 if api_key else 0.4)


def resolve_gene_id(gene: str, email: str, api_key: str) -> Optional[str]:
    """Return the NCBI Gene ID for *gene* in Homo sapiens.

    Returns ``None`` if the gene symbol can't be resolved.
    """
    params = {
        "db": "gene",
        "term": f"{gene}[Gene Symbol] AND human[Organism]",
        "retmode": "json",
        "retmax": 5,
        **_eutils_params(email, api_key),
    }
    try:
        r = requests.get(NCBI_ESEARCH, params=params, timeout=30)
    except Exception as exc:
        logger.warning("Gene ID lookup failed for %s: %s", gene, exc)
        return None
    if r.status_code != 200:
        return None
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    return ids[0]


def variants_for_pmid(
    pmid: str, gene: str, email: str, api_key: str, max_retries: int = 5
) -> Optional[list[tuple[Optional[str], Optional[str], str]]]:
    """Return ``[(cdna, protein, title), ...]`` for *pmid* in *gene*'s ClinVar."""
    for attempt in range(max_retries):
        params = {
            "db": "clinvar",
            "term": f"{pmid}[PMID] AND {gene}[gene]",
            "retmode": "json",
            "retmax": 500,
            **_eutils_params(email, api_key),
        }
        try:
            r = requests.get(NCBI_ESEARCH, params=params, timeout=30)
        except Exception:
            time.sleep(2 * (attempt + 1))
            continue
        if r.status_code == 429:
            time.sleep(2 * (attempt + 1))
            continue
        if r.status_code != 200:
            return None
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        _throttle(api_key)

        all_variants: list[tuple[Optional[str], Optional[str], str]] = []
        for batch_start in range(0, len(ids), 200):
            batch = ids[batch_start : batch_start + 200]
            params2 = {
                "db": "clinvar",
                "id": ",".join(batch),
                "retmode": "json",
                **_eutils_params(email, api_key),
            }
            r2 = requests.get(NCBI_ESUMMARY, params=params2, timeout=30)
            if r2.status_code != 200:
                continue
            data = r2.json().get("result", {})
            for vid in batch:
                info = data.get(vid)
                if not info:
                    continue
                title = info.get("title", "")
                cdna, protein = parse(title)
                for v in info.get("variation_set", []):
                    name = v.get("variation_name", "") or v.get("cdna_change", "")
                    c2, p2 = parse(name)
                    if c2 and not cdna:
                        cdna = c2
                    if p2 and not protein:
                        protein = p2
                if cdna or protein:
                    all_variants.append((cdna, protein, title))
            _throttle(api_key)
        return all_variants
    return None


# ---------------------------------------------------------------------------
# DB inserts
# ---------------------------------------------------------------------------


def ensure_paper(
    con: sqlite3.Connection, pmid: str, gene: str, *, allow_stub: bool
) -> bool:
    if con.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,)).fetchone():
        return True
    if not allow_stub:
        return False
    con.execute(
        "INSERT OR IGNORE INTO papers (pmid, gene_symbol, extraction_summary) "
        "VALUES (?, ?, 'ClinVar-only stub')",
        (pmid, gene),
    )
    return True


def ensure_variant(
    con: sqlite3.Connection, gene: str, cdna: Optional[str], protein: Optional[str]
) -> int:
    row = con.execute(
        """SELECT variant_id FROM variants
           WHERE gene_symbol=? AND cdna_notation IS ?
             AND protein_notation IS ? AND genomic_position IS NULL""",
        (gene, cdna, protein),
    ).fetchone()
    if row:
        return int(row[0])
    cur = con.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES (?, ?, ?)",
        (gene, cdna, protein),
    )
    return int(cur.lastrowid)


def load_gold_pmids(gold_path: Path) -> set[str]:
    gold_pmids: set[str] = set()
    with gold_path.open() as f:
        for r in csv.DictReader(f):
            if r.get("pmid"):
                gold_pmids.add(r["pmid"])
    return gold_pmids


def load_db_pmids(con: sqlite3.Connection, gene: str) -> set[str]:
    """Return PMIDs the extraction DB already knows about for this gene."""

    try:
        rows = con.execute(
            """SELECT DISTINCT pmid FROM papers
               WHERE pmid IS NOT NULL AND TRIM(pmid) != ''
                 AND (gene_symbol IS NULL OR UPPER(gene_symbol)=?)""",
            (gene,),
        ).fetchall()
    except sqlite3.OperationalError:
        rows = []
    pmids = {str(r[0]) for r in rows if str(r[0]).isdigit()}
    if pmids:
        return pmids
    try:
        rows = con.execute(
            "SELECT DISTINCT pmid FROM variant_papers WHERE pmid IS NOT NULL"
        ).fetchall()
    except sqlite3.OperationalError:
        rows = []
    return {str(r[0]) for r in rows if str(r[0]).isdigit()}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--gene", required=True, help="Gene symbol (e.g. KCNH2, KCNQ1)")
    p.add_argument("--db", required=True, help="Path to the recall SQLite DB")
    p.add_argument(
        "--gold",
        default=None,
        help="Path to the gold-standard recall CSV. Required only with --pmid-source gold.",
    )
    p.add_argument(
        "--pmid-source",
        choices=("db", "gold"),
        default="db",
        help=(
            "Which PMID set to query. Default 'db' enriches only papers already "
            "observed by the extraction DB. 'gold' is evaluation-aided and should "
            "not be used for cold-start recall claims."
        ),
    )
    p.add_argument(
        "--allow-stub-papers",
        action="store_true",
        help="Allow inserting paper rows for PMIDs absent from the DB.",
    )
    p.add_argument(
        "--email", default=None, help="Override NCBI email (else NCBI_EMAIL env var)"
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    load_dotenv()

    gene = args.gene.upper()
    db_path = Path(args.db)
    gold_path = Path(args.gold) if args.gold else None
    email = args.email or os.environ.get("NCBI_EMAIL")
    api_key = os.environ.get("NCBI_API_KEY", "")

    if not db_path.exists():
        sys.exit(f"DB not found: {db_path}")
    if args.pmid_source == "gold" and not gold_path:
        sys.exit("--gold is required when --pmid-source=gold")
    if gold_path and not gold_path.exists():
        sys.exit(f"Gold CSV not found: {gold_path}")
    if not email:
        sys.exit("NCBI_EMAIL not set (export it or pass --email)")

    # Resolve gene ID once — used for logging only (the query uses the symbol).
    gene_id = resolve_gene_id(gene, email, api_key)
    logger.info("Gene %s → NCBI Gene ID %s", gene, gene_id or "(unresolved)")

    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row

    pmids = (
        load_gold_pmids(gold_path)
        if args.pmid_source == "gold" and gold_path
        else load_db_pmids(con, gene)
    )
    logger.info("%s PMIDs: %d", args.pmid_source.upper(), len(pmids))

    added = 0
    con.execute("BEGIN")
    for i, pmid in enumerate(sorted(pmids), 1):
        if i % 30 == 0:
            logger.info("  progress: %d / %d", i, len(pmids))
        res = variants_for_pmid(pmid, gene, email, api_key)
        if not res:
            continue
        if not ensure_paper(con, pmid, gene, allow_stub=args.allow_stub_papers):
            continue
        for cdna, protein, title in res:
            vid = ensure_variant(con, gene, cdna, protein)
            if not con.execute(
                "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?",
                (vid, pmid),
            ).fetchone():
                con.execute(
                    """INSERT INTO variant_papers
                       (variant_id, pmid, source_location, additional_notes, key_quotes)
                       VALUES (?, ?, 'ClinVar (PMID citation)', ?, ?)""",
                    (vid, pmid, title, json.dumps([])),
                )
                added += 1
    con.commit()
    logger.info("Added %d variant_papers across %d PMIDs", added, len(pmids))
    print(
        json.dumps(
            {
                "gene": gene,
                "added": added,
                "pmids": len(pmids),
                "pmid_source": args.pmid_source,
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
