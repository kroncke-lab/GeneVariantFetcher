#!/usr/bin/env python3
"""Ingest PubTator3 text-mined mutations for PMIDs already present in a DB.

For each PMID already present in the SQLite DB, query NCBI PubTator3 for mutation
annotations, filter to the target gene by NCBI Gene ID, and insert any
new (PMID, variant) pairs into the recall DB as ``variant_papers`` rows
tagged ``PubTator3 (text-mined)``.

Gene-agnostic: pass ``--gene`` and the matching DB. Gold-PMID mode still exists
as an explicit evaluation/debug option, but it is intentionally not the default
because using the recall target to choose enrichment PMIDs leaks gold-set
knowledge into the scored DB.

Usage::

    python scripts/recall_recovery/ingest_pubtator.py \\
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
from collections import defaultdict
from pathlib import Path
from typing import Optional

import requests
from dotenv import load_dotenv

logger = logging.getLogger("ingest_pubtator")

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBTATOR_EXPORT = (
    "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson"
)


# ---------------------------------------------------------------------------
# Gene-ID resolution (reused pattern from ingest_clinvar.py)
# ---------------------------------------------------------------------------


def resolve_gene_id(gene: str, email: str, api_key: str) -> Optional[str]:
    params = {
        "db": "gene",
        "term": f"{gene}[Gene Symbol] AND human[Organism]",
        "retmode": "json",
        "retmax": 5,
        "tool": "GVF",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    try:
        r = requests.get(NCBI_ESEARCH, params=params, timeout=30)
    except Exception as exc:
        logger.warning("Gene ID lookup failed for %s: %s", gene, exc)
        return None
    if r.status_code != 200:
        return None
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else None


# ---------------------------------------------------------------------------
# PubTator client
# ---------------------------------------------------------------------------


def fetch_pubtator_mutations(
    pmids: list[str], gene_id: str, batch_size: int = 50
) -> dict[str, list[tuple[str, str]]]:
    """Return ``{pmid: [(text, hgvs), ...]}`` for *gene_id* across *pmids*."""
    out: dict[str, list[tuple[str, str]]] = defaultdict(list)
    gene_marker = f"CorrespondingGene:{gene_id}"
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        url = f"{PUBTATOR_EXPORT}?pmids={','.join(batch)}"
        try:
            r = requests.get(url, timeout=120)
        except Exception as exc:
            logger.warning("Batch %d: %s", i, exc)
            continue
        if r.status_code != 200:
            logger.warning("Batch %d: status %d", i, r.status_code)
            continue
        try:
            data = json.loads(r.text)
        except json.JSONDecodeError:
            continue
        for doc in data.get("PubTator3", []):
            pmid = doc.get("id")
            for passage in doc.get("passages", []):
                for ann in passage.get("annotations", []):
                    t = ann.get("infons", {}).get("type", "")
                    if t not in {
                        "Mutation",
                        "DNAMutation",
                        "ProteinMutation",
                        "SNP",
                        "Variant",
                    }:
                        continue
                    text = ann.get("text", "")
                    identifier = ann.get("infons", {}).get("identifier", "")
                    if gene_marker not in identifier:
                        continue
                    hgvs_m = re.search(r"HGVS:([pc]\.[^;]+)", identifier)
                    hgvs = hgvs_m.group(1) if hgvs_m else text
                    out[pmid].append((text, hgvs))
        time.sleep(0.5)
        if i % 200 == 0:
            logger.info("  progress: %d / %d", i, len(pmids))
    return out


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
        "VALUES (?, ?, 'PubTator stub')",
        (pmid, gene),
    )
    return True


def ensure_source_layer_column(con: sqlite3.Connection) -> None:
    columns = {row[1] for row in con.execute("PRAGMA table_info(variant_papers)")}
    if "source_layer" not in columns:
        con.execute("ALTER TABLE variant_papers ADD COLUMN source_layer TEXT")


def parse_hgvs(hgvs: str) -> tuple[Optional[str], Optional[str]]:
    if hgvs.startswith("p."):
        return None, hgvs
    if hgvs.startswith("c."):
        return hgvs, None
    return None, None


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
        "--gene-id",
        default=None,
        help="Override NCBI Gene ID (else auto-resolved from --gene)",
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

    gene_id = args.gene_id or resolve_gene_id(gene, email, api_key)
    if not gene_id:
        sys.exit(
            f"Could not resolve NCBI Gene ID for {gene}; pass --gene-id explicitly"
        )
    logger.info("Gene %s → NCBI Gene ID %s", gene, gene_id)

    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    ensure_source_layer_column(con)

    pmids = (
        load_gold_pmids(gold_path)
        if args.pmid_source == "gold" and gold_path
        else load_db_pmids(con, gene)
    )
    logger.info("%s PMIDs: %d", args.pmid_source.upper(), len(pmids))

    pmid_list = sorted(pmids)
    all_muts = fetch_pubtator_mutations(pmid_list, gene_id)
    logger.info(
        "PMIDs with %s mutations: %d (total annotations: %d)",
        gene,
        len(all_muts),
        sum(len(v) for v in all_muts.values()),
    )

    added = 0
    con.execute("BEGIN")
    for pmid, muts in all_muts.items():
        if pmid not in pmids:
            continue
        if not ensure_paper(con, pmid, gene, allow_stub=args.allow_stub_papers):
            continue
        for text, hgvs in muts:
            cdna, protein = parse_hgvs(hgvs)
            if not (cdna or protein):
                text_clean = text.replace("p.", "")
                if re.match(r"^[A-Z]\d+[A-Z*X]+$", text_clean):
                    protein = f"p.{text_clean}"
                elif text_clean.startswith("c."):
                    cdna = text_clean
            if not (cdna or protein):
                continue
            vid = ensure_variant(con, gene, cdna, protein)
            if not con.execute(
                "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?",
                (vid, pmid),
            ).fetchone():
                con.execute(
                    """INSERT INTO variant_papers
                       (variant_id, pmid, source_location, additional_notes, source_layer)
                       VALUES (?, ?, 'PubTator3 (text-mined)', ?, 'pubtator')""",
                    (vid, pmid, text),
                )
                added += 1
    con.commit()
    logger.info("Added %d variant_papers from PubTator3", added)
    print(
        json.dumps(
            {
                "gene": gene,
                "gene_id": gene_id,
                "added": added,
                "pmids": len(pmids),
                "pmid_source": args.pmid_source,
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
