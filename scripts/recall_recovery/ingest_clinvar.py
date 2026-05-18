import requests
import time
import sqlite3
import csv
import json
import re
from pathlib import Path
import os
from dotenv import load_dotenv

load_dotenv()

DB = "results/KCNH2/20260517_074737/KCNH2.db"
GOLD = "gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv"
NCBI_EMAIL = "brett.kroncke@vanderbilt.edu"
NCBI_KEY = os.environ.get("NCBI_API_KEY", "")

gold_pmids = set()
with open(GOLD) as f:
    for r in csv.DictReader(f):
        gold_pmids.add(r["pmid"])

con = sqlite3.connect(DB)
con.row_factory = sqlite3.Row


def parse(name):
    cdna_m = re.search(r"c\.([^\s)]+)", name)
    prot_m = re.search(r"p\.\(?([A-Za-z]+\d+[A-Za-z*?_=>\-]*)\)?", name)
    return (
        f"c.{cdna_m.group(1)}" if cdna_m else None,
        f"p.{prot_m.group(1)}" if prot_m else None,
    )


def variants_for_pmid(pmid, max_retries=5):
    for attempt in range(max_retries):
        params = {
            "db": "clinvar",
            "term": f"{pmid}[PMID] AND KCNH2[gene]",
            "retmode": "json",
            "retmax": 500,
            "tool": "GVF",
            "email": NCBI_EMAIL,
        }
        if NCBI_KEY:
            params["api_key"] = NCBI_KEY
        try:
            r = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params=params,
                timeout=30,
            )
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
        time.sleep(0.12 if NCBI_KEY else 0.4)

        # Get details
        all_variants = []
        for batch_start in range(0, len(ids), 200):
            batch = ids[batch_start : batch_start + 200]
            params2 = {
                "db": "clinvar",
                "id": ",".join(batch),
                "retmode": "json",
                "tool": "GVF",
                "email": NCBI_EMAIL,
            }
            if NCBI_KEY:
                params2["api_key"] = NCBI_KEY
            r2 = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params=params2,
                timeout=30,
            )
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
            time.sleep(0.12 if NCBI_KEY else 0.4)
        return all_variants
    return None


def ensure_paper(pmid):
    if not con.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,)).fetchone():
        con.execute(
            "INSERT OR IGNORE INTO papers (pmid, gene_symbol, extraction_summary) VALUES (?, 'KCNH2', 'ClinVar-only stub')",
            (pmid,),
        )


def ensure_variant(cdna, protein):
    e_row = con.execute(
        """SELECT variant_id FROM variants
                            WHERE gene_symbol='KCNH2' AND cdna_notation IS ?
                              AND protein_notation IS ? AND genomic_position IS NULL""",
        (cdna, protein),
    ).fetchone()
    if e_row:
        return e_row["variant_id"]
    cur = con.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES ('KCNH2', ?, ?)",
        (cdna, protein),
    )
    return cur.lastrowid


added = 0
con.execute("BEGIN")
for i, pmid in enumerate(sorted(gold_pmids), 1):
    if i % 30 == 0:
        print(f"  ({i}/{len(gold_pmids)})")
    res = variants_for_pmid(pmid)
    if not res:
        continue
    ensure_paper(pmid)
    for cdna, protein, title in res:
        vid = ensure_variant(cdna, protein)
        if not con.execute(
            "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?", (vid, pmid)
        ).fetchone():
            con.execute(
                """INSERT INTO variant_papers (variant_id, pmid, source_location, additional_notes, key_quotes)
                           VALUES (?, ?, 'ClinVar (PMID citation)', ?, ?)""",
                (vid, pmid, title, json.dumps([])),
            )
            added += 1
con.commit()
print(f"\nAdded {added} variant_papers across {len(gold_pmids)} PMIDs")
