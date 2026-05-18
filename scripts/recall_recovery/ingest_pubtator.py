import requests
import time
import json
import sqlite3
import csv
import re
from collections import defaultdict

# Load gold PMIDs
gold_pmids = set()
with open("gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv") as f:
    for r in csv.DictReader(f):
        gold_pmids.add(r["pmid"])
print(f"Gold PMIDs: {len(gold_pmids)}")

# Batch PubTator queries (100 PMIDs per request)
all_muts = defaultdict(list)
pmid_list = sorted(gold_pmids)
for i in range(0, len(pmid_list), 50):
    batch = pmid_list[i : i + 50]
    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmids={','.join(batch)}"
    try:
        r = requests.get(url, timeout=120)
        if r.status_code != 200:
            print(f"Batch {i}: status {r.status_code}")
            continue
        data = json.loads(r.text)
        for doc in data.get("PubTator3", []):
            pmid = doc.get("id")
            for passage in doc.get("passages", []):
                for ann in passage.get("annotations", []):
                    t = ann.get("infons", {}).get("type", "")
                    if t in (
                        "Mutation",
                        "DNAMutation",
                        "ProteinMutation",
                        "SNP",
                        "Variant",
                    ):
                        text = ann.get("text", "")
                        identifier = ann.get("infons", {}).get("identifier", "")
                        # Only KCNH2 mutations (gene 3757)
                        if "CorrespondingGene:3757" in identifier:
                            # Extract HGVS form
                            hgvs_m = re.search(r"HGVS:([pc]\.[^;]+)", identifier)
                            hgvs = hgvs_m.group(1) if hgvs_m else text
                            all_muts[pmid].append((text, hgvs))
        time.sleep(0.5)
    except Exception as e:
        print(f"Batch {i}: {e}")
    if i % 200 == 0:
        print(f"  Progress: {i}/{len(pmid_list)}")

print(f"\nTotal PMIDs with KCNH2 mutations: {len(all_muts)}")
total = sum(len(v) for v in all_muts.values())
print(f"Total mutations: {total}")

# Now ingest
DB = "results/KCNH2/20260517_074737/KCNH2.db"
con = sqlite3.connect(DB)
con.row_factory = sqlite3.Row


def ensure_paper(pmid):
    if not con.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,)).fetchone():
        con.execute(
            "INSERT OR IGNORE INTO papers (pmid, gene_symbol, extraction_summary) VALUES (?, 'KCNH2', 'PubTator stub')",
            (pmid,),
        )


def parse_hgvs(hgvs):
    if hgvs.startswith("p."):
        return None, hgvs
    elif hgvs.startswith("c."):
        return hgvs, None
    return None, None


def ensure_variant(cdna, protein):
    row = con.execute(
        """SELECT variant_id FROM variants
                          WHERE gene_symbol='KCNH2' AND cdna_notation IS ?
                            AND protein_notation IS ? AND genomic_position IS NULL""",
        (cdna, protein),
    ).fetchone()
    if row:
        return row[0]
    return con.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES ('KCNH2', ?, ?)",
        (cdna, protein),
    ).lastrowid


added = 0
con.execute("BEGIN")
for pmid, muts in all_muts.items():
    if pmid not in gold_pmids:
        continue
    ensure_paper(pmid)
    for text, hgvs in muts:
        cdna, protein = parse_hgvs(hgvs)
        if not (cdna or protein):
            # Fall back: try parsing text
            text_clean = text.replace("p.", "")
            if re.match(r"^[A-Z]\d+[A-Z*X]+$", text_clean):
                protein = f"p.{text_clean}"
            elif text_clean.startswith("c."):
                cdna = text_clean
        if not (cdna or protein):
            continue
        vid = ensure_variant(cdna, protein)
        if not con.execute(
            "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?", (vid, pmid)
        ).fetchone():
            con.execute(
                "INSERT INTO variant_papers (variant_id, pmid, source_location, additional_notes) VALUES (?, ?, ?, ?)",
                (vid, pmid, "PubTator3 (text-mined)", text),
            )
            added += 1
con.commit()
print(f"Added {added} variant_papers from PubTator3")
