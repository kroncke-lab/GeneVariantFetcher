"""Merge variant_papers from an older v12 DB into the current DB for gold PMIDs where v12 has more variants.

Usage:
    python scripts/recall_recovery/merge_v12_db.py \
      --v12-db results/KCNH2/.../KCNH2_v12_manual_recovery_20260515.db \
      --target-db results/KCNH2/.../KCNH2.db \
      --gold-pmids results/KCNH2/.../KCNH2_gold_pmids.txt
"""

import argparse
import sqlite3
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v12-db", required=True)
    ap.add_argument("--target-db", required=True)
    ap.add_argument("--gold-pmids", required=True)
    args = ap.parse_args()

    gold = set(Path(args.gold_pmids).read_text().split())
    print(f"Gold PMIDs: {len(gold)}")

    v12 = sqlite3.connect(args.v12_db)
    e2e = sqlite3.connect(args.target_db)
    v12.row_factory = sqlite3.Row
    e2e.row_factory = sqlite3.Row

    need_merge = []
    for pmid in gold:
        v = v12.execute(
            "SELECT COUNT(*) FROM variant_papers WHERE pmid=?", (pmid,)
        ).fetchone()[0]
        e = e2e.execute(
            "SELECT COUNT(*) FROM variant_papers WHERE pmid=?", (pmid,)
        ).fetchone()[0]
        if v > e:
            need_merge.append((pmid, v, e))
    print(f"PMIDs to merge: {len(need_merge)}")

    def ensure_paper(pmid):
        if e2e.execute("SELECT pmid FROM papers WHERE pmid=?", (pmid,)).fetchone():
            return
        p = v12.execute("SELECT * FROM papers WHERE pmid=?", (pmid,)).fetchone()
        if not p:
            return
        cols = list(p.keys())
        e2e.execute(
            f"INSERT OR REPLACE INTO papers ({','.join(cols)}) VALUES ({','.join('?' for _ in cols)})",
            tuple(p),
        )

    def ensure_variant(v12_row):
        key = (
            v12_row["gene_symbol"],
            v12_row["cdna_notation"],
            v12_row["protein_notation"],
            v12_row["genomic_position"],
        )
        e_row = e2e.execute(
            """SELECT variant_id FROM variants
                WHERE gene_symbol IS ? AND cdna_notation IS ?
                  AND protein_notation IS ? AND genomic_position IS ?""",
            key,
        ).fetchone()
        if e_row:
            return e_row["variant_id"]
        cur = e2e.execute(
            """INSERT INTO variants (gene_symbol, cdna_notation, protein_notation,
                                     genomic_position, clinical_significance, evidence_level)
               VALUES (?,?,?,?,?,?)""",
            (
                v12_row["gene_symbol"],
                v12_row["cdna_notation"],
                v12_row["protein_notation"],
                v12_row["genomic_position"],
                v12_row["clinical_significance"],
                v12_row["evidence_level"],
            ),
        )
        return cur.lastrowid

    e2e.execute("BEGIN")
    added = 0
    for pmid, v_count, e_count in need_merge:
        ensure_paper(pmid)
        rows = v12.execute(
            """SELECT v.*, vp.source_location, vp.additional_notes, vp.key_quotes
               FROM variants v JOIN variant_papers vp ON v.variant_id = vp.variant_id
               WHERE vp.pmid = ?""",
            (pmid,),
        ).fetchall()
        for r in rows:
            new_vid = ensure_variant(r)
            if not e2e.execute(
                "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?",
                (new_vid, pmid),
            ).fetchone():
                e2e.execute(
                    """INSERT INTO variant_papers
                       (variant_id, pmid, source_location, additional_notes, key_quotes)
                       VALUES (?,?,?,?,?)""",
                    (
                        new_vid,
                        pmid,
                        r["source_location"],
                        r["additional_notes"],
                        r["key_quotes"],
                    ),
                )
                added += 1
        # Also copy phenotypes / individual_records / penetrance_data
        for tbl in ("phenotypes", "individual_records", "penetrance_data"):
            rows_t = v12.execute(
                f"""SELECT t.*, v.gene_symbol as g, v.cdna_notation as c,
                           v.protein_notation as p, v.genomic_position as gp
                    FROM {tbl} t JOIN variants v ON v.variant_id = t.variant_id
                    WHERE t.pmid = ?""",
                (pmid,),
            ).fetchall()
            for row in rows_t:
                key = (row["g"], row["c"], row["p"], row["gp"])
                mapped = e2e.execute(
                    """SELECT variant_id FROM variants
                        WHERE gene_symbol IS ? AND cdna_notation IS ?
                          AND protein_notation IS ? AND genomic_position IS ?""",
                    key,
                ).fetchone()
                if not mapped:
                    continue
                d = dict(row)
                d["variant_id"] = mapped["variant_id"]
                for k in ("g", "c", "p", "gp"):
                    d.pop(k, None)
                for idcol in ("id", "phenotype_id", "record_id", "penetrance_id"):
                    d.pop(idcol, None)
                cols = list(d.keys())
                try:
                    e2e.execute(
                        f"INSERT INTO {tbl} ({','.join(cols)}) VALUES ({','.join('?' for _ in cols)})",
                        tuple(d[c] for c in cols),
                    )
                except sqlite3.IntegrityError:
                    pass
    e2e.commit()
    print(f"Added {added} variant_papers rows across {len(need_merge)} PMIDs")


if __name__ == "__main__":
    main()
