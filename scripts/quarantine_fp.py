#!/usr/bin/env python3
"""Quarantine variantFeatures-flagged false-positive variants from a GVF run DB.

Default target is `misparse_out_of_range` only — variants whose protein residue
exceeds the gene's length, i.e. wrong-gene contamination (BRCA2-in-BRCA1,
APP/MTHFR-in-APOE) or cDNA-as-protein misparses. These are unambiguous.

It deliberately does NOT touch `no_notation_suspect` by default: that bucket
holds real-but-unnotated variants (promoter SNPs by rsID, ε-alleles by name).

Reversible: backs up the DB, copies the removed variants into a
`quarantined_variants` table, then deletes them from every variant-keyed table.
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
from datetime import datetime
from pathlib import Path


def variant_keyed_tables(con) -> list[str]:
    out = []
    for (name,) in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall():
        if name in ("vf_enrichment", "quarantined_variants", "sqlite_sequence"):
            continue
        cols = {r[1] for r in con.execute(f"PRAGMA table_info({name})")}
        if "variant_id" in cols:
            out.append(name)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument(
        "--classes",
        default="misparse_out_of_range",
        help="comma-sep fp_class values to quarantine",
    )
    ap.add_argument("--stamp", required=True, help="timestamp for the backup filename")
    args = ap.parse_args()
    classes = [c.strip() for c in args.classes.split(",") if c.strip()]

    con = sqlite3.connect(str(args.db))
    if not con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='vf_enrichment'"
    ).fetchone():
        print(
            f"  SKIP {args.db}: no vf_enrichment table (run enrich_from_variantfeatures first)"
        )
        return 1
    ph = ",".join("?" * len(classes))
    fp_ids = [
        r[0]
        for r in con.execute(
            f"SELECT variant_id FROM vf_enrichment WHERE fp_class IN ({ph})", classes
        )
    ]
    if not fp_ids:
        print(f"  {args.db.name}: 0 variants to quarantine for {classes}")
        return 0

    backup = args.db.with_suffix(f".before_fp_quarantine_{args.stamp}.db")
    shutil.copy2(args.db, backup)

    cur = con.cursor()
    cur.execute("""CREATE TABLE IF NOT EXISTS quarantined_variants(
        variant_id INTEGER, gene_symbol TEXT, protein_notation TEXT, cdna_notation TEXT,
        fp_class TEXT, n_papers INTEGER, quarantined_at TEXT)""")
    idset = ",".join(str(i) for i in fp_ids)
    cur.execute(
        f"""INSERT INTO quarantined_variants
        SELECT v.variant_id, v.gene_symbol, v.protein_notation, v.cdna_notation, e.fp_class,
               (SELECT COUNT(DISTINCT pmid) FROM variant_papers vp WHERE vp.variant_id=v.variant_id),
               ? FROM variants v JOIN vf_enrichment e ON e.variant_id=v.variant_id
        WHERE v.variant_id IN ({idset})""",
        (datetime.now().isoformat(timespec="seconds"),),
    )

    removed = {}
    for t in variant_keyed_tables(con):
        cur.execute(f"DELETE FROM {t} WHERE variant_id IN ({idset})")
        if cur.rowcount:
            removed[t] = cur.rowcount
    con.commit()
    v_after = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    con.close()

    print(f"  {args.db.name}: quarantined {len(fp_ids)} variants ({','.join(classes)})")
    print(
        f"    rows removed: "
        + ", ".join(f"{t}={n}" for t, n in sorted(removed.items()))
    )
    print(f"    variants remaining: {v_after} | backup: {backup.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
