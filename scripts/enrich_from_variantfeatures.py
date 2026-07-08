#!/usr/bin/env python3
"""Match a GVF run DB's variants against the variantFeatures annotation warehouse.

Three outputs, all GVF-side (the browser already builds its identities from the
same warehouse, so matching here canonicalizes names the publish/matcher uses):

  1. Canonical name  — link each GVF variant to its variantFeatures identity
     (canonical hgvs_p + 1-letter aa key + hgvs_c), via protein change then cDNA.
  2. In silico scores — attach annotations_pathogenicity predictors per variant.
  3. False positives — GVF variants with NO variantFeatures match for the gene
     are flagged (wrong-gene / mis-extraction candidates).

Writes a `vf_enrichment` table into the run DB and a FP-report CSV. Read-only on
variantFeatures. Idempotent (drops+rebuilds vf_enrichment).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
from pathlib import Path

VF_DEFAULT = Path.home() / "GitRepos" / "variantFeatures" / "data" / "variants.db"
# Headline predictors get their own columns; all present ones go in scores_json.
PREDICTORS = (
    "alphamissense",
    "revel",
    "cadd_phred",
    "metasvm",
    "primateai",
    "provean",
    "polyphen2_hvar",
    "eve",
    "clinpred",
    "bayesdel_add_af",
)
HEADLINE = ("alphamissense", "revel", "cadd_phred")

AA3TO1 = {
    "ala": "A",
    "arg": "R",
    "asn": "N",
    "asp": "D",
    "cys": "C",
    "gln": "Q",
    "glu": "E",
    "gly": "G",
    "his": "H",
    "ile": "I",
    "leu": "L",
    "lys": "K",
    "met": "M",
    "phe": "F",
    "pro": "P",
    "ser": "S",
    "thr": "T",
    "trp": "W",
    "tyr": "Y",
    "val": "V",
    "ter": "*",
    "sec": "U",
    "pyl": "O",
    "xaa": "X",
}
_LEAD_P = re.compile(r"^p\.?", re.I)


def parse_protein_key(notation: str) -> str | None:
    """GVF protein notation -> canonical 1-letter key 'C112R' / 'K1194*'. None if unparseable."""
    if not notation:
        return None
    s = _LEAD_P.sub("", str(notation).strip())
    if not s:
        return None
    # 3-letter missense / stop: Cys112Arg, Lys1194Ter/Ter
    m = re.fullmatch(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", s)
    if m:
        wt = AA3TO1.get(m.group(1).lower())
        mut = AA3TO1.get(m.group(3).lower())
        if wt and mut:
            return f"{wt}{m.group(2)}{mut}"
    # 3-letter stop with * / X / Ter
    m = re.fullmatch(r"([A-Za-z]{3})(\d+)(?:\*|Ter|X)", s)
    if m and AA3TO1.get(m.group(1).lower()):
        return f"{AA3TO1[m.group(1).lower()]}{m.group(2)}*"
    # 1-letter missense / stop: V629I, K1194X/*
    m = re.fullmatch(r"([A-Za-z])(\d+)([A-Za-z*])", s)
    if m:
        wt, pos, mut = m.group(1).upper(), m.group(2), m.group(3).upper()
        if wt in "ACDEFGHIKLMNPQRSTVWY":
            if mut in ("*", "X"):
                return f"{wt}{pos}*"
            if mut in "ACDEFGHIKLMNPQRSTVWY":
                return f"{wt}{pos}{mut}"
    return None


def norm_cdna(hgvs_c: str) -> str | None:
    """Strip transcript prefix: 'ENST..:c.416C>T' -> 'c.416C>T'. Bare 'c.416C>T' passes through."""
    if not hgvs_c:
        return None
    s = str(hgvs_c).strip()
    if ":" in s:
        s = s.split(":", 1)[1]
    s = s.strip()
    return s if s.startswith("c.") else None


def load_vf(vf: Path, gene: str) -> tuple[dict, dict, dict]:
    con = sqlite3.connect(f"file:{vf}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        prot_map, cdna_map, meta = {}, {}, {}
        # Representative consequence per variant: prefer MANE select, then canonical.
        best: dict[int, tuple[int, sqlite3.Row]] = {}
        for r in con.execute(
            "SELECT variant_id, hgvs_c, hgvs_p, aa_pos, aa_ref, aa_alt, consequence, "
            "is_canonical, is_mane_select FROM variant_consequences WHERE upper(gene_symbol)=?",
            (gene.upper(),),
        ):
            rank = (r["is_mane_select"] or 0) * 2 + (r["is_canonical"] or 0)
            cur = best.get(r["variant_id"])
            if cur is None or rank > cur[0]:
                best[r["variant_id"]] = (rank, r)
        for vid, (_rank, r) in best.items():
            meta[vid] = {
                "hgvs_p": r["hgvs_p"],
                "hgvs_c": r["hgvs_c"],
                "consequence": r["consequence"],
                "aa_ref": r["aa_ref"],
                "aa_pos": r["aa_pos"],
                "aa_alt": r["aa_alt"],
            }
            if r["aa_ref"] and r["aa_pos"] and r["aa_alt"]:
                key = f"{r['aa_ref']}{r['aa_pos']}{r['aa_alt']}".upper()
                prot_map.setdefault(key, vid)
            ck = norm_cdna(r["hgvs_c"])
            if ck:
                cdna_map.setdefault(ck, vid)
        # in silico scores
        ph = ",".join("?" * len(PREDICTORS))
        scores: dict[int, dict] = {}
        for r in con.execute(
            f"SELECT t.variant_id, t.predictor, t.score FROM annotations_pathogenicity t "
            f"JOIN (SELECT DISTINCT variant_id FROM variant_consequences WHERE upper(gene_symbol)=?) c "
            f"ON t.variant_id=c.variant_id WHERE t.predictor IN ({ph})",
            (gene.upper(), *PREDICTORS),
        ):
            if r["score"] is not None:
                scores.setdefault(r["variant_id"], {})[r["predictor"]] = r["score"]
        for vid in meta:
            meta[vid]["scores"] = scores.get(vid, {})
        # Canonical protein length (max aa_pos) + the set of positions VF knows a
        # variant at — used to classify unmatched GVF variants (out-of-range = FP).
        max_pos = max((m["aa_pos"] for m in meta.values() if m["aa_pos"]), default=0)
        positions = {m["aa_pos"] for m in meta.values() if m["aa_pos"]}
        return prot_map, cdna_map, meta, max_pos, positions
    finally:
        con.close()


_RES_RE = re.compile(r"(\d+)")


def classify_unmatched(protein: str, cdna: str, max_pos: int) -> str:
    """High-confidence FP signal from variantFeatures non-coverage."""
    p = (protein or "").strip()
    c = (cdna or "").strip()
    if not p:
        return "cdna_only_unmatched" if "c." in c else "no_notation_suspect"
    m = _RES_RE.search(p.replace("p.", ""))
    if not m:
        return "no_notation_suspect"
    if max_pos and int(m.group(1)) > max_pos:
        return (
            "misparse_out_of_range"  # residue can't exist in this gene -> FP/wrong-gene
        )
    return "novel_in_range"  # valid position, just not in the observed warehouse


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--vf", default=VF_DEFAULT, type=Path)
    ap.add_argument(
        "--fp-out", default=None, help="false-positive report CSV (default next to db)"
    )
    args = ap.parse_args()

    prot_map, cdna_map, meta, max_pos, vpos = load_vf(args.vf, args.gene)
    print(
        f"variantFeatures[{args.gene}]: {len(meta)} variants  "
        f"(protein keys={len(prot_map)}, cDNA keys={len(cdna_map)}, "
        f"scored={sum(1 for m in meta.values() if m['scores'])})"
    )

    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS vf_enrichment")
    cols = ", ".join(f"{p} REAL" for p in HEADLINE)
    cur.execute(f"""CREATE TABLE vf_enrichment(
        variant_id INTEGER PRIMARY KEY, vf_variant_id INTEGER, matched INTEGER,
        match_method TEXT, canonical_hgvs_p TEXT, canonical_aa_key TEXT,
        canonical_hgvs_c TEXT, consequence TEXT, {cols}, scores_json TEXT, fp_class TEXT)""")

    rows = list(
        cur.execute(
            "SELECT variant_id, protein_notation, cdna_notation FROM variants WHERE upper(gene_symbol)=?",
            (args.gene.upper(),),
        )
    )
    from collections import Counter

    n_prot = n_cdna = n_unmatched = 0
    fp_class_counts: Counter = Counter()
    out = []
    for r in rows:
        vid_g = r["variant_id"]
        pkey = parse_protein_key(r["protein_notation"])
        method, vf_id = None, None
        if pkey and pkey in prot_map:
            vf_id, method = prot_map[pkey], "protein"
            n_prot += 1
        else:
            ck = norm_cdna(r["cdna_notation"])
            if ck and ck in cdna_map:
                vf_id, method = cdna_map[ck], "cdna"
                n_cdna += 1
        if vf_id is None:
            n_unmatched += 1
            cls = classify_unmatched(r["protein_notation"], r["cdna_notation"], max_pos)
            fp_class_counts[cls] += 1
            out.append(
                (
                    vid_g,
                    None,
                    0,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    cls,
                )
            )
            continue
        m = meta[vf_id]
        sc = m["scores"]
        out.append(
            (
                vid_g,
                vf_id,
                1,
                method,
                m["hgvs_p"],
                f"{m['aa_ref']}{m['aa_pos']}{m['aa_alt']}" if m["aa_ref"] else None,
                m["hgvs_c"],
                m["consequence"],
                sc.get("alphamissense"),
                sc.get("revel"),
                sc.get("cadd_phred"),
                json.dumps(sc) if sc else None,
                None,
            )
        )
    cur.executemany(
        f"INSERT INTO vf_enrichment VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", out
    )
    con.commit()

    scored = cur.execute(
        "SELECT COUNT(*) FROM vf_enrichment WHERE matched=1 AND scores_json IS NOT NULL"
    ).fetchone()[0]
    total = len(rows)
    matched = n_prot + n_cdna
    print(
        f"GVF[{args.gene}]: {total} variants | matched {matched} "
        f"({100 * matched / total:.1f}%) = protein {n_prot} + cDNA {n_cdna} | with in-silico score {scored}"
    )
    print(
        f"  unmatched {n_unmatched}: "
        + ", ".join(f"{k}={v}" for k, v in sorted(fp_class_counts.items()))
    )

    # High-confidence FP report: out-of-range misparses + junk notation, with #papers (prioritize review).
    fp_out = (
        Path(args.fp_out)
        if args.fp_out
        else args.db.parent / f"{args.gene}_vf_false_positives.csv"
    )
    rep = cur.execute("""
        SELECT v.variant_id, v.protein_notation, v.cdna_notation, e.fp_class,
               (SELECT COUNT(DISTINCT pmid) FROM variant_papers vp WHERE vp.variant_id=v.variant_id) AS papers
        FROM vf_enrichment e JOIN variants v ON v.variant_id=e.variant_id
        WHERE e.fp_class IN ('misparse_out_of_range','no_notation_suspect')
        ORDER BY papers DESC, v.variant_id""").fetchall()
    with open(fp_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "gvf_variant_id",
                "protein_notation",
                "cdna_notation",
                "fp_class",
                "n_papers",
            ]
        )
        w.writerows(rep)
    con.close()
    print(f"  high-confidence false positives ({len(rep)}) → {fp_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
