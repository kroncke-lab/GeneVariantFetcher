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
from typing import Optional
import sqlite3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utils.gene_metadata import (  # noqa: E402
    resolve_variantfeatures_gene_symbols,
)

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


def load_vf(vf: Path, gene: str) -> tuple[dict, dict, dict, int, set[int]]:
    con = sqlite3.connect(f"file:{vf}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        prot_map, cdna_map, meta = {}, {}, {}
        # Resolve the stored gene-symbol casing once, then bind it as an equality
        # test in both queries below. `upper(gene_symbol)=?` cannot use
        # idx_consequences_gene, so it degraded each query into a full scan --
        # ~5 minutes on the pathogenicity join for byte-identical output.
        symbols = resolve_variantfeatures_gene_symbols(con, gene)
        if not symbols:
            # Gene absent from variantFeatures -- the normal case for a gene with
            # no slice yet. Returning here skips both scans, including the join
            # over annotations_pathogenicity (43.5M rows). Empty maps and a max
            # position of 0 are what the callers already expect for this case:
            # everything reports unmatched and nothing is ever quarantined.
            return {}, {}, {}, 0, set(), {}
        gene_ph = ",".join("?" * len(symbols))
        # Representative consequence per variant: prefer MANE select, then canonical.
        best: dict[int, tuple[int, sqlite3.Row]] = {}
        for r in con.execute(
            "SELECT variant_id, hgvs_c, hgvs_p, aa_pos, aa_ref, aa_alt, consequence, "
            "is_canonical, is_mane_select FROM variant_consequences "
            f"WHERE gene_symbol IN ({gene_ph})",
            symbols,
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
            f"JOIN (SELECT DISTINCT variant_id FROM variant_consequences "
            f"WHERE gene_symbol IN ({gene_ph})) c "
            f"ON t.variant_id=c.variant_id WHERE t.predictor IN ({ph})",
            (*symbols, *PREDICTORS),
        ):
            if r["score"] is not None:
                scores.setdefault(r["variant_id"], {})[r["predictor"]] = r["score"]
        for vid in meta:
            meta[vid]["scores"] = scores.get(vid, {})
        # Canonical protein length (max aa_pos) + the set of positions VF knows a
        # variant at — used to classify unmatched GVF variants (out-of-range = FP).
        max_pos = max((m["aa_pos"] for m in meta.values() if m["aa_pos"]), default=0)
        positions = {m["aa_pos"] for m in meta.values() if m["aa_pos"]}
        # Reference residue by position. Out-of-range was the only wrong-gene
        # signal here, and it is the weaker half: a BRCA1 polymorphism such as
        # P871L sits comfortably inside BRCA2's 3,418 residues, so it was
        # classified "novel_in_range" and never flagged even though BRCA2's
        # residue 871 is not P. Residue identity catches that class.
        residues: dict = {}
        for m in meta.values():
            if m["aa_pos"] and m["aa_ref"]:
                residues.setdefault(int(m["aa_pos"]), set()).add(
                    str(m["aa_ref"]).upper()
                )
        return prot_map, cdna_map, meta, max_pos, positions, residues
    finally:
        con.close()


_RES_RE = re.compile(r"(\d+)")


_PROT_CALL_RE = re.compile(
    r"^([ACDEFGHIKLMNPQRSTVWY])(\d{1,5})(?![0-9])", re.IGNORECASE
)


def classify_unmatched(
    protein: str,
    cdna: str,
    max_pos: int,
    residues: Optional[dict] = None,
    other_gene_residues: Optional[dict] = None,
) -> str:
    """High-confidence FP signal from variantFeatures non-coverage."""
    p = (protein or "").strip()
    c = (cdna or "").strip()
    if not p:
        return "cdna_only_unmatched" if "c." in c else "no_notation_suspect"
    bare = p.replace("p.", "")
    m = _RES_RE.search(bare)
    if not m:
        return "no_notation_suspect"
    if max_pos and int(m.group(1)) > max_pos:
        return (
            "misparse_out_of_range"  # residue can't exist in this gene -> FP/wrong-gene
        )
    # In range, but is the REFERENCE residue right for this gene? A wrong-gene
    # row usually lands in range and is only betrayed by the residue: BRCA1's
    # P871/E1038/K1183 haplotype under a BRCA2 run scores in-range at every one.
    call = _PROT_CALL_RE.match(bare)
    if call and residues:
        aa = call.group(1).upper()
        pos = int(call.group(2))
        known = residues.get(pos)
        if known and aa not in known:
            # A residue that disagrees is not automatically the wrong gene.
            # Legacy BIC protein numbering is offset from RefSeq, so an older
            # paper's call can land a few residues off and still be this gene's
            # variant. Quarantining those would delete real data: 125 BRCA1,
            # 29 BRCA2 and 24 BMPR2 rows in the 150-paper re-extraction fit
            # their own gene once a small offset is allowed.
            for off in (-3, -2, -1, 1, 2, 3):
                if aa in residues.get(pos + off, set()):
                    return "residue_offset_suspect"
            if other_gene_residues:
                for g, rmap in other_gene_residues.items():
                    if aa in rmap.get(pos, set()):
                        # Positively belongs to a different gene: the confident
                        # wrong-gene class, and the only one safe to quarantine.
                        return "wrong_gene_residue_mismatch"
            # Disagrees with this gene and matches no other gene we know.
            # Suspect, but not demonstrably misattributed — report, don't remove.
            return "residue_unverified"
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

    prot_map, cdna_map, meta, max_pos, vpos, residues = load_vf(args.vf, args.gene)
    # Residue maps for the other genes GVF knows, so a mismatch can be told
    # apart from a positive match elsewhere. Bounded and indexed by gene_symbol.
    other_residues: dict = {}
    try:
        from utils.gene_metadata import known_gene_aliases

        others = [
            g
            for g in known_gene_aliases(include_query_aliases=False)
            if g.upper() != args.gene.upper()
        ]
        con = sqlite3.connect(f"file:{args.vf}?mode=ro", uri=True)
        for g in others:
            rmap: dict = {}
            for pos, aa in con.execute(
                "SELECT DISTINCT aa_pos, aa_ref FROM variant_consequences "
                "WHERE gene_symbol=? AND aa_pos IS NOT NULL AND aa_ref IS NOT NULL",
                (g,),
            ):
                rmap.setdefault(int(pos), set()).add(str(aa).upper())
            if rmap:
                other_residues[g] = rmap
        con.close()
    except Exception as e:  # noqa: BLE001 - QC must never fail the run
        print(f"  (other-gene residue maps unavailable: {e})")
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

    # `cur` is the GVF *run* DB, not variantFeatures -- a different `variants`
    # table that does have gene_symbol (variantFeatures.variants does not).
    # upper() stays here on purpose: migrate_to_sqlite stores gene_symbol exactly
    # as extraction emitted it, with no case normalization, so the fold is
    # load-bearing. It costs nothing to leave -- the run DB holds one gene's
    # variants (thousands of rows, not millions).
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
            cls = classify_unmatched(
                r["protein_notation"],
                r["cdna_notation"],
                max_pos,
                residues,
                other_residues,
            )
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
        WHERE e.fp_class IN ('misparse_out_of_range','no_notation_suspect','wrong_gene_residue_mismatch','residue_unverified')
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
