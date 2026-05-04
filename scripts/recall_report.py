"""Comprehensive recall metrics across all dimensions vs the gold standard."""

import sys
import re
import sqlite3
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.variant_normalizer import create_variant_key  # noqa: E402


def main(run_dir: Path, gold_xls: Path, gene: str = "KCNH2") -> dict:
    db_path = run_dir / f"{gene}.db"
    summary_json = run_dir / f"{gene}_penetrance_summary.json"

    # ---------- gold standard ----------
    df = pd.read_excel(gold_xls)
    df = df[df["EXCLUDE?"] == 0].copy()
    df["pmid_str"] = df["PMID"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else None
    )
    df["lqt_num"] = pd.to_numeric(df["LQT"], errors="coerce")
    df["car_num"] = pd.to_numeric(df["CAR"], errors="coerce")
    df["una_num"] = pd.to_numeric(df["UNA"], errors="coerce")
    df["var_key"] = df["Variant"].apply(
        lambda v: create_variant_key({"protein_notation": str(v).strip()}, gene)
        if pd.notna(v)
        else None
    )

    gold_pmids = set(df["pmid_str"].dropna())
    gold_variants = set(df["var_key"].dropna())
    gold_total_rows = len(df)
    gold_total_carriers = int(df["car_num"].fillna(0).sum())
    gold_lqt_affected_rows = int((df["lqt_num"] >= 1).sum())
    gold_lqt_unaffected_rows = int((df["lqt_num"] == 0).sum())
    gold_clinically_phenotyped = int(df["lqt_num"].notna().sum())
    gold_unaffected_sum = int(df["una_num"].fillna(0).sum())

    # Per-variant counts in gold
    gold_per_variant = (
        df.dropna(subset=["var_key"])
        .groupby("var_key")
        .agg(
            carriers=("car_num", lambda s: int(s.fillna(0).sum())),
            affected=("lqt_num", lambda s: int((s >= 1).sum())),
            unaffected_lqt=("lqt_num", lambda s: int((s == 0).sum())),
            una_sum=("una_num", lambda s: int(s.fillna(0).sum())),
        )
    )

    # ---------- extracted (SQLite) ----------
    conn = sqlite3.connect(db_path)
    extr_pmids = set(
        r[0]
        for r in conn.execute("SELECT DISTINCT pmid FROM variant_papers").fetchall()
    )

    variant_rows = conn.execute(
        "SELECT variant_id, gene_symbol, cdna_notation, protein_notation FROM variants"
    ).fetchall()
    variant_id_to_key: dict[int, str] = {}
    for vid, vgene, cdna, prot in variant_rows:
        key = create_variant_key(
            {"gene_symbol": vgene, "cdna_notation": cdna, "protein_notation": prot},
            vgene or gene,
        )
        if key:
            variant_id_to_key[vid] = key
    extr_variants = set(variant_id_to_key.values())

    pen_rows = conn.execute(
        "SELECT variant_id, COALESCE(total_carriers_observed,0), "
        "COALESCE(affected_count,0), COALESCE(unaffected_count,0) "
        "FROM penetrance_data"
    ).fetchall()

    def _safe_int(value) -> int:
        try:
            return int(float(str(value).strip()))
        except (TypeError, ValueError):
            return 0

    extr_carriers_by_variant: dict[str, dict[str, int]] = {}
    for vid, tc, aff, unaff in pen_rows:
        key = variant_id_to_key.get(vid)
        if not key:
            continue
        d = extr_carriers_by_variant.setdefault(
            key, {"carriers": 0, "affected": 0, "unaffected": 0}
        )
        d["carriers"] += _safe_int(tc)
        d["affected"] += _safe_int(aff)
        d["unaffected"] += _safe_int(unaff)

    extr_total_carriers = sum(d["carriers"] for d in extr_carriers_by_variant.values())
    extr_total_affected = sum(d["affected"] for d in extr_carriers_by_variant.values())
    extr_total_unaffected = sum(
        d["unaffected"] for d in extr_carriers_by_variant.values()
    )

    # ---------- recall by dimension ----------
    pmid_overlap = gold_pmids & extr_pmids
    variant_overlap = gold_variants & extr_variants

    # Per-variant comparison: among matched variants, how often we report
    # nonzero affected and unaffected
    matched_with_affected = sum(
        1
        for v in variant_overlap
        if extr_carriers_by_variant.get(v, {}).get("affected", 0) > 0
    )
    matched_with_unaffected = sum(
        1
        for v in variant_overlap
        if extr_carriers_by_variant.get(v, {}).get("unaffected", 0) > 0
    )

    gold_variants_with_affected = sum(
        1 for v in gold_variants if int(gold_per_variant.loc[v, "affected"]) > 0
    )
    gold_variants_with_unaffected = sum(
        1 for v in gold_variants if int(gold_per_variant.loc[v, "unaffected_lqt"]) > 0
    )

    metrics = {
        "gold_pmids": len(gold_pmids),
        "extr_pmids": len(extr_pmids),
        "pmid_recall_count": len(pmid_overlap),
        "pmid_recall_pct": round(100 * len(pmid_overlap) / len(gold_pmids), 2),
        "gold_variants": len(gold_variants),
        "extr_variants": len(extr_variants),
        "variant_recall_count": len(variant_overlap),
        "variant_recall_pct": round(100 * len(variant_overlap) / len(gold_variants), 2),
        "gold_total_carriers": gold_total_carriers,
        "extr_total_carriers": extr_total_carriers,
        "gold_clinically_phenotyped": gold_clinically_phenotyped,
        "gold_lqt_affected_rows": gold_lqt_affected_rows,
        "gold_lqt_unaffected_rows": gold_lqt_unaffected_rows,
        "gold_unaffected_sum": gold_unaffected_sum,
        "extr_total_affected": extr_total_affected,
        "extr_total_unaffected": extr_total_unaffected,
        "gold_variants_with_affected": gold_variants_with_affected,
        "gold_variants_with_unaffected": gold_variants_with_unaffected,
        "matched_variants_with_affected": matched_with_affected,
        "matched_variants_with_unaffected": matched_with_unaffected,
    }
    return metrics


if __name__ == "__main__":
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    gold_xls = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else Path(
            "/Users/kronckbm/Desktop/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
        )
    )
    m = main(run_dir, gold_xls)
    print("\n=== KCNH2 Recall Report (run", run_dir.name, ") ===\n")
    print(
        f"PMID recall:     {m['pmid_recall_count']}/{m['gold_pmids']} = {m['pmid_recall_pct']}%"
    )
    print(
        f"Variant recall:  {m['variant_recall_count']}/{m['gold_variants']} = {m['variant_recall_pct']}%"
    )
    print()
    print(f"Gold total carriers (CAR sum):       {m['gold_total_carriers']}")
    print(f"Extracted total carriers:            {m['extr_total_carriers']}")
    print()
    print(f"Gold clinically phenotyped (LQT col):{m['gold_clinically_phenotyped']}")
    print(f"Gold LQT-affected rows:              {m['gold_lqt_affected_rows']}")
    print(f"Gold LQT-unaffected rows:            {m['gold_lqt_unaffected_rows']}")
    print(f"Gold UNA column sum (unaffected):    {m['gold_unaffected_sum']}")
    print()
    print(f"Extracted total affected:            {m['extr_total_affected']}")
    print(f"Extracted total unaffected:          {m['extr_total_unaffected']}")
    print()
    print(
        f"Gold variants with ≥1 affected carrier:    {m['gold_variants_with_affected']}"
    )
    print(
        f"Matched variants reporting affected count: {m['matched_variants_with_affected']}"
    )
    print(
        f"Gold variants with ≥1 unaffected carrier:  {m['gold_variants_with_unaffected']}"
    )
    print(
        f"Matched variants reporting unaffected count:{m['matched_variants_with_unaffected']}"
    )
    print()
