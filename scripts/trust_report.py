#!/usr/bin/env python3
"""Trust-tier report: how much of a gene DB's count evidence is trusted vs
quarantined, and why.

Read-only inspector for the two-tier DB the trust gate (``pipeline/trust_gate.py``)
produces. Use it to watch the quarantine rate across genes at scale — a rate that
explodes on a new gene means the gate is "working" by hiding everything, which is
exactly the fleet-acceptance signal ``docs/AUTONOMY_ROADMAP.md`` calls for.

    python scripts/trust_report.py --db results/KCNH2/latest/KCNH2.db [--list 20]
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any


def _connect_ro(db_path: str | Path) -> sqlite3.Connection:
    return sqlite3.connect(f"{Path(db_path).resolve().as_uri()}?mode=ro", uri=True)


def summarize_trust(db_path: str | Path) -> dict[str, Any]:
    """Return trust-tier counts for a DB's ``penetrance_data`` facts.

    ``{tiered, total, trusted, quarantine, quarantine_rate, by_reason,
    rule_versions, study_design_distribution, quarantine_by_study_design}``.
    ``tiered=False`` for a pre-trust-gate DB (no column).
    """
    conn = _connect_ro(db_path)
    try:
        cols = {r[1] for r in conn.execute("PRAGMA table_info(penetrance_data)")}
        total = conn.execute("SELECT COUNT(*) FROM penetrance_data").fetchone()[0]
        if "trust_tier" not in cols:
            return {"tiered": False, "total": total}
        em_cols = {r[1] for r in conn.execute("PRAGMA table_info(extraction_metadata)")}
        has_study_design = "study_design" in em_cols
        if has_study_design:
            rows = conn.execute(
                """
                SELECT pd.trust_tier,
                       pd.trust_reasons,
                       pd.trust_rule_version,
                       em.study_design
                FROM penetrance_data pd
                LEFT JOIN extraction_metadata em
                       ON em.pmid = pd.pmid
                      AND em.extraction_id = (
                              SELECT MAX(em2.extraction_id)
                              FROM extraction_metadata em2
                              WHERE em2.pmid = pd.pmid
                          )
                """
            ).fetchall()
        else:
            rows = conn.execute(
                "SELECT trust_tier, trust_reasons, trust_rule_version, NULL "
                "FROM penetrance_data"
            ).fetchall()
    finally:
        conn.close()

    tiers: Counter = Counter()
    reasons: Counter = Counter()
    versions: Counter = Counter()
    study_designs: Counter = Counter()
    quarantine_by_design: Counter = Counter()
    for tier, reasons_json, version, study_design in rows:
        tiers[tier or "untiered"] += 1
        versions[version or "none"] += 1
        design_key = (study_design or "").strip() or "unknown"
        study_designs[design_key] += 1
        if tier == "quarantine":
            quarantine_by_design[design_key] += 1
        if reasons_json:
            try:
                for reason in json.loads(reasons_json):
                    reasons[reason] += 1
            except (TypeError, ValueError):
                pass

    quarantine = tiers.get("quarantine", 0)
    return {
        "tiered": True,
        "total": total,
        "trusted": tiers.get("trusted", 0),
        "quarantine": quarantine,
        "untiered": tiers.get("untiered", 0),
        "quarantine_rate": (quarantine / total) if total else 0.0,
        "by_reason": dict(reasons.most_common()),
        "rule_versions": dict(versions),
        "study_design_distribution": dict(study_designs.most_common()),
        "quarantine_by_study_design": dict(quarantine_by_design.most_common()),
    }


def list_quarantined(db_path: str | Path, limit: int = 20) -> list[dict[str, Any]]:
    """The quarantined facts (variant, pmid, counts, reasons), most-carriers first."""
    conn = _connect_ro(db_path)
    try:
        rows = conn.execute(
            """
            SELECT v.gene_symbol,
                   COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position)
                       AS variant,
                   pd.pmid,
                   pd.total_carriers_observed,
                   pd.affected_count,
                   pd.unaffected_count,
                   pd.trust_reasons
            FROM penetrance_data pd
            LEFT JOIN variants v ON v.variant_id = pd.variant_id
            WHERE pd.trust_tier = 'quarantine'
            ORDER BY pd.total_carriers_observed IS NULL,
                     pd.total_carriers_observed DESC
            LIMIT ?
            """,
            (limit,),
        ).fetchall()
    finally:
        conn.close()
    return [
        {
            "gene": r[0],
            "variant": r[1],
            "pmid": r[2],
            "carriers": r[3],
            "affected": r[4],
            "unaffected": r[5],
            "reasons": json.loads(r[6]) if r[6] else [],
        }
        for r in rows
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True, help="Path to a gene DB.")
    parser.add_argument(
        "--list",
        type=int,
        default=0,
        metavar="N",
        help="Also list the top N quarantined facts.",
    )
    args = parser.parse_args()

    summary = summarize_trust(args.db)
    print(json.dumps(summary, indent=2))
    if not summary.get("tiered"):
        print(
            "\n(This DB predates the trust gate — no trust_tier column. Re-run the "
            "pipeline, or apply pipeline.trust_gate.apply_trust_gate to tier it.)"
        )
        return 0
    if summary["total"]:
        print(
            f"\nTrusted {summary['trusted']}/{summary['total']} "
            f"({100 * (1 - summary['quarantine_rate']):.1f}%); "
            f"quarantined {summary['quarantine']} "
            f"({100 * summary['quarantine_rate']:.1f}%)."
        )
    if args.list:
        print(f"\nTop {args.list} quarantined facts:")
        for row in list_quarantined(args.db, args.list):
            print(
                f"  {row['gene']} {row['variant']} PMID {row['pmid']}: "
                f"carriers={row['carriers']} → {', '.join(row['reasons'])}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
