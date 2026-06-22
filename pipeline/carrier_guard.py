"""Carrier-count sanity guard.

Some papers report a cohort size or a gnomAD allele-number (e.g. ~565k) that the
extractor captures as a per-variant ``total_carriers_observed`` — inflating
carrier aggregates by millions. This neutralizes observations above a plausible
per-variant ceiling: it copies the offending penetrance rows into a
``carrier_guard_quarantine`` table (reversible record) and NULLs their counts in
place, keeping the variant↔paper link but removing the bogus number.

Non-destructive (rows stay; only the counts are nulled) and idempotent.
"""

from __future__ import annotations

import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, Optional


def _ensure_column(cur: sqlite3.Cursor, table: str, name: str, decl: str) -> None:
    columns = {row[1] for row in cur.execute(f"PRAGMA table_info({table})")}
    if name not in columns:
        cur.execute(f"ALTER TABLE {table} ADD COLUMN {name} {decl}")


def apply_carrier_guard(
    db,
    *,
    threshold: int = 100_000,
    logger: Optional[Any] = None,
) -> dict:
    db = Path(db)
    con = sqlite3.connect(str(db))
    cur = con.cursor()
    try:
        n = cur.execute(
            "SELECT COUNT(*) FROM penetrance_data WHERE total_carriers_observed > ?",
            (threshold,),
        ).fetchone()[0]
    except sqlite3.OperationalError:
        con.close()
        if logger:
            logger.info("carrier guard: no penetrance_data table; nothing to do")
        return {"guarded": 0}
    if not n:
        con.close()
        if logger:
            logger.info("carrier guard: 0 observations over %d carriers", threshold)
        return {"guarded": 0}

    cur.execute(
        """CREATE TABLE IF NOT EXISTS carrier_guard_quarantine(
            penetrance_id INTEGER, variant_id INTEGER, pmid TEXT,
            total_carriers_observed INTEGER, affected_count INTEGER,
            unaffected_count INTEGER, uncertain_count INTEGER,
            penetrance_percentage REAL, threshold INTEGER, reason TEXT,
            guarded_at TEXT
        )"""
    )
    _ensure_column(cur, "carrier_guard_quarantine", "uncertain_count", "INTEGER")
    _ensure_column(cur, "carrier_guard_quarantine", "penetrance_percentage", "REAL")
    _ensure_column(cur, "carrier_guard_quarantine", "threshold", "INTEGER")
    _ensure_column(cur, "carrier_guard_quarantine", "reason", "TEXT")
    cur.execute(
        """INSERT INTO carrier_guard_quarantine (
            penetrance_id, variant_id, pmid, total_carriers_observed,
            affected_count, unaffected_count, uncertain_count,
            penetrance_percentage, threshold, reason, guarded_at
        )
        SELECT penetrance_id, variant_id, pmid, total_carriers_observed,
               affected_count, unaffected_count, uncertain_count,
               penetrance_percentage, ?,
               'total_carriers_observed exceeds plausible per-variant ceiling',
               ?
        FROM penetrance_data WHERE total_carriers_observed > ?""",
        (threshold, datetime.now().isoformat(timespec="seconds"), threshold),
    )
    cur.execute(
        """UPDATE penetrance_data
           SET total_carriers_observed = NULL,
               affected_count = NULL,
               unaffected_count = NULL,
               uncertain_count = NULL,
               penetrance_percentage = NULL
           WHERE total_carriers_observed > ?""",
        (threshold,),
    )
    con.commit()
    guarded_total = cur.execute(
        "SELECT COALESCE(SUM(total_carriers_observed), 0) FROM penetrance_data"
    ).fetchone()[0]
    con.close()
    if logger:
        logger.info(
            "carrier guard: neutralized %d observations over %d; guarded carrier total=%d",
            n,
            threshold,
            guarded_total,
        )
    return {
        "guarded": n,
        "threshold": threshold,
        "guarded_carrier_total": guarded_total,
    }
