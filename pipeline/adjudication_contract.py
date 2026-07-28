"""Shared adjudication-overlay contract: variant keying, verdicts, gold tiers.

``cli/compare_variants.py`` needs these on the scorer's hot path — three of its
imports were unguarded ``from scripts.ingest_review_adjudications import ...``
calls, which raise ``ImportError`` in an installed wheel because ``scripts`` is
excluded from the package. The scorer's whole adjudication-overlay path died
there. The definitions live here so both the packaged scorer and the
``scripts/ingest_review_adjudications`` CLI share exactly one implementation
rather than one importing the other across a packaging boundary.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Optional

#: Canonical cardiac channelopathy gene set backing the built-in gold tiers.
CARDIAC_GENES = frozenset({"KCNH2", "KCNQ1", "RYR2", "SCN5A"})

#: Reviewer verdict -> overlay action. The scorer translates a blank ``action``
#: through this SAME table the ingest writes, so a verdict newly added to ingest
#: resolves identically for scoring instead of being silently ignored.
VERDICT_TO_ACTION = {
    "confirm": "gold_confirmed",
    "correct_counts": "count_override",
    "wrong_variant": "false_positive",
    "wrong_paper": "excluded",
    "missing": "followup_missing",
    "other": "followup_other",
}

BUILTIN_GOLD_TIERS: dict[str, dict[str, Any]] = {
    "all": {
        "description": "All current lead-approved records, including non-cardiac genes.",
        "gene_mode": "all",
        "genes": [],
    },
    "cardiac": {
        "description": "Canonical cardiac channelopathy genes only.",
        "gene_mode": "include",
        "genes": sorted(CARDIAC_GENES),
    },
    "noncardiac": {
        "description": "Current records outside the canonical cardiac gene set.",
        "gene_mode": "exclude",
        "genes": sorted(CARDIAC_GENES),
    },
}


class GoldSyncError(RuntimeError):
    """The live gold contract could not be authenticated or validated."""


def variant_key(notation: str) -> str:
    """Canonicalize a protein notation the same way the scorer aggregates the DB.

    Mirrors ``aggregate_sqlite_data`` (canonical form, falling back to the
    normalized form) so an export row keys to the same bucket as the extracted
    (variant, paper) row would. The canonicalizers live in the scorer, which
    imports this module, so the import is deferred to call time to keep the
    dependency one-directional at module load.
    """
    from cli.compare_variants import normalize_variant, to_canonical_form

    canon = to_canonical_form(notation)
    return canon if canon else normalize_variant(notation)


def tier_definition_from_connection(
    conn: sqlite3.Connection, tier: str
) -> dict[str, Any]:
    name = str(tier or "").strip().lower()
    row = conn.execute(
        "SELECT name, description, gene_mode, genes_json, builtin "
        "FROM review_gold_tiers WHERE name = ?",
        (name,),
    ).fetchone()
    if row is None:
        raise GoldSyncError(f"Unknown review-gold tier: {tier!r}.")
    genes = json.loads(row[3])
    if not isinstance(genes, list):
        raise GoldSyncError(f"Review-gold tier {name!r} has invalid genes JSON.")
    return {
        "name": row[0],
        "description": row[1],
        "gene_mode": row[2],
        "genes": sorted({str(g).strip().upper() for g in genes if str(g).strip()}),
        "builtin": bool(row[4]),
    }


def tier_includes_gene(definition: dict[str, Any], gene: str) -> bool:
    normalized = str(gene or "").strip().upper()
    genes = set(definition["genes"])
    if definition["gene_mode"] == "all":
        return True
    if definition["gene_mode"] == "include":
        return normalized in genes
    return normalized not in genes


def read_gold_tier_definition(path: Optional[Path], tier: str) -> dict[str, Any]:
    """Resolve a built-in or cache-defined tier without mutating the cache."""
    name = str(tier or "").strip().lower()
    builtin = BUILTIN_GOLD_TIERS.get(name)
    if path is None or not Path(path).exists():
        if builtin is None:
            raise GoldSyncError(
                f"Review-gold tier {name!r} needs an existing cache definition."
            )
        return {"name": name, **builtin, "builtin": True}
    try:
        conn = sqlite3.connect(f"{Path(path).resolve().as_uri()}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        raise GoldSyncError(f"Could not open review-gold cache: {exc}") from exc
    try:
        try:
            return tier_definition_from_connection(conn, name)
        except sqlite3.Error:
            if builtin is not None:
                return {"name": name, **builtin, "builtin": True}
            raise GoldSyncError(f"Unknown review-gold tier: {name!r}.") from None
    finally:
        conn.close()


def gold_tier_includes_gene(path: Optional[Path], tier: str, gene: str) -> bool:
    return tier_includes_gene(read_gold_tier_definition(path, tier), gene)
