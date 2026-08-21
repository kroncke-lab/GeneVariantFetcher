"""Persistent paper-level scope exclusions shared by enrichment layers.

Extraction can determine that an entire paper is outside the human target-gene
scope before any variant is accepted.  Recovery layers must honor that decision:
otherwise text-mined or citation-derived variants can silently repopulate a paper
that extraction deliberately rejected.
"""

from __future__ import annotations

import sqlite3


NONHUMAN_ORTHOLOG_MODEL = "skipped-nonhuman-ortholog"
NONHUMAN_ORTHOLOG_REASON = "nonhuman_target_gene_ortholog"


def metadata_paper_scope_exclusion_reason(metadata: object) -> str | None:
    """Read the paper-scope decision from extraction JSON metadata."""

    if not isinstance(metadata, dict):
        return None
    explicit = str(metadata.get("paper_scope_exclusion_reason") or "").strip()
    if explicit:
        return explicit
    if (
        metadata.get("model_used") == NONHUMAN_ORTHOLOG_MODEL
        or metadata.get("skipped_nonhuman_ortholog")
        or metadata.get("dropped_nonhuman_ortholog")
    ):
        return NONHUMAN_ORTHOLOG_REASON
    return None


def db_paper_scope_exclusions(con: sqlite3.Connection) -> dict[str, str]:
    """Return ``{pmid: reason}`` for durable paper-level DB exclusions.

    Older/minimal DBs may lack ``extraction_metadata`` or ``model_used``.  They
    remain readable and simply provide no durable exclusions.
    """

    try:
        tables = {
            str(row[0])
            for row in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        if "extraction_metadata" not in tables:
            return {}
        columns = {
            str(row[1])
            for row in con.execute("PRAGMA table_info(extraction_metadata)").fetchall()
        }
        if "pmid" not in columns:
            return {}
        predicates: list[str] = []
        params: list[str] = []
        if "model_used" in columns:
            predicates.append("model_used=?")
            params.append(NONHUMAN_ORTHOLOG_MODEL)
        if "paper_scope_exclusion_reason" in columns:
            predicates.append(
                "paper_scope_exclusion_reason IS NOT NULL "
                "AND TRIM(paper_scope_exclusion_reason) != ''"
            )
        if not predicates:
            return {}
        rows = con.execute(
            f"""SELECT DISTINCT pmid FROM extraction_metadata
                WHERE ({" OR ".join(predicates)})
                  AND pmid IS NOT NULL AND TRIM(pmid) != ''""",
            tuple(params),
        ).fetchall()
    except sqlite3.DatabaseError:
        return {}
    return {
        str(row[0]): NONHUMAN_ORTHOLOG_REASON for row in rows if str(row[0]).isdigit()
    }


def paper_scope_exclusion_reason(con: sqlite3.Connection, pmid: str) -> str | None:
    """Return a stable exclusion reason for *pmid*, if one is persisted."""

    return db_paper_scope_exclusions(con).get(str(pmid))


def purge_excluded_paper_evidence(con: sqlite3.Connection) -> int:
    """Delete clinical evidence rows for papers already marked out of scope.

    This repairs DBs contaminated by a recovery layer before the shared scope
    contract existed.  The paper and extraction metadata remain as an auditable
    rejected attempt; evidence that could reach a clinical browser is removed.
    """

    excluded = db_paper_scope_exclusions(con)
    if not excluded:
        return 0
    placeholders = ",".join("?" for _ in excluded)
    params = tuple(sorted(excluded))
    evidence_tables = (
        "count_repair_log",
        "fact_provenance",
        "functional_data",
        "individual_records",
        "penetrance_data",
        "phenotypes",
        "tables_processed",
        "variant_metadata",
        "variant_papers",
    )
    removed = 0
    for table in evidence_tables:
        try:
            columns = {
                str(row[1])
                for row in con.execute(f"PRAGMA table_info({table})").fetchall()
            }
            if "pmid" not in columns:
                continue
            cursor = con.execute(
                f"DELETE FROM {table} WHERE pmid IN ({placeholders})", params
            )
            removed += max(int(cursor.rowcount), 0)
        except sqlite3.DatabaseError:
            continue
    return removed
