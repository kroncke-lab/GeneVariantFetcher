"""Deterministic per-variant count repair over a finalized run database.

One defect in the stored counts is fixable with no model call, no network, and
no re-extraction, because the information is already on disk. It was measured
on the locked Tier-1 baseline (48 cardiac papers, 1,000 gold rows).

``adopt_figure_counts``
    The figure reader is prompted for carriers/affected/unaffected and packs
    them into ``variant_papers.additional_notes`` as JSON, but nothing ever
    deserializes them: 80 figure rows carried counts in that blob and exactly 1
    reached ``penetrance_data``. Fills only null fields.
    *Measured: 90 of 91 fills exactly match gold (the miss is 19 vs 20 read off
    a Kaplan-Meier at-risk table).*

The former ``refuse_all_unaffected`` rule was retired after review. It cleared
raw affected/unaffected observations whenever *any* structural rule quarantined
the row. Quarantine is not proof that an all-unaffected observation is false,
and deleting the raw value violated the trust gate's audit-preservation
contract. Suspicious observations remain visible in raw data and are excluded
through the field-level trusted projection instead.

**Deliberately not implemented: ``unaffected = carriers - affected``.** It looks
free and it scored 59/61 on cardiac gold, but every one of those 61 derivations
wrote ``0`` — they were all rows where ``carriers == affected``, 45 of them with
a single carrier, and not one had a cohort of 10 or more. So the arithmetic was
never actually exercised; the rule was a "default unaffected to 0" in disguise,
which is the fabrication ``pipeline/prompts.py:68`` forbids in as many words. On
a population or biobank cohort — BMPR2, BRCA1/2, APOE, TTN, the genes with no
gold standard and therefore the real target — ``carriers - affected`` silently
reclassifies people who were never phenotyped as confirmed unaffected, inflating
apparent non-penetrance in the one direction that would mislead a curator.
``pipeline/claim_verifier.py:600`` already refuses this inference by name ("124
genotyped carriers but phenotype follow-up for only 88"). Closing that gap needs
an explicit unassessed/unphenotyped partition so that
``carriers = affected + unaffected + unassessed`` can be represented honestly,
not an arithmetic shortcut.

Every adopted value is recorded in ``count_repair_log`` with the previous value,
so the pass is reversible and auditable. It is idempotent because it only fills
null fields.
"""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

COUNT_COLUMNS = {
    "carriers": "total_carriers_observed",
    "affected": "affected_count",
    "unaffected": "unaffected_count",
}
RULES = ("adopt_figure_counts",)
COUNT_REPAIR_VERSION = "count-repair-v2"


def _is_count(value: Any) -> bool:
    """A usable count: a non-negative integer, and not a bool."""
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


# A figure number is only a per-variant patient count if the figure was
# reporting people, not a curve. The reader records what it was looking at in
# the blob's `context`, and these shapes are the ones whose integers are axis
# annotations, at-risk tallies, or plotted sample sizes rather than carriers.
# KCNQ1 18713323 is the worked example: a Kaplan-Meier at-risk table yielded
# 13/8/25/19, and the one figure fill that missed gold on Tier 1 (19 vs 20) came
# from exactly there.
NON_PATIENT_FIGURE_CONTEXT = (
    "kaplan",
    "at risk",
    "at-risk",
    "survival",
    "cumulative",
    "free survival",
    "hazard",
    "axis",
    "x-axis",
    "y-axis",
    "scale bar",
)


def figure_context_is_patient_level(context: Optional[str]) -> bool:
    """Whether a figure's own description denotes people rather than a curve."""
    lowered = (context or "").strip().lower()
    if not lowered:
        # No description is not evidence of a patient-level count. Fail closed:
        # on a gene with no gold standard nobody would catch the difference.
        return False
    return not any(marker in lowered for marker in NON_PATIENT_FIGURE_CONTEXT)


def figure_counts(additional_notes: Optional[str]) -> dict:
    """Counts the figure reader serialized into ``additional_notes``.

    Returns only the fields that carry a usable integer, and only when the
    reader's own recorded context says the figure was counting patients. A
    malformed or non-object blob yields nothing rather than raising: this runs
    over every row of a production database and a bad note must never fail the
    step.
    """
    if not additional_notes:
        return {}
    try:
        blob = json.loads(additional_notes)
    except (ValueError, TypeError):
        return {}
    if not isinstance(blob, dict):
        return {}
    if not figure_context_is_patient_level(blob.get("context")):
        return {}
    return {f: blob[f] for f in COUNT_COLUMNS if _is_count(blob.get(f))}


def adopt_figure_counts(
    counts: dict, notes: Optional[str], layer: Optional[str]
) -> dict:
    """Fill null counts from the figure reader's own serialized reading."""
    if "figure" not in (layer or ""):
        return {}
    return {
        f: value for f, value in figure_counts(notes).items() if counts.get(f) is None
    }


def repair_counts(
    counts: dict,
    notes: Optional[str],
    layer: Optional[str],
    *,
    rules: Optional[set[str]] = None,
) -> dict:
    """Apply enabled rules; return only fields whose stored value should change."""
    enabled = set(RULES) if rules is None else {rule for rule in rules if rule in RULES}
    if "adopt_figure_counts" not in enabled:
        return {}
    return {
        field: (value, "adopt_figure_counts")
        for field, value in adopt_figure_counts(counts, notes, layer).items()
    }


def _json_value(raw: Any, expected_type: type, default: Any) -> Any:
    try:
        value = json.loads(raw) if raw else default
    except (TypeError, ValueError, json.JSONDecodeError):
        return default
    return value if isinstance(value, expected_type) else default


def _quarantine_figure_fields(
    cur: sqlite3.Cursor,
    penetrance_id: int,
    fields: set[str],
    columns: set[str],
    *,
    row_created: bool,
) -> None:
    """Keep newly adopted figure fields pending until the trust gate evaluates them."""
    available = {
        name
        for name in (
            "trust_tier",
            "trust_reasons",
            "trust_rule_version",
            "field_trust",
            "trust_sources",
        )
        if name in columns
    }
    field_scoped = "field_trust" in available
    if not available or (not row_created and not field_scoped):
        # Legacy rows without a field-level mask may already contain trusted
        # text counts. Quarantining the whole row would hide those unrelated
        # observations, so only a brand-new row may use the row-level fallback.
        return
    row = cur.execute(
        f"SELECT {', '.join(sorted(available))} FROM penetrance_data "
        "WHERE penetrance_id=?",
        (penetrance_id,),
    ).fetchone()
    current = dict(zip(sorted(available), row)) if row else {}
    assignments: list[str] = []
    values: list[Any] = []

    def assign(column: str, value: Any) -> None:
        if column in available:
            assignments.append(f"{column}=?")
            values.append(value)

    assign("trust_tier", "quarantine")
    reasons = _json_value(current.get("trust_reasons"), list, [])
    assign(
        "trust_reasons",
        json.dumps(
            sorted(
                {str(reason) for reason in reasons} | {"figure_count_pending_review"}
            )
        ),
    )
    field_trust = _json_value(current.get("field_trust"), dict, {})
    trust_fields = {
        "total_carriers" if field == "carriers" else field for field in fields
    }
    field_trust.update({field: "quarantine" for field in trust_fields})
    assign("field_trust", json.dumps(field_trust, sort_keys=True))
    sources = _json_value(current.get("trust_sources"), list, [])
    assign(
        "trust_sources",
        json.dumps(
            sorted({str(source) for source in sources} | {"figure_count_repair"})
        ),
    )
    assign("trust_rule_version", COUNT_REPAIR_VERSION)
    if assignments:
        cur.execute(
            f"UPDATE penetrance_data SET {', '.join(assignments)} "
            "WHERE penetrance_id=?",
            (*values, penetrance_id),
        )


def apply_count_repair(
    db,
    *,
    rules: Optional[set] = None,
    dry_run: bool = False,
    logger: Optional[Any] = None,
) -> dict:
    """Repair counts across a run database. Returns a per-rule change summary."""
    db = Path(db)
    enabled = set(RULES) if rules is None else {r for r in rules if r in RULES}
    summary = {
        "rows_examined": 0,
        "rows_changed": 0,
        "rows_created": 0,
        "ambiguous_variant_papers": 0,
        "dry_run": dry_run,
    }
    summary.update({rule: 0 for rule in RULES})
    if not enabled:
        return summary

    con = sqlite3.connect(str(db))
    try:
        con.row_factory = sqlite3.Row
        cur = con.cursor()
        # Driven from variant_papers, not penetrance_data. A penetrance row only
        # exists once some count was already non-null (harvesting/
        # migrate_to_sqlite.py), so the figure rows this pass exists to rescue
        # have no row to UPDATE — they need one inserted.
        try:
            rows = list(
                cur.execute(
                    """SELECT pd.penetrance_id AS penetrance_id,
                              vp.variant_id    AS variant_id,
                              vp.pmid          AS pmid,
                              pd.total_carriers_observed AS carriers,
                              pd.affected_count          AS affected,
                              pd.unaffected_count        AS unaffected,
                              vp.additional_notes AS additional_notes,
                              vp.source_layer     AS source_layer,
                              (SELECT COUNT(*) FROM penetrance_data siblings
                               WHERE siblings.variant_id = vp.variant_id
                                 AND siblings.pmid = vp.pmid) AS penetrance_rows
                       FROM variant_papers vp
                       LEFT JOIN penetrance_data pd
                              ON pd.penetrance_id = (
                                  SELECT MIN(parent.penetrance_id)
                                  FROM penetrance_data parent
                                  WHERE parent.variant_id = vp.variant_id
                                    AND parent.pmid = vp.pmid
                              )"""
                )
            )
        except sqlite3.OperationalError:
            if logger:
                logger.info("count repair: no variant_papers table; nothing to do")
            return summary

        if not dry_run:
            cur.execute(
                """CREATE TABLE IF NOT EXISTS count_repair_log(
                    penetrance_id INTEGER, variant_id INTEGER, pmid TEXT,
                    field TEXT, previous_value INTEGER, new_value INTEGER,
                    rule TEXT, repaired_at TEXT
                )"""
            )

        stamp = datetime.now(timezone.utc).isoformat()
        penetrance_columns = {
            info[1] for info in cur.execute("PRAGMA table_info(penetrance_data)")
        }
        for row in rows:
            summary["rows_examined"] += 1
            if row["penetrance_rows"] > 1:
                # A single figure observation cannot be assigned safely to an
                # unspecified cohort. Leave every sibling untouched.
                summary["ambiguous_variant_papers"] += 1
                continue
            counts = {f: row[f] for f in COUNT_COLUMNS}
            changes = {
                field: (value, rule)
                for field, (value, rule) in repair_counts(
                    counts,
                    row["additional_notes"],
                    row["source_layer"],
                    rules=enabled,
                ).items()
            }
            if not changes:
                continue
            summary["rows_changed"] += 1
            penetrance_id = row["penetrance_id"]
            row_created = penetrance_id is None
            if penetrance_id is None and not dry_run:
                cur.execute(
                    "INSERT INTO penetrance_data(variant_id, pmid) VALUES(?,?)",
                    (row["variant_id"], row["pmid"]),
                )
                penetrance_id = cur.lastrowid
                summary["rows_created"] += 1
            elif penetrance_id is None:
                summary["rows_created"] += 1
            for field, (value, rule) in changes.items():
                summary[rule] += 1
                if dry_run:
                    continue
                cur.execute(
                    f"UPDATE penetrance_data SET {COUNT_COLUMNS[field]} = ? "
                    "WHERE penetrance_id = ?",
                    (value, penetrance_id),
                )
                cur.execute(
                    "INSERT INTO count_repair_log(penetrance_id, variant_id, pmid, "
                    "field, previous_value, new_value, rule, repaired_at) "
                    "VALUES(?,?,?,?,?,?,?,?)",
                    (
                        penetrance_id,
                        row["variant_id"],
                        row["pmid"],
                        field,
                        counts[field],
                        value,
                        rule,
                        stamp,
                    ),
                )
            if not dry_run:
                _quarantine_figure_fields(
                    cur,
                    penetrance_id,
                    set(changes),
                    penetrance_columns,
                    row_created=row_created,
                )
        if dry_run:
            con.rollback()
        else:
            con.commit()
    finally:
        con.close()

    if logger:
        logger.info(
            "count repair%s: %d/%d rows changed (figure=%d)",
            " (dry run)" if dry_run else "",
            summary["rows_changed"],
            summary["rows_examined"],
            summary["adopt_figure_counts"],
        )
    return summary
