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

**This pass writes values and provenance, never trust.** It runs at Step 3.45,
ahead of the guards and the trust gate, so they judge the adopted numbers. It
marks each field it filled with ``{field}_source = "figure_count_repair"``, and
``pipeline/trust_gate.py`` turns that into a ``figure_count_unverified:<field>``
quarantine of its own accord at Step 3.7.

An earlier version set ``trust_tier``/``field_trust``/``trust_reasons`` here
directly. That accomplished nothing — the gate rebuilds those columns from
scratch and preserves only ``paper_final_check:`` reasons, so the quarantine was
erased inside the same run and the adopted counts landed *trusted*, which is the
opposite of what the code appeared to say. Trust columns belong to the gate; a
step that lands a value states where it came from and lets the gate rule on it.
"""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

# The trust gate owns the reason vocabulary, and this source tag is the half of
# it that this module writes. Imported rather than re-spelled so the producer and
# the consumer cannot drift apart — which is exactly how the previous quarantine
# came to be silently discarded.
from pipeline.trust_gate import FIGURE_REPAIR_SOURCE

COUNT_COLUMNS = {
    "carriers": "total_carriers_observed",
    "affected": "affected_count",
    "unaffected": "unaffected_count",
}
RULES = ("adopt_figure_counts",)
COUNT_REPAIR_VERSION = "count-repair-v3"

#: Per-variant role stamped on each adopted field, in the vocabulary
#: ``pipeline.trust_gate`` reads. Kept in sync with its ``_RECOVERY_ROLE_BY_FIELD``.
FIGURE_COUNT_ROLES = {
    "carriers": "per_variant_carriers",
    "affected": "per_variant_affected",
    "unaffected": "per_variant_unaffected",
}


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


#: Columns the driving query needs, per table.
_REQUIRED_INPUTS = {
    "variant_papers": ("variant_id", "pmid", "additional_notes", "source_layer"),
    "penetrance_data": (
        "penetrance_id",
        "variant_id",
        "pmid",
        "total_carriers_observed",
        "affected_count",
        "unaffected_count",
    ),
}


def _missing_inputs(cur: sqlite3.Cursor) -> list[str]:
    """Which tables or columns the driving query needs and this DB lacks."""
    missing: list[str] = []
    for table, columns in _REQUIRED_INPUTS.items():
        present = {info[1] for info in cur.execute(f"PRAGMA table_info({table})")}
        if not present:
            missing.append(f"no {table} table")
            continue
        absent = [name for name in columns if name not in present]
        if absent:
            missing.append(f"{table} is missing {', '.join(absent)}")
    return missing


def _stamp_figure_provenance(
    cur: sqlite3.Cursor,
    variant_id: Any,
    pmid: Any,
    fields: set[str],
    columns: set[str],
) -> None:
    """Record which fields this pass filled, so the trust gate can rule on them.

    This pass writes values and provenance; it does not write trust. An earlier
    version set ``trust_tier``/``field_trust``/``trust_reasons`` directly, which
    accomplished nothing: the trust gate (Step 3.7) rebuilds those columns from
    scratch a few steps later and keeps only ``paper_final_check:`` reasons, so
    the quarantine, its reason, and its source tag were all erased inside the
    same run and the adopted counts landed trusted. Stamping provenance instead
    lets the gate re-derive ``figure_count_unverified`` itself, which also makes
    a standalone re-run of the gate idempotent.

    Marks only the fields actually adopted here, not every field the figure
    reader happened to see, so a count that some other layer supplied in the
    meantime is not mistaken for a vision read.
    """
    if "count_provenance" not in columns:
        return
    row = cur.execute(
        "SELECT count_provenance FROM variant_papers WHERE variant_id=? AND pmid=?",
        (variant_id, pmid),
    ).fetchone()
    provenance = _json_value(row[0] if row else None, dict, {})
    for field in sorted(fields):
        provenance[f"{field}_source"] = FIGURE_REPAIR_SOURCE
        provenance.setdefault(f"{field}_count_type", FIGURE_COUNT_ROLES[field])
    provenance["figure_count_repair_version"] = COUNT_REPAIR_VERSION
    cur.execute(
        "UPDATE variant_papers SET count_provenance=? WHERE variant_id=? AND pmid=?",
        (json.dumps(provenance, sort_keys=True), variant_id, pmid),
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
        except sqlite3.OperationalError as exc:
            # Name what is actually missing. A bare handler here reported "no
            # variant_papers table" for a missing *column* too, so a DB whose
            # variant_papers predates the `source_layer` migration silently
            # produced rows_examined=0 and the run looked like it had nothing to
            # repair rather than like it could not look.
            missing = _missing_inputs(cur)
            message = (
                f"count repair: {', '.join(missing)}; nothing to do"
                if missing
                else f"count repair: cannot read the count tables ({exc}); nothing to do"
            )
            if logger:
                logger.warning(message)
            summary["skipped_reason"] = message
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
        variant_paper_columns = {
            info[1] for info in cur.execute("PRAGMA table_info(variant_papers)")
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
                _stamp_figure_provenance(
                    cur,
                    row["variant_id"],
                    row["pmid"],
                    set(changes),
                    variant_paper_columns,
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
