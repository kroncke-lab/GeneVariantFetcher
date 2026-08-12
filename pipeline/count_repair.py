"""Deterministic per-variant count repair over a finalized run database.

Two defects in the stored counts are fixable with no model call, no network, and
no re-extraction, because the information is already on disk. Each was measured
on the locked Tier-1 baseline (48 cardiac papers, 1,000 gold rows) before being
written here; the per-fill accuracy is recorded against each rule.

Both rules only ever move a count toward what some part of the pipeline already
recorded. Neither computes a patient number that no lane asserted — on a gene
with no gold standard that is the difference between a resource a curator can
trust and one that quietly invents denominators.

``adopt_figure_counts``
    The figure reader is prompted for carriers/affected/unaffected and packs
    them into ``variant_papers.additional_notes`` as JSON, but nothing ever
    deserializes them: 80 figure rows carried counts in that blob and exactly 1
    reached ``penetrance_data``. Fills only null fields.
    *Measured: 90 of 91 fills exactly match gold (the miss is 19 vs 20 read off
    a Kaplan-Meier at-risk table).*

``refuse_all_unaffected``
    ``carriers = N, affected = 0, unaffected = N`` **on a row the trust gate has
    independently quarantined**. The shape alone is not a fabrication signature:
    it is also exactly how a control-cohort or benign-variant paper reports a
    real negative finding, and wiping those would delete the BS2-style evidence
    this resource exists to carry. So the rule never fires on its own judgement
    — it requires corroboration from a gate that reached its conclusion by a
    different route. Only the phenotype split is cleared; the carrier total is
    kept, because in this pattern it is the one number that is usually right.
    *Measured on Tier 1: the two rows with this shape separate perfectly.
    KCNQ1 33141630 T224M (124/0/124 against gold 124/34/54) is
    ``quarantine``/``population_count`` and is repaired, removing 104 units of
    absolute error while keeping the exact carrier total of 124. SCN5A 26746457
    p.Asp1243Asn (1/0/1, correct) is ``trusted`` with no reasons and is left
    alone.*

Order is load-bearing: figure counts are adopted first so a real reading is
available before any contradiction is judged.

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

Every change is recorded in ``count_repair_log`` with the previous value, so the
pass is reversible and auditable, and it is idempotent: re-running finds nothing
to do because each rule only writes into a null (or clears an exact pattern).
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
RULES = ("adopt_figure_counts", "refuse_all_unaffected")


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


def quarantined(trust_tier: Optional[str], trust_reasons: Optional[str]) -> bool:
    """Whether the trust gate independently flagged this row.

    Corroboration, not a second opinion from the same evidence: the gate reaches
    ``quarantine`` from study design and count plausibility, which is a
    different route than the arithmetic shape below.
    """
    del trust_reasons  # advisory reasons can sit on a trusted row; tier decides
    return (trust_tier or "").strip().lower() == "quarantine"


def refuse_all_unaffected(counts: dict, is_quarantined: bool = False) -> dict:
    """Clear a 100%-non-penetrant split that the trust gate also distrusts.

    Requires corroboration. ``carriers=N, affected=0, unaffected=N`` is equally
    the shape of a fabricated split and of a real control-cohort finding, and on
    a gene with no gold standard there is no way to tell them apart from the
    numbers alone — so the numbers alone are not allowed to decide.
    """
    if not is_quarantined:
        return {}
    carriers, affected, unaffected = (
        counts.get("carriers"),
        counts.get("affected"),
        counts.get("unaffected"),
    )
    if not _is_count(carriers) or carriers == 0:
        return {}
    if affected == 0 and unaffected == carriers:
        return {"affected": None, "unaffected": None}
    return {}


def repair_counts(
    counts: dict,
    notes: Optional[str],
    layer: Optional[str],
    *,
    trust_tier: Optional[str] = None,
    trust_reasons: Optional[str] = None,
) -> dict:
    """Apply every rule in order; return only the fields whose value changed."""
    working = dict(counts)
    changes: dict = {}
    flagged = quarantined(trust_tier, trust_reasons)
    # Each rule must see the previous rule's writes, so the deltas are computed
    # one at a time against the mutated `working` -- not eagerly in a literal.
    for rule, compute in (
        ("adopt_figure_counts", lambda w: adopt_figure_counts(w, notes, layer)),
        ("refuse_all_unaffected", lambda w: refuse_all_unaffected(w, flagged)),
    ):
        for field, value in compute(working).items():
            if working.get(field) == value:
                continue
            working[field] = value
            changes[field] = (value, rule)
    return changes


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
                              pd.trust_tier    AS trust_tier,
                              pd.trust_reasons AS trust_reasons,
                              vp.additional_notes AS additional_notes,
                              vp.source_layer     AS source_layer
                       FROM variant_papers vp
                       LEFT JOIN penetrance_data pd
                              ON pd.variant_id = vp.variant_id
                             AND pd.pmid = vp.pmid"""
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
        for row in rows:
            summary["rows_examined"] += 1
            counts = {f: row[f] for f in COUNT_COLUMNS}
            changes = {
                field: (value, rule)
                for field, (value, rule) in repair_counts(
                    counts,
                    row["additional_notes"],
                    row["source_layer"],
                    trust_tier=row["trust_tier"],
                    trust_reasons=row["trust_reasons"],
                ).items()
                if rule in enabled
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
        if dry_run:
            con.rollback()
        else:
            con.commit()
    finally:
        con.close()

    if logger:
        logger.info(
            "count repair%s: %d/%d rows changed (figure=%d refuse=%d)",
            " (dry run)" if dry_run else "",
            summary["rows_changed"],
            summary["rows_examined"],
            summary["adopt_figure_counts"],
            summary["refuse_all_unaffected"],
        )
    return summary
