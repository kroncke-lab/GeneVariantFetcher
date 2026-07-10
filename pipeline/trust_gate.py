"""Per-fact confidence/trust gate → two-tier (trusted / quarantine) evidence.

The autonomy target (see ``docs/AUTONOMY_ROADMAP.md``) is unattended operation
across many genes with automated quality gates as the primary control. This
module is the keystone: it decides, **per extracted count fact and without any
gold standard**, whether to trust it. Facts that fail a gold-free structural
check are **soft-quarantined** — ``trust_tier='quarantine'`` with a list of
reason codes — never deleted and never NULLed, so the evidence stays auditable
and can be re-tiered after a rule change.

The rules are deliberately **gene-class-agnostic** so they transfer from cardiac
missense to BRCA truncating/case-control without per-gene tuning:

- ``arith_inconsistent``  affected + unaffected + uncertain exceeds the total
  carrier count (the parts cannot exceed the whole).
- ``count_is_total``      the extracted "carrier" count is really a cohort
  denominator (count role is a total/screened-N), not per-variant carriers.
- ``population_count``     the count is a population allele magnitude (gnomAD /
  MAF / allele-frequency label, or a carrier count above a population ceiling).
- ``paper_outlier``       the carrier count is a wild outlier vs the other
  variants in the same paper (median x k).

The pure core (``evaluate_fact`` / ``decide_tier``) has no DB dependency and is
unit-tested directly; ``apply_trust_gate`` is the thin SQLite adapter.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import sqlite3
import statistics
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

# Count roles (from pipeline/prompts.py) that mean the number is a denominator /
# aggregate, not a per-variant carrier count.
DENOMINATOR_COUNT_TYPES = frozenset(
    {"cohort_total", "screened_n", "study_total", "total"}
)

# Column labels that signal a population allele count/frequency rather than a
# clinical carrier count (kept in sync in spirit with
# pipeline/count_classifier._SUSPICIOUS_LABEL_PATTERNS).
POPULATION_LABEL_RE = re.compile(
    r"gnomad|exac|topmed|1000\s*genomes|allele\s*freq|allele\s*count|"
    r"\bmaf\b|\ba[cn]\b|minor\s*allele|population\s*frequency",
    re.IGNORECASE,
)

# A per-variant clinical carrier count above this is almost certainly a
# population allele number (gnomAD AN-scale), not a cohort of carriers. Used as a
# REASON, not an in-place NULL (that is carrier_guard's legacy behavior).
POPULATION_CARRIER_CEILING = 100_000

# Within-paper outlier: a carrier count more than K x the paper's median and
# above an absolute floor (mirrors pipeline/count_outlier_guard defaults).
OUTLIER_K = 10
OUTLIER_ABS_FLOOR = 50

# Study designs where a clinical per-variant carrier/penetrance count is a trap
# (review, pure functional, GWAS allele tables). Soft-quarantine via
# study_type_mismatch; never delete the row.
STUDY_TYPE_MISMATCH_DESIGNS = frozenset({"review_meta", "functional_invitro", "gwas"})
# Population/biobank designs: strengthen population_count for large carriers.
POPULATION_STUDY_DESIGNS = frozenset({"cohort_population", "cohort_biobank", "gwas"})
# Soft floor for study-context population_count (below the hard gnomAD ceiling).
POPULATION_STUDY_CARRIER_FLOOR = 50

# Rule identifiers, ordered. The rule-set version below is derived from this so a
# stored fact records exactly which rule generation tiered it.
RULE_IDS = (
    "arith_inconsistent",
    "count_is_total",
    "population_count",
    "paper_outlier",
    "study_type_mismatch",
)


def rule_version() -> str:
    """Short content hash of the rule set + thresholds, so a stored fact's
    ``trust_rule_version`` pins the exact generation of rules that tiered it."""
    payload = json.dumps(
        {
            "rules": RULE_IDS,
            "denominator_types": sorted(DENOMINATOR_COUNT_TYPES),
            "population_ceiling": POPULATION_CARRIER_CEILING,
            "outlier_k": OUTLIER_K,
            "outlier_floor": OUTLIER_ABS_FLOOR,
            "population_label_re": POPULATION_LABEL_RE.pattern,
            "study_type_mismatch_designs": sorted(STUDY_TYPE_MISMATCH_DESIGNS),
            "population_study_designs": sorted(POPULATION_STUDY_DESIGNS),
            "population_study_carrier_floor": POPULATION_STUDY_CARRIER_FLOOR,
        },
        sort_keys=True,
    )
    # tg2: study-record aware generation (study_type_mismatch + population study).
    return "tg2-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _as_int(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def evaluate_fact(
    counts: dict[str, Any],
    provenance: Optional[dict[str, Any]] = None,
    paper_stats: Optional[dict[str, Any]] = None,
    study_context: Optional[dict[str, Any]] = None,
) -> list[str]:
    """Return the sorted reason codes a single count fact trips (empty = clean).

    ``counts``: ``carriers`` / ``affected`` / ``unaffected`` / ``uncertain``
    (any may be None). ``provenance``: the fact's ``count_provenance`` dict
    (``carriers_count_type``, ``carriers_column_label``). ``paper_stats``:
    ``{"carriers_median": float, "carriers_n": int}`` for the fact's paper.
    ``study_context``: optional paper-level fields from extraction_metadata
    (``study_design``, ``ascertainment``, ``study_type``).
    """
    provenance = provenance or {}
    paper_stats = paper_stats or {}
    study_context = study_context or {}
    reasons: set[str] = set()

    carriers = _as_int(counts.get("carriers"))
    affected = _as_int(counts.get("affected"))
    unaffected = _as_int(counts.get("unaffected"))
    uncertain = _as_int(counts.get("uncertain"))

    # arith_inconsistent: the phenotype breakdown cannot exceed the total.
    if carriers is not None and (affected is not None or unaffected is not None):
        parts = (affected or 0) + (unaffected or 0) + (uncertain or 0)
        if parts > carriers:
            reasons.add("arith_inconsistent")

    # count_is_total: the carrier number is really a cohort denominator.
    ctype = str(provenance.get("carriers_count_type") or "").strip().lower()
    if ctype in DENOMINATOR_COUNT_TYPES:
        reasons.add("count_is_total")

    # population_count: population allele magnitude, by label or by ceiling.
    label = str(provenance.get("carriers_column_label") or "")
    if label and POPULATION_LABEL_RE.search(label):
        reasons.add("population_count")
    if carriers is not None and carriers > POPULATION_CARRIER_CEILING:
        reasons.add("population_count")

    design = str(study_context.get("study_design") or "").strip().lower()
    # Strengthen population_count when the paper is itself a population/biobank/GWAS study.
    if (
        design in POPULATION_STUDY_DESIGNS
        and carriers is not None
        and carriers >= POPULATION_STUDY_CARRIER_FLOOR
    ):
        reasons.add("population_count")

    # study_type_mismatch: clinical carrier/penetrance counts from review /
    # functional / GWAS papers are soft-quarantined (still auditable).
    has_clinical_count = any(
        v is not None and v > 0 for v in (carriers, affected, unaffected)
    )
    if design in STUDY_TYPE_MISMATCH_DESIGNS and has_clinical_count:
        reasons.add("study_type_mismatch")

    # paper_outlier: carrier count is a wild outlier vs the paper's other rows.
    median = paper_stats.get("carriers_median")
    n = paper_stats.get("carriers_n") or 0
    if (
        carriers is not None
        and median is not None
        and n >= 3
        and median > 0
        and carriers > median * OUTLIER_K
        and carriers > OUTLIER_ABS_FLOOR
    ):
        reasons.add("paper_outlier")

    return sorted(reasons)


def decide_tier(reasons: list[str]) -> str:
    """Map reason codes to a tier. Soft: any reason quarantines (keeps the row);
    a clean fact stays trusted. (Severity weighting can come later.)"""
    return "quarantine" if reasons else "trusted"


def _paper_carrier_stats(
    rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Per-PMID median + count of non-null carrier values, for the outlier rule."""
    by_pmid: dict[str, list[int]] = {}
    for r in rows:
        c = _as_int(r.get("carriers"))
        if c is not None:
            by_pmid.setdefault(str(r.get("pmid")), []).append(c)
    stats: dict[str, dict[str, Any]] = {}
    for pmid, values in by_pmid.items():
        stats[pmid] = {
            "carriers_median": statistics.median(values) if values else None,
            "carriers_n": len(values),
        }
    return stats


def apply_trust_gate(db_path: str | Path) -> dict[str, Any]:
    """Tier every ``penetrance_data`` fact in the DB. Soft-quarantine only: sets
    ``trust_tier`` / ``trust_reasons`` / ``trust_rule_version``; never NULLs a
    count or deletes a row. Idempotent — re-running re-tiers against current
    rules. Returns ``{trusted, quarantine, by_reason, rule_version}``.
    """
    db_path = str(db_path)
    version = rule_version()
    conn = sqlite3.connect(db_path)
    try:
        conn.row_factory = sqlite3.Row
        # Pull each count fact with its fact-level count provenance (variant ↔
        # paper carries the count_provenance JSON).
        # Optional study columns (PR-1 A2); tolerate pre-study-record DBs.
        em_cols = {r[1] for r in conn.execute("PRAGMA table_info(extraction_metadata)")}
        study_select = []
        for col in (
            "study_design",
            "ascertainment",
            "study_type",
            "cohort_source",
            "population",
        ):
            if col in em_cols:
                study_select.append(f"em.{col} AS {col}")
            else:
                study_select.append(f"NULL AS {col}")
        study_sql = ",\n                       ".join(study_select)

        rows = [
            dict(r)
            for r in conn.execute(
                f"""
                SELECT pd.penetrance_id       AS penetrance_id,
                       pd.pmid                 AS pmid,
                       pd.total_carriers_observed AS carriers,
                       pd.affected_count       AS affected,
                       pd.unaffected_count     AS unaffected,
                       pd.uncertain_count      AS uncertain,
                       vp.count_provenance     AS count_provenance,
                       {study_sql}
                FROM penetrance_data pd
                LEFT JOIN variant_papers vp
                       ON vp.variant_id = pd.variant_id AND vp.pmid = pd.pmid
                LEFT JOIN extraction_metadata em
                       ON em.pmid = pd.pmid
                """
            )
        ]

        paper_stats = _paper_carrier_stats(rows)
        stats: dict[str, Any] = {
            "trusted": 0,
            "quarantine": 0,
            "by_reason": {},
            "rule_version": version,
        }
        updates: list[tuple[str, str, str, int]] = []
        for r in rows:
            provenance = {}
            raw = r.get("count_provenance")
            if raw:
                try:
                    provenance = json.loads(raw)
                except (TypeError, ValueError, json.JSONDecodeError):
                    provenance = {}
            study_context = {
                "study_design": r.get("study_design"),
                "ascertainment": r.get("ascertainment"),
                "study_type": r.get("study_type"),
            }
            reasons = evaluate_fact(
                r,
                provenance,
                paper_stats.get(str(r.get("pmid"))),
                study_context=study_context,
            )
            tier = decide_tier(reasons)
            stats[tier] += 1
            for reason in reasons:
                stats["by_reason"][reason] = stats["by_reason"].get(reason, 0) + 1
            updates.append((tier, json.dumps(reasons), version, r["penetrance_id"]))

        conn.executemany(
            "UPDATE penetrance_data SET trust_tier=?, trust_reasons=?, "
            "trust_rule_version=? WHERE penetrance_id=?",
            updates,
        )
        conn.commit()
        logger.info(
            "trust gate: %d trusted, %d quarantined (%s) [%s]",
            stats["trusted"],
            stats["quarantine"],
            stats["by_reason"],
            version,
        )
        return stats
    finally:
        conn.close()
