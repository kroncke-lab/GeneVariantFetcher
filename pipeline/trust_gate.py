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

- ``arith_inconsistent``  the phenotype breakdown is impossible: parts exceed
  the total, a fully-specified partition does not equal the total, or any count
  is negative.
- ``count_is_total``      the extracted "carrier" count is really a cohort
  denominator (count role is a total/screened-N), not per-variant carriers.
- ``population_count``     the count is a population allele magnitude (gnomAD /
  MAF / allele-frequency label, or a carrier count above a population ceiling).
- ``paper_outlier``       the carrier count is a wild outlier vs the other
  variants in the same paper (median x k).
- ``negative_count``      any carrier/affected/unaffected/uncertain count is
  negative (never valid); folded into ``arith_inconsistent`` semantics but kept
  as a distinct reason for auditability.
- ``implied_unaffected_zero``  a *derived* ``unaffected == 0`` (affected equals
  the carrier total and no unaffected column/quote backs it) in a study whose
  design deliberately includes potentially-unaffected carriers (population /
  biobank screening, case-control, family cascade/segregation). This is an
  unsupported 100%-penetrance assertion, so **only the unaffected field** is
  soft-quarantined — affected/total stay trusted. Proband case reports/series
  and unknown-design papers never trip it (the cardiac default is preserved).

The pure core (``evaluate_fact`` / ``decide_tier``) has no DB dependency and is
unit-tested directly; ``apply_trust_gate`` is the thin SQLite adapter. The
adapter preserves any already-composed ``paper_final_check:`` field findings;
the paper gate still follows this step in ``gvf-run`` to validate freshness and
perform a full deterministic recomposition.
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

# Designs / ascertainments that deliberately enroll carriers who may be
# UNAFFECTED (population/biobank screening, case-control, family cascade or
# segregation). In these, a *derived* unaffected=0 (all carriers labelled
# affected, no unaffected column/quote) is an unsupported 100%-penetrance claim.
# The set is deliberately AFFIRMATIVE: an unknown/null design (the cardiac legacy
# default) is NOT a member, so the rule stays dormant there and only bites the
# cohort/screening/cascade studies (common in BRCA) where it is actually wrong.
UNAFFECTED_EXPECTED_DESIGNS = frozenset(
    {"cohort_population", "cohort_biobank", "case_control", "family_segregation"}
)
UNAFFECTED_EXPECTED_ASCERTAINMENTS = frozenset(
    {"population_screening", "biobank", "family_cascade"}
)
# Count-provenance values that mean the unaffected number was NOT read from a
# real unaffected column (so a 0 there is inferred, not sourced).
_UNSOURCED_UNAFFECTED_TYPES = frozenset({"", "unknown"})

# Rule identifiers, ordered. The rule-set version below is derived from this so a
# stored fact records exactly which rule generation tiered it.
RULE_IDS = (
    "arith_inconsistent",
    "count_is_total",
    "population_count",
    "paper_outlier",
    "study_type_mismatch",
    "negative_count",
    "implied_unaffected_zero",
)

COUNT_TRUST_FIELDS = ("total_carriers", "affected", "unaffected", "uncertain")

# Which count fields a reason masks from the trusted projection. A reason not
# listed here masks the WHOLE fact (all count fields) — the historical
# all-or-nothing behavior. Field-specific reasons let us quarantine only the
# unsupported field (e.g. a bad derived unaffected) while keeping the
# well-supported affected/total counts trusted.
REASON_FIELDS: dict[str, tuple[str, ...]] = {
    "implied_unaffected_zero": ("unaffected",),
}


def _fields_masked_by(reason: str) -> tuple[str, ...]:
    """Count fields a single reason quarantines (defaults to the whole fact)."""
    return REASON_FIELDS.get(reason, COUNT_TRUST_FIELDS)


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
            "unaffected_expected_designs": sorted(UNAFFECTED_EXPECTED_DESIGNS),
            "unaffected_expected_ascertainments": sorted(
                UNAFFECTED_EXPECTED_ASCERTAINMENTS
            ),
            "reason_fields": {k: sorted(v) for k, v in sorted(REASON_FIELDS.items())},
        },
        sort_keys=True,
    )
    # tg3: negative_count + full-partition arith + field-scoped implied_unaffected_zero.
    return "tg3-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


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

    # negative_count: a count can never be negative. Kept distinct for audit but
    # also implies arith_inconsistent.
    if any(
        v is not None and v < 0 for v in (carriers, affected, unaffected, uncertain)
    ):
        reasons.add("negative_count")
        reasons.add("arith_inconsistent")

    # arith_inconsistent: the phenotype breakdown cannot exceed the total, and a
    # fully-specified partition (affected + unaffected + uncertain all reported)
    # must equal the total. A partial partition summing to LESS than the total is
    # left alone (the residual is simply unphenotyped carriers).
    if carriers is not None and (affected is not None or unaffected is not None):
        parts = (affected or 0) + (unaffected or 0) + (uncertain or 0)
        full_partition = (
            affected is not None and unaffected is not None and uncertain is not None
        )
        if parts > carriers or (full_partition and parts != carriers):
            reasons.add("arith_inconsistent")

    # implied_unaffected_zero: a *derived* unaffected==0 (every carrier labelled
    # affected, nothing sourced the zero) in a study design that deliberately
    # enrolls potentially-unaffected carriers is an unsupported 100%-penetrance
    # claim. Masks only the unaffected field; affected/total stay trusted.
    design = str(study_context.get("study_design") or "").strip().lower()
    ascertainment = str(study_context.get("ascertainment") or "").strip().lower()
    if (
        unaffected == 0
        and carriers is not None
        and carriers > 0
        and affected is not None
        and affected == carriers
    ):
        unaff_label = str(provenance.get("unaffected_column_label") or "").strip()
        unaff_type = str(provenance.get("unaffected_count_type") or "").strip().lower()
        unaffected_sourced = bool(unaff_label) or (
            unaff_type not in _UNSOURCED_UNAFFECTED_TYPES
        )
        design_expects_unaffected = (
            design in UNAFFECTED_EXPECTED_DESIGNS
            or ascertainment in UNAFFECTED_EXPECTED_ASCERTAINMENTS
        )
        if not unaffected_sourced and design_expects_unaffected:
            reasons.add("implied_unaffected_zero")

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
    structural rules while retaining already-composed paper-level field masks.
    The paper gate must still run afterward to validate their payload hash.
    Returns ``{trusted, quarantine, by_reason, rule_version}``.
    """
    db_path = str(db_path)
    version = rule_version()
    conn = sqlite3.connect(db_path)
    try:
        conn.row_factory = sqlite3.Row
        pd_cols = {r[1] for r in conn.execute("PRAGMA table_info(penetrance_data)")}
        for name, declaration in (
            ("trust_tier", "TEXT DEFAULT 'trusted'"),
            ("trust_reasons", "TEXT"),
            ("trust_rule_version", "TEXT"),
            ("field_trust", "TEXT"),
            ("trust_sources", "TEXT"),
        ):
            if name not in pd_cols:
                conn.execute(
                    f"ALTER TABLE penetrance_data ADD COLUMN {name} {declaration}"
                )
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
                       pd.trust_tier            AS existing_trust_tier,
                       pd.trust_reasons         AS existing_trust_reasons,
                       pd.trust_rule_version    AS existing_trust_rule_version,
                       pd.field_trust           AS existing_field_trust,
                       vp.count_provenance     AS count_provenance,
                       {study_sql}
                FROM penetrance_data pd
                LEFT JOIN variant_papers vp
                       ON vp.variant_id = pd.variant_id AND vp.pmid = pd.pmid
                LEFT JOIN extraction_metadata em
                       ON em.pmid = pd.pmid
                      AND em.extraction_id = (
                              SELECT MAX(em2.extraction_id)
                              FROM extraction_metadata em2
                              WHERE em2.pmid = pd.pmid
                          )
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
        updates: list[tuple[str, str, str, str, str, int]] = []
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
            existing_reasons: list[str] = []
            try:
                parsed_reasons = json.loads(r.get("existing_trust_reasons") or "[]")
                if isinstance(parsed_reasons, list):
                    existing_reasons = [str(value) for value in parsed_reasons]
            except (TypeError, ValueError, json.JSONDecodeError):
                pass
            paper_reasons = [
                reason
                for reason in existing_reasons
                if reason.startswith("paper_final_check:")
            ]
            combined_reasons = sorted(set(reasons) | set(paper_reasons))

            # Per-reason field masking: most reasons quarantine the whole fact
            # (all count fields); field-scoped reasons (e.g.
            # implied_unaffected_zero) quarantine only the unsupported field so a
            # well-sourced affected/total count on the same row stays trusted.
            field_trust = {field: "trusted" for field in COUNT_TRUST_FIELDS}
            for reason in reasons:
                for field in _fields_masked_by(reason):
                    field_trust[field] = "quarantine"
            if paper_reasons:
                try:
                    existing_fields = json.loads(r.get("existing_field_trust") or "{}")
                except (TypeError, ValueError, json.JSONDecodeError):
                    existing_fields = {}
                if isinstance(existing_fields, dict):
                    for field in COUNT_TRUST_FIELDS:
                        if existing_fields.get(field) == "quarantine":
                            field_trust[field] = "quarantine"
                # Recover the mask from reason metadata if a legacy row lacks a
                # valid field_trust object.
                for reason in paper_reasons:
                    field = reason.rsplit(":", 1)[-1]
                    if field in COUNT_TRUST_FIELDS:
                        field_trust[field] = "quarantine"

            tier = (
                "quarantine"
                if any(value == "quarantine" for value in field_trust.values())
                else "trusted"
            )
            stats[tier] += 1
            for reason in combined_reasons:
                stats["by_reason"][reason] = stats["by_reason"].get(reason, 0) + 1
            sources = []
            if reasons:
                sources.append("structural")
            if paper_reasons:
                sources.append("paper_final_check")
            prior_composer_versions = [
                part
                for part in str(r.get("existing_trust_rule_version") or "").split("|")
                if part.startswith("pfcg")
            ]
            composed_version = "|".join([version, *prior_composer_versions])
            updates.append(
                (
                    tier,
                    json.dumps(combined_reasons),
                    composed_version,
                    json.dumps(field_trust, sort_keys=True),
                    json.dumps(sources),
                    r["penetrance_id"],
                )
            )

        conn.executemany(
            "UPDATE penetrance_data SET trust_tier=?, trust_reasons=?, "
            "trust_rule_version=?, field_trust=?, trust_sources=? "
            "WHERE penetrance_id=?",
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
