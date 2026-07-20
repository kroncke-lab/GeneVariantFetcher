"""Risk-triage for the expensive per-paper final check (Step 3.8).

The source-grounded final check (``azure_ai/gpt-5.6-sol`` at ``xhigh``) is the
dominant run cost and today runs on EVERY paper. This module decides, per paper
and **without any LLM call**, which of three lanes a paper takes:

* ``full``  — the full ``@xhigh`` source-grounded audit (current behavior).
* ``cheap`` — the SAME source-grounded prompt shape on a cheaper model / lower
  reasoning effort, escalated to ``full`` on any flag / miss / low confidence.
* ``skip``  — no LLM call (nothing to check: no counts AND no usable source).

Design rationale (three-way review — codex/grok/agy, 2026-07-20; see
``docs/PROTOCOL_COST_EVAL.md``):

* The final check does TWO jobs: a **count-trust audit** and a **completeness**
  check (which reported carrier groups the pipeline missed). The free
  deterministic trust gate (:mod:`pipeline.trust_gate`) already covers most of
  the count-trust audit, but **completeness has no free substitute and does NOT
  correlate with extraction difficulty** — a clean, high-confidence extraction
  can still miss an entire supplemental table. Therefore we NEVER map
  "clean" → ``skip``; a clean paper with usable source demotes to a
  completeness-preserving ``cheap`` pass.
* ``skip`` is reserved for the existing ``empty_no_source`` case (no counts and
  no usable source), plus deterministically-clean DB-only papers where no
  completeness claim is possible anyway.
* Zero counts WITH source is the HIGHEST completeness value → ``full``.

Signals are limited to what is available in the migrated DB (trust-gate reasons,
``extraction_confidence``, fact/quote/role counts). The paper census risk flags
are NOT yet persisted to the DB (they live only in the extraction JSON); wiring
them in is a follow-up that would sharpen the completeness lane. Until then a
missing/unknown signal is treated as *risky*, never clean.

Pure and unit-tested; mirrors :mod:`pipeline.trust_gate`. The thin SQLite/LLM
orchestration lives in :func:`pipeline.paper_final_check.apply_paper_final_check`.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

# ---- Thresholds (v1; calibrate on the 101-paper shadow set before enforcing) ----
# A count surface this large is worth the full audit regardless of trust tier.
MAX_SIMPLE_COUNT_FACTS = 24
# This many distinct count roles on one paper signals a mixed/complex table.
MAX_SIMPLE_COUNT_ROLES = 3
# Below this fraction of count facts carrying a source quote, the trust gate
# cannot have source-grounded the residual count-trust errors → full.
MIN_QUOTE_COVERAGE = 0.25
# extraction_confidence values that are NOT a clean signal. "medium" is the
# hardcoded default on several extraction paths, so it is noise, not signal;
# only an explicit "low" (or a missing value) counts as risk.
RISKY_CONFIDENCE = frozenset({"low"})
# Cheap-tier self-reported confidence at/above which we do NOT escalate on
# confidence alone (escalation still fires on flags / misses / coverage gaps).
CHEAP_MIN_CONFIDENCE = 0.90
# Stable random audit rate: this fraction of otherwise-cheap papers are forced
# to ``full`` so silent cheap-tier false negatives stay measurable.
AUDIT_SAMPLE_RATE = 0.10

TRIAGE_POLICY_VERSION = "tri1"


def policy_hash() -> str:
    """Short content hash of the triage policy + thresholds, so a stored decision
    records exactly which policy generation produced it (folded into the
    final-check version string when triage is enforced)."""
    payload = json.dumps(
        {
            "policy": TRIAGE_POLICY_VERSION,
            "max_simple_count_facts": MAX_SIMPLE_COUNT_FACTS,
            "max_simple_count_roles": MAX_SIMPLE_COUNT_ROLES,
            "min_quote_coverage": MIN_QUOTE_COVERAGE,
            "risky_confidence": sorted(RISKY_CONFIDENCE),
            "cheap_min_confidence": CHEAP_MIN_CONFIDENCE,
            "audit_sample_rate": AUDIT_SAMPLE_RATE,
        },
        sort_keys=True,
    )
    return (
        f"{TRIAGE_POLICY_VERSION}-" + hashlib.sha256(payload.encode()).hexdigest()[:8]
    )


@dataclass(frozen=True)
class PaperRiskView:
    """Everything the triage predicate needs about one paper — pure inputs."""

    pmid: str
    n_count_facts: int  # extracted rows carrying at least one non-null count
    any_quarantine: bool  # trust gate quarantined any field on this paper
    any_trust_reason: bool  # trust gate raised any reason on this paper
    frac_facts_with_quote: float  # fraction of count facts backed by a quote (0..1)
    n_distinct_count_roles: int  # distinct carriers_count_type values seen
    extraction_confidence: Optional[str]  # 'high'|'medium'|'low'|None
    has_usable_source: bool  # scouted source on disk usable for completeness
    empty_no_source: bool  # no count facts AND no usable source (existing skip)


@dataclass(frozen=True)
class TriageDecision:
    tier: str  # 'full' | 'cheap' | 'skip'
    reasons: tuple[str, ...] = field(default_factory=tuple)


def _stable_audit_pick(pmid: str, rate: float) -> bool:
    """Deterministic per-PMID audit sampling (no RNG, resume-safe): keeps a fixed
    ~``rate`` fraction of PMIDs so a stable random slice of cheap papers is still
    audited at full effort."""
    if rate <= 0:
        return False
    if rate >= 1:
        return True
    bucket = int(hashlib.sha256(pmid.encode()).hexdigest()[:8], 16) % 10_000
    return bucket < int(rate * 10_000)


def decide_final_check_tier(
    view: PaperRiskView,
    *,
    audit_sample_rate: float = AUDIT_SAMPLE_RATE,
    zero_count_source_tier: str = "full",
) -> TriageDecision:
    """Route one paper to ``full`` / ``cheap`` / ``skip`` from DB-only signals.

    Never returns ``skip`` for a paper with usable source — completeness is not
    predictable from clean count signals, so clean source-grounded papers demote
    to ``cheap`` (a completeness-preserving pass), not ``skip``.

    ``zero_count_source_tier`` controls the biggest cost lever: papers that
    extracted NO counts but have source. The conservative default ``"full"``
    (the three-way review's recommendation — "highest completeness value") sends
    them to the full ``@xhigh`` audit; on a real cardiac DB these are ~55% of
    papers, so triage then only saves ~30%. Setting it to ``"cheap"`` routes them
    through the cheap completeness lane instead (still completeness-preserving:
    the cheap pass enumerates carrier groups and escalates on any found miss, and
    the audit sample catches systematic false negatives) for a much larger saving.
    Which is right is a calibration question — measure cheap-vs-full completeness
    agreement on the zero-count bucket before enforcing ``"cheap"``.
    """
    reasons: list[str] = []

    # 1) Nothing to check: no counts and no usable source (existing behavior).
    if view.empty_no_source:
        return TriageDecision("skip", ("empty_no_source",))

    # 2) Zero counts WITH source: completeness-critical. Routed per the knob.
    if view.n_count_facts == 0 and view.has_usable_source:
        if zero_count_source_tier == "cheap" and not _stable_audit_pick(
            view.pmid, audit_sample_rate
        ):
            return TriageDecision("cheap", ("zero_counts_with_source",))
        return TriageDecision("full", ("zero_counts_with_source",))

    # 3) Hard risk → full.
    if view.any_quarantine or view.any_trust_reason:
        reasons.append("trust_gate_hit")
    conf = (view.extraction_confidence or "").strip().lower()
    if conf in RISKY_CONFIDENCE or view.extraction_confidence is None:
        reasons.append("low_or_unknown_confidence")
    if view.n_count_facts >= 1 and view.frac_facts_with_quote < MIN_QUOTE_COVERAGE:
        reasons.append("thin_provenance")
    if view.n_count_facts >= MAX_SIMPLE_COUNT_FACTS:
        reasons.append("complex_count_surface")
    if view.n_distinct_count_roles >= MAX_SIMPLE_COUNT_ROLES:
        reasons.append("mixed_count_roles")
    if reasons:
        return TriageDecision("full", tuple(reasons))

    # 4) Stable random audit slice of otherwise-cheap papers → full (measures
    #    silent cheap-tier false negatives).
    if _stable_audit_pick(view.pmid, audit_sample_rate):
        return TriageDecision("full", ("audit_sample",))

    # 5) Clean. With usable source → cheap completeness-preserving pass. DB-only
    #    clean (no source) has no completeness lane and its counts are already
    #    deterministically vetted by the trust gate → skip.
    if view.has_usable_source:
        return TriageDecision("cheap", ("clean_with_source",))
    return TriageDecision("skip", ("deterministic_clean_db_only",))


def _table_cols(conn: sqlite3.Connection, table: str) -> set[str]:
    try:
        return {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
    except sqlite3.Error:
        return set()


def collect_risk_views(
    db_path: str | Path,
    *,
    has_source: Callable[[str], bool],
) -> dict[str, PaperRiskView]:
    """Build a :class:`PaperRiskView` per paper from a migrated DB (no LLM).

    ``has_source(pmid) -> bool`` supplies whether usable source text is on disk
    (the real adapter uses the final-check source resolver; the offline report
    approximates it from the corpus). Papers with zero count rows are still
    included — completeness is most valuable there.
    """
    conn = sqlite3.connect(str(db_path))
    try:
        conn.row_factory = sqlite3.Row
        pd_cols = _table_cols(conn, "penetrance_data")
        has_trust = "trust_tier" in pd_cols
        # Seed every paper (even zero-count), like gather_paper_payloads.
        views: dict[str, dict[str, Any]] = {}
        for r in conn.execute("SELECT pmid FROM papers"):
            views[str(r["pmid"])] = {
                "n_count_facts": 0,
                "any_quarantine": False,
                "any_trust_reason": False,
                "quoted": 0,
                "roles": set(),
            }
        # count_provenance / key_quotes per (variant_id, pmid) for quotes + roles.
        prov: dict[tuple[int, str], dict[str, Any]] = {}
        if "count_provenance" in _table_cols(conn, "variant_papers"):
            for r in conn.execute(
                "SELECT variant_id, pmid, count_provenance, key_quotes FROM variant_papers"
            ):
                role = None
                raw = r["count_provenance"]
                if raw:
                    try:
                        role = (json.loads(raw) or {}).get("carriers_count_type")
                    except (TypeError, ValueError, json.JSONDecodeError):
                        role = None
                quoted = bool(
                    r["key_quotes"] and str(r["key_quotes"]).strip() not in ("", "[]")
                )
                prov[(int(r["variant_id"]), str(r["pmid"]))] = {
                    "role": role,
                    "quoted": quoted,
                }

        trust_sel = ", trust_tier, trust_reasons" if has_trust else ""
        for r in conn.execute(
            f"SELECT pmid, variant_id, total_carriers_observed, affected_count, "
            f"unaffected_count, uncertain_count{trust_sel} FROM penetrance_data"
        ):
            pmid = str(r["pmid"])
            v = views.setdefault(
                pmid,
                {
                    "n_count_facts": 0,
                    "any_quarantine": False,
                    "any_trust_reason": False,
                    "quoted": 0,
                    "roles": set(),
                },
            )
            has_count = any(
                r[c] is not None
                for c in (
                    "total_carriers_observed",
                    "affected_count",
                    "unaffected_count",
                    "uncertain_count",
                )
            )
            if not has_count:
                continue
            v["n_count_facts"] += 1
            p = prov.get((int(r["variant_id"]), pmid), {})
            if p.get("quoted"):
                v["quoted"] += 1
            if p.get("role"):
                v["roles"].add(p["role"])
            if has_trust:
                if str(r["trust_tier"] or "") == "quarantine":
                    v["any_quarantine"] = True
                tr = r["trust_reasons"]
                if tr and str(tr).strip() not in ("", "[]"):
                    v["any_trust_reason"] = True

        # extraction_confidence: latest row per pmid.
        conf: dict[str, Optional[str]] = {}
        if "extraction_confidence" in _table_cols(conn, "extraction_metadata"):
            for r in conn.execute(
                "SELECT pmid, extraction_confidence FROM extraction_metadata "
                "ORDER BY extraction_id"
            ):
                conf[str(r["pmid"])] = r["extraction_confidence"]

        out: dict[str, PaperRiskView] = {}
        for pmid, v in views.items():
            n = v["n_count_facts"]
            src = bool(has_source(pmid))
            out[pmid] = PaperRiskView(
                pmid=pmid,
                n_count_facts=n,
                any_quarantine=v["any_quarantine"],
                any_trust_reason=v["any_trust_reason"],
                frac_facts_with_quote=(v["quoted"] / n) if n else 0.0,
                n_distinct_count_roles=len(v["roles"]),
                extraction_confidence=conf.get(pmid),
                has_usable_source=src,
                empty_no_source=(n == 0 and not src),
            )
        return out
    finally:
        conn.close()


def should_escalate_cheap(
    result: dict[str, Any], *, min_confidence: float = CHEAP_MIN_CONFIDENCE
) -> tuple[bool, tuple[str, ...]]:
    """Whether a cheap-tier result must be re-run at full ``@xhigh``.

    Never trusts the cheap model's self-reported confidence alone: escalate on
    any flag, any derived miss, a non-complete completeness status, a truncated
    view, an invalid response, or low confidence.
    """
    reasons: list[str] = []
    if not isinstance(result, dict):
        return True, ("invalid_response",)
    if (
        str(result.get("verdict") or result.get("paper_verdict") or "").lower()
        == "flag"
    ):
        reasons.append("cheap_flagged")
    if result.get("flags"):
        reasons.append("cheap_flags")
    comp = result.get("completeness") or {}
    if isinstance(comp, dict):
        if int(comp.get("n_missing") or 0) > 0:
            reasons.append("cheap_missing_carriers")
        status = str(comp.get("status") or "").lower()
        if status and status not in {"complete", "ok"}:
            reasons.append("completeness_not_complete")
        if comp.get("truncated") or result.get("source_truncated"):
            reasons.append("source_truncated")
    conf = result.get("confidence")
    try:
        if conf is not None and float(conf) < min_confidence:
            reasons.append("low_cheap_confidence")
    except (TypeError, ValueError):
        reasons.append("unparseable_confidence")
    return (bool(reasons), tuple(reasons))
