"""Compose source-grounded paper-final-check findings into count trust.

The final checker is intentionally an auditor: it preserves raw extracted
counts and records structured findings.  This module is the enforcement
boundary.  It binds only explicit ``fact_id`` + field findings to
``penetrance_data``, composes them with the structural trust reasons already on
the row, and writes a field-level trusted projection.  It never NULLs or
deletes raw evidence and never applies a paper-level verdict to every fact.

Completeness gaps remain a re-extraction signal.  They are reported in the
returned statistics but do not make otherwise-supported captured facts
untrusted.
"""

from __future__ import annotations

import hashlib
import json
import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

COUNT_FIELDS = ("total_carriers", "affected", "unaffected", "uncertain")
FIELD_ALIASES = {
    "carriers": "total_carriers",
    "carrier": "total_carriers",
    "patient_count": "total_carriers",
    "total_carriers_observed": "total_carriers",
    "affected_count": "affected",
    "unaffected_count": "unaffected",
    "uncertain_count": "uncertain",
}
PFC_REASON_PREFIX = "paper_final_check:"
COMPOSE_POLICY_VERSION = "pfcg2"
SEVERITY_RANK = {"low": 0, "medium": 1, "high": 2}
ENFORCEABLE_REASON_CODES = frozenset(
    {
        "count_is_total",
        "population_count",
        "wrong_column",
        "arith_inconsistent",
        "phenotype_misclassified",
        "wrong_gene",
    }
)


def compose_version(*, min_severity: str, require_source_grounded: bool) -> str:
    """Pin the exact enforcement policy written into trust metadata."""
    payload = json.dumps(
        {
            "policy": COMPOSE_POLICY_VERSION,
            "fields": COUNT_FIELDS,
            "min_severity": min_severity,
            "require_source_grounded": require_source_grounded,
            "enforceable_reason_codes": sorted(ENFORCEABLE_REASON_CODES),
        },
        sort_keys=True,
    )
    return (
        f"{COMPOSE_POLICY_VERSION}-"
        + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
    )


def _load_json(raw: Any, default: Any) -> Any:
    if raw in (None, ""):
        return default
    if isinstance(raw, (dict, list)):
        return raw
    try:
        return json.loads(str(raw))
    except (TypeError, ValueError, json.JSONDecodeError):
        return default


def _reason_list(raw: Any) -> list[str]:
    values = _load_json(raw, [])
    if not isinstance(values, list):
        return []
    return sorted({str(value) for value in values if str(value).strip()})


def _strip_compose_version(raw: Any) -> str:
    return "|".join(
        part
        for part in str(raw or "").split("|")
        if part and not part.startswith("pfcg")
    )


def _structural_payload_signature(payload: dict[str, Any]) -> str:
    """Hash exactly the structural payload shape used by the reviewer.

    ``gather_paper_payloads`` already removes prior paper-final-check reasons
    while retaining a reasonless legacy quarantine. Re-deriving the tier here
    would make that valid reviewer payload look stale.
    """
    from pipeline.paper_final_check import payload_content_hash

    return payload_content_hash(payload)[:12]


def _ensure_column(
    conn: sqlite3.Connection, table: str, name: str, declaration: str
) -> None:
    columns = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
    if name not in columns:
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {name} {declaration}")


def ensure_gate_schema(conn: sqlite3.Connection) -> None:
    """Add the backward-compatible trust projection and action audit schema."""
    _ensure_column(conn, "penetrance_data", "trust_tier", "TEXT DEFAULT 'trusted'")
    _ensure_column(conn, "penetrance_data", "trust_reasons", "TEXT")
    _ensure_column(conn, "penetrance_data", "trust_rule_version", "TEXT")
    _ensure_column(conn, "penetrance_data", "field_trust", "TEXT")
    _ensure_column(conn, "penetrance_data", "trust_sources", "TEXT")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS paper_final_check_actions (
            action_id       INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid            TEXT,
            penetrance_id   INTEGER,
            variant         TEXT,
            severity        TEXT,
            fields_json     TEXT,
            reason_code     TEXT,
            issue           TEXT,
            match_status    TEXT,
            action          TEXT,
            check_version   TEXT,
            applied_at      TEXT
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pfca_pmid ON paper_final_check_actions(pmid)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pfca_penetrance "
        "ON paper_final_check_actions(penetrance_id)"
    )


def _normalize_flag(raw: Any) -> dict[str, Any]:
    item = raw if isinstance(raw, dict) else {}
    severity = str(item.get("severity") or "medium").strip().lower()
    if severity not in SEVERITY_RANK:
        severity = "medium"

    raw_fact_ids = item.get("fact_ids")
    if raw_fact_ids is None and item.get("fact_id") is not None:
        raw_fact_ids = [item.get("fact_id")]
    if not isinstance(raw_fact_ids, list):
        raw_fact_ids = []
    fact_ids: list[int] = []
    for value in raw_fact_ids:
        try:
            fact_id = int(value)
        except (TypeError, ValueError):
            continue
        if fact_id > 0 and fact_id not in fact_ids:
            fact_ids.append(fact_id)

    raw_variant_ids = item.get("variant_ids")
    if not isinstance(raw_variant_ids, list):
        raw_variant_ids = []
    variant_ids: list[int] = []
    for value in raw_variant_ids:
        try:
            variant_id = int(value)
        except (TypeError, ValueError):
            continue
        if variant_id > 0:
            variant_ids.append(variant_id)

    raw_fields = item.get("fields")
    if isinstance(raw_fields, str):
        raw_fields = [raw_fields]
    fields = []
    for value in raw_fields or []:
        field = str(value or "").strip().lower()
        field = FIELD_ALIASES.get(field, field)
        if field in COUNT_FIELDS and field not in fields:
            fields.append(field)

    reason_code = str(item.get("reason_code") or "other").strip().lower()
    reason_code = (
        "".join(
            char if char.isalnum() or char == "_" else "_" for char in reason_code
        ).strip("_")
        or "other"
    )
    return {
        "fact_ids": fact_ids,
        "variant_ids": variant_ids,
        "variant": str(item.get("variant") or "").strip() or "(unspecified)",
        "fields": fields,
        "reason_code": reason_code,
        "evidence_quote": str(item.get("evidence_quote") or "").strip(),
        "evidence_quote_verified": bool(item.get("evidence_quote_verified")),
        "issue": str(item.get("issue") or "").strip(),
        "severity": severity,
    }


def apply_paper_final_check_gate(
    db_path: str | Path,
    *,
    min_severity: str = "high",
    require_source_grounded: bool = True,
) -> dict[str, Any]:
    """Apply exact fact/field final-check findings to the trusted projection.

    The operation is a deterministic full recomposition.  Existing final-check
    reasons are removed first, structural/legacy reasons are preserved, and the
    current ``paper_final_check`` rows are then applied.  Low-severity,
    ungrounded, unbound, or ambiguous findings remain auditable but cannot
    quarantine a fact.
    """
    min_severity = str(min_severity or "high").strip().lower()
    if min_severity not in SEVERITY_RANK:
        raise ValueError("min_severity must be low, medium, or high")

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    try:
        tables = {
            row[0]
            for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }
        if "penetrance_data" not in tables:
            return {"skipped": "penetrance_data table absent"}

        ensure_gate_schema(conn)
        policy_version = compose_version(
            min_severity=min_severity,
            require_source_grounded=require_source_grounded,
        )
        applied_at = datetime.now(timezone.utc).isoformat()
        stats: dict[str, Any] = {
            "policy_version": policy_version,
            "min_severity": min_severity,
            "require_source_grounded": require_source_grounded,
            "papers": 0,
            "source_grounded_papers": 0,
            "missing_carriers": 0,
            "flags_total": 0,
            "actionable_flags": 0,
            "advisory_flags": 0,
            "unresolved_actionable": 0,
            "ungrounded_actionable": 0,
            "unverified_actionable": 0,
            "advisory_reason_flags": 0,
            "stale_actionable": 0,
            "stale_missing_carriers": 0,
            "stale_papers": 0,
            "stale_reviewer_version": 0,
            "applied_bindings": 0,
            "applied_facts": 0,
            "fields_quarantined": {field: 0 for field in COUNT_FIELDS},
            "trusted": 0,
            "quarantine": 0,
        }

        rows = [
            dict(row)
            for row in conn.execute(
                """
                SELECT penetrance_id, variant_id, pmid, trust_tier, trust_reasons,
                       trust_rule_version
                FROM penetrance_data
                ORDER BY penetrance_id
                """
            )
        ]
        states: dict[int, dict[str, Any]] = {}
        for row in rows:
            all_reasons = _reason_list(row.get("trust_reasons"))
            prior_pfc = [
                reason for reason in all_reasons if reason.startswith(PFC_REASON_PREFIX)
            ]
            base_reasons = [
                reason
                for reason in all_reasons
                if not reason.startswith(PFC_REASON_PREFIX)
            ]
            legacy_quarantine = (
                str(row.get("trust_tier") or "").lower() == "quarantine"
                and not prior_pfc
                and not base_reasons
            )
            base_quarantine = bool(base_reasons or legacy_quarantine)
            states[int(row["penetrance_id"])] = {
                "variant_id": int(row["variant_id"]),
                "pmid": str(row.get("pmid")),
                "reasons": set(base_reasons),
                "field_trust": {
                    field: "quarantine" if base_quarantine else "trusted"
                    for field in COUNT_FIELDS
                },
                "sources": {"legacy" if legacy_quarantine else "structural"}
                if base_quarantine
                else set(),
                "base_version": _strip_compose_version(row.get("trust_rule_version")),
            }

        conn.execute("DELETE FROM paper_final_check_actions")

        pfc_columns = (
            {row[1] for row in conn.execute("PRAGMA table_info(paper_final_check)")}
            if "paper_final_check" in tables
            else set()
        )
        action_rows: list[tuple[Any, ...]] = []
        applied_fact_ids: set[int] = set()
        quarantined_fields: set[tuple[int, str]] = set()

        if pfc_columns:
            # The stored check_version ends in the semantic payload hash that
            # produced its fact IDs. Durable findings may be reused during a
            # reviewer outage only while that hash still matches the current DB.
            from pipeline.paper_final_check import (
                PROMPT_VERSION,
                SUMMARY_PROMPT_VERSION,
                gather_paper_payloads,
            )

            current_payload_sigs = {
                str(payload.get("pmid")): _structural_payload_signature(payload)
                for payload in gather_paper_payloads(conn)
            }
            source_expr = "source_grounded" if "source_grounded" in pfc_columns else "0"
            missing_expr = "n_missing" if "n_missing" in pfc_columns else "0"
            version_expr = "check_version" if "check_version" in pfc_columns else "NULL"
            prompt_expr = (
                "prompt_version" if "prompt_version" in pfc_columns else "NULL"
            )
            for paper in conn.execute(
                f"""
                SELECT pmid, verdict, flags_json,
                       {source_expr} AS source_grounded,
                       {missing_expr} AS n_missing,
                       {version_expr} AS check_version,
                       {prompt_expr} AS prompt_version
                FROM paper_final_check
                ORDER BY pmid
                """
            ):
                stats["papers"] += 1
                grounded = bool(paper["source_grounded"])
                if grounded:
                    stats["source_grounded_papers"] += 1
                flags = _load_json(paper["flags_json"], [])
                if not isinstance(flags, list):
                    flags = []
                check_version = str(paper["check_version"] or "")
                stored_payload_sig = check_version.rsplit("-", 1)[-1]
                payload_current = len(
                    stored_payload_sig
                ) == 12 and stored_payload_sig == current_payload_sigs.get(
                    str(paper["pmid"]), ""
                )
                prompt_version = str(paper["prompt_version"] or "")
                reviewer_current = prompt_version in {
                    PROMPT_VERSION,
                    SUMMARY_PROMPT_VERSION,
                }
                if not reviewer_current:
                    stats["stale_reviewer_version"] += 1
                if not payload_current or not reviewer_current:
                    stats["stale_papers"] += 1
                    if grounded:
                        stats["stale_missing_carriers"] += int(paper["n_missing"] or 0)
                    for raw_flag in flags:
                        flag = _normalize_flag(raw_flag)
                        stats["flags_total"] += 1
                        if (
                            SEVERITY_RANK[flag["severity"]]
                            >= SEVERITY_RANK[min_severity]
                            and flag["reason_code"] in ENFORCEABLE_REASON_CODES
                        ):
                            stats["actionable_flags"] += 1
                            stats["stale_actionable"] += 1
                        else:
                            stats["advisory_flags"] += 1
                            if flag["reason_code"] not in ENFORCEABLE_REASON_CODES:
                                stats["advisory_reason_flags"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                None,
                                flag["variant"],
                                flag["severity"],
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "stale_payload_hash"
                                if not payload_current
                                else "stale_reviewer_version",
                                "not_applied",
                                paper["check_version"],
                                applied_at,
                            )
                        )
                    continue

                if grounded:
                    stats["missing_carriers"] += int(paper["n_missing"] or 0)
                for raw_flag in flags:
                    flag = _normalize_flag(raw_flag)
                    stats["flags_total"] += 1
                    severity = flag["severity"]
                    actionable = SEVERITY_RANK[severity] >= SEVERITY_RANK[min_severity]
                    if not actionable:
                        stats["advisory_flags"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                None,
                                flag["variant"],
                                severity,
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "not_evaluated",
                                "advisory_below_threshold",
                                paper["check_version"],
                                applied_at,
                            )
                        )
                        continue

                    if flag["reason_code"] not in ENFORCEABLE_REASON_CODES:
                        stats["advisory_reason_flags"] += 1
                        stats["advisory_flags"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                None,
                                flag["variant"],
                                severity,
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "reason_not_enforceable",
                                "advisory_only",
                                paper["check_version"],
                                applied_at,
                            )
                        )
                        continue

                    stats["actionable_flags"] += 1
                    if require_source_grounded and not grounded:
                        stats["ungrounded_actionable"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                None,
                                flag["variant"],
                                severity,
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "ungrounded",
                                "not_applied",
                                paper["check_version"],
                                applied_at,
                            )
                        )
                        continue

                    if grounded and not flag["evidence_quote_verified"]:
                        stats["unverified_actionable"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                None,
                                flag["variant"],
                                severity,
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "source_quote_unverified",
                                "not_applied",
                                paper["check_version"],
                                applied_at,
                            )
                        )
                        continue

                    if (
                        not flag["fact_ids"]
                        or not flag["fields"]
                        or len(flag["variant_ids"]) != len(flag["fact_ids"])
                    ):
                        stats["unresolved_actionable"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                None,
                                flag["variant"],
                                severity,
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "missing_fact_variant_or_field_binding",
                                "not_applied",
                                paper["check_version"],
                                applied_at,
                            )
                        )
                        continue

                    for binding_index, fact_id in enumerate(flag["fact_ids"]):
                        state = states.get(fact_id)
                        expected_variant_id = flag["variant_ids"][binding_index]
                        if (
                            state is None
                            or state["pmid"] != str(paper["pmid"])
                            or state["variant_id"] != expected_variant_id
                        ):
                            stats["unresolved_actionable"] += 1
                            action_rows.append(
                                (
                                    str(paper["pmid"]),
                                    fact_id,
                                    flag["variant"],
                                    severity,
                                    json.dumps(flag["fields"]),
                                    flag["reason_code"],
                                    flag["issue"],
                                    "fact_variant_binding_mismatch",
                                    "not_applied",
                                    paper["check_version"],
                                    applied_at,
                                )
                            )
                            continue

                        for field in flag["fields"]:
                            reason = (
                                f"{PFC_REASON_PREFIX}{flag['reason_code']}:"
                                f"{severity}:{field}"
                            )
                            state["reasons"].add(reason)
                            state["field_trust"][field] = "quarantine"
                            quarantined_fields.add((fact_id, field))
                        state["sources"].add("paper_final_check")
                        applied_fact_ids.add(fact_id)
                        stats["applied_bindings"] += 1
                        action_rows.append(
                            (
                                str(paper["pmid"]),
                                fact_id,
                                flag["variant"],
                                severity,
                                json.dumps(flag["fields"]),
                                flag["reason_code"],
                                flag["issue"],
                                "exact_fact_id",
                                "quarantine_fields",
                                paper["check_version"],
                                applied_at,
                            )
                        )
        stats["applied_facts"] = len(applied_fact_ids)
        for _, field in quarantined_fields:
            stats["fields_quarantined"][field] += 1

        updates: list[tuple[str, str, str, str, str, int]] = []
        for fact_id, state in states.items():
            tier = (
                "quarantine"
                if any(value == "quarantine" for value in state["field_trust"].values())
                else "trusted"
            )
            stats[tier] += 1
            base_version = state["base_version"]
            version = (
                f"{base_version}|{policy_version}" if base_version else policy_version
            )
            updates.append(
                (
                    tier,
                    json.dumps(sorted(state["reasons"])),
                    version,
                    json.dumps(state["field_trust"], sort_keys=True),
                    json.dumps(sorted(state["sources"])),
                    fact_id,
                )
            )
        conn.executemany(
            """
            UPDATE penetrance_data
            SET trust_tier=?, trust_reasons=?, trust_rule_version=?,
                field_trust=?, trust_sources=?
            WHERE penetrance_id=?
            """,
            updates,
        )
        if action_rows:
            conn.executemany(
                """
                INSERT INTO paper_final_check_actions (
                    pmid, penetrance_id, variant, severity, fields_json,
                    reason_code, issue, match_status, action, check_version,
                    applied_at
                ) VALUES (?,?,?,?,?,?,?,?,?,?,?)
                """,
                action_rows,
            )
        conn.commit()
        logger.info(
            "paper final check gate: %d fact(s), %d binding(s), "
            "%d unresolved, %d grounded missing carrier group(s) [%s]",
            stats["applied_facts"],
            stats["applied_bindings"],
            stats["unresolved_actionable"],
            stats["missing_carriers"],
            policy_version,
        )
        return stats
    finally:
        conn.close()
