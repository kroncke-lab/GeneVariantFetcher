"""Per-paper "final check" — a gpt-5.6-sol sniff test over each paper's
extracted variant counts and their provenance.

Autonomy QC. After extraction, migration, and the gold-free structural trust
gate (see :mod:`pipeline.trust_gate`), this step asks a strong reasoning model
(default ``azure_ai/gpt-5.6-sol`` at ``xhigh`` effort) to eyeball **each paper's**
extracted variant/count rows against whatever provenance was captured (source
location, verbatim quote, and the count's role/label) and flag anything that
fails the smell test — count-role confusion, a count the quote does not support,
a population/allele figure masquerading as a carrier count, an implausible
magnitude.

It is a **judgment layer, not a mutation**: it RECORDS a soft per-paper verdict
in ``paper_final_check`` and never edits or deletes a raw count. The pure helpers
(:func:`gather_paper_payloads`, :func:`build_paper_check_prompt`,
:func:`normalize_result`) have no LLM/DB coupling and are unit-tested directly;
:func:`apply_paper_final_check` is the thin SQLite + LLM adapter and accepts an
injected ``checker`` so tests never hit the network.
"""

from __future__ import annotations

import hashlib
import json
import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Protocol

logger = logging.getLogger(__name__)

# Defaults. The mainline (gvf-run) resolves the real model/effort from
# ``config.settings`` (PAPER_FINAL_CHECK_MODEL / _REASONING_EFFORT); these are
# the standalone fallbacks so the module is usable on its own.
DEFAULT_MODEL = "azure_ai/gpt-5.6-sol"
DEFAULT_REASONING_EFFORT = "xhigh"
# xhigh spends hidden reasoning tokens before the visible JSON, so give the
# completion generous headroom (clamp_max_tokens caps it per model).
DEFAULT_MAX_TOKENS = 16000
# Cap facts per paper so a 70-row paper stays a "quick" sniff test.
MAX_FACTS_PER_PAPER = 60
# Bound each provenance field so the prompt stays compact.
MAX_QUOTE_CHARS = 800

PROMPT_VERSION = "pfc1"
VALID_VERDICTS = ("ok", "flag")
VALID_SEVERITIES = ("low", "medium", "high")

_PROMPT_TEMPLATE = """You are the final quality-control reviewer (a "sniff test") for an automated
pipeline that extracts variant carrier counts from biomedical papers. You are
NOT re-extracting the paper. For ONE paper you are given the variant/count rows
the pipeline extracted and whatever provenance it captured (a source location,
verbatim quote(s), and the count's declared role/label). Judge quickly but
carefully whether each row looks trustworthy, and flag the ones that do not.

A row is SUSPECT if any of these apply:
- The count looks like a cohort/study denominator (everyone enrolled, screened,
  sampled, or sequenced), not the number of people who carry THIS variant.
- The count is a study-wide, family-set, or mutation-class total copied onto a
  single variant.
- The number is a population/allele figure (gnomAD, ExAC, TOPMed, 1000 Genomes,
  MAF, AC/AN, allele frequency), not a clinical carrier count.
- affected + unaffected exceeds total_carriers.
- The number is an assay replicate/cell count, or a table value that is really a
  prediction score, allele frequency, exon/domain number, or genomic position.
- The count has no supporting quote or source location (treat as low confidence).
- The magnitude is implausible for the carriers of a single variant.

Paper (JSON):
{payload}

Return STRICT JSON only, no prose outside the object:
{{
  "paper_verdict": "ok" | "flag",
  "confidence": <number between 0 and 1: your confidence the extraction for this
                 paper is correct overall>,
  "flags": [
    {{"variant": "<the variant field from a suspect row>",
      "issue": "<short reason it is suspect>",
      "severity": "low" | "medium" | "high"}}
  ],
  "summary": "<one or two sentences: what this study appears to be, and whether
              its extracted counts pass the smell test>"
}}
Only include a flag for rows that are actually suspect. If every row looks fine,
return "ok", an empty flags list, and a brief summary."""


class _Checker(Protocol):
    def check(self, payload: dict[str, Any]) -> dict[str, Any]:  # pragma: no cover
        ...


def check_version(model: str, reasoning_effort: Optional[str]) -> str:
    """Content hash of the prompt generation + model + effort, so a stored
    ``paper_final_check`` row pins exactly which reviewer produced it."""
    payload = json.dumps(
        {
            "prompt_version": PROMPT_VERSION,
            "template": _PROMPT_TEMPLATE,
            "model": model,
            "reasoning_effort": reasoning_effort,
        },
        sort_keys=True,
    )
    return (
        f"{PROMPT_VERSION}-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
    )


def _as_int(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _clip(value: Any, limit: int = MAX_QUOTE_CHARS) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return text[:limit]


def _load_json(raw: Any) -> Any:
    if not raw:
        return None
    if isinstance(raw, (dict, list)):
        return raw
    try:
        return json.loads(raw)
    except (TypeError, ValueError, json.JSONDecodeError):
        return None


def _penetrance_has_trust(conn: sqlite3.Connection) -> bool:
    cols = {row[1] for row in conn.execute("PRAGMA table_info(penetrance_data)")}
    return "trust_tier" in cols


def gather_paper_payloads(
    conn: sqlite3.Connection,
    pmids: Optional[list[str]] = None,
    max_facts: int = MAX_FACTS_PER_PAPER,
) -> list[dict[str, Any]]:
    """Assemble one sniff-test payload per paper that has extracted counts.

    Provenance is read from ``variant_papers`` (``key_quotes`` / ``source_location``
    / ``count_provenance``), which is present across old and new DBs; the trust
    tier is included only when the trust-gate columns exist. Papers with no
    ``penetrance_data`` rows are skipped — there is nothing to sniff-test.
    """
    conn.row_factory = sqlite3.Row
    sel_trust = (
        ", pd.trust_tier AS trust_tier, pd.trust_reasons AS trust_reasons"
        if _penetrance_has_trust(conn)
        else ""
    )
    where = ""
    params: list[Any] = []
    if pmids:
        placeholders = ",".join("?" for _ in pmids)
        where = f" WHERE pd.pmid IN ({placeholders})"
        params = list(pmids)
    sql = f"""
        SELECT pd.pmid AS pmid,
               COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position)
                   AS variant,
               pd.total_carriers_observed AS total_carriers,
               pd.affected_count          AS affected,
               pd.unaffected_count        AS unaffected,
               pd.uncertain_count         AS uncertain,
               vp.source_location         AS source_location,
               vp.key_quotes              AS key_quotes,
               vp.count_provenance        AS count_provenance,
               p.title                    AS title,
               COALESCE(p.gene_symbol, v.gene_symbol) AS gene_symbol
               {sel_trust}
        FROM penetrance_data pd
        JOIN variants v ON v.variant_id = pd.variant_id
        LEFT JOIN variant_papers vp
               ON vp.variant_id = pd.variant_id AND vp.pmid = pd.pmid
        LEFT JOIN papers p ON p.pmid = pd.pmid
        {where}
        ORDER BY pd.pmid
    """
    grouped: dict[str, dict[str, Any]] = {}
    for row in conn.execute(sql, params):
        pmid = str(row["pmid"])
        bucket = grouped.setdefault(
            pmid,
            {
                "pmid": pmid,
                "gene": row["gene_symbol"],
                "title": row["title"],
                "facts": [],
                "n_facts_total": 0,
            },
        )
        bucket["n_facts_total"] += 1
        if len(bucket["facts"]) >= max_facts:
            continue

        counts = {
            "total_carriers": _as_int(row["total_carriers"]),
            "affected": _as_int(row["affected"]),
            "unaffected": _as_int(row["unaffected"]),
            "uncertain": _as_int(row["uncertain"]),
        }
        fact: dict[str, Any] = {
            "n": len(bucket["facts"]) + 1,
            "variant": row["variant"] or "(unnamed)",
            "counts": {k: v for k, v in counts.items() if v is not None},
        }
        prov = _load_json(row["count_provenance"]) or {}
        if isinstance(prov, dict):
            role = _clip(prov.get("carriers_count_type"), 80)
            label = _clip(prov.get("carriers_column_label"), 120)
            if role:
                fact["count_role"] = role
            if label:
                fact["count_label"] = label
        loc = _clip(row["source_location"], 200)
        if loc:
            fact["source_location"] = loc
        quotes = _clip(row["key_quotes"])
        if quotes:
            fact["quote"] = quotes
        if sel_trust:
            tier = row["trust_tier"]
            if tier:
                fact["trust_tier"] = tier
            reasons = _load_json(row["trust_reasons"])
            if reasons:
                fact["trust_reasons"] = reasons
        bucket["facts"].append(fact)

    payloads = []
    for payload in grouped.values():
        payload["truncated"] = payload["n_facts_total"] > len(payload["facts"])
        payloads.append(payload)
    return payloads


def build_paper_check_prompt(payload: dict[str, Any]) -> str:
    return _PROMPT_TEMPLATE.format(
        payload=json.dumps(payload, ensure_ascii=False, indent=2)
    )


def _coerce_confidence(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        conf = float(value)
    except (TypeError, ValueError):
        return None
    return max(0.0, min(1.0, conf))


def normalize_result(raw: Any, payload: dict[str, Any]) -> dict[str, Any]:
    """Coerce the model's JSON into a stable record. Defensive: a missing or
    malformed verdict falls back to whether any flags were returned."""
    raw = raw if isinstance(raw, dict) else {}

    flags: list[dict[str, Any]] = []
    for item in raw.get("flags") or []:
        if not isinstance(item, dict):
            continue
        severity = str(item.get("severity") or "").strip().lower()
        flags.append(
            {
                "variant": str(item.get("variant") or "").strip() or "(unspecified)",
                "issue": str(item.get("issue") or "").strip(),
                "severity": severity if severity in VALID_SEVERITIES else "medium",
            }
        )

    verdict = str(raw.get("paper_verdict") or "").strip().lower()
    if verdict not in VALID_VERDICTS:
        verdict = "flag" if flags else "ok"

    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "verdict": verdict,
        "confidence": _coerce_confidence(raw.get("confidence")),
        "n_facts": len(payload.get("facts") or []),
        "n_flagged": len(flags),
        "flags": flags,
        "summary": str(raw.get("summary") or "").strip(),
    }


def ensure_schema(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS paper_final_check (
            pmid             TEXT PRIMARY KEY,
            gene_symbol      TEXT,
            model            TEXT,
            reasoning_effort TEXT,
            prompt_version   TEXT,
            check_version    TEXT,
            verdict          TEXT,
            confidence       REAL,
            n_facts          INTEGER,
            n_flagged        INTEGER,
            flags_json       TEXT,
            summary          TEXT,
            checked_at       TEXT
        )
        """
    )


def _record_row(
    conn: sqlite3.Connection,
    *,
    result: dict[str, Any],
    model: str,
    reasoning_effort: Optional[str],
    version: str,
    checked_at: str,
) -> None:
    conn.execute(
        """
        INSERT OR REPLACE INTO paper_final_check (
            pmid, gene_symbol, model, reasoning_effort, prompt_version,
            check_version, verdict, confidence, n_facts, n_flagged, flags_json,
            summary, checked_at
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        (
            result.get("pmid"),
            result.get("gene_symbol"),
            model,
            reasoning_effort,
            PROMPT_VERSION,
            version,
            result.get("verdict"),
            result.get("confidence"),
            result.get("n_facts"),
            result.get("n_flagged"),
            json.dumps(result.get("flags") or []),
            result.get("summary"),
            checked_at,
        ),
    )


class PaperFinalChecker:
    """LLM wrapper for the per-paper sniff test (thin BaseLLMCaller adapter)."""

    def __init__(
        self,
        model: str = DEFAULT_MODEL,
        reasoning_effort: Optional[str] = DEFAULT_REASONING_EFFORT,
        max_tokens: int = DEFAULT_MAX_TOKENS,
        temperature: float = 0.0,
    ):
        # Imported lazily so the pure helpers above don't drag in litellm.
        from utils.llm_utils import BaseLLMCaller, clamp_max_tokens

        self._caller = BaseLLMCaller(
            model=model,
            temperature=temperature,
            max_tokens=clamp_max_tokens(model, max_tokens),
            reasoning_effort=reasoning_effort,
        )

    def check(self, payload: dict[str, Any]) -> dict[str, Any]:
        raw = self._caller.call_llm_json(build_paper_check_prompt(payload))
        return normalize_result(raw, payload)


def apply_paper_final_check(
    db_path: str | Path,
    *,
    model: Optional[str] = None,
    reasoning_effort: Optional[str] = None,
    max_tokens: int = DEFAULT_MAX_TOKENS,
    pmids: Optional[list[str]] = None,
    limit: Optional[int] = None,
    checker: Optional[_Checker] = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Run the sniff test over every paper with extracted counts and record a
    soft verdict in ``paper_final_check``. Idempotent (keyed by pmid). Never
    edits or deletes a raw count.

    ``model`` / ``reasoning_effort`` default to ``config.settings`` when unset.
    Pass ``checker`` to inject a stub (tests); pass ``dry_run=True`` to build and
    return payloads without any LLM call.
    """
    db_path = str(db_path)
    if model is None or reasoning_effort is None:
        from config.settings import get_settings

        settings = get_settings()
        model = model or settings.paper_final_check_model
        if reasoning_effort is None:
            reasoning_effort = settings.paper_final_check_reasoning_effort

    version = check_version(model, reasoning_effort)
    checked_at = datetime.now(timezone.utc).isoformat()
    stats: dict[str, Any] = {
        "papers": 0,
        "checked": 0,
        "ok": 0,
        "flag": 0,
        "error": 0,
        "flagged_facts": 0,
        "model": model,
        "reasoning_effort": reasoning_effort,
        "version": version,
        "dry_run": dry_run,
    }

    conn = sqlite3.connect(db_path)
    try:
        ensure_schema(conn)
        payloads = gather_paper_payloads(conn, pmids=pmids)
        if limit is not None:
            payloads = payloads[:limit]
        stats["papers"] = len(payloads)

        if dry_run:
            stats["payloads"] = payloads
            return stats

        if checker is None:
            checker = PaperFinalChecker(
                model=model, reasoning_effort=reasoning_effort, max_tokens=max_tokens
            )

        for payload in payloads:
            try:
                result = checker.check(payload)
                verdict = result.get("verdict") or "ok"
                stats[verdict] = stats.get(verdict, 0) + 1
                stats["flagged_facts"] += int(result.get("n_flagged") or 0)
            except Exception as exc:  # noqa: BLE001 - one bad paper must not abort
                logger.warning(
                    "paper final check failed for pmid=%s: %s", payload.get("pmid"), exc
                )
                stats["error"] += 1
                result = {
                    "pmid": payload.get("pmid"),
                    "gene_symbol": payload.get("gene"),
                    "verdict": "error",
                    "confidence": None,
                    "n_facts": len(payload.get("facts") or []),
                    "n_flagged": 0,
                    "flags": [],
                    "summary": f"final check error: {exc}",
                }
            stats["checked"] += 1
            _record_row(
                conn,
                result=result,
                model=model,
                reasoning_effort=reasoning_effort,
                version=version,
                checked_at=checked_at,
            )
        conn.commit()
        logger.info(
            "paper final check: %d papers, ok=%d flag=%d error=%d (%s, effort=%s) [%s]",
            stats["papers"],
            stats["ok"],
            stats["flag"],
            stats["error"],
            model,
            reasoning_effort,
            version,
        )
        return stats
    finally:
        conn.close()
