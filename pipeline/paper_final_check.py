"""Per-paper "final check" — a gpt-5.6-sol pass over each paper.

Autonomy QC. After extraction, migration, and the gold-free structural trust
gate (see :mod:`pipeline.trust_gate`), this step asks a strong reasoning model
(default ``azure_ai/gpt-5.6-sol`` at ``xhigh`` effort) to review **each paper**.

Two modes, one step:

* **Source-grounded (default when the paper's source text is on disk):** the
  model reads the paper's actual source text PLUS the pipeline's extracted rows
  and returns, in one call, (a) a task-focused study/cohort summary, (b) an
  enumerated ``carrier_groups`` superset — every carrier group the paper reports,
  each tagged ``in_extraction`` and grounded by a verbatim quote — from which we
  DERIVE the carriers the pipeline MISSED, and (c) the trust-verdict/flags. This
  is the completeness signal: "what carriers does the paper report that we did
  not capture." Anti-parrot by construction — misses are derived from a required
  per-group flag, never asked as "what did you miss?".
* **DB-only fallback (no source text available):** the original "sniff test" —
  the model eyeballs the extracted rows + their captured provenance and flags
  count-role errors. No completeness signal is possible without the source.

It is a **judgment layer, not a mutation**: it RECORDS soft results in
``paper_final_check`` (+ the ``paper_carrier_groups`` child table) and never
edits or deletes a raw count. Pure helpers (``resolve_source_text``,
``_select_source_for_summary``, ``build_*_prompt``, ``normalize_*``) have no
LLM/DB coupling and are unit-tested directly; ``apply_paper_final_check`` is the
thin SQLite + LLM adapter and accepts injected checkers so tests never hit the
network.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Protocol

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]

# Defaults. The mainline (gvf-run) resolves the real model/effort from
# ``config.settings`` (PAPER_FINAL_CHECK_MODEL / _REASONING_EFFORT); these are
# the standalone fallbacks so the module is usable on its own.
DEFAULT_MODEL = "azure_ai/gpt-5.6-sol"
DEFAULT_REASONING_EFFORT = "xhigh"
# xhigh spends hidden reasoning tokens before the visible JSON, so give the
# completion generous headroom (clamp_max_tokens caps it per model).
DEFAULT_MAX_TOKENS = 16000
# The source-grounded summary reasons over full text at xhigh, so it needs far
# more output budget (reasoning tokens count against it) than the DB-only sniff
# test — otherwise the JSON truncates to empty on long papers.
SUMMARY_MAX_TOKENS = 64000
# Cap facts per paper so a 70-row paper stays a "quick" pass.
MAX_FACTS_PER_PAPER = 60
# Bound each provenance/quote field so the prompt stays compact.
MAX_QUOTE_CHARS = 800
# Cap enumerated carrier groups so output tokens stay bounded.
MAX_CARRIER_GROUPS = 80
# Default source-text budget fed to the summary prompt (~15K input tokens).
DEFAULT_MAX_SOURCE_CHARS = 60000
# FULL_CONTEXT files above this drop to the CLEANED sibling before truncation,
# so a supplement-folded multi-MB monster never loads whole.
OVERSIZE_SOURCE_BYTES = 3_000_000
# Minimum usable source length (chars); below this, treat as no source.
MIN_USABLE_SOURCE_CHARS = 200

PROMPT_VERSION = "pfc1"  # DB-only sniff test (stable)
SUMMARY_PROMPT_VERSION = "pfs1"  # source-grounded summary + completeness
SOURCE_RECIPE_VERSION = "src1"  # bump when resolve/select logic changes

VALID_VERDICTS = ("ok", "flag")
VALID_SEVERITIES = ("low", "medium", "high")
VALID_COMPLETENESS = ("complete", "gaps", "unsure")

# Preferred full-text file suffixes, most-complete first (FULL_CONTEXT carries
# folded supplements; CLEANED is noise-stripped; DATA_ZONES is scout-narrowed).
_FULLTEXT_SUFFIXES = ("_FULL_CONTEXT.md", "_CLEANED.md", "_DATA_ZONES.md")

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

_SUMMARY_PROMPT_TEMPLATE = """You are a task-focused COMPLETENESS auditor AND quality reviewer for an
automated pipeline whose job is to capture EVERY carrier of {gene} variants
described in a paper, each carrier's phenotype (affected vs unaffected, and which
disease), and where that carrier's cohort came from. You are AUDITING against the
paper's real text — not re-typing the whole paper.

You are given TWO blocks:

=== SOURCE TEXT (the paper; may include folded supplement tables) ===
{source}
=== END SOURCE TEXT ===

=== ALREADY CAPTURED BY THE PIPELINE (rows we already extracted) ===
{payload}
=== END ALREADY CAPTURED ===

Do ALL of the following and return ONE strict JSON object.

TASK A — Enumerate carrier_groups: list EVERY distinct carrier group present in
the SOURCE TEXT above. A "group" is a variant together with the people reported
to carry it: its total carriers, the affected vs unaffected split, the disease/
phenotype, and the cohort it came from. Scan the body text, ALL tables, folded
supplement tables (look for markers like "# SUPPLEMENTAL FILE" or
"GVF_FOLDED_SUPPLEMENTS"), and pedigrees — supplement and per-variant multi-row
tables are where carriers are most often reported. For EACH group:
  - copy a VERBATIM evidence_quote from the SOURCE above (never paraphrase),
  - set in_extraction=true ONLY if that group clearly matches one of the
    ALREADY CAPTURED rows (notation may differ — you decide the match),
  - otherwise in_extraction=false. Groups with in_extraction=false are carriers
    the pipeline MISSED — that is the point of this audit.
Never invent a carrier you cannot quote from the SOURCE. If the source is
truncated or is only an abstract, do not claim coverage beyond the text shown.

TASK B — study_summary (2-4 sentences: what the study is and which disease(s))
and cohort_sources (where the patients come from: a registry, a single family,
a biobank, a referral clinic, a cited prior cohort, etc.).

TASK C — trust flags: flag any ALREADY CAPTURED row that is SUSPECT (a cohort/
study denominator, a study-wide or mutation-class total copied onto one variant,
a gnomAD/allele-frequency figure, affected+unaffected exceeding the total, or an
implausible magnitude), cross-checked against the real text.

Return STRICT JSON only, no prose outside the object:
{{
  "paper_verdict": "ok" | "flag",
  "confidence": <0..1 confidence that the pipeline's captured rows are correct>,
  "study_summary": "<2-4 sentences>",
  "cohort_sources": ["<where carriers came from>", "..."],
  "carrier_groups": [
    {{
      "variant": "<as written in the paper>",
      "total_carriers": <int or null>,
      "affected": <int or null>,
      "unaffected": <int or null>,
      "disease": "<phenotype/disease for THIS group>",
      "cohort_source": "<which cohort this group came from>",
      "location_in_paper": "<Table 2 | Suppl Table S3 | Results text | Figure 1 pedigree>",
      "evidence_quote": "<VERBATIM snippet copied from the SOURCE above>",
      "in_extraction": true | false
    }}
  ],
  "completeness": {{
    "status": "complete" | "gaps" | "unsure",
    "notes": "<why; call out supplement/table gaps>"
  }},
  "flags": [
    {{"variant": "<row>", "issue": "<why suspect>", "severity": "low"|"medium"|"high"}}
  ]
}}"""


class _Checker(Protocol):
    def check(self, payload: dict[str, Any]) -> dict[str, Any]:  # pragma: no cover
        ...


class _SummaryChecker(Protocol):
    def check(  # pragma: no cover
        self, payload: dict[str, Any], source_text: str
    ) -> dict[str, Any]: ...


# ---------------------------------------------------------------------------
# Versioning
# ---------------------------------------------------------------------------


def check_version(model: str, reasoning_effort: Optional[str]) -> str:
    """Content hash pinning the DB-only sniff-test generation for a stored row."""
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


def summary_check_version(
    model: str, reasoning_effort: Optional[str], max_source_chars: int
) -> str:
    """Content hash pinning the source-grounded summary generation (folds the
    prompt, the source-recipe version, and the source budget), so a fleet re-run
    only re-pays for papers whose reviewer generation actually changed."""
    payload = json.dumps(
        {
            "prompt_version": SUMMARY_PROMPT_VERSION,
            "recipe_version": SOURCE_RECIPE_VERSION,
            "template": _SUMMARY_PROMPT_TEMPLATE,
            "model": model,
            "reasoning_effort": reasoning_effort,
            "max_source_chars": max_source_chars,
        },
        sort_keys=True,
    )
    return (
        f"{SUMMARY_PROMPT_VERSION}-"
        + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
    )


# ---------------------------------------------------------------------------
# Small pure helpers
# ---------------------------------------------------------------------------


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


def _norm_ws(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip().casefold()


def _normalized_contains(haystack: str, needle: str) -> bool:
    """Whitespace/case-insensitive containment — robust to markdown pipes, OCR
    spacing, and casing so a real quote isn't rejected on cosmetics."""
    needle_n = _norm_ws(needle)
    if len(needle_n) < 8:  # too short to verify meaningfully
        return False
    return needle_n in _norm_ws(haystack)


# ---------------------------------------------------------------------------
# Source-text resolution (filesystem only; unit-tested with tmp_path)
# ---------------------------------------------------------------------------


def _is_usable_source(text: Optional[str]) -> bool:
    return bool(text and len(text.strip()) >= MIN_USABLE_SOURCE_CHARS)


def _resolve_path(raw: str) -> Optional[Path]:
    path = Path(raw)
    if path.is_absolute():
        return path if path.exists() else None
    for base in (Path.cwd(), REPO_ROOT):
        cand = base / path
        if cand.exists():
            return cand
    return None


def _read_text(path: Path) -> Optional[str]:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _pick_fulltext(dirpath: Path, pmid: str, oversize_bytes: int) -> Optional[Path]:
    """Pick the most-complete usable full-text sibling for a pmid in a dir.

    Prefer FULL_CONTEXT (folded supplements) unless it is oversized and a smaller
    sibling exists; among non-oversized candidates pick the largest.
    """
    cands: list[tuple[str, Path, int]] = []
    for suffix in _FULLTEXT_SUFFIXES:
        f = dirpath / f"{pmid}{suffix}"
        try:
            if f.exists():
                cands.append((suffix, f, f.stat().st_size))
        except OSError:
            continue
    if not cands:
        return None
    full = next((c for c in cands if c[0] == "_FULL_CONTEXT.md"), None)
    if full and full[2] <= oversize_bytes:
        return full[1]
    under = [c for c in cands if c[2] <= oversize_bytes]
    if under:
        return max(under, key=lambda c: c[2])[1]
    return min(cands, key=lambda c: c[2])[1]


def _abstract_from_json(path: Path) -> Optional[str]:
    data = _load_json(_read_text(path))
    if isinstance(data, dict):
        abstract = str(data.get("abstract") or "").strip()
        return abstract or None
    return None


def resolve_source_text(
    pmid: str,
    gene: Optional[str],
    run_dir: Optional[str | Path] = None,
    source_file_hint: Optional[str] = None,
    oversize_bytes: int = OVERSIZE_SOURCE_BYTES,
) -> tuple[Optional[str], str, Optional[str]]:
    """Return (text, source_kind, source_path) for a paper's on-disk source.

    ``source_kind`` is ``"fulltext"`` | ``"abstract_only"`` | ``"none"``.
    Precedence: (a) the extraction_metadata.source_file pointer (the exact file
    extraction saw); (b) run_dir/pmc_fulltext; (c) corpus/<GENE>/<PMID>/;
    (d) run_dir/abstract_json. Never raises.
    """
    pmid = str(pmid)
    gene = (gene or "").strip()
    run_path = Path(run_dir) if run_dir else None

    # (a) extraction_metadata.source_file hint
    if source_file_hint:
        resolved = _resolve_path(source_file_hint)
        if resolved is not None:
            if resolved.suffix == ".json":
                abstract = _abstract_from_json(resolved)
                if abstract:
                    return abstract, "abstract_only", str(resolved)
            else:
                best = _pick_fulltext(resolved.parent, pmid, oversize_bytes)
                if best is not None:
                    text = _read_text(best)
                    if _is_usable_source(text):
                        return text, "fulltext", str(best)

    # (b) run_dir/pmc_fulltext
    if run_path is not None:
        best = _pick_fulltext(run_path / "pmc_fulltext", pmid, oversize_bytes)
        if best is not None:
            text = _read_text(best)
            if _is_usable_source(text):
                return text, "fulltext", str(best)

    # (c) corpus fallback (cross-run superset)
    if gene:
        best = _pick_fulltext(REPO_ROOT / "corpus" / gene / pmid, pmid, oversize_bytes)
        if best is not None:
            text = _read_text(best)
            if _is_usable_source(text):
                return text, "fulltext", str(best)

    # (d) abstract_json in run_dir
    if run_path is not None:
        aj = run_path / "abstract_json" / f"{pmid}.json"
        if aj.exists():
            abstract = _abstract_from_json(aj)
            if abstract:
                return abstract, "abstract_only", str(aj)

    return None, "none", None


# Sections that never contain carriers — dropped before the char budget.
_DEAD_SECTION_RE = re.compile(
    r"^#{1,6}\s*(references|bibliography|acknowledg|author\s+contributions?|"
    r"funding|conflicts?\s+of\s+interest|competing\s+interests?|"
    r"declaration|data\s+availability)",
    re.IGNORECASE,
)
_HEADER_RE = re.compile(r"^#{1,6}\s+\S")
_TABLEISH_RE = re.compile(
    r"^\s*\||supplemental?\s+file|supplementary\s+table|GVF_FOLDED", re.IGNORECASE
)


def _select_source_for_summary(
    text: str, max_chars: int = DEFAULT_MAX_SOURCE_CHARS
) -> tuple[str, bool]:
    """Trim source to the completeness-relevant text within a char budget.

    Drops dead sections (references/acknowledgments/funding/...); then, if still
    over budget, keeps a prose head plus ALL table/supplement-bearing lines
    (that is where carriers concentrate) and truncates the long prose middle.
    Returns (selected_text, truncated) — truncated is True only when the char
    budget forced dropping potentially-relevant content, so a gap under
    truncation is auditable and re-runnable at a higher cap.
    """
    lines = text.split("\n")
    kept: list[str] = []
    skipping = False
    for ln in lines:
        if _HEADER_RE.match(ln):
            skipping = bool(_DEAD_SECTION_RE.match(ln))
        if not skipping:
            kept.append(ln)
    body = "\n".join(kept)
    if len(body) <= max_chars:
        return body, False

    table_lines = [ln for ln in kept if _TABLEISH_RE.search(ln)]
    head_budget = int(max_chars * 0.6)
    head = body[:head_budget]
    tail_budget = max_chars - len(head)
    tail = "\n".join(table_lines)[: max(0, tail_budget)]
    selected = (
        head
        + "\n\n[... prose truncated; table/supplement rows preserved below ...]\n\n"
        + tail
    )
    return selected[: max_chars + 256], True


# ---------------------------------------------------------------------------
# Payload gathering
# ---------------------------------------------------------------------------


def _penetrance_has_trust(conn: sqlite3.Connection) -> bool:
    cols = {row[1] for row in conn.execute("PRAGMA table_info(penetrance_data)")}
    return "trust_tier" in cols


def gather_paper_payloads(
    conn: sqlite3.Connection,
    pmids: Optional[list[str]] = None,
    max_facts: int = MAX_FACTS_PER_PAPER,
) -> list[dict[str, Any]]:
    """Assemble one payload per paper that has extracted counts (the pipeline's
    already-captured rows + their provenance). Papers with no ``penetrance_data``
    rows are skipped — there is nothing to check."""
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


def _source_hints(conn: sqlite3.Connection) -> dict[str, dict[str, Any]]:
    """Map pmid -> {source_file, source_type, abstract_only} from the DB's
    extraction_metadata pointer (latest row per pmid). Empty when unavailable."""
    tables = {
        r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
    }
    if "extraction_metadata" not in tables:
        return {}
    cols = {r[1] for r in conn.execute("PRAGMA table_info(extraction_metadata)")}
    if "source_file" not in cols:
        return {}
    order_col = "extraction_id" if "extraction_id" in cols else "rowid"
    select_cols = ["pmid", "source_file"]
    for c in ("source_type", "abstract_only"):
        if c in cols:
            select_cols.append(c)
    sql = (
        f"SELECT {', '.join(select_cols)}, MAX({order_col}) "
        f"FROM extraction_metadata GROUP BY pmid"
    )
    hints: dict[str, dict[str, Any]] = {}
    for row in conn.execute(sql):
        record = dict(zip(select_cols, row))
        pmid = record.get("pmid")
        if pmid is not None:
            hints[str(pmid)] = record
    return hints


# ---------------------------------------------------------------------------
# Prompt building + result normalization
# ---------------------------------------------------------------------------


def build_paper_check_prompt(payload: dict[str, Any]) -> str:
    return _PROMPT_TEMPLATE.format(
        payload=json.dumps(payload, ensure_ascii=False, indent=2)
    )


def build_paper_summary_prompt(payload: dict[str, Any], source_text: str) -> str:
    return _SUMMARY_PROMPT_TEMPLATE.format(
        gene=payload.get("gene") or "the target gene",
        source=source_text,
        payload=json.dumps(payload, ensure_ascii=False, indent=2),
    )


def _coerce_confidence(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        conf = float(value)
    except (TypeError, ValueError):
        return None
    return max(0.0, min(1.0, conf))


def _norm_flags(raw: dict[str, Any]) -> tuple[list[dict[str, Any]], str]:
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
    return flags, verdict


def normalize_result(raw: Any, payload: dict[str, Any]) -> dict[str, Any]:
    """Coerce the DB-only sniff-test JSON into a stable record."""
    raw = raw if isinstance(raw, dict) else {}
    flags, verdict = _norm_flags(raw)
    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "source_grounded": False,
        "verdict": verdict,
        "confidence": _coerce_confidence(raw.get("confidence")),
        "n_facts": len(payload.get("facts") or []),
        "n_flagged": len(flags),
        "flags": flags,
        "summary": str(raw.get("summary") or "").strip(),
    }


def normalize_summary_result(
    raw: Any, payload: dict[str, Any], source_text: str
) -> dict[str, Any]:
    """Coerce the source-grounded JSON into a stable record + derive the
    missed-carrier signal (never asked directly; derived from in_extraction)."""
    raw = raw if isinstance(raw, dict) else {}
    flags, verdict = _norm_flags(raw)

    groups: list[dict[str, Any]] = []
    for g in (raw.get("carrier_groups") or [])[:MAX_CARRIER_GROUPS]:
        if not isinstance(g, dict):
            continue
        in_extraction = bool(g.get("in_extraction"))
        quote = str(g.get("evidence_quote") or "").strip()
        groups.append(
            {
                "variant": str(g.get("variant") or "").strip() or "(unspecified)",
                "total_carriers": _as_int(g.get("total_carriers")),
                "affected": _as_int(g.get("affected")),
                "unaffected": _as_int(g.get("unaffected")),
                "disease": str(g.get("disease") or "").strip(),
                "cohort_source": str(g.get("cohort_source") or "").strip(),
                "location_in_paper": str(g.get("location_in_paper") or "").strip(),
                "evidence_quote": quote[:MAX_QUOTE_CHARS],
                "in_extraction": in_extraction,
                # Soft: keep-with-flag, never drop, so a real miss isn't lost to
                # markdown/whitespace/OCR cosmetics.
                "quote_verified": 1 if _normalized_contains(source_text, quote) else 0,
                "status": "reported_extracted" if in_extraction else "reported_missing",
            }
        )

    n_reported = len(groups)
    n_in = sum(1 for g in groups if g["in_extraction"])
    n_missing = n_reported - n_in

    comp = raw.get("completeness")
    comp = comp if isinstance(comp, dict) else {}
    status = str(comp.get("status") or "").strip().lower()
    if status not in VALID_COMPLETENESS:
        status = "gaps" if n_missing > 0 else ("complete" if n_reported else "unsure")

    cohort_sources = [
        str(c).strip() for c in (raw.get("cohort_sources") or []) if str(c).strip()
    ][:20]
    summary = str(raw.get("study_summary") or raw.get("summary") or "").strip()

    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "source_grounded": True,
        "verdict": verdict,
        "confidence": _coerce_confidence(raw.get("confidence")),
        "n_facts": len(payload.get("facts") or []),
        "n_flagged": len(flags),
        "flags": flags,
        "summary": summary,
        "cohort_sources": cohort_sources,
        "carrier_groups": groups,
        "completeness": {
            "status": status,
            "n_reported": n_reported,
            "n_in_extraction": n_in,
            "n_missing": n_missing,
            "notes": str(comp.get("notes") or "").strip(),
        },
    }


# ---------------------------------------------------------------------------
# Schema + persistence
# ---------------------------------------------------------------------------

_SUMMARY_COLUMNS = (
    ("source_grounded", "INTEGER"),
    ("source_kind", "TEXT"),
    ("source_file", "TEXT"),
    ("source_chars", "INTEGER"),
    ("source_truncated", "INTEGER"),
    ("summary_json", "TEXT"),
    ("cohort_sources", "TEXT"),
    ("completeness_status", "TEXT"),
    ("n_reported", "INTEGER"),
    ("n_in_extraction", "INTEGER"),
    ("n_missing", "INTEGER"),
)


def ensure_schema(conn: sqlite3.Connection) -> None:
    """Create the base paper_final_check table (DB-only sniff test)."""
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


def ensure_summary_schema(conn: sqlite3.Connection) -> None:
    """Widen paper_final_check with the summary columns (idempotent ALTERs) and
    create the paper_carrier_groups child table + indexes."""
    ensure_schema(conn)
    existing = {row[1] for row in conn.execute("PRAGMA table_info(paper_final_check)")}
    for name, coltype in _SUMMARY_COLUMNS:
        if name not in existing:
            conn.execute(f"ALTER TABLE paper_final_check ADD COLUMN {name} {coltype}")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS paper_carrier_groups (
            id                INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid              TEXT,
            gene_symbol       TEXT,
            variant           TEXT,
            total_carriers    INTEGER,
            affected          INTEGER,
            unaffected        INTEGER,
            disease           TEXT,
            cohort_source     TEXT,
            location_in_paper TEXT,
            evidence_quote    TEXT,
            quote_verified    INTEGER,
            in_extraction     INTEGER,
            status            TEXT,
            check_version     TEXT,
            checked_at        TEXT
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pcg_pmid ON paper_carrier_groups(pmid)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pcg_status ON paper_carrier_groups(status)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pcg_gene_status "
        "ON paper_carrier_groups(gene_symbol, status)"
    )


def _record_row(
    conn: sqlite3.Connection,
    *,
    result: dict[str, Any],
    model: str,
    reasoning_effort: Optional[str],
    prompt_version: str,
    version: str,
    checked_at: str,
    source_kind: str = "none",
    source_file: Optional[str] = None,
    source_chars: Optional[int] = None,
    source_truncated: Optional[int] = None,
) -> None:
    """Upsert one paper_final_check row (works for DB-only and source-grounded)."""
    comp = result.get("completeness") or {}
    conn.execute(
        """
        INSERT OR REPLACE INTO paper_final_check (
            pmid, gene_symbol, model, reasoning_effort, prompt_version,
            check_version, verdict, confidence, n_facts, n_flagged, flags_json,
            summary, checked_at, source_grounded, source_kind, source_file,
            source_chars, source_truncated, summary_json, cohort_sources,
            completeness_status, n_reported, n_in_extraction, n_missing
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        (
            result.get("pmid"),
            result.get("gene_symbol"),
            model,
            reasoning_effort,
            prompt_version,
            version,
            result.get("verdict"),
            result.get("confidence"),
            result.get("n_facts"),
            result.get("n_flagged"),
            json.dumps(result.get("flags") or []),
            result.get("summary"),
            checked_at,
            1 if result.get("source_grounded") else 0,
            source_kind,
            source_file,
            source_chars,
            source_truncated,
            json.dumps(result) if result.get("source_grounded") else None,
            json.dumps(result.get("cohort_sources") or [])
            if result.get("source_grounded")
            else None,
            comp.get("status") if result.get("source_grounded") else None,
            comp.get("n_reported") if result.get("source_grounded") else None,
            comp.get("n_in_extraction") if result.get("source_grounded") else None,
            comp.get("n_missing") if result.get("source_grounded") else None,
        ),
    )


def _record_carrier_groups(
    conn: sqlite3.Connection,
    *,
    pmid: Any,
    gene_symbol: Any,
    groups: list[dict[str, Any]],
    version: str,
    checked_at: str,
) -> None:
    """Replace this paper's carrier-group rows (DELETE-before-insert idempotency,
    in the caller's transaction)."""
    conn.execute("DELETE FROM paper_carrier_groups WHERE pmid = ?", (pmid,))
    conn.executemany(
        """
        INSERT INTO paper_carrier_groups (
            pmid, gene_symbol, variant, total_carriers, affected, unaffected,
            disease, cohort_source, location_in_paper, evidence_quote,
            quote_verified, in_extraction, status, check_version, checked_at
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        [
            (
                pmid,
                gene_symbol,
                g.get("variant"),
                g.get("total_carriers"),
                g.get("affected"),
                g.get("unaffected"),
                g.get("disease"),
                g.get("cohort_source"),
                g.get("location_in_paper"),
                g.get("evidence_quote"),
                g.get("quote_verified"),
                1 if g.get("in_extraction") else 0,
                g.get("status"),
                version,
                checked_at,
            )
            for g in groups
        ],
    )


def _existing_check_version(conn: sqlite3.Connection, pmid: str) -> Optional[str]:
    row = conn.execute(
        "SELECT check_version FROM paper_final_check WHERE pmid = ?", (pmid,)
    ).fetchone()
    return row[0] if row else None


# ---------------------------------------------------------------------------
# LLM adapters
# ---------------------------------------------------------------------------


class PaperFinalChecker:
    """DB-only sniff test (thin BaseLLMCaller adapter)."""

    def __init__(
        self,
        model: str = DEFAULT_MODEL,
        reasoning_effort: Optional[str] = DEFAULT_REASONING_EFFORT,
        max_tokens: int = DEFAULT_MAX_TOKENS,
        temperature: float = 0.0,
    ):
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


class PaperSummaryChecker:
    """Source-grounded summary + completeness (thin BaseLLMCaller adapter)."""

    def __init__(
        self,
        model: str = DEFAULT_MODEL,
        reasoning_effort: Optional[str] = DEFAULT_REASONING_EFFORT,
        max_tokens: int = SUMMARY_MAX_TOKENS,
        temperature: float = 0.0,
    ):
        from utils.llm_utils import BaseLLMCaller, clamp_max_tokens

        self._caller = BaseLLMCaller(
            model=model,
            temperature=temperature,
            max_tokens=clamp_max_tokens(model, max_tokens),
            reasoning_effort=reasoning_effort,
        )

    def check(self, payload: dict[str, Any], source_text: str) -> dict[str, Any]:
        raw = self._caller.call_llm_json(
            build_paper_summary_prompt(payload, source_text)
        )
        return normalize_summary_result(raw, payload, source_text)


def _error_result(payload: dict[str, Any], exc: Exception) -> dict[str, Any]:
    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "source_grounded": False,
        "verdict": "error",
        "confidence": None,
        "n_facts": len(payload.get("facts") or []),
        "n_flagged": 0,
        "flags": [],
        "summary": f"final check error: {exc}",
    }


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def apply_paper_final_check(
    db_path: str | Path,
    *,
    model: Optional[str] = None,
    reasoning_effort: Optional[str] = None,
    max_tokens: int = DEFAULT_MAX_TOKENS,
    run_dir: Optional[str | Path] = None,
    gene: Optional[str] = None,
    source_grounded: Optional[bool] = None,
    max_source_chars: Optional[int] = None,
    pmids: Optional[list[str]] = None,
    limit: Optional[int] = None,
    checker: Optional[_Checker] = None,
    summary_checker: Optional[_SummaryChecker] = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Review every paper with extracted counts and record a soft result.

    Source-grounded when the paper's source text is on disk (``source_grounded``
    default from settings): one gpt-5.6-sol/xhigh call produces the sniff-test
    verdict + a study/cohort summary + an enumerated ``carrier_groups`` superset
    (each tagged ``in_extraction``, grounded by a verbatim quote) from which the
    MISSED carriers are derived into ``paper_carrier_groups`` (status=
    ``reported_missing``). Falls back to the DB-only sniff test when no source is
    available. Idempotent (keyed by pmid + version-skip). Never mutates counts.
    """
    db_path = str(db_path)
    if (
        model is None
        or reasoning_effort is None
        or source_grounded is None
        or max_source_chars is None
    ):
        from config.settings import get_settings

        settings = get_settings()
        model = model or settings.paper_final_check_model
        if reasoning_effort is None:
            reasoning_effort = settings.paper_final_check_reasoning_effort
        if source_grounded is None:
            source_grounded = getattr(settings, "paper_summary_source_grounded", True)
        if max_source_chars is None:
            max_source_chars = getattr(
                settings, "paper_summary_max_source_chars", DEFAULT_MAX_SOURCE_CHARS
            )

    db_version = check_version(model, reasoning_effort)
    sum_version = summary_check_version(model, reasoning_effort, max_source_chars)
    checked_at = datetime.now(timezone.utc).isoformat()
    stats: dict[str, Any] = {
        "papers": 0,
        "checked": 0,
        "skipped": 0,
        "ok": 0,
        "flag": 0,
        "error": 0,
        "source_grounded": 0,
        "source_absent": 0,
        "flagged_facts": 0,
        "missing_carriers": 0,
        "model": model,
        "reasoning_effort": reasoning_effort,
        "db_version": db_version,
        "summary_version": sum_version,
        "dry_run": dry_run,
    }

    conn = sqlite3.connect(db_path)
    try:
        ensure_summary_schema(conn)
        payloads = gather_paper_payloads(conn, pmids=pmids)
        if limit is not None:
            payloads = payloads[:limit]
        stats["papers"] = len(payloads)
        hints = _source_hints(conn)

        if dry_run:
            stats["payloads"] = payloads
            return stats

        for payload in payloads:
            pmid = str(payload.get("pmid"))
            paper_gene = gene or payload.get("gene")

            # Resolve source (cheap: path stat + a read) to decide the mode.
            source_text, source_kind, source_path = (None, "none", None)
            if source_grounded:
                source_text, source_kind, source_path = resolve_source_text(
                    pmid,
                    paper_gene,
                    run_dir=run_dir,
                    source_file_hint=(hints.get(pmid) or {}).get("source_file"),
                )
            use_summary = bool(
                source_grounded and source_kind == "fulltext" and source_text
            )
            intended_version = sum_version if use_summary else db_version

            # Version-skip: identical reviewer generation already recorded.
            if _existing_check_version(conn, pmid) == intended_version:
                stats["skipped"] += 1
                continue

            if use_summary:
                selected, truncated = _select_source_for_summary(
                    source_text, max_source_chars
                )
                if summary_checker is None:
                    summary_checker = PaperSummaryChecker(
                        model=model,
                        reasoning_effort=reasoning_effort,
                        max_tokens=SUMMARY_MAX_TOKENS,
                    )
                try:
                    result = summary_checker.check(payload, selected)
                    stats["source_grounded"] += 1
                    comp = result.get("completeness") or {}
                    stats["missing_carriers"] += int(comp.get("n_missing") or 0)
                except Exception as exc:  # noqa: BLE001 - one paper must not abort
                    logger.warning("paper summary failed for pmid=%s: %s", pmid, exc)
                    result = _error_result(payload, exc)
                    use_summary = False  # record as an error row, no groups
                _record_row(
                    conn,
                    result=result,
                    model=model,
                    reasoning_effort=reasoning_effort,
                    prompt_version=SUMMARY_PROMPT_VERSION,
                    version=sum_version,
                    checked_at=checked_at,
                    source_kind=source_kind,
                    source_file=source_path,
                    source_chars=len(selected),
                    source_truncated=1 if truncated else 0,
                )
                if result.get("source_grounded"):
                    _record_carrier_groups(
                        conn,
                        pmid=pmid,
                        gene_symbol=result.get("gene_symbol"),
                        groups=result.get("carrier_groups") or [],
                        version=sum_version,
                        checked_at=checked_at,
                    )
            else:
                if source_kind == "none":
                    stats["source_absent"] += 1
                if checker is None:
                    checker = PaperFinalChecker(
                        model=model,
                        reasoning_effort=reasoning_effort,
                        max_tokens=max_tokens,
                    )
                try:
                    result = checker.check(payload)
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "paper final check failed for pmid=%s: %s", pmid, exc
                    )
                    result = _error_result(payload, exc)
                _record_row(
                    conn,
                    result=result,
                    model=model,
                    reasoning_effort=reasoning_effort,
                    prompt_version=PROMPT_VERSION,
                    version=db_version,
                    checked_at=checked_at,
                    source_kind=source_kind,
                    source_file=source_path,
                )

            verdict = result.get("verdict") or "ok"
            stats[verdict] = stats.get(verdict, 0) + 1
            stats["flagged_facts"] += int(result.get("n_flagged") or 0)
            stats["checked"] += 1

        conn.commit()
        logger.info(
            "paper final check: %d papers, checked=%d skipped=%d "
            "(grounded=%d, missing_carriers=%d) ok=%d flag=%d error=%d [%s / %s]",
            stats["papers"],
            stats["checked"],
            stats["skipped"],
            stats["source_grounded"],
            stats["missing_carriers"],
            stats["ok"],
            stats["flag"],
            stats["error"],
            db_version,
            sum_version,
        )
        return stats
    finally:
        conn.close()
