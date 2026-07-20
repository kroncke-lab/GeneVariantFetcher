"""Per-paper "final check" — a gpt-5.6-sol pass over each paper.

Autonomy QC. After extraction, migration, and the gold-free structural trust
gate (see :mod:`pipeline.trust_gate`), this step asks a strong reasoning model
(default ``azure_ai/gpt-5.6-sol`` at ``xhigh`` effort) to review **each paper**.

Two modes, one step:

* **Source-grounded (default when the paper's scouted source is on disk):** the
  model reads the paper's ABSTRACT + the Data-Scout ``DATA_ZONES`` (the
  scout-narrowed high-value regions, NOT the full text) PLUS the extracted rows
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

It is an **auditor, not a raw-data mutation layer**: it records structured
findings in ``paper_final_check`` (+ the ``paper_carrier_groups`` child table)
and never edits or deletes a raw count. The separate
``pipeline.paper_final_check_gate`` composer can bind source-grounded findings
to exact count facts and quarantine only the affected trusted projection. Pure
helpers (``resolve_summary_source``,
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
# Cap provenance-rich facts per paper so a 70-row paper stays a "quick" pass.
# Every captured row still appears in ``captured_fact_index``; this cap only
# bounds quotes/locations/trust detail.
MAX_FACTS_PER_PAPER = 60
# Bound each provenance/quote field so the prompt stays compact.
MAX_QUOTE_CHARS = 800
# Cap enumerated carrier groups so output tokens stay bounded.
MAX_CARRIER_GROUPS = 80
# Char budget for the abstract + scouted zones fed to the summary prompt.
DEFAULT_MAX_SOURCE_CHARS = 60000
# Minimum usable source length (chars); below this, treat as no source.
MIN_USABLE_SOURCE_CHARS = 200

PROMPT_VERSION = "pfc7"  # source-quoted fact/field flags + phenotype reason
SUMMARY_PROMPT_VERSION = "pfs12"  # attached table-title evidence + phenotype reason
SOURCE_RECIPE_VERSION = "src2"  # abstract + Data-Scout DATA_ZONES (no full text)

VALID_VERDICTS = ("ok", "flag")
VALID_SEVERITIES = ("low", "medium", "high")
VALID_COMPLETENESS = ("complete", "gaps", "unsure")
VALID_FLAG_FIELDS = ("total_carriers", "affected", "unaffected", "uncertain")
FLAG_FIELD_ALIASES = {
    "carriers": "total_carriers",
    "carrier": "total_carriers",
    "patient_count": "total_carriers",
    "total_carriers_observed": "total_carriers",
    "affected_count": "affected",
    "unaffected_count": "unaffected",
    "uncertain_count": "uncertain",
}
VALID_FLAG_REASONS = (
    "count_is_total",
    "population_count",
    "wrong_column",
    "arith_inconsistent",
    "phenotype_misclassified",
    "unsupported_count",
    "wrong_gene",
    "other",
)

_PROMPT_TEMPLATE = """You are the final quality-control reviewer (a "sniff test") for an automated
pipeline that extracts variant carrier counts from biomedical papers. You are
NOT re-extracting the paper. For ONE paper you are given a COMPLETE compact
``captured_fact_index`` containing every extracted identity/count row, plus a
possibly capped ``facts`` list with richer provenance (source location, verbatim
quote(s), and declared count role/label). Use the complete index when deciding
which rows exist; use the detailed facts when judging provenance. Judge quickly
but carefully whether each row looks trustworthy, and flag the ones that do not.

A row is SUSPECT if any of these apply:
- The count looks like a cohort/study denominator (everyone enrolled, screened,
  sampled, or sequenced), not the number of people who carry THIS variant.
- The count is a study-wide, family-set, or mutation-class total copied onto a
  single variant.
- The number is a population/allele figure (gnomAD, ExAC, TOPMed, 1000 Genomes,
  MAF, AC/AN, allele frequency), not a clinical carrier count.
- affected + unaffected exceeds total_carriers.
- The number is an assay replicate/cell count, or a table value that is really a
  prediction score, allele frequency, exon/domain number, genomic position, age
  at diagnosis, mean age, or a family-history disease/case count.
- A claimed count column does not line up with the target row in the supplied
  table header/row context (watch especially for blank leading/group columns).
- The count has no supporting quote or source location (treat as low confidence).
- The magnitude is implausible for the carriers of a single variant.

Every flag MUST name the exact ``fact_id`` value(s) copied from
``captured_fact_index`` and the specific count field(s) that are suspect. Never
blanket-flag every fact in a paper or rely on variant text alone.

Paper (JSON):
{payload}

Return STRICT JSON only, no prose outside the object:
{{
  "paper_verdict": "ok" | "flag",
  "confidence": <number between 0 and 1: your confidence the extraction for this
                 paper is correct overall>,
  "flags": [
    {{"fact_ids": [<one or more integer fact_id values copied exactly from the
                     suspect captured_fact_index rows>],
      "variant": "<the variant field from the suspect row(s)>",
      "fields": ["total_carriers"|"affected"|"unaffected"|"uncertain"],
      "reason_code": "count_is_total"|"population_count"|"wrong_column"|
                     "arith_inconsistent"|"phenotype_misclassified"|
                     "unsupported_count"|"wrong_gene"|"other",
      "evidence_quote": "<verbatim header/row text that demonstrates the problem>",
      "issue": "<short reason the named fields are suspect>",
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

WHY THIS PROJECT EXISTS (judge relevance against this): GVF builds an auditable
database of per-variant PATIENT carrier counts and their phenotypes so the
Kroncke Lab can estimate variant PENETRANCE for clinical interpretation. For each
{gene} variant we want exactly three things: how many INDIVIDUAL humans carry it,
how many are AFFECTED vs UNAFFECTED (and with which disease), and the cohort/
provenance. We do NOT want, and you should treat as the WRONG kind of data:
allele frequencies (gnomAD/ExAC/MAF), in-vitro / electrophysiology / functional
measurements, in-silico predictions (SIFT/PolyPhen/CADD), or study-wide/cohort
denominators — none of these are per-variant human carrier counts. A table is the
"right" table for us only if it reports individual human carriers of specific
variants; if a captured row came from a "wrong kind" table, flag it.

You are given TWO blocks:

=== SOURCE TEXT (the paper; may include folded supplement tables) ===
{source}
=== END SOURCE TEXT ===

=== ALREADY CAPTURED BY THE PIPELINE ===
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
  - set in_extraction=true ONLY if that group clearly matches an entry in the
    COMPLETE ``captured_fact_index`` (notation may differ — use every identifier
    provided and decide the match). The richer ``facts`` list may be capped and
    MUST NOT be treated as the complete inventory,
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
implausible magnitude). Also flag age-at-diagnosis, mean-age, and family-history
disease/case values that were misread as carrier or affected counts. Cross-check
the claimed column against any compact header + target-row evidence supplied.
Every actionable flag MUST name the exact ``fact_id`` value(s) from
``captured_fact_index`` and the specific count field(s) that are unsupported.
Never use a paper-wide or variant-text-only flag as a substitute for enumerating
the affected fact IDs. Do not flag a field merely because its evidence was
omitted from the capped rich ``facts`` list; the complete index is authoritative
for what was captured.

Important phenotype rule: a cohort-level inclusion criterion can support a
per-variant affected status. If the counted subjects are explicitly patients or
cases with the target disease, their per-variant patient count may legitimately
equal affected and unaffected may be zero even when the table does not repeat a
separate phenotype column. Flag this only when the source shows a mixed cohort,
controls/asymptomatic carriers, or another direct contradiction. Lack of a
repeated row-level affected/unaffected split is not by itself a high-severity
error.

Every trust flag MUST copy a concise verbatim ``evidence_quote`` from the SOURCE
TEXT that demonstrates the mismatch. For a table problem, include enough header
and target-row text to establish the column role. A flag may bind multiple facts
from one wrong-gene table; for that table-level problem, quote its title, header,
and one representative target row (do not append a second distant row). Use
``population_count`` for a
gnomAD/ExAC/allele-count column and ``wrong_column`` for age, score, or shifted
column values. Use ``phenotype_misclassified`` only when the source directly
contradicts the captured affected/unaffected status or split; absent row-level
phenotype detail alone is not enough. ``unsupported_count`` means weak/absent
support rather than a demonstrated contradiction; it is advisory and must not
be presented as a hard column error.

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
    {{"fact_ids": [<one or more integer fact_id values copied exactly from the
                     suspect captured_fact_index rows>],
      "variant": "<row>",
      "fields": ["total_carriers"|"affected"|"unaffected"|"uncertain"],
      "reason_code": "count_is_total"|"population_count"|"wrong_column"|
                     "arith_inconsistent"|"phenotype_misclassified"|
                     "unsupported_count"|"wrong_gene"|"other",
      "evidence_quote": "<VERBATIM source header/row demonstrating the mismatch>",
      "issue": "<why the named fields are suspect>",
      "severity": "low"|"medium"|"high"}}
  ]
}}"""


class _Checker(Protocol):
    def check(self, payload: dict[str, Any]) -> dict[str, Any]:  # pragma: no cover
        ...


class _SummaryChecker(Protocol):
    def check(  # pragma: no cover
        self, payload: dict[str, Any], source_text: str, truncated: bool = False
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


def _source_evidence_verified(source_text: str, quote: str) -> bool:
    """Verify a contiguous quote or fragments from one logical table.

    Table headers and target rows are often non-adjacent. Reviewers may return
    them in one multi-line evidence block. A table title may also be included,
    but it must be in the short preamble immediately attached to that same table.
    This permits title + sampled head/tail rows while preventing a model from
    stitching a header or title from one table to a target row in another.
    Markdown separator rows carry no independent evidence and are ignored.
    """
    if _normalized_contains(source_text, quote):
        return True
    fragments = []
    for raw_line in str(quote or "").splitlines():
        line = raw_line.strip()
        if not line or re.fullmatch(r"[|:\-\s]+", line):
            continue
        if len(_norm_ws(line)) >= 8:
            fragments.append(line)
    if len(fragments) < 2:
        return False

    table_fragments = [line for line in fragments if line.lstrip().startswith("|")]
    context_fragments = [
        line for line in fragments if not line.lstrip().startswith("|")
    ]
    if len(table_fragments) < 2:
        return False

    # A header and target row may be non-adjacent, but they must belong to one
    # contiguous markdown table. Searching the whole paper independently lets a
    # model stitch a population header to an unrelated age-table row and turn
    # that false evidence into an enforceable action.
    source_lines = str(source_text or "").splitlines()
    table_blocks: list[tuple[int, list[str]]] = []
    current: list[str] = []
    current_start = 0
    for line_no, source_line in enumerate(source_lines):
        if source_line.lstrip().startswith("|"):
            if not current:
                current_start = line_no
            current.append(source_line)
        elif current:
            table_blocks.append((current_start, current))
            current = []
    if current:
        table_blocks.append((current_start, current))

    for start, block in table_blocks:
        block_text = "\n".join(block)
        if not all(
            _normalized_contains(block_text, fragment) for fragment in table_fragments
        ):
            continue
        if not context_fragments:
            return True

        # Titles/captions belong immediately before a table. Bound the preamble
        # so unrelated prose elsewhere in the paper cannot validate the row.
        preamble_lines: list[str] = []
        cursor = start - 1
        while cursor >= 0 and len(preamble_lines) < 8:
            candidate = source_lines[cursor]
            if candidate.lstrip().startswith("|") or candidate.strip() == "---":
                break
            if candidate.strip():
                preamble_lines.append(candidate)
            cursor -= 1
        preamble_text = "\n".join(reversed(preamble_lines))
        if all(
            _normalized_contains(preamble_text, fragment)
            for fragment in context_fragments
        ):
            return True
    return False


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


def _load_scout_zones(
    pmid: str, gene: str, run_path: Optional[Path]
) -> tuple[Optional[str], Optional[str]]:
    """Data-Scout ``DATA_ZONES.md`` (scout-narrowed high-value regions) for a pmid.

    Looks in ``run_dir/scout_output`` first, then the corpus copy. Returns
    ``(text, path)`` or ``(None, None)``. Never raises."""
    dirs: list[Path] = []
    if run_path is not None:
        dirs.append(run_path / "scout_output")
    if gene:
        dirs.append(REPO_ROOT / "corpus" / gene / pmid)
    for d in dirs:
        f = d / f"{pmid}_DATA_ZONES.md"
        try:
            if f.exists():
                text = _read_text(f)
                if _is_usable_source(text):
                    return text, str(f)
        except OSError:
            continue
    return None, None


def _abstract_from_json(path: Path) -> Optional[str]:
    data = _load_json(_read_text(path))
    if isinstance(data, dict):
        abstract = str(data.get("abstract") or "").strip()
        return abstract or None
    return None


def _load_abstract(
    pmid: str, run_path: Optional[Path], hint_path: Optional[Path]
) -> tuple[Optional[str], Optional[str]]:
    """Abstract text for a pmid: a JSON source-file hint, else run_dir/abstract_json."""
    if hint_path is not None and hint_path.suffix == ".json":
        abstract = _abstract_from_json(hint_path)
        if _is_usable_source(abstract):
            return abstract, str(hint_path)
    if run_path is not None:
        aj = run_path / "abstract_json" / f"{pmid}.json"
        if aj.exists():
            abstract = _abstract_from_json(aj)
            if _is_usable_source(abstract):
                return abstract, str(aj)
    return None, None


def resolve_summary_source(
    pmid: str,
    gene: Optional[str],
    run_dir: Optional[str | Path] = None,
    source_file_hint: Optional[str] = None,
) -> tuple[Optional[str], str, Optional[str]]:
    """Return (text, source_kind, source_path) for the sniff test's source-grounded
    summary: the paper's ABSTRACT plus the Data-Scout ``DATA_ZONES`` (scout-narrowed
    high-value regions where carriers concentrate).

    Full text is DELIBERATELY excluded — the completeness check judges against the
    abstract + scouted data zones, not tens of thousands of characters of prose.
    ``source_kind`` is ``"abstract_scout"`` | ``"scout_only"`` | ``"abstract_only"``
    | ``"none"``. Never raises.
    """
    pmid = str(pmid)
    gene = (gene or "").strip()
    run_path = Path(run_dir) if run_dir else None
    hint_path = _resolve_path(source_file_hint) if source_file_hint else None

    abstract, a_path = _load_abstract(pmid, run_path, hint_path)
    zones, z_path = _load_scout_zones(pmid, gene, run_path)

    parts: list[str] = []
    if abstract:
        parts.append("## Abstract\n\n" + abstract.strip())
    if zones:
        parts.append("## Data-scouted zones\n\n" + zones.strip())
    combined = "\n\n".join(parts).strip()
    if not _is_usable_source(combined):
        return None, "none", None

    if abstract and zones:
        kind = "abstract_scout"
    elif zones:
        kind = "scout_only"
    else:
        kind = "abstract_only"
    return combined, kind, (z_path or a_path)


# Sections that never contain carriers — dropped before the char budget.
_DEAD_SECTION_RE = re.compile(
    r"^#{1,6}\s*(references|bibliography|acknowledg|author\s+contributions?|"
    r"funding|conflicts?\s+of\s+interest|competing\s+interests?|"
    r"declarations\b|data\s+availability)",
    re.IGNORECASE,
)
_HEADER_RE = re.compile(r"^#{1,6}\s+\S")
_TABLEISH_RE = re.compile(
    r"^\s*\||supplemental?\s+file|supplementary\s+table|GVF_FOLDED", re.IGNORECASE
)


# A long variant table is repetitive: the sniff test needs the columns + a sample
# of rows to judge whether it is the right KIND of table (per-variant human
# carriers) for our purposes, not every row. Sample long runs to head + tail.
_PIPE_ROW_RE = re.compile(r"^\s*\|")
TABLE_ROW_SAMPLE_HEAD = 20  # rows kept from the top (header + separator + data)
TABLE_ROW_SAMPLE_TAIL = 3  # rows kept from the bottom
MAX_TABLE_ROWS = (
    TABLE_ROW_SAMPLE_HEAD + TABLE_ROW_SAMPLE_TAIL + 2
)  # only cap clearly-long


def _cap_repetitive_tables(text: str) -> str:
    """Sample any markdown table longer than ``MAX_TABLE_ROWS`` down to its first
    ``TABLE_ROW_SAMPLE_HEAD`` and last ``TABLE_ROW_SAMPLE_TAIL`` rows, with a count
    marker in between. Preserves the table's columns and content shape at a
    fraction of the tokens; short tables pass through untouched."""
    lines = text.split("\n")
    out: list[str] = []
    i, n = 0, len(lines)
    while i < n:
        if not _PIPE_ROW_RE.match(lines[i]):
            out.append(lines[i])
            i += 1
            continue
        j = i
        while j < n and _PIPE_ROW_RE.match(lines[j]):
            j += 1
        run = lines[i:j]
        if len(run) > MAX_TABLE_ROWS:
            omitted = len(run) - TABLE_ROW_SAMPLE_HEAD - TABLE_ROW_SAMPLE_TAIL
            out.extend(run[:TABLE_ROW_SAMPLE_HEAD])
            out.append(
                f"| … {omitted} more rows (same columns) omitted for the sniff test … |"
            )
            out.extend(run[-TABLE_ROW_SAMPLE_TAIL:])
        else:
            out.extend(run)
        i = j
    return "\n".join(out)


def _select_source_for_summary(
    text: str, max_chars: int = DEFAULT_MAX_SOURCE_CHARS
) -> tuple[str, bool]:
    """Trim source to the completeness-relevant text within a char budget.

    Drops dead sections (references/acknowledgments/funding/...); samples long
    repetitive tables to a representative head+tail; then, if still over budget,
    keeps a prose head plus ALL table/supplement-bearing lines (that is where
    carriers concentrate) and truncates the long prose middle. Returns
    (selected_text, truncated) — truncated is True only when the char budget forced
    dropping content, so a gap under truncation is auditable and re-runnable.
    """
    lines = text.split("\n")
    kept: list[str] = []
    skipping = False
    for ln in lines:
        if _HEADER_RE.match(ln):
            skipping = bool(_DEAD_SECTION_RE.match(ln))
        if not skipping:
            kept.append(ln)
    # Sample long repetitive tables first, so a 200-row variant table never
    # dominates the budget or the model's attention.
    body = _cap_repetitive_tables("\n".join(kept))
    if len(body) <= max_chars:
        return body, False

    body_lines = body.split("\n")
    table_lines = [ln for ln in body_lines if _TABLEISH_RE.search(ln)]
    head_budget = int(max_chars * 0.6)
    head = body[:head_budget]
    nl = head.rfind("\n")  # cut at a line boundary so a variant notation isn't split
    if nl > 0:
        head = head[:nl]
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
    """Assemble one payload per paper, including zero-count-row papers.

    ``facts`` retains provenance-rich detail for at most ``max_facts`` rows.
    ``captured_fact_index`` is a complete compact identity/count inventory and
    is never capped; source-grounded matching must use it so large papers cannot
    turn already-captured rows into false misses.
    """
    # Local cursor with Row factory — don't mutate the caller's connection.
    cur = conn.cursor()
    cur.row_factory = sqlite3.Row
    sel_trust = (
        ", pd.trust_tier AS trust_tier, pd.trust_reasons AS trust_reasons"
        if _penetrance_has_trust(conn)
        else ""
    )
    where = ""
    paper_where = ""
    params: list[Any] = []
    if pmids:
        placeholders = ",".join("?" for _ in pmids)
        where = f" WHERE pd.pmid IN ({placeholders})"
        paper_where = f" WHERE p.pmid IN ({placeholders})"
        params = list(pmids)

    # Seed from papers, not penetrance_data. A source-grounded audit is most
    # valuable when extraction produced zero count rows.
    grouped: dict[str, dict[str, Any]] = {}
    for row in cur.execute(
        f"""
        SELECT p.pmid AS pmid, p.title AS title, p.gene_symbol AS gene_symbol
        FROM papers p
        {paper_where}
        ORDER BY p.pmid
        """,
        params,
    ):
        pmid = str(row["pmid"])
        grouped[pmid] = {
            "pmid": pmid,
            "gene": row["gene_symbol"],
            "title": row["title"],
            "facts": [],
            "captured_fact_index": [],
            "n_facts_total": 0,
        }

    variant_cols = {row[1] for row in conn.execute("PRAGMA table_info(variants)")}
    provenance_cols = {
        row[1] for row in conn.execute("PRAGMA table_info(fact_provenance)")
    }
    evidence_by_variant: dict[tuple[int, str], list[dict[str, Any]]] = {}
    if provenance_cols:
        provenance_where = ""
        provenance_params: list[Any] = []
        if pmids:
            placeholders = ",".join("?" for _ in pmids)
            provenance_where = f" AND pmid IN ({placeholders})"
            provenance_params = list(pmids)
        for ev in cur.execute(
            f"""
            SELECT variant_id, pmid, fact_type, fact_value, source_table,
                   source_row, source_column, evidence_quote, count_type,
                   source_layer
            FROM fact_provenance
            WHERE fact_type IN (
                'patient_count', 'total_carriers_observed', 'affected_count',
                'unaffected_count', 'uncertain_count'
            )
            {provenance_where}
            ORDER BY variant_id, pmid, provenance_id
            """,
            provenance_params,
        ):
            evidence = {
                key: ev[key]
                for key in (
                    "fact_type",
                    "fact_value",
                    "source_table",
                    "source_row",
                    "source_column",
                    "evidence_quote",
                    "count_type",
                    "source_layer",
                )
                if ev[key] not in (None, "")
            }
            key = (int(ev["variant_id"]), str(ev["pmid"]))
            if evidence and evidence not in evidence_by_variant.setdefault(key, []):
                evidence_by_variant[key].append(evidence)
    structural_expr = (
        "v.structural_description"
        if "structural_description" in variant_cols
        else "NULL"
    )
    variant_class_expr = (
        "v.variant_class" if "variant_class" in variant_cols else "NULL"
    )
    variant_identity_expr = (
        "COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position, "
        f"{structural_expr})"
    )
    sql = f"""
        SELECT pd.pmid AS pmid,
               pd.penetrance_id         AS penetrance_id,
               v.variant_id              AS variant_id,
               {variant_identity_expr} AS variant,
               v.protein_notation      AS protein_notation,
               v.cdna_notation         AS cdna_notation,
               v.genomic_position      AS genomic_position,
               {structural_expr}       AS structural_description,
               {variant_class_expr}    AS variant_class,
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
        ORDER BY pd.pmid, pd.penetrance_id
    """
    for row in cur.execute(sql, params):
        pmid = str(row["pmid"])
        bucket = grouped.setdefault(
            pmid,
            {
                "pmid": pmid,
                "gene": row["gene_symbol"],
                "title": row["title"],
                "facts": [],
                "captured_fact_index": [],
                "n_facts_total": 0,
            },
        )
        counts = {
            "total_carriers": _as_int(row["total_carriers"]),
            "affected": _as_int(row["affected"]),
            "unaffected": _as_int(row["unaffected"]),
            "uncertain": _as_int(row["uncertain"]),
        }
        compact_counts = {k: v for k, v in counts.items() if v is not None}
        identifiers = {
            key: row[key]
            for key in (
                "protein_notation",
                "cdna_notation",
                "genomic_position",
                "structural_description",
                "variant_class",
            )
            if row[key] not in (None, "")
        }
        fact_number = bucket["n_facts_total"] + 1
        compact_fact: dict[str, Any] = {
            "fact_id": int(row["penetrance_id"]),
            "variant_id": int(row["variant_id"]),
            "variant": row["variant"] or "(unnamed)",
            "counts": compact_counts,
        }
        if identifiers:
            compact_fact["identifiers"] = identifiers
        bucket["captured_fact_index"].append(compact_fact)
        bucket["n_facts_total"] = fact_number

        fact: dict[str, Any] = {
            "fact_id": int(row["penetrance_id"]),
            "variant_id": int(row["variant_id"]),
            "variant": row["variant"] or "(unnamed)",
            "counts": compact_counts,
        }
        if identifiers:
            fact["identifiers"] = identifiers
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
        raw_quotes = _load_json(row["key_quotes"])
        if isinstance(raw_quotes, list):
            quotes = _clip("\n".join(str(value) for value in raw_quotes if value))
        else:
            quotes = _clip(row["key_quotes"])
        if quotes:
            fact["quote"] = quotes
        evidence_context = []
        count_values = {
            "patient_count": counts["total_carriers"],
            "total_carriers_observed": counts["total_carriers"],
            "affected_count": counts["affected"],
            "unaffected_count": counts["unaffected"],
            "uncertain_count": counts["uncertain"],
        }
        for evidence in evidence_by_variant.get((int(row["variant_id"]), pmid), []):
            fact_type = str(evidence.get("fact_type") or "")
            expected_value = count_values.get(fact_type)
            if (
                expected_value is None
                or _as_int(evidence.get("fact_value")) != expected_value
            ):
                continue
            clipped = {
                key: _clip(value, MAX_QUOTE_CHARS if key == "evidence_quote" else 160)
                for key, value in evidence.items()
            }
            clipped = {key: value for key, value in clipped.items() if value}
            if clipped and clipped not in evidence_context:
                evidence_context.append(clipped)
            if len(evidence_context) >= 4:
                break
        if evidence_context:
            fact["table_evidence"] = evidence_context
        if sel_trust:
            all_reasons = _load_json(row["trust_reasons"])
            all_reasons = all_reasons if isinstance(all_reasons, list) else []
            structural_reasons = [
                str(reason)
                for reason in all_reasons
                if not str(reason).startswith("paper_final_check:")
            ]
            # The reviewer must see/hash the structural extraction state, not
            # its own prior decision. This keeps a durable PFC composition from
            # invalidating the next reviewer payload or creating a feedback loop.
            legacy_quarantine = (
                str(row["trust_tier"] or "").lower() == "quarantine" and not all_reasons
            )
            if structural_reasons:
                fact["trust_tier"] = "quarantine"
                fact["trust_reasons"] = structural_reasons
            elif legacy_quarantine:
                fact["trust_tier"] = "quarantine"
            else:
                fact["trust_tier"] = "trusted"
        bucket["facts"].append(fact)

    payloads = []
    for payload in grouped.values():
        if len(payload["facts"]) > max_facts:

            def review_priority(fact: dict[str, Any]) -> tuple[int, int, int, int]:
                values = [
                    abs(value)
                    for value in (fact.get("counts") or {}).values()
                    if isinstance(value, int)
                ]
                max_count = max(values, default=0)
                label = str(fact.get("count_label") or "").lower()
                ambiguous_label = int(
                    any(
                        token in label
                        for token in ("age", "family", "case", "cohort", "screen")
                    )
                )
                return (
                    ambiguous_label,
                    int(bool(fact.get("table_evidence"))),
                    max_count,
                    -int(fact.get("fact_id") or 0),
                )

            selected = sorted(payload["facts"], key=review_priority, reverse=True)[
                :max_facts
            ]
            payload["facts"] = sorted(
                selected, key=lambda fact: fact.get("fact_id") or 0
            )
        payload["truncated"] = payload["n_facts_total"] > len(payload["facts"])
        payloads.append(payload)
    return payloads


def payload_content_hash(payload: dict[str, Any]) -> str:
    """Hash the semantic DB payload independently of row/insertion order.

    A changed count, identity, fact/variant ID, provenance field, title, or trust
    decision must invalidate version-skip even when the source text and reviewer
    model did not change.
    """

    def stable_items(items: Any) -> list[Any]:
        values = list(items or [])
        return sorted(
            values,
            key=lambda item: json.dumps(
                item, ensure_ascii=False, sort_keys=True, separators=(",", ":")
            ),
        )

    canonical = {
        "pmid": payload.get("pmid"),
        "gene": payload.get("gene"),
        "title": payload.get("title"),
        "facts": stable_items(payload.get("facts")),
        "captured_fact_index": stable_items(payload.get("captured_fact_index")),
        "n_facts_total": payload.get("n_facts_total", 0),
        "truncated": bool(payload.get("truncated")),
    }
    encoded = json.dumps(
        canonical, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


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


def _as_bool(value: Any) -> bool:
    """Robust truthiness for LLM JSON: the string "false"/"no"/"0" must be False
    (Python's bool("false") is True, which would silently hide a real miss)."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes", "y", "t")
    return False


def _norm_flags(
    raw: dict[str, Any], payload: Optional[dict[str, Any]] = None
) -> tuple[list[dict[str, Any]], str]:
    """Normalize findings and seal their fact-to-variant bindings.

    The reviewer chooses exact ``fact_ids`` from the payload.  ``variant_ids``
    are copied from that same payload rather than trusted from model output, so
    the enforcement step can reject a fact-id collision after a re-migration.
    """
    fact_variants: dict[int, int] = {}
    for fact in (payload or {}).get("captured_fact_index") or []:
        if not isinstance(fact, dict):
            continue
        fact_id = _as_int(fact.get("fact_id"))
        variant_id = _as_int(fact.get("variant_id"))
        if fact_id is not None and variant_id is not None:
            fact_variants[fact_id] = variant_id

    flags: list[dict[str, Any]] = []
    raw_flags = raw.get("flags")
    if isinstance(raw_flags, dict):  # model sometimes returns a single object
        raw_flags = [raw_flags]
    for item in raw_flags or []:
        if not isinstance(item, dict):
            continue
        severity = str(item.get("severity") or "").strip().lower()
        raw_fact_ids = item.get("fact_ids")
        if raw_fact_ids is None and item.get("fact_id") is not None:
            raw_fact_ids = [item.get("fact_id")]
        if not isinstance(raw_fact_ids, list):
            raw_fact_ids = []
        fact_ids: list[int] = []
        for value in raw_fact_ids:
            fact_id = _as_int(value)
            if fact_id is not None and fact_id > 0 and fact_id not in fact_ids:
                fact_ids.append(fact_id)

        raw_fields = item.get("fields")
        if isinstance(raw_fields, str):
            raw_fields = [raw_fields]
        fields = []
        for value in raw_fields or []:
            field = str(value or "").strip().lower()
            field = FLAG_FIELD_ALIASES.get(field, field)
            if field in VALID_FLAG_FIELDS and field not in fields:
                fields.append(field)

        reason_code = str(item.get("reason_code") or "other").strip().lower()
        if reason_code not in VALID_FLAG_REASONS:
            reason_code = "other"
        flags.append(
            {
                "fact_ids": fact_ids,
                "variant_ids": [
                    fact_variants[fact_id]
                    for fact_id in fact_ids
                    if fact_id in fact_variants
                ],
                "variant": str(item.get("variant") or "").strip() or "(unspecified)",
                "fields": fields,
                "reason_code": reason_code,
                "evidence_quote": str(item.get("evidence_quote") or "").strip()[
                    :MAX_QUOTE_CHARS
                ],
                "issue": str(item.get("issue") or "").strip(),
                "severity": severity if severity in VALID_SEVERITIES else "medium",
            }
        )
    verdict = str(raw.get("paper_verdict") or "").strip().lower()
    if verdict not in VALID_VERDICTS:
        verdict = "flag" if flags else "ok"
    return flags, verdict


class FinalCheckResponseError(ValueError):
    """The model returned JSON that cannot support an honest final verdict."""


def _validate_common_response(raw: Any, *, mode: str) -> dict[str, Any]:
    """Validate fields shared by both final-check response schemas.

    Empty/partial JSON is an outage signal, not an implicit ``ok``. Raising here
    lets the orchestration persist ``verdict=error``; error rows are deliberately
    excluded from version-skip and therefore remain retryable.
    """
    if not isinstance(raw, dict) or not raw:
        raise FinalCheckResponseError(f"{mode} response must be a non-empty object")

    required = {"paper_verdict", "confidence", "flags"}
    missing = sorted(required - raw.keys())
    if missing:
        raise FinalCheckResponseError(
            f"{mode} response missing required fields: {', '.join(missing)}"
        )

    verdict = str(raw.get("paper_verdict") or "").strip().lower()
    if verdict not in VALID_VERDICTS:
        raise FinalCheckResponseError(
            f"{mode} response has invalid paper_verdict: {verdict or '<empty>'}"
        )

    confidence = raw.get("confidence")
    if isinstance(confidence, bool):
        raise FinalCheckResponseError(f"{mode} confidence must be a number in [0, 1]")
    try:
        confidence_value = float(confidence)
    except (TypeError, ValueError) as exc:
        raise FinalCheckResponseError(
            f"{mode} confidence must be a number in [0, 1]"
        ) from exc
    if not 0.0 <= confidence_value <= 1.0:
        raise FinalCheckResponseError(f"{mode} confidence must be a number in [0, 1]")

    if not isinstance(raw.get("flags"), list):
        raise FinalCheckResponseError(f"{mode} flags must be a list")
    return raw


def _validate_db_response(raw: Any) -> dict[str, Any]:
    raw = _validate_common_response(raw, mode="DB-only final-check")
    if "summary" not in raw or not isinstance(raw.get("summary"), str):
        raise FinalCheckResponseError(
            "DB-only final-check response missing required string field: summary"
        )
    if not raw["summary"].strip():
        raise FinalCheckResponseError("DB-only final-check summary must not be empty")
    return raw


def _validate_summary_response(raw: Any) -> dict[str, Any]:
    raw = _validate_common_response(raw, mode="source-grounded final-check")
    required = {
        "study_summary",
        "cohort_sources",
        "carrier_groups",
        "completeness",
    }
    missing = sorted(required - raw.keys())
    if missing:
        raise FinalCheckResponseError(
            "source-grounded final-check response missing required fields: "
            + ", ".join(missing)
        )
    if (
        not isinstance(raw.get("study_summary"), str)
        or not raw["study_summary"].strip()
    ):
        raise FinalCheckResponseError(
            "source-grounded final-check study_summary must be a non-empty string"
        )
    if not isinstance(raw.get("cohort_sources"), list):
        raise FinalCheckResponseError(
            "source-grounded final-check cohort_sources must be a list"
        )
    groups = raw.get("carrier_groups")
    if not isinstance(groups, list):
        raise FinalCheckResponseError(
            "source-grounded final-check carrier_groups must be a list"
        )
    for index, group in enumerate(groups):
        if not isinstance(group, dict):
            raise FinalCheckResponseError(
                f"source-grounded carrier_groups[{index}] must be an object"
            )
        missing_group = {
            "variant",
            "evidence_quote",
            "in_extraction",
        } - group.keys()
        if missing_group:
            raise FinalCheckResponseError(
                f"source-grounded carrier_groups[{index}] missing required fields: "
                + ", ".join(sorted(missing_group))
            )
    completeness = raw.get("completeness")
    if not isinstance(completeness, dict):
        raise FinalCheckResponseError(
            "source-grounded final-check completeness must be an object"
        )
    status = str(completeness.get("status") or "").strip().lower()
    if status not in VALID_COMPLETENESS:
        raise FinalCheckResponseError(
            "source-grounded final-check completeness.status must be one of "
            + ", ".join(VALID_COMPLETENESS)
        )
    return raw


def _validate_checker_result(result: Any, *, source_grounded: bool) -> dict[str, Any]:
    """Defense-in-depth for injected/custom checkers after normalization."""
    mode = "source-grounded" if source_grounded else "DB-only"
    if not isinstance(result, dict) or not result:
        raise FinalCheckResponseError(f"{mode} checker returned an empty result")
    if result.get("verdict") not in VALID_VERDICTS:
        raise FinalCheckResponseError(
            f"{mode} checker returned invalid verdict: {result.get('verdict')!r}"
        )
    required = {"confidence", "flags", "summary"}
    missing = sorted(required - result.keys())
    if missing:
        raise FinalCheckResponseError(
            f"{mode} checker result missing required fields: {', '.join(missing)}"
        )
    if source_grounded:
        if not isinstance(result.get("carrier_groups"), list) or not isinstance(
            result.get("completeness"), dict
        ):
            raise FinalCheckResponseError(
                "source-grounded checker result requires carrier_groups and completeness"
            )
    return result


def normalize_result(raw: Any, payload: dict[str, Any]) -> dict[str, Any]:
    """Coerce the DB-only sniff-test JSON into a stable record."""
    raw = _validate_db_response(raw)
    flags, verdict = _norm_flags(raw, payload)
    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "source_grounded": False,
        "verdict": verdict,
        "confidence": _coerce_confidence(raw.get("confidence")),
        "n_facts": int(
            payload.get("n_facts_total", len(payload.get("facts") or [])) or 0
        ),
        "n_flagged": len(flags),
        "flags": flags,
        "summary": str(raw.get("summary") or "").strip(),
    }


def normalize_summary_result(
    raw: Any, payload: dict[str, Any], source_text: str, truncated: bool = False
) -> dict[str, Any]:
    """Coerce the source-grounded JSON into a stable record + derive the
    missed-carrier signal (never asked directly; derived from in_extraction).
    ``truncated`` = the fed source was truncated, so "complete" is not claimable.
    """
    raw = _validate_summary_response(raw)
    flags, _model_verdict = _norm_flags(
        raw, payload
    )  # verdict re-derived from grounded signal
    for flag in flags:
        quote = flag.get("evidence_quote") or ""
        flag["evidence_quote_verified"] = bool(
            quote and _source_evidence_verified(source_text, quote)
        )

    raw_groups = raw.get("carrier_groups")
    if isinstance(raw_groups, dict):  # tolerate a single object
        raw_groups = [raw_groups]
    groups: list[dict[str, Any]] = []
    for g in raw_groups or []:
        if not isinstance(g, dict):
            continue
        in_extraction = _as_bool(g.get("in_extraction"))
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
                "quote_verified": 1
                if _source_evidence_verified(source_text, quote)
                else 0,
                "status": "reported_extracted" if in_extraction else "reported_missing",
            }
        )

    # Count across ALL groups the model returned (before any persistence cap).
    # A reported-missing group is only a GROUNDED miss when its evidence quote is
    # actually present in the (scouted) source; ungrounded "missing" claims are the
    # model's over-report and must not drive the gaps/verdict signal. This is the
    # over-sensitivity fix: n_missing counts quote-verified misses only.
    n_reported = len(groups)
    n_in = sum(1 for g in groups if g["in_extraction"])
    missing_groups = [g for g in groups if g["status"] == "reported_missing"]
    n_missing = sum(1 for g in missing_groups if g["quote_verified"])
    n_missing_unverified = len(missing_groups) - n_missing

    # Persist grounded misses first, then ungrounded misses, then already-extracted
    # rows, so the cap never drops a quote-verified miss.
    def _group_priority(g: dict[str, Any]) -> int:
        if g["status"] != "reported_missing":
            return 2
        return 0 if g["quote_verified"] else 1

    groups.sort(key=_group_priority)
    capped = n_reported > MAX_CARRIER_GROUPS
    kept_groups = groups[:MAX_CARRIER_GROUPS]

    # Completeness is DERIVED from the grounded miss signal, never the model's
    # (over-eager) self-report: "gaps" REQUIRES a quote-verified miss; ungrounded
    # miss claims are only "unsure"; a model "complete" is trusted just when no
    # misses were claimed at all (and the source was not truncated).
    comp = raw.get("completeness")
    comp = comp if isinstance(comp, dict) else {}
    if n_missing > 0:
        status = "gaps"
    elif n_missing_unverified > 0:
        status = "unsure"
    else:
        model_status = str(comp.get("status") or "").strip().lower()
        complete_ok = (n_reported > 0 or model_status == "complete") and not truncated
        status = "complete" if complete_ok else "unsure"

    notes = str(comp.get("notes") or "").strip()
    if n_missing_unverified:
        notes = (
            f"{notes} [{n_missing_unverified} ungrounded missing-carrier "
            "claim(s) not counted]"
        ).strip()
    if capped:
        notes = (
            f"{notes} [carrier_groups capped at {MAX_CARRIER_GROUPS} of {n_reported}]"
        ).strip()

    # Verdict = TRUST only, orthogonal to completeness. Flag ONLY on a real
    # per-variant trust flag (wrong-gene, count-role, denominator, ...). A paper
    # that merely LOOKS incomplete is NOT flagged — the missed carriers live in the
    # separate ``completeness`` signal (status=gaps + carrier_groups worklist), so
    # "flag" stays a precise, actionable "a captured count is suspect" queue rather
    # than firing on every paper the model thinks it under-read. Also reconciles a
    # model 'ok' that shipped real flags up to 'flag'.
    verdict = "flag" if flags else "ok"

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
        "n_facts": int(
            payload.get("n_facts_total", len(payload.get("facts") or [])) or 0
        ),
        "n_flagged": len(flags),
        "flags": flags,
        "summary": summary,
        "cohort_sources": cohort_sources,
        "carrier_groups": kept_groups,
        "completeness": {
            "status": status,
            "n_reported": n_reported,
            "n_in_extraction": n_in,
            "n_missing": n_missing,
            "n_missing_unverified": n_missing_unverified,
            "notes": notes,
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
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS paper_final_check_attempts (
            attempt_id      INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid            TEXT,
            check_version   TEXT,
            verdict         TEXT,
            summary         TEXT,
            attempted_at    TEXT
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pfcat_pmid ON paper_final_check_attempts(pmid)"
    )


def _record_attempt(
    conn: sqlite3.Connection,
    *,
    result: dict[str, Any],
    version: str,
    attempted_at: str,
) -> None:
    """Persist every reviewer attempt, including transient errors.

    ``paper_final_check`` remains the last successful/current decision row;
    this append-only table keeps an outage auditable without erasing durable
    grounded findings that still match the current payload.
    """
    conn.execute(
        """
        INSERT INTO paper_final_check_attempts (
            pmid, check_version, verdict, summary, attempted_at
        ) VALUES (?,?,?,?,?)
        """,
        (
            result.get("pmid"),
            version,
            result.get("verdict"),
            result.get("summary"),
            attempted_at,
        ),
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


def _existing_check(
    conn: sqlite3.Connection, pmid: str
) -> tuple[Optional[str], Optional[str]]:
    """(check_version, verdict) of a paper's recorded row, or (None, None)."""
    row = conn.execute(
        "SELECT check_version, verdict FROM paper_final_check WHERE pmid = ?", (pmid,)
    ).fetchone()
    return (row[0], row[1]) if row else (None, None)


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

    def check(
        self, payload: dict[str, Any], source_text: str, truncated: bool = False
    ) -> dict[str, Any]:
        raw = self._caller.call_llm_json(
            build_paper_summary_prompt(payload, source_text)
        )
        return normalize_summary_result(raw, payload, source_text, truncated)


def _error_result(payload: dict[str, Any], exc: Exception) -> dict[str, Any]:
    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "source_grounded": False,
        "verdict": "error",
        "confidence": None,
        "n_facts": int(
            payload.get("n_facts_total", len(payload.get("facts") or [])) or 0
        ),
        "n_flagged": 0,
        "flags": [],
        "summary": f"final check error: {exc}",
    }


def _skipped_result(payload: dict[str, Any], reason: str) -> dict[str, Any]:
    """Durable marker for a paper that cannot be reviewed without source/facts."""
    return {
        "pmid": payload.get("pmid"),
        "gene_symbol": payload.get("gene"),
        "source_grounded": False,
        "verdict": "skipped",
        "confidence": None,
        "n_facts": 0,
        "n_flagged": 0,
        "flags": [],
        "summary": reason,
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
        model = model or settings.get_paper_final_check_model()
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
        "skipped_empty_no_source": 0,
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
        # Review and enforcement must hash the same payload shape even on a DB
        # created before field-level trust columns existed. The import is local
        # to avoid a module-load cycle with the deterministic composer.
        from pipeline.paper_final_check_gate import ensure_gate_schema

        ensure_gate_schema(conn)
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
                source_text, source_kind, source_path = resolve_summary_source(
                    pmid,
                    paper_gene,
                    run_dir=run_dir,
                    source_file_hint=(hints.get(pmid) or {}).get("source_file"),
                )
            use_summary = bool(
                source_grounded
                and source_text
                and source_kind in {"abstract_scout", "scout_only", "abstract_only"}
            )

            selected, truncated = "", False
            if use_summary:
                selected, truncated = _select_source_for_summary(
                    source_text, max_source_chars
                )
                # An abstract can ground reported groups, but it can never prove
                # full-paper completeness. Reuse the truncation guard so a model
                # cannot persist "complete" from an abstract-only view.
                truncated = truncated or source_kind == "abstract_only"
                # Fold a signature of the SELECTED source into the version so a
                # later-enriched source (folded supplements) re-runs instead of
                # being version-skipped.
                source_sig = hashlib.sha256(selected.encode("utf-8")).hexdigest()[:12]
                payload_sig = payload_content_hash(payload)[:12]
                intended_version = f"{sum_version}-{source_sig}-{payload_sig}"
            else:
                payload_sig = payload_content_hash(payload)[:12]
                intended_version = f"{db_version}-{payload_sig}"

            empty_no_source = not payload.get("captured_fact_index") and not use_summary

            # Version-skip: same reviewer generation AND same source already
            # recorded, same DB payload, and NOT a prior error (transient
            # failures are retryable).
            prev_version, prev_verdict = _existing_check(conn, pmid)
            if prev_version == intended_version and prev_verdict != "error":
                stats["skipped"] += 1
                if empty_no_source:
                    stats["skipped_empty_no_source"] += 1
                continue

            if empty_no_source:
                reason = "paper has no extracted count facts and no usable source text"
                logger.info("paper final check skipped pmid=%s: %s", pmid, reason)
                result = _skipped_result(payload, reason)
                _record_row(
                    conn,
                    result=result,
                    model=model,
                    reasoning_effort=reasoning_effort,
                    prompt_version=PROMPT_VERSION,
                    version=intended_version,
                    checked_at=checked_at,
                    source_kind=source_kind,
                    source_file=source_path,
                )
                _record_carrier_groups(
                    conn,
                    pmid=pmid,
                    gene_symbol=result.get("gene_symbol"),
                    groups=[],
                    version=intended_version,
                    checked_at=checked_at,
                )
                stats["skipped"] += 1
                stats["skipped_empty_no_source"] += 1
                if source_kind == "none":
                    stats["source_absent"] += 1
                conn.commit()
                continue

            if use_summary:
                record_prompt_version = SUMMARY_PROMPT_VERSION
                record_source_chars: Optional[int] = len(selected)
                record_source_truncated: Optional[int] = 1 if truncated else 0
                if summary_checker is None:
                    summary_checker = PaperSummaryChecker(
                        model=model,
                        reasoning_effort=reasoning_effort,
                        max_tokens=SUMMARY_MAX_TOKENS,
                    )
                try:
                    result = summary_checker.check(payload, selected, truncated)
                    result = _validate_checker_result(result, source_grounded=True)
                    stats["source_grounded"] += 1
                    comp = result.get("completeness") or {}
                    stats["missing_carriers"] += int(comp.get("n_missing") or 0)
                except Exception as exc:  # noqa: BLE001 - one paper must not abort
                    logger.warning("paper summary failed for pmid=%s: %s", pmid, exc)
                    result = _error_result(payload, exc)
            else:
                record_prompt_version = PROMPT_VERSION
                record_source_chars = None
                record_source_truncated = None
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
                    result = _validate_checker_result(result, source_grounded=False)
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "paper final check failed for pmid=%s: %s", pmid, exc
                    )
                    result = _error_result(payload, exc)

            _record_attempt(
                conn,
                result=result,
                version=intended_version,
                attempted_at=checked_at,
            )
            preserve_previous_success = result.get(
                "verdict"
            ) == "error" and prev_verdict in {"ok", "flag"}
            if preserve_previous_success:
                logger.warning(
                    "paper final check retained prior successful finding for "
                    "pmid=%s after transient reviewer error",
                    pmid,
                )
            else:
                _record_row(
                    conn,
                    result=result,
                    model=model,
                    reasoning_effort=reasoning_effort,
                    prompt_version=record_prompt_version,
                    version=intended_version,
                    checked_at=checked_at,
                    source_kind=source_kind,
                    source_file=source_path,
                    source_chars=record_source_chars,
                    source_truncated=record_source_truncated,
                )
                # A successful re-check replaces the paper's carrier groups. An
                # error must not erase prior quote-verified misses; the failed
                # attempt is already durable in paper_final_check_attempts.
                _record_carrier_groups(
                    conn,
                    pmid=pmid,
                    gene_symbol=result.get("gene_symbol"),
                    groups=(result.get("carrier_groups") or [])
                    if result.get("source_grounded")
                    else [],
                    version=intended_version,
                    checked_at=checked_at,
                )

            verdict = result.get("verdict") or "ok"
            stats[verdict] = stats.get(verdict, 0) + 1
            stats["flagged_facts"] += int(result.get("n_flagged") or 0)
            stats["checked"] += 1
            # Incremental commit: a crash mid-run keeps completed papers.
            conn.commit()

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
