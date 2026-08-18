"""Table-routing extraction.

Two-phase replacement for the "load whole paper into the LLM" extraction path:

  1. Enumerate every markdown table in the paper, capture its caption + header +
     a small data-row preview.
  2. Send only that compact inventory to the Tier-2 LLM (Kimi). The router
     returns which tables hold variant data and a column→field mapping for
     each one (e.g. "col 2 = cdna, col 3 = protein, col 4 = patient_count").
  3. Run a deterministic parser over the routed tables — no LLM cell reads.

For papers with no usable variant tables (case reports, narrative-only studies)
the caller falls back to the existing full-text extraction path, so we never
*lose* recall, we just stop spending tokens on the easy cases.

Design notes:
  - The router prompt is intentionally small (~1-2k tokens in, <500 tokens out)
    so 8-way concurrency at Azure barely registers against TPM ceilings.
  - Column mapping comes back as integer indices into the header row's
    pipe-split cells, so the parser doesn't have to do any header-string
    fuzzy-matching at extract time.
  - The parser preserves the existing variant schema produced by
    `ExpertExtractor._parse_markdown_table_variants`, so downstream
    aggregation / SQLite migration code is unchanged.
"""

from __future__ import annotations

import json
import logging
import re
from functools import lru_cache
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from utils.gene_metadata import gene_alias_regex, get_gene_aliases, known_gene_aliases

logger = logging.getLogger(__name__)


# Headers we accept in the column mapping — also used to validate router output.
_KNOWN_FIELDS = {
    "gene",
    "cdna",
    "protein",
    "patient_count",
    "affected",
    "unaffected",
    "uncertain",
    "phenotype",
    "clinical_significance",
}

_INFER_ROW_PATIENT_COUNT = -1


@dataclass
class MarkdownTable:
    """One markdown table sliced out of a paper's full-text markdown."""

    table_id: str
    caption: Optional[str]
    header_line: str
    header_cells: List[str]
    data_lines: List[str]
    char_start: int
    char_end: int
    section: Optional[str] = None

    def header_preview(self, max_chars: int = 240) -> str:
        # Render parsed cells with `|` separators and per-cell trimming so the
        # router sees useful content even when the original markdown was padded
        # with hundreds of spaces per cell (PMC-converted supplementary tables).
        cells = [c.strip()[:60] for c in self.header_cells if c.strip()]
        h = "| " + " | ".join(cells) + " |" if cells else self.header_line.strip()
        return h if len(h) <= max_chars else h[: max_chars - 1] + "…"

    def data_preview(self, n_rows: int = 2, max_chars: int = 240) -> List[str]:
        out: List[str] = []
        for row in self.data_lines[:n_rows]:
            cells = [c.strip()[:60] for c in _split_pipe_row(row)]
            r = "| " + " | ".join(cells) + " |" if cells else row.strip()
            if len(r) > max_chars:
                r = r[: max_chars - 1] + "…"
            out.append(r)
        return out


@dataclass
class RoutedTable:
    """A table the router flagged as containing variant data, plus its column mapping."""

    table_id: str
    column_mapping: Dict[str, int]
    confidence: Optional[float] = None
    notes: Optional[str] = None
    table: Optional[MarkdownTable] = None


@dataclass
class RouterResult:
    """Outcome of the routing phase."""

    routed_tables: List[RoutedTable] = field(default_factory=list)
    raw_response: Optional[str] = None
    error: Optional[str] = None
    used_fallback: bool = False
    # Every table id `enumerate_markdown_tables` produced for this paper, routed
    # or not. Recorded so a decline can be told apart from "there was nothing to
    # route" when auditing the router after the fact.
    enumerated_table_ids: List[str] = field(default_factory=list)


# Outcome codes written to extraction_metadata["table_router_outcome"].
ROUTER_OUTCOME_APPROVED = "approved"
# Router approved tables, but the caller found too few variants to short-circuit
# the full-text path, so the router output was demoted to extraction hints.
ROUTER_OUTCOME_APPROVED_BELOW_THRESHOLD = "approved_below_threshold"
ROUTER_OUTCOME_NO_USABLE_TABLES = "no_usable_tables"
ROUTER_OUTCOME_LLM_ERROR = "llm_error"
ROUTER_OUTCOME_CRASHED = "crashed"
ROUTER_OUTCOME_DISABLED = "disabled"
# Feature flag on, but a precondition (gene symbol, source text) was missing.
ROUTER_OUTCOME_SKIPPED = "skipped"

ROUTER_OUTCOMES = (
    ROUTER_OUTCOME_APPROVED,
    ROUTER_OUTCOME_APPROVED_BELOW_THRESHOLD,
    ROUTER_OUTCOME_NO_USABLE_TABLES,
    ROUTER_OUTCOME_LLM_ERROR,
    ROUTER_OUTCOME_CRASHED,
    ROUTER_OUTCOME_DISABLED,
    ROUTER_OUTCOME_SKIPPED,
)

# Caps so papers with dozens of supplementary tables cannot bloat the
# extraction JSON. Counts stay exact; only the per-table detail is capped.
MAX_ROUTER_DETAIL_TABLES = 12
MAX_ROUTER_NOTE_CHARS = 200
MAX_ROUTER_ERROR_CHARS = 300


def build_router_decision_metadata(
    *,
    model: Optional[str],
    outcome: str,
    routed: Optional[Sequence[RoutedTable]] = None,
    enumerated_table_ids: Optional[Sequence[str]] = None,
    variants_parsed: int = 0,
    error: Optional[str] = None,
    max_detail_tables: int = MAX_ROUTER_DETAIL_TABLES,
) -> Dict[str, Any]:
    """Build a compact, JSON-safe record of one router decision.

    The router's column->field mapping is applied deterministically to every row
    of a routed table, so a single mapping error scales to hundreds of rows.
    Persisting the mapping (plus what was enumerated but declined) is what makes
    that judgment call auditable and scorable after the run.

    Observability only: nothing here feeds back into routing.
    """
    routed_tables = list(routed or [])
    enumerated = [str(t) for t in (enumerated_table_ids or [])]
    routed_ids = {r.table_id for r in routed_tables}
    declined = [t for t in enumerated if t not in routed_ids]

    mappings: Dict[str, Dict[str, int]] = {}
    confidences: Dict[str, float] = {}
    notes: Dict[str, str] = {}
    for r in routed_tables[:max_detail_tables]:
        mappings[r.table_id] = dict(r.column_mapping or {})
        if r.confidence is not None:
            confidences[r.table_id] = r.confidence
        if r.notes:
            notes[r.table_id] = str(r.notes)[:MAX_ROUTER_NOTE_CHARS]

    meta: Dict[str, Any] = {
        "table_router_model": model,
        "table_router_outcome": outcome,
        "table_router_tables_enumerated": len(enumerated),
        "table_router_tables_routed": len(routed_tables),
        "table_router_variants_parsed": int(variants_parsed),
        "table_router_mappings": mappings,
        "table_router_confidences": confidences,
        "table_router_notes": notes,
    }
    if declined:
        meta["table_router_tables_declined"] = declined[:max_detail_tables]
    if len(routed_tables) > max_detail_tables or len(declined) > max_detail_tables:
        meta["table_router_detail_truncated"] = True
    if error:
        meta["table_router_error"] = str(error)[:MAX_ROUTER_ERROR_CHARS]
    return meta


_TABLE_CAPTION_RE = re.compile(
    r"^\s*(e?table\s+\d+[a-z]?[.:]|tbl\.?\s*\d+[.:])", re.IGNORECASE
)

# A markdown heading that names a table ("### Table 1", "## eTable 2A"). PMC's
# XML→markdown conversion emits the label as a heading and the descriptive
# caption as a separate paragraph beneath it.
_TABLE_HEADING_RE = re.compile(
    r"^\s*#{1,6}\s*(e?table\s+\d+[a-z]?|tbl\.?\s*\d+[a-z]?)\s*[.:]?\s*$",
    re.IGNORECASE,
)


def _split_pipe_row(line: str) -> List[str]:
    """Split a markdown row on `|`, dropping empty leading/trailing cells."""
    parts = [c.strip() for c in line.split("|")]
    if parts and parts[0] == "":
        parts = parts[1:]
    if parts and parts[-1] == "":
        parts = parts[:-1]
    return parts


def _is_separator(line: str) -> bool:
    """Detect markdown table separator rows like `|---|---|` or `|:--|--:|`.

    Linear-time character scan — replaced an earlier regex that exhibited
    catastrophic backtracking on long non-separator lines.
    """
    s = line.strip()
    if not s:
        return False
    if "|" not in s:
        return False
    if "-" not in s:
        return False
    # Every non-pipe character must be one of -, :, or whitespace
    for ch in s:
        if ch not in "|-: \t":
            return False
    # Must contain at least one stretch of >=3 dashes, e.g. ---
    return "---" in s


# Header content cues that suggest a candidate is actually a variant / patient
# data table rather than artifact text rendered as a pipe-table.
_TABLE_CONTENT_KEYWORDS = (
    "variant",
    "mutation",
    "cdna",
    "protein",
    "amino acid",
    "aachange",
    "patient",
    "carrier",
    "proband",
    "family",
    "subject",
    "case",
    "kindred",
    "affected",
    "unaffected",
    "genotype",
    "phenotype",
    "lqt",
    "syndrome",
    "snv",
    "snp",
    "allele",
    "exon",
    "intron",
    "nucleotide",
    "missense",
    "nonsense",
    "frameshift",
    "splice",
    "deletion",
    "insertion",
    "lqts",
    "p.",
    "c.",
)

_GWAS_OR_ASSAY_HEADERS = {
    "chr",
    "chromosome",
    "bp",
    "position",
    "snp",
    "snv",
    "rsid",
    "locus",
    "ea",
    "eaf",
    "beta",
    "se",
    "p",
    "pvalue",
    "construct",
    "current",
    "activation",
    "tailcurrent",
}

_HEADER_FIELD_KEYWORDS = {
    "gene": ("gene", "symbol"),
    "cdna": ("cdna", "codingdna", "nucleotide", "dna", "ntchange", "cdnachange"),
    "protein": (
        "protein",
        "aminoacid",
        "aachange",
        "proteinchange",
        "mutation",
        "variant",
    ),
    "patient_count": (
        "carrier",
        "patient",
        "proband",
        "subject",
        "individual",
        "family",
        "kindred",
        "ncarrier",
        "number",
        "count",
        "total",
    ),
    "affected": (
        "affected",
        "symptomatic",
        "case",
        "lqt",
        "cpvt",
        "brs",
        "disease",
    ),
    "unaffected": (
        "unaffected",
        "asymptomatic",
        "control",
        "healthy",
        "normal",
        "nonaffected",
    ),
    "uncertain": ("uncertain", "unknown", "equivocal", "borderline", "ambiguous"),
    "phenotype": ("phenotype", "diagnosis", "clinical", "symptom"),
    "clinical_significance": (
        "pathogenic",
        "classification",
        "significance",
        "interpretation",
    ),
}

# Subsets of the affected/unaffected keyword lists that are ambiguous between a
# per-variant clinical count and a case-control *cohort* total. On assay /
# case-control tables, a column whose header matches ONLY these (e.g. "Cases",
# "Controls") is a cohort column, not a per-variant affected/unaffected count;
# mapping it inflates counts (the worst MAE axis). Gated by strict_cohort_labels.
_AMBIGUOUS_COHORT_KEYWORDS = {
    "affected": ("case", "disease"),
    "unaffected": ("control", "healthy", "normal"),
}


def _field_header_match(
    header: str,
    field: str,
    *,
    strict_cohort_labels: bool,
    has_assay_or_gwas_cue: bool,
) -> bool:
    """True when ``header`` maps to ``field``, honoring the cohort-label guard.

    Equivalent to ``any(kw in header for kw in _HEADER_FIELD_KEYWORDS[field])``
    except that, when ``strict_cohort_labels`` is set and the table looks like a
    case-control / assay table, an affected/unaffected header matching ONLY
    ambiguous cohort terms (case/disease, control/healthy/normal) is rejected so
    a cohort total is not assigned as a per-variant count.
    """
    matches = [kw for kw in _HEADER_FIELD_KEYWORDS[field] if kw in header]
    if not matches:
        return False
    if (
        strict_cohort_labels
        and has_assay_or_gwas_cue
        and field in _AMBIGUOUS_COHORT_KEYWORDS
        and all(kw in _AMBIGUOUS_COHORT_KEYWORDS[field] for kw in matches)
    ):
        return False
    return True


def _strict_cohort_labels_enabled() -> bool:
    """Read the strict_cohort_labels flag from Settings (default False)."""
    try:
        from config.settings import get_settings

        return bool(getattr(get_settings(), "strict_cohort_labels", False))
    except Exception:
        return False


_ROW_LEVEL_SUBJECT_KEYWORDS = (
    "patient",
    "patientnumber",
    "patientno",
    "patientid",
    "proband",
    "probandnumber",
    "probandno",
    "probandid",
    "subject",
    "subjectnumber",
    "subjectno",
    "subjectid",
    "caseid",
    "case",
    "casenumber",
    "caseno",
    "family",
    "kindred",
    "individual",
    "adultnumber",
    "childnumber",
    "participant",
    "pedigree",
    "member",
)

_CLINICAL_CONTEXT_KEYWORDS = (
    "age",
    "sex",
    "gender",
    "phenotype",
    "diagnosis",
    "clinical",
    "symptom",
    "syncope",
    "qtc",
    "qtinterval",
    "ecg",
    "arrhythmia",
    "onset",
    "sudden",
    "death",
    "therapy",
    "treatment",
    "schwartz",
    "score",
)

_FAMILY_HISTORY_CONTEXT_TOKENS = (
    "infamily",
    "familyhistory",
    "familial",
    "familymember",
    "affectedrelative",
    "relativewith",
    "pedigreehistory",
    "kindredhistory",
)

_CLINICAL_TIMING_TOKENS = (
    "age",
    "onset",
    "diagnosis",
    "duration",
    "followup",
)

_CLINICAL_MEASURE_TOKENS = (
    "age",
    "case",
    "cancer",
    "diagnosis",
    "duration",
    "event",
    "family",
    "followup",
    "month",
    "onset",
    "relative",
    "year",
)

# Clinical status read out of free text, split by how much the wording can carry.
#
# The unambiguous set states the carrier's phenotype and nothing else, so it may
# match anywhere inside a clinical cell — "asymptomatic LQT2 carrier" is how
# papers write genotype-positive/phenotype-negative, and that row is the
# non-penetrance signal this database exists to capture, so the disease token
# must not outrank it.
_UNAFFECTED_CLINICAL_RE = re.compile(
    r"\b(unaffected|asymptomatic|no symptoms?)\b",
    re.IGNORECASE,
)
# The ambiguous set is ordinary assay and cohort English: "normal trafficking",
# "control-like kinetics", "normal current density". A substring match on these
# used to be enough to record a confirmed unaffected carrier, so a proband whose
# phenotype cell said "aborted cardiac arrest" was stored as affected=0,
# unaffected=1. They are only believed as the *entire* clinical cell.
_UNAFFECTED_WHOLE_CELL_RE = re.compile(
    r"^(unaffected|asymptomatic|control|controls|healthy|normal|"
    r"no symptoms?|none|negative)$",
    re.IGNORECASE,
)
_UNINFORMATIVE_PHENOTYPE_RE = re.compile(
    r"\b(unknown|uncertain|not reported|n/?a)\b", re.IGNORECASE
)
# A negated finding ("normal ECG, no syncope", "without arrhythmia") describes a
# carrier who was assessed and found well, but the wording is a denied symptom
# rather than a status term, so it does not reach the unaffected lexicon above.
# Naming it affected because the cell mentions a condition would invent the
# opposite partition, and enumerating every symptom that can be negated is not a
# rule anyone can keep correct. Such a cell yields null, which is honest.
# The negator must be followed by a word, so an ordinal ("Patient No. 3") is not
# read as a denial.
_NEGATED_FINDING_RE = re.compile(
    r"\b(?:no|without|free\s+of|denies|negative\s+for)\s+[a-z]|\bnone\b",
    re.IGNORECASE,
)


def _phenotype_says_unaffected(phenotype: Optional[str]) -> bool:
    """Whether the phenotype cell states the carrier had no disease phenotype.

    Checked before the affirmative-condition branch on purpose: an unambiguous
    status wins over a disease token in the same cell, because "asymptomatic
    LQT2 carrier" is a genotype-positive/phenotype-negative row, not an affected
    one. Silence is never read as unaffected — an empty or unmapped cell leaves
    both partitions null.
    """
    text = (phenotype or "").strip()
    if not text:
        return False
    return bool(
        _UNAFFECTED_CLINICAL_RE.search(text) or _UNAFFECTED_WHOLE_CELL_RE.match(text)
    )


def _phenotype_names_a_condition(phenotype: Optional[str]) -> bool:
    """Whether the phenotype cell affirmatively names a condition.

    Runs only after ``_phenotype_says_unaffected`` has declined, so a cell that
    both denies findings and names a condition ("normal ECG, no syncope") is
    left null rather than counted as an affected carrier.
    """
    text = (phenotype or "").strip()
    if not text or not re.search(r"[A-Za-z]", text):
        return False
    if _UNINFORMATIVE_PHENOTYPE_RE.search(text) or _NEGATED_FINDING_RE.search(text):
        return False
    return True


# Caption can also be embedded as the first cell of the header row, e.g.
#   | Table S3: rare variants ... | ... |
_EMBEDDED_TABLE_LABEL_RE = re.compile(
    r"\btable\s*(?:s|supp|supplementary)?\s*\d", re.IGNORECASE
)


def _looks_like_variant_table(table: "MarkdownTable") -> bool:
    """Heuristic: is this candidate worth shipping to the LLM router?

    The router only judges *which* table holds variant data, but many markdown
    docs contain dozens of pseudo-tables (text formatted with pipes). Cheap
    Python-side filtering keeps the router prompt small and the LLM call fast.
    """
    if len(table.header_cells) < 2 or len(table.data_lines) < 1:
        return False
    if table.caption and _TABLE_CAPTION_RE.match(table.caption):
        return True
    header_text = " ".join(table.header_cells).lower()
    # Embedded "Table N:" label inside the header itself counts as a real table
    if _EMBEDDED_TABLE_LABEL_RE.search(header_text):
        return True
    return any(kw in header_text for kw in _TABLE_CONTENT_KEYWORDS)


def _normalize_header(value: str) -> str:
    """Normalize a header cell for deterministic field matching."""
    return re.sub(r"[^a-z0-9]+", "", (value or "").lower())


def _column_values(table: MarkdownTable, idx: int, limit: int = 8) -> List[str]:
    values: List[str] = []
    for row in table.data_lines[:limit]:
        cells = _split_pipe_row(row)
        if 0 <= idx < len(cells):
            values.append(cells[idx].strip())
    return values


def _looks_numeric_column(values: List[str]) -> bool:
    non_empty = [v for v in values if v and v not in {"-", "–", "—"}]
    if not non_empty:
        return False
    numeric = sum(1 for value in non_empty if _coerce_int(value) is not None)
    return numeric >= max(1, len(non_empty) // 2)


_EXPLICIT_CARRIER_COUNT_HEADERS = (
    "carrier",
    "carriers",
    "ncarrier",
    "ncarriers",
    "numberofcarriers",
    "noofcarriers",
)

_ROW_IDENTIFIER_HEADERS = (
    "adultnumber",
    "childnumber",
    "patientnumber",
    "patientno",
    "patientid",
    "probandnumber",
    "probandno",
    "probandid",
    "subjectnumber",
    "subjectno",
    "subjectid",
    "individualnumber",
    "individualno",
    "individualid",
    "participantnumber",
    "participantno",
    "participantid",
    "casenumber",
    "caseno",
    "caseid",
    "familynumber",
    "familyno",
    "familyid",
    "kindrednumber",
    "kindredno",
    "kindredid",
    "pedigreenumber",
    "pedigreeno",
    "pedigreeid",
    "sampleid",
)

_POPULATION_FREQUENCY_HEADERS = (
    "maf",
    "minorallelefrequency",
    "allelefrequency",
    "alternateallelefrequency",
    "effectallelefrequency",
    "eaf",
    "frequency",
    "freq",
    "occurrence",
    "occurrences",
    "numberofoccurrences",
    "noofoccurrences",
    "allelecount",
    "alleles",
    "exac",
    "gnomad",
    "1000genomes",
    "topmed",
    "clinvar",
    "dbsnp",
)

# In-silico pathogenicity predictors. A numeric column under one of these is a
# computational score, never a human carrier/person count (extraction rule 3).
# Tokens are matched as substrings of normalized (lowercased, punctuation-
# stripped) headers, so keep them long enough to avoid colliding with real
# count headers (e.g. "sift", "cadd", "revel" are safe; bare "score" is not and
# stays out — it is a clinical-context cue, e.g. Schwartz score).
_PREDICTION_SCORE_HEADERS = (
    "sift",
    "polyphen",
    "cadd",
    "revel",
    "provean",
    "mutationtaster",
    "mutationassessor",
    "alphamissense",
    "fathmm",
    "phylop",
    "phastcons",
    "gerp",
    "metasvm",
    "metalr",
    "primateai",
    "condel",
    "insilico",
)

_DENOMINATOR_HEADERS_WHEN_CARRIER_PRESENT = (
    "totalcase",
    "totalcases",
    "case",
    "cases",
    "totalcontrol",
    "totalcontrols",
    "control",
    "controls",
    "totalsample",
    "samplesize",
    "cohortsize",
    "screened",
    "screenedn",
    "ntested",
    "numbertested",
)


def _has_explicit_carrier_header(headers: List[str]) -> bool:
    return any(
        any(token in header for token in _EXPLICIT_CARRIER_COUNT_HEADERS)
        for header in headers
    )


def _is_row_identifier_header(header: str) -> bool:
    return any(token in header for token in _ROW_IDENTIFIER_HEADERS)


def _is_population_frequency_header(header: str) -> bool:
    return any(token in header for token in _POPULATION_FREQUENCY_HEADERS)


def _is_prediction_score_header(header: str) -> bool:
    return any(token in header for token in _PREDICTION_SCORE_HEADERS)


def _is_variant_annotation_header(header: str) -> bool:
    """True for columns that annotate a variant rather than count people.

    Population-frequency allele counts (gnomAD/ExAC/TOPMed…) and in-silico
    pathogenicity predictors (SIFT/PolyPhen/CADD/REVEL…) are variant metadata,
    never carrier/person counts (extraction rules 2-3).
    """
    return _is_population_frequency_header(header) or _is_prediction_score_header(
        header
    )


def _is_denominator_header_when_carrier_present(
    header: str, normalized_headers: List[str]
) -> bool:
    if not _has_explicit_carrier_header(normalized_headers):
        return False
    return any(token in header for token in _DENOMINATOR_HEADERS_WHEN_CARRIER_PRESENT)


def _is_non_variant_count_header(header: str, normalized_headers: List[str]) -> bool:
    """Reject numeric columns that are visibly not per-variant carrier counts."""
    row_identifier = _is_row_identifier_header(header)
    # ``Family (No.)`` is usually an identifier, but mutation-list tables can
    # also use it as a family count when a separate serial-number column is the
    # actual row identifier.  Resolve that role from the table schema rather
    # than from any paper, gene, or caption identity.
    if (
        row_identifier
        and header in {"familyno", "familynumber"}
        and any(
            candidate in {"no", "number", "rowno", "serialno"}
            for candidate in normalized_headers
        )
    ):
        row_identifier = False
    # Clinical-characteristics tables often put a mutation next to numeric
    # phenotype summaries such as age at diagnosis, affected relatives, and
    # early-onset events.  A fuzzy ``cases`` match must not turn those measures
    # into variant carriers.
    # Require the table-level combination here rather than blacklisting every
    # ``cases`` column: ordinary per-variant case tables remain valid.
    clinical_characteristics_table = any(
        any(token in candidate for token in _FAMILY_HISTORY_CONTEXT_TOKENS)
        for candidate in normalized_headers
    ) and any(
        any(token in candidate for token in _CLINICAL_TIMING_TOKENS)
        for candidate in normalized_headers
    )
    clinical_measure = clinical_characteristics_table and any(
        token in header for token in _CLINICAL_MEASURE_TOKENS
    )
    return (
        row_identifier
        or _is_population_frequency_header(header)
        or _is_denominator_header_when_carrier_present(header, normalized_headers)
        or clinical_measure
    )


def _has_header_keyword(headers: List[str], keywords: tuple[str, ...]) -> bool:
    return any(any(kw in header for kw in keywords) for header in headers)


def _looks_like_row_level_clinical_list(
    table: MarkdownTable, mapping: Dict[str, int], has_assay_or_gwas_cue: bool
) -> bool:
    """Detect clinical mutation-list tables where each row is one carrier.

    Many curated clinical papers list one proband/patient/family per row with
    mutation + phenotype columns but no explicit ``n`` column. In that shape,
    the row itself is the carrier count. Keep the rule conservative so variant
    definition tables, assay tables, and GWAS summaries do not become fake
    patient counts.
    """
    if not any(k in mapping for k in ("cdna", "protein")):
        return False
    if any(k in mapping for k in ("patient_count", "affected", "unaffected")):
        return False

    headers = [_normalize_header(c) for c in table.header_cells]
    caption = _normalize_header(table.caption or "")
    subject_cue = _has_header_keyword(headers, _ROW_LEVEL_SUBJECT_KEYWORDS)
    clinical_cue = _has_header_keyword(headers, _CLINICAL_CONTEXT_KEYWORDS) or any(
        kw in caption for kw in _CLINICAL_CONTEXT_KEYWORDS
    )

    if has_assay_or_gwas_cue and not subject_cue:
        return False
    # An annotation/summary table keyed by variant-annotation columns
    # (population-frequency allele counts like gnomAD/ExAC, or in-silico
    # predictors like SIFT/PolyPhen/CADD/REVEL) with no subject column is
    # variant annotation, not a one-carrier-per-row clinical list. A bare
    # clinical caption/"score" cue must not mint a carrier for every annotated
    # row (e.g. PMID 33013630 Table 1: gnomAD allele count + SIFT/PolyPhen,
    # zero patient columns; or a "REVEL score | CADD score" prediction table).
    has_annotation_header = any(_is_variant_annotation_header(h) for h in headers)
    if has_annotation_header and not subject_cue:
        return False
    return subject_cue or clinical_cue


def _infer_column_mapping_from_headers(
    table: MarkdownTable,
    *,
    strict_cohort_labels: bool = False,
    target_gene: Optional[str] = None,
) -> Optional[Dict[str, int]]:
    """Infer a table mapping without an LLM when headers/data are unambiguous."""
    mapping: Dict[str, int] = {}
    normalized_headers = [_normalize_header(c) for c in table.header_cells]

    # Avoid obvious GWAS/assay tables unless their data columns contain actual
    # HGVS-like variant notation. This keeps "AA/n" allele-frequency tables out.
    has_assay_or_gwas_cue = any(h in _GWAS_OR_ASSAY_HEADERS for h in normalized_headers)

    for idx, header in enumerate(normalized_headers):
        values = _column_values(table, idx)
        protein_by_data = any(_normalize_protein(v) for v in values)
        cdna_by_data = any(_normalize_cdna(v) for v in values)
        non_variant_count_header = _is_non_variant_count_header(
            header, normalized_headers
        )

        if "gene" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["gene"]
        ):
            # Substring matching on the normalized header alone is not enough:
            # "Generation" and "General comments" both contain "gene". A false
            # gene column is worse than none, because it (a) suppresses the
            # whole-table caption reject and (b) makes the per-row gene filter
            # compare variant rows against values like "F1"/"F2" and drop every
            # one — a table with a Generation column lost ALL of its own rows.
            if _column_values_look_like_genes(values):
                mapping["gene"] = idx

        if (
            not non_variant_count_header
            and "unaffected" not in mapping
            and _field_header_match(
                header,
                "unaffected",
                strict_cohort_labels=strict_cohort_labels,
                has_assay_or_gwas_cue=has_assay_or_gwas_cue,
            )
        ):
            mapping["unaffected"] = idx
            continue

        if (
            not non_variant_count_header
            and "affected" not in mapping
            and _field_header_match(
                header,
                "affected",
                strict_cohort_labels=strict_cohort_labels,
                has_assay_or_gwas_cue=has_assay_or_gwas_cue,
            )
        ):
            mapping["affected"] = idx
            continue

        if "cdna" not in mapping and (
            any(kw in header for kw in _HEADER_FIELD_KEYWORDS["cdna"]) or cdna_by_data
        ):
            mapping["cdna"] = idx
            continue

        if "protein" not in mapping and (
            any(kw in header for kw in _HEADER_FIELD_KEYWORDS["protein"])
            or protein_by_data
        ):
            mapping["protein"] = idx
            continue

        if "uncertain" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["uncertain"]
        ):
            mapping["uncertain"] = idx
            continue

        if "phenotype" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["phenotype"]
        ):
            mapping["phenotype"] = idx
            continue

        if "clinical_significance" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["clinical_significance"]
        ):
            mapping["clinical_significance"] = idx
            continue

        if (
            not non_variant_count_header
            and "patient_count" not in mapping
            and any(kw in header for kw in _HEADER_FIELD_KEYWORDS["patient_count"])
        ):
            if _looks_numeric_column(values):
                mapping["patient_count"] = idx

    # Converters can erase the header of a rowspan-like gene column. Recover it
    # from column structure so parse_routed_table can forward-fill groups and
    # keep off-target rows out of any target-gene extraction.
    if "gene" not in mapping:
        inferred_gene_idx = _infer_unnamed_gene_column(
            table.header_cells,
            table.data_lines,
            target_gene=target_gene,
        )
        if inferred_gene_idx is not None:
            mapping["gene"] = inferred_gene_idx

    # If a notation column has a generic "mutation/variant" header, data decides
    # whether it is cDNA or protein. Prefer explicit cDNA/protein when present.
    for idx, _header in enumerate(normalized_headers):
        if "cdna" in mapping and "protein" in mapping:
            break
        values = _column_values(table, idx)
        if "cdna" not in mapping and any(_normalize_cdna(v) for v in values):
            mapping["cdna"] = idx
        if "protein" not in mapping and any(_normalize_protein(v) for v in values):
            mapping["protein"] = idx

    if _looks_like_row_level_clinical_list(table, mapping, has_assay_or_gwas_cue):
        mapping["patient_count"] = _INFER_ROW_PATIENT_COUNT

    has_notation = any(k in mapping for k in ("cdna", "protein"))
    has_count = any(k in mapping for k in ("patient_count", "affected", "unaffected"))
    if not (has_notation and has_count):
        return None
    if has_assay_or_gwas_cue and not any(
        any(
            _normalize_cdna(v) or _normalize_protein(v)
            for v in _column_values(table, i)
        )
        for i in range(len(table.header_cells))
    ):
        return None
    return mapping


# A separator-less pipe block needs at least this many consecutive `|` rows
# (header + >=2 data rows) before we treat it as a table. Office-format
# supplements (.doc/.docx) converted to text emit pipe-delimited rows with no
# `|---|` separator at all, so requiring one made every such table invisible.
_MIN_BORDERLESS_ROWS = 3


def _is_header_continuation(cells: List[str]) -> bool:
    """Is this row a wrapped fragment of the header above it, not a data row?

    Office-to-text conversion wraps a long header cell across physical lines,
    leaving every *other* cell blank::

        |Region |Nucleotide |Variant |Mutation Type |Location |No. of |
        |       |           |        |              |         |patient|
        |       |           |        |              |         |s      |

    A data row fills most of its cells; a wrapped header fragment fills one or
    two. Callers only consult this *before* the first real data row, so a
    mostly-blank continuation of a data row can never reach it.
    """
    if not cells:
        return False
    non_blank = sum(1 for c in cells if c.strip())
    if non_blank == 0:
        return False
    # An in-band gene divider is, by definition, a row with exactly one
    # populated cell — the same shape as a wrapped header fragment. In the
    # separator-less branch it was therefore concatenated into the header
    # ("Nucleotide" + "BRCA1" -> "NucleotideBRCA1") and the table never
    # partitioned, so the second gene's rows leaked into the first gene's run.
    # A cell that reads as a bare gene label is a divider, not a header.
    if _gene_section_divider(cells) is not None:
        return False
    return non_blank <= max(1, len(cells) // 3)


def _join_wrapped_header(
    lines: List[str], header_idx: int, stop_idx: int
) -> Tuple[List[str], int]:
    """Rebuild a header split across physical lines.

    Returns ``(header_cells, next_index)`` where ``next_index`` is the first
    line that is *not* part of the header.

    Fragments are concatenated with no separator. Whether the break was
    mid-word (``patient`` + ``s``) or between words (``No. of`` + ``patients``)
    is undecidable from the text alone, but the choice is immaterial for
    matching: ``_normalize_header`` keeps only alphanumerics, so both
    ``"No. ofpatients"`` and ``"No. of patient s"`` normalize to
    ``noofpatients``. No-space join is preferred so a mid-word break — the
    common case, since conversion wraps on width — stays one readable token in
    the preview shipped to the router.
    """
    header_cells = _split_pipe_row(lines[header_idx])
    idx = header_idx + 1
    while idx < stop_idx:
        cells = _split_pipe_row(lines[idx])
        if not cells or not _is_header_continuation(cells):
            break
        for col, frag in enumerate(cells):
            frag = frag.strip()
            if not frag or col >= len(header_cells):
                continue
            header_cells[col] = (header_cells[col].rstrip() + frag).strip()
        idx += 1
    return header_cells, idx


def _lookback_caption(lines: List[str], header_index: int) -> Optional[str]:
    """Return the caption above a table header, or None when unanchored.

    Two shapes occur in practice and both must be captured, because the caption
    is the ONLY gene-scoping signal for a table with no gene column:

        (a) inline      "Table 1. Mutations in BRCA1 gene"
        (b) PMC/XML→md  "### Table 1" / blank / "Mutations in BRCA1 gene"

    Shape (b) is what the PMC converter emits and it was invisible to the
    original lookback, which required `_TABLE_CAPTION_RE` to match the first
    non-blank line and then `break`ed unconditionally. Across a 400-paper corpus
    sample that scanner captured 0 captions out of 330 tables, which silently
    disabled every caption-derived guard — see the cross-gene reject in
    `parse_routed_table`.

    Deliberately anchored, not greedy: a descriptive line is accepted only when
    it sits under a "Table N" label. Treating arbitrary preceding prose as a
    caption lets a body paragraph naming another gene reject a legitimate table,
    which is a measured recall loss, not a hypothetical one.
    """
    preceding: List[tuple[int, str]] = []
    for back in range(1, 8):
        idx = header_index - back
        if idx < 0:
            break
        candidate = lines[idx].strip()
        if not candidate:
            continue
        if candidate.startswith("|"):
            break
        preceding.append((idx, candidate))
        if len(preceding) >= 3:
            break

    for pos, (_, candidate) in enumerate(preceding):
        if _TABLE_CAPTION_RE.match(candidate):
            return candidate
        heading = _TABLE_HEADING_RE.match(candidate)
        if heading:
            label = heading.group(1).strip()
            descriptive = preceding[pos - 1][1] if pos > 0 else ""
            if descriptive and not _TABLE_HEADING_RE.match(descriptive):
                return f"{label}. {descriptive}"
            return label
    return None


def enumerate_markdown_tables(
    text: str, *, only_variant_like: bool = True
) -> List[MarkdownTable]:
    """Walk a markdown document and slice out every pipe-table.

    A "table" here is a stretch of consecutive `|`-rows containing a separator
    row (`|---|---|`). We capture the caption (immediately preceding line that
    matches "Table N." or starts with "Table" + a digit) plus all data rows.

    With ``only_variant_like=True`` (default) we drop candidates that have no
    "Table N." caption *and* no variant-ish keyword in their header — these are
    almost always paragraph text that happens to contain `|` separators in the
    PMC XML→markdown conversion.
    """
    lines = text.splitlines()
    cum_offsets = [0]
    for line in lines:
        cum_offsets.append(cum_offsets[-1] + len(line) + 1)  # +1 for newline

    tables: List[MarkdownTable] = []
    i = 0
    counter = 0
    while i < len(lines):
        line = lines[i]
        # Detect a separator → step back one line to find the header
        if _is_separator(line) and i > 0:
            header_line = lines[i - 1]
            if not header_line.strip().startswith("|"):
                i += 1
                continue
            header_cells = _split_pipe_row(header_line)
            if len(header_cells) < 2:
                i += 1
                continue
            # Walk forward collecting data rows until we hit a non-pipe line
            data_lines: List[str] = []
            j = i + 1
            while j < len(lines):
                row = lines[j]
                if not row.strip().startswith("|") and not row.strip().startswith("|"):
                    if not row.strip():
                        # blank line ends the table (markdown convention)
                        break
                    if not row.lstrip().startswith("|"):
                        break
                if _is_separator(row):
                    j += 1
                    continue
                if row.strip().startswith("|"):
                    data_lines.append(row)
                    j += 1
                    continue
                break

            caption = _lookback_caption(lines, i - 1)

            char_start = cum_offsets[max(0, i - 1)]
            char_end = cum_offsets[min(len(lines), j)]
            counter += 1
            tables.append(
                MarkdownTable(
                    table_id=f"T{counter}",
                    caption=caption,
                    header_line=header_line,
                    header_cells=header_cells,
                    data_lines=data_lines,
                    char_start=char_start,
                    char_end=char_end,
                )
            )
            i = j
            continue

        # Separator-less pipe block. Office-format supplements (.doc/.docx)
        # converted to text emit pipe-delimited rows with no `|---|` row at
        # all, so the branch above never fired and the table was invisible to
        # the router — the whole document scored as "no usable variant tables".
        # Only claim a run that contains NO separator, so a bordered table is
        # still handled by the branch above (which needs to see the separator
        # first) and can never be emitted twice.
        if line.strip().startswith("|") and not _is_separator(line):
            j = i
            saw_separator = False
            while j < len(lines) and lines[j].strip().startswith("|"):
                if _is_separator(lines[j]):
                    saw_separator = True
                    break
                j += 1
            run_len = j - i
            if saw_separator or run_len < _MIN_BORDERLESS_ROWS:
                i += 1
                continue

            header_cells, data_start = _join_wrapped_header(lines, i, j)
            data_lines = [
                row for row in lines[data_start:j] if row.strip().startswith("|")
            ]
            if len(header_cells) < 2 or not data_lines:
                i += 1
                continue

            # Same anchored lookback as the bordered branch above — this
            # branch was added in the same series and kept the pre-fix logic
            # (first non-blank line, caption regex only, unconditional break),
            # so borderless tables stayed captionless and unscoped.
            caption = _lookback_caption(lines, i)

            counter += 1
            tables.append(
                MarkdownTable(
                    table_id=f"T{counter}",
                    caption=caption,
                    header_line=lines[i],
                    header_cells=header_cells,
                    data_lines=data_lines,
                    char_start=cum_offsets[i],
                    char_end=cum_offsets[min(len(lines), j)],
                )
            )
            i = j
            continue
        i += 1

    if only_variant_like:
        tables = [t for t in tables if _looks_like_variant_table(t)]
    return tables


def build_router_prompt(
    tables: List[MarkdownTable],
    gene_symbol: str,
    *,
    sample_rows: int = 3,
) -> str:
    """Compact JSON-friendly inventory of every detected table.

    The router only sees: caption, header, and a few sample rows per table — plus
    the gene we are extracting. It never sees the full paper text. Sample rows
    are critical for tables exported from Excel/CSV where the markdown header is
    just `Unnamed: 1 | Unnamed: 2 | ...` (a common PMC supplementary artifact);
    the router then infers column types from the data values.
    """
    lines: List[str] = []
    for t in tables:
        block = [f"### {t.table_id}"]
        if t.caption:
            block.append(f"caption: {t.caption.strip()[:240]}")
        block.append(f"header: {t.header_preview(max_chars=320)}")
        previews = t.data_preview(n_rows=sample_rows, max_chars=320)
        for k, row in enumerate(previews, 1):
            block.append(f"row{k}: {row}")
        lines.append("\n".join(block))
    inventory = "\n\n".join(lines)

    instructions = f"""You are routing markdown tables for downstream deterministic
extraction of variant data. You will NOT extract any variants yourself. The
target gene of interest is {gene_symbol}, but multi-gene tables (LQTS panels,
arrhythmia panels, supplementary variant lists across many genes) ALSO qualify
as long as they have a "Gene" column — the parser will filter rows by gene.

Output each qualifying table's column→field mapping using 0-based column
indices into the header row.

Allowed field names (omit a field if no column matches):
  - gene              — column listing the gene symbol (KCNH2, SCN5A, etc.)
  - cdna              — column with c. or HGVS coding-DNA notation
  - protein           — column with p. or short-form amino acid notation
  - patient_count     — explicit variant-specific carriers / patients / probands
                        / "n detected" counts
  - affected          — affected/symptomatic carrier count
  - unaffected        — unaffected/asymptomatic carrier count
  - uncertain         — equivocal/borderline carrier count
  - phenotype         — clinical phenotype text
  - clinical_significance — pathogenicity classification

A table qualifies if it has at least one of {{cdna, protein}} AND at least one
of {{patient_count, affected, unaffected}}. A clinical mutation-list table with
one patient/proband/family/subject/case per row also qualifies even if it has no
explicit count column; in that case set "patient_count": -1 to tell the parser
to count one carrier per row. Multi-gene panel tables qualify even if the
preview rows show non-{gene_symbol} variants. Tables that list functional
assays, drug screens, primer sequences, in silico predictions,
GWAS/association statistics, allele-frequency summaries, lead SNPs, or columns
such as Locus/SNV/CHR/BP/EA/AA/EAF/beta/SE/p/n do NOT qualify — skip them.
In those tables, AA usually means alternate allele and n is cohort size, not a
patient/carrier count.
Do NOT map denominator, row-ID, or population-frequency columns to patient_count:
"Total case(s)" is a screened/case denominator when a "Carrier(s)" column is
present; use "Carrier(s)" instead. "Adult number", "child number", "patient ID",
and "case no." are row identifiers; for one-patient-per-row clinical lists set
"patient_count": -1. "MAF", "allele frequency", "No. of occurrences",
ExAC/gnomAD/1000 Genomes counts, and allele denominators are population data,
not clinical carrier counts.

Return strict JSON. No prose. No markdown fences. Schema:

{{
  "variant_tables": [
    {{
      "table_id": "T1",
      "column_mapping": {{
        "cdna": 0,
        "protein": 1,
        "patient_count": 2,
        "affected": 3,
        "unaffected": 4
      }},
      "confidence": 0.0_to_1.0,
      "notes": "optional short reason"
    }}
  ]
}}

If no table qualifies, return {{"variant_tables": []}} exactly.

Tables to route:

{inventory}
"""
    return instructions.strip()


def parse_router_response(raw: str) -> List[RoutedTable]:
    """Parse the router's JSON output into RoutedTable objects.

    Tolerates fenced JSON, leading prose, and trailing chatter. Drops mappings
    that violate the field allowlist or use non-integer indices.
    """
    if not raw:
        return []

    text = raw.strip()
    # Strip ```json fences if present
    if text.startswith("```"):
        text = re.sub(r"^```[a-zA-Z0-9_]*\s*", "", text)
        text = re.sub(r"\s*```\s*$", "", text)
    # Locate the first decodable JSON object so prose-prefixed responses work
    # even when the prose itself contains braces. raw_decode ignores trailing
    # prose after the object; rfind('}') was unsafe when trailing analysis
    # embedded braces (e.g. "I picked table {ID: 5}").
    candidate_starts = [idx for idx, char in enumerate(text) if char == "{"]
    if not candidate_starts:
        return []

    data = None
    decoder = json.JSONDecoder()
    last_error = None
    for start in candidate_starts:
        try:
            parsed, _end = decoder.raw_decode(text[start:])
        except json.JSONDecodeError as e:
            last_error = e
            continue
        if isinstance(parsed, dict):
            data = parsed
            break

    if data is None:
        if last_error:
            logger.warning("table_router: JSON parse failed: %s", last_error)
        return []

    items = data.get("variant_tables") or []
    out: List[RoutedTable] = []
    for item in items:
        tid = item.get("table_id")
        mapping = item.get("column_mapping") or {}
        if not tid or not isinstance(mapping, dict):
            continue
        clean: Dict[str, int] = {}
        for field_name, idx in mapping.items():
            if field_name not in _KNOWN_FIELDS:
                continue
            try:
                parsed_idx = int(idx)
            except (TypeError, ValueError):
                continue
            if parsed_idx < 0 and not (
                field_name == "patient_count" and parsed_idx == _INFER_ROW_PATIENT_COUNT
            ):
                continue
            clean[field_name] = parsed_idx
        if not clean:
            continue
        # Must have at least one notation column AND one count column
        has_notation = any(k in clean for k in ("cdna", "protein"))
        has_count = any(k in clean for k in ("patient_count", "affected", "unaffected"))
        if not (has_notation and has_count):
            continue
        out.append(
            RoutedTable(
                table_id=str(tid),
                column_mapping=clean,
                confidence=item.get("confidence"),
                notes=item.get("notes"),
            )
        )
    return out


def _coerce_int(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    s = str(value).strip()
    # Strip footnote markers, percent signs, and trailing parenthetical notes
    s = re.sub(r"[  ]", "", s)
    s = re.sub(r"\([^)]*\)", "", s).strip()
    s = s.rstrip("*†‡§¶")
    try:
        return int(s)
    except ValueError:
        try:
            return int(float(s))
        except ValueError:
            return None


_AA_TOKEN = r"(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
_PROTEIN_VARIANT_RE = re.compile(
    r"^(?:p\.)?(?:"
    rf"{_AA_TOKEN}\d{{1,4}}"
    rf"(?:_{_AA_TOKEN}\d{{1,4}})?"
    r"(?:"
    rf"(?:{_AA_TOKEN}|[ACDEFGHIKLMNPQRSTVWY])?fs(?:X|\*)?\d*"
    rf"|{_AA_TOKEN}"
    r"|[ACDEFGHIKLMNPQRSTVWY*X?]"
    r"|del\d*"
    r"|dup"
    rf"|ins(?:{_AA_TOKEN}|[ACDEFGHIKLMNPQRSTVWY]+)?"
    r")"
    r"|\d{1,4}_\d{1,4}ins(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY]+)"
    r")$",
    re.IGNORECASE,
)

_CDNA_VARIANT_RE = re.compile(
    r"^c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
    r"|^c\.\d+(?:[+-]\d+)?(?:del|dup|ins)[ACGT]*$"
    r"|^c\.\d+(?:_\d+)?(?:del|dup|ins)[ACGT]*$",
    re.IGNORECASE,
)


def _normalize_cdna(value: str) -> Optional[str]:
    s = value.strip().replace(" ", "")
    if not s or s.lower() in {"-", "na", "n/a", "none", "."}:
        return None
    if not s.lower().startswith("c."):
        s = "c." + s
    if not _CDNA_VARIANT_RE.match(s):
        return None
    return s


def _normalize_protein(value: str) -> Optional[str]:
    s = value.strip().replace(" ", "")
    if not s or s.lower() in {"-", "na", "n/a", "none", "."}:
        return None
    if s.endswith("*"):
        without_footnote = s[:-1]
        if without_footnote and _PROTEIN_VARIANT_RE.match(without_footnote):
            s = without_footnote
    if not _PROTEIN_VARIANT_RE.match(s):
        return None
    return s


def _split_variant_cell(value: Optional[str]) -> List[str]:
    """Split cells that list parallel cDNA/protein pairs."""
    if not value:
        return []
    return [
        part.strip()
        for part in re.split(r"\s*(?:,|;)\s*", value)
        if part and part.strip()
    ]


def _normalize_cdna_values(value: Optional[str]) -> List[str]:
    return [
        normalized
        for part in _split_variant_cell(value)
        if (normalized := _normalize_cdna(part))
    ]


def _normalize_protein_values(value: Optional[str]) -> List[str]:
    return [
        normalized
        for part in _split_variant_cell(value)
        if (normalized := _normalize_protein(part))
    ]


_GENE_SYMBOL_CELL_RE = re.compile(r"^[A-Z][A-Z0-9-]{2,11}$")
_GENE_TOKEN_PROTEIN_LIKE_RE = re.compile(
    r"^(?:p\.)?(?:(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|"
    r"Phe|Pro|Ser|Thr|Trp|Tyr|Val)|[ACDEFGHIKLMNPQRSTVWY])\d{1,5}"
    r"(?:[A-Z*?]|[a-z]{2}|fs|del|dup|ins)",
    re.IGNORECASE,
)
_GENE_SYMBOL_IGNORE = {
    "ALL",
    "DNA",
    "RNA",
    "ECG",
    "QT",
    "QTC",
    "LQT",
    "LQTS",
    "CPVT",
    "SCD",
    "WT",
    "NA",
    "PMID",
    "EXON",
    "INTRON",
    "DELETION",
    "DUPLICATION",
    "INSERTION",
    "GRCH37",
    "GRCH38",
    "HG19",
    "HG38",
    # Mutation-type and protein-domain COLUMN VALUES. `_GENE_SYMBOL_CELL_RE`
    # matches any 3-12 character uppercase token, so a `Mutation Type` cell
    # reading "Missense" and a `Location` cell reading "N-Terminal" both
    # registered as "some other gene" — and because
    # `_row_has_off_target_gene_without_target` then saw an off-target token
    # with no target token, it discarded EVERY row of the table. Any variant
    # table carrying a mutation-type or domain column but no mapped Gene column
    # lost all of its rows; PMID 26496715 (99 count-bearing cardiac gold rows)
    # is the type specimen.
    #
    # This stays a denylist rather than a "must look like a gene symbol"
    # allowlist because digit-free real symbols (LMNA, TTN, ENG, MYH7's
    # all-alpha peers) must keep tripping the guard — that is the cross-gene
    # contamination defense pinned by
    # tests/unit/test_extraction_table_parser.py. The vocabulary below is a
    # closed class (the value space of two column kinds), not open English.
    "MISSENSE",
    "NONSENSE",
    "FRAMESHIFT",
    "SYNONYMOUS",
    "SILENT",
    "TRUNCATING",
    "SPLICE",
    "SPLICING",
    "INFRAME",
    "INDEL",
    "STOPGAIN",
    "STOPLOSS",
    "REGULATORY",
    "INTRONIC",
    "EXONIC",
    "UTR",
    "PROMOTER",
    "TERMINAL",
    "N-TERMINAL",
    "C-TERMINAL",
    "DOMAIN",
    "LINKER",
    "PORE",
    "HELIX",
    "TRANSMEMBRANE",
    "CYTOPLASMIC",
    "EXTRACELLULAR",
    "TYPE",
    "REGION",
    "LOCATION",
}


def _target_gene_tokens(gene_symbol: str) -> set[str]:
    gene = (gene_symbol or "").strip().upper()
    aliases = get_gene_aliases(gene, include_query_aliases=True)
    return {alias.upper() for alias in aliases if alias} | {gene}


# Query aliases that also read as ordinary clinical prose. "FH" is "family
# history" far more often than it is familial hypercholesterolaemia, and this
# scope drives a whole-table reject, so an ambiguous alias alone must not scope
# a caption.
_CAPTION_AMBIGUOUS_ALIASES = {"FH", "IRIS"}


def _query_only_aliases(gene: str) -> tuple:
    """Aliases a gene has ONLY as discovery/query terms, not as core synonyms."""
    from utils.gene_metadata import BUILTIN_GENE_METADATA

    meta = BUILTIN_GENE_METADATA.get(str(gene).strip().upper())
    if not meta:
        return ()
    core = {a.upper() for a in (meta.aliases or ())}
    return tuple(a for a in (meta.query_aliases or ()) if a.upper() not in core)


def _caption_gene_scope(caption: Optional[str]) -> set[str]:
    """Return target genes explicitly implied by a table caption."""
    if not caption:
        return set()
    # Query aliases mostly EARN their place in a caption: "LQT1" reliably means
    # KCNQ1, "CPVT1" means RYR2, and those captions are the only gene signal a
    # table without a gene column has.
    #
    # A few are ordinary prose, though, and this scope drives a whole-table
    # reject. LDLR's query alias "FH" matches "probands with FH of sudden
    # cardiac death" (family history), which scoped a cardiac mutation table to
    # {LDLR} and discarded all of it. BRCA1's "IRIS" is the same shape. So a
    # gene is admitted on a query alias only when at least one NON-ambiguous
    # alias matched.
    scope: set[str] = set()
    for gene in known_gene_aliases(include_query_aliases=True):
        if gene_alias_regex(gene, include_query_aliases=False).search(caption):
            scope.add(gene)
            continue
        if not gene_alias_regex(gene, include_query_aliases=True).search(caption):
            continue
        informative = [
            alias
            for alias in _query_only_aliases(gene)
            if alias.upper() not in _CAPTION_AMBIGUOUS_ALIASES
            and re.search(
                rf"(?<![A-Za-z0-9]){re.escape(alias)}(?![A-Za-z0-9])",
                caption,
                re.IGNORECASE,
            )
        ]
        if informative:
            scope.add(gene)
    return scope


# A divider cell may carry benign qualifiers around the gene symbol; real
# corpus dividers include "BRCA1 Gene" and "BRCA1 gene mutations", so exact
# match is too strict and reopens the cross-gene leak.
_DIVIDER_BENIGN_WORDS = {
    "gene",
    "genes",
    "mutation",
    "mutations",
    "variant",
    "variants",
    "alteration",
    "alterations",
    "sequence",
    "carriers",
}

# ...but a cell that NEGATES or subsets the gene is not a divider at all.
# "BRCA1 negative subgroup" describes patients WITHOUT the gene's mutations;
# reading it as a divider flips the section gene for the rest of the table and
# silently deletes the target gene's remaining rows.
_DIVIDER_NEGATION_WORDS = {
    "negative",
    "without",
    "non",
    "excluding",
    "excluded",
    "wild",
    "wt",
    "vs",
    "versus",
    "free",
    "absent",
    "no",
    "not",
    "control",
    "controls",
    "unrelated",
}


def _looks_like_gene_symbol(token: str) -> bool:
    """Open-vocabulary gene-symbol shape test for a single token."""
    token = token.strip().strip("*: ").upper()
    if not (2 < len(token) <= 15):
        return False
    if not token[0].isalpha() or not token.replace("-", "").isalnum():
        return False
    if re.fullmatch(r"[A-Z]{1,2}\d{1,3}", token):  # F1, P040 — sample IDs
        return False
    return sum(c.isalpha() for c in token) >= 3


def _gene_section_divider(
    cells: List[str],
    *,
    target_gene: Optional[str] = None,
    allow_unknown: bool = False,
) -> Optional[str]:
    """Return the gene named by an in-band section-divider row, else None.

    Papers that report two genes in one table partition it with a row whose only
    populated cell is a bare gene symbol::

        | BRCA1 |  |  |  |
        | 710 C > T | c.591C > T | C197C | ... |
        | BRCA2 |  |  |  |
        | 1093A > C | c.865A > C | N289H | ... |

    Such a table has no gene *column*, so without recognising the divider every
    row inherits the run's target gene and half the table is attributed to the
    wrong one. PMID 21232165 Table 3 is the reference case.
    """
    populated = [c.strip() for c in cells if c.strip()]
    if len(populated) != 1:
        return None
    token = populated[0].strip().strip("*: ")
    if not token or len(token) > 40:
        return None

    words = [w for w in re.split(r"[\s,;/()\-]+", token.lower()) if w]
    if any(word in _DIVIDER_NEGATION_WORDS for word in words):
        return None

    genes = {
        gene
        for gene in known_gene_aliases(include_query_aliases=False)
        if gene_alias_regex(gene, include_query_aliases=False).search(token)
    }
    if len(genes) > 1:
        return None
    if not genes:
        # Off-roster genes must be able to partition a table too. Only 14 genes
        # are registered, so without this a `| TTN |` or `| LMNA |` row after a
        # `| BRCA1 |` row is not a divider, `section_gene` stays BRCA1, and the
        # TTN run loses every one of its rows while BRCA1 absorbs them.
        #
        # Two anchors keep this from becoming an open-vocabulary guess: the
        # token is the gene we are extracting for, or this table has ALREADY
        # demonstrated the bare-gene-cell partition idiom with a known gene.
        # A gene-shaped token with neither anchor is left alone.
        candidate = token.strip().strip("*: ")
        if not _looks_like_gene_symbol(candidate):
            return None
        if target_gene and candidate.upper() == target_gene.strip().upper():
            return candidate.upper()
        if allow_unknown and len(words) == 1:
            return candidate.upper()
        return None
    gene = next(iter(genes))

    gene_lower = gene.lower()
    leftover = [
        word
        for word in words
        if word not in _DIVIDER_BENIGN_WORDS
        and word not in gene_lower
        and gene_lower not in word
    ]
    if leftover:
        return None
    return gene


def _cell_mentions_target_gene(value: Optional[str], gene_symbol: str) -> bool:
    if not value or not gene_symbol:
        return False
    return bool(gene_alias_regex(gene_symbol, include_query_aliases=True).search(value))


def _gene_symbol_tokens(value: Optional[str]) -> set[str]:
    """Return gene-looking tokens from a cell, excluding variant notation."""
    if not value:
        return set()
    out: set[str] = set()
    for part in re.split(r"[\s,;/()]+", str(value)):
        token = part.strip().strip("[]{}:.'\"").upper()
        if not token or token in _GENE_SYMBOL_IGNORE:
            continue
        if (
            re.match(r"^[A-Z]\d", token)
            or _GENE_TOKEN_PROTEIN_LIKE_RE.match(token)
            or _normalize_cdna(token)
        ):
            continue
        if _GENE_SYMBOL_CELL_RE.match(token) and sum(ch.isalpha() for ch in token) >= 2:
            out.add(token)
    return out


def _column_values_look_like_genes(values: List[str]) -> bool:
    """True when a column's cells read as gene symbols rather than IDs/prose.

    Deliberately open-vocabulary: a rostered symbol is strong evidence, but a
    panel column naming genes GVF has never seen must still qualify, so an
    all-caps alphanumeric token carrying a letter and not looking like a short
    family/sample ID also counts.
    """
    seen = [v.strip() for v in values if v and v.strip()]
    if not seen:
        return False
    recognized = {
        str(gene).strip().upper()
        for gene in known_gene_aliases(include_query_aliases=False)
    }
    for value in seen:
        token = value.strip().strip("*: ").upper()
        if not token or len(token) > 15:
            continue
        if token in recognized:
            return True
        # Gene-shaped: >=3 chars, starts with a letter, alphanumeric, and not a
        # bare enumeration like F1/P040/S3.
        if (
            len(token) >= 3
            and token[0].isalpha()
            and token.replace("-", "").isalnum()
            and not re.fullmatch(r"[A-Z]{1,2}\d{1,3}", token)
            and any(c.isalpha() for c in token)
            and sum(c.isalpha() for c in token) >= 3
        ):
            return True
    return False


def _infer_unnamed_gene_column(
    header_cells: List[str],
    data_lines: List[str],
    *,
    target_gene: Optional[str] = None,
) -> Optional[int]:
    """Infer a blank-header gene column from table structure.

    Known symbols and aliases remain strong evidence, but are not required.
    Two distinct cells that consist solely of gene-looking tokens are enough to
    identify an open-vocabulary multi-gene grouping column.  Requiring a blank
    header plus repeated token-only cells avoids treating prose, phenotype, or
    family-ID columns as genes while allowing previously unseen panel genes.
    """
    recognized_genes = {
        str(gene).strip().upper()
        for gene in known_gene_aliases(include_query_aliases=True)
        if str(gene).strip()
    }
    if target_gene:
        recognized_genes.update(_target_gene_tokens(target_gene))

    for idx, header in enumerate(header_cells):
        if _normalize_header(header):
            continue

        observed_tokens: set[str] = set()
        token_only_values: set[str] = set()
        saw_blank_continuation = False
        for row in data_lines:
            cells = _split_pipe_row(row)
            if idx >= len(cells):
                continue
            value = cells[idx].strip()
            if not value:
                saw_blank_continuation = True
                continue
            tokens = _gene_symbol_tokens(value)
            observed_tokens.update(tokens)
            cleaned = re.sub(r"[*_`\\]", "", value).strip().upper()
            if len(tokens) == 1 and cleaned in tokens:
                token_only_values.update(tokens)

        if observed_tokens & recognized_genes or (
            saw_blank_continuation and len(token_only_values) >= 2
        ):
            return idx
    return None


def _row_has_off_target_gene_without_target(cells: List[str], gene_symbol: str) -> bool:
    """Detect rows with an explicit non-target gene token but no target token.

    This is deliberately narrower than "row must mention the target gene".
    Single-gene clinical tables often list only variant notations per row, so a
    hard positive target mention would destroy true positives. The guard only
    drops rows that visibly name another gene, covering misrouted panel tables
    whose gene column was not mapped by the router.
    """
    target_tokens = _target_gene_tokens(gene_symbol)
    if _cell_mentions_target_gene(" ".join(cells), gene_symbol):
        return False
    row_tokens: set[str] = set()
    for cell in cells:
        row_tokens.update(_gene_symbol_tokens(cell))
    if not row_tokens:
        return False
    if row_tokens & target_tokens:
        return False
    return bool(row_tokens - target_tokens)


def _sum_optional_int(a: Optional[int], b: Optional[int]) -> Optional[int]:
    if a is None and b is None:
        return None
    return (a or 0) + (b or 0)


def route_tables(
    text: str,
    gene_symbol: str,
    *,
    model: str,
    llm_caller: Optional[Any] = None,
    max_tokens: int = 8192,
    temperature: float = 0.0,
    reasoning_effort: Optional[str] = None,
    pmid: Optional[str] = None,
) -> RouterResult:
    """Run the full enumerate → LLM-route flow and return parsed routes.

    `llm_caller` is any callable matching the `litellm.completion` shape — it
    receives `model`, `messages`, `temperature`, and `max_tokens` kwargs and
    returns an object whose `.choices[0].message.content` is the model
    response. Plumbed as a parameter so tests can stub it out.

    Routing runs under its own trace stage and emits a
    ``table_routing_decision`` event. Previously the router's prompt/response
    record carried no stage or component at all, so it was indistinguishable
    from primary extraction, and the routing outcome existed only inside
    ``extraction_metadata`` — not as a normalized decision a curator could
    follow. ``pmid`` is accepted for standalone callers; inside extraction the
    enclosing paper scope already supplies gene and PMID.
    """
    from utils.llm_trace import (
        attempt_link_summary,
        llm_attempt_ledger,
        llm_trace_scope,
        record_trace_event,
    )

    with (
        llm_trace_scope(
            gene=gene_symbol,
            pmid=pmid,
            stage="clinical_table_routing",
            component="table_router",
        ),
        llm_attempt_ledger(),
    ):
        result = _route_tables_impl(
            text,
            gene_symbol,
            model=model,
            llm_caller=llm_caller,
            max_tokens=max_tokens,
            temperature=temperature,
            reasoning_effort=reasoning_effort,
        )
        routed_ids = [r.table_id for r in result.routed_tables]
        links = attempt_link_summary()
        if links["provider_attempts"] == 0:
            # Every table mapped deterministically (or there were none), so no
            # model was consulted. Say so instead of implying a model call.
            decision_source = "deterministic"
        elif result.error:
            decision_source = "router_call_failed"
        elif links["accepted_response_trace_id"] is None:
            decision_source = "router_response_unusable"
        else:
            decision_source = "model_routed"
        record_trace_event(
            "table_routing_decision",
            {
                "model": model if links["provider_attempts"] else None,
                "enumerated_table_ids": list(result.enumerated_table_ids or []),
                "routed_table_ids": routed_ids,
                "declined_table_ids": [
                    table_id
                    for table_id in (result.enumerated_table_ids or [])
                    if table_id not in set(routed_ids)
                ],
                "routed_table_count": len(routed_ids),
                "deterministic_route_count": len(routed_ids)
                - len(
                    [
                        r
                        for r in result.routed_tables
                        if r.notes != "deterministic header/data mapping"
                    ]
                ),
                "used_fallback": result.used_fallback,
                "error": result.error,
                "decision_source": decision_source,
                "selection_policy": (
                    "Deterministic header/data mapping first; only unmapped "
                    "tables are offered to the router model. A router failure "
                    "or unusable response falls back to the deterministic routes."
                ),
                **links,
            },
        )
        return result


def _route_tables_impl(
    text: str,
    gene_symbol: str,
    *,
    model: str,
    llm_caller: Optional[Any] = None,
    max_tokens: int = 8192,
    temperature: float = 0.0,
    reasoning_effort: Optional[str] = None,
) -> RouterResult:
    tables = enumerate_markdown_tables(text)
    if not tables:
        return RouterResult(routed_tables=[], used_fallback=True)

    enumerated_ids = [t.table_id for t in tables]
    deterministic: List[RoutedTable] = []
    llm_candidates: List[MarkdownTable] = []
    strict_cohort_labels = _strict_cohort_labels_enabled()
    for table in tables:
        mapping = _infer_column_mapping_from_headers(
            table,
            strict_cohort_labels=strict_cohort_labels,
            target_gene=gene_symbol,
        )
        if mapping:
            deterministic.append(
                RoutedTable(
                    table_id=table.table_id,
                    column_mapping=mapping,
                    confidence=1.0,
                    notes="deterministic header/data mapping",
                    table=table,
                )
            )
        else:
            llm_candidates.append(table)

    if deterministic and not llm_candidates:
        return RouterResult(
            routed_tables=deterministic,
            used_fallback=False,
            enumerated_table_ids=enumerated_ids,
        )

    if not llm_candidates:
        return RouterResult(
            routed_tables=deterministic,
            used_fallback=not deterministic,
            enumerated_table_ids=enumerated_ids,
        )

    if llm_caller is None:
        from utils.llm_utils import litellm_completion as llm_caller  # type: ignore

    prompt = build_router_prompt(llm_candidates, gene_symbol)

    try:
        from utils.llm_utils import (
            build_reasoning_effort_kwargs,
            wait_for_llm_rate_limit,
        )
        from utils.llm_trace import (
            OUTCOME_DISCARDED,
            OUTCOME_ERROR,
            OUTCOME_PARSE_FAILED,
            OUTCOME_PARSED,
            last_llm_trace,
            note_llm_attempt,
            note_llm_outcome,
        )
        from utils.retry_utils import llm_retry

        # The router calls litellm_completion DIRECTLY (it is not a
        # BaseLLMCaller), so nothing registered its calls in the ledger: its
        # decision event carried no accepted link and its empty-response retry
        # and parse failures were invisible. Register each call here.
        @llm_retry
        def _call(role: str) -> tuple[Any, Optional[str]]:
            wait_for_llm_rate_limit(model)
            try:
                response = llm_caller(
                    model=model,
                    messages=[
                        {
                            "role": "system",
                            "content": (
                                "You are a careful router for biomedical table "
                                "extraction. Output strict JSON. No prose."
                            ),
                        },
                        {"role": "user", "content": prompt},
                    ],
                    temperature=temperature,
                    max_tokens=max_tokens,
                    **build_reasoning_effort_kwargs(model, reasoning_effort),
                )
            except BaseException:
                note_llm_outcome(
                    note_llm_attempt(last_llm_trace(), role=role), OUTCOME_ERROR
                )
                raise
            return response, note_llm_attempt(last_llm_trace(), role=role)

        response, router_trace_id = _call("table_router")
        raw = response.choices[0].message.content or ""
        if not raw.strip():
            logger.warning(
                "table_router: empty LLM response from %s; retrying once", model
            )
            note_llm_outcome(router_trace_id, OUTCOME_DISCARDED)
            response, router_trace_id = _call("empty_content_retry")
            raw = response.choices[0].message.content or ""
    except Exception as e:  # noqa: BLE001
        logger.warning("table_router: LLM call failed: %s", e)
        return RouterResult(
            routed_tables=deterministic,
            error=str(e),
            used_fallback=not deterministic,
            enumerated_table_ids=enumerated_ids,
        )

    if not raw.strip():
        error = f"empty LLM response from table router model {model}"
        logger.warning("table_router: %s", error)
        note_llm_outcome(router_trace_id, OUTCOME_DISCARDED)
        return RouterResult(
            routed_tables=deterministic,
            error=error,
            used_fallback=not deterministic,
            enumerated_table_ids=enumerated_ids,
        )

    routed = parse_router_response(raw)
    if not routed:
        # Non-empty content that yielded no usable route is a parse failure, not
        # an accepted response — the deterministic routes are what get used.
        note_llm_outcome(router_trace_id, OUTCOME_PARSE_FAILED)
    else:
        note_llm_outcome(router_trace_id, OUTCOME_PARSED)
    # Attach the corresponding MarkdownTable for downstream parsing
    by_id = {t.table_id: t for t in llm_candidates}
    for r in routed:
        r.table = by_id.get(r.table_id)
    routed = [r for r in routed if r.table is not None]
    return RouterResult(
        routed_tables=deterministic + routed,
        raw_response=raw,
        enumerated_table_ids=enumerated_ids,
    )


def extract_via_router(
    text: str,
    gene_symbol: str,
    *,
    model: str,
    llm_caller: Optional[Any] = None,
    max_tokens: int = 8192,
    reasoning_effort: Optional[str] = None,
    pmid: Optional[str] = None,
) -> Dict[str, Any]:
    """End-to-end: route, then deterministically parse, returning a variant list.

    Returns a dict with keys:
      - ``routed`` (list[RoutedTable]): tables the router selected
      - ``variants`` (list[dict]): variant records produced by the parser
      - ``used_fallback`` (bool): True if no usable tables were found
      - ``error`` (str|None): set when the router LLM call failed
      - ``enumerated_table_ids`` (list[str]): every table id offered to the
        router, so callers can record what was declined, not just what was kept
    """
    result = route_tables(
        text,
        gene_symbol,
        model=model,
        llm_caller=llm_caller,
        max_tokens=max_tokens,
        reasoning_effort=reasoning_effort,
        pmid=pmid,
    )
    variants: List[Dict[str, Any]] = []
    for routed in result.routed_tables:
        if routed.table is None:
            continue
        variants.extend(
            parse_routed_table(routed.table, routed.column_mapping, gene_symbol)
        )
    return {
        "routed": result.routed_tables,
        "variants": variants,
        "used_fallback": result.used_fallback or not result.routed_tables,
        "error": result.error,
        "enumerated_table_ids": result.enumerated_table_ids,
    }


def _header_label(table: MarkdownTable, idx: Optional[int]) -> Optional[str]:
    if idx is None or idx < 0 or idx >= len(table.header_cells):
        return None
    label = table.header_cells[idx].strip()
    return label or None


def _usable_count_index(
    table: MarkdownTable,
    idx: Optional[int],
    normalized_headers: List[str],
) -> Optional[int]:
    """Drop routed count columns that are clearly IDs, denominators, or AFs."""
    if idx is None or idx == _INFER_ROW_PATIENT_COUNT:
        return idx
    if idx < 0 or idx >= len(normalized_headers):
        return None
    if _is_non_variant_count_header(normalized_headers[idx], normalized_headers):
        logger.debug(
            "table_router: ignoring non-carrier count column %r in %s",
            _header_label(table, idx),
            table.table_id,
        )
        return None
    return idx


def _router_count_type_for_label(label: Optional[str], role: str) -> Optional[str]:
    normalized = _normalize_header(label or "")
    if not normalized:
        return None
    if role == "unaffected" and any(
        token in normalized
        for token in (
            "control",
            "healthy",
            "normal",
            "unaffected",
            "asymptomatic",
        )
    ):
        return "unaffected_control"
    if role == "affected" and normalized in {
        "case",
        "cases",
        "brs",
        "lqt",
        "brslqt",
        "disease",
        "diseased",
    }:
        return "case"
    if "family" in normalized or "families" in normalized or "kindred" in normalized:
        return "family_count"
    if "proband" in normalized:
        return "proband_count"
    if any(
        token in normalized
        for token in ("screened", "ntested", "numbertested", "screenedn")
    ):
        return "screened_N"
    if any(
        token in normalized
        for token in (
            "totalcase",
            "totalcontrol",
            "totalsample",
            "samplesize",
            "cohortsize",
        )
    ):
        return "cohort_total"
    return "per_variant_carrier"


def _router_count_provenance(
    table: MarkdownTable,
    count_idx: Optional[int],
    aff_idx: Optional[int],
    unaff_idx: Optional[int],
) -> Dict[str, Optional[str]]:
    carrier_label = (
        "implicit one carrier per clinical row"
        if count_idx == _INFER_ROW_PATIENT_COUNT
        else _header_label(table, count_idx)
    )
    affected_label = _header_label(table, aff_idx)
    unaffected_label = _header_label(table, unaff_idx)
    return {
        "carriers_column_label": carrier_label,
        "carriers_count_type": _router_count_type_for_label(carrier_label, "carriers"),
        "affected_column_label": affected_label,
        "affected_count_type": _router_count_type_for_label(affected_label, "affected"),
        "unaffected_column_label": unaffected_label,
        "unaffected_count_type": _router_count_type_for_label(
            unaffected_label, "unaffected"
        ),
    }


def _router_fact_rows(
    *,
    table: MarkdownTable,
    row_number: int,
    row_text: str,
    cdna: Optional[str],
    protein: Optional[str],
    total: Optional[int],
    affected: Optional[int],
    unaffected: Optional[int],
    uncertain: Optional[int],
    count_idx: Optional[int],
    aff_idx: Optional[int],
    unaff_idx: Optional[int],
    unc_idx: Optional[int],
    source_location: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    source_table = table.caption or f"Table {table.table_id}"
    base = {
        "source_location": source_location,
        "source_table": source_table,
        "source_row": str(row_number),
        "evidence_quote": row_text.strip(),
        "source_layer": "llm_table",
    }
    for value in (protein, cdna):
        if value:
            rows.append(
                {
                    **base,
                    "fact_type": "variant_identity",
                    "fact_value": value,
                }
            )

    def add_count(
        fact_type: str, value: Optional[int], idx: Optional[int], count_type: str
    ) -> None:
        if value is None:
            return
        rows.append(
            {
                **base,
                "fact_type": fact_type,
                "fact_value": value,
                "source_column": (
                    "implicit one carrier per clinical row"
                    if idx == _INFER_ROW_PATIENT_COUNT
                    else _header_label(table, idx)
                ),
                "count_type": count_type,
            }
        )

    add_count(
        "patient_count",
        total,
        count_idx,
        _router_count_type_for_label(
            "implicit one carrier per clinical row"
            if count_idx == _INFER_ROW_PATIENT_COUNT
            else _header_label(table, count_idx),
            "carriers",
        )
        or "per_variant_carrier",
    )
    add_count(
        "total_carriers_observed",
        total,
        count_idx,
        _router_count_type_for_label(
            "implicit one carrier per clinical row"
            if count_idx == _INFER_ROW_PATIENT_COUNT
            else _header_label(table, count_idx),
            "carriers",
        )
        or "per_variant_carrier",
    )
    add_count(
        "affected_count",
        affected,
        aff_idx if aff_idx is not None else count_idx,
        _router_count_type_for_label(
            _header_label(table, aff_idx)
            or (
                "implicit one carrier per clinical row"
                if count_idx == _INFER_ROW_PATIENT_COUNT
                else _header_label(table, count_idx)
            ),
            "affected",
        )
        or "per_variant_carrier",
    )
    add_count(
        "unaffected_count",
        unaffected,
        unaff_idx if unaff_idx is not None else count_idx,
        _router_count_type_for_label(
            _header_label(table, unaff_idx)
            or (
                "implicit one carrier per clinical row"
                if count_idx == _INFER_ROW_PATIENT_COUNT
                else _header_label(table, count_idx)
            ),
            "unaffected",
        )
        or "per_variant_carrier",
    )
    add_count("uncertain_count", uncertain, unc_idx, "per_variant_carrier")
    return rows


def _router_observation_provenance(
    table: MarkdownTable,
    row_number: int,
    count_idx: Optional[int],
) -> Dict[str, Any]:
    source_ref = table.caption or f"Table {table.table_id}"
    source_container = (
        "supplement"
        if re.search(
            r"\bsupp(?:lement(?:ary|al)?)?\b|etable|e-table", source_ref, re.IGNORECASE
        )
        else "main"
    )
    column_ref = (
        "implicit one carrier per clinical row"
        if count_idx == _INFER_ROW_PATIENT_COUNT
        else _header_label(table, count_idx)
    )
    return {
        "source_container": source_container,
        "source_kind": "table",
        "source_ref": source_ref,
        "row_ordinal": row_number,
        "column_ref": column_ref,
        "locator_extra": {
            "parser": "table_router",
            "table_id": table.table_id,
        },
    }


def parse_routed_table(
    table: MarkdownTable, mapping: Dict[str, int], gene_symbol: str
) -> List[Dict[str, Any]]:
    """Run the deterministic cell parser over a router-approved table.

    The variant dict shape matches what `ExpertExtractor._parse_markdown_table_variants`
    produces, so downstream aggregation/SQLite paths don't need to change.

    When the router supplied a `gene` column, rows whose gene cell does not
    match `gene_symbol` (case-insensitive substring match against canonical
    gene symbol) are skipped — multi-gene panel tables therefore contribute
    only the relevant variants.
    """
    gene_idx = mapping.get("gene")
    cdna_idx = mapping.get("cdna")
    protein_idx = mapping.get("protein")
    normalized_headers = [_normalize_header(c) for c in table.header_cells]
    count_idx = _usable_count_index(
        table, mapping.get("patient_count"), normalized_headers
    )
    aff_idx = _usable_count_index(table, mapping.get("affected"), normalized_headers)
    unaff_idx = _usable_count_index(
        table, mapping.get("unaffected"), normalized_headers
    )
    unc_idx = _usable_count_index(table, mapping.get("uncertain"), normalized_headers)
    pheno_idx = mapping.get("phenotype")
    clin_idx = mapping.get("clinical_significance")

    if count_idx is None and aff_idx is None and unaff_idx is None:
        sanitized_mapping = {
            key: value
            for key, value in mapping.items()
            if key
            not in {
                "patient_count",
                "affected",
                "unaffected",
                "uncertain",
            }
        }
        has_assay_or_gwas_cue = any(
            h in _GWAS_OR_ASSAY_HEADERS for h in normalized_headers
        )
        if _looks_like_row_level_clinical_list(
            table, sanitized_mapping, has_assay_or_gwas_cue
        ):
            count_idx = _INFER_ROW_PATIENT_COUNT
        else:
            return []

    variants: List[Dict[str, Any]] = []
    by_key: Dict[tuple[str, str], Dict[str, Any]] = {}
    target_gene_lower = gene_symbol.strip().lower() if gene_symbol else ""
    caption_scope = _caption_gene_scope(table.caption)
    # Reject a table whose caption names a DIFFERENT gene — but only when that
    # absence is actually informative.
    #
    # `_caption_gene_scope` can only resolve genes in the built-in roster (14).
    # For an off-roster TARGET (TTN, LMNA, MYH7, PKP2, DSP — LMNA alone has
    # ~1,149 corpus papers) the target can never appear in any caption scope, so
    # this test degenerates into "reject every captioned table" and the gene
    # yields nothing at all. That breaks the standing requirement that GVF work
    # turnkey on genes with no gold standard.
    #
    # A mapped Gene column is also strictly better evidence than the caption: it
    # adjudicates per row, so a two-gene table keeps the target's rows instead of
    # losing all of them. Let the row filter below handle those.
    target_upper = gene_symbol.strip().upper()
    target_is_rostered = target_upper in known_gene_aliases(include_query_aliases=True)
    if (
        caption_scope
        and target_is_rostered
        and gene_idx is None
        and target_upper not in caption_scope
    ):
        return []
    # A caption that affirmatively scopes the table to the target gene ("Summary
    # of putative LQT1-associated mutations in KCNQ1") has already answered the
    # question the row-level off-target guard exists to answer, and answered it
    # from a more reliable signal. Running the row guard anyway is not merely
    # redundant, it is destructive: `_gene_symbol_tokens` accepts any 3-12
    # character uppercase token, so protein-domain shorthand (CNBD, DII, DI-S6,
    # PAC, SAR), mutation-type fragments (FRAME, SHIFT, DEL, INS, DUP, SITE) and
    # even raw nucleotide runs (CGGGGCGAC) all read as "some other gene" and
    # delete the row. That vocabulary is open-ended — the ignore-list cannot
    # converge on it — but it only ever appears in tables whose gene is already
    # established. PMID 26496715 lost 33 of 99 count-bearing gold rows this way
    # even after the common mutation-type words were excluded by name.
    #
    # The guard still runs on captions that name no gene, which is the misrouted
    # panel case it was written for, and a caption naming a DIFFERENT gene is
    # rejected outright above (the cross-gene defense pinned by
    # tests/unit/test_extraction_table_parser.py).
    #
    # UNAMBIGUOUS captions only. A caption naming SEVERAL genes ("Variant
    # Calling of Germline Mutations in BRCA1 and BRCA2 in FFPE Samples") does
    # NOT answer the question the row guard exists to answer — it says the
    # table is about both, so the rows still need separating. Treating that as
    # "scoped to the target" skipped the guard for the whole table and returned
    # an identical row set for either gene. This was latent while captions were
    # never captured; capturing them made it live, and it regressed
    # PMID 27767231 (a BRCA1 review-queue paper whose gene lives in a
    # "Germline mutation" column that `_infer_column_mapping_from_headers`
    # maps to `protein`, not `gene`, so there is no gene-column guard either).
    caption_scopes_to_target = (
        len(caption_scope) == 1 and gene_symbol.strip().upper() in caption_scope
    )
    # Markdown rowspan: a blank Gene cell inherits the gene from the row above
    # (gene-grouped tables list the gene once, then leave continuation rows
    # blank). Forward-fill so the gene-filter below sees the true gene of every
    # row; otherwise off-target continuation rows (e.g. HCN4 Val759Ile under a
    # KCNH2 extraction) leak through and are mis-stamped with the target gene.
    last_gene_cell = ""
    # In-band gene partition: set by a section-divider row, applies to every
    # row until the next divider. Independent of caption scope, because a
    # caption naming BOTH genes ("Polymorphisms in BRCA1 and BRCA2 genes")
    # legitimately passes the cross-gene reject yet still needs its rows split.
    section_gene: Optional[str] = None
    seen_divider = False

    for row_number, row in enumerate(table.data_lines, start=1):
        cells = _split_pipe_row(row)
        if not cells or len(cells) < 2:
            continue
        if all(not c.strip() for c in cells):
            continue

        divider_gene = _gene_section_divider(
            cells, target_gene=gene_symbol, allow_unknown=seen_divider
        )
        if divider_gene is not None:
            section_gene = divider_gene
            seen_divider = True
            continue
        if (
            section_gene is not None
            and gene_symbol
            and section_gene.upper() != gene_symbol.strip().upper()
        ):
            continue

        def cell(idx: Optional[int]) -> Optional[str]:
            if idx is None or idx < 0 or idx >= len(cells):
                return None
            return cells[idx]

        # Gene-filter: skip rows that explicitly belong to a different gene.
        if gene_idx is not None and target_gene_lower:
            gene_cell = (cell(gene_idx) or "").strip()
            if gene_cell:
                last_gene_cell = gene_cell
            else:
                gene_cell = last_gene_cell  # inherit markdown rowspan
            gene_tokens = _gene_symbol_tokens(gene_cell.upper())
            if (
                gene_cell
                and not _cell_mentions_target_gene(gene_cell, gene_symbol)
                and not (gene_tokens & _target_gene_tokens(gene_symbol))
            ):
                continue
        elif target_gene_lower and not caption_scopes_to_target:
            non_gene_context_indices = {
                idx
                for idx in (count_idx, aff_idx, unaff_idx, unc_idx, pheno_idx, clin_idx)
                if idx is not None and idx >= 0
            }
            cells_for_gene_guard = [
                value
                for idx, value in enumerate(cells)
                if idx not in non_gene_context_indices
            ]
            if _row_has_off_target_gene_without_target(
                cells_for_gene_guard, gene_symbol
            ):
                continue

        cdna_values = _normalize_cdna_values(cell(cdna_idx))
        protein_values = _normalize_protein_values(cell(protein_idx))

        if not cdna_values and not protein_values:
            continue

        infer_one_carrier = count_idx == _INFER_ROW_PATIENT_COUNT
        total = 1 if infer_one_carrier else _coerce_int(cell(count_idx))
        affected = _coerce_int(cell(aff_idx))
        unaffected = _coerce_int(cell(unaff_idx))
        uncertain = _coerce_int(cell(unc_idx))
        phenotype = cell(pheno_idx)
        clinical = cell(clin_idx)

        # A carrier/patient total does not establish either phenotype partition.
        # Only the mapped phenotype cell can populate one. Reading clinical
        # status out of the whole row let an unrelated functional-effect cell
        # ("normal trafficking", "control-like kinetics") record a symptomatic
        # proband as a confirmed unaffected carrier. `clinical` is deliberately
        # not consulted: `clinical_significance` is the variant's pathogenicity
        # call, not the carrier's phenotype.
        if total is not None and affected is None and unaffected is None:
            if infer_one_carrier and _phenotype_says_unaffected(phenotype):
                affected = 0
                unaffected = total
            elif infer_one_carrier and _phenotype_names_a_condition(phenotype):
                affected = total

        # If we have affected + unaffected but no total, derive it.
        if total is None and affected is not None and unaffected is not None:
            total = affected + unaffected + (uncertain or 0)
            if total == 0:
                total = None

        n_variants = max(len(cdna_values), len(protein_values), 1)
        for variant_idx in range(n_variants):
            cdna = (
                cdna_values[variant_idx]
                if variant_idx < len(cdna_values)
                else cdna_values[0]
                if len(cdna_values) == 1
                else None
            )
            protein = (
                protein_values[variant_idx]
                if variant_idx < len(protein_values)
                else protein_values[0]
                if len(protein_values) == 1
                else None
            )

            if not cdna and not protein:
                continue

            source_table = table.caption or f"Table {table.table_id}"
            source_location = f"{source_table}, row {row_number} (router+deterministic)"
            observation_provenance = _router_observation_provenance(
                table, row_number, count_idx
            )
            fact_rows = _router_fact_rows(
                table=table,
                row_number=row_number,
                row_text=row,
                cdna=cdna,
                protein=protein,
                total=total,
                affected=affected,
                unaffected=unaffected,
                uncertain=uncertain,
                count_idx=count_idx,
                aff_idx=aff_idx,
                unaff_idx=unaff_idx,
                unc_idx=unc_idx,
                source_location=source_location,
            )
            dedup_key = ((cdna or "").lower(), (protein or "").lower())
            if dedup_key in by_key:
                existing = by_key[dedup_key]
                existing.setdefault("fact_provenance", []).extend(fact_rows)
                existing["patients"]["count"] = _sum_optional_int(
                    existing["patients"].get("count"), total
                )
                pen = existing["penetrance_data"]
                pen["total_carriers_observed"] = _sum_optional_int(
                    pen.get("total_carriers_observed"), total
                )
                pen["affected_count"] = _sum_optional_int(
                    pen.get("affected_count"), affected
                )
                pen["unaffected_count"] = _sum_optional_int(
                    pen.get("unaffected_count"), unaffected
                )
                pen["uncertain_count"] = _sum_optional_int(
                    pen.get("uncertain_count"), uncertain
                )
                continue

            variant = {
                "gene_symbol": gene_symbol,
                "cdna_notation": cdna,
                "protein_notation": protein,
                "clinical_significance": (clinical or "pathogenic").strip().lower()
                if clinical
                else "pathogenic",
                "patients": {
                    "count": total,
                    "phenotype": phenotype.strip() if phenotype else None,
                    **observation_provenance,
                },
                "penetrance_data": {
                    "total_carriers_observed": total,
                    "affected_count": affected,
                    "unaffected_count": unaffected,
                    "uncertain_count": uncertain,
                    "penetrance_percentage": None,
                    "age_dependent_penetrance": [],
                },
                "individual_records": [],
                "functional_data": {"summary": "", "assays": []},
                "segregation_data": None,
                "population_frequency": None,
                "evidence_level": "medium",
                "source_location": source_location,
                **observation_provenance,
                "additional_notes": (
                    f"Parsed deterministically from {table.table_id}; "
                    f"caption={table.caption!r}"
                    if table.caption
                    else f"Parsed deterministically from {table.table_id}"
                ),
                "key_quotes": [],
                "count_provenance": _router_count_provenance(
                    table, count_idx, aff_idx, unaff_idx
                ),
                "fact_provenance": fact_rows,
            }
            by_key[dedup_key] = variant
            variants.append(variant)

    return variants
