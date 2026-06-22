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
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

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


_TABLE_CAPTION_RE = re.compile(
    r"^\s*(e?table\s+\d+[a-z]?[.:]|tbl\.?\s*\d+[.:])", re.IGNORECASE
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

_UNAFFECTED_TEXT_RE = re.compile(
    r"\b(unaffected|asymptomatic|control|healthy|normal|no symptoms?)\b",
    re.IGNORECASE,
)

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
    return (
        _is_row_identifier_header(header)
        or _is_population_frequency_header(header)
        or _is_denominator_header_when_carrier_present(header, normalized_headers)
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

            # Look back for a caption — up to 3 lines above the header
            caption: Optional[str] = None
            for back in range(2, 5):
                idx = i - back
                if idx < 0:
                    break
                candidate = lines[idx].strip()
                if not candidate:
                    continue
                if candidate.startswith("|"):
                    break
                if _TABLE_CAPTION_RE.match(candidate):
                    caption = candidate
                break

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
        i += 1

    if only_variant_like:
        tables = [t for t in tables if _looks_like_variant_table(t)]
    return tables


def build_router_prompt(
    tables: List[MarkdownTable],
    gene_symbol: str,
    *,
    sample_rows: int = 3,
    cell_chars: int = 60,
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
}


def _target_gene_tokens(gene_symbol: str) -> set[str]:
    gene = (gene_symbol or "").strip().upper()
    aliases = get_gene_aliases(gene, include_query_aliases=True)
    return {alias.upper() for alias in aliases if alias} | {gene}


def _caption_gene_scope(caption: Optional[str]) -> set[str]:
    """Return target genes explicitly implied by a table caption."""
    if not caption:
        return set()
    scope: set[str] = set()
    for gene in known_gene_aliases(include_query_aliases=True):
        if gene_alias_regex(gene, include_query_aliases=True).search(caption):
            scope.add(gene)
    return scope


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
) -> RouterResult:
    """Run the full enumerate → LLM-route flow and return parsed routes.

    `llm_caller` is any callable matching the `litellm.completion` shape — it
    receives `model`, `messages`, `temperature`, and `max_tokens` kwargs and
    returns an object whose `.choices[0].message.content` is the model
    response. Plumbed as a parameter so tests can stub it out.
    """
    tables = enumerate_markdown_tables(text)
    if not tables:
        return RouterResult(routed_tables=[], used_fallback=True)

    deterministic: List[RoutedTable] = []
    llm_candidates: List[MarkdownTable] = []
    strict_cohort_labels = _strict_cohort_labels_enabled()
    for table in tables:
        mapping = _infer_column_mapping_from_headers(
            table, strict_cohort_labels=strict_cohort_labels
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
        return RouterResult(routed_tables=deterministic, used_fallback=False)

    if not llm_candidates:
        return RouterResult(
            routed_tables=deterministic, used_fallback=not deterministic
        )

    if llm_caller is None:
        from litellm import completion as llm_caller  # type: ignore

    prompt = build_router_prompt(llm_candidates, gene_symbol)

    try:
        from utils.llm_utils import (
            build_reasoning_effort_kwargs,
            wait_for_llm_rate_limit,
        )
        from utils.retry_utils import llm_retry

        @llm_retry
        def _call() -> Any:
            wait_for_llm_rate_limit(model)
            return llm_caller(
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

        response = _call()
        raw = response.choices[0].message.content or ""
    except Exception as e:  # noqa: BLE001
        logger.warning("table_router: LLM call failed: %s", e)
        return RouterResult(
            routed_tables=deterministic,
            error=str(e),
            used_fallback=not deterministic,
        )

    routed = parse_router_response(raw)
    # Attach the corresponding MarkdownTable for downstream parsing
    by_id = {t.table_id: t for t in llm_candidates}
    for r in routed:
        r.table = by_id.get(r.table_id)
    routed = [r for r in routed if r.table is not None]
    return RouterResult(routed_tables=deterministic + routed, raw_response=raw)


def extract_via_router(
    text: str,
    gene_symbol: str,
    *,
    model: str,
    llm_caller: Optional[Any] = None,
    max_tokens: int = 8192,
    reasoning_effort: Optional[str] = None,
) -> Dict[str, Any]:
    """End-to-end: route, then deterministically parse, returning a variant list.

    Returns a dict with keys:
      - ``routed`` (list[RoutedTable]): tables the router selected
      - ``variants`` (list[dict]): variant records produced by the parser
      - ``used_fallback`` (bool): True if no usable tables were found
      - ``error`` (str|None): set when the router LLM call failed
    """
    result = route_tables(
        text,
        gene_symbol,
        model=model,
        llm_caller=llm_caller,
        max_tokens=max_tokens,
        reasoning_effort=reasoning_effort,
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
    affected_label = _header_label(table, aff_idx) or carrier_label
    unaffected_label = _header_label(table, unaff_idx)
    return {
        "carriers_column_label": carrier_label,
        "carriers_count_type": "per_variant_carrier" if carrier_label else None,
        "affected_column_label": affected_label,
        "affected_count_type": "per_variant_carrier" if affected_label else None,
        "unaffected_column_label": unaffected_label,
        "unaffected_count_type": "per_variant_carrier" if unaffected_label else None,
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
        "per_variant_carrier",
    )
    add_count(
        "total_carriers_observed",
        total,
        count_idx,
        "per_variant_carrier",
    )
    add_count(
        "affected_count",
        affected,
        aff_idx if aff_idx is not None else count_idx,
        "per_variant_carrier",
    )
    add_count(
        "unaffected_count",
        unaffected,
        unaff_idx if unaff_idx is not None else count_idx,
        "per_variant_carrier",
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
    if caption_scope and gene_symbol.strip().upper() not in caption_scope:
        return []
    # Markdown rowspan: a blank Gene cell inherits the gene from the row above
    # (gene-grouped tables list the gene once, then leave continuation rows
    # blank). Forward-fill so the gene-filter below sees the true gene of every
    # row; otherwise off-target continuation rows (e.g. HCN4 Val759Ile under a
    # KCNH2 extraction) leak through and are mis-stamped with the target gene.
    last_gene_cell = ""

    for row_number, row in enumerate(table.data_lines, start=1):
        cells = _split_pipe_row(row)
        if not cells or len(cells) < 2:
            continue
        if all(not c.strip() for c in cells):
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
        elif target_gene_lower:
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

        # If only patient_count exists, treat all carriers as affected (matches
        # the assumption already used by _parse_markdown_table_variants).
        if total is not None and affected is None and unaffected is None:
            row_text = " ".join(cells)
            if infer_one_carrier and _UNAFFECTED_TEXT_RE.search(row_text):
                affected = 0
                unaffected = total
            else:
                affected = total
                unaffected = 0

        # If we have affected + unaffected but no total, derive it.
        if total is None and (affected is not None or unaffected is not None):
            total = (affected or 0) + (unaffected or 0) + (uncertain or 0)
            if total == 0:
                total = None

        phenotype = cell(pheno_idx)
        clinical = cell(clin_idx)

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
