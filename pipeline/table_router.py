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
    r"^\s*(table\s+\d+[a-z]?[.:]|tbl\.?\s*\d+[.:])", re.IGNORECASE
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

_ROW_LEVEL_SUBJECT_KEYWORDS = (
    "patient",
    "proband",
    "subject",
    "caseid",
    "case",
    "family",
    "kindred",
    "individual",
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
    return subject_cue or clinical_cue


def _infer_column_mapping_from_headers(
    table: MarkdownTable,
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

        if "gene" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["gene"]
        ):
            mapping["gene"] = idx

        if "unaffected" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["unaffected"]
        ):
            mapping["unaffected"] = idx
            continue

        if "affected" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["affected"]
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

        if "patient_count" not in mapping and any(
            kw in header for kw in _HEADER_FIELD_KEYWORDS["patient_count"]
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
  - patient_count     — total carriers / "n detected" / "number of times" / "patients"
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
    for table in tables:
        mapping = _infer_column_mapping_from_headers(table)
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
        from utils.llm_utils import wait_for_llm_rate_limit
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
) -> Dict[str, Any]:
    """End-to-end: route, then deterministically parse, returning a variant list.

    Returns a dict with keys:
      - ``routed`` (list[RoutedTable]): tables the router selected
      - ``variants`` (list[dict]): variant records produced by the parser
      - ``used_fallback`` (bool): True if no usable tables were found
      - ``error`` (str|None): set when the router LLM call failed
    """
    result = route_tables(
        text, gene_symbol, model=model, llm_caller=llm_caller, max_tokens=max_tokens
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
    count_idx = mapping.get("patient_count")
    aff_idx = mapping.get("affected")
    unaff_idx = mapping.get("unaffected")
    unc_idx = mapping.get("uncertain")
    pheno_idx = mapping.get("phenotype")
    clin_idx = mapping.get("clinical_significance")

    variants: List[Dict[str, Any]] = []
    by_key: Dict[tuple[str, str], Dict[str, Any]] = {}
    target_gene_lower = gene_symbol.strip().lower() if gene_symbol else ""

    for row in table.data_lines:
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
            gene_cell = (cell(gene_idx) or "").strip().lower()
            if gene_cell and target_gene_lower not in gene_cell:
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

            dedup_key = ((cdna or "").lower(), (protein or "").lower())
            if dedup_key in by_key:
                existing = by_key[dedup_key]
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
                "source_location": f"Table {table.table_id} (router+deterministic)",
                "additional_notes": (
                    f"Parsed deterministically from {table.table_id}; "
                    f"caption={table.caption!r}"
                    if table.caption
                    else f"Parsed deterministically from {table.table_id}"
                ),
                "key_quotes": [],
            }
            by_key[dedup_key] = variant
            variants.append(variant)

    return variants
