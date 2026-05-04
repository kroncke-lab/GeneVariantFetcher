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
of {{patient_count, affected, unaffected}}. Multi-gene panel tables qualify
even if the preview rows show non-{gene_symbol} variants. Tables that list
functional assays, drug screens, primer sequences, or in silico predictions
do NOT qualify — skip them.

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
    # Find the first { and last }
    first = text.find("{")
    last = text.rfind("}")
    if first == -1 or last == -1 or last < first:
        return []
    text = text[first : last + 1]

    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        logger.warning("table_router: JSON parse failed: %s", e)
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
                clean[field_name] = int(idx)
            except (TypeError, ValueError):
                continue
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


def _normalize_cdna(value: str) -> Optional[str]:
    s = value.strip().replace(" ", "")
    if not s or s.lower() in {"-", "na", "n/a", "none", "."}:
        return None
    if not s.lower().startswith("c."):
        s = "c." + s
    return s


def _normalize_protein(value: str) -> Optional[str]:
    s = value.strip().replace(" ", "")
    if not s or s.lower() in {"-", "na", "n/a", "none", "."}:
        return None
    return s


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

    if llm_caller is None:
        from litellm import completion as llm_caller  # type: ignore

    prompt = build_router_prompt(tables, gene_symbol)

    try:
        from utils.retry_utils import llm_retry

        @llm_retry
        def _call() -> Any:
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
        return RouterResult(error=str(e), used_fallback=True)

    routed = parse_router_response(raw)
    # Attach the corresponding MarkdownTable for downstream parsing
    by_id = {t.table_id: t for t in tables}
    for r in routed:
        r.table = by_id.get(r.table_id)
    routed = [r for r in routed if r.table is not None]
    return RouterResult(routed_tables=routed, raw_response=raw)


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
    seen_keys: set[str] = set()
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

        cdna_raw = cell(cdna_idx)
        protein_raw = cell(protein_idx)
        cdna = _normalize_cdna(cdna_raw) if cdna_raw else None
        protein = _normalize_protein(protein_raw) if protein_raw else None

        if not cdna and not protein:
            continue

        # Dedup by the most stable key available
        dedup_key = (cdna or protein or "").lower()
        if not dedup_key or dedup_key in seen_keys:
            continue
        seen_keys.add(dedup_key)

        total = _coerce_int(cell(count_idx))
        affected = _coerce_int(cell(aff_idx))
        unaffected = _coerce_int(cell(unaff_idx))
        uncertain = _coerce_int(cell(unc_idx))

        # If only patient_count exists, treat all carriers as affected (matches
        # the assumption already used by _parse_markdown_table_variants).
        if total is not None and affected is None and unaffected is None:
            affected = total
            unaffected = 0

        # If we have affected + unaffected but no total, derive it.
        if total is None and (affected is not None or unaffected is not None):
            total = (affected or 0) + (unaffected or 0) + (uncertain or 0)
            if total == 0:
                total = None

        phenotype = cell(pheno_idx)
        clinical = cell(clin_idx)

        variants.append(
            {
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
        )

    return variants
