"""A small, source-grounded GPT-5.6 Sol extraction protocol.

This module is deliberately independent from the production extraction cascade.
One reasoning model receives the complete cleaned paper directly when it is a
safe size, plus bounded read-only evidence tools.  Oversized papers stay behind
the same tools.  A deterministic parsing tool lets the model ask the host to
parse a table *after* it has identified human, per-variant count semantics.

The model never gets to make an unsupported fact authoritative.  Direct facts
must point to an exact source quote, accepted deterministic table batches are
revalidated, and duplicate restatements are merged without summing conflicting
counts.  The result uses the existing GVF extraction JSON shape so the normal
SQLite migration and scorer can consume it unchanged.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import asdict, dataclass, is_dataclass, replace
from typing import Any, Callable, Dict, List, Mapping, Optional, Protocol, Sequence

from pipeline.table_router import (
    MarkdownTable,
    enumerate_markdown_tables,
    parse_routed_table,
)


SOL_PROTOCOL_VERSION = "sol-simple-v2-full-context"
SOL_CLEANER_VERSION = "markdown-references-v1"
FULL_CONTEXT_MAX_CHARS = 500_000

_REFERENCE_HEADING_RE = re.compile(
    r"^(?:references?|bibliograph(?:y|ies)|works\s+cited|literature\s+cited)"
    r"(?:\s+and\s+notes?)?\s*$",
    re.IGNORECASE,
)
_SUPPLEMENT_HEADING_RE = re.compile(
    r"^(?:supplement(?:ary|al)?|supporting)\b", re.IGNORECASE
)
_HEADING_RE = re.compile(r"^(#{1,6})[ \t]+(.+?)[ \t]*#*[ \t]*$")
_SOURCE_LOCATOR_RE = re.compile(r"^S:(\d+):(\d+)$")
_TABLE_ROW_LOCATOR_RE = re.compile(r"^(T\d+):R(\d+)$", re.IGNORECASE)
_INTEGER_RE = re.compile(r"(?<!\d)(\d{1,3}(?:,\d{3})*|\d+)(?!\d)")

_AA3_TO_1 = {
    "ala": "A",
    "arg": "R",
    "asn": "N",
    "asp": "D",
    "cys": "C",
    "gln": "Q",
    "glu": "E",
    "gly": "G",
    "his": "H",
    "ile": "I",
    "leu": "L",
    "lys": "K",
    "met": "M",
    "phe": "F",
    "pro": "P",
    "ser": "S",
    "thr": "T",
    "trp": "W",
    "tyr": "Y",
    "val": "V",
    "ter": "*",
    "stop": "*",
}

_COUNT_KEYS = (
    "total_carriers_observed",
    "affected_count",
    "unaffected_count",
    "uncertain_count",
)
_DIRECT_TO_GVF_COUNT = {
    "total_carriers": "total_carriers_observed",
    "affected": "affected_count",
    "unaffected": "unaffected_count",
    "uncertain": "uncertain_count",
}
_DIRECT_COUNT_ROLES = {
    "total_carriers": "human_variant_carriers",
    "affected": "human_variant_affected",
    "unaffected": "human_variant_unaffected",
    "uncertain": "human_variant_uncertain",
}
_COUNT_MAPPING_ROLES = {
    "patient_count": {"human_variant_carriers", "human_variant_affected"},
    "affected": {"human_variant_affected"},
    "unaffected": {"human_variant_unaffected"},
    "uncertain": {"human_variant_uncertain"},
}
_ALLOWED_COLUMN_FIELDS = {
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


class ResponsesRunner(Protocol):
    """Narrow interface implemented by :class:`ResponsesAPIClient`."""

    def run_tool_loop(self, **kwargs: Any) -> Any:
        """Run one Responses API tool loop and return a ResponsesResult."""


def strip_irrelevant_context(text: str) -> str:
    """Remove Markdown reference sections without deleting later supplements.

    A reference heading starts a removable region.  The region ends at the next
    heading of the same or a higher level.  A supplement/supporting-information
    heading always ends it, even when a converter nested that heading one level
    below ``References``.  Tables and all non-reference prose are byte-for-byte
    preserved apart from the removed region and surrounding blank-line trim.
    """

    if not text:
        return ""
    lines = text.splitlines(keepends=True)
    kept: List[str] = []
    removing_level: Optional[int] = None

    for line in lines:
        heading = _HEADING_RE.match(line.rstrip("\r\n"))
        if heading:
            level = len(heading.group(1))
            heading_text = heading.group(2).strip()
            if removing_level is not None:
                if level <= removing_level or _SUPPLEMENT_HEADING_RE.match(
                    heading_text
                ):
                    removing_level = None
                else:
                    continue
            if _REFERENCE_HEADING_RE.fullmatch(heading_text):
                removing_level = level
                continue
        if removing_level is None:
            kept.append(line)

    return "".join(kept).strip()


# Short alias for callers that treat this as the protocol's cleaner.
clean_context = strip_irrelevant_context


@dataclass(frozen=True)
class _RowEvidence:
    table_id: str
    row_number: int
    text: str
    start: int
    end: int

    @property
    def locator(self) -> str:
        return f"{self.table_id}:R{self.row_number}"

    @property
    def source_locator(self) -> str:
        return f"S:{self.start}:{self.end}"


class EvidenceIndex:
    """Bounded, stable read-only access to one cleaned Markdown paper."""

    MAX_TABLE_LIST = 20
    MAX_TABLE_ROWS = 100
    MAX_ROW_CHARS = 2_000
    MAX_SOURCE_SPAN = 6_000
    MAX_SEARCH_HITS = 12
    SEARCH_CONTEXT_CHARS = 280

    def __init__(self, source_text: str) -> None:
        self.source_text = source_text or ""
        self.source_hash = hashlib.sha256(self.source_text.encode("utf-8")).hexdigest()
        self.tables: List[MarkdownTable] = enumerate_markdown_tables(
            self.source_text, only_variant_like=False
        )
        self._tables = {table.table_id.upper(): table for table in self.tables}
        self._rows: Dict[tuple[str, int], _RowEvidence] = {}
        self._index_rows()

    def _index_rows(self) -> None:
        for table in self.tables:
            cursor = max(0, table.char_start)
            table_end = min(len(self.source_text), max(cursor, table.char_end))
            for row_number, row in enumerate(table.data_lines, start=1):
                start = self.source_text.find(row, cursor, table_end)
                if start < 0:
                    # Converter offsets can be off by one at EOF; keep a stable
                    # table-row locator even if the auxiliary char span is not
                    # recoverable from the enumerator's range.
                    start = self.source_text.find(row, cursor)
                if start < 0:
                    start = cursor
                end = min(len(self.source_text), start + len(row))
                self._rows[(table.table_id.upper(), row_number)] = _RowEvidence(
                    table_id=table.table_id,
                    row_number=row_number,
                    text=row,
                    start=start,
                    end=end,
                )
                cursor = end

    def get_table(self, table_id: str) -> MarkdownTable:
        table = self._tables.get(str(table_id).strip().upper())
        if table is None:
            raise ValueError(f"unknown table_id: {table_id!r}")
        return table

    def list_tables(self, *, offset: int = 0, limit: int = 20) -> Dict[str, Any]:
        offset = _bounded_int(offset, "offset", minimum=0)
        limit = min(
            self.MAX_TABLE_LIST,
            _bounded_int(limit, "limit", minimum=1),
        )
        selected = self.tables[offset : offset + limit]
        return {
            "total_tables": len(self.tables),
            "offset": offset,
            "returned": len(selected),
            "has_more": offset + len(selected) < len(self.tables),
            "tables": [self._table_summary(table) for table in selected],
        }

    def _table_summary(self, table: MarkdownTable) -> Dict[str, Any]:
        return {
            "table_id": table.table_id,
            "caption": _trim_text(table.caption, 300),
            "header": _trim_text(table.header_preview(max_chars=500), 500),
            "row_count": len(table.data_lines),
            "source_locator": (
                f"S:{max(0, table.char_start)}:"
                f"{min(len(self.source_text), table.char_end)}"
            ),
        }

    def read_table(
        self, *, table_id: str, row_start: int = 1, row_end: int = 20
    ) -> Dict[str, Any]:
        table = self.get_table(table_id)
        row_start = _bounded_int(row_start, "row_start", minimum=1)
        row_end = _bounded_int(row_end, "row_end", minimum=row_start)
        if row_start > len(table.data_lines):
            raise ValueError(
                f"row_start {row_start} exceeds {table.table_id} row count "
                f"{len(table.data_lines)}"
            )
        bounded_end = min(
            len(table.data_lines), row_end, row_start + self.MAX_TABLE_ROWS - 1
        )
        rows: List[Dict[str, Any]] = []
        for row_number in range(row_start, bounded_end + 1):
            row = self._rows[(table.table_id.upper(), row_number)]
            rows.append(
                {
                    "locator": row.locator,
                    "source_locator": row.source_locator,
                    "text": _trim_text(row.text, self.MAX_ROW_CHARS),
                    "truncated": len(row.text) > self.MAX_ROW_CHARS,
                }
            )
        return {
            "table_id": table.table_id,
            "caption": _trim_text(table.caption, 500),
            "header": _trim_text(table.header_line, self.MAX_ROW_CHARS),
            "total_rows": len(table.data_lines),
            "row_start": row_start,
            "row_end": bounded_end,
            "requested_row_end": row_end,
            "truncated": bounded_end < min(row_end, len(table.data_lines)),
            "batch_locator": f"{table.table_id}:R{row_start}-R{bounded_end}",
            "rows": rows,
        }

    def read_source_span(self, *, start: int, end: int) -> Dict[str, Any]:
        start = _bounded_int(start, "start", minimum=0)
        end = _bounded_int(end, "end", minimum=start + 1)
        if start >= len(self.source_text):
            raise ValueError(
                f"start {start} exceeds source length {len(self.source_text)}"
            )
        requested_end = min(end, len(self.source_text))
        bounded_end = min(requested_end, start + self.MAX_SOURCE_SPAN)
        return {
            "locator": f"S:{start}:{bounded_end}",
            "requested_end": end,
            "source_length": len(self.source_text),
            "truncated": bounded_end < requested_end,
            "text": self.source_text[start:bounded_end],
        }

    def search_source(self, *, query: str, max_hits: int = 8) -> Dict[str, Any]:
        query = str(query or "").strip()
        if len(query) < 2:
            raise ValueError("query must contain at least two characters")
        if len(query) > 200:
            raise ValueError("query must be at most 200 characters")
        max_hits = min(
            self.MAX_SEARCH_HITS,
            _bounded_int(max_hits, "max_hits", minimum=1),
        )
        matches = list(re.finditer(re.escape(query), self.source_text, re.IGNORECASE))
        hits: List[Dict[str, Any]] = []
        for match in matches[:max_hits]:
            start = max(0, match.start() - self.SEARCH_CONTEXT_CHARS)
            end = min(len(self.source_text), match.end() + self.SEARCH_CONTEXT_CHARS)
            hits.append(
                {
                    "match_locator": f"S:{match.start()}:{match.end()}",
                    "context_locator": f"S:{start}:{end}",
                    "text": self.source_text[start:end],
                }
            )
        return {
            "query": query,
            "total_hits": len(matches),
            "returned": len(hits),
            "truncated": len(matches) > len(hits),
            "hits": hits,
        }

    def resolve_locator(self, locator: str) -> str:
        source_match = _SOURCE_LOCATOR_RE.fullmatch(str(locator).strip())
        if source_match:
            start, end = map(int, source_match.groups())
            if start < 0 or end <= start or end > len(self.source_text):
                raise ValueError(f"invalid source locator bounds: {locator!r}")
            if end - start > self.MAX_SOURCE_SPAN:
                raise ValueError(
                    f"source locator exceeds {self.MAX_SOURCE_SPAN} characters"
                )
            return self.source_text[start:end]

        row_match = _TABLE_ROW_LOCATOR_RE.fullmatch(str(locator).strip())
        if row_match:
            table_id, row_number_text = row_match.groups()
            row = self._rows.get((table_id.upper(), int(row_number_text)))
            if row is None:
                raise ValueError(f"unknown table row locator: {locator!r}")
            return row.text
        raise ValueError(
            "locator must be an exact source span S:start:end or table row Tn:Rn"
        )

    def row_evidence(self, table_id: str, row_number: int) -> _RowEvidence:
        row = self._rows.get((str(table_id).upper(), int(row_number)))
        if row is None:
            raise ValueError(f"unknown table row: {table_id}:R{row_number}")
        return row


@dataclass
class _ParsedBatch:
    batch_id: str
    table_id: str
    row_start: int
    row_end: int
    variants: List[Dict[str, Any]]
    request: Dict[str, Any]


class _ParsedBatchStore:
    """Host-side cache for table parses that Sol can accept by short ID."""

    MAX_PREVIEW_VARIANTS = 8

    def __init__(self, index: EvidenceIndex, gene_symbol: str) -> None:
        self.index = index
        self.gene_symbol = gene_symbol
        self.batches: Dict[str, _ParsedBatch] = {}

    def parse(self, args: Mapping[str, Any]) -> Dict[str, Any]:
        table_id = str(args.get("table_id") or "").strip().upper()
        table = self.index.get_table(table_id)
        row_start = _bounded_int(args.get("row_start"), "row_start", minimum=1)
        row_end = _bounded_int(args.get("row_end"), "row_end", minimum=row_start)
        if row_start > len(table.data_lines):
            raise ValueError(
                f"row_start {row_start} exceeds {table.table_id} row count "
                f"{len(table.data_lines)}"
            )
        bounded_end = min(
            row_end,
            len(table.data_lines),
            row_start + EvidenceIndex.MAX_TABLE_ROWS - 1,
        )

        mapping = _validate_column_mapping(args.get("column_mapping"), table)
        human_subjects = args.get("human_subjects")
        per_variant = args.get("per_variant")
        if human_subjects is not True or per_variant is not True:
            raise ValueError(
                "parse_variant_table requires human_subjects=true and per_variant=true"
            )
        semantics = args.get("count_semantics")
        if not isinstance(semantics, dict):
            raise ValueError("count_semantics must be an object")
        normalized_semantics = _validate_count_semantics(mapping, semantics)

        request = {
            "table_id": table.table_id,
            "row_start": row_start,
            "row_end": bounded_end,
            "column_mapping": mapping,
            "human_subjects": True,
            "per_variant": True,
            "count_semantics": normalized_semantics,
        }
        fingerprint = json.dumps(request, sort_keys=True, separators=(",", ":"))
        batch_id = (
            "PB-"
            + hashlib.sha256(
                f"{self.index.source_hash}:{fingerprint}".encode("utf-8")
            ).hexdigest()[:14]
        )
        if batch_id in self.batches:
            return self._response(self.batches[batch_id], requested_end=row_end)

        variants: List[Dict[str, Any]] = []
        # Parse each row separately.  parse_routed_table normally sums duplicate
        # identities inside one table; row-wise parsing preserves observations
        # so our conflict-aware deduper, not an implicit sum, makes the decision.
        last_gene_cell = self._last_gene_cell_before(
            table, mapping.get("gene"), row_start
        )
        for row_number in range(row_start, bounded_end + 1):
            raw_row = table.data_lines[row_number - 1]
            gene_idx = mapping.get("gene")
            if gene_idx is not None:
                cells = _split_table_row(raw_row)
                gene_cell = cells[gene_idx].strip() if gene_idx < len(cells) else ""
                if gene_cell:
                    last_gene_cell = gene_cell
                else:
                    gene_cell = last_gene_cell
                if not _cell_has_exact_gene(gene_cell, self.gene_symbol):
                    continue

            # Avoid table_router's corpus-wide gene-alias lookup here.  Sol has
            # already mapped the gene column and the host just applied an exact
            # target-cell check.  A caption-less, gene-less one-row table lets
            # the existing deterministic cell parser do notation/count parsing
            # without scanning the local VariantFeatures metadata database.
            parser_mapping = {
                key: value for key, value in mapping.items() if key != "gene"
            }
            one_row_table = replace(
                table,
                caption=None,
                data_lines=[raw_row],
            )
            parsed = parse_routed_table(one_row_table, parser_mapping, "")
            for variant in parsed:
                variant["gene_symbol"] = self.gene_symbol
                self._restore_row_locator(
                    variant,
                    table=table,
                    row_number=row_number,
                    semantics=normalized_semantics,
                )
                variants.append(variant)

        batch = _ParsedBatch(
            batch_id=batch_id,
            table_id=table.table_id,
            row_start=row_start,
            row_end=bounded_end,
            variants=variants,
            request=request,
        )
        self.batches[batch_id] = batch
        return self._response(batch, requested_end=row_end)

    @staticmethod
    def _last_gene_cell_before(
        table: MarkdownTable, gene_idx: Optional[int], row_start: int
    ) -> str:
        if gene_idx is None or row_start <= 1:
            return ""
        last = ""
        for raw_row in table.data_lines[: row_start - 1]:
            cells = _split_table_row(raw_row)
            if gene_idx < len(cells) and cells[gene_idx].strip():
                last = cells[gene_idx].strip()
        return last

    def _restore_row_locator(
        self,
        variant: Dict[str, Any],
        *,
        table: MarkdownTable,
        row_number: int,
        semantics: Mapping[str, Any],
    ) -> None:
        row = self.index.row_evidence(table.table_id, row_number)
        locator = row.locator
        variant["source_location"] = locator
        variant["source_table"] = table.caption or table.table_id
        variant["source_row"] = str(row_number)
        variant["evidence_quote"] = row.text
        variant["key_quotes"] = [row.text]
        variant["source_layer"] = "llm_table"
        variant["row_ordinal"] = row_number
        variant.setdefault("locator_extra", {})["stable_locator"] = locator

        patients = variant.setdefault("patients", {})
        patients["source_location"] = locator
        patients["row_ordinal"] = row_number
        patients.setdefault("locator_extra", {})["stable_locator"] = locator

        for fact in variant.get("fact_provenance", []):
            if not isinstance(fact, dict):
                continue
            fact["source_location"] = locator
            fact["source_row"] = str(row_number)
            fact["evidence_quote"] = row.text

        # parse_routed_table treats a lone patient-count column as affected.
        # Keep that only when Sol explicitly identified the column as affected;
        # a generic human-carrier column has unknown phenotype split.
        if (
            "patient_count" in semantics
            and semantics["patient_count"] == "human_variant_carriers"
            and semantics.get("affected") is None
            and semantics.get("unaffected") is None
        ):
            penetrance = variant.setdefault("penetrance_data", {})
            penetrance["affected_count"] = None
            penetrance["unaffected_count"] = None
            for fact in list(variant.get("fact_provenance", [])):
                if fact.get("fact_type") in {"affected_count", "unaffected_count"}:
                    variant["fact_provenance"].remove(fact)

    def _response(self, batch: _ParsedBatch, *, requested_end: int) -> Dict[str, Any]:
        preview: List[Dict[str, Any]] = []
        for variant in batch.variants[: self.MAX_PREVIEW_VARIANTS]:
            penetrance = variant.get("penetrance_data") or {}
            preview.append(
                {
                    "gene_symbol": variant.get("gene_symbol"),
                    "cdna_notation": variant.get("cdna_notation"),
                    "protein_notation": variant.get("protein_notation"),
                    "total_carriers": penetrance.get("total_carriers_observed"),
                    "affected": penetrance.get("affected_count"),
                    "unaffected": penetrance.get("unaffected_count"),
                    "locator": variant.get("source_location"),
                }
            )
        return {
            "batch_id": batch.batch_id,
            "table_id": batch.table_id,
            "batch_locator": f"{batch.table_id}:R{batch.row_start}-R{batch.row_end}",
            "requested_row_end": requested_end,
            "truncated": batch.row_end
            < min(
                requested_end,
                len(self.index.get_table(batch.table_id).data_lines),
            ),
            "candidate_count": len(batch.variants),
            "preview": preview,
            "preview_omitted": max(0, len(batch.variants) - self.MAX_PREVIEW_VARIANTS),
            "instruction": (
                "If these semantics are correct, include only this batch_id in "
                "accepted_table_batches; do not echo its rows as evidence_records."
            ),
        }


def _split_table_row(line: str) -> List[str]:
    parts = [cell.strip() for cell in str(line).split("|")]
    if parts and not parts[0]:
        parts = parts[1:]
    if parts and not parts[-1]:
        parts = parts[:-1]
    return parts


def _cell_has_exact_gene(value: str, gene_symbol: str) -> bool:
    tokens = {
        token.upper() for token in re.findall(r"[A-Za-z][A-Za-z0-9-]*", value or "")
    }
    return gene_symbol.strip().upper() in tokens


def _validate_column_mapping(raw_mapping: Any, table: MarkdownTable) -> Dict[str, int]:
    if not isinstance(raw_mapping, dict):
        raise ValueError("column_mapping must be an object")
    mapping: Dict[str, int] = {}
    for raw_key, raw_value in raw_mapping.items():
        key = str(raw_key).strip()
        if key not in _ALLOWED_COLUMN_FIELDS:
            raise ValueError(f"unsupported column mapping field: {key!r}")
        if raw_value is None:
            continue
        if isinstance(raw_value, bool) or not isinstance(raw_value, int):
            raise ValueError(f"column_mapping[{key!r}] must be an integer")
        if raw_value < 0 or raw_value >= len(table.header_cells):
            raise ValueError(
                f"column_mapping[{key!r}]={raw_value} is outside the "
                f"{len(table.header_cells)}-column table"
            )
        mapping[key] = raw_value
    if not any(key in mapping for key in ("cdna", "protein")):
        raise ValueError("column_mapping needs a cdna or protein column")
    if not any(key in mapping for key in _COUNT_MAPPING_ROLES):
        raise ValueError("column_mapping needs an explicit human count column")
    return mapping


def _validate_count_semantics(
    mapping: Mapping[str, int], raw_semantics: Mapping[str, Any]
) -> Dict[str, str]:
    normalized: Dict[str, str] = {}
    mapped_count_keys = {key for key in mapping if key in _COUNT_MAPPING_ROLES}
    supplied_keys = {
        str(key) for key, value in raw_semantics.items() if value is not None
    }
    if supplied_keys != mapped_count_keys:
        raise ValueError(
            "count_semantics keys must exactly match mapped count columns: "
            f"expected {sorted(mapped_count_keys)}, got {sorted(supplied_keys)}"
        )
    for key in mapped_count_keys:
        role = str(raw_semantics[key]).strip()
        if role not in _COUNT_MAPPING_ROLES[key]:
            raise ValueError(
                f"invalid semantics for {key}: {role!r}; expected one of "
                f"{sorted(_COUNT_MAPPING_ROLES[key])}"
            )
        normalized[key] = role
    return normalized


def _bounded_int(value: Any, name: str, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _trim_text(value: Any, max_chars: int) -> Optional[str]:
    if value is None:
        return None
    text = str(value)
    return text if len(text) <= max_chars else text[: max_chars - 1] + "…"


def _extract_title(source_text: str, explicit_title: Optional[str]) -> str:
    if explicit_title and explicit_title.strip():
        return explicit_title.strip()
    for line in source_text.splitlines()[:80]:
        match = _HEADING_RE.match(line.strip())
        if match and not re.search(
            r"\b(?:abstract|full text|methods?|results?|introduction)\b",
            match.group(2),
            re.IGNORECASE,
        ):
            return match.group(2).strip()
    return ""


def _extract_abstract(source_text: str, max_chars: int = 4_000) -> str:
    lines = source_text.splitlines()
    start: Optional[int] = None
    level = 6
    for idx, line in enumerate(lines):
        heading = _HEADING_RE.match(line.strip())
        if heading and re.search(r"\babstract\b", heading.group(2), re.IGNORECASE):
            start = idx + 1
            level = len(heading.group(1))
            break
    if start is None:
        return ""
    collected: List[str] = []
    for line in lines[start:]:
        heading = _HEADING_RE.match(line.strip())
        if heading and len(heading.group(1)) <= level:
            break
        collected.append(line)
    return _trim_text("\n".join(collected).strip(), max_chars) or ""


def build_compact_context(
    *,
    index: EvidenceIndex,
    gene_symbol: str,
    pmid: str,
    title: Optional[str] = None,
) -> str:
    """Build the initial paper input, embedding normal-size sources in full.

    Each embedded source chunk is paired with an exact locator already accepted
    by :meth:`EvidenceIndex.resolve_locator`.  Sol can therefore quote directly
    from the supplied paper without spending tool calls merely to rediscover
    text that fit safely in the first request.  Only oversized papers use the
    index-first path.
    """

    inventory = index.list_tables(offset=0, limit=EvidenceIndex.MAX_TABLE_LIST)
    full_source = len(index.source_text) <= FULL_CONTEXT_MAX_CHARS
    source_spans = []
    if full_source:
        for start in range(0, len(index.source_text), index.MAX_SOURCE_SPAN):
            end = min(len(index.source_text), start + index.MAX_SOURCE_SPAN)
            source_spans.append(
                {
                    "locator": f"S:{start}:{end}",
                    "text": index.source_text[start:end],
                }
            )
    payload = {
        "paper": {
            "pmid": str(pmid),
            "target_gene": gene_symbol,
            "title": _extract_title(index.source_text, title),
            "abstract": "" if full_source else _extract_abstract(index.source_text),
            "cleaned_source_chars": len(index.source_text),
            "cleaned_source_sha256": index.source_hash,
        },
        "context_delivery": {
            "mode": "full_cleaned_source" if full_source else "indexed_tools",
            "complete": full_source,
            "source_span_count": len(source_spans),
        },
        "source_spans": source_spans,
        "table_inventory": inventory,
        "note": (
            "The source_spans above are the complete cleaned paper. Read them "
            "directly; their S:start:end labels are valid evidence locators. "
            "Do not spend search/read calls on content already supplied. Use "
            "the deterministic table parser when appropriate."
            if full_source
            else "This paper exceeds the direct-context threshold. Search/read "
            "the cleaned source and tables with bounded tools before returning "
            "evidence."
        ),
    }
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def _tool_definitions() -> List[Dict[str, Any]]:
    nullable_column = {"type": ["integer", "null"], "minimum": 0}
    column_properties = {
        key: dict(nullable_column) for key in sorted(_ALLOWED_COLUMN_FIELDS)
    }
    semantic_values = sorted(
        {role for roles in _COUNT_MAPPING_ROLES.values() for role in roles}
    )
    nullable_semantic = {
        "type": ["string", "null"],
        "enum": [*semantic_values, None],
    }
    semantic_properties = {
        key: dict(nullable_semantic) for key in sorted(_COUNT_MAPPING_ROLES)
    }
    return [
        {
            "type": "function",
            "name": "list_tables",
            "description": "List a bounded page of table captions and headers.",
            "parameters": {
                "type": "object",
                "properties": {
                    "offset": {"type": "integer", "minimum": 0},
                    "limit": {"type": "integer", "minimum": 1, "maximum": 20},
                },
                "required": ["offset", "limit"],
                "additionalProperties": False,
            },
            "strict": True,
        },
        {
            "type": "function",
            "name": "search_source",
            "description": (
                "Case-insensitive literal search of the cleaned source. Returns "
                "bounded contexts with exact S:start:end locators."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {"type": "string", "minLength": 2, "maxLength": 200},
                    "max_hits": {
                        "type": "integer",
                        "minimum": 1,
                        "maximum": 12,
                    },
                },
                "required": ["query", "max_hits"],
                "additionalProperties": False,
            },
            "strict": True,
        },
        {
            "type": "function",
            "name": "read_source_span",
            "description": (
                "Read at most 6000 characters from the cleaned source by "
                "zero-based character offsets."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "start": {"type": "integer", "minimum": 0},
                    "end": {"type": "integer", "minimum": 1},
                },
                "required": ["start", "end"],
                "additionalProperties": False,
            },
            "strict": True,
        },
        {
            "type": "function",
            "name": "read_table",
            "description": (
                "Read a bounded inclusive table-row range with stable Tn:Rn locators."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "table_id": {"type": "string", "pattern": "^T[0-9]+$"},
                    "row_start": {"type": "integer", "minimum": 1},
                    "row_end": {"type": "integer", "minimum": 1},
                },
                "required": ["table_id", "row_start", "row_end"],
                "additionalProperties": False,
            },
            "strict": True,
        },
        {
            "type": "function",
            "name": "parse_variant_table",
            "description": (
                "Host-parse a bounded table batch only after identifying "
                "explicit HUMAN PER-VARIANT count-column semantics. The result "
                "is cached under a short batch_id for final acceptance."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "table_id": {"type": "string", "pattern": "^T[0-9]+$"},
                    "row_start": {"type": "integer", "minimum": 1},
                    "row_end": {"type": "integer", "minimum": 1},
                    "column_mapping": {
                        "type": "object",
                        "properties": column_properties,
                        "required": list(column_properties),
                        "additionalProperties": False,
                    },
                    "human_subjects": {"type": "boolean"},
                    "per_variant": {"type": "boolean"},
                    "count_semantics": {
                        "type": "object",
                        "properties": semantic_properties,
                        "required": list(semantic_properties),
                        "additionalProperties": False,
                    },
                },
                "required": [
                    "table_id",
                    "row_start",
                    "row_end",
                    "column_mapping",
                    "human_subjects",
                    "per_variant",
                    "count_semantics",
                ],
                "additionalProperties": False,
            },
            "strict": True,
        },
    ]


def _final_text_format() -> Dict[str, Any]:
    nullable_string = {"type": ["string", "null"]}
    nullable_count = {"type": ["integer", "null"], "minimum": 0}
    direct_semantic_properties = {
        key: {
            "type": ["string", "null"],
            "enum": [expected_role, None],
        }
        for key, expected_role in _DIRECT_COUNT_ROLES.items()
    }
    record_properties = {
        "gene_symbol": {"type": "string"},
        "cdna_notation": nullable_string,
        "protein_notation": nullable_string,
        "genomic_position": nullable_string,
        "clinical_significance": nullable_string,
        "phenotype": nullable_string,
        "total_carriers": nullable_count,
        "affected": nullable_count,
        "unaffected": nullable_count,
        "uncertain": nullable_count,
        "count_semantics": {
            "type": "object",
            "properties": direct_semantic_properties,
            "required": list(direct_semantic_properties),
            "additionalProperties": False,
        },
        "locator": {"type": "string"},
        "evidence_quote": {"type": "string", "minLength": 1},
    }
    schema = {
        "type": "object",
        "properties": {
            "eligibility": {
                "type": "object",
                "properties": {
                    "human_clinical": {"type": "boolean"},
                    "reason": {"type": "string"},
                    "species": nullable_string,
                },
                "required": ["human_clinical", "reason", "species"],
                "additionalProperties": False,
            },
            "accepted_table_batches": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {"batch_id": {"type": "string"}},
                    "required": ["batch_id"],
                    "additionalProperties": False,
                },
            },
            "evidence_records": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": record_properties,
                    "required": list(record_properties),
                    "additionalProperties": False,
                },
            },
        },
        "required": [
            "eligibility",
            "accepted_table_batches",
            "evidence_records",
        ],
        "additionalProperties": False,
    }
    return {
        "format": {
            "type": "json_schema",
            "name": "gvf_sol_extraction",
            "strict": True,
            "schema": schema,
        }
    }


_SOL_INSTRUCTIONS = """You extract target-gene HUMAN clinical variant evidence.

Use only the supplied paper and read-only tools. References/bibliography have
already been removed; the abstract, body, tables, and supplements remain.

Rules:
1. If subjects are animal/model-organism only, set human_clinical=false and
   return no batches or records. Never turn canine/feline/murine samples into
   human patients.
2. Counts must be explicit human counts tied to one target-gene variant. Never
   use study N, families, alleles, chromosomes, controls, tumors, cells,
   replicates, assay n, percentages, or population frequency as carriers.
3. For regular Markdown tables, inspect headers/rows, then prefer
   parse_variant_table. Call it only when human_subjects and per_variant are
   definitely true and declare the exact count role for every mapped count
   column. Accept its batch_id once; do not echo parsed rows.
4. For narrative evidence or tables the parser cannot represent, return one
   evidence_record per variant observation. Copy a verbatim quote and give its
   exact Tn:Rn or S:start:end locator. Every notation and every numeric count in
   a direct record must occur in that quote/span. For every non-null count, set
   its count_semantics field to the matching fixed human_variant_* role; both the
   count and semantic must be null when unsupported. Use null rather than inference.
5. Do not add counts from repeated mentions or apparently different cohorts.
   Preserve separate direct observations; the host resolves duplicates safely.
6. Extract only the requested gene. When context_delivery.mode is
   full_cleaned_source, the supplied source_spans are complete: inspect them
   directly and do not search/read the same content again. Tools remain useful
   for deterministic table parsing. In indexed_tools mode, search/read until all
   relevant tables and narrative evidence have been checked. Then emit the
   strict final JSON.
"""


def _notation_key(value: Any) -> str:
    text = str(value or "")

    def replace_aa(match: re.Match[str]) -> str:
        return _AA3_TO_1.get(match.group(0).lower(), match.group(0))

    text = re.sub(
        r"(?i)Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|"
        r"Ser|Thr|Trp|Tyr|Val|Ter|Stop",
        replace_aa,
        text,
    )
    return re.sub(r"[^A-Za-z0-9*]", "", text).lower()


def _notation_supported(value: str, evidence: str) -> bool:
    if value.casefold() in evidence.casefold():
        return True
    key = _notation_key(value)
    return bool(key and key in _notation_key(evidence))


def _quote_has_integer(quote: str, value: int) -> bool:
    for match in _INTEGER_RE.finditer(quote):
        try:
            if int(match.group(1).replace(",", "")) == value:
                return True
        except ValueError:
            continue
    return False


def _validate_direct_record(
    raw: Any,
    *,
    index: EvidenceIndex,
    gene_symbol: str,
) -> tuple[Optional[Dict[str, Any]], List[str]]:
    fatal_reasons: List[str] = []
    count_reasons: List[str] = []
    if not isinstance(raw, dict):
        return None, ["record is not an object"]
    raw_gene = str(raw.get("gene_symbol") or "").strip()
    if raw_gene.casefold() != gene_symbol.casefold():
        fatal_reasons.append(
            f"gene mismatch: expected {gene_symbol!r}, got {raw_gene!r}"
        )

    locator = str(raw.get("locator") or "").strip()
    quote = str(raw.get("evidence_quote") or "")
    try:
        located_text = index.resolve_locator(locator)
    except ValueError as exc:
        located_text = ""
        fatal_reasons.append(str(exc))
    if not quote.strip():
        fatal_reasons.append("evidence_quote is empty")
    elif quote not in located_text:
        fatal_reasons.append("evidence_quote is not an exact substring of locator")

    identities = [
        (key, str(raw.get(key) or "").strip())
        for key in ("cdna_notation", "protein_notation", "genomic_position")
        if raw.get(key)
    ]
    if not identities:
        fatal_reasons.append("record has no variant notation")
    for key, value in identities:
        if quote and not _notation_supported(value, quote):
            fatal_reasons.append(f"{key} is not supported by evidence_quote")

    counts: Dict[str, Optional[int]] = {}
    semantics = raw.get("count_semantics")
    if not isinstance(semantics, dict):
        semantics = {}
    for direct_key, gvf_key in _DIRECT_TO_GVF_COUNT.items():
        value = raw.get(direct_key)
        semantic = semantics.get(direct_key)
        if value is None:
            counts[gvf_key] = None
            if semantic is not None:
                count_reasons.append(
                    f"{direct_key} semantic was supplied without a count and was ignored"
                )
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            count_reasons.append(
                f"{direct_key} cleared: must be a non-negative integer or null"
            )
            counts[gvf_key] = None
            continue
        if semantic != _DIRECT_COUNT_ROLES[direct_key]:
            count_reasons.append(
                f"{direct_key}={value} cleared: count_semantics must be "
                f"{_DIRECT_COUNT_ROLES[direct_key]!r}"
            )
            counts[gvf_key] = None
            continue
        if not quote or not _quote_has_integer(quote, value):
            count_reasons.append(
                f"{direct_key}={value} cleared: value does not occur in evidence_quote"
            )
            counts[gvf_key] = None
            continue
        counts[gvf_key] = value

    total = counts.get("total_carriers_observed")
    affected = counts.get("affected_count")
    unaffected = counts.get("unaffected_count")
    uncertain = counts.get("uncertain_count")
    known_parts = [v for v in (affected, unaffected, uncertain) if v is not None]
    if total is not None and sum(known_parts) > total:
        count_reasons.append(
            "all counts cleared: affected + unaffected + uncertain exceeds "
            "total_carriers"
        )
        counts = {key: None for key in _COUNT_KEYS}
    if fatal_reasons:
        return None, [*fatal_reasons, *count_reasons]
    return (
        _direct_record_to_variant(
            raw,
            gene_symbol=gene_symbol,
            locator=locator,
            quote=quote,
            counts=counts,
        ),
        count_reasons,
    )


def _direct_record_to_variant(
    raw: Mapping[str, Any],
    *,
    gene_symbol: str,
    locator: str,
    quote: str,
    counts: Mapping[str, Optional[int]],
) -> Dict[str, Any]:
    total = counts.get("total_carriers_observed")
    fact_rows: List[Dict[str, Any]] = []
    for notation_key in ("protein_notation", "cdna_notation", "genomic_position"):
        value = raw.get(notation_key)
        if value:
            fact_rows.append(
                {
                    "fact_type": "variant_identity",
                    "fact_value": value,
                    "source_location": locator,
                    "evidence_quote": quote,
                    "source_layer": (
                        "llm_table" if locator.upper().startswith("T") else "llm_text"
                    ),
                    "provenance_kind": "source_grounded_sol",
                }
            )
    for fact_type, value in counts.items():
        if value is not None:
            fact_rows.append(
                {
                    "fact_type": fact_type,
                    "fact_value": value,
                    "source_location": locator,
                    "evidence_quote": quote,
                    "count_type": "per_variant_carrier",
                    "source_layer": (
                        "llm_table" if locator.upper().startswith("T") else "llm_text"
                    ),
                    "provenance_kind": "source_grounded_sol",
                }
            )
    clinical = str(raw.get("clinical_significance") or "").strip() or None
    phenotype = str(raw.get("phenotype") or "").strip() or None
    return {
        "gene_symbol": gene_symbol,
        "cdna_notation": raw.get("cdna_notation"),
        "protein_notation": raw.get("protein_notation"),
        "genomic_position": raw.get("genomic_position"),
        "clinical_significance": clinical,
        "patients": {
            "count": total,
            "phenotype": phenotype,
            "source_kind": "table" if locator.upper().startswith("T") else "text",
            "source_ref": locator,
            "locator_extra": {"stable_locator": locator},
        },
        "penetrance_data": {
            **counts,
            "penetrance_percentage": None,
            "age_dependent_penetrance": [],
        },
        "individual_records": [],
        "functional_data": {"summary": "", "assays": []},
        "segregation_data": None,
        "population_frequency": None,
        "evidence_level": "source_grounded",
        "source_location": locator,
        "evidence_quote": quote,
        "source_layer": "llm_table" if locator.upper().startswith("T") else "llm_text",
        "additional_notes": "GPT-5.6 Sol direct source-grounded evidence",
        "key_quotes": [quote],
        "count_provenance": {
            "carriers_count_type": (
                "per_variant_carrier" if total is not None else None
            ),
            "affected_count_type": (
                "per_variant_carrier"
                if counts.get("affected_count") is not None
                else None
            ),
            "unaffected_count_type": (
                "per_variant_carrier"
                if counts.get("unaffected_count") is not None
                else None
            ),
        },
        "fact_provenance": fact_rows,
        "locator_extra": {"stable_locator": locator},
    }


def _parsed_variant_as_direct_record(variant: Mapping[str, Any]) -> Dict[str, Any]:
    penetrance = variant.get("penetrance_data") or {}
    count_values = {
        "total_carriers": penetrance.get("total_carriers_observed"),
        "affected": penetrance.get("affected_count"),
        "unaffected": penetrance.get("unaffected_count"),
        "uncertain": penetrance.get("uncertain_count"),
    }
    return {
        "gene_symbol": variant.get("gene_symbol"),
        "cdna_notation": variant.get("cdna_notation"),
        "protein_notation": variant.get("protein_notation"),
        "genomic_position": variant.get("genomic_position"),
        "clinical_significance": variant.get("clinical_significance"),
        "phenotype": (variant.get("patients") or {}).get("phenotype"),
        **count_values,
        "count_semantics": {
            key: (_DIRECT_COUNT_ROLES[key] if value is not None else None)
            for key, value in count_values.items()
        },
        "locator": variant.get("source_location"),
        "evidence_quote": variant.get("evidence_quote")
        or ((variant.get("key_quotes") or [""])[0]),
    }


def _identity_key(variant: Mapping[str, Any]) -> tuple[str, str]:
    for field in ("cdna_notation", "protein_notation", "genomic_position"):
        value = variant.get(field)
        if value:
            return field, _notation_key(value)
    return "unknown", hashlib.sha256(
        json.dumps(dict(variant), sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()


def _merge_duplicate_variants(
    variants: Sequence[Dict[str, Any]], audit: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Merge repeated evidence without assuming observations are additive."""

    merged: Dict[tuple[str, str], Dict[str, Any]] = {}
    order: List[tuple[str, str]] = []
    for incoming in variants:
        key = _identity_key(incoming)
        if key not in merged:
            merged[key] = incoming
            order.append(key)
            continue
        existing = merged[key]
        existing_pen = existing.setdefault("penetrance_data", {})
        incoming_pen = incoming.get("penetrance_data") or {}
        conflicting: List[str] = []
        for count_key in _COUNT_KEYS:
            old = existing_pen.get(count_key)
            new = incoming_pen.get(count_key)
            if old is None:
                existing_pen[count_key] = new
            elif new is None or old == new:
                continue
            else:
                conflicting.append(count_key)
                existing_pen[count_key] = None
        existing.setdefault("patients", {})["count"] = existing_pen.get(
            "total_carriers_observed"
        )

        for list_key in ("key_quotes", "fact_provenance"):
            target = existing.setdefault(list_key, [])
            for item in incoming.get(list_key, []) or []:
                if item not in target:
                    target.append(item)
        if conflicting:
            audit.append(
                {
                    "reason": "conflicting_duplicate_counts_not_summed",
                    "identity": {"field": key[0], "normalized": key[1]},
                    "fields_cleared": conflicting,
                    "first_locator": existing.get("source_location"),
                    "other_locator": incoming.get("source_location"),
                }
            )
            note = existing.get("additional_notes") or ""
            existing["additional_notes"] = (
                note
                + "; conflicting duplicate counts retained in provenance and cleared"
            ).strip("; ")
    return [merged[key] for key in order]


def _serializable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if is_dataclass(value):
        return {key: _serializable(item) for key, item in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(key): _serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serializable(item) for item in value]
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return _serializable(to_dict())
    return str(value)


def _parse_response_value(result: Any) -> Optional[Dict[str, Any]]:
    value = getattr(result, "value", None)
    if isinstance(value, dict):
        return value
    candidates = [value, getattr(result, "output_text", None)]
    for candidate in candidates:
        if not isinstance(candidate, str) or not candidate.strip():
            continue
        text = candidate.strip()
        if text.startswith("```"):
            text = re.sub(r"^```(?:json)?\s*|\s*```$", "", text, flags=re.I)
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            continue
        if isinstance(parsed, dict):
            return parsed
    return None


def _empty_output(
    *,
    pmid: str,
    title: str,
    gene_symbol: str,
    source_hash: str,
    cleaned_hash: str,
    model: str,
    reasoning_effort: str,
    eligibility: Mapping[str, Any],
    audit: Optional[List[Dict[str, Any]]] = None,
    response_meta: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    return {
        "paper_metadata": {
            "pmid": str(pmid),
            "title": title or f"Paper {pmid}",
            "gene_symbol": gene_symbol,
            "extraction_summary": str(eligibility.get("reason") or "No evidence"),
        },
        "variants": [],
        "tables_processed": [],
        "extraction_metadata": {
            "total_variants_found": 0,
            "model_used": model,
            "reasoning_effort": reasoning_effort,
            "protocol_version": SOL_PROTOCOL_VERSION,
            "cleaner_version": SOL_CLEANER_VERSION,
            "source_sha256": source_hash,
            "cleaned_source_sha256": cleaned_hash,
            "eligibility": dict(eligibility),
            "sol_audit": {"rejected_raw_facts": audit or []},
            "responses": dict(response_meta or {}),
        },
    }


class SolExtractor:
    """One-model, evidence-tool extraction lane for GPT-5.6 Sol."""

    def __init__(
        self,
        runner: ResponsesRunner,
        *,
        model: str = "gpt-5.6-sol",
        reasoning_effort: str = "xhigh",
        limits: Any = None,
        max_output_tokens: int = 100_000,
    ) -> None:
        self.runner = runner
        self.model = model
        self.reasoning_effort = reasoning_effort
        self.limits = limits
        self.max_output_tokens = max_output_tokens

    def extract(
        self,
        *,
        source_text: str,
        gene_symbol: str,
        pmid: str,
        title: Optional[str] = None,
    ) -> Dict[str, Any]:
        source_text = source_text or ""
        gene_symbol = str(gene_symbol or "").strip().upper()
        if not gene_symbol:
            raise ValueError("gene_symbol is required")
        pmid = str(pmid or "").strip()
        if not pmid:
            raise ValueError("pmid is required")

        raw_hash = hashlib.sha256(source_text.encode("utf-8")).hexdigest()
        cleaned = strip_irrelevant_context(source_text)
        cleaned_hash = hashlib.sha256(cleaned.encode("utf-8")).hexdigest()
        resolved_title = _extract_title(cleaned, title)

        # This invariant is cheaper and more reliable than asking any model.
        from pipeline.steps import (
            _apply_nonhuman_clinical_count_guard,
            _nonhuman_veterinary_source_reason,
        )

        nonhuman_reason = _nonhuman_veterinary_source_reason(cleaned)
        if nonhuman_reason:
            output = _empty_output(
                pmid=pmid,
                title=resolved_title,
                gene_symbol=gene_symbol,
                source_hash=raw_hash,
                cleaned_hash=cleaned_hash,
                model=self.model,
                reasoning_effort=self.reasoning_effort,
                eligibility={
                    "human_clinical": False,
                    "reason": nonhuman_reason,
                    "species": "non-human/veterinary",
                    "deterministic_short_circuit": True,
                },
            )
            _apply_nonhuman_clinical_count_guard(output, cleaned)
            return output

        index = EvidenceIndex(cleaned)
        batches = _ParsedBatchStore(index, gene_symbol)
        tool_handlers: Dict[str, Callable[[Dict[str, Any]], Any]] = {
            "list_tables": lambda args: index.list_tables(**args),
            "search_source": lambda args: index.search_source(**args),
            "read_source_span": lambda args: index.read_source_span(**args),
            "read_table": lambda args: index.read_table(**args),
            "parse_variant_table": batches.parse,
        }
        initial_input = build_compact_context(
            index=index,
            gene_symbol=gene_symbol,
            pmid=pmid,
            title=resolved_title,
        )
        context_mode = (
            "full_cleaned_source"
            if len(index.source_text) <= FULL_CONTEXT_MAX_CHARS
            else "indexed_tools"
        )

        run_kwargs: Dict[str, Any] = {
            "model": self.model,
            "reasoning_effort": self.reasoning_effort,
            "instructions": _SOL_INSTRUCTIONS,
            "initial_input": initial_input,
            "tools": _tool_definitions(),
            "tool_handlers": tool_handlers,
            "text": _final_text_format(),
            "max_output_tokens": self.max_output_tokens,
        }
        if self.limits is not None:
            run_kwargs["limits"] = self.limits
        result = self.runner.run_tool_loop(**run_kwargs)
        response_meta = {
            "ok": bool(getattr(result, "ok", False)),
            "status": getattr(result, "status", None),
            "response_id": getattr(result, "response_id", None),
            "error": getattr(result, "error", None),
            "usage": _serializable(getattr(result, "usage", None)),
            "telemetry": _serializable(getattr(result, "telemetry", None)),
        }
        final = _parse_response_value(result)
        if not bool(getattr(result, "ok", False)) or final is None:
            audit = [
                {
                    "reason": "responses_tool_loop_failed_or_unparseable",
                    "status": response_meta["status"],
                    "error": response_meta["error"],
                }
            ]
            output = _empty_output(
                pmid=pmid,
                title=resolved_title,
                gene_symbol=gene_symbol,
                source_hash=raw_hash,
                cleaned_hash=cleaned_hash,
                model=self.model,
                reasoning_effort=self.reasoning_effort,
                eligibility={
                    "human_clinical": False,
                    "reason": "model/tool loop failed",
                    "species": None,
                    "deterministic_short_circuit": False,
                },
                audit=audit,
                response_meta=response_meta,
            )
            _apply_nonhuman_clinical_count_guard(output, cleaned)
            return output

        eligibility = final.get("eligibility")
        if not isinstance(eligibility, dict):
            eligibility = {
                "human_clinical": False,
                "reason": "invalid or missing eligibility object",
                "species": None,
            }
        rejected: List[Dict[str, Any]] = []
        candidates: List[Dict[str, Any]] = []
        accepted_ids: set[str] = set()
        seen_batch_ids: set[str] = set()

        if eligibility.get("human_clinical") is not True:
            for raw in final.get("evidence_records") or []:
                rejected.append(
                    {
                        "reason": "model classified paper as ineligible",
                        "raw": _serializable(raw),
                    }
                )
            for raw_batch in final.get("accepted_table_batches") or []:
                rejected.append(
                    {
                        "reason": "model classified paper as ineligible",
                        "raw": _serializable(raw_batch),
                    }
                )
        else:
            for raw_batch in final.get("accepted_table_batches") or []:
                batch_id = (
                    str(raw_batch.get("batch_id") or "").strip()
                    if isinstance(raw_batch, dict)
                    else ""
                )
                if not batch_id or batch_id in seen_batch_ids:
                    if batch_id in seen_batch_ids:
                        rejected.append(
                            {
                                "reason": "duplicate accepted batch reference",
                                "raw": _serializable(raw_batch),
                            }
                        )
                    continue
                seen_batch_ids.add(batch_id)
                batch = batches.batches.get(batch_id)
                if batch is None:
                    rejected.append(
                        {
                            "reason": "unknown or unexecuted table batch",
                            "raw": _serializable(raw_batch),
                        }
                    )
                    continue
                accepted_ids.add(batch_id)
                for parsed_variant in batch.variants:
                    direct_raw = _parsed_variant_as_direct_record(parsed_variant)
                    validated, reasons = _validate_direct_record(
                        direct_raw, index=index, gene_symbol=gene_symbol
                    )
                    if validated is None:
                        rejected.append(
                            {
                                "reason": "; ".join(reasons),
                                "origin": batch_id,
                                "raw": _serializable(direct_raw),
                            }
                        )
                    else:
                        if reasons:
                            rejected.append(
                                {
                                    "reason": "unsupported_count_fields_cleared",
                                    "details": reasons,
                                    "origin": batch_id,
                                    "raw": _serializable(direct_raw),
                                }
                            )
                        # The compact adapter is intentionally common to direct
                        # and deterministic candidates, which gives both paths
                        # identical source/count validation and migration shape.
                        candidates.append(validated)

            for raw in final.get("evidence_records") or []:
                validated, reasons = _validate_direct_record(
                    raw, index=index, gene_symbol=gene_symbol
                )
                if validated is None:
                    rejected.append(
                        {
                            "reason": "; ".join(reasons),
                            "origin": "direct",
                            "raw": _serializable(raw),
                        }
                    )
                else:
                    if reasons:
                        rejected.append(
                            {
                                "reason": "unsupported_count_fields_cleared",
                                "details": reasons,
                                "origin": "direct",
                                "raw": _serializable(raw),
                            }
                        )
                    candidates.append(validated)

        conflict_audit: List[Dict[str, Any]] = []
        variants = _merge_duplicate_variants(candidates, conflict_audit)
        rejected.extend(conflict_audit)
        output: Dict[str, Any] = {
            "paper_metadata": {
                "pmid": pmid,
                "title": resolved_title or f"Paper {pmid}",
                "gene_symbol": gene_symbol,
                "extraction_summary": (
                    f"{len(variants)} source-grounded target-gene variants; "
                    f"{len(rejected)} rejected/audited facts"
                ),
            },
            "variants": variants,
            "tables_processed": [],
            "extraction_metadata": {
                "total_variants_found": len(variants),
                "model_used": self.model,
                "reasoning_effort": self.reasoning_effort,
                "protocol_version": SOL_PROTOCOL_VERSION,
                "cleaner_version": SOL_CLEANER_VERSION,
                "source_sha256": raw_hash,
                "cleaned_source_sha256": cleaned_hash,
                "eligibility": _serializable(eligibility),
                "sol_audit": {
                    "context_delivery": {
                        "mode": context_mode,
                        "full_context_max_chars": FULL_CONTEXT_MAX_CHARS,
                        "source_spans_supplied": (
                            (len(index.source_text) + index.MAX_SOURCE_SPAN - 1)
                            // index.MAX_SOURCE_SPAN
                            if context_mode == "full_cleaned_source"
                            else 0
                        ),
                    },
                    "table_inventory": [
                        {
                            "table_id": table.table_id,
                            "row_count": len(table.data_lines),
                        }
                        for table in index.tables
                    ],
                    "accepted_table_batch_ids": sorted(accepted_ids),
                    "executed_table_batches": [
                        {
                            "batch_id": batch.batch_id,
                            "table_id": batch.table_id,
                            "row_start": batch.row_start,
                            "row_end": batch.row_end,
                            "candidate_count": len(batch.variants),
                            "accepted": batch.batch_id in accepted_ids,
                        }
                        for batch in batches.batches.values()
                    ],
                    "accepted_candidate_facts": len(candidates),
                    "rejected_raw_facts": rejected,
                },
                "responses": response_meta,
            },
        }
        # Defense in depth: if future deterministic detection improves, the
        # established guard still clears any clinical count that slipped here.
        _apply_nonhuman_clinical_count_guard(output, cleaned)
        return output


def extract_paper(
    *,
    gene: str,
    pmid: str,
    source_text: str,
    source_sha256: Optional[str] = None,
    reasoning_effort: str,
    model: str = "azure_ai/gpt-5.6-sol",
    title: Optional[str] = None,
    runner: Optional[ResponsesRunner] = None,
) -> tuple[Dict[str, Any], Dict[str, Any]]:
    """Driver-facing extraction entry point with a frozen identity receipt.

    ``model`` is retained verbatim in telemetry; ``ResponsesAPIClient`` removes
    provider prefixes only in the wire request.  A supplied source hash is an
    assertion, never trusted metadata, so a changed fixture fails before spend.
    """

    from utils.responses_api import (
        REASONING_EFFORTS,
        ResponsesAPIClient,
        ResponsesLoopLimits,
    )

    if reasoning_effort not in REASONING_EFFORTS:
        raise ValueError(
            "reasoning_effort must be one of: " + ", ".join(REASONING_EFFORTS)
        )
    observed_source_hash = hashlib.sha256(
        (source_text or "").encode("utf-8")
    ).hexdigest()
    if source_sha256 is not None and (
        str(source_sha256).strip().lower() != observed_source_hash
    ):
        raise ValueError(
            "source_sha256 mismatch: "
            f"expected {source_sha256}, observed {observed_source_hash}"
        )
    cleaned = strip_irrelevant_context(source_text or "")
    cleaned_hash = hashlib.sha256(cleaned.encode("utf-8")).hexdigest()

    active_runner = runner
    if active_runner is None:
        active_runner = ResponsesAPIClient.azure(
            os.environ["AZURE_AI_API_BASE"],
            os.environ["AZURE_AI_API_KEY"],
        )
    extraction = SolExtractor(
        active_runner,
        model=model,
        reasoning_effort=reasoning_effort,
        limits=ResponsesLoopLimits(
            max_rounds=32,
            max_tool_calls=128,
            max_tool_output_chars=30_000,
            max_total_tool_output_chars=640_000,
        ),
        max_output_tokens=100_000,
    ).extract(
        source_text=source_text,
        gene_symbol=gene,
        pmid=pmid,
        title=title,
    )

    extraction_meta = extraction.get("extraction_metadata") or {}
    response = extraction_meta.get("responses") or {}
    eligibility = extraction_meta.get("eligibility") or {}
    response_ok = response.get("ok") is True
    deterministic_short_circuit = eligibility.get("deterministic_short_circuit") is True
    if deterministic_short_circuit:
        status = "success"
        model_status = "deterministic_nonhuman_short_circuit"
        effective_effort: Optional[str] = None
    elif response_ok:
        status = "success"
        model_status = str(response.get("status") or "completed")
        effective_effort = reasoning_effort
    else:
        status = "model_error"
        model_status = str(response.get("status") or "failed")
        effective_effort = None

    telemetry: Dict[str, Any] = {
        "status": status,
        "model_status": model_status,
        "model": model,
        "reasoning_effort": reasoning_effort,
        "requested_reasoning_effort": reasoning_effort,
        # This is populated only after the API actually accepted the effort.
        "effective_reasoning_effort": effective_effort,
        "source_sha256": observed_source_hash,
        "cleaned_source_sha256": cleaned_hash,
        "source_chars": len(source_text or ""),
        "cleaned_source_chars": len(cleaned),
        "protocol_version": SOL_PROTOCOL_VERSION,
        "cleaner_version": SOL_CLEANER_VERSION,
        "usage": response.get("usage"),
        "response_telemetry": response.get("telemetry"),
        "response_id": response.get("response_id"),
        "response_error": response.get("error"),
        "response_status": response.get("status"),
    }
    return extraction, telemetry


__all__ = [
    "FULL_CONTEXT_MAX_CHARS",
    "SOL_CLEANER_VERSION",
    "SOL_PROTOCOL_VERSION",
    "EvidenceIndex",
    "ResponsesRunner",
    "SolExtractor",
    "build_compact_context",
    "clean_context",
    "extract_paper",
    "strip_irrelevant_context",
]
