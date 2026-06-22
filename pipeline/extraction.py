"""
Expert Extractor module (Tier 3) for the Tiered Biomedical Extraction Pipeline.

The heavy lifter that processes full-text papers and extracts structured
genetic variant data using advanced LLM prompting.
"""

import copy
import json
import logging
import os
import re
from pathlib import Path
from typing import Any, List, Optional

from config.constants import (
    ADAPTIVE_TABLE_THRESHOLD,
    CONTINUATION_VARIANT_GAP,
    DETERMINISTIC_PARSER_MIN_VARIANTS,
    GENE_CONTEXT_WINDOW_LINES,
    LARGE_TABLE_ROW_THRESHOLD,
    MIN_ALPHANUMERIC_RATIO,
    MIN_CONDENSED_SIZE,
    MIN_EXTRACTION_INPUT_SIZE,
    SCANNER_MAX_HINTS,
    SCANNER_MERGE_MIN_CONFIDENCE,
    TABLE_MIN_COLUMNS,
    TEXT_TRUNCATION_MAX_CHARS,
)
from config.settings import get_settings
from pipeline.prompts import (
    COMPACT_EXTRACTION_PROMPT,
    CONTINUATION_PROMPT,
    EXTRACTION_PROMPT,
    HIGH_VARIANT_THRESHOLD,
)
from utils.gene_metadata import gene_alias_regex, known_gene_aliases
from utils.llm_utils import BaseLLMCaller, clamp_max_tokens
from utils.models import ExtractionResult, Paper
from utils.source_layers import infer_source_layer_from_text
from utils.variant_scanner import (
    VariantScanner,
    merge_scanner_results,
    scan_document_for_variants,
)

logger = logging.getLogger(__name__)

TABLE_HINT_MAX_VARIANTS = SCANNER_MAX_HINTS
TABLE_REGEX_MERGE_MAX_VARIANTS = 500
SCANNER_REGEX_MERGE_MAX_VARIANTS = 150
TABLE_REGEX_OVERFLOW_MERGE_MAX_VARIANTS = int(
    os.environ.get("GVF_TABLE_OVERFLOW_MERGE_MAX_VARIANTS", "2000")
)
TABLE_REGEX_OVERFLOW_CHUNK_SIZE = int(
    os.environ.get("GVF_TABLE_OVERFLOW_CHUNK_SIZE", "250")
)


def _find_data_zones_file(
    pmid: str, search_dirs: Optional[List[str]] = None
) -> Optional[Path]:
    """
    Search for a DATA_ZONES.md file for the given PMID.

    .. deprecated::
        This function uses directory scanning which is inefficient and fragile.
        Prefer using manifest-based file discovery via `cli.extract` which reads
        file paths directly from scout_manifest.json.

    Args:
        pmid: PubMed ID to search for
        search_dirs: Optional list of directories to search in. If provided with
                     a single directory, searches that directory directly first.

    Returns:
        Path to DATA_ZONES.md if found, None otherwise
    """
    import warnings

    filename = f"{pmid}_DATA_ZONES.md"

    # If explicit search directories provided, check them first
    if search_dirs:
        for search_dir in search_dirs:
            if search_dir is None:
                continue
            path = Path(search_dir)
            if path.exists() and path.is_dir():
                # Direct match in the specified directory
                zones_file = path / filename
                if zones_file.exists():
                    logger.debug(f"Found DATA_ZONES.md at {zones_file}")
                    return zones_file

                # Search subdirectories (one level deep)
                for subdir in path.iterdir():
                    if subdir.is_dir():
                        zones_file = subdir / filename
                        if zones_file.exists():
                            logger.debug(f"Found DATA_ZONES.md at {zones_file}")
                            return zones_file

    # Fallback: search common output directory patterns (DEPRECATED behavior)
    fallback_dirs = [".", "pmc_fulltext", "output"]
    used_fallback = False
    for search_dir in fallback_dirs:
        path = Path(search_dir)
        if path.exists() and path.is_dir():
            zones_file = path / filename
            if zones_file.exists():
                if not used_fallback:
                    warnings.warn(
                        f"Using deprecated directory scanning fallback for {pmid}. "
                        "Use cli.extract with --manifest for better reliability.",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                    used_fallback = True
                logger.debug(f"Found DATA_ZONES.md at {zones_file} (fallback)")
                return zones_file

    # Also search subdirectories matching test/output patterns (DEPRECATED)
    cwd = Path(".")
    for subdir in cwd.iterdir():
        if subdir.is_dir() and (
            subdir.name.startswith("test_")
            or subdir.name.startswith("output_")
            or pmid in subdir.name
        ):
            zones_file = subdir / filename
            if zones_file.exists():
                if not used_fallback:
                    warnings.warn(
                        f"Using deprecated directory scanning fallback for {pmid}. "
                        "Use cli.extract with --manifest for better reliability.",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                logger.debug(f"Found DATA_ZONES.md at {zones_file} (subdir fallback)")
                return zones_file

    return None


class ExpertExtractor(BaseLLMCaller):
    """
    Tier 3: Expert-level extraction using advanced LLM (GPT-4, etc.).
    Handles full-text papers and markdown tables to extract structured variant data.

    Prompts are defined in pipeline/prompts.py for easier maintenance.
    """

    def __init__(
        self,
        models: Optional[List[str]] = None,
        temperature: Optional[float] = None,
        max_tokens: Optional[int] = None,
        reasoning_effort: Optional[str] = None,
        tier_threshold: int = 1,
        fulltext_dir: Optional[str] = None,
    ):
        """
        Initialize the Expert Extractor.

        Args:
            models: List of LiteLLM model identifiers. If None, uses config (TIER3_MODELS).
            temperature: Model temperature. If None, uses config.
            max_tokens: Maximum tokens for response. If None, uses config.
            tier_threshold: If the first model finds fewer variants than this, the next model is tried.
            fulltext_dir: Directory where full-text/DATA_ZONES files are stored.
        """
        settings = get_settings()

        self.models = models or settings.get_tier3_models()
        self.temperature = (
            temperature if temperature is not None else settings.tier3_temperature
        )
        # Store the requested max tokens, but clamp per-model to avoid API errors
        self.requested_max_tokens = (
            max_tokens if max_tokens is not None else settings.tier3_max_tokens
        )
        self.max_tokens = self._clamp_max_tokens(
            self.models[0], self.requested_max_tokens
        )
        self.tier_threshold = tier_threshold
        self.fulltext_dir = fulltext_dir
        self.use_condensed = settings.scout_use_condensed
        self.enable_ensemble_qa = settings.enable_tier3_ensemble_qa
        self.adjudicator_models = settings.get_tier3_adjudicator_models()
        self.adjudication_risk_threshold = settings.tier3_adjudication_risk_threshold
        self.evidence_packet_max_chars = settings.tier3_evidence_packet_max_chars
        self.adjudication_max_tokens = settings.tier3_adjudication_max_tokens
        self.max_verifier_cards = settings.tier3_max_verifier_cards

        super().__init__(
            model=self.models[0],
            temperature=self.temperature,
            max_tokens=self.max_tokens,
            reasoning_effort=(
                reasoning_effort
                if reasoning_effort is not None
                else settings.tier3_reasoning_effort
            ),
        )
        logger.debug(
            f"ExpertExtractor initialized with models={self.models}, temp={self.temperature}, max_tokens={self.max_tokens}"
        )

    def _clamp_max_tokens(self, model: str, requested: int) -> int:
        """Clamp max_tokens to model-safe limits to avoid provider errors."""
        return clamp_max_tokens(model, requested)

    # Use centralized constant for minimum condensed size threshold
    MIN_CONDENSED_SIZE = MIN_CONDENSED_SIZE

    def _prepare_full_text(self, paper: Paper) -> str:
        """
        Prepare full text for extraction.

        If scout_use_condensed is enabled and a DATA_ZONES.md file exists,
        prefer the condensed version for more focused extraction.

        IMPORTANT: Falls back to full text if DATA_ZONES.md is too small,
        as this indicates Data Scout failed to identify useful zones.
        """
        # Try to use condensed DATA_ZONES.md if enabled
        if self.use_condensed and paper.pmid:
            search_dirs = [self.fulltext_dir] if self.fulltext_dir else None
            zones_file = _find_data_zones_file(paper.pmid, search_dirs)

            if zones_file:
                try:
                    condensed_text = zones_file.read_text(encoding="utf-8")

                    # Check if condensed text is useful (not just headers/metadata)
                    # Look for actual data content, not just "No high-value data zones identified"
                    has_no_zones = (
                        "No high-value data zones identified" in condensed_text
                    )
                    is_too_small = len(condensed_text) < self.MIN_CONDENSED_SIZE

                    # Content-aware safety: reject DATA_ZONES that look "big
                    # enough" but don't even mention the target gene. Without
                    # this we can ship 600-byte zone files identifying
                    # off-target receptors as "high-value" and never see the
                    # real text. Skips when we can't determine the gene.
                    gene = (paper.gene_symbol or "").strip()
                    missing_gene = bool(
                        gene and gene.lower() not in condensed_text.lower()
                    )

                    if has_no_zones or is_too_small or missing_gene:
                        reason = (
                            "no zones identified"
                            if has_no_zones
                            else "too small"
                            if is_too_small
                            else f"target gene {gene!r} not mentioned"
                        )
                        logger.warning(
                            f"PMID {paper.pmid} - DATA_ZONES.md unusable "
                            f"({len(condensed_text)} chars, {reason}); "
                            f"falling back to full text"
                        )
                        print(
                            f"⚠ DATA_ZONES.md insufficient ({len(condensed_text)} chars, {reason}) - using full text instead"
                        )
                        # Fall through to use full text
                    elif (
                        condensed_text and len(condensed_text) > self.MIN_CONDENSED_SIZE
                    ):
                        lines = len(condensed_text.splitlines())
                        logger.info(
                            f"PMID {paper.pmid} - Using condensed {paper.pmid}_DATA_ZONES.md ({len(condensed_text)} chars)"
                        )
                        print(
                            f"Using {paper.pmid}_DATA_ZONES.md for extraction: {len(condensed_text):,} chars, {lines:,} lines"
                        )
                        return condensed_text
                except Exception as e:
                    logger.warning(
                        f"PMID {paper.pmid} - Failed to read {paper.pmid}_DATA_ZONES.md: {e}"
                    )

        # Fall back to paper.full_text
        if paper.full_text:
            lines = len(paper.full_text.splitlines())
            print(
                f"Using {paper.pmid}_FULL_CONTEXT.md for extraction: {len(paper.full_text):,} chars, {lines:,} lines"
            )
            return paper.full_text
        elif paper.abstract:
            logger.warning(
                f"PMID {paper.pmid} - Full text not available, using abstract only"
            )
            print(f"Using ABSTRACT ONLY for extraction: {len(paper.abstract):,} chars")
            return f"[ABSTRACT ONLY - FULL TEXT NOT AVAILABLE]\n\n{paper.abstract}"
        else:
            return "[NO TEXT AVAILABLE]"

    def _merge_segments(self, segments: List[tuple], total_lines: int) -> List[tuple]:
        """Merge overlapping (start, end) line index segments."""
        if not segments:
            return []
        merged = []
        for start, end in sorted(segments, key=lambda x: x[0]):
            start = max(0, start)
            end = min(total_lines, end)
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        return [(s, e) for s, e in merged]

    def _gene_focused_truncation(
        self,
        full_text: str,
        gene_symbol: Optional[str],
        max_chars: int,
        context_window: int = GENE_CONTEXT_WINDOW_LINES,
    ) -> Optional[str]:
        """
        Prioritize slices that mention the target gene or contain gene-specific tables.
        Returns concatenated segments if they fit within max_chars; otherwise None.
        """
        if not gene_symbol:
            return None

        lines = full_text.splitlines()
        lower_gene = gene_symbol.lower()
        segments: List[tuple] = []

        for idx, line in enumerate(lines):
            lower_line = line.lower()

            # Capture any table caption that pairs the gene with a table label
            if "table" in lower_line and lower_gene in lower_line:
                end = idx + 1
                while end < len(lines) and (
                    lines[end].strip() == ""
                    or lines[end].lstrip().startswith("|")
                    or lines[end].strip().startswith("*")
                ):
                    end += 1
                segments.append((max(0, idx - 5), end))
                continue

            # Generic gene mentions get a surrounding window
            if lower_gene in lower_line:
                segments.append((max(0, idx - context_window), idx + context_window))

        merged = self._merge_segments(segments, len(lines))
        if not merged:
            return None

        pieces = []
        total_len = 0
        for start, end in merged:
            block = "\n".join(lines[start:end])
            block_len = len(block)
            if total_len + block_len > max_chars:
                # Add as much as possible from this block
                remaining = max_chars - total_len
                if remaining <= 0:
                    break
                block = block[:remaining]
                block_len = len(block)
            pieces.append(block)
            total_len += block_len
            if total_len >= max_chars:
                break

        focused = (
            f"[GENE-FOCUSED TRUNCATION for {gene_symbol}]\n\n"
            + "\n\n---\n\n".join(pieces)
        )
        return focused[:max_chars] if focused else None

    def _truncate_text_for_prompt(
        self,
        full_text: str,
        gene_symbol: Optional[str] = None,
        max_chars: int = TEXT_TRUNCATION_MAX_CHARS,
    ) -> str:
        """
        Keep three slices (head/mid/tail) so large supplemental tables aren't dropped.
        The middle slice is centered on the last table/supplement mention when found.
        """
        if len(full_text) <= max_chars:
            return full_text

        focused = self._gene_focused_truncation(full_text, gene_symbol, max_chars)
        if focused:
            return focused

        marker = "\n\n[TRUNCATED FOR PROMPT - non-adjacent segments preserved to keep tables]\n\n"
        mid_header = "### MIDDLE SEGMENT (around tables/supplements) ###\n"
        tail_header = "### TRAILING SEGMENT ###\n"

        # Allocate space across three segments
        fixed_overhead = len(marker) * 2 + len(mid_header) + len(tail_header)
        available = max_chars - fixed_overhead
        if available < 0:
            available = max_chars

        # Distribute roughly evenly
        head_len = available // 3
        mid_len = available // 3
        tail_len = available - head_len - mid_len

        lower_text = full_text.lower()

        # Anchor the middle slice near the last table/supplement mention
        anchors = [
            lower_text.rfind("table 2"),
            lower_text.rfind("table"),
            lower_text.rfind("supplement"),
        ]
        anchor_idx = max(anchors)
        if anchor_idx == -1:
            anchor_idx = len(full_text) // 2

        mid_start = max(0, anchor_idx - mid_len // 2)
        mid_chunk = full_text[mid_start : mid_start + mid_len]

        truncated = (
            f"{full_text[:head_len]}"
            f"{marker}"
            f"{mid_header}"
            f"{mid_chunk}"
            f"{marker}"
            f"{tail_header}"
            f"{full_text[-tail_len:]}"
        )

        return truncated[:max_chars]

    def _estimate_table_rows(self, full_text: str) -> int:
        """
        Estimate how many rows of tabular data are present.

        This is used to detect when we're looking at a large variant table so
        we can be more aggressive about rerunning extraction with stronger
        models if the first attempt undercounts variants.
        """
        count = 0
        for line in full_text.splitlines():
            stripped = line.strip()
            # Count lines that look like real table rows (multiple columns)
            if stripped.startswith("|") and stripped.count("|") >= TABLE_MIN_COLUMNS:
                count += 1
        return count

    def _extract_variants_from_tables(
        self, full_text: str, gene_symbol: Optional[str]
    ) -> List[dict]:
        """
        Extract variants directly from markdown tables using regex patterns.

        This is a supplemental extraction pass that catches variants the LLM might miss,
        especially in large tables with many rows.

        Returns a list of minimal variant dicts.

        RE-ENABLED (2026-02-10): Now uses VariantNormalizer.is_non_target_variant()
        to filter out non-target gene variants (TP53, KRAS, BRAF hotspots etc.) and
        validates position against gene-specific protein length.
        """
        from utils.variant_normalizer import VariantNormalizer, normalize_variant

        if not gene_symbol:
            return []

        # Create normalizer for validation
        normalizer = VariantNormalizer(gene_symbol)
        target_gene = gene_symbol.upper()

        variants = []
        seen_variants = set()
        filtered_count = 0
        row_gene_filtered_count = 0

        # Pre-process: normalize Unicode arrows (same as variant_scanner)
        full_text = (
            full_text.replace("\u2192", ">")
            .replace("\u2190", "<")
            .replace("\u21d2", ">")
            .replace("\\_", "_")
            .replace("\\*", "*")
        )

        # Variant patterns to look for in table cells
        protein_pattern = re.compile(
            r"p\.\(?("
            r"(?:[A-Z][a-z]{2}|[A-Z])\d{1,4}"
            r"(?:_(?:[A-Z][a-z]{2}|[A-Z])\d{1,4})?"
            r"(?:"
            r"[A-Z][a-z]{2}|[A-Z*]|"
            r"fs(?:Ter|X|\*)?\d*|"
            r"del|dup|ins(?:[A-Z][a-z]{2}|[A-Z])*"
            r")"
            r")\)?",
            re.IGNORECASE,
        )
        cdna_pattern = re.compile(
            r"c\.(\d+[+-]?\d*[ACGT]>[ACGT]|\d+(?:_\d+)?(?:del|dup|ins)[ACGT]*|\d+[ACGT]>[ACGT])",
            re.IGNORECASE,
        )
        short_protein_pattern = re.compile(r"\b([A-Z]\d{2,4}(?:[A-Z*]|fs|del|dup))\b")

        def genes_mentioned_in_cells(cells: list[str]) -> set[str]:
            genes = set()
            for cell in cells:
                searchable = re.sub(r"[*_`\\]", " ", cell)
                searchable_upper = searchable.upper()
                for alias, official in self.TABLE_LABEL_GENE_ALIASES.items():
                    if re.search(
                        rf"(?<![A-Z0-9]){re.escape(alias)}(?![A-Z0-9])",
                        searchable_upper,
                    ):
                        genes.add(official)
                for gene in known_gene_aliases(include_query_aliases=False):
                    if gene_alias_regex(gene, include_query_aliases=False).search(
                        searchable
                    ):
                        genes.add(gene)
                for explicit_gene in self.TABLE_LABEL_EXPLICIT_GENES:
                    if re.search(
                        rf"(?<![A-Z0-9]){re.escape(explicit_gene)}(?![A-Z0-9])",
                        searchable_upper,
                    ):
                        genes.add(explicit_gene)
            return genes

        active_table_genes: set[str] = set()

        # Find all markdown table rows
        for line in full_text.splitlines():
            line = line.strip()
            if not line.startswith("|"):
                active_table_genes = set()
                continue
            if set(line) <= {"|", "-", " ", ":"}:  # Skip separator rows
                continue

            cells = [cell.strip() for cell in line.strip("|").split("|")]
            mentioned_genes = genes_mentioned_in_cells(cells)
            if mentioned_genes:
                active_table_genes = mentioned_genes

            row_genes = []
            for cell in cells:
                normalized_cell = re.sub(r"[^A-Za-z0-9]", "", cell).upper()
                row_gene = self.TABLE_LABEL_GENE_ALIASES.get(normalized_cell)
                if not row_gene and normalized_cell in self.TABLE_LABEL_EXPLICIT_GENES:
                    row_gene = normalized_cell
                if row_gene:
                    row_genes.append(row_gene)
                row_genes.extend(genes_mentioned_in_cells([cell]))
            context_genes = set(row_genes) or active_table_genes
            if context_genes and target_gene not in context_genes:
                row_gene_filtered_count += 1
                continue

            # Extract variants from this row
            # Look for protein notation
            for match in protein_pattern.finditer(line):
                notation = f"p.{match.group(1)}"
                normalized = normalize_variant(notation, gene_symbol)

                # Validate: skip non-target gene variants
                is_non_target, reason = normalizer.is_non_target_variant(normalized)
                if is_non_target:
                    filtered_count += 1
                    logger.debug(f"Table regex: filtered {notation} - {reason}")
                    continue

                if normalized not in seen_variants:
                    seen_variants.add(normalized)
                    variants.append(
                        {
                            "gene_symbol": gene_symbol,
                            "protein_notation": normalized,
                            "source_location": "Table (regex extraction)",
                        }
                    )

            # Look for cDNA notation
            for match in cdna_pattern.finditer(line):
                notation = f"c.{match.group(1)}"
                if notation not in seen_variants:
                    seen_variants.add(notation)
                    variants.append(
                        {
                            "gene_symbol": gene_symbol,
                            "cdna_notation": notation,
                            "source_location": "Table (regex extraction)",
                        }
                    )

            # Look for short protein notation (A561V)
            for match in short_protein_pattern.finditer(line):
                notation = match.group(1)
                # Skip if it looks like a figure/table reference (e.g., "S1", "T1")
                if len(notation) < 4:
                    continue

                normalized = normalize_variant(notation, gene_symbol)

                # Validate: skip non-target gene variants
                is_non_target, reason = normalizer.is_non_target_variant(normalized)
                if is_non_target:
                    filtered_count += 1
                    logger.debug(f"Table regex: filtered {notation} - {reason}")
                    continue

                if normalized not in seen_variants:
                    seen_variants.add(normalized)
                    variants.append(
                        {
                            "gene_symbol": gene_symbol,
                            "protein_notation": normalized,
                            "source_location": "Table (regex extraction)",
                        }
                    )

        if variants:
            logger.info(
                f"Table regex extraction found {len(variants)} variants "
                f"(filtered {filtered_count} non-target, "
                f"{row_gene_filtered_count} off-target gene rows)"
            )
        return variants

    def _merge_table_variants(
        self, extracted_data: dict, table_variants: List[dict]
    ) -> dict:
        """
        Merge variants found via table regex extraction with LLM-extracted variants.

        Avoids duplicates by checking normalized notation.
        """
        from utils.variant_normalizer import normalize_variant

        existing_variants = extracted_data.get("variants", [])

        # Build set of existing variant keys
        existing_keys = set()
        for v in existing_variants:
            protein = normalize_variant(v.get("protein_notation", "") or "")
            cdna = normalize_variant(v.get("cdna_notation", "") or "")
            if protein:
                existing_keys.add(protein)
            if cdna:
                existing_keys.add(cdna)

        # Add new variants not already present
        added_count = 0
        for tv in table_variants:
            protein = normalize_variant(tv.get("protein_notation", "") or "")
            cdna = normalize_variant(tv.get("cdna_notation", "") or "")

            # Check if already exists
            if protein and protein in existing_keys:
                continue
            if cdna and cdna in existing_keys:
                continue

            # Add to existing variants. Router-produced table variants already
            # carry patient / penetrance fields, while regex-only hints are
            # intentionally sparse. Preserve rich fields when present and fill
            # only the missing schema defaults.
            new_variant = {
                "gene_symbol": tv.get("gene_symbol"),
                "cdna_notation": tv.get("cdna_notation"),
                "protein_notation": tv.get("protein_notation"),
                "clinical_significance": "unknown",
                "evidence_level": "low",
                "source_location": tv.get("source_location", "Table"),
                "additional_notes": "Added via table regex extraction",
                "patients": {},
                "penetrance_data": {},
                "individual_records": [],
                "functional_data": {"summary": "", "assays": []},
                "key_quotes": [],
            }
            for key, value in tv.items():
                if value is not None:
                    new_variant[key] = value
            new_variant.setdefault("patients", {})
            new_variant.setdefault("penetrance_data", {})
            new_variant.setdefault("individual_records", [])
            new_variant.setdefault("functional_data", {"summary": "", "assays": []})
            new_variant.setdefault("key_quotes", [])
            existing_variants.append(new_variant)
            added_count += 1

            # Update keys
            if protein:
                existing_keys.add(protein)
            if cdna:
                existing_keys.add(cdna)

        if added_count > 0:
            logger.info(
                f"Merged {added_count} additional variants from table extraction"
            )
            # Update metadata
            if "extraction_metadata" in extracted_data:
                extracted_data["extraction_metadata"]["total_variants_found"] = len(
                    existing_variants
                )
                extracted_data["extraction_metadata"]["table_regex_added"] = added_count

        extracted_data["variants"] = existing_variants
        return extracted_data

    def _table_variant_key(self, variant: dict) -> tuple[str, str]:
        """Return a normalized dedupe key for a table-derived variant."""
        from utils.variant_normalizer import normalize_variant

        protein = normalize_variant(variant.get("protein_notation", "") or "")
        cdna = normalize_variant(variant.get("cdna_notation", "") or "")
        return (cdna.lower(), protein.lower())

    def _is_structured_table_variant(self, variant: dict) -> bool:
        """True when a table variant carries row/count provenance, not just regex text."""
        source_location = str(variant.get("source_location") or "").lower()
        patients = variant.get("patients") or {}
        penetrance = variant.get("penetrance_data") or {}
        return (
            "router+deterministic" in source_location
            or bool(variant.get("fact_provenance"))
            or bool(variant.get("count_provenance"))
            or patients.get("count") is not None
            or penetrance.get("total_carriers_observed") is not None
            or penetrance.get("affected_count") is not None
        )

    def _dedupe_table_variants(self, table_variants: List[dict]) -> List[dict]:
        """Dedupe table candidates while preferring richer row-level variants."""
        selected: List[dict] = []
        index_by_key: dict[tuple[str, str], int] = {}
        for variant in table_variants:
            key = self._table_variant_key(variant)
            if not any(key):
                continue
            if key not in index_by_key:
                index_by_key[key] = len(selected)
                selected.append(variant)
                continue
            existing_index = index_by_key[key]
            existing = selected[existing_index]
            if self._is_structured_table_variant(
                variant
            ) and not self._is_structured_table_variant(existing):
                selected[existing_index] = variant
        return selected

    def _select_table_variants_for_merge(
        self, table_variants: List[dict]
    ) -> tuple[List[dict], Optional[dict]]:
        """Select candidates for table merge and describe overflow decisions."""
        deduped = self._dedupe_table_variants(table_variants)
        if len(table_variants) <= TABLE_REGEX_MERGE_MAX_VARIANTS:
            return deduped, None

        overflow_max = max(
            TABLE_REGEX_MERGE_MAX_VARIANTS, TABLE_REGEX_OVERFLOW_MERGE_MAX_VARIANTS
        )
        structured = [v for v in deduped if self._is_structured_table_variant(v)]
        regex_only = [v for v in deduped if not self._is_structured_table_variant(v)]

        selected: List[dict] = []
        selected.extend(structured[:overflow_max])
        remaining_slots = max(0, overflow_max - len(selected))
        if remaining_slots:
            selected.extend(regex_only[:remaining_slots])

        overflow = {
            "candidate_count": len(table_variants),
            "deduped_count": len(deduped),
            "normal_safety_cap": TABLE_REGEX_MERGE_MAX_VARIANTS,
            "overflow_merge_cap": overflow_max,
            "chunk_size": max(1, TABLE_REGEX_OVERFLOW_CHUNK_SIZE),
            "structured_candidate_count": len(structured),
            "regex_only_candidate_count": len(regex_only),
            "selected_for_merge": len(selected),
            "omitted_after_dedupe": max(0, len(deduped) - len(selected)),
            "reason": "candidate_count_exceeds_safety_cap",
        }
        return selected, overflow

    def _merge_table_variants_with_overflow_qc(
        self, extracted_data: dict, table_variants: List[dict], pmid: Optional[str]
    ) -> dict:
        """Merge table candidates; when over cap, merge in bounded chunks with QC."""
        selected, overflow = self._select_table_variants_for_merge(table_variants)
        if overflow:
            logger.warning(
                "PMID %s - Table regex/router merge overflow: %d candidates "
                "(deduped=%d), merging %d with overflow cap %d",
                pmid,
                overflow["candidate_count"],
                overflow["deduped_count"],
                overflow["selected_for_merge"],
                overflow["overflow_merge_cap"],
            )

        before_count = len(extracted_data.get("variants", []) or [])
        if selected:
            chunk_size = (
                max(1, TABLE_REGEX_OVERFLOW_CHUNK_SIZE) if overflow else len(selected)
            )
            for start in range(0, len(selected), chunk_size):
                extracted_data = self._merge_table_variants(
                    extracted_data, selected[start : start + chunk_size]
                )

        after_count = len(extracted_data.get("variants", []) or [])
        metadata = extracted_data.setdefault("extraction_metadata", {})
        metadata["table_merge_candidate_count"] = len(table_variants)
        if overflow:
            overflow["merged_added"] = max(0, after_count - before_count)
            metadata["table_merge_overflow"] = overflow
            metadata["table_regex_added"] = overflow["merged_added"]
        return extracted_data

    def _format_table_hints(
        self, table_variants: List[dict], max_hints: int = TABLE_HINT_MAX_VARIANTS
    ) -> str:
        """
        Format pre-extracted table variants as structured hints for the LLM prompt.

        This gives the LLM a head start by showing variants already detected via regex,
        helping it focus on extracting clinical details and catching any missed variants.
        """
        if not table_variants:
            return ""

        hints = [
            "\n\n--- PRE-EXTRACTED TABLE VARIANTS (regex detection) ---",
            f"The following {len(table_variants)} variant(s) were detected in tables via pattern matching.",
            "Use these as hints - verify and enrich with clinical details from the text:",
            "",
        ]

        shown_variants = table_variants[:max_hints]
        omitted = max(0, len(table_variants) - len(shown_variants))

        for i, v in enumerate(shown_variants, 1):
            parts = [f"{i}."]
            if v.get("cdna_notation"):
                parts.append(v["cdna_notation"])
            if v.get("protein_notation"):
                parts.append(f"/ {v['protein_notation']}")
            if v.get("source_location"):
                parts.append(f"[{v['source_location']}]")
            hints.append(" ".join(parts))

        if omitted:
            hints.append(
                f"... {omitted} additional table-detected variants omitted from prompt hints."
            )

        hints.append("\n--- END PRE-EXTRACTED HINTS ---\n")
        return "\n".join(hints)

    def _filter_by_gene(self, extracted_data: dict, target_gene: str) -> dict:
        """
        Filter variants to only keep those matching the target gene.

        This removes variants that were extracted for other genes (common in papers
        that study multiple genes). Matching is case-insensitive.

        Args:
            extracted_data: The extraction result dictionary
            target_gene: The target gene symbol to filter for (e.g., "KCNH2")

        Returns:
            Updated extracted_data with filtered variants
        """
        variants = extracted_data.get("variants", [])
        if not variants:
            return extracted_data

        target_upper = target_gene.upper()
        original_count = len(variants)

        filtered_variants = []
        removed_genes = set()

        for v in variants:
            gene = v.get("gene_symbol", "") or ""
            gene_upper = gene.upper()

            # Keep if gene matches target or is empty/unknown
            if not gene or gene_upper == target_upper or gene_upper in ("UNKNOWN", ""):
                filtered_variants.append(v)
            else:
                removed_genes.add(gene)

        removed_count = original_count - len(filtered_variants)

        if removed_count > 0:
            logger.info(
                f"Filtered out {removed_count} variants from non-target genes: {removed_genes}"
            )
            # Update metadata
            if "extraction_metadata" in extracted_data:
                extracted_data["extraction_metadata"]["total_variants_found"] = len(
                    filtered_variants
                )
                extracted_data["extraction_metadata"]["filtered_out_by_gene"] = (
                    removed_count
                )
                extracted_data["extraction_metadata"]["filtered_genes"] = list(
                    removed_genes
                )

        extracted_data["variants"] = filtered_variants
        return extracted_data

    # Artifact patterns that indicate failed/invalid extraction
    ARTIFACT_PATTERNS = [
        r"^p\.XXX$",  # Placeholder notation
        r"^p\.unknown$",  # Unknown notation
        r"^p\.null$",  # Null notation
        r"^null$",  # Raw null
        r"^unknown$",  # Raw unknown
        r"^Splicesite$",  # Generic splice label (not actual notation)
        r"^splice$",  # Generic splice label
        r"^N/A$",  # Not applicable
        r"^NA$",  # Not applicable variant
        r"^-$",  # Dash placeholder
        r"^\?$",  # Question mark placeholder
    ]
    AA3_RE = (
        r"Ala|Cys|Asp|Glu|Phe|Gly|His|Ile|Lys|Leu|Met|Asn|Pro|Gln|"
        r"Arg|Ser|Thr|Val|Trp|Tyr|Ter|Stop|Xaa"
    )
    PROTEIN_NOTATION_RE = re.compile(
        rf"^(?:p\.)?(?:{AA3_RE}|[ACDEFGHIKLMNPQRSTVWY])"
        r"\d{1,4}"
        rf"(?:[_-](?:{AA3_RE}|[ACDEFGHIKLMNPQRSTVWY])\d{{1,4}})?"
        rf"(?:{AA3_RE}|[ACDEFGHIKLMNPQRSTVWY*X?]|fs(?:X|\*)?\d*|del|dup|ins)",
        re.IGNORECASE,
    )
    CDNA_NOTATION_RE = re.compile(
        r"^c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
        r"|^c\.\d+(?:[+-]\d+)?(?:del|dup|ins)[ACGT]*$"
        r"|^c\.\d+(?:_\d+)?(?:del|dup|ins)[ACGT]*$",
        re.IGNORECASE,
    )

    def _filter_extraction_artifacts(
        self, extracted_data: dict, target_gene: str
    ) -> dict:
        """
        Filter out extraction artifacts and invalid variants.

        Removes:
        - Placeholder notations: p.XXX, p.unknown, null, Splicesite, etc.
        - Variants with amino acid position exceeding protein length
          (e.g., position > 1159 for KCNH2)

        This is a post-extraction cleanup step that removes noise introduced
        by LLM extraction errors or malformed source data.

        Args:
            extracted_data: The extraction result dictionary
            target_gene: The target gene symbol for position validation

        Returns:
            Updated extracted_data with artifacts removed
        """
        from utils.variant_normalizer import PROTEIN_LENGTHS, VariantNormalizer

        variants = extracted_data.get("variants", [])
        if not variants:
            return extracted_data

        # Compile artifact patterns
        artifact_re = re.compile("|".join(self.ARTIFACT_PATTERNS), re.IGNORECASE)

        # Get protein length for position validation
        protein_length = PROTEIN_LENGTHS.get(target_gene.upper())
        normalizer = VariantNormalizer(target_gene)

        original_count = len(variants)
        filtered_variants = []
        artifact_count = 0
        malformed_count = 0
        position_invalid_count = 0
        artifact_examples = []
        malformed_examples = []
        position_examples = []

        for v in variants:
            protein = v.get("protein_notation", "") or ""
            cdna = v.get("cdna_notation", "") or ""

            # Check for artifact patterns
            is_artifact = False
            if protein and artifact_re.match(protein.strip()):
                is_artifact = True
                if len(artifact_examples) < 5:
                    artifact_examples.append(protein)
            if cdna and artifact_re.match(cdna.strip()):
                is_artifact = True
                if len(artifact_examples) < 5:
                    artifact_examples.append(cdna)

            if is_artifact:
                artifact_count += 1
                continue

            # Malformed notation guard. Protein entries need an amino-acid
            # position; cDNA entries need HGVS-like c. notation. This catches
            # table-parser artifacts such as A/T/G/C, p-values, allele
            # frequencies, and header fragments before they reach SQLite.
            protein_clean = protein.strip().replace(" ", "")
            cdna_clean = cdna.strip().replace(" ", "")
            malformed = False
            if protein_clean and not self.PROTEIN_NOTATION_RE.match(protein_clean):
                malformed = True
                if len(malformed_examples) < 5:
                    malformed_examples.append(protein)
            if cdna_clean and not self.CDNA_NOTATION_RE.match(cdna_clean):
                malformed = True
                if len(malformed_examples) < 5:
                    malformed_examples.append(cdna)
            if malformed:
                malformed_count += 1
                continue

            # Check position validity (for protein variants)
            if protein and protein_length:
                position = normalizer.extract_position(protein)
                if position and position > protein_length:
                    position_invalid_count += 1
                    if len(position_examples) < 5:
                        position_examples.append(f"{protein} (pos={position})")
                    continue

            filtered_variants.append(v)

        removed_count = original_count - len(filtered_variants)

        if removed_count > 0:
            logger.info(
                f"Filtered out {removed_count} artifacts: "
                f"{artifact_count} artifact patterns, "
                f"{malformed_count} malformed notations, "
                f"{position_invalid_count} invalid positions"
            )
            if artifact_examples:
                logger.debug(f"Artifact examples: {artifact_examples}")
            if malformed_examples:
                logger.debug(f"Malformed examples: {malformed_examples}")
            if position_examples:
                logger.debug(f"Position invalid examples: {position_examples}")

            # Update metadata
            if "extraction_metadata" in extracted_data:
                extracted_data["extraction_metadata"]["total_variants_found"] = len(
                    filtered_variants
                )
                extracted_data["extraction_metadata"]["artifacts_filtered"] = (
                    artifact_count
                )
                extracted_data["extraction_metadata"]["malformed_filtered"] = (
                    malformed_count
                )
                extracted_data["extraction_metadata"]["position_invalid_filtered"] = (
                    position_invalid_count
                )

        extracted_data["variants"] = filtered_variants
        return extracted_data

    # Header patterns for variant table detection - broadened to catch more table formats
    # Groups: cdna_headers, protein_headers, count_headers
    VARIANT_CDNA_HEADERS = {
        "nucleotide",
        "nucleotide change",
        "cds",
        "cdna",
        "c.",
        "cdna change",
        "cdna_change",
        "hgvs cdna",
        "hgvs_cdna",
        "dna change",
        "dna mutation",
        "coding change",
        "nt change",
        "base change",
    }
    VARIANT_PROTEIN_HEADERS = {
        "variant",
        "amino acid",
        "protein",
        "p.",
        "aa change",
        "aachange",
        "protein change",
        "mutation",
        "hgvs protein",
        "hgvs_protein",
        "amino acid change",
        "coding effect",
        "protein mutation",
        "missense",
    }
    VARIANT_GENE_HEADERS = {
        "gene",
        "gene symbol",
        "genesymbol",
        "symbol",
    }
    VARIANT_COUNT_HEADERS = {
        "patient",
        "patients",
        "no. of patients",
        "no. of patient",
        "carrier",
        "carriers",
        "probands",
        "families",
        "subjects",
        "count",
        "cases",
        "individuals",
        "number",
    }
    VARIANT_AFFECTED_HEADERS = {
        "affected",
        "affected carriers",
        "affected_count",
        "naffected",
        "n affected",
        "no. affected",
        "symptomatic",
        "symptomatic carriers",
        "cases",
    }
    VARIANT_UNAFFECTED_HEADERS = {
        "unaffected",
        "unaffected carriers",
        "unaffected_count",
        "nunaffected",
        "n unaffected",
        "no. unaffected",
        "asymptomatic",
        "asymptomatic carriers",
        "controls",
    }
    VARIANT_ROW_LEVEL_HEADERS = {
        "patient",
        "patient id",
        "proband",
        "subject",
        "case",
        "family",
        "kindred",
        "individual",
        "participant",
        "pedigree",
    }
    VARIANT_CLINICAL_CONTEXT_HEADERS = {
        "age",
        "sex",
        "gender",
        "phenotype",
        "diagnosis",
        "clinical",
        "symptom",
        "syncope",
        "qtc",
        "qt interval",
        "ecg",
        "arrhythmia",
        "onset",
        "sudden death",
        "therapy",
        "treatment",
        "schwartz score",
    }
    GWAS_ASSOCIATION_HEADERS = {
        "locus",
        "snv",
        "snp",
        "rsid",
        "rs id",
        "chr",
        "chromosome",
        "bp",
        "position",
        "ea",
        "effect allele",
        "aa",
        "alternate allele",
        "eaf",
        "maf",
        "beta",
        "se",
        "p",
        "p value",
        "p-value",
        "or",
        "odds ratio",
    }
    TABLE_LABEL_GENE_RE = re.compile(r"\b[A-Z][A-Z0-9]{1,11}\b")
    TABLE_LABEL_GENE_STOPWORDS = {
        "AA",
        "AF",
        "ALL",
        "AND",
        "BASIC",
        "BRS",
        "CASES",
        "CANCER",
        "CARRIER",
        "CARRIERS",
        "CDNA",
        "CHR",
        "CLINICAL",
        "CONTROL",
        "CONTROLS",
        "CPVT",
        "DNA",
        "EFFECT",
        "ETABLE",
        "ETABLES",
        "EXON",
        "EXONS",
        "FIG",
        "FIGURE",
        "FUNCTIONAL",
        "GENE",
        "GENES",
        "HGVS",
        "HUMAN",
        "HEREDITARY",
        "IDENTIFIED",
        "IN",
        "INDEL",
        "INCLUDED",
        "LIST",
        "LQTS",
        "MAF",
        "MISSENSE",
        "MUTATION",
        "MUTATIONS",
        "NO",
        "NON",
        "NONSYNONYMOUS",
        "NSSNV",
        "NSSNVS",
        "OF",
        "OBSERVED",
        "PATIENT",
        "PATIENTS",
        "PMID",
        "PROBAND",
        "PROBANDS",
        "PROTEIN",
        "RARE",
        "RNA",
        "SNV",
        "SNVS",
        "SNP",
        "SNPS",
        "SUPP",
        "SUPPLEMENT",
        "SUPPLEMENTAL",
        "SUPPLEMENTARY",
        "TABLE",
        "THE",
        "VARIANT",
        "VARIANTS",
        "WES",
        "WT",
    }
    TABLE_LABEL_GENE_ALIASES = {
        "LQT1": "KCNQ1",
        "LQT2": "KCNH2",
        "LQT3": "SCN5A",
    }
    TABLE_LABEL_CONTEXTUAL_GENE_RE = re.compile(
        r"\b([A-Z][A-Z0-9-]{2,11})\s+"
        r"(?:gene|genes|variant|variants|mutation|mutations)\b"
        r"|"
        r"\b(?:gene|genes|variant|variants|mutation|mutations)\s+"
        r"(?:(?:in|of|for|from)\s+|:\s*)?"
        r"([A-Z][A-Z0-9-]{2,11})\b"
        r"|"
        r"\b([A-Z][A-Z0-9-]{2,11})(?=-(?:associated|related)\b)",
        re.IGNORECASE,
    )
    TABLE_LABEL_EXPLICIT_GENES = {
        "ANK2",
        "CACNA1C",
        "CALM1",
        "CALM2",
        "CALM3",
        "CASQ2",
        "CAV3",
        "KCNE1",
        "KCNE2",
        "KCNH2",
        "KCNJ2",
        "KCNQ1",
        "RYR2",
        "SCN5A",
        "SNTA1",
        "TECRL",
        "TRDN",
    }
    LQTS_COMPENDIUM_MUTATION_TYPES = (
        "Inframe insertion",
        "Inframe deletion",
        "Splice site",
        "Frameshift",
        "Missense",
        "Nonsense",
    )

    PROTEIN_TABLE_VALUE_RE = re.compile(
        r"^(?:p\.)?(?:"
        r"(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
        r"\d{1,4}"
        r"(?:[_-](?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])\d{1,4})?"
        r"(?:"
        r"(?:(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY]))?fs(?:Ter|X|\*)?\d*"
        r"|(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
        r"|[ACDEFGHIKLMNPQRSTVWY*X?]"
        r"|del(?:\d*|[ACDEFGHIKLMNPQRSTVWY]+)"
        r"|dup"
        r"|ins(?:(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])+)?"
        r")"
        r"|\d{1,4}_\d{1,4}ins(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY]+)"
        r")$",
        re.IGNORECASE,
    )
    CDNA_TABLE_VALUE_RE = re.compile(
        r"^c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
        r"|^c\.\d+(?:[+-]\d+)?(?:del|dup|ins)[ACGT]*$"
        r"|^c\.\d+(?:_\d+)?(?:del|dup|ins)[ACGT]*$",
        re.IGNORECASE,
    )

    def _header_matches(self, header: str, candidates: set[str]) -> bool:
        """Match table headers without broad substring hits like `n` in `Unnamed`."""
        normalized = re.sub(r"[^a-z0-9.+/-]+", " ", header.lower()).strip()
        compact = normalized.replace(" ", "")
        tokens = set(normalized.split())

        for candidate in candidates:
            cand = candidate.lower()
            cand_compact = re.sub(r"[^a-z0-9.+/-]+", "", cand)

            # Very short headers are dangerous as substrings. Require exact
            # token/cell matches so GWAS columns like AA/n do not become
            # protein/count columns.
            if len(cand_compact) <= 2:
                if normalized == cand or compact == cand_compact or cand in tokens:
                    return True
                continue

            if cand in normalized or cand_compact in compact:
                return True

        return False

    def _looks_like_gwas_header(self, cells: list[str]) -> bool:
        """Identify association/GWAS summary tables, not clinical carrier tables."""
        normalized_cells = [
            re.sub(r"[^a-z0-9.+/-]+", " ", c.lower()).strip() for c in cells
        ]
        marker_count = 0
        for cell in normalized_cells:
            compact = cell.replace(" ", "")
            if any(
                cell == marker or compact == marker.replace(" ", "")
                for marker in self.GWAS_ASSOCIATION_HEADERS
            ):
                marker_count += 1

        clinical_count = sum(
            1
            for cell in normalized_cells
            if any(
                term in cell
                for term in (
                    "patient",
                    "carrier",
                    "affected",
                    "unaffected",
                    "proband",
                    "case",
                    "family",
                )
            )
        )
        return marker_count >= 4 and clinical_count == 0

    def _looks_like_row_level_clinical_header(self, cells: list[str]) -> bool:
        """Detect patient/proband-level mutation lists without count columns."""
        normalized_cells = [
            re.sub(r"[^a-z0-9.+/-]+", " ", c.lower()).strip() for c in cells
        ]
        has_row_subject = any(
            any(term in cell for term in self.VARIANT_ROW_LEVEL_HEADERS)
            for cell in normalized_cells
        )
        has_clinical_context = any(
            any(term in cell for term in self.VARIANT_CLINICAL_CONTEXT_HEADERS)
            for cell in normalized_cells
        )
        return has_row_subject or has_clinical_context

    def _clean_table_cell(self, value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        cleaned = (
            value.strip()
            .replace(" ", "")
            .replace("\u00a0", "")
            .replace("\u2009", "")
            .replace("\u200a", "")
            .replace("→", ">")
            .replace("–", "-")
            .replace("−", "-")
            .rstrip(",;")
        )
        if not cleaned or cleaned.lower() in {"nan", "na", "n/a", "none", "-", "."}:
            return None
        return cleaned

    def _valid_table_protein(self, value: Optional[str]) -> Optional[str]:
        cleaned = self._clean_table_cell(value)
        if not cleaned:
            return None
        if (
            cleaned.lower().startswith("p")
            and not cleaned.lower().startswith("p.")
            and len(cleaned) > 1
            and cleaned[1].isalpha()
        ):
            cleaned = f"p.{cleaned[1:]}"
        deletion_prefix = re.match(
            r"^del(\d{1,4})([ACDEFGHIKLMNPQRSTVWY])$", cleaned, re.IGNORECASE
        )
        if deletion_prefix:
            cleaned = f"{deletion_prefix.group(2).upper()}{deletion_prefix.group(1)}del"
        deletion_infix = re.match(
            r"^(?:p\.)?(\d{1,4})del([ACDEFGHIKLMNPQRSTVWY])$",
            cleaned,
            re.IGNORECASE,
        )
        if deletion_infix:
            prefix = "p." if cleaned.lower().startswith("p.") else ""
            cleaned = (
                f"{prefix}{deletion_infix.group(2).upper()}{deletion_infix.group(1)}del"
            )
        cleaned = cleaned.rstrip("⁎†‡§")
        residue_range_deletion = re.match(
            r"^(p\.)?([ACDEFGHIKLMNPQRSTVWY]{2,})(\d{1,4})-(\d{1,4})del$",
            cleaned,
            re.IGNORECASE,
        )
        if residue_range_deletion:
            residues = residue_range_deletion.group(2).upper()
            cleaned = (
                f"{residue_range_deletion.group(1) or ''}"
                f"{residues[0]}{residue_range_deletion.group(3)}_"
                f"{residues[-1]}{residue_range_deletion.group(4)}del"
            )
        ocr_terminal_i = re.match(
            r"^(p\.)?([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})1$",
            cleaned,
            re.IGNORECASE,
        )
        if ocr_terminal_i:
            cleaned = (
                f"{ocr_terminal_i.group(1) or ''}"
                f"{ocr_terminal_i.group(2).upper()}{ocr_terminal_i.group(3)}I"
            )
        frameshift_slash = re.match(
            r"^(p\.)?([ACDEFGHIKLMNPQRSTVWY]\d+)fs/\d+$",
            cleaned,
            re.IGNORECASE,
        )
        if frameshift_slash:
            cleaned = f"{frameshift_slash.group(1) or ''}{frameshift_slash.group(2)}fsX"
        splice = re.match(
            r"^(p\.)?([ACDEFGHIKLMNPQRSTVWY])(\d+)/?sp$",
            cleaned,
            re.IGNORECASE,
        )
        if splice:
            cleaned = (
                f"{splice.group(1) or ''}{splice.group(2).upper()}{splice.group(3)}sp"
            )
            return cleaned
        if cleaned.endswith("*"):
            without_footnote = cleaned[:-1]
            if without_footnote and self.PROTEIN_TABLE_VALUE_RE.match(without_footnote):
                cleaned = without_footnote
        if self.PROTEIN_TABLE_VALUE_RE.match(cleaned):
            return cleaned
        return None

    def _valid_table_cdna(self, value: Optional[str]) -> Optional[str]:
        cleaned = self._clean_table_cell(value)
        if not cleaned:
            return None
        if not cleaned.lower().startswith("c."):
            if cleaned.lower().startswith("c"):
                cleaned = "c." + cleaned[1:]
            else:
                cleaned = "c." + cleaned
        body = cleaned[2:].upper()
        body = re.sub(
            r"(DEL|DUP|INS)([ACGT]*)",
            lambda match: match.group(1).lower() + match.group(2),
            body,
        )
        cleaned = "c." + body
        if self.CDNA_TABLE_VALUE_RE.match(cleaned):
            return cleaned
        return None

    def _is_variant_table_header(self, line: str) -> bool:
        """
        Check if a line looks like a variant table header row.

        Requires at least one cDNA/protein-like header AND one count-like header,
        OR two different variant notation headers (cDNA + protein).
        """
        if "|" not in line:
            return False

        parts = [p.strip() for p in line.split("|") if p.strip()]
        if self._looks_like_gwas_header(parts):
            return False

        normalized_parts = [
            re.sub(r"[^a-z0-9.+/-]+", " ", p.lower()).strip() for p in parts
        ]
        unnamed_count = sum(p.startswith("unnamed") for p in normalized_parts)
        if unnamed_count and (
            unnamed_count >= max(2, len(parts) // 2)
            or any(p.startswith("table ") for p in normalized_parts)
        ):
            return False

        has_cdna_header = any(
            self._header_matches(part, self.VARIANT_CDNA_HEADERS) for part in parts
        )
        has_protein_header = any(
            self._header_matches(part, self.VARIANT_PROTEIN_HEADERS) for part in parts
        )
        has_count_header = any(
            self._header_matches(part, self.VARIANT_COUNT_HEADERS) for part in parts
        )
        has_unaffected_header = any(
            self._header_matches(part, self.VARIANT_UNAFFECTED_HEADERS)
            for part in parts
        )
        has_affected_header = any(
            self._header_matches(part, self.VARIANT_AFFECTED_HEADERS)
            and not self._header_matches(part, self.VARIANT_UNAFFECTED_HEADERS)
            for part in parts
        )
        has_clinical_count_header = (
            has_count_header or has_affected_header or has_unaffected_header
        )
        has_row_level_clinical_header = self._looks_like_row_level_clinical_header(
            parts
        )

        # Accept: (cDNA OR protein) AND count
        # OR: row-level clinical mutation list where each row is one carrier
        # OR: cDNA AND protein (variant mapping table)
        return (has_cdna_header or has_protein_header) and (
            has_clinical_count_header
            or has_row_level_clinical_header
            or (has_cdna_header and has_protein_header)
        )

    def _gene_symbols_in_table_label(self, label: str) -> set[str]:
        """Return likely gene symbols mentioned in a table title/caption."""
        if not label:
            return set()

        symbols: set[str] = set()
        for match in self.TABLE_LABEL_GENE_RE.finditer(label.upper()):
            token = match.group(0)
            if token in self.TABLE_LABEL_GENE_STOPWORDS:
                continue
            if re.fullmatch(r"[A-Z]\d+", token):
                continue
            if len(token) <= 3 and token not in {
                "TTN",
                "APC",
                "RET",
                "VHL",
                "NF1",
                "NF2",
            }:
                continue
            if token.startswith("LQT") or token.startswith("SQT"):
                continue
            symbols.add(token)
        return symbols

    def _gene_scope_from_table_label(self, label: str) -> set[str]:
        """Return explicit target genes implied by a table title."""
        symbols = {
            symbol
            for symbol in self._gene_symbols_in_table_label(label)
            if symbol in self.TABLE_LABEL_EXPLICIT_GENES
        }
        label_upper = label.upper()
        for alias, gene in self.TABLE_LABEL_GENE_ALIASES.items():
            if re.search(rf"\b{re.escape(alias)}\b", label_upper):
                symbols.add(gene)
        for gene in known_gene_aliases(include_query_aliases=False):
            if gene_alias_regex(gene, include_query_aliases=False).search(label):
                symbols.add(gene)
        return symbols

    # Rare all-alpha gene symbols worth scoping even without a digit. Mirrors the
    # short-token whitelist in `_gene_symbols_in_table_label` so caption prose
    # (COMPENDIUM/PROPERTIES/TESTING/...) is not mistaken for a gene.
    _OPEN_VOCAB_ALPHA_GENES = {"TTN", "APC", "RET", "VHL", "NF1", "NF2"}

    def _contextual_gene_symbols_in_table_label(self, label: str) -> set[str]:
        """Return all-letter gene tokens when the caption labels them as genes."""
        if not label:
            return set()

        from pipeline.table_router import _gene_symbol_tokens

        symbols: set[str] = set()
        for match in self.TABLE_LABEL_CONTEXTUAL_GENE_RE.finditer(label.upper()):
            for raw_token in match.groups():
                token = (raw_token or "").strip()
                if (
                    not token
                    or token in self.TABLE_LABEL_GENE_STOPWORDS
                    or token.startswith(("LQT", "SQT"))
                    or token not in _gene_symbol_tokens(token)
                ):
                    continue
                if token.isalpha() and len(token) > 6:
                    continue
                symbols.add(token)
        return symbols

    def _fixed_width_gene_scope_from_table_label(self, label: str) -> set[str]:
        """Open-vocabulary gene scope for the fixed-width off-target guard.

        The markdown parser scopes arbitrary genes via
        `_gene_symbols_in_table_label`, but that helper is too noisy for the
        fixed-width guard's `bool(scope) and target not in scope` logic: it keeps
        ALL-CAPS non-genes (ETABLE/COMPENDIUM/PROPERTIES/...) that would make
        scope non-empty for nearly every caption and over-suppress claimable
        tables. Compose a tighter scope instead:

          explicit cardiac genes + LQT alias resolution
          ∪ confidently gene-shaped tokens from the careful `_gene_symbol_tokens`
            extractor (it excludes protein-/cDNA-like tokens), restricted to
            tokens that carry a digit, are a known all-alpha gene symbol, or are
            all-alpha tokens in a caption context that explicitly labels them as
            genes/variants/mutations.

        `_gene_symbol_tokens` ignores bare LQT/LQTS but still emits LQT1/LQT2/...,
        so drop LQT*/SQT* here and let the alias map resolve them to a gene.
        """
        from pipeline.table_router import _gene_symbol_tokens

        scope = set(self._gene_scope_from_table_label(label))
        for token in _gene_symbol_tokens(label):
            if (
                token in self.TABLE_LABEL_GENE_STOPWORDS
                or token.startswith("LQT")
                or token.startswith("SQT")
            ):
                continue
            if (
                any(ch.isdigit() for ch in token)
                or token in self._OPEN_VOCAB_ALPHA_GENES
            ):
                scope.add(token)
        scope.update(self._contextual_gene_symbols_in_table_label(label))
        return scope

    def _parse_lqt_fixed_width_variant_row(
        self,
        stripped: str,
        gene_symbol: str,
        current_table_label: str,
    ) -> Optional[dict]:
        """Parse LQT registry supplement rows with cDNA/protein/count columns."""
        if not re.search(r"(?:^|[\s*])(?:[cｃ]\.?\s*\d|p\.)", stripped):
            return None

        cdna = None
        cdna_match = re.search(
            r"(?:[cｃ]\.?\s*\d+(?:[+-]\d+|[-_]\d+)?\s*"
            r"(?:[ACGTacgt]\s*>\s*[ACGTacgt]|del\s*[ACGTacgt]*|"
            r"dup\s*[ACGTacgt]*|ins\s*[ACGTacgt]*))",
            stripped,
        )
        if cdna_match:
            cdna_raw = (
                cdna_match.group(0)
                .replace("ｃ", "c")
                .replace(" ", "")
                .replace("\u3000", "")
            )
            cdna = self._valid_table_cdna(cdna_raw)

        protein = None
        protein_match = re.search(r"\bp\.\s*", stripped, re.IGNORECASE)
        if protein_match:
            protein_tokens: list[str] = []
            for token in stripped[protein_match.end() :].split():
                clean_token = token.strip(",;")
                if (
                    not clean_token
                    or clean_token.startswith("+")
                    or clean_token.isdigit()
                    or clean_token.upper()
                    in {"N", "C", "MS", "TM", "DI", "DII", "DIII", "DIV"}
                    or clean_token.lower() in {"c-loop", "s5", "s6", "diii-div"}
                    or "term" in clean_token.lower()
                    or "pore" in clean_token.lower()
                    or clean_token.lower().startswith("splicing")
                ):
                    break
                protein_tokens.append(clean_token)
            if protein_tokens:
                protein_raw = "p." + "".join(protein_tokens)
                protein = self._valid_table_protein(protein_raw)

        if not cdna and not protein:
            return None

        numeric_tokens = [int(token) for token in stripped.split() if token.isdigit()]
        if len(numeric_tokens) < 3:
            return None
        count = numeric_tokens[0]
        affected_count = numeric_tokens[-2]
        unaffected_count = max(count - affected_count, 0)

        return {
            "gene_symbol": gene_symbol,
            "cdna_notation": cdna,
            "protein_notation": protein,
            "clinical_significance": "pathogenic",
            "patients": {
                "count": count,
                "phenotype": f"{gene_symbol}-associated disease",
            },
            "penetrance_data": {
                "total_carriers_observed": count,
                "affected_count": affected_count,
                "unaffected_count": unaffected_count,
            },
            "individual_records": [],
            "functional_data": {"summary": "", "assays": []},
            "segregation_data": None,
            "population_frequency": None,
            "evidence_level": "medium",
            "source_location": current_table_label or "Fixed-width LQT mutation table",
            "additional_notes": "Parsed via deterministic fixed-width LQT mutation table parser",
            "key_quotes": [],
        }

    def _normalize_lqts_compendium_protein(
        self, protein_raw: Optional[str]
    ) -> Optional[str]:
        if not protein_raw:
            return None
        compact = re.sub(r"\s+", " ", protein_raw.strip())

        insertion_dup = re.match(
            r"^([ACDEFGHIKLMNPQRSTVWY]+)\s+(\d+)-(\d+)\s+dup$",
            compact,
            re.IGNORECASE,
        )
        if insertion_dup:
            compact = f"{insertion_dup.group(1)[0].upper()}{insertion_dup.group(2)}ins"

        range_deletion = re.match(
            r"^(\d+)-(\d+)\s+del\s+([ACDEFGHIKLMNPQRSTVWY]+)$",
            compact,
            re.IGNORECASE,
        )
        if range_deletion:
            compact = (
                f"{range_deletion.group(3)[-1].upper()}{range_deletion.group(2)}del"
            )

        numeric_deletion = re.match(
            r"^(\d+)\s+del\s+([ACDEFGHIKLMNPQRSTVWY])$",
            compact,
            re.IGNORECASE,
        )
        if numeric_deletion:
            compact = (
                f"{numeric_deletion.group(2).upper()}{numeric_deletion.group(1)}del"
            )

        numeric_insertion = re.match(
            r"^(\d+)\s+ins\s+([ACDEFGHIKLMNPQRSTVWY]+)$",
            compact,
            re.IGNORECASE,
        )
        if numeric_insertion:
            compact = (
                f"{numeric_insertion.group(2)[0].upper()}"
                f"{numeric_insertion.group(1)}ins"
            )

        frameshift = re.match(
            r"^([ACDEFGHIKLMNPQRSTVWY]\d+)fs/\d+$",
            compact,
            re.IGNORECASE,
        )
        if frameshift:
            compact = f"{frameshift.group(1).upper()}fsX"

        splice = re.match(
            r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)sp$",
            compact,
            re.IGNORECASE,
        )
        if splice:
            compact = f"{splice.group(1).upper()}{splice.group(2)}X"

        return self._valid_table_protein(compact)

    def _extract_lqts_compendium_protein(
        self, prefix: str, mutation_type: str
    ) -> Optional[str]:
        patterns = [
            r"[ACDEFGHIKLMNPQRSTVWY]+\s+\d+-\d+\s+dup",
            r"\d+-\d+\s+del\s+[ACDEFGHIKLMNPQRSTVWY]+",
            r"\d+\s+del\s+[ACDEFGHIKLMNPQRSTVWY]",
            r"\d+\s+ins\s+[ACDEFGHIKLMNPQRSTVWY]+",
            r"[ACDEFGHIKLMNPQRSTVWY]\d+fs/\d+",
            r"[ACDEFGHIKLMNPQRSTVWY]\d+sp",
            r"[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWYX]",
        ]
        candidates: list[str] = []
        for pattern in patterns:
            candidates.extend(
                match.group(0)
                for match in re.finditer(pattern, prefix, flags=re.IGNORECASE)
            )
        if not candidates:
            return None
        return max(candidates, key=len)

    def _parse_lqts_compendium_status(
        self, count_ethnicity: str, status: str
    ) -> Optional[tuple[int, int, int, bool]]:
        count_match = re.match(r"\s*(\d+)", count_ethnicity)
        if not count_match:
            return None
        count = int(count_match.group(1))
        control_like = bool(
            re.search(r"\b(rare control|polymorphism)\b", status, re.IGNORECASE)
        )
        if control_like:
            return count, 0, count, True
        return count, count, 0, False

    def _lqts_compendium_variant(
        self,
        *,
        gene_symbol: str,
        protein_raw: Optional[str],
        count_ethnicity: str,
        status: str,
        source_location: str,
    ) -> Optional[dict]:
        protein = self._normalize_lqts_compendium_protein(protein_raw)
        if not protein:
            return None
        parsed_status = self._parse_lqts_compendium_status(count_ethnicity, status)
        if not parsed_status:
            return None
        count, affected_count, unaffected_count, control_like = parsed_status

        return {
            "gene_symbol": gene_symbol,
            "cdna_notation": None,
            "protein_notation": protein,
            "clinical_significance": "benign" if control_like else "pathogenic",
            "patients": {
                "count": count,
                "phenotype": "unaffected control"
                if control_like
                else f"{gene_symbol}-associated disease",
            },
            "penetrance_data": {
                "total_carriers_observed": count,
                "affected_count": affected_count,
                "unaffected_count": unaffected_count,
            },
            "individual_records": [],
            "functional_data": {"summary": "", "assays": []},
            "segregation_data": None,
            "population_frequency": None,
            "evidence_level": "medium",
            "source_location": source_location
            or "LQTS compendium summary of all mutations",
            "additional_notes": "Parsed via deterministic LQTS compendium table parser",
            "key_quotes": [],
        }

    def _parse_lqts_compendium_variant_row(
        self,
        stripped: str,
        gene_symbol: str,
        current_table_label: str,
    ) -> Optional[dict]:
        """Parse PMID 19841300-style mixed-gene LQTS compendium rows."""
        if stripped.startswith("|"):
            cells = [cell.strip() for cell in stripped.strip("|").split("|")]
            cells = [cell for cell in cells if cell]
            if (
                len(cells) < 7
                or set("".join(cells)) <= {"-", ":", " "}
                or not cells[0].upper().startswith(gene_symbol.upper())
            ):
                return None
            mutation_type_idx = next(
                (
                    idx
                    for idx, cell in enumerate(cells)
                    if cell.lower()
                    in {kind.lower() for kind in self.LQTS_COMPENDIUM_MUTATION_TYPES}
                ),
                None,
            )
            if mutation_type_idx is None or mutation_type_idx < 2:
                return None
            return self._lqts_compendium_variant(
                gene_symbol=gene_symbol,
                protein_raw=cells[mutation_type_idx - 1],
                count_ethnicity=cells[-2],
                status=cells[-1],
                source_location=current_table_label,
            )

        if not stripped.upper().startswith(gene_symbol.upper()):
            return None
        status_match = re.search(
            r"\s(?P<count_ethnicity>\d+\s*[A-Za-z]+(?:[>=][A-Za-z]+)*)\s+"
            r"(?P<status>Rare control|Polymorphism|Case)\s*$",
            stripped,
            re.IGNORECASE,
        )
        if not status_match:
            return None
        left = stripped[: status_match.start()].rstrip()
        type_match = re.search(
            r"\s(?P<mutation_type>"
            + "|".join(re.escape(kind) for kind in self.LQTS_COMPENDIUM_MUTATION_TYPES)
            + r")\s+",
            left,
            re.IGNORECASE,
        )
        if not type_match:
            return None
        protein_raw = self._extract_lqts_compendium_protein(
            left[: type_match.start()], type_match.group("mutation_type")
        )
        return self._lqts_compendium_variant(
            gene_symbol=gene_symbol,
            protein_raw=protein_raw,
            count_ethnicity=status_match.group("count_ethnicity"),
            status=status_match.group("status"),
            source_location=current_table_label,
        )

    def _parse_lqts_compendium_vertical_variants(
        self, full_text: str, gene_symbol: str
    ) -> list[dict]:
        """Parse PDF-converted compendium rows emitted one cell per line."""
        lines = [line.strip() for line in full_text.splitlines()]
        variants: list[dict] = []
        current_table_label = ""
        in_compendium = False
        gene_line_re = re.compile(r"^(KCNQ1|KCNH2|SCN5A)(?:\s+intron\s+\d+)?$")
        count_re = re.compile(r"^\d+\s*[A-Za-z]+(?:[>=][A-Za-z]+)*$")
        status_re = re.compile(r"^(Rare control|Polymorphism|Case)$", re.IGNORECASE)
        mutation_type_re = re.compile(
            r"^("
            + "|".join(re.escape(kind) for kind in self.LQTS_COMPENDIUM_MUTATION_TYPES)
            + r")\b",
            re.IGNORECASE,
        )

        i = 0
        while i < len(lines):
            stripped = lines[i]
            if re.match(
                r"^(?:eTable|(?:Supplementary|Supplemental)?\s*Table)\s+\w+",
                stripped,
                re.IGNORECASE,
            ):
                current_table_label = stripped
                in_compendium = (
                    "compendium summary of all mutations" in stripped.lower()
                )
                i += 1
                continue
            if not in_compendium or not gene_line_re.match(stripped):
                i += 1
                continue

            row_gene = stripped.split()[0].upper()
            j = i + 1
            while j < len(lines) and not gene_line_re.match(lines[j]):
                if re.match(r"^#{1,6}\s+", lines[j]):
                    break
                j += 1

            if row_gene != gene_symbol.upper():
                i = j
                continue

            cells = [cell for cell in lines[i + 1 : j] if cell]
            status_idx = next(
                (idx for idx, cell in enumerate(cells) if status_re.match(cell)),
                None,
            )
            if status_idx is None:
                i = j
                continue
            count_idx = next(
                (
                    idx
                    for idx in range(status_idx - 1, -1, -1)
                    if count_re.match(cells[idx])
                ),
                None,
            )
            mutation_type_idx = next(
                (
                    idx
                    for idx, cell in enumerate(cells[:status_idx])
                    if mutation_type_re.match(cell)
                ),
                None,
            )
            if count_idx is None or mutation_type_idx is None or mutation_type_idx == 0:
                i = j
                continue

            variant = self._lqts_compendium_variant(
                gene_symbol=gene_symbol,
                protein_raw=cells[mutation_type_idx - 1],
                count_ethnicity=cells[count_idx],
                status=cells[status_idx],
                source_location=current_table_label,
            )
            if variant:
                variants.append(variant)
            i = j

        return variants

    LINEARIZED_TABLE_CAPTION_RE = re.compile(
        r"^(?:eTable|(?:Supplementary|Supplemental)?\s*Table)\s+\w+.*"
        r"(?:mutation|variant|carrier|patient|proband)",
        re.IGNORECASE,
    )
    LINEARIZED_TABLE_STOP_RE = re.compile(
        r"^(?:eTable|eFigure|(?:Supplementary|Supplemental)?\s*"
        r"(?:Table|Figure))\s+\w+",
        re.IGNORECASE,
    )

    def _augment_pdf_linearized_tables(self, full_text: str) -> str:
        """Append markdown reconstructions of one-cell-per-line PDF tables."""
        # Reset the per-paper cDNA<->protein pairing map. The reconstructed
        # tables carry an explicit cDNA+protein pairing per row; downstream
        # parsers/scanners sometimes emit a row with only one notation, so we
        # stash the pairing here and backfill the missing side after extraction
        # (see _backfill_variant_notation_pairs).
        self._linearized_variant_pairs = {}
        blocks = self._reconstruct_pdf_linearized_tables(full_text)
        if not blocks:
            return full_text
        self._linearized_variant_pairs = self._pairs_from_reconstructed_blocks(blocks)
        logger.info("Reconstructed %d PDF-linearized table(s)", len(blocks))
        appendix = "\n\n".join(blocks)
        return (
            f"{full_text.rstrip()}\n\n"
            "## Reconstructed PDF-linearized tables\n\n"
            f"{appendix}"
        )

    @staticmethod
    def _notation_pair_key(notation: str) -> str:
        """Whitespace-insensitive, case-insensitive key for pairing notations."""
        return re.sub(r"\s+", "", (notation or "")).lower()

    def _pairs_from_reconstructed_blocks(self, blocks: list[str]) -> dict[str, str]:
        """Map each reconstructed cDNA<->protein cell to its row partner.

        Keyed by ``_notation_pair_key`` so a cDNA-only extracted row can recover
        its protein (and vice versa). cDNA and protein notations never collide
        (``c.`` vs ``p.``), so a single dict is unambiguous. First write wins.
        """
        pairs: dict[str, str] = {}
        for block in blocks:
            rows = [line for line in block.splitlines() if line.strip().startswith("|")]
            if len(rows) < 3:  # header, separator, >=1 data row
                continue
            header = [c.strip().lower() for c in rows[0].strip().strip("|").split("|")]
            if "cdna" not in header or "protein" not in header:
                continue
            cdna_idx = header.index("cdna")
            protein_idx = header.index("protein")
            for line in rows[2:]:  # skip header + separator row
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if max(cdna_idx, protein_idx) >= len(cells):
                    continue
                cdna, protein = cells[cdna_idx], cells[protein_idx]
                if not cdna or not protein:
                    continue
                pairs.setdefault(self._notation_pair_key(cdna), protein)
                pairs.setdefault(self._notation_pair_key(protein), cdna)
        return pairs

    def _backfill_variant_notation_pairs(
        self, extracted_data: Optional[dict]
    ) -> Optional[dict]:
        """Fill a missing cDNA or protein notation from the reconstructed pairing.

        Only fills the empty side; never overwrites a notation the extractor
        already produced. No-op when no PDF-linearized table was reconstructed.
        """
        pairs = getattr(self, "_linearized_variant_pairs", None)
        if not pairs or not extracted_data:
            return extracted_data
        for variant in extracted_data.get("variants", []) or []:
            if not isinstance(variant, dict):
                continue
            cdna = (variant.get("cdna_notation") or "").strip()
            protein = (variant.get("protein_notation") or "").strip()
            if cdna and not protein:
                partner = pairs.get(self._notation_pair_key(cdna))
                if partner and partner.lower().startswith("p."):
                    variant["protein_notation"] = partner
            elif protein and not cdna:
                partner = pairs.get(self._notation_pair_key(protein))
                if partner and partner.lower().startswith("c"):
                    variant["cdna_notation"] = partner
        return extracted_data

    def _reconstruct_pdf_linearized_tables(self, full_text: str) -> list[str]:
        """Return pipe-table blocks reconstructed from PDF-linearized cell runs.

        These supplements often preserve every table cell on its own line. We
        detect a variant/count header, then rebuild rows by variant-line
        boundaries and header-derived column positions so the existing
        deterministic markdown parser can consume the table.
        """
        lines = full_text.splitlines()
        blocks: list[str] = []
        i = 0
        while i < len(lines):
            header = self._find_linearized_table_header(lines, i)
            if header is None:
                i += 1
                continue

            caption, headers, data_start = header
            rows, end_idx = self._collect_linearized_table_rows(
                lines, data_start, headers
            )
            if rows:
                blocks.append(self._render_reconstructed_table(caption, headers, rows))
                i = max(end_idx, data_start + 1)
            else:
                i += 1
        return blocks

    def _find_linearized_table_header(
        self, lines: list[str], start_idx: int
    ) -> Optional[tuple[str, list[str], int]]:
        first = lines[start_idx].strip()
        if (
            not first
            or "|" in first
            or not self._looks_like_linearized_header_cell(first)
        ):
            return None

        raw_cells: list[str] = []
        j = start_idx
        max_header_lines = min(len(lines), start_idx + 24)
        while j < max_header_lines:
            cell = lines[j].strip()
            if not cell:
                j += 1
                continue
            if "|" in cell:
                return None
            if self._line_has_variant_notation(cell):
                break
            if j != start_idx and self.LINEARIZED_TABLE_STOP_RE.match(cell):
                return None
            if not self._looks_like_linearized_header_cell(cell):
                if raw_cells and re.fullmatch(r"[A-Za-z/]{1,12}", cell):
                    raw_cells.append(cell)
                    j += 1
                    continue
                if len(raw_cells) < 3:
                    return None
                break
            raw_cells.extend(self._split_linearized_header_line(cell))
            j += 1

        if j >= len(lines) or not self._line_has_variant_notation(lines[j].strip()):
            return None

        header_cells = self._normalize_linearized_header_cells(raw_cells)
        if not self._linearized_header_qualifies(header_cells):
            return None

        output_headers = self._linearized_output_headers(header_cells)
        caption = self._nearest_linearized_table_caption(lines, start_idx)
        if not caption and len(output_headers) < 5:
            return None
        return caption, output_headers, j

    def _nearest_linearized_table_caption(
        self, lines: list[str], header_idx: int
    ) -> str:
        for idx in range(header_idx - 1, max(-1, header_idx - 8), -1):
            candidate = lines[idx].strip()
            if not candidate:
                continue
            if self.LINEARIZED_TABLE_CAPTION_RE.match(candidate):
                return candidate
            if candidate.startswith("|"):
                break
        return ""

    def _split_linearized_header_line(self, line: str) -> list[str]:
        return [
            part.strip()
            for part in re.split(r"\t+|\s{2,}", line.strip())
            if part.strip()
        ]

    def _normalize_linearized_header_cells(self, cells: list[str]) -> list[str]:
        normalized: list[str] = []
        join_previous = {
            "acid",
            "change",
            "effect",
            "count",
            "counts",
            "number",
            "patients",
            "carriers",
            "observed",
        }
        previous_prefix = {
            "amino",
            "amino acid",
            "coding",
            "dna",
            "cdna",
            "nucleotide",
            "protein",
            "mean",
            "total",
            "number of",
            "no.",
        }
        for raw in cells:
            cell = re.sub(r"\s+", " ", raw.strip().strip("*")).strip()
            if not cell:
                continue
            lower = cell.lower()
            if normalized and re.fullmatch(r"[A-Za-z]{1,4}\s+\([^)]{1,30}\)", cell):
                normalized[-1] = f"{normalized[-1]}{cell}"
                continue
            if normalized and re.fullmatch(r"\([^)]{1,30}\)", cell):
                normalized[-1] = f"{normalized[-1]} {cell}"
                continue
            if normalized and (
                lower in join_previous or normalized[-1].lower() in previous_prefix
            ):
                normalized[-1] = f"{normalized[-1]} {cell}"
                continue
            normalized.append(cell)
        return normalized

    def _linearized_header_qualifies(self, cells: list[str]) -> bool:
        if len(cells) < 4 or len(cells) > 18:
            return False
        normalized = [re.sub(r"[^a-z0-9./+-]+", " ", c.lower()).strip() for c in cells]
        has_variant = any(
            any(term in cell for term in ("mutation", "variant", "hgvs", "protein"))
            or cell in {"cdna", "c.", "nucleotide"}
            for cell in normalized
        )
        has_count = any(
            cell in {"n", "no", "number"}
            or "(n)" in cell
            or any(
                term in cell
                for term in (
                    "patient",
                    "carrier",
                    "proband",
                    "affected",
                    "unaffected",
                    "syncope",
                    "case",
                )
            )
            for cell in normalized
        )
        return has_variant and has_count

    def _linearized_output_headers(self, cells: list[str]) -> list[str]:
        headers: list[str] = []
        inserted_variant_columns = False
        for cell in cells:
            lower = re.sub(r"[^a-z0-9./+-]+", " ", cell.lower()).strip()
            if (
                not inserted_variant_columns
                and any(term in lower for term in ("mutation", "variant", "hgvs"))
                and "type" not in lower
            ):
                headers.extend(["cDNA", "Protein"])
                inserted_variant_columns = True
                continue
            if lower == "n":
                headers.append("No. of patients")
            elif "syncope" in lower:
                headers.append("affected")
            else:
                headers.append(cell)
        if not inserted_variant_columns:
            headers = ["cDNA", "Protein", *headers]
        return headers

    def _collect_linearized_table_rows(
        self, lines: list[str], data_start: int, headers: list[str]
    ) -> tuple[list[list[str]], int]:
        rows: list[list[str]] = []
        current: Optional[list[str]] = None
        i = data_start
        while i < len(lines):
            cell = lines[i].strip()
            if not cell or self._is_linearized_page_noise(cell):
                i += 1
                continue
            if rows and self.LINEARIZED_TABLE_STOP_RE.match(cell):
                break
            if rows and self._is_linearized_table_footnote_or_prose(cell):
                break
            if "|" in cell and rows:
                break

            if self._line_has_variant_notation(cell):
                if current:
                    aligned = self._align_linearized_row_cells(current, headers)
                    if aligned:
                        rows.append(aligned)
                current = self._split_linearized_variant_cell(cell)
            elif current and self._looks_like_linearized_data_cell(cell):
                current.append(cell)
            elif rows:
                break
            else:
                return [], i + 1
            i += 1

        if current:
            aligned = self._align_linearized_row_cells(current, headers)
            if aligned:
                rows.append(aligned)
        return rows, i

    def _line_has_variant_notation(self, line: str) -> bool:
        cleaned = line.strip()
        if not cleaned or len(cleaned) > 140:
            return False
        return bool(
            re.search(
                r"(?:^|[\s*])(?:[cｃ]\.?\s*\d|p\.\s*"
                r"(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY0-9]))",
                cleaned,
                re.IGNORECASE,
            )
        )

    def _split_linearized_variant_cell(self, cell: str) -> list[str]:
        cdna = ""
        protein = ""
        consumed_end = 0
        cdna_match = re.search(
            r"(?:[cｃ]\.?\s*\d+(?:[+-]\d+|[-_]\d+)?\s*"
            r"(?:[ACGTacgt]\s*>\s*[ACGTacgt]|del\s*[ACGTacgt]*|"
            r"dup\s*[ACGTacgt]*|ins\s*[ACGTacgt]*))",
            cell,
        )
        if cdna_match:
            cdna = (
                self._valid_table_cdna(
                    cdna_match.group(0).replace("ｃ", "c").replace("\u3000", "")
                )
                or ""
            )
            consumed_end = max(consumed_end, cdna_match.end())

        protein_match = re.search(r"\bp\.\s*", cell, re.IGNORECASE)
        if protein_match:
            protein_tokens: list[str] = []
            cursor = protein_match.end()
            for match in re.finditer(r"\S+", cell[protein_match.end() :]):
                token_start = protein_match.end() + match.start()
                token_end = protein_match.end() + match.end()
                token = match.group(0).strip(",;")
                if self._linearized_variant_token_stops_protein(token):
                    break
                protein_tokens.append(token)
                cursor = token_end
            if protein_tokens:
                protein = (
                    self._valid_table_protein("p." + "".join(protein_tokens)) or ""
                )
                consumed_end = max(consumed_end, cursor)

        cells = [cdna, protein]
        remainder = cell[consumed_end:].strip(" ,;")
        if remainder and self._looks_like_linearized_variant_remainder(remainder):
            cells.append(remainder)
        return cells

    def _linearized_variant_token_stops_protein(self, token: str) -> bool:
        clean = token.strip().strip(",;")
        if not clean or clean.startswith("+") or clean.isdigit():
            return True
        upper = clean.upper()
        lower = clean.lower()
        return (
            upper in {"N", "C", "MS", "TM", "DI", "DII", "DIII", "DIV"}
            or lower in {"c-loop", "s5", "s6", "diii-div"}
            or "term" in lower
            or "pore" in lower
            or lower.startswith("splicing")
            or lower == "error"
        )

    def _looks_like_linearized_variant_remainder(self, value: str) -> bool:
        lower = value.lower()
        if "splicing" in lower or "error" in lower:
            return False
        if len(value) > 35:
            return False
        return bool(
            re.fullmatch(r"[A-Za-z0-9/().+-]+(?:\s+[A-Za-z0-9/().+-]+)?", value)
        )

    def _align_linearized_row_cells(
        self, cells: list[str], headers: list[str]
    ) -> Optional[list[str]]:
        if not any(cells[:2]) or not any(self._is_integer_cell(c) for c in cells[2:]):
            return None

        count_idx = self._linearized_count_header_index(headers)
        if count_idx is None or count_idx < 2:
            row = cells[: len(headers)]
            return row + [""] * (len(headers) - len(row))

        row = cells[:2]
        rest = cells[2:]
        if (
            row[1]
            and re.search(r"fs(?:x|ter|\*)$", row[1], re.IGNORECASE)
            and rest
            and re.fullmatch(r"\d{2,4}", rest[0])
            and len(rest) > 1
            and not self._is_integer_cell(rest[1])
        ):
            row[1] = f"{row[1]}{rest.pop(0)}"
        descriptor_slots = max(0, count_idx - 2)
        descriptors: list[str] = []
        while rest and len(descriptors) < descriptor_slots:
            if self._is_integer_cell(rest[0]):
                break
            descriptors.append(rest.pop(0))
        descriptors.extend([""] * (descriptor_slots - len(descriptors)))

        numeric_slots = len(headers) - count_idx
        numeric = rest[:]
        mean_qtc_idx = next(
            (
                idx
                for idx, header in enumerate(headers)
                if "qtc" in header.lower() or "qt interval" in header.lower()
            ),
            None,
        )
        if len(numeric) == numeric_slots - 1 and mean_qtc_idx is not None:
            relative_qtc_idx = mean_qtc_idx - count_idx
            if 0 <= relative_qtc_idx <= len(numeric):
                qtc_value = (
                    int(numeric[relative_qtc_idx])
                    if relative_qtc_idx < len(numeric)
                    and numeric[relative_qtc_idx].isdigit()
                    else None
                )
                if qtc_value is None or qtc_value < 100:
                    numeric.insert(relative_qtc_idx, "")
        numeric = numeric[:numeric_slots]
        numeric.extend([""] * (numeric_slots - len(numeric)))

        row.extend(descriptors)
        row.extend(numeric)
        row = row[: len(headers)]
        row.extend([""] * (len(headers) - len(row)))
        return row

    def _linearized_count_header_index(self, headers: list[str]) -> Optional[int]:
        for idx, header in enumerate(headers):
            normalized = re.sub(r"[^a-z0-9]+", " ", header.lower()).strip()
            if normalized in {"n", "no of patients", "number of patients"}:
                return idx
            if "patient" in normalized or "carrier" in normalized:
                return idx
        return None

    def _is_integer_cell(self, value: str) -> bool:
        return bool(re.fullmatch(r"\d+", value.strip()))

    def _looks_like_linearized_header_cell(self, cell: str) -> bool:
        if len(cell) > 80 or cell.startswith("#"):
            return False
        lower = cell.lower().strip()
        if re.fullmatch(r"\([^)]{1,30}\)", cell) or re.fullmatch(
            r"[A-Za-z]{1,4}\s+\([^)]{1,30}\)", cell
        ):
            return True
        return bool(
            re.search(
                r"\b(mutation|variant|hgvs|protein|nucleotide|cdna|site|"
                r"female|male|proband|patient|carrier|mean|qtc|syncope|"
                r"ca/vf|affected|unaffected|number|count|n)\b",
                lower,
            )
        )

    def _looks_like_linearized_data_cell(self, cell: str) -> bool:
        if len(cell) > 90 or cell.startswith("#"):
            return False
        if self._is_integer_cell(cell):
            return True
        if re.fullmatch(r"[A-Za-z0-9/().+-]+(?:\s+[A-Za-z0-9/().+-]+){0,4}", cell):
            return True
        return False

    def _is_linearized_page_noise(self, cell: str) -> bool:
        return bool(
            re.match(r"^#{1,6}\s+Page\s+\d+", cell, re.IGNORECASE)
            or cell.startswith("©")
            or "all rights reserved" in cell.lower()
        )

    def _is_linearized_table_footnote_or_prose(self, cell: str) -> bool:
        if re.match(r"^[A-Z]:\s", cell):
            return True
        if len(cell) > 90 and not self._line_has_variant_notation(cell):
            return True
        if cell.endswith(".") and len(cell.split()) > 7:
            return True
        return False

    def _render_reconstructed_table(
        self, caption: str, headers: list[str], rows: list[list[str]]
    ) -> str:
        def clean(value: str) -> str:
            cleaned = re.sub(r"\s+", " ", (value or "").replace("|", "/")).strip()
            return cleaned or "-"

        block: list[str] = []
        if caption:
            block.append(clean(caption))
        else:
            block.append("Reconstructed PDF-linearized variant table")
        block.append("| " + " | ".join(clean(header) for header in headers) + " |")
        block.append("| " + " | ".join("---" for _ in headers) + " |")
        for row in rows:
            block.append("| " + " | ".join(clean(cell) for cell in row) + " |")
        return "\n".join(block)

    def _parse_markdown_table_variants(
        self, full_text: str, gene_symbol: Optional[str]
    ) -> List[dict]:
        """
        Best-effort parser for simple markdown tables (fast path for very large tables).

        Returns a minimal variant list without calling the LLM. Intended for papers
        like PMID 19716085 where a single giant table lists hundreds of variants.

        BROADENED (2026-02-10): Now recognizes many more header patterns:
        - cDNA: nucleotide, cDNA, c., HGVS cdna, etc.
        - Protein: variant, amino acid, protein, p., AAChange, mutation, etc.
        - Counts: patient(s), N, affected, carriers, probands, families, subjects, etc.
        """
        if not gene_symbol:
            return []

        lines = full_text.splitlines()
        table_started = False
        variants = []
        header_idx = {}
        header_mapping = {}  # Maps our standard keys to actual column indices
        header_multi = {"count": [], "affected": [], "unaffected": []}
        active_headers = []
        table_row_ordinal = 0
        row_level_clinical_table = False
        current_table_label = ""
        active_table_label = ""
        normalized_active_headers = []

        from pipeline.table_router import (
            _is_non_variant_count_header,
            _normalize_header,
        )

        def count_type_for_header(label: Optional[str]) -> Optional[str]:
            if not label:
                return None
            normalized = _normalize_header(label)
            if not normalized:
                return None
            if "family" in normalized or "kindred" in normalized:
                return "family_count"
            if "proband" in normalized:
                return "proband_count"
            if "control" in normalized:
                return "unaffected_control"
            if normalized in {"case", "cases"}:
                return "case"
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
            if any(
                token in normalized for token in ("screened", "ntested", "numbertested")
            ):
                return "screened_N"
            return "per_variant_carrier"

        def update_table_label(text: str) -> None:
            nonlocal current_table_label
            candidate = text
            # Captions are frequently emitted as a single-cell table row, e.g.
            # `| **Supplemental Table 1. KCNQ1 variants ...** | | | |`. Pull the
            # first non-empty cell so a gene named only in the caption (and not
            # in a per-row Gene column) still scopes the table.
            if candidate.lstrip().startswith("|"):
                row_cells = [c.strip() for c in candidate.split("|") if c.strip()]
                candidate = row_cells[0] if row_cells else ""
            heading = candidate.lstrip("#").strip().replace("*", "").strip()
            if re.match(
                r"^(?:etable|(?:supplement\w*\s+)?table)\s+\w+(?:[\.:]|\s|$)",
                heading,
                re.IGNORECASE,
            ):
                current_table_label = heading

        def table_label_matches_target(label: str) -> bool:
            if not label:
                return True
            mentioned = set(self._gene_scope_from_table_label(label))
            mentioned.update(self._contextual_gene_symbols_in_table_label(label))
            if not mentioned:
                mentioned = self._gene_symbols_in_table_label(label)
            return not mentioned or gene_symbol.upper() in mentioned

        for line_number, line in enumerate(lines, start=1):
            stripped_line = line.strip()
            if not table_started:
                update_table_label(stripped_line)

                # Detect header row using broadened patterns
                if self._is_variant_table_header(line):
                    parts = [p.strip() for p in line.split("|") if p.strip()]
                    if self._looks_like_gwas_header(parts):
                        table_started = False
                        header_mapping = {}
                        header_idx = {}
                        header_multi = {"count": [], "affected": [], "unaffected": []}
                        row_level_clinical_table = False
                        continue
                    row_level_clinical_table = (
                        self._looks_like_row_level_clinical_header(parts)
                    )
                    active_table_label = current_table_label
                    active_headers = parts
                    normalized_active_headers = [
                        _normalize_header(header) for header in active_headers
                    ]
                    table_row_ordinal = 0
                    for idx, name in enumerate(parts):
                        name_lower = name.lower().strip()
                        header_idx[name_lower] = idx

                        # Map to standard keys for retrieval
                        if self._header_matches(name, self.VARIANT_CDNA_HEADERS):
                            header_mapping["cdna"] = idx
                        if (
                            self._header_matches(name, self.VARIANT_PROTEIN_HEADERS)
                            and "type" not in name_lower
                        ):
                            header_mapping["protein"] = idx
                        if self._header_matches(name, self.VARIANT_GENE_HEADERS):
                            header_mapping["gene"] = idx
                        if self._header_matches(name, self.VARIANT_COUNT_HEADERS):
                            header_mapping["count"] = idx
                            header_multi["count"].append(idx)
                        if self._header_matches(name, self.VARIANT_ROW_LEVEL_HEADERS):
                            header_mapping["row_subject"] = idx
                        is_unaffected_header = self._header_matches(
                            name, self.VARIANT_UNAFFECTED_HEADERS
                        )
                        if is_unaffected_header:
                            header_mapping["unaffected"] = idx
                            header_multi["unaffected"].append(idx)
                        elif self._header_matches(name, self.VARIANT_AFFECTED_HEADERS):
                            header_mapping["affected"] = idx
                            header_multi["affected"].append(idx)

                    table_started = True
                continue

            # Stop when table ends
            if not stripped_line.startswith("|"):
                table_started = False
                header_mapping = {}
                header_idx = {}
                header_multi = {"count": [], "affected": [], "unaffected": []}
                active_headers = []
                normalized_active_headers = []
                table_row_ordinal = 0
                row_level_clinical_table = False
                active_table_label = ""
                continue

            # Skip separator rows
            if set(line.strip()) <= {"|", "-", " ", ":"}:
                continue

            cells = [c.strip() for c in line.split("|") if c.strip()]
            if len(cells) < 2:  # Relaxed from 3 to 2 for simpler tables
                continue
            table_row_ordinal += 1

            def usable_indices(key: str) -> List[int]:
                raw_indices = list(header_multi.get(key, []))
                mapped_idx = header_mapping.get(key)
                if mapped_idx is not None and mapped_idx not in raw_indices:
                    raw_indices.append(mapped_idx)
                if not raw_indices:
                    direct_idx = header_idx.get(key)
                    if direct_idx is not None:
                        raw_indices.append(direct_idx)

                indices: List[int] = []
                for idx in raw_indices:
                    if idx is None or idx >= len(cells):
                        continue
                    if (
                        key in {"count", "affected", "unaffected"}
                        and idx < len(normalized_active_headers)
                        and _is_non_variant_count_header(
                            normalized_active_headers[idx], normalized_active_headers
                        )
                    ):
                        continue
                    indices.append(idx)
                return indices

            def get_col(key: str) -> Optional[str]:
                indices = usable_indices(key)
                if not indices:
                    return None
                idx = indices[0]
                return cells[idx] or None

            def get_cols(key: str) -> List[str]:
                values = []
                for idx in usable_indices(key):
                    if idx < len(cells) and cells[idx]:
                        values.append(cells[idx])
                return values

            def header_labels(key: str) -> List[str]:
                labels = []
                for idx in usable_indices(key):
                    if idx < len(active_headers):
                        label = active_headers[idx].strip()
                        if label:
                            labels.append(label)
                return labels

            def header_label(key: str) -> Optional[str]:
                labels = header_labels(key)
                return "; ".join(labels) if labels else None

            cdna = get_col("cdna")
            protein = get_col("protein")
            row_gene = get_col("gene")
            patient_count_raw = get_col("count")
            affected_raw = get_col("affected")
            unaffected_raw = get_col("unaffected")
            row_subject_raw = get_col("row_subject")

            if row_gene and row_gene.strip().upper() != gene_symbol.upper():
                continue

            row_text = " ".join(cells)
            if not table_label_matches_target(active_table_label):
                continue
            row_control_like = bool(
                re.search(
                    r"\b(rare control|common polymorphism|polymorphism|control variant)\b",
                    row_text,
                    re.IGNORECASE,
                )
            )

            # Validate actual cell values. Header cues alone are not enough:
            # GWAS tables use AA=alternate allele and n=participants, which
            # previously produced fake variants like A/G/T/C with huge counts.
            cdna = self._valid_table_cdna(cdna)
            protein = self._valid_table_protein(protein)

            if not cdna and not protein:
                continue

            def parse_count(value: Optional[str]) -> Optional[int]:
                if not value:
                    return None
                cleaned = value.strip().replace(",", "")
                if cleaned.isdigit():
                    return int(cleaned)
                return None

            def parse_count_sum(values: List[str]) -> Optional[int]:
                parsed = [parse_count(value) for value in values]
                parsed = [value for value in parsed if value is not None]
                if not parsed:
                    return None
                return sum(parsed)

            patient_count = parse_count_sum(get_cols("count")) or parse_count(
                patient_count_raw
            )
            affected_count = parse_count_sum(get_cols("affected")) or parse_count(
                affected_raw
            )
            unaffected_count = parse_count_sum(get_cols("unaffected")) or parse_count(
                unaffected_raw
            )

            # Fallback: use last cell if it looks numeric
            if (
                patient_count is None
                and affected_count is None
                and unaffected_count is None
                and patient_count_raw is None
                and row_subject_raw is None
                and cells
            ):
                tail = cells[-1]
                if tail.isdigit():
                    patient_count = int(tail)

            if (
                patient_count is None
                and affected_count is None
                and unaffected_count is None
                and row_level_clinical_table
            ):
                patient_count = 1
                row_text = " ".join(cells)
                if re.search(
                    r"\b(unaffected|asymptomatic|control|healthy|normal|no symptoms?)\b",
                    row_text,
                    re.IGNORECASE,
                ):
                    affected_count = 0
                    unaffected_count = 1
                else:
                    affected_count = 1
                    unaffected_count = 0

            if patient_count is None and (
                affected_count is not None or unaffected_count is not None
            ):
                patient_count = (affected_count or 0) + (unaffected_count or 0)

            if row_control_like and patient_count is not None:
                affected_count = 0
                unaffected_count = patient_count
            elif patient_count is not None:
                if affected_count is None and unaffected_count is None:
                    affected_count = patient_count
                    unaffected_count = 0
                elif affected_count is None and unaffected_count is not None:
                    remainder = patient_count - unaffected_count
                    affected_count = remainder if remainder >= 0 else None
                elif unaffected_count is None and affected_count is not None:
                    remainder = patient_count - affected_count
                    unaffected_count = remainder if remainder >= 0 else None

            source_ref = active_table_label or "Markdown table"
            observation_provenance = {
                "source_container": "supplement"
                if re.search(
                    r"\bsupp(?:lement(?:ary|al)?)?\b|etable|e-table",
                    source_ref,
                    re.IGNORECASE,
                )
                else "main",
                "source_kind": "table",
                "source_ref": source_ref,
                "row_ordinal": table_row_ordinal,
                "column_ref": header_label("count"),
                "locator_extra": {
                    "parser": "markdown_table",
                    "line_number": line_number,
                },
            }
            carrier_label = header_label("count")
            affected_label = header_label("affected") or carrier_label
            unaffected_label = header_label("unaffected")
            variant = {
                "gene_symbol": gene_symbol,
                "cdna_notation": f"c.{cdna}"
                if cdna and not cdna.startswith("c.")
                else cdna,
                "protein_notation": protein,
                "clinical_significance": "benign" if row_control_like else "pathogenic",
                "patients": {
                    "count": patient_count,
                    "phenotype": "unaffected control"
                    if row_control_like
                    else f"{gene_symbol}-associated disease",
                    **observation_provenance,
                },
                "penetrance_data": {
                    "total_carriers_observed": patient_count,
                    "affected_count": affected_count,
                    "unaffected_count": unaffected_count,
                },
                "individual_records": [],
                "functional_data": {"summary": "", "assays": []},
                "segregation_data": None,
                "population_frequency": None,
                "evidence_level": "medium",
                "source_location": active_table_label or "Markdown table",
                **observation_provenance,
                "additional_notes": "Parsed via deterministic table parser",
                "key_quotes": [],
                "count_provenance": {
                    "carriers_column_label": carrier_label,
                    "carriers_count_type": count_type_for_header(carrier_label),
                    "affected_column_label": affected_label,
                    "affected_count_type": count_type_for_header(affected_label),
                    "unaffected_column_label": unaffected_label,
                    "unaffected_count_type": count_type_for_header(unaffected_label),
                },
            }
            variants.append(variant)

        if variants:
            logger.info(
                f"Parsed {len(variants)} variants via deterministic markdown table parser"
            )
        return variants

    def _parse_vertical_gene_table_variants(
        self, full_text: str, gene_symbol: Optional[str]
    ) -> List[dict]:
        """
        Parse stacked PDF/HTML supplement tables where each cell is on its own line.

        Example shape after text extraction:

            Gene
            Nucleotide
            Amino acids
            ...
            RYR2
            1258c>t
            R420W
            NT

        This is common in publisher/PDF supplements that are not valid markdown
        tables or fixed-width rows. We only emit rows where the target gene is an
        isolated cell and a valid cDNA or protein token appears immediately after.
        """
        if not gene_symbol:
            return []

        variants: List[dict] = []
        by_key: dict[tuple[str, str], dict] = {}
        current_table_label = ""
        target_gene = gene_symbol.strip().upper()
        lines = full_text.splitlines()

        def add_variant(cdna: Optional[str], protein: Optional[str]) -> None:
            key = ((cdna or "").lower(), (protein or "").lower())
            if not key[0] and not key[1]:
                return
            if key in by_key:
                existing = by_key[key]
                patients = existing.setdefault("patients", {})
                pdata = existing.setdefault("penetrance_data", {})
                for obj, field in (
                    (patients, "count"),
                    (pdata, "total_carriers_observed"),
                    (pdata, "affected_count"),
                ):
                    value = obj.get(field)
                    obj[field] = 1 if value is None else value + 1
                if pdata.get("unaffected_count") is None:
                    pdata["unaffected_count"] = 0
                return

            variant = {
                "gene_symbol": gene_symbol,
                "cdna_notation": cdna,
                "protein_notation": protein,
                "clinical_significance": "pathogenic",
                "patients": {
                    "count": 1,
                    "phenotype": f"{gene_symbol}-associated disease",
                },
                "penetrance_data": {
                    "total_carriers_observed": 1,
                    "affected_count": 1,
                    "unaffected_count": 0,
                },
                "individual_records": [],
                "functional_data": {"summary": "", "assays": []},
                "segregation_data": None,
                "population_frequency": None,
                "evidence_level": "medium",
                "source_location": current_table_label or "Vertical gene table",
                "additional_notes": "Parsed via deterministic vertical gene-table parser",
                "key_quotes": [],
            }
            by_key[key] = variant
            variants.append(variant)

        for idx, line in enumerate(lines):
            stripped = line.strip()
            if re.match(
                r"^(?:Supplementary|Supplemental)?\s*Table\s+\w+",
                stripped,
                re.IGNORECASE,
            ):
                current_table_label = stripped[:240]

            if stripped.upper() != target_gene:
                continue

            window: list[str] = []
            for raw in lines[idx + 1 : idx + 9]:
                cell = raw.strip()
                if not cell:
                    continue
                window.append(cell)
                if len(window) >= 6:
                    break

            cdna = None
            protein = None
            for value in window:
                # Cells sometimes contain parallel nucleotide changes separated
                # by commas; trying the whole cell and then pieces keeps the
                # common simple cases while preserving at least the protein
                # notation for complex multi-base substitutions.
                pieces = [value] + [
                    part for part in re.split(r"\s*[,;]\s*", value) if part
                ]
                for piece in pieces:
                    if cdna is None:
                        cdna = self._valid_table_cdna(piece)
                    if protein is None:
                        protein = self._valid_table_protein(piece)
                if cdna and protein:
                    break

            if cdna or protein:
                add_variant(cdna, protein)

        if variants:
            logger.info(
                f"Parsed {len(variants)} variants via deterministic vertical gene-table parser"
            )
        return variants

    def _parse_fixed_width_table_variants(
        self, full_text: str, gene_symbol: Optional[str]
    ) -> List[dict]:
        """
        Parse pdftotext -layout style variant tables.

        Some recovered publisher/manuscript PDFs preserve tables as fixed-width
        rows rather than markdown pipes, e.g.:

            c5350G>A  28  pGlu1784Lys  Gain of function  69

        This keeps large supplement tables on the deterministic path instead of
        asking the LLM to reconstruct hundreds of row-aligned variants.
        """
        if not gene_symbol:
            return []

        variants: List[dict] = []
        seen = set()
        current_table_label = ""
        protein_count_table = False
        nssnv_case_control_table = False
        functional_status_table = False
        clinical_mutation_table = False
        target_gene_lower = gene_symbol.lower()
        target_gene_upper = gene_symbol.upper()
        current_table_gene_scope: set[str] = set()
        effect_phrases = [
            ("Gain", "and", "loss"),
            ("Loss", "of", "function"),
            ("Gain", "of", "function"),
        ]

        for line in full_text.splitlines():
            stripped = line.strip()
            if re.match(
                r"^(?:eTable|(?:Supplementary|Supplemental)?\s*Table)\s+\w+",
                stripped,
                re.IGNORECASE,
            ):
                current_table_label = stripped
                label_lower = current_table_label.lower()
                if "continued" not in label_lower:
                    current_table_gene_scope = (
                        self._fixed_width_gene_scope_from_table_label(
                            current_table_label
                        )
                    )
                    protein_count_table = (
                        "list of mutations by coding effect" in label_lower
                        and "frequency" in label_lower
                    )
                    label_mentions_target = (
                        target_gene_lower in label_lower
                        or target_gene_upper in current_table_gene_scope
                    )
                    nssnv_case_control_table = (
                        label_mentions_target
                        and "nssnv" in label_lower
                        and ("properties" in label_lower or "variant" in label_lower)
                    )
                    functional_status_table = (
                        label_mentions_target
                        and "functionally characterized" in label_lower
                    )
                    clinical_mutation_table = (
                        label_mentions_target
                        and "mutation" in label_lower
                        and "variant" in label_lower
                    )

            lower = stripped.lower()
            table_has_off_target_gene_scope = (
                bool(current_table_gene_scope)
                and target_gene_upper not in current_table_gene_scope
            )
            if protein_count_table and lower.startswith("total"):
                protein_count_table = False
                continue
            if (
                current_table_label
                and "list of mutations by coding effect" in current_table_label.lower()
                and "coding effect" in lower
                and "count" in lower
            ):
                protein_count_table = True
                continue
            if "nucleotide" in lower and "coding effect" in lower:
                clinical_mutation_table = True
                continue
            if (
                "status" in lower
                and "variant" in lower
                and "control" in lower
                and ("case" in lower or "brs" in lower or "lqt" in lower)
            ):
                nssnv_case_control_table = True
                continue
            if table_has_off_target_gene_scope:
                continue

            tokens = line.strip().split()
            lqt_mutation_table = (
                target_gene_upper in current_table_gene_scope
                and "mutations or rare variants" in current_table_label.lower()
            )
            if lqt_mutation_table:
                if stripped.startswith("|"):
                    continue
                variant = self._parse_lqt_fixed_width_variant_row(
                    stripped,
                    gene_symbol,
                    current_table_label,
                )
                if variant:
                    key = (
                        variant.get("cdna_notation") or "",
                        variant.get("protein_notation") or "",
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    variants.append(variant)
                    continue

            lqts_compendium_table = (
                "compendium summary of all mutations" in current_table_label.lower()
            )
            if lqts_compendium_table:
                variant = self._parse_lqts_compendium_variant_row(
                    stripped,
                    gene_symbol,
                    current_table_label,
                )
                if variant:
                    key = (
                        variant.get("cdna_notation") or "",
                        variant.get("protein_notation") or "",
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    variants.append(variant)
                    continue

            if protein_count_table and len(tokens) >= 3 and tokens[-1].isdigit():
                protein = self._valid_table_protein(tokens[0])
                if not protein:
                    continue
                count = int(tokens[-1])
                key = ("", protein)
                if key in seen:
                    continue
                seen.add(key)
                variants.append(
                    {
                        "gene_symbol": gene_symbol,
                        "cdna_notation": None,
                        "protein_notation": protein,
                        "clinical_significance": "pathogenic",
                        "patients": {
                            "count": count,
                            "phenotype": f"{gene_symbol}-associated disease",
                        },
                        "penetrance_data": {
                            "total_carriers_observed": count,
                            "affected_count": count,
                            "unaffected_count": 0,
                        },
                        "individual_records": [],
                        "functional_data": {"summary": "", "assays": []},
                        "segregation_data": None,
                        "population_frequency": None,
                        "evidence_level": "medium",
                        "source_location": current_table_label
                        or "Supplemental coding-effect count table",
                        "additional_notes": "Parsed via deterministic fixed-width coding-effect count table parser",
                        "key_quotes": [],
                    }
                )
                continue

            if clinical_mutation_table:
                match = re.match(
                    r"^(?P<cdna>(?:\d+(?:_\d+)?\s*(?:[ACGT]>[ACGT]|del[ACGT]+|dup[ACGT]+|ins[ACGT]+)|IVS\d+[+-]\d+\s*(?:[ACGT]>[ACGT]|del[ACGT]+)))"
                    r"\s+(?P<protein>\S+)?"
                    r"(?:\s+\((?P<count>\d+)\))?\b",
                    stripped,
                    re.IGNORECASE,
                )
                if match:
                    cdna_raw = re.sub(r"\s+", "", match.group("cdna"))
                    protein_raw = match.group("protein")
                    count_match = None
                    if protein_raw:
                        count_match = re.match(
                            r"\s+\((\d+)\)", stripped[match.end("protein") :]
                        )
                    count = int(
                        count_match.group(1)
                        if count_match
                        else match.group("count") or "1"
                    )
                    cdna = (
                        f"c.{cdna_raw}"
                        if not cdna_raw.upper().startswith("IVS")
                        else cdna_raw
                    )
                    if not cdna.upper().startswith("IVS"):
                        cdna = self._valid_table_cdna(cdna)
                    protein = self._valid_table_protein(protein_raw)
                    if not cdna and not protein:
                        continue
                    key = (cdna or "", protein or "")
                    if key in seen:
                        continue
                    seen.add(key)
                    variants.append(
                        {
                            "gene_symbol": gene_symbol,
                            "cdna_notation": cdna,
                            "protein_notation": protein,
                            "clinical_significance": "pathogenic",
                            "patients": {
                                "count": count,
                                "phenotype": f"{gene_symbol}-associated disease",
                            },
                            "penetrance_data": {
                                "total_carriers_observed": count,
                                "affected_count": count,
                                "unaffected_count": 0,
                            },
                            "individual_records": [],
                            "functional_data": {"summary": "", "assays": []},
                            "segregation_data": None,
                            "population_frequency": None,
                            "evidence_level": "medium",
                            "source_location": current_table_label
                            or "Fixed-width clinical mutation table",
                            "additional_notes": "Parsed via deterministic fixed-width clinical mutation table parser",
                            "key_quotes": [],
                        }
                    )
                    continue

            if nssnv_case_control_table:
                match = re.match(
                    r"^(case\s+nsSNV|control\s+nsSNV|Polymorphism)\s+"
                    r"(\d+)\s+([ACGT]>[ACGT])\s+(\S+)\s+\S+\s+"
                    r"(\d+)\s+(\d+)\s+(\d+)\b",
                    stripped,
                    re.IGNORECASE,
                )
                if match:
                    status, pos, change, protein_raw, brs, lqt, control = match.groups()
                    protein = self._valid_table_protein(protein_raw)
                    if not protein:
                        continue
                    cdna = self._valid_table_cdna(f"{pos}{change}")
                    affected_count = int(brs) + int(lqt)
                    unaffected_count = int(control)
                    patient_count = affected_count + unaffected_count
                    key = (cdna or "", protein)
                    if key in seen:
                        continue
                    seen.add(key)
                    variants.append(
                        {
                            "gene_symbol": gene_symbol,
                            "cdna_notation": cdna,
                            "protein_notation": protein,
                            "clinical_significance": "pathogenic"
                            if affected_count
                            else "benign",
                            "patients": {
                                "count": patient_count,
                                "phenotype": f"{gene_symbol} nsSNV carrier",
                            },
                            "penetrance_data": {
                                "total_carriers_observed": patient_count,
                                "affected_count": affected_count,
                                "unaffected_count": unaffected_count,
                            },
                            "individual_records": [],
                            "functional_data": {"summary": "", "assays": []},
                            "segregation_data": None,
                            "population_frequency": None,
                            "evidence_level": "medium",
                            "source_location": current_table_label
                            or f"Supplemental {gene_symbol} nsSNV table",
                            "additional_notes": f"Parsed via deterministic fixed-width nsSNV table parser; status={status}",
                            "key_quotes": [],
                        }
                    )
                    continue

            if functional_status_table:
                match = re.match(
                    r"^((?:p\.)?[A-Z][a-z]{0,2}\d+[A-Z][a-z]{0,2})\s+"
                    r"\d+(?:\.\d+)?%\s+"
                    r"(Wildtype|Abnormal|Conflicting(?:\s+Data)?)\b",
                    stripped,
                    re.IGNORECASE,
                )
                if match:
                    protein_raw, ep_status = match.groups()
                    protein = self._valid_table_protein(protein_raw)
                    if not protein:
                        continue
                    if any(existing_protein == protein for _, existing_protein in seen):
                        continue
                    key = ("", protein)
                    seen.add(key)
                    ep_status_clean = re.sub(r"\s+", " ", ep_status).lower()
                    significance = {
                        "abnormal": "pathogenic",
                        "wildtype": "benign",
                    }.get(ep_status_clean, "unknown")
                    variants.append(
                        {
                            "gene_symbol": gene_symbol,
                            "cdna_notation": None,
                            "protein_notation": protein,
                            "clinical_significance": significance,
                            "patients": {
                                "count": None,
                                "phenotype": f"{gene_symbol} functionally characterized nsSNV",
                            },
                            "penetrance_data": {
                                "total_carriers_observed": None,
                                "affected_count": None,
                                "unaffected_count": None,
                            },
                            "individual_records": [],
                            "functional_data": {
                                "summary": ep_status_clean,
                                "assays": [],
                            },
                            "segregation_data": None,
                            "population_frequency": None,
                            "evidence_level": "low",
                            "source_location": current_table_label
                            or f"Supplemental {gene_symbol} functional variant table",
                            "additional_notes": "Parsed via deterministic fixed-width functional-status table parser",
                            "key_quotes": [],
                        }
                    )
                    continue

            if (
                current_table_label
                and gene_symbol.upper() in current_table_label.upper()
                and len(tokens) >= 3
                and tokens[0].startswith("#")
            ):
                protein = self._valid_table_protein(tokens[1])
                if not protein:
                    continue
                key = ("", protein)
                if key in seen:
                    continue
                seen.add(key)
                variants.append(
                    {
                        "gene_symbol": gene_symbol,
                        "cdna_notation": None,
                        "protein_notation": protein,
                        "clinical_significance": "pathogenic",
                        "patients": {
                            "count": 1,
                            "phenotype": f"{gene_symbol}-associated disease",
                        },
                        "penetrance_data": {
                            "total_carriers_observed": 1,
                            "affected_count": 1,
                            "unaffected_count": 0,
                        },
                        "individual_records": [],
                        "functional_data": {"summary": "", "assays": []},
                        "segregation_data": None,
                        "population_frequency": None,
                        "evidence_level": "medium",
                        "source_location": current_table_label,
                        "additional_notes": "Parsed via deterministic fixed-width patient table parser",
                        "key_quotes": [],
                    }
                )
                continue

            if len(tokens) >= 8 and tokens[0].startswith("Chr") and ":" in tokens[0]:
                if tokens[3].upper() != gene_symbol.upper():
                    continue
                location = tokens[4]
                cdna_idx = 6 if location.lower() == "exon" else 5
                if cdna_idx >= len(tokens):
                    continue
                cdna_raw = tokens[cdna_idx]
                protein_raw = (
                    tokens[cdna_idx + 1] if cdna_idx + 1 < len(tokens) else None
                )
                cdna = self._valid_table_cdna(cdna_raw)
                protein = self._valid_table_protein(protein_raw)
                if not cdna and not protein:
                    continue

                key = (cdna or "", protein or "")
                if key in seen:
                    continue
                seen.add(key)
                variants.append(
                    {
                        "gene_symbol": gene_symbol,
                        "cdna_notation": cdna,
                        "protein_notation": protein,
                        "clinical_significance": "unknown",
                        "patients": {
                            "count": 1,
                            "phenotype": f"{gene_symbol} variant carrier",
                        },
                        "penetrance_data": {
                            "total_carriers_observed": 1,
                            "affected_count": 1,
                            "unaffected_count": 0,
                        },
                        "individual_records": [],
                        "functional_data": {"summary": "", "assays": []},
                        "segregation_data": None,
                        "population_frequency": None,
                        "evidence_level": "medium",
                        "source_location": "Supplemental Table 1",
                        "additional_notes": "Parsed via deterministic fixed-width WES table parser",
                        "key_quotes": [],
                    }
                )
                continue

            if len(tokens) < 4 or not tokens[0].lower().startswith("c"):
                continue
            if not tokens[-1].isdigit():
                continue

            exon_idx = None
            for idx in range(1, min(len(tokens), 5)):
                if tokens[idx].isdigit():
                    exon_idx = idx
                    break
            if exon_idx is None:
                continue

            cdna_raw = " ".join(tokens[:exon_idx])
            count = int(tokens[-1])
            middle = tokens[exon_idx + 1 : -1]

            for phrase in effect_phrases:
                if (
                    len(middle) >= len(phrase)
                    and tuple(middle[-len(phrase) :]) == phrase
                ):
                    middle = middle[: -len(phrase)]
                    break

            protein_raw = middle[0] if middle else None
            cdna = self._valid_table_cdna(cdna_raw) or self._valid_table_cdna(tokens[0])
            protein = self._valid_table_protein(protein_raw)
            if not cdna and not protein:
                continue

            key = (cdna or "", protein or "")
            if key in seen:
                continue
            seen.add(key)

            variants.append(
                {
                    "gene_symbol": gene_symbol,
                    "cdna_notation": cdna,
                    "protein_notation": protein,
                    "clinical_significance": "pathogenic",
                    "patients": {
                        "count": count,
                        "phenotype": f"{gene_symbol}-associated disease",
                    },
                    "penetrance_data": {
                        "total_carriers_observed": count,
                        "affected_count": count,
                        "unaffected_count": 0,
                    },
                    "individual_records": [],
                    "functional_data": {"summary": "", "assays": []},
                    "segregation_data": None,
                    "population_frequency": None,
                    "evidence_level": "medium",
                    "source_location": "Fixed-width table",
                    "additional_notes": "Parsed via deterministic fixed-width table parser",
                    "key_quotes": [],
                }
            )

        for variant in self._parse_lqts_compendium_vertical_variants(
            full_text, gene_symbol
        ):
            key = (
                variant.get("cdna_notation") or "",
                variant.get("protein_notation") or "",
            )
            if key in seen:
                continue
            seen.add(key)
            variants.append(variant)

        if variants:
            logger.info(
                f"Parsed {len(variants)} variants via deterministic fixed-width table parser"
            )
        return variants

    def _get_extracted_variants_summary(self, variants: list) -> str:
        """Create a compact summary of extracted variants for continuation prompts."""
        summaries = []
        for v in variants:
            cdna = v.get("cdna_notation", "") or ""
            protein = v.get("protein_notation", "") or ""
            summaries.append(f"- {cdna} / {protein}")
        return "\n".join(summaries)

    def _merge_continuation_results(
        self, base_data: dict, continuation_data: dict
    ) -> dict:
        """Merge continuation extraction results into base results."""
        # Add continuation variants
        base_variants = base_data.get("variants", [])
        continuation_variants = continuation_data.get("variants", [])

        # Deduplicate by cdna_notation + protein_notation
        existing_keys = set()
        for v in base_variants:
            key = (v.get("cdna_notation", ""), v.get("protein_notation", ""))
            existing_keys.add(key)

        new_variants = []
        for v in continuation_variants:
            key = (v.get("cdna_notation", ""), v.get("protein_notation", ""))
            if key not in existing_keys:
                new_variants.append(v)
                existing_keys.add(key)

        base_variants.extend(new_variants)
        base_data["variants"] = base_variants

        # Update metadata
        if "extraction_metadata" in base_data:
            base_data["extraction_metadata"]["total_variants_found"] = len(
                base_variants
            )
            continuation_count = continuation_data.get("extraction_metadata", {}).get(
                "continuation_variants_found", len(new_variants)
            )
            base_data["extraction_metadata"]["notes"] = (
                base_data["extraction_metadata"].get("notes", "")
                + f" [Continuation added {continuation_count} variants]"
            ).strip()

        return base_data

    def _attempt_continuation(
        self, paper: Paper, model: str, base_data: dict, full_text: str
    ) -> dict:
        """Attempt to extract remaining variants after truncation."""
        variants = base_data.get("variants", [])
        expected_count = base_data.get("extraction_metadata", {}).get(
            "total_variants_found", len(variants)
        )

        # Only continue if we're missing a significant number of variants
        if (
            len(variants) >= expected_count
            or expected_count - len(variants) < CONTINUATION_VARIANT_GAP
        ):
            return base_data

        logger.info(
            f"PMID {paper.pmid} - Truncation detected: extracted {len(variants)} of {expected_count} variants. "
            f"Attempting continuation extraction."
        )

        prompt = CONTINUATION_PROMPT.format(
            extracted_count=len(variants),
            expected_count=expected_count,
            extracted_variants_list=self._get_extracted_variants_summary(variants),
            gene_symbol=paper.gene_symbol or "UNKNOWN",
            title=paper.title or "Unknown Title",
            full_text=self._truncate_text_for_prompt(
                full_text, gene_symbol=paper.gene_symbol
            ),
        )

        try:
            continuation_data = self.call_llm_json(prompt)
            merged = self._merge_continuation_results(base_data, continuation_data)
            logger.info(
                f"PMID {paper.pmid} - Continuation successful. Total variants now: {len(merged.get('variants', []))}"
            )
            return merged
        except Exception as e:
            logger.warning(f"PMID {paper.pmid} - Continuation extraction failed: {e}")
            return base_data

    def _normalize_stop_codon_notation(self, extracted_data: dict) -> dict:
        """
        Normalize stop codon notation in extracted variants.

        Converts various stop codon representations to the standard '*' notation:
        - 'X' -> '*' (e.g., p.Arg412X -> p.Arg412*)
        - 'Ter' -> '*' (e.g., p.Arg412Ter -> p.Arg412*)
        - 'Stop' -> '*'
        - Handles frameshift notation

        Also handles novel mutation markers (asterisk AFTER the mutation).
        """
        for variant in extracted_data.get("variants", []):
            protein = variant.get("protein_notation", "")
            if not protein:
                continue

            original = protein

            # Convert 'Ter' (termination) at the end of protein notation to '*'
            # e.g., p.Arg412Ter -> p.Arg412*, p.R412Ter -> p.R412*
            protein = re.sub(
                r"([A-Za-z]{1,3}\d+)Ter$", r"\1*", protein, flags=re.IGNORECASE
            )

            # Convert 'Stop' at the end to '*'
            protein = re.sub(
                r"([A-Za-z]{1,3}\d+)Stop$", r"\1*", protein, flags=re.IGNORECASE
            )

            # Convert X at the end of protein notation to '*' (stop codon)
            # e.g., p.Arg412X -> p.Arg412*, p.R412X -> p.R412*
            protein = re.sub(r"([A-Za-z]{1,3}\d+)X$", r"\1*", protein)

            # Pattern for frameshift with X or Ter: fs + digits + X/Ter
            # e.g., p.Gly24fs+34X -> p.Gly24fs*34, p.Gly24fsTer58 -> p.Gly24fs*58
            protein = re.sub(r"(fs\+?)(\d*)X$", r"\1*\2", protein)
            protein = re.sub(
                r"(fs\+?)(\d*)Ter$", r"\1*\2", protein, flags=re.IGNORECASE
            )
            protein = re.sub(r"(fs)\*(\d+)X$", r"\1*\2", protein)

            # Handle cases where frameshift shows as fs*NUMBER or fsX NUMBER
            # Normalize fs*58 format (already correct but might have extra chars)
            protein = re.sub(r"(fs)\s*\*\s*(\d+)", r"\1*\2", protein)

            # Handle truncating mutations that end with unusual patterns
            # Sometimes 'stop gained' is written as the amino acid that replaces (incorrectly)
            # This is harder to detect without context, but we can flag suspicious patterns

            # If notation changed, update it
            if protein != original:
                variant["protein_notation"] = protein
                notes = variant.get("additional_notes", "") or ""
                if notes:
                    notes += " "
                notes += f"[Stop codon notation normalized: {original} -> {protein}]"
                variant["additional_notes"] = notes

        return extracted_data

    def _populate_penetrance_from_patient_count(self, extracted_data: dict) -> dict:
        """
        Populate penetrance_data fields from patients.count when appropriate.

        For disease-associated mutation tables (like "LQT2-associated mutations"),
        the patient count represents affected carriers. This function maps:
        - patients.count -> penetrance_data.total_carriers_observed
        - patients.count -> penetrance_data.affected_count (for pathogenic variants)

        Only populates when penetrance_data fields are null/missing but patients.count exists.
        """
        for variant in extracted_data.get("variants", []):
            patients = variant.get("patients", {})
            patient_count = patients.get("count") if patients else None

            # Skip if no patient count
            if patient_count is None or patient_count == 0:
                continue

            # Get or create penetrance_data
            pdata = variant.get("penetrance_data", {})
            if pdata is None:
                pdata = {}
                variant["penetrance_data"] = pdata

            # Only populate if fields are missing/null
            total_carriers = pdata.get("total_carriers_observed")
            affected_count = pdata.get("affected_count")

            if total_carriers is None and affected_count is None:
                # Determine if this is a pathogenic/disease-associated variant
                significance = (variant.get("clinical_significance") or "").lower()
                phenotype = (patients.get("phenotype") or "").lower()
                source_location = (variant.get("source_location") or "").lower()

                # Check if from a disease-associated context
                is_disease_associated = any(
                    [
                        "pathogenic" in significance,
                        "lqt" in phenotype,  # Long QT syndrome
                        "brugada" in phenotype,
                        "arrhythmia" in phenotype,
                        "syndrome" in phenotype,
                        "disease" in phenotype,
                        "affected" in phenotype,
                        "mutation" in source_location,
                        "lqt" in source_location,
                    ]
                )

                # For disease-associated variants, patient count = affected count
                if is_disease_associated:
                    pdata["total_carriers_observed"] = patient_count
                    pdata["affected_count"] = patient_count
                    pdata["unaffected_count"] = 0

                    # Add note about the mapping
                    notes = variant.get("additional_notes", "") or ""
                    if notes:
                        notes += " "
                    notes += f"[Penetrance data inferred from patient count: {patient_count} affected carriers]"
                    variant["additional_notes"] = notes
                else:
                    # For uncertain significance, just set total carriers
                    pdata["total_carriers_observed"] = patient_count

        return extracted_data

    @staticmethod
    def _coerce_count(value) -> int | None:
        if value is None or value == "":
            return None
        if isinstance(value, bool):
            return None
        if isinstance(value, (int, float)):
            return int(value)
        match = re.search(r"\d+", str(value))
        return int(match.group()) if match else None

    @staticmethod
    def _count_source_is_deterministic(variant: dict) -> bool:
        text = " ".join(
            str(variant.get(key) or "")
            for key in ("additional_notes", "source_location", "extraction_source")
        ).lower()
        return any(
            marker in text
            for marker in (
                "deterministic",
                "table regex",
                "fixed-width",
                "router-first",
            )
        )

    @staticmethod
    def _infer_source_layer(variant: dict) -> str:
        return infer_source_layer_from_text(
            source_location=variant.get("source_location"),
            additional_notes=variant.get("additional_notes"),
            extraction_source=variant.get("extraction_source"),
            source_layer=variant.get("source_layer"),
        )

    def _annotate_source_layers(self, extracted_data: dict | None) -> dict | None:
        if not extracted_data:
            return extracted_data
        for variant in extracted_data.get("variants", []) or []:
            if isinstance(variant, dict):
                variant["source_layer"] = self._infer_source_layer(variant)
        return extracted_data

    def _suppress_repeated_study_wide_counts(self, extracted_data: dict) -> dict:
        """
        Clear likely cohort-wide totals copied onto every variant by the LLM.

        A recurring recall failure is an abstract/table statement such as
        "43 carriers: 28 affected, 15 unaffected" getting repeated on each
        variant in the paper. Row-level deterministic parsers are protected;
        this only acts on non-deterministic variants with the same moderate/large
        count tuple repeated across several variants.
        """
        groups: dict[tuple[int, int | None, int | None], list[dict]] = {}
        for variant in extracted_data.get("variants", []):
            if self._count_source_is_deterministic(variant):
                continue
            patients = variant.get("patients") or {}
            pdata = variant.get("penetrance_data") or {}
            total = self._coerce_count(
                pdata.get("total_carriers_observed", patients.get("count"))
            )
            affected = self._coerce_count(pdata.get("affected_count"))
            unaffected = self._coerce_count(pdata.get("unaffected_count"))
            if total is None or total < 10:
                continue
            groups.setdefault((total, affected, unaffected), []).append(variant)

        suppressed = 0
        for (total, affected, unaffected), variants in groups.items():
            if len(variants) < 3 or total < len(variants) * 2:
                continue
            for variant in variants:
                patients = variant.setdefault("patients", {})
                pdata = variant.setdefault("penetrance_data", {})
                if self._coerce_count(patients.get("count")) == total:
                    patients["count"] = None
                if self._coerce_count(pdata.get("total_carriers_observed")) == total:
                    pdata["total_carriers_observed"] = None
                if (
                    affected is not None
                    and self._coerce_count(pdata.get("affected_count")) == affected
                ):
                    pdata["affected_count"] = None
                if (
                    unaffected is not None
                    and self._coerce_count(pdata.get("unaffected_count")) == unaffected
                ):
                    pdata["unaffected_count"] = None
                if self._coerce_count(pdata.get("uncertain_count")) == total:
                    pdata["uncertain_count"] = None
                notes = variant.get("additional_notes", "") or ""
                if notes:
                    notes += " "
                notes += (
                    "[Cleared repeated cohort-wide count tuple copied across "
                    f"{len(variants)} variants: total={total}, affected={affected}, "
                    f"unaffected={unaffected}]"
                )
                variant["additional_notes"] = notes
                suppressed += 1

        if suppressed and "extraction_metadata" in extracted_data:
            extracted_data["extraction_metadata"]["study_wide_count_suppressed"] = (
                suppressed
            )
        return extracted_data

    def _suppress_repeated_study_wide_context_fields(
        self, extracted_data: dict
    ) -> dict:
        """Clear cohort-summary phenotype/demographics broadcast to many rows."""
        study_wide_markers = re.compile(
            r"\b(overall|cohort|study population|entire|all carriers|"
            r"designated variant carriers|median|mean|average|percent)\b|%",
            re.IGNORECASE,
        )
        groups: dict[tuple[str, str], list[dict]] = {}
        for variant in extracted_data.get("variants", []):
            if self._count_source_is_deterministic(variant):
                continue
            patients = variant.get("patients") or {}
            for field in ("phenotype", "demographics"):
                value = patients.get(field)
                if not isinstance(value, str):
                    continue
                normalized = re.sub(r"\s+", " ", value).strip()
                if len(normalized) < 16 or not study_wide_markers.search(normalized):
                    continue
                groups.setdefault((field, normalized.lower()), []).append(variant)

        suppressed = 0
        for (field, value), variants in groups.items():
            if len(variants) < 5:
                continue
            for variant in variants:
                patients = variant.setdefault("patients", {})
                current = patients.get(field)
                if not isinstance(current, str):
                    continue
                if re.sub(r"\s+", " ", current).strip().lower() != value:
                    continue
                patients[field] = None
                notes = variant.get("additional_notes", "") or ""
                if notes:
                    notes += " "
                notes += (
                    f"[Cleared repeated cohort-wide {field} copied across "
                    f"{len(variants)} variants]"
                )
                variant["additional_notes"] = notes
                suppressed += 1

        if suppressed and "extraction_metadata" in extracted_data:
            extracted_data["extraction_metadata"]["study_wide_context_suppressed"] = (
                suppressed
            )
        return extracted_data

    @staticmethod
    def _source_has_missing_table_bodies(text: str) -> bool:
        table_refs = len(re.findall(r"\bTables?\s+\d+(?:\s*(?:and|,|-)\s*\d+)?", text))
        table_body_markers = len(
            re.findall(r"(?m)^\s*\|.+\|\s*$", text)
            + re.findall(r"(?m)^\s*(?:#{1,6}\s*)?(?:Table|TABLE)\s+\d+[.:]", text)
            + re.findall(r"(?i)<table\b", text)
        )
        return table_refs >= 3 and table_body_markers == 0 and len(text) > 5_000

    @staticmethod
    def _compact_variant_for_evidence(variant: dict) -> dict[str, Any]:
        return {
            "gene_symbol": variant.get("gene_symbol"),
            "cdna_notation": variant.get("cdna_notation"),
            "protein_notation": variant.get("protein_notation"),
            "genomic_position": variant.get("genomic_position"),
            "clinical_significance": variant.get("clinical_significance"),
            "patients": variant.get("patients"),
            "penetrance_data": variant.get("penetrance_data"),
            "source_location": variant.get("source_location"),
            "additional_notes": variant.get("additional_notes"),
        }

    def _count_tuple(
        self, variant: dict
    ) -> tuple[int | None, int | None, int | None, int | None]:
        patients = variant.get("patients") or {}
        pdata = variant.get("penetrance_data") or {}
        total = self._coerce_count(
            pdata.get("total_carriers_observed", patients.get("count"))
        )
        affected = self._coerce_count(pdata.get("affected_count"))
        unaffected = self._coerce_count(pdata.get("unaffected_count"))
        uncertain = self._coerce_count(pdata.get("uncertain_count"))
        return total, affected, unaffected, uncertain

    def _assess_extraction_risk(
        self,
        *,
        extracted_data: dict,
        source_text: str,
        estimated_variants: int | None = None,
        scanner_variant_count: int = 0,
    ) -> dict[str, Any]:
        """Score whether an extraction should get second-model adjudication."""
        variants = extracted_data.get("variants") or []
        if not isinstance(variants, list):
            variants = []

        reasons: list[str] = []
        source_blockers: list[str] = []
        score = 0

        if self._source_has_missing_table_bodies(source_text):
            source_blockers.append("source_mentions_tables_but_no_table_bodies")

        if estimated_variants and estimated_variants >= 8:
            if len(variants) < max(2, estimated_variants // 4):
                score += 2
                reasons.append(
                    f"low_variant_yield_vs_table_hint:{len(variants)}/{estimated_variants}"
                )

        if scanner_variant_count >= 5 and len(variants) < max(
            2, scanner_variant_count // 2
        ):
            score += 1
            reasons.append(
                f"low_variant_yield_vs_scanner:{len(variants)}/{scanner_variant_count}"
            )

        repeated_groups: dict[tuple[int, int | None, int | None], int] = {}
        arithmetic_mismatches = 0
        variants_with_counts = 0
        for variant in variants:
            total, affected, unaffected, uncertain = self._count_tuple(variant)
            if any(
                value is not None for value in (total, affected, unaffected, uncertain)
            ):
                variants_with_counts += 1
            if (
                total is not None
                and affected is not None
                and unaffected is not None
                and affected + unaffected + (uncertain or 0) != total
            ):
                arithmetic_mismatches += 1
            if (
                total is not None
                and total >= 10
                and not self._count_source_is_deterministic(variant)
            ):
                repeated_groups[(total, affected, unaffected)] = (
                    repeated_groups.get((total, affected, unaffected), 0) + 1
                )

        repeated = [
            (count_tuple, count)
            for count_tuple, count in repeated_groups.items()
            if count >= 3 and count_tuple[0] >= count * 2
        ]
        if repeated:
            score += 2
            reasons.append(
                "repeated_large_count_tuple:"
                + ";".join(
                    f"{count_tuple}x{count}" for count_tuple, count in repeated[:3]
                )
            )

        if arithmetic_mismatches:
            score += 2
            reasons.append(f"count_arithmetic_mismatch:{arithmetic_mismatches}")

        meta = extracted_data.get("extraction_metadata") or {}
        if meta.get("study_wide_count_suppressed"):
            score += 2
            reasons.append(
                f"study_wide_count_suppressed:{meta['study_wide_count_suppressed']}"
            )

        if len(variants) >= 5 and variants_with_counts <= len(variants) // 3:
            score += 1
            reasons.append(
                f"many_variants_missing_counts:{variants_with_counts}/{len(variants)}"
            )

        return {
            "score": score,
            "reasons": reasons,
            "source_blockers": source_blockers,
            "estimated_variants": estimated_variants,
            "scanner_variant_count": scanner_variant_count,
            "variant_count": len(variants),
            "requires_adjudication": score >= self.adjudication_risk_threshold,
        }

    def _select_evidence_lines(
        self,
        *,
        source_text: str,
        gene_symbol: str | None,
        variants: list[dict],
        max_chars: int,
    ) -> str:
        variant_terms = set()
        for variant in variants:
            for key in ("protein_notation", "cdna_notation", "genomic_position"):
                value = str(variant.get(key) or "").strip()
                if value:
                    variant_terms.add(value)
                    if value.startswith("p."):
                        variant_terms.add(value[2:])

        gene = (gene_symbol or "").strip().lower()
        count_terms = (
            "carrier",
            "carriers",
            "affected",
            "unaffected",
            "asymptomatic",
            "symptomatic",
            "patient",
            "patients",
            "proband",
            "probands",
            "control",
            "controls",
            "table",
            "mutation",
            "variant",
        )
        selected: list[str] = []
        seen: set[int] = set()
        lines = source_text.splitlines()
        for idx, line in enumerate(lines):
            lower = line.lower()
            hit = (gene and gene in lower) or any(
                term.lower() in lower for term in variant_terms
            )
            hit = hit or any(term in lower for term in count_terms)
            if not hit:
                continue
            for nearby in range(max(0, idx - 1), min(len(lines), idx + 2)):
                if nearby in seen:
                    continue
                seen.add(nearby)
                selected.append(f"L{nearby + 1}: {lines[nearby][:500]}")
            if sum(len(item) + 1 for item in selected) >= max_chars:
                break
        return "\n".join(selected)[:max_chars]

    def _build_adjudication_prompt(
        self,
        *,
        paper: Paper,
        primary_model: str,
        extracted_data: dict,
        source_text: str,
        risk: dict[str, Any],
    ) -> str:
        variants = extracted_data.get("variants") or []
        compact_variants = [
            self._compact_variant_for_evidence(variant)
            for variant in variants[:150]
            if isinstance(variant, dict)
        ]
        primary_json = json.dumps(compact_variants, ensure_ascii=False, indent=2)
        evidence_budget = max(
            4_000,
            self.evidence_packet_max_chars - min(len(primary_json), 12_000) - 4_000,
        )
        evidence_lines = self._select_evidence_lines(
            source_text=source_text,
            gene_symbol=paper.gene_symbol,
            variants=variants,
            max_chars=evidence_budget,
        )

        prompt = f"""You are a careful biomedical variant-extraction adjudicator.

You receive a compact evidence packet, not the full paper. Your job is to correct
the primary extraction only when the provided evidence supports the correction.
Do not invent counts. Do not copy study-wide totals onto each variant. If a count
is aggregate for the whole study, family set, table, domain, or mutation class,
set the variant-level count fields to null and explain that in additional_notes.

Target gene: {paper.gene_symbol or "UNKNOWN"}
PMID: {paper.pmid}
Title: {paper.title or "Unknown Title"}
Primary extraction model: {primary_model}
Risk assessment: {json.dumps(risk, ensure_ascii=False)}

Primary extracted variants:
{primary_json[:12000]}

Evidence snippets:
{evidence_lines}

Return strict JSON with this schema:
{{
  "variants": [
    {{
      "gene_symbol": "{paper.gene_symbol or "UNKNOWN"}",
      "cdna_notation": "string or null",
      "protein_notation": "string or null",
      "genomic_position": "string or null",
      "clinical_significance": "string or null",
      "patients": {{"count": "integer or null", "phenotype": "string or null", "demographics": "string or null"}},
      "penetrance_data": {{
        "total_carriers_observed": "integer or null",
        "affected_count": "integer or null",
        "unaffected_count": "integer or null",
        "uncertain_count": "integer or null",
        "penetrance_percentage": "float or null",
        "age_dependent_penetrance": []
      }},
      "source_location": "specific table/line/section if supported",
      "additional_notes": "brief evidence rationale, including why counts are null"
    }}
  ],
  "extraction_metadata": {{
    "adjudication_confidence": "high/medium/low",
    "adjudication_notes": "brief summary of changes and unresolved evidence gaps"
  }}
}}
"""
        return prompt[: self.evidence_packet_max_chars]

    def _adjudicator_candidates(self, primary_model: str) -> list[str]:
        seen = {primary_model}
        out: list[str] = []
        for model in [*self.adjudicator_models, *self.models]:
            if not model or model in seen:
                continue
            out.append(model)
            seen.add(model)
        return out

    def _call_adjudicator(self, model: str, prompt: str) -> dict[str, Any]:
        saved_model, saved_max_tokens, saved_temperature = (
            self.model,
            self.max_tokens,
            self.temperature,
        )
        try:
            self.model = model
            self.max_tokens = self._clamp_max_tokens(
                model, self.adjudication_max_tokens
            )
            self.temperature = 0.0
            return self.call_llm_json(prompt)
        finally:
            self.model = saved_model
            self.max_tokens = saved_max_tokens
            self.temperature = saved_temperature

    def _merge_adjudicated_extraction(
        self,
        *,
        extracted_data: dict,
        adjudicated_data: dict[str, Any],
        adjudicator_model: str,
        risk: dict[str, Any],
    ) -> dict:
        adjudicated_variants = adjudicated_data.get("variants")
        if not isinstance(adjudicated_variants, list) or not adjudicated_variants:
            raise ValueError("Adjudicator returned no variants")

        merged = copy.deepcopy(extracted_data)
        merged["variants"] = adjudicated_variants
        metadata = merged.setdefault("extraction_metadata", {})
        metadata.update(adjudicated_data.get("extraction_metadata") or {})
        metadata["total_variants_found"] = len(adjudicated_variants)
        metadata["adjudication_applied"] = True
        metadata["adjudication_model"] = adjudicator_model
        metadata["adjudication_risk"] = risk
        return merged

    @staticmethod
    def _result_model_label(primary_model: str, extracted_data: dict) -> str:
        metadata = extracted_data.get("extraction_metadata", {})
        if metadata.get("claim_verification_applied"):
            return (
                f"{primary_model}+verified:{metadata.get('claim_verification_model')}"
            )
        if not metadata.get("adjudication_applied"):
            return primary_model
        return f"{primary_model}+adjudicated:{metadata.get('adjudication_model')}"

    def _verify_claim_cards_for_extraction(
        self,
        *,
        paper: Paper,
        primary_model: str,
        extracted_data: dict,
        source_text: str,
        verifier_model: str,
    ) -> dict:
        from pipeline.claim_verifier import (
            VariantClaimVerifier,
            apply_verification_to_variant,
            build_claim_card,
        )

        variants = extracted_data.get("variants") or []
        if not isinstance(variants, list) or not variants:
            raise ValueError("No variants available for claim verification")

        verifier = VariantClaimVerifier(
            model=verifier_model,
            max_tokens=min(self.adjudication_max_tokens, 2500),
            reasoning_effort=self.reasoning_effort,
        )
        updated_variants: list[dict] = []
        verification_results: list[dict[str, Any]] = []
        field_changes: list[dict[str, Any]] = []
        verified_count = 0

        for idx, variant in enumerate(variants):
            if verified_count >= self.max_verifier_cards:
                updated_variants.append(variant)
                continue
            card = build_claim_card(
                source_text=source_text,
                gene=paper.gene_symbol or variant.get("gene_symbol") or "UNKNOWN",
                disease=paper.disease,
                pmid=paper.pmid,
                title=paper.title,
                variant=variant,
                max_evidence_chars=max(2_000, self.evidence_packet_max_chars // 6),
            )
            if card is None or not card.evidence.strip():
                updated = dict(variant)
                updated["claim_verification"] = {
                    "verdict": "source_missing",
                    "field_verdicts": {
                        "variant": "ambiguous",
                        "total_carriers": "source_missing",
                        "affected": "source_missing",
                        "unaffected": "source_missing",
                    },
                    "corrected_values": {
                        "total_carriers": None,
                        "affected": None,
                        "unaffected": None,
                    },
                    "reason": "No compact evidence lines found for this claim.",
                    "evidence_quote": "",
                }
                updated_variants.append(updated)
                continue

            verification = verifier.verify(card)
            updated, changes = apply_verification_to_variant(variant, verification)
            updated_variants.append(updated)
            verified_count += 1
            verification_results.append(
                {
                    "variant_index": idx,
                    "variant": card.variant,
                    "extracted": card.extracted,
                    "verification": verification,
                }
            )
            if changes:
                field_changes.append(
                    {
                        "variant_index": idx,
                        "variant": card.variant,
                        "changes": changes,
                    }
                )

        merged = copy.deepcopy(extracted_data)
        merged["variants"] = updated_variants
        metadata = merged.setdefault("extraction_metadata", {})
        metadata["claim_verification_applied"] = True
        metadata["claim_verification_model"] = verifier_model
        metadata["claim_verification_primary_model"] = primary_model
        metadata["claim_verification_cards"] = verified_count
        metadata["claim_verification_field_changes"] = field_changes
        metadata["claim_verification_results"] = verification_results[
            : self.max_verifier_cards
        ]
        return merged

    def _maybe_adjudicate_extraction(
        self,
        *,
        paper: Paper,
        primary_model: str,
        extracted_data: dict,
        source_text: str,
        estimated_variants: int | None,
        scanner_variant_count: int,
    ) -> dict:
        risk = self._assess_extraction_risk(
            extracted_data=extracted_data,
            source_text=source_text,
            estimated_variants=estimated_variants,
            scanner_variant_count=scanner_variant_count,
        )
        metadata = extracted_data.setdefault("extraction_metadata", {})
        metadata["extraction_risk"] = risk

        if not self.enable_ensemble_qa:
            metadata["adjudication_skipped_reason"] = "ensemble_qa_disabled"
            return extracted_data
        if not risk["requires_adjudication"]:
            metadata["adjudication_skipped_reason"] = "risk_below_threshold"
            return extracted_data
        if risk["source_blockers"] and not risk["reasons"]:
            metadata["adjudication_skipped_reason"] = (
                "source_blocker_without_extraction_risk"
            )
            return extracted_data

        candidates = self._adjudicator_candidates(primary_model)
        if not candidates:
            metadata["adjudication_skipped_reason"] = "no_adjudicator_model_available"
            return extracted_data

        claim_errors: list[str] = []
        for verifier_model in candidates[:1]:
            try:
                return self._verify_claim_cards_for_extraction(
                    paper=paper,
                    primary_model=primary_model,
                    extracted_data=extracted_data,
                    source_text=source_text,
                    verifier_model=verifier_model,
                )
            except Exception as exc:  # noqa: BLE001
                claim_errors.append(f"{verifier_model}: {exc}")
                logger.warning(
                    "PMID %s - claim verification with %s failed: %s",
                    paper.pmid,
                    verifier_model,
                    exc,
                )

        prompt = self._build_adjudication_prompt(
            paper=paper,
            primary_model=primary_model,
            extracted_data=extracted_data,
            source_text=source_text,
            risk=risk,
        )
        errors: list[str] = []
        for adjudicator_model in candidates[:1]:
            try:
                adjudicated = self._call_adjudicator(adjudicator_model, prompt)
                return self._merge_adjudicated_extraction(
                    extracted_data=extracted_data,
                    adjudicated_data=adjudicated,
                    adjudicator_model=adjudicator_model,
                    risk=risk,
                )
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{adjudicator_model}: {exc}")
                logger.warning(
                    "PMID %s - adjudication with %s failed: %s",
                    paper.pmid,
                    adjudicator_model,
                    exc,
                )

        metadata["adjudication_skipped_reason"] = "adjudication_failed"
        metadata["adjudication_errors"] = [*claim_errors, *errors]
        return extracted_data

    # Patterns that indicate failed content extraction
    FAILED_EXTRACTION_PATTERNS = [
        "[PDF file available at:",
        "[Error converting",
        "[Error reading",
        "[Legacy .doc file available at:",
        "text extraction failed",
        "manual review required",
        "[NO TEXT AVAILABLE]",
        "[ABSTRACT ONLY",
        "Access Denied",
        "403 Forbidden",
        "404 Not Found",
        "Page not found",
        "Unable to retrieve",
        "Subscription required",
        "Please login",
        "Session expired",
    ]

    # Patterns indicating HTML/markup garbage
    HTML_GARBAGE_PATTERNS = [
        r"<(?:div|span|script|style|html|body|head)[^>]*>",
        r"</(?:div|span|script|style|html|body|head)>",
        r"\{[\s]*[\"\']?[a-zA-Z_]+[\"\']?\s*:",  # JSON-like fragments
        r"class=[\"\'][^\"\']+[\"\']",
        r"style=[\"\'][^\"\']+[\"\']",
    ]

    def _assess_input_quality(
        self, text: str, gene_symbol: Optional[str]
    ) -> tuple[bool, str]:
        """
        Assess whether the input text is of sufficient quality for extraction.

        This is a circuit breaker that prevents wasting LLM calls on garbage input.
        Checks for:
        - Minimum content length
        - Failed extraction placeholders
        - HTML/markup garbage
        - Sufficient alphanumeric content ratio
        - Relevant variant or gene content

        Returns:
            (is_usable, reason) - whether text is usable and why/why not
        """
        # Check 1: Minimum length threshold
        if not text:
            return False, "No text provided"

        text_len = len(text)
        if text_len < MIN_EXTRACTION_INPUT_SIZE:
            return (
                False,
                f"Text too short ({text_len} chars, min={MIN_EXTRACTION_INPUT_SIZE})",
            )

        # Check 2: Failed extraction placeholders
        failed_count = sum(
            1 for pattern in self.FAILED_EXTRACTION_PATTERNS if pattern in text
        )
        if failed_count >= 2:
            return (
                False,
                f"Contains {failed_count} failed extraction markers",
            )

        # Check 3: HTML/markup garbage detection
        html_matches = sum(
            1 for pattern in self.HTML_GARBAGE_PATTERNS if re.search(pattern, text)
        )
        if html_matches >= 3:
            return False, f"Contains HTML/markup garbage ({html_matches} patterns)"

        # Check 4: Alphanumeric content ratio
        # Filters out text that's mostly punctuation, whitespace, or special chars
        alphanumeric_chars = sum(1 for c in text if c.isalnum() or c.isspace())
        ratio = alphanumeric_chars / text_len if text_len > 0 else 0
        if ratio < MIN_ALPHANUMERIC_RATIO:
            return (
                False,
                f"Low alphanumeric ratio ({ratio:.2f}, min={MIN_ALPHANUMERIC_RATIO})",
            )

        # Check 5: Variant-like or gene content
        variant_patterns = [
            r"c\.\d+",
            r"p\.[A-Z]",
            r"[A-Z]\d+[A-Z]",  # Basic variant patterns
            r"mutation",
            r"variant",
            r"carrier",  # Clinical terms
        ]
        has_variant_content = any(
            re.search(p, text, re.IGNORECASE) for p in variant_patterns
        )

        if not has_variant_content:
            # Check if gene is at least mentioned
            if gene_symbol and gene_symbol.upper() in text.upper():
                return True, "Gene mentioned but no variant patterns found"
            return False, "No variant patterns or gene mentions in text"

        return True, "Text appears usable"

    def _try_table_router(
        self, paper: Paper, scanner_text: str
    ) -> Optional[ExtractionResult]:
        """Router-first extraction. Returns None to fall through to full-text Tier 3.

        Caller is responsible for the feature flag + gene-symbol guard. We only
        return a successful ExtractionResult if the router found ≥ 1 variant
        across the routed tables — anything less means the paper is either
        narrative-only or has tables the router rejected, and the existing
        full-text path is a better answer.
        """
        from pipeline.table_router import extract_via_router

        settings = get_settings()
        router_model = settings.get_table_router_model()
        try:
            outcome = extract_via_router(
                scanner_text,
                paper.gene_symbol or "UNKNOWN",
                model=router_model,
                max_tokens=settings.table_router_max_tokens,
                reasoning_effort=settings.table_router_reasoning_effort,
            )
        except Exception as e:  # noqa: BLE001
            logger.warning(
                "PMID %s - Table router crashed (%s); falling back to full-text path",
                paper.pmid,
                e,
            )
            return None

        variants = outcome.get("variants") or []
        if outcome.get("error"):
            logger.info(
                "PMID %s - Router LLM error (%s); falling back",
                paper.pmid,
                outcome["error"],
            )
            return None
        if not variants:
            logger.info(
                "PMID %s - Router found no usable variant tables (%d candidates routed); "
                "falling back to full-text Tier 3",
                paper.pmid,
                len(outcome.get("routed", []) or []),
            )
            return None

        routed = outcome.get("routed", []) or []
        table_ids = ", ".join(r.table_id for r in routed)
        logger.info(
            "PMID %s - Router approved tables [%s]; deterministically parsed %d variants",
            paper.pmid,
            table_ids,
            len(variants),
        )

        extracted_data = {
            "paper_metadata": {
                "pmid": paper.pmid,
                "title": paper.title or "Unknown Title",
                "first_author": (paper.authors or [None])[0],
                "journal": paper.journal,
                "publication_date": paper.publication_date,
                "doi": paper.doi,
                "pmc_id": paper.pmc_id,
                "extraction_summary": (
                    f"Router+deterministic parse of {len(variants)} variants from "
                    f"{len(routed)} routed table(s)"
                ),
            },
            "variants": variants,
            "extraction_metadata": {
                "total_variants_found": len(variants),
                "extraction_confidence": "medium",
                "compact_mode": True,
                "notes": (
                    f"Router-first path: tables {table_ids} parsed deterministically "
                    f"using {router_model}"
                ),
            },
        }
        return ExtractionResult(
            pmid=paper.pmid,
            success=True,
            extracted_data=extracted_data,
            model_used=f"router+{router_model}",
        )

    def _attempt_extraction(
        self,
        paper: Paper,
        model: str,
        prepared_full_text: Optional[str] = None,
        estimated_variants: Optional[int] = None,
    ) -> ExtractionResult:
        """Attempt extraction with a single model, with continuation for truncated responses.

        Uses compact extraction mode for papers with many variants to avoid output truncation.
        """
        logger.info(f"PMID {paper.pmid} - Starting expert extraction with {model}")
        saved_model, saved_max_tokens = self.model, self.max_tokens
        try:
            self.model = model
            self.max_tokens = self._clamp_max_tokens(model, self.requested_max_tokens)

            result = self._do_attempt_extraction(
                paper, model, prepared_full_text, estimated_variants
            )
            if result.success:
                result.extracted_data = self._backfill_variant_notation_pairs(
                    result.extracted_data
                )
                result.extracted_data = self._annotate_source_layers(
                    result.extracted_data
                )
            return result
        finally:
            self.model, self.max_tokens = saved_model, saved_max_tokens

    def _do_attempt_extraction(
        self,
        paper: Paper,
        model: str,
        prepared_full_text: Optional[str] = None,
        estimated_variants: Optional[int] = None,
    ) -> ExtractionResult:
        """Inner extraction logic called with self.model/self.max_tokens already set."""
        full_text = (
            prepared_full_text
            if prepared_full_text is not None
            else self._prepare_full_text(paper)
        )
        if full_text == "[NO TEXT AVAILABLE]":
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error="No text available",
                model_used=model,
            )

        # SCANNER TEXT: Always use the original full text for regex scanning
        # and table extraction, even when the LLM gets condensed DATA_ZONES.
        # The scanner is pure regex (no API cost), so it should see everything.
        # This prevents table data from being lost during Scout condensation.
        scanner_text = paper.full_text if paper.full_text else full_text
        raw_scanner_text = scanner_text
        scanner_text = self._augment_pdf_linearized_tables(scanner_text)
        if scanner_text != raw_scanner_text:
            augmented_table_rows = self._estimate_table_rows(scanner_text)
            if estimated_variants is None or augmented_table_rows > estimated_variants:
                estimated_variants = augmented_table_rows
            reconstructed_appendix = scanner_text[len(raw_scanner_text) :].lstrip()
            if full_text == raw_scanner_text:
                full_text = scanner_text
            elif (
                reconstructed_appendix
                and reconstructed_appendix.strip() not in full_text
            ):
                full_text = f"{full_text.rstrip()}\n\n{reconstructed_appendix}"

        # Assess input quality before sending to LLM
        is_usable, quality_reason = self._assess_input_quality(
            full_text, paper.gene_symbol
        )
        if not is_usable:
            logger.warning(
                f"PMID {paper.pmid} - Input quality check failed: {quality_reason}"
            )
            print(f"⚠ PMID {paper.pmid}: Skipping LLM extraction - {quality_reason}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"Input quality insufficient: {quality_reason}",
                model_used=model,
            )

        # Estimate variant count if not provided
        # Use scanner_text (original full text) for table row estimation
        if estimated_variants is None:
            estimated_variants = self._estimate_table_rows(scanner_text)

        # Fast path: deterministic table parser for very large tables to avoid slow LLM calls
        deterministic_min_variants = (
            1 if self.tier_threshold <= 0 else DETERMINISTIC_PARSER_MIN_VARIANTS
        )
        if estimated_variants >= LARGE_TABLE_ROW_THRESHOLD:
            parsed_variants = self._parse_markdown_table_variants(
                scanner_text, paper.gene_symbol
            )
            if len(parsed_variants) >= deterministic_min_variants:
                extracted_data = {
                    "paper_metadata": {
                        "pmid": paper.pmid,
                        "title": paper.title or "Unknown Title",
                        "first_author": (paper.authors or [None])[0],
                        "journal": paper.journal,
                        "publication_date": paper.publication_date,
                        "doi": paper.doi,
                        "pmc_id": paper.pmc_id,
                        "extraction_summary": f"Deterministic table parse of {len(parsed_variants)} variants",
                    },
                    "variants": parsed_variants,
                    "extraction_metadata": {
                        "total_variants_found": len(parsed_variants),
                        "extraction_confidence": "medium",
                        "compact_mode": True,
                        "notes": "Bypassed LLM using markdown table parser for large table",
                    },
                }
                return ExtractionResult(
                    pmid=paper.pmid,
                    success=True,
                    extracted_data=extracted_data,
                    model_used="deterministic-table-parser",
                )

        fixed_width_min_variants = (
            1
            if self.tier_threshold <= 0
            else min(DETERMINISTIC_PARSER_MIN_VARIANTS, 20)
        )
        fixed_width_variants = self._parse_fixed_width_table_variants(
            scanner_text, paper.gene_symbol
        )
        if fixed_width_variants:
            logger.info(
                "deterministic_fixed_width_parser_hit pmid=%s gene=%s variants=%d",
                paper.pmid,
                paper.gene_symbol,
                len(fixed_width_variants),
            )
        fixed_width_has_lqts_compendium = any(
            "deterministic LQTS compendium table parser"
            in (variant.get("additional_notes") or "")
            for variant in fixed_width_variants
        )
        vertical_variants = (
            []
            if fixed_width_has_lqts_compendium
            else self._parse_vertical_gene_table_variants(
                scanner_text, paper.gene_symbol
            )
        )
        deterministic_variants = list(fixed_width_variants)
        deterministic_keys = {
            (
                (variant.get("cdna_notation") or "").lower(),
                (variant.get("protein_notation") or "").lower(),
            )
            for variant in deterministic_variants
        }
        for variant in vertical_variants:
            key = (
                (variant.get("cdna_notation") or "").lower(),
                (variant.get("protein_notation") or "").lower(),
            )
            if key not in deterministic_keys:
                deterministic_variants.append(variant)
                deterministic_keys.add(key)

        if len(deterministic_variants) >= fixed_width_min_variants:
            extracted_data = {
                "paper_metadata": {
                    "pmid": paper.pmid,
                    "title": paper.title or "Unknown Title",
                    "first_author": (paper.authors or [None])[0],
                    "journal": paper.journal,
                    "publication_date": paper.publication_date,
                    "doi": paper.doi,
                    "pmc_id": paper.pmc_id,
                    "extraction_summary": f"Deterministic table parse of {len(deterministic_variants)} variants",
                },
                "variants": deterministic_variants,
                "extraction_metadata": {
                    "total_variants_found": len(deterministic_variants),
                    "extraction_confidence": "medium",
                    "compact_mode": True,
                    "notes": "Bypassed LLM using deterministic parser for large PDF/vertical table",
                    "deterministic_parser_counts": {
                        "fixed_width": len(fixed_width_variants),
                        "vertical": len(vertical_variants),
                    },
                },
            }
            parser_model = (
                "deterministic-fixed-width-table-parser"
                if fixed_width_variants and not vertical_variants
                else "deterministic-table-layout-parser"
            )
            return ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=extracted_data,
                model_used=parser_model,
            )

        # Router-first extraction: ask a cheap LLM to classify tables, then
        # parse them deterministically. Cuts ~10× tokens per typical paper
        # vs sending 60k chars of full text. Falls through to full-text Tier 3
        # below if the router finds no usable tables OR the router LLM fails.
        settings = get_settings()
        router_variants: List[dict] = []
        if settings.enable_table_router and paper.gene_symbol and scanner_text:
            router_outcome = self._try_table_router(paper, scanner_text)
            if router_outcome is not None:
                routed_data = router_outcome.extracted_data or {}
                router_variants = routed_data.get("variants", []) or []
                if len(router_variants) >= deterministic_min_variants:
                    return router_outcome
                logger.info(
                    "PMID %s - Router parsed only %d variants; continuing with "
                    "full-text extraction and using router output as merge hints",
                    paper.pmid,
                    len(router_variants),
                )

        truncated_text = self._truncate_text_for_prompt(
            full_text, gene_symbol=paper.gene_symbol
        )

        # Pre-extract variants from tables BEFORE LLM call (provides hints to improve extraction)
        # Use scanner_text (original full text) so we see tables even if condensed
        pre_extracted_variants = self._extract_variants_from_tables(
            scanner_text, paper.gene_symbol
        )
        table_hint_variants = (
            router_variants + deterministic_variants + pre_extracted_variants
        )
        table_hints = self._format_table_hints(table_hint_variants)
        if table_hint_variants:
            logger.info(
                f"PMID {paper.pmid} - Pre-extracted {len(table_hint_variants)} variants from tables as LLM hints"
            )
            print(f"Pre-extracted {len(table_hint_variants)} variant hints from tables")

        # Run comprehensive variant scanner on ORIGINAL full text (catches narrative mentions + tables)
        # Uses scanner_text to bypass condensation and see all content
        scanner_result = scan_document_for_variants(
            scanner_text,
            gene_symbol=paper.gene_symbol or "UNKNOWN",
            source=f"PMID_{paper.pmid}",
        )
        scanner_variant_count = len(scanner_result.variants)
        scanner_merge_enabled = (
            scanner_variant_count <= SCANNER_REGEX_MERGE_MAX_VARIANTS
        )
        effective_scanner_variant_count = scanner_variant_count
        if scanner_variant_count and not scanner_merge_enabled:
            logger.warning(
                "PMID %s - Skipping scanner hints/merge for %d candidates; "
                "candidate count exceeds safety cap %d",
                paper.pmid,
                scanner_variant_count,
                SCANNER_REGEX_MERGE_MAX_VARIANTS,
            )
            print(
                f"Variant scanner found {scanner_variant_count} potential variants; "
                "skipping scanner hints/merge due safety cap"
            )
            scanner_hints = ""
            effective_scanner_variant_count = 0
        else:
            scanner_hints = scanner_result.get_hints_for_prompt(
                max_hints=SCANNER_MAX_HINTS
            )

        if scanner_variant_count:
            logger.info(
                f"PMID {paper.pmid} - Variant scanner found {scanner_variant_count} candidates"
            )
            if scanner_merge_enabled:
                print(
                    f"Variant scanner found {scanner_variant_count} potential variants in text"
                )

        # Combine all hints (table + scanner)
        all_hints = table_hints + scanner_hints

        # Use compact mode for high-variant papers to avoid output truncation
        use_compact = estimated_variants >= HIGH_VARIANT_THRESHOLD
        if use_compact:
            logger.info(
                f"PMID {paper.pmid} - Using COMPACT extraction mode ({estimated_variants} estimated variants)"
            )
            print(
                f"High-variant paper detected ({estimated_variants}+ rows) - using compact extraction mode"
            )
            prompt = COMPACT_EXTRACTION_PROMPT.format(
                gene_symbol=paper.gene_symbol or "UNKNOWN",
                title=paper.title or "Unknown Title",
                full_text=truncated_text + all_hints,
                pmid=paper.pmid,
                estimated_variants=estimated_variants,
            )
        else:
            prompt = EXTRACTION_PROMPT.format(
                gene_symbol=paper.gene_symbol or "UNKNOWN",
                title=paper.title or "Unknown Title",
                full_text=truncated_text + all_hints,
                pmid=paper.pmid,
            )

        try:
            # Use the new method that tracks truncation status
            extracted_data, was_truncated, raw_text = self.call_llm_json_with_status(
                prompt
            )

            # Normalize stop codon notation
            extracted_data = self._normalize_stop_codon_notation(extracted_data)
            # Populate penetrance data from patient counts
            extracted_data = self._populate_penetrance_from_patient_count(
                extracted_data
            )

            # Check if we need continuation extraction
            variants = extracted_data.get("variants", [])
            expected_count = extracted_data.get("extraction_metadata", {}).get(
                "total_variants_found", len(variants)
            )

            if was_truncated and len(variants) < expected_count:
                extracted_data = self._attempt_continuation(
                    paper, model, extracted_data, full_text
                )
                # Normalize again after continuation
                extracted_data = self._normalize_stop_codon_notation(extracted_data)
                extracted_data = self._populate_penetrance_from_patient_count(
                    extracted_data
                )

            # Merge pre-extracted table variants to catch any the LLM missed.
            # Above the normal cap, keep a bounded deterministic overflow path
            # instead of dropping the entire table safety-net.
            if table_hint_variants:
                extracted_data = self._merge_table_variants_with_overflow_qc(
                    extracted_data, table_hint_variants, paper.pmid
                )

            # NEW: Merge scanner-found variants (catches narrative mentions the LLM missed)
            if scanner_result.variants and scanner_merge_enabled:
                extracted_data = merge_scanner_results(
                    extracted_data,
                    scanner_result,
                    paper.gene_symbol or "UNKNOWN",
                    min_confidence=SCANNER_MERGE_MIN_CONFIDENCE,
                )
            elif scanner_result.variants:
                metadata = extracted_data.setdefault("extraction_metadata", {})
                metadata["scanner_merge_skipped"] = {
                    "candidate_count": scanner_variant_count,
                    "safety_cap": SCANNER_REGEX_MERGE_MAX_VARIANTS,
                    "reason": "candidate_count_exceeds_safety_cap",
                }

            extracted_data = self._suppress_repeated_study_wide_counts(extracted_data)
            extracted_data = self._suppress_repeated_study_wide_context_fields(
                extracted_data
            )

            # Filter variants to only keep those matching the target gene
            if paper.gene_symbol:
                extracted_data = self._filter_by_gene(extracted_data, paper.gene_symbol)
                # Filter out extraction artifacts (p.XXX, invalid positions, etc.)
                extracted_data = self._filter_extraction_artifacts(
                    extracted_data, paper.gene_symbol
                )

            extracted_data = self._maybe_adjudicate_extraction(
                paper=paper,
                primary_model=model,
                extracted_data=extracted_data,
                source_text=scanner_text,
                estimated_variants=estimated_variants,
                scanner_variant_count=effective_scanner_variant_count,
            )
            extracted_data = self._suppress_repeated_study_wide_counts(extracted_data)
            extracted_data = self._suppress_repeated_study_wide_context_fields(
                extracted_data
            )

            num_variants = len(extracted_data.get("variants", []))
            logger.info(
                f"PMID {paper.pmid} - Extraction with {model} successful. Found {num_variants} variants."
            )
            return ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=extracted_data,
                model_used=self._result_model_label(model, extracted_data),
            )
        except json.JSONDecodeError as e:
            logger.error(f"PMID {paper.pmid} - JSON parsing error with {model}: {e}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"JSON error: {e}",
                model_used=model,
            )
        except Exception as e:
            logger.error(f"PMID {paper.pmid} - Extraction failed with {model}: {e}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"Extraction error: {e}",
                model_used=model,
            )

    def extract(self, paper: Paper) -> ExtractionResult:
        """
        Extract structured variant data from a paper using a tiered model approach.
        """
        if not self.models:
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error="No models configured for extraction",
            )

        best_successful_result: Optional[ExtractionResult] = None
        best_successful_count = -1
        prepared_full_text = self._prepare_full_text(paper)
        # Keep original full text for scanner (bypasses condensation)
        scanner_full_text = paper.full_text if paper.full_text else prepared_full_text

        # Circuit breaker: Check input quality BEFORE attempting any model
        # This prevents wasting LLM calls on garbage/unusable input
        is_usable, quality_reason = self._assess_input_quality(
            prepared_full_text, paper.gene_symbol
        )
        if not is_usable:
            logger.warning(
                f"PMID {paper.pmid} - Circuit breaker triggered: {quality_reason}"
            )
            print(f"⚡ PMID {paper.pmid}: SKIPPED (circuit breaker) - {quality_reason}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"SKIPPED: {quality_reason}",
                model_used=None,
            )

        table_row_hint = self._estimate_table_rows(scanner_full_text)

        # If we detect a large table, raise the bar so we try the next (stronger) model
        adaptive_threshold = self.tier_threshold
        if table_row_hint >= ADAPTIVE_TABLE_THRESHOLD:
            table_based_threshold = min(50, max(5, table_row_hint // 3))
            adaptive_threshold = max(self.tier_threshold, table_based_threshold)
            logger.info(
                f"PMID {paper.pmid} - Detected {table_row_hint} table-like rows; "
                f"using adaptive variant threshold {adaptive_threshold}"
            )

        for idx, model in enumerate(self.models):
            result = self._attempt_extraction(
                paper,
                model,
                prepared_full_text=prepared_full_text,
                estimated_variants=table_row_hint,
            )

            if not result.success:
                if idx + 1 < len(self.models):
                    logger.info(
                        f"PMID {paper.pmid} - Extraction with {model} failed ({result.error}). Trying next model."
                    )
                    continue
                # If we have a previous successful result, return it instead of the failure
                if best_successful_result is not None:
                    logger.info(
                        f"PMID {paper.pmid} - Later model failed, falling back to previous successful result "
                        f"from {best_successful_result.model_used}."
                    )
                    return best_successful_result
                return result

            extracted_data = result.extracted_data or {}
            variants = extracted_data.get("variants", [])
            num_variants = len(variants) if isinstance(variants, list) else 0
            if num_variants > best_successful_count:
                best_successful_result = result
                best_successful_count = num_variants
            threshold = adaptive_threshold
            if num_variants < threshold and idx + 1 < len(self.models):
                next_model = self.models[idx + 1]
                logger.info(
                    f"PMID {paper.pmid} - Found {num_variants} variants with {model} "
                    f"(threshold: {threshold}). Retrying with {next_model}."
                )
                continue

            if best_successful_result is not result:
                logger.info(
                    f"PMID {paper.pmid} - Keeping prior extraction from "
                    f"{best_successful_result.model_used} with {best_successful_count} variants; "
                    f"{model} found {num_variants}."
                )
                return best_successful_result
            return result

        # If we looped through all models and have a successful result, return it
        if best_successful_result is not None:
            return best_successful_result

        return ExtractionResult(
            pmid=paper.pmid, success=False, error="Extraction failed with all models"
        )

    def extract_batch(self, papers: List[Paper]) -> List[ExtractionResult]:
        """Extract data from multiple papers."""
        return [self.extract(paper) for paper in papers]


def extract_variants_from_paper(
    paper: Paper,
    models: Optional[List[str]] = None,
    fulltext_dir: Optional[str] = None,
) -> ExtractionResult:
    """
    Convenience function to extract variants from a single paper.

    Args:
        paper: Paper object with text to extract from
        models: Optional list of model identifiers to use
        fulltext_dir: Optional directory where DATA_ZONES.md files are stored

    Returns:
        ExtractionResult with extracted variant data
    """
    extractor = ExpertExtractor(models=models, fulltext_dir=fulltext_dir)
    return extractor.extract(paper)
