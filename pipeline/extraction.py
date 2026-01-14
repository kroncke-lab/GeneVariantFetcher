"""
Expert Extractor module (Tier 3) for the Tiered Biomedical Extraction Pipeline.

The heavy lifter that processes full-text papers and extracts structured
genetic variant data using advanced LLM prompting.
"""

import json
import logging
import re
from pathlib import Path
from typing import List, Optional

from config.settings import get_settings
from pipeline.prompts import (
    COMPACT_EXTRACTION_PROMPT,
    CONTINUATION_PROMPT,
    EXTRACTION_PROMPT,
    HIGH_VARIANT_THRESHOLD,
)
from utils.llm_utils import BaseLLMCaller
from utils.models import ExtractionResult, Paper

logger = logging.getLogger(__name__)


def _find_data_zones_file(pmid: str, search_dirs: Optional[List[str]] = None) -> Optional[Path]:
    """
    Search for a DATA_ZONES.md file for the given PMID.

    Args:
        pmid: PubMed ID to search for
        search_dirs: Optional list of directories to search in. If provided with
                     a single directory, searches that directory directly first.

    Returns:
        Path to DATA_ZONES.md if found, None otherwise
    """
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

    # Fallback: search common output directory patterns
    fallback_dirs = ['.', 'pmc_fulltext', 'output']
    for search_dir in fallback_dirs:
        path = Path(search_dir)
        if path.exists() and path.is_dir():
            zones_file = path / filename
            if zones_file.exists():
                logger.debug(f"Found DATA_ZONES.md at {zones_file} (fallback)")
                return zones_file

    # Also search subdirectories matching test/output patterns (e.g., test_gene_pmid/)
    cwd = Path('.')
    for subdir in cwd.iterdir():
        if subdir.is_dir() and (subdir.name.startswith('test_') or subdir.name.startswith('output_') or pmid in subdir.name):
            zones_file = subdir / filename
            if zones_file.exists():
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

        self.models = models or settings.tier3_models
        self.temperature = temperature if temperature is not None else settings.tier3_temperature
        # Store the requested max tokens, but clamp per-model to avoid API errors
        self.requested_max_tokens = max_tokens if max_tokens is not None else settings.tier3_max_tokens
        self.max_tokens = self._clamp_max_tokens(self.models[0], self.requested_max_tokens)
        self.tier_threshold = tier_threshold
        self.fulltext_dir = fulltext_dir
        self.use_condensed = settings.scout_use_condensed

        super().__init__(model=self.models[0], temperature=self.temperature, max_tokens=self.max_tokens)
        logger.debug(f"ExpertExtractor initialized with models={self.models}, temp={self.temperature}, max_tokens={self.max_tokens}")

    def _clamp_max_tokens(self, model: str, requested: int) -> int:
        """
        Clamp max_tokens to model-safe limits to avoid provider errors.

        OpenAI 4o/4o-mini currently cap at 16384 completion tokens; use 15000
        as a safety margin.
        """
        limit = None
        m = model.lower() if model else ""
        if "gpt-4o" in m:
            limit = 15000
        elif "gpt-3.5" in m:
            limit = 4000

        if limit is None:
            return requested

        if requested > limit:
            logger.warning(
                f"Requested max_tokens={requested} exceeds safe limit {limit} for model {model}; clamping."
            )
        return min(requested, limit)

    # Minimum size for DATA_ZONES.md to be considered useful
    # If smaller than this, fall back to full text
    MIN_CONDENSED_SIZE = 500

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
                    condensed_text = zones_file.read_text(encoding='utf-8')

                    # Check if condensed text is useful (not just headers/metadata)
                    # Look for actual data content, not just "No high-value data zones identified"
                    has_no_zones = "No high-value data zones identified" in condensed_text
                    is_too_small = len(condensed_text) < self.MIN_CONDENSED_SIZE

                    if has_no_zones or is_too_small:
                        logger.warning(
                            f"PMID {paper.pmid} - DATA_ZONES.md too small or empty "
                            f"({len(condensed_text)} chars), falling back to full text"
                        )
                        print(f"⚠ DATA_ZONES.md insufficient ({len(condensed_text)} chars) - using full text instead")
                        # Fall through to use full text
                    elif condensed_text and len(condensed_text) > self.MIN_CONDENSED_SIZE:
                        lines = len(condensed_text.splitlines())
                        logger.info(f"PMID {paper.pmid} - Using condensed {paper.pmid}_DATA_ZONES.md ({len(condensed_text)} chars)")
                        print(f"Using {paper.pmid}_DATA_ZONES.md for extraction: {len(condensed_text):,} chars, {lines:,} lines")
                        return condensed_text
                except Exception as e:
                    logger.warning(f"PMID {paper.pmid} - Failed to read {paper.pmid}_DATA_ZONES.md: {e}")

        # Fall back to paper.full_text
        if paper.full_text:
            lines = len(paper.full_text.splitlines())
            print(f"Using {paper.pmid}_FULL_CONTEXT.md for extraction: {len(paper.full_text):,} chars, {lines:,} lines")
            return paper.full_text
        elif paper.abstract:
            logger.warning(f"PMID {paper.pmid} - Full text not available, using abstract only")
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
        context_window: int = 80,
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
        max_chars: int = 60000,
    ) -> str:
        """
        Keep three slices (head/mid/tail) so large supplemental tables aren’t dropped.
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
            if stripped.startswith("|") and stripped.count("|") >= 3:
                count += 1
        return count

    def _parse_markdown_table_variants(self, full_text: str, gene_symbol: Optional[str]) -> List[dict]:
        """
        Best-effort parser for simple markdown tables (fast path for very large tables).

        Returns a minimal variant list without calling the LLM. Intended for papers
        like PMID 19716085 where a single giant table lists hundreds of variants.
        """
        if not gene_symbol:
            return []

        lines = full_text.splitlines()
        table_started = False
        variants = []
        header_idx = {}

        for line in lines:
            if not table_started:
                # Detect header row
                if "Nucleotide" in line and "Variant" in line and "patient" in line:
                    parts = [p.strip() for p in line.split("|") if p.strip()]
                    for idx, name in enumerate(parts):
                        header_idx[name.lower()] = idx
                    table_started = True
                continue

            # Stop when table ends
            if not line.strip().startswith("|"):
                if variants:
                    break
                else:
                    table_started = False
                    continue

            # Skip separator rows
            if set(line.strip()) <= {"|", "-", " "}:
                continue

            cells = [c.strip() for c in line.split("|") if c.strip()]
            if len(cells) < 3:
                continue

            def get_col(key: str) -> Optional[str]:
                idx = header_idx.get(key)
                if idx is None or idx >= len(cells):
                    return None
                return cells[idx] or None

            cdna = get_col("nucleotide")
            protein = get_col("variant") or get_col("amino acid")  # safety
            patient_count_raw = get_col("no. of patients") or get_col("no. of patient") or get_col("patients")

            if not cdna and not protein:
                continue

            # Clean cdna/protein formatting
            cdna = cdna.replace(" ", "") if cdna else None
            protein = protein.replace(" ", "") if protein else None

            # Parse patient count
            patient_count = None
            if patient_count_raw:
                try:
                    patient_count = int(patient_count_raw)
                except ValueError:
                    patient_count = None
            # Fallback: use last cell if it looks numeric
            if patient_count is None and cells:
                tail = cells[-1]
                if tail.isdigit():
                    patient_count = int(tail)

            variant = {
                "gene_symbol": gene_symbol,
                "cdna_notation": f"c.{cdna}" if cdna and not cdna.startswith("c.") else cdna,
                "protein_notation": protein,
                "clinical_significance": "pathogenic",  # table is disease-associated
                "patients": {"count": patient_count, "phenotype": "LQT2"},
                "penetrance_data": {
                    "total_carriers_observed": patient_count,
                    "affected_count": patient_count,
                    "unaffected_count": 0 if patient_count is not None else None,
                },
                "individual_records": [],
                "functional_data": {"summary": "", "assays": []},
                "segregation_data": None,
                "population_frequency": None,
                "evidence_level": "medium",
                "source_location": "Table 2",
                "additional_notes": "Parsed via deterministic table parser",
                "key_quotes": [],
            }
            variants.append(variant)

        if variants:
            logger.info(f"Parsed {len(variants)} variants via deterministic markdown table parser")
        return variants

    def _get_extracted_variants_summary(self, variants: list) -> str:
        """Create a compact summary of extracted variants for continuation prompts."""
        summaries = []
        for v in variants:
            cdna = v.get('cdna_notation', '') or ''
            protein = v.get('protein_notation', '') or ''
            summaries.append(f"- {cdna} / {protein}")
        return "\n".join(summaries)

    def _merge_continuation_results(self, base_data: dict, continuation_data: dict) -> dict:
        """Merge continuation extraction results into base results."""
        # Add continuation variants
        base_variants = base_data.get('variants', [])
        continuation_variants = continuation_data.get('variants', [])

        # Deduplicate by cdna_notation + protein_notation
        existing_keys = set()
        for v in base_variants:
            key = (v.get('cdna_notation', ''), v.get('protein_notation', ''))
            existing_keys.add(key)

        new_variants = []
        for v in continuation_variants:
            key = (v.get('cdna_notation', ''), v.get('protein_notation', ''))
            if key not in existing_keys:
                new_variants.append(v)
                existing_keys.add(key)

        base_variants.extend(new_variants)
        base_data['variants'] = base_variants

        # Update metadata
        if 'extraction_metadata' in base_data:
            base_data['extraction_metadata']['total_variants_found'] = len(base_variants)
            continuation_count = continuation_data.get('extraction_metadata', {}).get('continuation_variants_found', len(new_variants))
            base_data['extraction_metadata']['notes'] = (
                base_data['extraction_metadata'].get('notes', '') +
                f" [Continuation added {continuation_count} variants]"
            ).strip()

        return base_data

    def _attempt_continuation(self, paper: Paper, model: str, base_data: dict, full_text: str) -> dict:
        """Attempt to extract remaining variants after truncation."""
        variants = base_data.get('variants', [])
        expected_count = base_data.get('extraction_metadata', {}).get('total_variants_found', len(variants))

        # Only continue if we're missing a significant number of variants
        if len(variants) >= expected_count or expected_count - len(variants) < 5:
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
            full_text=self._truncate_text_for_prompt(full_text, gene_symbol=paper.gene_symbol)
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
        for variant in extracted_data.get('variants', []):
            protein = variant.get('protein_notation', '')
            if not protein:
                continue

            original = protein

            # Convert 'Ter' (termination) at the end of protein notation to '*'
            # e.g., p.Arg412Ter -> p.Arg412*, p.R412Ter -> p.R412*
            protein = re.sub(r'([A-Za-z]{1,3}\d+)Ter$', r'\1*', protein, flags=re.IGNORECASE)

            # Convert 'Stop' at the end to '*'
            protein = re.sub(r'([A-Za-z]{1,3}\d+)Stop$', r'\1*', protein, flags=re.IGNORECASE)

            # Convert X at the end of protein notation to '*' (stop codon)
            # e.g., p.Arg412X -> p.Arg412*, p.R412X -> p.R412*
            protein = re.sub(r'([A-Za-z]{1,3}\d+)X$', r'\1*', protein)

            # Pattern for frameshift with X or Ter: fs + digits + X/Ter
            # e.g., p.Gly24fs+34X -> p.Gly24fs*34, p.Gly24fsTer58 -> p.Gly24fs*58
            protein = re.sub(r'(fs\+?)(\d*)X$', r'\1*\2', protein)
            protein = re.sub(r'(fs\+?)(\d*)Ter$', r'\1*\2', protein, flags=re.IGNORECASE)
            protein = re.sub(r'(fs)\*(\d+)X$', r'\1*\2', protein)

            # Handle cases where frameshift shows as fs*NUMBER or fsX NUMBER
            # Normalize fs*58 format (already correct but might have extra chars)
            protein = re.sub(r'(fs)\s*\*\s*(\d+)', r'\1*\2', protein)

            # Handle truncating mutations that end with unusual patterns
            # Sometimes 'stop gained' is written as the amino acid that replaces (incorrectly)
            # This is harder to detect without context, but we can flag suspicious patterns

            # If notation changed, update it
            if protein != original:
                variant['protein_notation'] = protein
                notes = variant.get('additional_notes', '') or ''
                if notes:
                    notes += ' '
                notes += f'[Stop codon notation normalized: {original} -> {protein}]'
                variant['additional_notes'] = notes

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
        for variant in extracted_data.get('variants', []):
            patients = variant.get('patients', {})
            patient_count = patients.get('count') if patients else None

            # Skip if no patient count
            if patient_count is None or patient_count == 0:
                continue

            # Get or create penetrance_data
            pdata = variant.get('penetrance_data', {})
            if pdata is None:
                pdata = {}
                variant['penetrance_data'] = pdata

            # Only populate if fields are missing/null
            total_carriers = pdata.get('total_carriers_observed')
            affected_count = pdata.get('affected_count')

            if total_carriers is None and affected_count is None:
                # Determine if this is a pathogenic/disease-associated variant
                significance = (variant.get('clinical_significance') or '').lower()
                phenotype = (patients.get('phenotype') or '').lower()
                source_location = (variant.get('source_location') or '').lower()

                # Check if from a disease-associated context
                is_disease_associated = any([
                    'pathogenic' in significance,
                    'lqt' in phenotype,  # Long QT syndrome
                    'brugada' in phenotype,
                    'arrhythmia' in phenotype,
                    'syndrome' in phenotype,
                    'disease' in phenotype,
                    'affected' in phenotype,
                    'mutation' in source_location,
                    'lqt' in source_location,
                ])

                # For disease-associated variants, patient count = affected count
                if is_disease_associated:
                    pdata['total_carriers_observed'] = patient_count
                    pdata['affected_count'] = patient_count
                    pdata['unaffected_count'] = 0

                    # Add note about the mapping
                    notes = variant.get('additional_notes', '') or ''
                    if notes:
                        notes += ' '
                    notes += f'[Penetrance data inferred from patient count: {patient_count} affected carriers]'
                    variant['additional_notes'] = notes
                else:
                    # For uncertain significance, just set total carriers
                    pdata['total_carriers_observed'] = patient_count

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
    ]

    def _assess_input_quality(self, text: str, gene_symbol: Optional[str]) -> tuple[bool, str]:
        """
        Assess whether the input text is of sufficient quality for extraction.

        Returns:
            (is_usable, reason) - whether text is usable and why/why not
        """
        if not text or len(text) < 200:
            return False, f"Text too short ({len(text) if text else 0} chars)"

        # Count failed extraction placeholders
        failed_count = sum(1 for pattern in self.FAILED_EXTRACTION_PATTERNS if pattern in text)

        # If most of the content is failure placeholders, skip
        if failed_count >= 3:
            return False, f"Multiple failed extraction placeholders ({failed_count} found)"

        # Check for actual variant-like content
        variant_patterns = [
            r'c\.\d+', r'p\.[A-Z]', r'[A-Z]\d+[A-Z]',  # Basic variant patterns
            r'mutation', r'variant', r'carrier',  # Clinical terms
        ]
        has_variant_content = any(re.search(p, text, re.IGNORECASE) for p in variant_patterns)

        if not has_variant_content:
            # Check if gene is at least mentioned
            if gene_symbol and gene_symbol.upper() in text.upper():
                return True, "Gene mentioned but no variant patterns found"
            return False, "No variant patterns or gene mentions in text"

        return True, "Text appears usable"

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
        self.model = model
        self.max_tokens = self._clamp_max_tokens(model, self.requested_max_tokens)

        full_text = prepared_full_text if prepared_full_text is not None else self._prepare_full_text(paper)
        if full_text == "[NO TEXT AVAILABLE]":
            return ExtractionResult(pmid=paper.pmid, success=False, error="No text available", model_used=model)

        # Assess input quality before sending to LLM
        is_usable, quality_reason = self._assess_input_quality(full_text, paper.gene_symbol)
        if not is_usable:
            logger.warning(f"PMID {paper.pmid} - Input quality check failed: {quality_reason}")
            print(f"⚠ PMID {paper.pmid}: Skipping LLM extraction - {quality_reason}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"Input quality insufficient: {quality_reason}",
                model_used=model
            )

        # Estimate variant count if not provided
        if estimated_variants is None:
            estimated_variants = self._estimate_table_rows(full_text)

        # Fast path: deterministic table parser for very large tables to avoid slow LLM calls
        if estimated_variants >= 100:
            parsed_variants = self._parse_markdown_table_variants(full_text, paper.gene_symbol)
            if len(parsed_variants) >= 50:
                extracted_data = {
                    "paper_metadata": {
                        "pmid": paper.pmid,
                        "title": paper.title or "Unknown Title",
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

        truncated_text = self._truncate_text_for_prompt(full_text, gene_symbol=paper.gene_symbol)

        # Use compact mode for high-variant papers to avoid output truncation
        use_compact = estimated_variants >= HIGH_VARIANT_THRESHOLD
        if use_compact:
            logger.info(f"PMID {paper.pmid} - Using COMPACT extraction mode ({estimated_variants} estimated variants)")
            print(f"High-variant paper detected ({estimated_variants}+ rows) - using compact extraction mode")
            prompt = COMPACT_EXTRACTION_PROMPT.format(
                gene_symbol=paper.gene_symbol or "UNKNOWN",
                title=paper.title or "Unknown Title",
                full_text=truncated_text,
                pmid=paper.pmid,
                estimated_variants=estimated_variants,
            )
        else:
            prompt = EXTRACTION_PROMPT.format(
                gene_symbol=paper.gene_symbol or "UNKNOWN",
                title=paper.title or "Unknown Title",
                full_text=truncated_text,
                pmid=paper.pmid
            )

        try:
            # Use the new method that tracks truncation status
            extracted_data, was_truncated, raw_text = self.call_llm_json_with_status(prompt)

            # Normalize stop codon notation
            extracted_data = self._normalize_stop_codon_notation(extracted_data)
            # Populate penetrance data from patient counts
            extracted_data = self._populate_penetrance_from_patient_count(extracted_data)

            # Check if we need continuation extraction
            variants = extracted_data.get('variants', [])
            expected_count = extracted_data.get('extraction_metadata', {}).get('total_variants_found', len(variants))

            if was_truncated and len(variants) < expected_count:
                extracted_data = self._attempt_continuation(paper, model, extracted_data, full_text)
                # Normalize again after continuation
                extracted_data = self._normalize_stop_codon_notation(extracted_data)
                extracted_data = self._populate_penetrance_from_patient_count(extracted_data)

            num_variants = len(extracted_data.get('variants', []))
            logger.info(f"PMID {paper.pmid} - Extraction with {model} successful. Found {num_variants} variants.")
            return ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=extracted_data,
                model_used=model
            )
        except json.JSONDecodeError as e:
            logger.error(f"PMID {paper.pmid} - JSON parsing error with {model}: {e}")
            return ExtractionResult(pmid=paper.pmid, success=False, error=f"JSON error: {e}", model_used=model)
        except Exception as e:
            logger.error(f"PMID {paper.pmid} - Extraction failed with {model}: {e}")
            return ExtractionResult(pmid=paper.pmid, success=False, error=f"Extraction error: {e}", model_used=model)

    def extract(self, paper: Paper) -> ExtractionResult:
        """
        Extract structured variant data from a paper using a tiered model approach.
        """
        if not self.models:
            return ExtractionResult(pmid=paper.pmid, success=False, error="No models configured for extraction")

        best_successful_result: Optional[ExtractionResult] = None
        prepared_full_text = self._prepare_full_text(paper)
        table_row_hint = self._estimate_table_rows(prepared_full_text)

        # If we detect a large table, raise the bar so we try the next (stronger) model
        adaptive_threshold = self.tier_threshold
        if table_row_hint >= 50:
            table_based_threshold = min(50, max(5, table_row_hint // 3))
            adaptive_threshold = max(self.tier_threshold, table_based_threshold)
            logger.info(
                f"PMID {paper.pmid} - Detected {table_row_hint} table-like rows; "
                f"using adaptive variant threshold {adaptive_threshold}"
            )

        for idx, model in enumerate(self.models):
            result = self._attempt_extraction(
                paper, model,
                prepared_full_text=prepared_full_text,
                estimated_variants=table_row_hint,
            )

            if not result.success:
                if idx + 1 < len(self.models):
                    logger.info(f"PMID {paper.pmid} - Extraction with {model} failed ({result.error}). Trying next model.")
                    continue
                # If we have a previous successful result, return it instead of the failure
                if best_successful_result is not None:
                    logger.info(
                        f"PMID {paper.pmid} - Later model failed, falling back to previous successful result "
                        f"from {best_successful_result.model_used}."
                    )
                    return best_successful_result
                return result

            num_variants = result.extracted_data.get('extraction_metadata', {}).get('total_variants_found', 0)
            threshold = adaptive_threshold
            if num_variants < threshold and idx + 1 < len(self.models):
                # Store this successful result as a fallback in case next model fails
                # Always update to the most recent successful result
                best_successful_result = result
                next_model = self.models[idx + 1]
                logger.info(
                    f"PMID {paper.pmid} - Found {num_variants} variants with {model} "
                    f"(threshold: {threshold}). Retrying with {next_model}."
                )
                continue

            return result

        # If we looped through all models and have a successful result, return it
        if best_successful_result is not None:
            return best_successful_result

        return ExtractionResult(pmid=paper.pmid, success=False, error="Extraction failed with all models")

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
