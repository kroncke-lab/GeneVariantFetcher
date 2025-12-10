"""
Expert Extractor module (Tier 3) for the Tiered Biomedical Extraction Pipeline.

The heavy lifter that processes full-text papers and extracts structured
genetic variant data using advanced LLM prompting.
"""

import logging
import json
import re
from pathlib import Path
from typing import Optional, List
from utils.models import Paper, ExtractionResult
from utils.llm_utils import BaseLLMCaller
from config.settings import get_settings

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
    """

    CONTINUATION_PROMPT = """You previously extracted variants from this paper but the response was truncated.
You extracted {extracted_count} variants so far. The paper contains approximately {expected_count} variants total.

Previously extracted variants (DO NOT re-extract these):
{extracted_variants_list}

Please continue extracting the REMAINING variants starting AFTER the last one listed above.
Return ONLY the variants you haven't extracted yet in the same JSON format.

TARGET GENE: {gene_symbol}
Paper Title: {title}

Full Text:
{full_text}

Return a JSON object with this structure:
{{
    "continuation": true,
    "variants": [
        ... remaining variants only ...
    ],
    "extraction_metadata": {{
        "continuation_variants_found": integer,
        "notes": "any notes about this continuation"
    }}
}}
"""

    EXTRACTION_PROMPT = """You are an expert medical geneticist and data extraction specialist. Your task is to extract genetic variant information from the provided scientific paper, with special emphasis on penetrance data (affected vs unaffected carriers).

RESPONSE SIZE GUIDANCE:
- You MUST extract ALL variants found in the paper - do not stop early or truncate the list
- For papers with many variants (large tables), use a simplified format per variant:
  * Include: gene_symbol, cdna_notation, protein_notation, clinical_significance, patients.count, source_location
  * Use empty arrays for: individual_records, key_quotes, functional_data.assays, age_dependent_penetrance
  * Use null for: genomic_position, penetrance_data counts (unless explicitly stated), segregation_data, population_frequency
  * Use brief strings for: patients.phenotype, additional_notes
- If there are >50 variants, limit individual_records to at most 1 per variant (or empty list)
- Prioritize completeness over detail: extract ALL variants with minimal fields rather than few variants with full detail

TARGET GENE: {gene_symbol}

CRITICAL: Only extract variants in the gene "{gene_symbol}". Ignore all variants in other genes.

Paper Title: {title}

Full Text (including tables):
{full_text}

EXTRACTION INSTRUCTIONS:
Extract ALL variants in the {gene_symbol} gene mentioned in this paper with the following structured information:

For each variant, provide:
1. Gene Symbol (e.g., "BRCA1", "TP53")
2. Variant Notation:
   - cDNA notation (e.g., "c.1234G>A") - convert from nucleotide position if needed (e.g., "47 A>C" → "c.47A>C")
   - Protein notation (e.g., "p.Arg412His") - use full HGVS notation:
     * Single-letter amino acids should be converted to three-letter (e.g., "D16A" → "p.Asp16Ala")
     * Nonsense/stop mutations: "X" or "*" in protein position means stop codon (e.g., "R412X" → "p.Arg412*")
     * Frameshift: include fs and stop position (e.g., "G24fs+34X" → "p.Gly24fs*58" or "p.Gly24fs")
   - Genomic position (if available)

   IMPORTANT - Novel mutation markers:
   - An asterisk (*) AFTER a mutation name (e.g., "D16A*", "R20G*") typically indicates "novel mutation"
     (first reported in this study). This is NOT a stop codon - note it in additional_notes as "novel mutation".
   - An asterisk (*) IN the protein position (e.g., "R412*", "p.Arg412*") indicates a stop codon/nonsense mutation.

3. Clinical Significance: pathogenic, likely pathogenic, benign, likely benign, VUS, etc.
4. Patient Information:
   - Number of patients/cases
   - Demographics (age, sex, ethnicity if mentioned)
   - Phenotype/disease presentation
5. Functional Data:
   - In vitro studies
   - Functional assays
   - Protein effects
6. Segregation Data: Does it segregate with disease in families?
7. Population Frequency: Allele frequency in databases (gnomAD, ExAC, etc.)
8. Evidence Level: How strong is the evidence for pathogenicity?
9. Additional Notes: Any other relevant clinical or functional information

PENETRANCE DATA EXTRACTION (CRITICAL):
For calibrating disease prediction models, extract detailed penetrance information:

A. Individual-Level Records:
   - Extract EVERY individual person mentioned who carries the variant
   - Look for: proband, case, patient, subject, individual, family member, sibling, parent, or designations like II-1, III-2, P1, Case 2, etc.
   - For each individual, capture:
     * Age at evaluation/assessment (if mentioned)
     * Age at disease onset (if mentioned)
     * Age at diagnosis (if mentioned)
     * Sex
     * Affected status: "affected", "unaffected", or "uncertain"
     * Phenotype details specific to that individual
     * Exact sentence where this information appears
   - Use age-appropriate penetrance logic: a young person may be "unaffected" but not past the risk window

B. Cohort/Aggregate Penetrance Data:
   - Extract study-level statistics when provided (e.g., "10 carriers, 4 affected, 6 unaffected")
   - Total carriers observed
   - Affected count
   - Unaffected count
   - Uncertain/unclear cases
   - Age-dependent penetrance data if stratified by age groups
   - Penetrance percentages if explicitly stated

C. Age-Dependent Penetrance:
   - Capture penetrance stratified by age ranges if mentioned
   - Note age at which penetrance was assessed
   - For age-dependent diseases, track penetrance by age group

CRITICAL REQUIREMENTS:
- Extract data from BOTH the main text AND all tables
- Pay special attention to supplementary table references
- For tables: int ALL rows with variant data AND individual person data
- If a variant is mentioned multiple times, consolidate the information but preserve individual-level detail
- Include exact quotes for key clinical descriptions
- Note the specific section/table where each variant was found
- Distinguish between individual case reports and cohort study statistics
- When paper states "10 people with variant X, 4 with disease" → extract: total_carriers=10, affected=4, unaffected=6
- When individual cases are described, create individual_records entries for each person
- IMPORTANT: For disease-associated mutation tables (e.g., "LQT2-associated mutations", "pathogenic variants"):
  * Patient counts in these tables represent AFFECTED carriers
  * Extract patient count to BOTH patients.count AND penetrance_data.total_carriers_observed/affected_count
  * Example: Table shows "No. of patients: 1" for a pathogenic variant → patients.count=1, total_carriers_observed=1, affected_count=1

OUTPUT FORMAT:
Return a JSON object with this structure:
{{
    "paper_metadata": {{
        "pmid": "{pmid}",
        "title": "{title}",
        "extraction_summary": "Brief summary of what was extracted"
    }},
    "variants": [
        {{
            "gene_symbol": "string",
            "cdna_notation": "string or null",
            "protein_notation": "string or null",
            "genomic_position": "string or null",
            "clinical_significance": "string",
            "patients": {{
                "count": integer,
                "demographics": "string",
                "phenotype": "string"
            }},
            "penetrance_data": {{
                "total_carriers_observed": "integer or null (total people with variant)",
                "affected_count": "integer or null (number with disease)",
                "unaffected_count": "integer or null (number without disease, past risk window)",
                "uncertain_count": "integer or null (unclear status or too young)",
                "penetrance_percentage": "float or null (calculated: affected/total_carriers * 100)",
                "age_dependent_penetrance": [
                    {{
                        "age_range": "string (e.g., '40-50 years')",
                        "penetrance_percentage": "float",
                        "carriers_in_range": "integer",
                        "affected_in_range": "integer"
                    }}
                ]
            }},
            "individual_records": [
                {{
                    "individual_id": "string (e.g., 'II-1', 'P1', 'Case_2', or generate unique ID)",
                    "age_at_evaluation": "integer or null",
                    "age_at_onset": "integer or null (age when symptoms started)",
                    "age_at_diagnosis": "integer or null (age when diagnosed)",
                    "sex": "string or null (male/female/other)",
                    "affected_status": "string (affected/unaffected/uncertain)",
                    "phenotype_details": "string (disease manifestations for this person)",
                    "evidence_sentence": "string (exact sentence from paper)"
                }}
            ],
            "functional_data": {{
                "summary": "string",
                "assays": ["list of assays performed"]
            }},
            "segregation_data": "string or null",
            "population_frequency": "string or null",
            "evidence_level": "string",
            "source_location": "e.g., 'Table 2, Row 3' or 'Results, paragraph 4'",
            "additional_notes": "string",
            "key_quotes": ["relevant quotes from paper"]
        }}
    ],
    "tables_processed": [
        {{
            "table_name": "string (e.g., 'Table 1', 'Supplementary Table 3')",
            "table_caption": "string",
            "variants_extracted": integer
        }}
    ],
    "extraction_metadata": {{
        "total_variants_found": integer,
        "extraction_confidence": "high/medium/low",
        "challenges": ["any issues during extraction"],
        "notes": "any additional notes about the extraction process"
    }}
}}

IMPORTANT NOTES:
- ONLY include variants in the {gene_symbol} gene. Do NOT include variants from other genes even if they are mentioned in the paper.
- If the paper mentions the target gene {gene_symbol} but does not report any variants, return an empty variants list.
- If full text is not available and only abstract is provided, note this limitation
- Be thorough but accurate - don't invent data not present in the paper
- If a field is not available, use null or an empty string as appropriate
- Preserve exact nomenclature from the paper
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
        self.max_tokens = max_tokens if max_tokens is not None else settings.tier3_max_tokens
        self.tier_threshold = tier_threshold
        self.fulltext_dir = fulltext_dir
        self.use_condensed = settings.scout_use_condensed

        super().__init__(model=self.models[0], temperature=self.temperature, max_tokens=self.max_tokens)
        logger.debug(f"ExpertExtractor initialized with models={self.models}, temp={self.temperature}, max_tokens={self.max_tokens}")

    def _prepare_full_text(self, paper: Paper) -> str:
        """
        Prepare full text for extraction.

        If scout_use_condensed is enabled and a DATA_ZONES.md file exists,
        prefer the condensed version for more focused extraction.
        """
        # Try to use condensed DATA_ZONES.md if enabled
        if self.use_condensed and paper.pmid:
            search_dirs = [self.fulltext_dir] if self.fulltext_dir else None
            zones_file = _find_data_zones_file(paper.pmid, search_dirs)

            if zones_file:
                try:
                    condensed_text = zones_file.read_text(encoding='utf-8')
                    if condensed_text and len(condensed_text) > 100:
                        lines = len(condensed_text.splitlines())
                        logger.info(f"PMID {paper.pmid} - Using condensed DATA_ZONES.md ({len(condensed_text)} chars)")
                        print(f"Using DATA_ZONES.md for extraction: {len(condensed_text):,} chars, {lines:,} lines")
                        return condensed_text
                except Exception as e:
                    logger.warning(f"PMID {paper.pmid} - Failed to read DATA_ZONES.md: {e}")

        # Fall back to paper.full_text
        if paper.full_text:
            lines = len(paper.full_text.splitlines())
            print(f"Using FULL_CONTEXT.md for extraction: {len(paper.full_text):,} chars, {lines:,} lines")
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

        prompt = self.CONTINUATION_PROMPT.format(
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

    def _attempt_extraction(self, paper: Paper, model: str, prepared_full_text: Optional[str] = None) -> ExtractionResult:
        """Attempt extraction with a single model, with continuation for truncated responses."""
        logger.info(f"PMID {paper.pmid} - Starting expert extraction with {model}")
        self.model = model

        full_text = prepared_full_text if prepared_full_text is not None else self._prepare_full_text(paper)
        if full_text == "[NO TEXT AVAILABLE]":
            return ExtractionResult(pmid=paper.pmid, success=False, error="No text available", model_used=model)

        truncated_text = self._truncate_text_for_prompt(full_text, gene_symbol=paper.gene_symbol)
        prompt = self.EXTRACTION_PROMPT.format(
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
            result = self._attempt_extraction(paper, model, prepared_full_text=prepared_full_text)

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
