"""
Supplement-Specific Extraction for Genetic Variants.

Handles extraction from supplementary materials which have different characteristics
than main text - dense tables, appendix sections, figure legends, and methods.
"""

import logging
import re
from typing import List, Optional, Dict, Any
from pathlib import Path

from utils.models import Paper
from utils.llm_utils import BaseLLMCaller
from pipeline.prompts import (
    SUPPLEMENT_TABLE_PROMPT,
    SUPPLEMENT_APPENDIX_PROMPT,
)

logger = logging.getLogger(__name__)


class SupplementExtractor:
    """Specialized extraction for supplementary materials"""

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        # Higher token budget for supplements (tables are dense)
        self.max_tokens = 8000
        self.llm_caller = BaseLLMCaller(
            model=config.get("model", "gpt-4"),
            temperature=0.2,
            max_tokens=8000,
        )

    def identify_supplement_type(self, content: str) -> str:
        """Classify content type from supplement files"""
        content_lower = content.lower()

        # Look for markdown table patterns
        table_indicators = [
            "|",
            "\t",  # Tab-separated
            ",",  # CSV
        ]
        table_lines = sum(
            1
            for line in content.splitlines()
            if any(indicator in line for indicator in table_indicators)
            and len(line.strip()) > 10
        )

        # Table content with headers
        if content.count("|") > 10 and any(
            keyword in content_lower
            for keyword in ["variant", "mutation", "patient", "phenotype"]
        ):
            return "table"

        # Appendix patterns
        if any(
            keyword in content_lower
            for keyword in ["appendix", "supplementary", "supplemental"]
        ):
            return "appendix"

        # Methods/protocol descriptions
        methods_keywords = [
            "methods",
            "protocol",
            "procedure",
            "experimental",
            "assay",
            "functional study",
            "in vitro",
        ]
        if any(keyword in content_lower for keyword in methods_keywords):
            return "methods"

        # Figure legends/image descriptions
        figure_keywords = ["figure", "table", "legend", "panel", "supplementary figure"]
        if any(keyword in content_lower for keyword in figure_keywords):
            return "figure_legend"

        return "other"

    def extract_from_table(self, table_text: str, gene_symbol: str) -> List[dict]:
        """Extract variants from supplement tables (typically dense, structured)"""
        supplement_type = self.identify_supplement_type(table_text)

        if supplement_type != "table":
            logger.info(f"Content type {supplement_type} requested table extraction")
            return []

        # Attempt deterministic extraction for known table formats
        variants = self._parse_standard_supplement_table(table_text, gene_symbol)
        if variants:
            logger.info(f"Deterministic extraction found {len(variants)} variants")
            return variants

        # Use LLM for complex tables
        prompt = SUPPLEMENT_TABLE_PROMPT.format(
            gene_symbol=gene_symbol, table_content=table_text
        )

        try:
            result = self.llm_caller.call_llm_json(prompt)
            variants = result.get("variants", [])
            logger.info(f"LLM extraction found {len(variants)} variants from table")
            return variants
        except Exception as e:
            logger.error(f"Failed to extract from supplement table: {e}")
            return []

    def extract_from_appendix(self, appendix_text: str, gene_symbol: str) -> List[dict]:
        """Extract variants from appendix sections (typically narrative)"""
        supplement_type = self.identify_supplement_type(appendix_text)

        if supplement_type != "appendix":
            logger.info(f"Content type {supplement_type} requested appendix extraction")
            return []

        # Use appendix-specific prompt for narrative content
        prompt = SUPPLEMENT_APPENDIX_PROMPT.format(
            gene_symbol=gene_symbol, appendix_content=appendix_text
        )

        try:
            result = self.llm_caller.call_llm_json(prompt)
            variants = result.get("variants", [])
            logger.info(f"Found {len(variants)} variants from appendix")
            return variants
        except Exception as e:
            logger.error(f"Failed to extract from appendix: {e}")
            return []

    def extract_from_methods(self, methods_text: str, gene_symbol: str) -> List[dict]:
        """Extract variants from methods sections (protocol details)"""
        # Methods sections typically contain variant creation/engineering info
        # but rarely contain patient data. Extract functional data instead.
        variants = []

        # Parse functional assay data
        functional_patterns = [
            r"(p\.\w+\d+\w+|c\.\w+)",  # Protein/cDNA notation
            r"(?:functional|assay|current|trafficking)",  # Functional keywords
        ]

        lines = methods_text.splitlines()
        for line in lines:
            variants.extend(self._extract_functional_variants(line, gene_symbol))

        return variants

    def extract_from_figure_legend(
        self, legend_text: str, gene_symbol: str
    ) -> List[dict]:
        """Extract variants from figure legends (often contains supplementary data)"""
        # Check for supplementary table references in figure legends
        return []

    def _parse_standard_supplement_table(
        self, table_text: str, gene_symbol: str
    ) -> List[dict]:
        """Parse standard supplement table formats"""
        variants = []
        lines = table_text.splitlines()

        # Skip empty lines at start
        lines = [line for line in lines if line.strip()]
        if not lines:
            return variants

        # Detect header row
        header_row = None
        for i, line in enumerate(lines):
            line_lower = line.lower()
            if any(
                keyword in line_lower
                for keyword in ["variant", "mutation", "cDNA", "protein", "patient"]
            ):
                header_row = i
                break

        if header_row is None:
            return variants

        # Parse header and get column mapping
        header_line = lines[header_row]
        header_cols = self._parse_columns(header_line)

        if not header_cols:
            return variants

        # Map column names to standard fields
        column_mapping = self._map_columns_to_fields(header_cols)

        # Parse data rows
        for line in lines[header_row + 1 :]:
            line = line.strip()
            if not line or line.startswith("---") or line == "":
                continue

            variant = self._parse_table_row(line, column_mapping, gene_symbol)
            if variant and variant.get("gene_symbol") == gene_symbol:
                variants.append(variant)

        return variants

    def _parse_columns(self, header_line: str) -> List[str]:
        """Parse columns from headers (handles various delimiters)"""
        # Handle markdown table headers
        if "|" in header_line:
            parts = [part.strip() for part in header_line.split("|")]
            return [part for part in parts if part and not re.match(r"^-+$", part)]

        # Handle tab-separated
        if "\t" in header_line:
            return [part.strip() for part in header_line.split("\t")]

        # Handle comma-separated CSV
        if "," in header_line:
            parts = []
            for part in header_line.split(","):
                part = part.strip().strip("\"'")
                parts.append(part)
            return parts

        # Clean up column names
        return [col.strip() for col in header_line.split() if col.strip()]

    def _map_columns_to_fields(
        self, column_names: List[str]
    ) -> Dict[str, Optional[int]]:
        """Map detected column names to standard fields"""
        mapping = {
            "cdna_notation": None,
            "protein_notation": None,
            "patient_count": None,
            "phenotype": None,
            "clinical_significance": None,
            "functional_data": None,
            "segregation": None,
        }

        for idx, col_name in enumerate(column_names):
            col_lower = col_name.lower()

            if any(
                kw in col_lower for kw in ["nucleotide", "cDNA", "variant", "mutation"]
            ):
                mapping["cdna_notation"] = idx
            elif any(kw in col_lower for kw in ["protein", "amino acid", "aa"]):
                mapping["protein_notation"] = idx
            elif any(kw in col_lower for kw in ["patient", "case", "n", "count"]):
                mapping["patient_count"] = idx
            elif any(
                kw in col_lower for kw in ["phenotype", "clinical", "presentation"]
            ):
                mapping["phenotype"] = idx
            elif any(
                kw in col_lower for kw in ["significance", "pathogenic", "benign"]
            ):
                mapping["clinical_significance"] = idx
            elif any(kw in col_lower for kw in ["functional", "assay", "current"]):
                mapping["functional_data"] = idx
            elif any(
                kw in col_lower for kw in ["segregation", "family", "inheritance"]
            ):
                mapping["segregation"] = idx

        return mapping

    def _parse_table_row(
        self, row_line: str, mapping: Dict[str, Optional[int]], gene_symbol: str
    ) -> Optional[Dict]:
        """Parse a single table row into variant structure"""
        # Parse row based on delimiter
        if "|" in row_line:
            cells = [cell.strip() for cell in row_line.split("|")]
        elif "\t" in row_line:
            cells = [cell.strip() for cell in row_line.split("\t")]
        else:
            cells = [cell.strip() for cell in self._split_csv_row(row_line)]

        if len(cells) <= max(pos for pos in mapping.values() if pos is not None):
            return None

        # Extract variant data
        cdna = (
            None
            if mapping["cdna_notation"] is None
            else cells[mapping["cdna_notation"]]
        )
        protein = (
            None
            if mapping["protein_notation"] is None
            else cells[mapping["protein_notation"]]
        )
        patient_count = (
            None
            if mapping["patient_count"] is None
            else cells[mapping["patient_count"]]
        )
        phenotype = (
            None if mapping["phenotype"] is None else cells[mapping["phenotype"]]
        )

        # Clean and normalize data
        variant = {
            "gene_symbol": gene_symbol,
            "cdna_notation": self._normalize_notation(cdna, "cdna"),
            "protein_notation": self._normalize_notation(protein, "protein"),
            "clinical_significance": "unknown",  # Default
            "patients": {},
            "penetrance_data": {},
            "individual_records": [],
            "functional_data": {"summary": "", "assays": []},
            "evidence_level": "medium",
            "source_location": "Supplementary Table",
            "additional_notes": "",
            "key_quotes": [],
        }

        # Parse patient count
        try:
            if patient_count:
                count = int(patient_count.strip())
                variant["patients"]["count"] = count
                variant["penetrance_data"]["total_carriers_observed"] = count
                variant["penetrance_data"]["affected_count"] = (
                    count  # Disease-associated
                )
        except (ValueError, TypeError):
            pass

        if phenotype:
            variant["patients"]["phenotype"] = phenotype

        return variant

    def _normalize_notation(
        self, notation: Optional[str], notation_type: str
    ) -> Optional[str]:
        """Normalize variant notation"""
        if not notation or notation in ["—", "-", "NA", "n/a"]:
            return None

        notation = notation.strip()

        if notation_type == "cdna" and not notation.startswith("c."):
            return f"c.{notation}"

        if notation_type == "protein" and notation and not notation.startswith("p."):
            # Handle protein shorthand like "R123W"
            match = re.match(r"([A-Za-z]\d+[A-Za-z*]+)", notation)
            if match:
                return f"p.{match.group(1)}"
            return f"p.{notation}"

        return notation

    def _split_csv_row(self, row: str) -> List[str]:
        """Simple CSV parsing"""
        # Handle quoted values
        parts = []
        current = ""
        in_quotes = False

        for char in row + ",":
            if char == '"':
                in_quotes = not in_quotes
            elif char == "," and not in_quotes:
                parts.append(current.strip())
                current = ""
            else:
                current += char

        return [part.strip() for part in parts if part.strip()]

    def _extract_functional_variants(
        self, text_line: str, gene_symbol: str
    ) -> List[Dict]:
        """Extract functional assay variants from methods text"""
        variants = []

        # Find protein notation patterns
        protein_pattern = re.compile(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}")
        cdna_pattern = re.compile(r"c\.\d+\w+")

        for match in protein_pattern.finditer(text_line):
            variants.append(
                {
                    "gene_symbol": gene_symbol,
                    "protein_notation": match.group(0),
                    "source_location": "Methods/Protocol",
                }
            )

        for match in cdna_pattern.finditer(text_line):
            variants.append(
                {
                    "gene_symbol": gene_symbol,
                    "cdna_notation": match.group(0),
                    "source_location": "Methods/Protocol",
                }
            )

        return variants

    def merge_main_and_supplement(
        self, main_variants: List[Dict], supplement_variants: List[Dict]
    ) -> List[Dict]:
        """Merge main text and supplement variants with deduplication"""
        from utils.variant_normalizer import normalize_variant

        # Create normalized keys for deduplication
        all_variants = []
        seen_keys = set()

        def get_variant_key(v):
            """Generate unique key for deduplication"""
            cdna = normalize_variant(v.get("cdna_notation", "") or "")
            protein = normalize_variant(v.get("protein_notation", "") or "")
            gene_symbol = v.get("gene_symbol", "").upper()
            return f"{gene_symbol}_{cdna}_{protein}"

        # Add main variants first (higher priority)
        for variant in main_variants:
            key = get_variant_key(variant)
            if key not in seen_keys:
                variant["source_priority"] = "main_text"
                all_variants.append(variant)
                seen_keys.add(key)

        # Add supplement variants, enriching main variants when possible
        for variant in supplement_variants:
            key = get_variant_key(variant)

            if key in seen_keys:
                # Enrich existing variant with supplement data
                existing = next(v for v in all_variants if get_variant_key(v) == key)
                self._enrich_variant_with_supplement(existing, variant)
            else:
                # Add new supplement variant
                variant["source_priority"] = "supplement"
                all_variants.append(variant)
                seen_keys.add(key)

        return all_variants

    def _enrich_variant_with_supplement(self, existing: Dict, supplement: Dict):
        """Enrich main variant with supplement data"""
        # Add supplement notes
        notes = existing.get("additional_notes", "") or ""
        if notes:
            notes += " | "
        notes += (
            f"Supplement: {supplement.get('source_location', 'supplementary data')}"
        )
        existing["additional_notes"] = notes

        # Merge patient counts if supplement has higher resolution
        existing_patients = existing.get("patients", {})
        supplement_patients = supplement.get("patients", {})

        if supplement_patients.get("count") and not existing_patients.get("count"):
            existing["patients"] = supplement_patients
