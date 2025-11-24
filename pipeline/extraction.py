"""
Expert Extractor module (Tier 3) for the Tiered Biomedical Extraction Pipeline.

The heavy lifter that processes full-text papers and extracts structured
genetic variant data using advanced LLM prompting.
"""

import logging
import json
from typing import Optional, Dict, Any
from models import Paper, ExtractionResult
from pipeline.utils.llm_utils import BaseLLMCaller
from config.settings import get_settings

logger = logging.getLogger(__name__)


class ExpertExtractor(BaseLLMCaller):
    """
    Tier 3: Expert-level extraction using advanced LLM (GPT-4, etc.).
    Handles full-text papers and markdown tables to extract structured variant data.
    """

    # Complex JSON extraction prompt
    EXTRACTION_PROMPT = """You are an expert medical geneticist and data extraction specialist. Your task is to extract ALL genetic variant information from the provided scientific paper, with special emphasis on penetrance data (affected vs unaffected carriers).

Paper Title: {title}

Full Text (including tables):
{full_text}

EXTRACTION INSTRUCTIONS:
Extract ALL genetic variants mentioned in this paper with the following structured information:

For each variant, provide:
1. Gene Symbol (e.g., "BRCA1", "TP53")
2. Variant Notation:
   - cDNA notation (e.g., "c.1234G>A")
   - Protein notation (e.g., "p.Arg412His")
   - Genomic position (if available)
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
- For tables: extract ALL rows with variant data AND individual person data
- If a variant is mentioned multiple times, consolidate the information but preserve individual-level detail
- Include exact quotes for key clinical descriptions
- Note the specific section/table where each variant was found
- Distinguish between individual case reports and cohort study statistics
- When paper states "10 people with variant X, 4 with disease" → extract: total_carriers=10, affected=4, unaffected=6
- When individual cases are described, create individual_records entries for each person

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
- If full text is not available and only abstract is provided, note this limitation
- Be thorough but accurate - don't invent data not present in the paper
- If a field is not available, use null or an empty string as appropriate
- Preserve exact nomenclature from the paper
"""

    def __init__(
        self,
        model: Optional[str] = None,
        temperature: float = 0.0,
        max_tokens: int = 8000  # Increased to accommodate penetrance data
    ):
        """
        Initialize the Expert Extractor.

        Args:
            model: LiteLLM model identifier (e.g., 'gpt-4o', 'claude-3-opus-20240229'). If None, uses config.
            temperature: Model temperature (0.0 for most deterministic).
            max_tokens: Maximum tokens for response.
        """
        settings = get_settings()
        model = model or settings.extractor_model
        super().__init__(model=model, temperature=temperature, max_tokens=max_tokens)

    def _prepare_full_text(self, paper: Paper) -> str:
        """
        Prepare full text for extraction, including abstract if full text not available.

        Args:
            paper: Paper object.

        Returns:
            Formatted text for extraction.
        """
        if paper.full_text:
            return paper.full_text
        elif paper.abstract:
            logger.warning(
                f"PMID {paper.pmid} - Full text not available, using abstract only"
            )
            return f"[ABSTRACT ONLY - FULL TEXT NOT AVAILABLE]\n\n{paper.abstract}"
        else:
            return "[NO TEXT AVAILABLE]"

    def extract(self, paper: Paper) -> ExtractionResult:
        """
        Extract structured variant data from a paper using expert LLM.

        Args:
            paper: Paper object with full text.

        Returns:
            ExtractionResult with extracted data or error.
        """
        logger.info(f"PMID {paper.pmid} - Starting expert extraction with {self.model}")

        # Prepare full text
        full_text = self._prepare_full_text(paper)

        if full_text == "[NO TEXT AVAILABLE]":
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error="No text available for extraction",
                model_used=self.model
            )

        # Construct prompt
        prompt = self.EXTRACTION_PROMPT.format(
            title=paper.title or "Unknown Title",
            full_text=full_text[:30000],  # Truncate extremely long texts
            pmid=paper.pmid
        )

        try:
            # Use shared LLM utility (includes retry logic)
            extracted_data = self.call_llm_json(prompt)

            logger.info(
                f"PMID {paper.pmid} - Extraction successful. "
                f"Found {extracted_data.get('extraction_metadata', {}).get('total_variants_found', 0)} variants."
            )

            return ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=extracted_data,
                model_used=self.model
            )

        except json.JSONDecodeError as e:
            logger.error(f"PMID {paper.pmid} - JSON parsing error: {e}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"JSON parsing error: {str(e)}",
                model_used=self.model
            )

        except Exception as e:
            logger.error(f"PMID {paper.pmid} - Extraction failed: {e}")
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error=f"Extraction error: {str(e)}",
                model_used=self.model
            )

    def extract_batch(self, papers: list[Paper]) -> list[ExtractionResult]:
        """
        Extract data from multiple papers.

        Args:
            papers: List of Paper objects.

        Returns:
            List of ExtractionResults.
        """
        results = []
        for paper in papers:
            result = self.extract(paper)
            results.append(result)

        return results


def extract_variants_from_paper(paper: Paper, model: Optional[str] = None) -> ExtractionResult:
    """
    Convenience function to extract variants from a single paper.

    Args:
        paper: Paper object with full text.
        model: LiteLLM model identifier. If None, uses config.

    Returns:
        ExtractionResult.
    """
    extractor = ExpertExtractor(model=model)
    return extractor.extract(paper)