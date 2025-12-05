"""
Expert Extractor module (Tier 3) for the Tiered Biomedical Extraction Pipeline.

The heavy lifter that processes full-text papers and extracts structured
genetic variant data using advanced LLM prompting.
"""

import logging
import json
from typing import Optional, List
from utils.models import Paper, ExtractionResult
from utils.llm_utils import BaseLLMCaller
from config.settings import get_settings

logger = logging.getLogger(__name__)


class ExpertExtractor(BaseLLMCaller):
    """
    Tier 3: Expert-level extraction using advanced LLM (GPT-4, etc.).
    Handles full-text papers and markdown tables to extract structured variant data.
    """

    EXTRACTION_PROMPT = """You are an expert medical geneticist and data extraction specialist. Your task is to extract genetic variant information from the provided scientific paper, with special emphasis on penetrance data (affected vs unaffected carriers).

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
- For tables: int ALL rows with variant data AND individual person data
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
    ):
        """
        Initialize the Expert Extractor.

        Args:
            models: List of LiteLLM model identifiers. If None, uses config (TIER3_MODELS).
            temperature: Model temperature. If None, uses config.
            max_tokens: Maximum tokens for response. If None, uses config.
            tier_threshold: If the first model finds fewer variants than this, the next model is tried.
        """
        settings = get_settings()

        self.models = models or settings.tier3_models
        self.temperature = temperature if temperature is not None else settings.tier3_temperature
        self.max_tokens = max_tokens if max_tokens is not None else settings.tier3_max_tokens
        self.tier_threshold = tier_threshold

        super().__init__(model=self.models[0], temperature=self.temperature, max_tokens=self.max_tokens)
        logger.debug(f"ExpertExtractor initialized with models={self.models}, temp={self.temperature}, max_tokens={self.max_tokens}")

    def _prepare_full_text(self, paper: Paper) -> str:
        """Prepare full text for extraction."""
        if paper.full_text:
            return paper.full_text
        elif paper.abstract:
            logger.warning(f"PMID {paper.pmid} - Full text not available, using abstract only")
            return f"[ABSTRACT ONLY - FULL TEXT NOT AVAILABLE]\n\n{paper.abstract}"
        else:
            return "[NO TEXT AVAILABLE]"

    def _attempt_extraction(self, paper: Paper, model: str) -> ExtractionResult:
        """Attempt extraction with a single model."""
        logger.info(f"PMID {paper.pmid} - Starting expert extraction with {model}")
        self.model = model

        full_text = self._prepare_full_text(paper)
        if full_text == "[NO TEXT AVAILABLE]":
            return ExtractionResult(pmid=paper.pmid, success=False, error="No text available", model_used=model)

        prompt = self.EXTRACTION_PROMPT.format(
            gene_symbol=paper.gene_symbol or "UNKNOWN",
            title=paper.title or "Unknown Title",
            full_text=full_text[:30000],
            pmid=paper.pmid
        )

        try:
            extracted_data = self.call_llm_json(prompt)
            num_variants = extracted_data.get('extraction_metadata', {}).get('total_variants_found', 0)
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

        first_model = self.models[0]
        result = self._attempt_extraction(paper, first_model)

        if result.success:
            num_variants = result.extracted_data.get('extraction_metadata', {}).get('total_variants_found', 0)
            if num_variants < self.tier_threshold and len(self.models) > 1:
                next_model = self.models[1]
                logger.info(f"PMID {paper.pmid} - Found {num_variants} variants with {first_model} (threshold: {self.tier_threshold}). Retrying with {next_model}.")
                return self._attempt_extraction(paper, next_model)

        return result

    def extract_batch(self, papers: List[Paper]) -> List[ExtractionResult]:
        """Extract data from multiple papers."""
        return [self.extract(paper) for paper in papers]


def extract_variants_from_paper(paper: Paper, models: Optional[List[str]] = None) -> ExtractionResult:
    """
    Convenience function to extract variants from a single paper.
    """
    extractor = ExpertExtractor(models=models)
    return extractor.extract(paper)
