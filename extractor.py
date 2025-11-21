"""
Expert Extractor module (Tier 3) for the Tiered Biomedical Extraction Pipeline.

The heavy lifter that processes full-text papers and extracts structured
genetic variant data using advanced LLM prompting.
"""

import logging
import json
from typing import Optional, Dict, Any
from litellm import completion
from tenacity import retry, stop_after_attempt, wait_exponential
from models import Paper, ExtractionResult

logger = logging.getLogger(__name__)


class ExpertExtractor:
    """
    Tier 3: Expert-level extraction using advanced LLM (GPT-4, etc.).
    Handles full-text papers and markdown tables to extract structured variant data.
    """

    # Complex JSON extraction prompt
    EXTRACTION_PROMPT = """You are an expert medical geneticist and data extraction specialist. Your task is to extract ALL genetic variant information from the provided scientific paper.

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

CRITICAL REQUIREMENTS:
- Extract data from BOTH the main text AND all tables
- Pay special attention to supplementary table references
- For tables: extract ALL rows with variant data
- If a variant is mentioned multiple times, consolidate the information
- Include exact quotes for key clinical descriptions
- Note the specific section/table where each variant was found

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
            "cdna_notation": "string",
            "protein_notation": "string",
            "genomic_position": "string or null",
            "clinical_significance": "string",
            "patients": {{
                "count": integer,
                "demographics": "string",
                "phenotype": "string"
            }},
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
        model: str = "gpt-4o",  # Use more capable model for extraction
        temperature: float = 0.0,
        max_tokens: int = 4000
    ):
        """
        Initialize the Expert Extractor.

        Args:
            model: LiteLLM model identifier (e.g., 'gpt-4o', 'claude-3-opus-20240229').
            temperature: Model temperature (0.0 for most deterministic).
            max_tokens: Maximum tokens for response.
        """
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens

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

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=2, min=4, max=30))
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
            # Call LiteLLM with large context model
            response = completion(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                response_format={"type": "json_object"}
            )

            # Parse response
            result_text = response.choices[0].message.content
            extracted_data = json.loads(result_text)

            # Get token usage
            tokens_used = None
            if hasattr(response, 'usage'):
                tokens_used = response.usage.total_tokens

            logger.info(
                f"PMID {paper.pmid} - Extraction successful. "
                f"Found {extracted_data.get('extraction_metadata', {}).get('total_variants_found', 0)} variants. "
                f"Tokens used: {tokens_used}"
            )

            return ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=extracted_data,
                model_used=self.model,
                tokens_used=tokens_used
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


def extract_variants_from_paper(paper: Paper, model: str = "gpt-4o") -> ExtractionResult:
    """
    Convenience function to extract variants from a single paper.

    Args:
        paper: Paper object with full text.
        model: LiteLLM model identifier.

    Returns:
        ExtractionResult.
    """
    extractor = ExpertExtractor(model=model)
    return extractor.extract(paper)
