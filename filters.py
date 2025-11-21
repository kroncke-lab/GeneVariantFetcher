"""
Filter modules (The Sieves) for the Tiered Biomedical Extraction Pipeline.

Tier 1: KeywordFilter - Fast, cheap keyword-based filtering
Tier 2: InternFilter - LLM-based classification (gpt-4o-mini)
"""

import logging
import re
from typing import List, Optional
from litellm import completion
from tenacity import retry, stop_after_attempt, wait_exponential
from models import Paper, FilterResult, FilterDecision, FilterTier

logger = logging.getLogger(__name__)


class KeywordFilter:
    """
    Tier 1 Filter: Fast keyword-based filtering of abstracts.
    Checks for clinical/variant-related keywords to reduce obviously irrelevant papers.
    """

    DEFAULT_CLINICAL_KEYWORDS = [
        # Variant/mutation terms
        "variant", "mutation", "polymorphism", "SNP", "deletion", "insertion",
        "substitution", "missense", "nonsense", "frameshift", "splice",
        "indel", "copy number", "CNV",

        # Clinical terms
        "patient", "patients", "clinical", "disease", "syndrome", "phenotype",
        "diagnosis", "treatment", "therapy", "outcome", "prognosis",
        "symptom", "symptoms", "manifestation", "pathogenic", "benign",

        # Study types
        "case report", "case series", "cohort", "clinical trial", "study",
        "analysis", "association", "genotype", "phenotype",

        # Medical/genetic terms
        "pathology", "molecular", "genetic", "genomic", "exome", "sequencing",
        "gene", "chromosome", "allele", "heterozygous", "homozygous",
        "carrier", "inheritance", "familial", "sporadic"
    ]

    def __init__(self, keywords: Optional[List[str]] = None, min_keyword_matches: int = 2):
        """
        Initialize the keyword filter.

        Args:
            keywords: List of keywords to search for. If None, uses DEFAULT_CLINICAL_KEYWORDS.
            min_keyword_matches: Minimum number of keyword matches required to pass.
        """
        self.keywords = keywords or self.DEFAULT_CLINICAL_KEYWORDS
        self.min_keyword_matches = min_keyword_matches

        # Compile regex patterns for efficiency
        self.patterns = [
            re.compile(r'\b' + re.escape(keyword) + r'\b', re.IGNORECASE)
            for keyword in self.keywords
        ]

    def filter(self, paper: Paper) -> FilterResult:
        """
        Apply keyword filter to a paper's abstract.

        Args:
            paper: Paper object with abstract.

        Returns:
            FilterResult with decision and reason.
        """
        if not paper.abstract:
            logger.warning(f"PMID {paper.pmid} has no abstract, failing keyword filter")
            return FilterResult(
                decision=FilterDecision.FAIL,
                tier=FilterTier.TIER_1_KEYWORD,
                reason="No abstract available",
                pmid=paper.pmid,
                confidence=1.0,
                metadata={"matched_keywords": []}
            )

        # Find matching keywords
        matched_keywords = []
        for keyword, pattern in zip(self.keywords, self.patterns):
            if pattern.search(paper.abstract):
                matched_keywords.append(keyword)

        num_matches = len(matched_keywords)

        # Determine pass/fail
        if num_matches >= self.min_keyword_matches:
            decision = FilterDecision.PASS
            reason = f"Found {num_matches} clinical keywords: {', '.join(matched_keywords[:5])}"
            if num_matches > 5:
                reason += f" (and {num_matches - 5} more)"
        else:
            decision = FilterDecision.FAIL
            reason = f"Only {num_matches} keyword match(es), minimum {self.min_keyword_matches} required"

        logger.debug(f"PMID {paper.pmid} - Keyword filter: {decision.value} ({num_matches} matches)")

        return FilterResult(
            decision=decision,
            tier=FilterTier.TIER_1_KEYWORD,
            reason=reason,
            pmid=paper.pmid,
            confidence=min(num_matches / 10.0, 1.0),  # Simple confidence score
            metadata={
                "matched_keywords": matched_keywords,
                "num_matches": num_matches
            }
        )


class InternFilter:
    """
    Tier 2 Filter: LLM-based classification using gpt-4o-mini.
    Classifies papers as "Original Clinical Data" vs "Review/Irrelevant".
    """

    CLASSIFICATION_PROMPT = """You are a medical research classifier. Your job is to determine if a scientific paper contains ORIGINAL CLINICAL DATA about genetic variants.

Classify the paper as ONE of:
1. PASS - Contains original clinical data (case reports, patient cohorts, clinical studies with genetic variant data)
2. FAIL - Review article, meta-analysis, methodology paper, or not clinically relevant

Title: {title}

Abstract: {abstract}

Instructions:
- Look for original patient data, case reports, or clinical studies
- PASS if the paper reports new genetic variant findings in patients
- FAIL if it's a review, meta-analysis, or purely methodological
- FAIL if it's basic science without clinical application

Respond with a JSON object:
{{
    "decision": "PASS" or "FAIL",
    "reason": "Brief explanation (1-2 sentences)",
    "confidence": 0.0-1.0
}}"""

    def __init__(
        self,
        model: str = "gpt-4o-mini",
        temperature: float = 0.1,
        max_tokens: int = 150
    ):
        """
        Initialize the Intern filter.

        Args:
            model: LiteLLM model identifier.
            temperature: Model temperature (lower = more deterministic).
            max_tokens: Maximum tokens for response.
        """
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def filter(self, paper: Paper) -> FilterResult:
        """
        Apply LLM-based filter to classify the paper.

        Args:
            paper: Paper object with title and abstract.

        Returns:
            FilterResult with LLM-based decision.
        """
        if not paper.abstract or not paper.title:
            logger.warning(f"PMID {paper.pmid} missing title/abstract for Intern filter")
            return FilterResult(
                decision=FilterDecision.FAIL,
                tier=FilterTier.TIER_2_INTERN,
                reason="Missing title or abstract for LLM classification",
                pmid=paper.pmid,
                confidence=1.0
            )

        # Construct prompt
        prompt = self.CLASSIFICATION_PROMPT.format(
            title=paper.title,
            abstract=paper.abstract[:2000]  # Truncate very long abstracts
        )

        try:
            logger.debug(f"PMID {paper.pmid} - Calling LLM for classification")

            # Call LiteLLM
            response = completion(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                response_format={"type": "json_object"}
            )

            # Parse response
            import json
            result_text = response.choices[0].message.content
            result_data = json.loads(result_text)

            decision_str = result_data.get("decision", "FAIL").upper()
            decision = FilterDecision.PASS if decision_str == "PASS" else FilterDecision.FAIL
            reason = result_data.get("reason", "No reason provided")
            confidence = float(result_data.get("confidence", 0.5))

            logger.info(
                f"PMID {paper.pmid} - Intern filter: {decision.value} "
                f"(confidence: {confidence:.2f}) - {reason}"
            )

            return FilterResult(
                decision=decision,
                tier=FilterTier.TIER_2_INTERN,
                reason=reason,
                pmid=paper.pmid,
                confidence=confidence,
                metadata={
                    "model": self.model,
                    "tokens_used": response.usage.total_tokens if hasattr(response, 'usage') else None
                }
            )

        except Exception as e:
            logger.error(f"PMID {paper.pmid} - Intern filter failed: {e}")
            # On error, fail the paper to be conservative
            return FilterResult(
                decision=FilterDecision.FAIL,
                tier=FilterTier.TIER_2_INTERN,
                reason=f"LLM classification error: {str(e)}",
                pmid=paper.pmid,
                confidence=0.0,
                metadata={"error": str(e)}
            )
