"""
Filter modules (The Sieves) for the Tiered Biomedical Extraction Pipeline.

Tier 1: KeywordFilter - Fast, cheap keyword-based filtering
Tier 2: InternFilter - LLM-based classification (gpt-4o-mini)
"""

import logging
import re
from typing import List, Optional
from models import Paper, FilterResult, FilterDecision, FilterTier
from pipeline.utils.llm_utils import BaseLLMCaller
from config.settings import get_settings

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


class InternFilter(BaseLLMCaller):
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
        model: Optional[str] = None,
        temperature: float = 0.1,
        max_tokens: int = 150
    ):
        """
        Initialize the Intern filter.

        Args:
            model: LiteLLM model identifier. If None, uses config settings.
            temperature: Model temperature (lower = more deterministic).
            max_tokens: Maximum tokens for response.
        """
        settings = get_settings()
        model = model or settings.intern_model
        super().__init__(model=model, temperature=temperature, max_tokens=max_tokens)

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

            # Use shared LLM utility (includes retry logic)
            result_data = self.call_llm_json(prompt)

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
                    "model": self.model
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


class ClinicalDataTriageFilter(BaseLLMCaller):
    """
    Clinical Data Triage Filter: Determines if text contains ORIGINAL clinical data.

    This is a specialized filter that triages papers based on whether they contain
    new clinical case data vs. reviews, animal studies, or other non-clinical research.

    Rules:
    - REJECT: Review articles, Meta-analyses (unless individual patient data attached),
              Animal studies, Cell studies only, Variant interpretation guidelines without new cases
    - ACCEPT: Case reports, Case series, Clinical cohorts, Functional studies that also describe phenotype

    Output: JSON with "KEEP" or "DROP" decision, reason, and confidence score.
    """

    TRIAGE_PROMPT = """You are a Triage Assistant for clinical genetic research papers.

Task: Determine if the text contains ORIGINAL clinical data (case reports, cohort studies) for {gene}.
Input: Abstract/Introduction.

Rules:
1. REJECT if:
   - Review article
   - Meta-analysis (unless individual patient data attached)
   - Animal study (mouse, rat, zebrafish, etc.)
   - Cell study only (in vitro, cell lines, no patients)
   - Variant interpretation guidelines (without new cases)
   - Purely computational/bioinformatics study
   - Basic science without patient phenotypes

2. ACCEPT if:
   - Case report (even single patient)
   - Case series (multiple patients)
   - Clinical cohort study
   - Functional study that *also* describes the patient phenotype
   - Clinical trial with patient-level data
   - Family study with clinical data

Key Question: Does this paper report NEW patient-level clinical data?

Title: {title}

Abstract/Introduction: {abstract}

Output format: JSON only.
{{
  "decision": "KEEP" | "DROP",
  "reason": "Brief explanation (e.g., 'Review article only', 'Mouse model only', 'New case report detected')",
  "confidence": 0.0-1.0
}}

Respond ONLY with valid JSON. Be conservative - when in doubt about borderline cases, use confidence < 0.5."""

    def __init__(
        self,
        model: Optional[str] = None,
        temperature: float = 0.1,
        max_tokens: int = 200
    ):
        """
        Initialize the Clinical Data Triage filter.

        Args:
            model: LiteLLM model identifier (default: gpt-4o-mini for cost efficiency). If None, uses config.
            temperature: Model temperature (lower = more deterministic).
            max_tokens: Maximum tokens for response.
        """
        settings = get_settings()
        model = model or settings.intern_model
        super().__init__(model=model, temperature=temperature, max_tokens=max_tokens)

    def triage(
        self,
        title: str,
        abstract: str,
        gene: str = "the gene of interest",
        pmid: Optional[str] = None
    ) -> dict:
        """
        Triage a paper based on title and abstract.

        Args:
            title: Paper title.
            abstract: Paper abstract or introduction text.
            gene: Gene symbol being studied.
            pmid: Optional PMID for logging.

        Returns:
            Dictionary with:
                - decision: "KEEP" or "DROP"
                - reason: Explanation for the decision
                - confidence: Confidence score (0.0-1.0)
                - pmid: PMID if provided
        """
        if not abstract or not title:
            logger.warning(f"Missing title or abstract for triage{f' (PMID: {pmid})' if pmid else ''}")
            return {
                "decision": "DROP",
                "reason": "Missing title or abstract",
                "confidence": 1.0,
                "pmid": pmid
            }

        # Construct prompt
        prompt = self.TRIAGE_PROMPT.format(
            gene=gene,
            title=title,
            abstract=abstract[:2500]  # Truncate very long abstracts
        )

        try:
            logger.debug(f"Triaging paper{f' PMID {pmid}' if pmid else ''} for gene {gene}")

            # Use shared LLM utility (includes retry logic)
            result_data = self.call_llm_json(prompt)

            # Validate and normalize decision
            decision = result_data.get("decision", "DROP").upper()
            if decision not in ["KEEP", "DROP"]:
                logger.warning(f"Invalid decision '{decision}', defaulting to DROP")
                decision = "DROP"

            reason = result_data.get("reason", "No reason provided")
            confidence = float(result_data.get("confidence", 0.5))
            confidence = max(0.0, min(1.0, confidence))  # Clamp to [0, 1]

            result = {
                "decision": decision,
                "reason": reason,
                "confidence": confidence,
                "pmid": pmid,
                "model": self.model
            }

            logger.info(
                f"Triage result{f' for PMID {pmid}' if pmid else ''}: {decision} "
                f"(confidence: {confidence:.2f}) - {reason}"
            )

            return result

        except Exception as e:
            logger.error(f"Triage failed{f' for PMID {pmid}' if pmid else ''}: {e}")
            # On error, drop the paper to be conservative
            return {
                "decision": "DROP",
                "reason": f"Triage error: {str(e)}",
                "confidence": 0.0,
                "pmid": pmid,
                "error": str(e)
            }

    def triage_paper(self, paper: Paper, gene: Optional[str] = None) -> dict:
        """
        Triage a Paper object.

        Args:
            paper: Paper object with title and abstract.
            gene: Optional gene symbol (uses paper.gene_symbol if not provided).

        Returns:
            Dictionary with triage results.
        """
        gene_symbol = gene or paper.gene_symbol or "the gene of interest"

        return self.triage(
            title=paper.title or "",
            abstract=paper.abstract or "",
            gene=gene_symbol,
            pmid=paper.pmid
        )