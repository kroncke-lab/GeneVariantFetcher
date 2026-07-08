"""
Filter modules (The Sieves) for the Tiered Biomedical Extraction Pipeline.

Tier 1: KeywordFilter - Fast, cheap keyword-based filtering
Tier 2: InternFilter - LLM-based classification (gpt-4o-mini)
"""

import logging
import re
from typing import List, Optional

from config.constants import FILTER_CLINICAL_KEYWORDS
from config.settings import get_settings
from utils.gene_metadata import get_gene_aliases
from utils.llm_utils import BaseLLMCaller
from utils.models import FilterDecision, FilterResult, FilterTier, Paper

logger = logging.getLogger(__name__)


# Substrings that flag a model as a reasoning model (large hidden reasoning
# pass before the visible response). These need a higher max_tokens floor or
# they emit empty/truncated content. Match against the lower-cased model id.
_REASONING_MODEL_HINTS = (
    "kimi",
    "grok-4",
    "grok-reasoning",
    "gpt-5",
    "gpt5",
    "o1",
    "o3",
)

GENOTYPED_COHORT_SIGNALS = (
    "MUTATION",
    "MUTATIONS",
    "VARIANT",
    "VARIANTS",
    "CARRIER",
    "CARRIERS",
    "GENOTYP",
    "SEQUENC",
    "SCREEN",
    "SCREENED",
    "SCREENING",
    "IDENTIFIED",
    "FOUND IN",
    "PROBAND",
    "PROBANDS",
    "FAMILY",
    "FAMILIES",
    "PATIENT",
    "PATIENTS",
    "COHORT",
    "REGISTRY",
)

REVIEW_ONLY_SIGNALS = (
    "REVIEW",
    "META-ANALYSIS",
    "SYSTEMATIC REVIEW",
    "GUIDELINE",
    "CONSENSUS",
)


def _is_reasoning_model(model: Optional[str]) -> bool:
    if not model:
        return False
    m = model.lower()
    return any(hint in m for hint in _REASONING_MODEL_HINTS)


class KeywordFilter:
    """
    Tier 1 Filter: Fast keyword-based filtering of abstracts.
    Checks for clinical/variant-related keywords to reduce obviously irrelevant papers.
    """

    # Use centralized constant for keywords
    DEFAULT_CLINICAL_KEYWORDS = FILTER_CLINICAL_KEYWORDS

    def __init__(
        self,
        keywords: Optional[List[str]] = None,
        min_keyword_matches: Optional[int] = None,
    ):
        """
        Initialize the keyword filter.

        Args:
            keywords: List of keywords to search for. If None, uses FILTER_CLINICAL_KEYWORDS from config.
            min_keyword_matches: Minimum number of keyword matches required to pass. If None, uses config default.
        """
        settings = get_settings()

        self.keywords = keywords or FILTER_CLINICAL_KEYWORDS
        self.min_keyword_matches = (
            min_keyword_matches
            if min_keyword_matches is not None
            else settings.tier1_min_keywords
        )

        # Compile single combined regex pattern for 10x+ speedup
        # Uses alternation (|) to match any keyword in one pass
        escaped_keywords = [re.escape(keyword) for keyword in self.keywords]
        combined_pattern = r"\b(" + "|".join(escaped_keywords) + r")\b"
        self.pattern = re.compile(combined_pattern, re.IGNORECASE)

        logger.debug(
            f"KeywordFilter initialized with {len(self.keywords)} keywords, min_matches={self.min_keyword_matches}"
        )

    def filter(self, paper: Paper) -> FilterResult:
        """
        Apply keyword filter to a paper's abstract.

        Args:
            paper: Paper object with abstract.

        Returns:
            FilterResult with decision and reason.
        """
        if not paper.abstract:
            logger.warning(
                f"PMID {paper.pmid} has no abstract, passing through (fail-open)"
            )
            return FilterResult(
                decision=FilterDecision.PASS,
                tier=FilterTier.TIER_1_KEYWORD,
                reason="No abstract available — fail-open to Tier 2",
                pmid=paper.pmid,
                confidence=0.0,
                metadata={"matched_keywords": [], "fail_open": True},
            )

        # Find all matching keywords in one pass (10x faster than individual searches)
        matches = self.pattern.findall(paper.abstract.lower())
        # Get unique matches (case-insensitive)
        matched_keywords = list(
            dict.fromkeys(matches)
        )  # Preserves order, removes duplicates
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

        logger.debug(
            f"PMID {paper.pmid} - Keyword filter: {decision.value} ({num_matches} matches)"
        )

        return FilterResult(
            decision=decision,
            tier=FilterTier.TIER_1_KEYWORD,
            reason=reason,
            pmid=paper.pmid,
            confidence=min(num_matches / 10.0, 1.0),  # Simple confidence score
            metadata={"matched_keywords": matched_keywords, "num_matches": num_matches},
        )


class InternFilter(BaseLLMCaller):
    """
    Tier 2 Filter: LLM-based classification using gpt-4o-mini.
    Classifies papers as "Original Clinical Data" vs "Review/Irrelevant".
    """

    CLASSIFICATION_PROMPT = """You are a medical research classifier. Your job is to decide whether downstream variant curation could plausibly extract one or more genetic variants for {gene_clause} from this paper.

Remember: you only see the abstract. The full text and tables are available downstream and almost always contain the specific HGVS identifiers when the abstract describes a clinical genotyped cohort, family study, or case report. Lean toward PASS when the abstract makes it likely the body has variants.

Classify the paper as ONE of:
1. PASS - Likely contains extractable variants. This includes:
   - Case reports, case series, family/pedigree studies, or clinical cohorts of genotyped patients (e.g. "we describe a family with a mutation in the target gene", "26 carriers of a target-gene mutation", "patients with variants in the target gene")
   - Studies that screen patients with the target disease or phenotype for variants — even if the abstract only names the gene/disease/phenotype and not the specific HGVS
   - Functional / in-vitro / iPSC / heterologous-expression studies that name one or more specific variants (p.Arg123Trp / R123W / c.1234A>G) OR reference patient-derived variants
   - Pharmacogenomics or drug-induced phenotype papers in genotyped patients
   - NGS panel / diagnostic-yield papers that report results from a patient cohort
2. FAIL only if at least one of:
   - Pure review or meta-analysis with no new patient/functional data
   - Animal-only or non-human cell study with no specific human variant
   - Bioinformatics/computational paper that does not report specific variants
   - Methodology paper that explicitly does NOT include patient findings (e.g. assay-development only)
   - Unrelated topic (gene/disease not relevant)

Title: {title}

Abstract: {abstract}

Decision rules:
- "A mutation in the target gene", "carriers of a target-gene mutation", "target-gene positive patients" — PASS. The variant identifier is often in the body or supplement.
- "Patients with the target disease/phenotype were genotyped" — PASS. Specific variants may be in tables.
- Family/pedigree, twins, pregnancy/perinatal, pediatric/neonatal, adult-onset, tumor, population, diagnostic-yield, or cohort studies — PASS when genotyped patients are mentioned and the topic is relevant to the target gene/disease.
- "In vitro" or "functional study" alone is NOT a reason to fail — check whether named variants or patient-derived variants appear.
- When uncertain, prefer PASS. Cost of false-negative (lost relevant variant-bearing papers) >> cost of false-positive (wasted Tier-3 call).
{disease_clause}
Respond with a JSON object containing exactly these keys:
{{
    "decision": "PASS" or "FAIL",
    "reason": "Brief explanation (1-2 sentences). If FAIL, state which exclusion rule applied.",
    "confidence": 0.0-1.0
}}

Output the JSON object FIRST, before any reasoning. Keep the response under 200 tokens."""

    DISEASE_PROMPT_ADDENDUM = (
        "- PRIORITIZE papers reporting original patient or functional data for {disease}.\n"
        "- FAIL papers that are pure reviews of {disease} genetics without new patient/functional data.\n"
    )

    def __init__(
        self,
        model: Optional[str] = None,
        temperature: Optional[float] = None,
        max_tokens: Optional[int] = None,
        reasoning_effort: Optional[str] = None,
        confidence_threshold: Optional[float] = None,
        disease: Optional[str] = None,
    ):
        """
        Initialize the Intern filter.

        Args:
            model: LiteLLM model identifier. If None, uses config settings (TIER2_MODEL).
            temperature: Model temperature (lower = more deterministic). If None, uses config.
            max_tokens: Maximum tokens for response. If None, uses config.
            confidence_threshold: Minimum confidence to pass (0.0-1.0). If None, uses config.
            disease: Optional disease term (e.g. "atrial fibrillation"). When set,
                an addendum is added to the classification prompt prioritizing
                original patient/functional data for that disease and rejecting
                pure reviews. When None (default), the prompt is unchanged.
        """
        settings = get_settings()

        # Use config defaults if not specified
        model = model or settings.get_tier2_model() or settings.intern_model
        temperature = (
            temperature if temperature is not None else settings.tier2_temperature
        )
        max_tokens = max_tokens if max_tokens is not None else settings.tier2_max_tokens
        reasoning_effort = (
            reasoning_effort
            if reasoning_effort is not None
            else settings.tier2_reasoning_effort
        )

        # Reasoning models (Kimi-K2.6, grok-reasoning, gpt-5-codex) burn most
        # of the token budget on hidden reasoning before emitting visible JSON.
        # 2.5k is borderline and produces empty/truncated output ~10% of the
        # time — the JSON-repair path saves the parse but yields an empty `{}`
        # that bypasses the schema and silently fails papers. Lift the floor
        # to 8192 for any reasoning model (matches the table-router baseline).
        if max_tokens < 8192 and _is_reasoning_model(model):
            logger.info(
                f"InternFilter: bumping max_tokens {max_tokens} -> 8192 "
                f"for reasoning model {model!r} (avoids empty-JSON truncation)"
            )
            max_tokens = 8192

        self.confidence_threshold = (
            confidence_threshold
            if confidence_threshold is not None
            else settings.tier2_confidence_threshold
        )
        self.disease = disease.strip() if disease and disease.strip() else None

        super().__init__(
            model=model,
            temperature=temperature,
            max_tokens=max_tokens,
            reasoning_effort=reasoning_effort,
        )

        logger.debug(
            f"InternFilter initialized with model={model}, temp={temperature}, confidence_threshold={self.confidence_threshold}, disease={self.disease!r}"
        )

    @staticmethod
    def _contains_alias(text: str, alias: str) -> bool:
        pattern = rf"(?<![A-Z0-9]){re.escape(alias.upper())}(?![A-Z0-9])"
        return re.search(pattern, text) is not None

    def _target_gene_or_alias_mentioned(self, paper: Paper) -> bool:
        gene = (paper.gene_symbol or "").strip().upper()
        if not gene:
            return False
        text = f"{paper.title or ''}\n{paper.abstract or ''}".upper()
        aliases = get_gene_aliases(gene, include_query_aliases=True)
        return any(self._contains_alias(text, alias) for alias in aliases)

    def _should_fail_open_target_gene_abstract(self, paper: Paper) -> bool:
        """Recover papers the LLM labels FAIL despite target-gene cohort signals."""
        text = f"{paper.title or ''}\n{paper.abstract or ''}".upper()
        if not self._target_gene_or_alias_mentioned(paper):
            return False
        if not any(signal in text for signal in GENOTYPED_COHORT_SIGNALS):
            return False
        if (
            any(signal in text for signal in REVIEW_ONLY_SIGNALS)
            and "CASE" not in text
            and "IDENTIFIED" not in text
            and "PATIENT" not in text
        ):
            return False
        return True

    def filter(self, paper: Paper) -> FilterResult:
        """
        Apply LLM-based filter to classify the paper.

        Args:
            paper: Paper object with title and abstract.

        Returns:
            FilterResult with LLM-based decision.
        """
        if not paper.abstract or not paper.title:
            logger.warning(
                f"PMID {paper.pmid} missing title/abstract for Intern filter"
            )
            return FilterResult(
                decision=FilterDecision.PASS,
                tier=FilterTier.TIER_2_INTERN,
                reason="Missing title or abstract for LLM classification; fail-open for recall",
                pmid=paper.pmid,
                confidence=0.0,
                metadata={
                    "model": self.model,
                    "fail_open": True,
                    "missing_title_or_abstract": True,
                },
            )

        # Construct prompt — disease clause is empty unless --disease was set
        disease_clause = (
            self.DISEASE_PROMPT_ADDENDUM.format(disease=self.disease)
            if self.disease
            else ""
        )
        prompt = self.CLASSIFICATION_PROMPT.format(
            gene_clause=paper.gene_symbol or "the target gene",
            title=paper.title,
            abstract=paper.abstract[:2000],  # Truncate very long abstracts
            disease_clause=disease_clause,
        )

        try:
            logger.debug(f"PMID {paper.pmid} - Calling LLM for classification")

            # Use shared LLM utility (includes retry logic)
            result_data = self.call_llm_json(prompt)

            raw_decision = result_data.get("decision")
            decision_str = (
                raw_decision.strip().upper() if isinstance(raw_decision, str) else ""
            )

            # Reasoning models (e.g. Kimi-K2.6) sometimes burn the token budget on
            # hidden reasoning and emit empty/keyless JSON that parses fine but is
            # missing the schema. The previous default-to-FAIL behavior silently
            # killed ~320 papers (incl. ~30 gold-standard variants) per KCNH2 run.
            # Treat unrecognized output as inconclusive and fail-OPEN, matching the
            # exception path's policy.
            if decision_str not in {"PASS", "FAIL"}:
                logger.warning(
                    f"PMID {paper.pmid} - Intern filter: malformed/empty LLM "
                    f"output (keys={list(result_data.keys())}); failing OPEN."
                )
                return FilterResult(
                    decision=FilterDecision.PASS,
                    tier=FilterTier.TIER_2_INTERN,
                    reason="Inconclusive LLM output (missing/invalid decision key); fail-open",
                    pmid=paper.pmid,
                    confidence=0.0,
                    metadata={
                        "model": self.model,
                        "fail_open": True,
                        "raw_keys": list(result_data.keys()),
                    },
                )

            reason = result_data.get("reason") or "No reason provided"
            try:
                confidence = float(result_data.get("confidence", 0.5))
            except (TypeError, ValueError):
                confidence = 0.5

            # Apply confidence threshold - if confidence is below threshold, fail the paper
            if decision_str == "PASS" and confidence < self.confidence_threshold:
                decision = FilterDecision.FAIL
                reason = f"Low confidence ({confidence:.2f} < {self.confidence_threshold}): {reason}"
                logger.info(
                    f"PMID {paper.pmid} - Failed due to low confidence: {confidence:.2f}"
                )
            else:
                decision = (
                    FilterDecision.PASS
                    if decision_str == "PASS"
                    else FilterDecision.FAIL
                )

            fail_open_target_gene_signal = False
            if (
                decision == FilterDecision.FAIL
                and self._should_fail_open_target_gene_abstract(paper)
            ):
                fail_open_target_gene_signal = True
                decision = FilterDecision.PASS
                reason = (
                    "Fail-open override: abstract/title explicitly mention the "
                    f"target gene or alias for {paper.gene_symbol} plus genotyped "
                    f"cohort/variant signals. Original classifier reason: {reason}"
                )

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
                    "confidence_threshold": self.confidence_threshold,
                    "fail_open_target_gene_signal": fail_open_target_gene_signal,
                },
            )

        except Exception as e:
            logger.error(f"PMID {paper.pmid} - Intern filter failed: {e}")
            # On error, fail-open (PASS) to avoid losing papers due to transient LLM issues
            return FilterResult(
                decision=FilterDecision.PASS,
                tier=FilterTier.TIER_2_INTERN,
                reason=f"LLM error (fail-open): {str(e)}",
                pmid=paper.pmid,
                confidence=0.0,
                metadata={"error": str(e), "fail_open": True},
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

    TRIAGE_PROMPT = """You are a Triage Assistant for genetic variant curation. Your job is to decide whether the paper could plausibly contain extractable patient, family, cohort, functional, or table/supplement variant data for {gene} (or its protein).

Input: Abstract/Introduction.

Rules:
1. REJECT if:
   - Pure review or meta-analysis with NO new variant findings
   - Animal-only study (mouse, rat, zebrafish) with NO specific human variant identified
   - Variant interpretation guidelines / methodology paper with NO new cases, cohorts, functional variant assays, or patient findings
   - Purely computational/bioinformatics study with NO specific variants reported
   - Off-topic (unrelated to {gene})

2. ACCEPT if:
   - Case report, even a single patient, with a named variant or clear genotyping
   - Case series, clinical cohort, screening study, diagnostic-yield study, or panel-sequencing study likely to report variants in full text/tables
   - Family study with clinical data and variant/genotype data
   - Clinical trial or pharmacogenomic study with patient-level variant data
   - **Functional/in-vitro study (cell lines, iPSC, heterologous expression, electrophysiology) that names one or more SPECIFIC variants** (e.g. p.Arg123Trp, R123W, c.1234A>G). These are still primary sources for variant data.
   - Re-analysis of public sequencing data (gnomAD, ESP, exome cohorts) reporting specific variants
{disease_clause}
Key Question: Should downstream full-text/supplement extraction spend effort on this paper?
Do NOT require the abstract itself to name specific variants. Cohort and screening abstracts often hide variant lists in full-text tables or supplements.
"In vitro" or "functional study" alone is NOT a reason to drop — only drop if no patient-derived or specific human variant data is likely.

Title: {title}

Abstract/Introduction: {abstract}

Output format: JSON only.
{{
  "decision": "KEEP" | "DROP",
  "reason": "Brief explanation (e.g., 'Review article only', 'Mouse model only', 'New case report detected')",
  "confidence": 0.0-1.0
}}

Respond ONLY with valid JSON. Be recall-biased: when in doubt, KEEP with low confidence rather than DROP."""

    DISEASE_PROMPT_ADDENDUM = (
        "\n3. DISEASE-AWARE PRIORITIZATION (this run is filtering for {disease}):\n"
        "   - PRIORITIZE original patient or functional data for {disease}.\n"
        "   - REJECT papers that are pure reviews of {disease} genetics without new cases.\n"
    )

    def __init__(
        self,
        model: Optional[str] = None,
        temperature: Optional[float] = None,
        max_tokens: Optional[int] = None,
        reasoning_effort: Optional[str] = None,
        disease: Optional[str] = None,
    ):
        """
        Initialize the Clinical Data Triage filter.

        Args:
            model: LiteLLM model identifier (default: uses TIER2_MODEL from config). If None, uses config.
            temperature: Model temperature (lower = more deterministic). If None, uses config.
            max_tokens: Maximum tokens for response. If None, uses config default (200).
            disease: Optional disease term. When set, an addendum is added to
                the triage prompt prioritizing original patient/functional data
                for that disease. When None (default), the prompt is unchanged.
        """
        settings = get_settings()

        model = model or settings.get_tier2_model() or settings.intern_model
        temperature = (
            temperature if temperature is not None else settings.tier2_temperature
        )
        max_tokens = (
            max_tokens if max_tokens is not None else 200
        )  # Slightly higher than default for triage
        reasoning_effort = (
            reasoning_effort
            if reasoning_effort is not None
            else settings.tier2_reasoning_effort
        )

        self.disease = disease.strip() if disease and disease.strip() else None

        super().__init__(
            model=model,
            temperature=temperature,
            max_tokens=max_tokens,
            reasoning_effort=reasoning_effort,
        )

        logger.debug(
            f"ClinicalDataTriageFilter initialized with model={model}, disease={self.disease!r}"
        )

    def triage(
        self,
        title: str,
        abstract: str,
        gene: str = "the gene of interest",
        pmid: Optional[str] = None,
    ) -> str:
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
            logger.warning(
                f"Missing title or abstract for triage{f' (PMID: {pmid})' if pmid else ''}"
            )
            return {
                "decision": "KEEP",
                "reason": "Missing title or abstract; fail-open for recall",
                "confidence": 0.0,
                "pmid": pmid,
                "fail_open": True,
            }

        # Construct prompt — disease clause is empty unless --disease was set
        disease_clause = (
            self.DISEASE_PROMPT_ADDENDUM.format(disease=self.disease)
            if self.disease
            else ""
        )
        prompt = self.TRIAGE_PROMPT.format(
            gene=gene,
            title=title,
            abstract=abstract[:2500],  # Truncate very long abstracts
            disease_clause=disease_clause,
        )

        try:
            logger.debug(
                f"Triaging paper{f' PMID {pmid}' if pmid else ''} for gene {gene}"
            )

            # Use shared LLM utility (includes retry logic)
            result_data = self.call_llm_json(prompt)

            # Validate and normalize decision
            decision = result_data.get("decision", "DROP").upper()
            if decision not in ["KEEP", "DROP"]:
                logger.warning(f"Invalid decision '{decision}', fail-opening to KEEP")
                decision = "KEEP"

            reason = result_data.get("reason", "No reason provided")
            confidence = float(result_data.get("confidence", 0.5))
            confidence = max(0.0, min(1.0, confidence))  # Clamp to [0, 1]

            result = {
                "decision": decision,
                "reason": reason,
                "confidence": confidence,
                "pmid": pmid,
                "model": self.model,
            }

            logger.info(
                f"Triage result{f' for PMID {pmid}' if pmid else ''}: {decision} "
                f"(confidence: {confidence:.2f}) - {reason}"
            )

            return result

        except Exception as e:
            logger.error(f"Triage failed{f' for PMID {pmid}' if pmid else ''}: {e}")
            # On error, fail-open (KEEP) to avoid losing papers due to transient LLM issues
            return {
                "decision": "KEEP",
                "reason": f"Triage error (fail-open): {str(e)}",
                "confidence": 0.0,
                "pmid": pmid,
                "error": str(e),
                "fail_open": True,
            }

    def triage_paper(self, paper: Paper, gene: Optional[str] = None) -> str:
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
            pmid=paper.pmid,
        )
