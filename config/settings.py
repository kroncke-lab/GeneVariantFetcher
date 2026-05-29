"""Validated configuration for the Gene Variant Fetcher pipeline."""

import logging
import os
from functools import lru_cache
from pathlib import Path
from typing import ClassVar, List, Optional, Union

from pydantic import AliasChoices, Field, field_validator, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_ENV_PATH = Path(__file__).resolve().parent.parent / ".env"
logger = logging.getLogger(__name__)


# Default Anthropic model strings. Used by the get_*_model() helpers when
# MODEL_PROVIDER=anthropic and the user has not set the per-tier env var.
# All routed through LiteLLM, which understands the "anthropic/" prefix.
ANTHROPIC_TIER2_DEFAULT = "anthropic/claude-haiku-4-5-20251001"
ANTHROPIC_TIER3_DEFAULT = "anthropic/claude-sonnet-4-6,anthropic/claude-opus-4-7"
ANTHROPIC_TABLE_ROUTER_DEFAULT = "anthropic/claude-haiku-4-5-20251001"
ANTHROPIC_VISION_DEFAULT = "anthropic/claude-sonnet-4-6"

# Tier-model env var names. Used to decide whether the user explicitly set
# a tier model (in which case --model-provider does not override it) vs.
# left it at the field default.
_TIER_ENV_VARS = {
    "tier2_model": "TIER2_MODEL",
    "tier3_models": "TIER3_MODELS",
    "tier3_adjudicator_models": "TIER3_ADJUDICATOR_MODELS",
    "table_router_model": "TABLE_ROUTER_MODEL",
    "vision_model": "VISION_MODEL",
}


class Settings(BaseSettings):
    """Centralized application settings loaded from environment variables."""

    # API Keys - all loaded from .env file
    openai_api_key: Optional[str] = Field(default=None, env="OPENAI_API_KEY")
    anthropic_api_key: Optional[str] = Field(default=None, env="ANTHROPIC_API_KEY")
    gemini_api_key: Optional[str] = Field(default=None, env="GEMINI_API_KEY")
    ncbi_email: Optional[str] = Field(default=None, env="NCBI_EMAIL")
    ncbi_api_key: Optional[str] = Field(default=None, env="NCBI_API_KEY")
    elsevier_api_key: Optional[str] = Field(default=None, env="ELSEVIER_API_KEY")
    elsevier_insttoken: Optional[str] = Field(
        default=None,
        env="ELSEVIER_INSTTOKEN",
        description="Elsevier institution token for subscription access",
    )
    wiley_api_key: Optional[str] = Field(default=None, env="WILEY_API_KEY")
    springer_api_key: Optional[str] = Field(default=None, env="SPRINGER_API_KEY")

    # Azure AI Foundry — primary LLM provider. Single resource hosts multiple
    # deployments accessed via LiteLLM model strings like "azure_ai/<deployment>".
    # LiteLLM reads AZURE_AI_API_KEY / AZURE_AI_API_BASE from the process
    # environment automatically; we expose them here so other modules can read
    # them via Settings and so missing-key validation can include them.
    azure_ai_api_key: Optional[str] = Field(default=None, env="AZURE_AI_API_KEY")
    azure_ai_api_base: Optional[str] = Field(
        default=None,
        env="AZURE_AI_API_BASE",
        description="Azure AI Foundry endpoint, e.g. https://<resource>.services.ai.azure.com",
    )
    azure_ai_api_version: Optional[str] = Field(
        default="2024-08-01-preview", env="AZURE_AI_API_VERSION"
    )
    azure_deployment_gpt5_codex: Optional[str] = Field(
        default=None, env="AZURE_DEPLOYMENT_GPT5_CODEX"
    )
    azure_deployment_kimi: Optional[str] = Field(
        default=None, env="AZURE_DEPLOYMENT_KIMI"
    )
    azure_deployment_grok: Optional[str] = Field(
        default=None, env="AZURE_DEPLOYMENT_GROK"
    )
    azure_deployment_gpt54: Optional[str] = Field(
        default=None,
        validation_alias=AliasChoices("GPT54_DEPLOYMENT", "AZURE_DEPLOYMENT_GPT54"),
        description="Azure AI Foundry GPT-5.4 deployment name",
    )
    azure_deployment_deepseek: Optional[str] = Field(
        default=None,
        validation_alias=AliasChoices(
            "DEEPSEEK_DEPLOYMENT", "AZURE_DEPLOYMENT_DEEPSEEK"
        ),
        description="Azure AI Foundry DeepSeek deployment name",
    )

    # Model provider selector. When set to "anthropic" (default), per-tier
    # model strings default to the ANTHROPIC_* values below — but only for
    # tiers where the corresponding TIER*_MODEL env var is unset. Setting an
    # explicit env var always wins, so each tier remains independently
    # configurable. Set MODEL_PROVIDER=azure (or pass --model-provider azure)
    # to fall back to the Azure AI Foundry deployments.
    model_provider: str = Field(
        default="anthropic",
        env="MODEL_PROVIDER",
        description=(
            "LLM provider selector: 'anthropic' (default), 'azure', or 'openai'."
            " 'anthropic' switches every tier to Claude defaults unless the"
            " specific TIER*_MODEL env var is also set."
        ),
    )

    # Model Configuration
    tier1_model: Optional[str] = Field(
        default=None,
        env="TIER1_MODEL",
        description="Optional LLM for Tier 1 (if using LLM-based Tier 1)",
    )
    tier2_model: str = Field(
        default="gpt-4o-mini",
        env="TIER2_MODEL",
        description="Model for Tier 2 classification (cheap)",
    )
    tier3_models: Union[str, List[str]] = Field(
        default="gpt-4o-mini,gpt-4o",
        env="TIER3_MODELS",
        description="Comma-separated list of models for Tier 3 extraction (e.g., 'gpt-4o-mini,gpt-4o')",
    )

    # Anthropic per-tier overrides. These take effect when MODEL_PROVIDER=anthropic
    # AND the matching TIER*_MODEL env var is unset. Configurable independently
    # via env so users can mix-and-match (e.g. Anthropic for Tier 3 only).
    anthropic_tier2_model: str = Field(
        default=ANTHROPIC_TIER2_DEFAULT,
        env="ANTHROPIC_TIER2_MODEL",
        description="Anthropic Tier 2 classifier model (used when MODEL_PROVIDER=anthropic)",
    )
    anthropic_tier3_models: str = Field(
        default=ANTHROPIC_TIER3_DEFAULT,
        env="ANTHROPIC_TIER3_MODELS",
        description="Anthropic Tier 3 extraction cascade (used when MODEL_PROVIDER=anthropic)",
    )
    anthropic_table_router_model: str = Field(
        default=ANTHROPIC_TABLE_ROUTER_DEFAULT,
        env="ANTHROPIC_TABLE_ROUTER_MODEL",
        description="Anthropic table-router model (used when MODEL_PROVIDER=anthropic)",
    )
    anthropic_vision_model: str = Field(
        default=ANTHROPIC_VISION_DEFAULT,
        env="ANTHROPIC_VISION_MODEL",
        description="Anthropic vision model for figure/pedigree analysis (used when MODEL_PROVIDER=anthropic)",
    )

    # Legacy alias for backward compatibility
    intern_model: str = Field(default="gpt-4o-mini", env="INTERN_MODEL")

    # Tiered Classification Configuration
    enable_tier1: bool = Field(
        default=True,
        env="ENABLE_TIER1",
        description="Enable Tier 1 keyword/heuristic filtering",
    )
    enable_tier2: bool = Field(
        default=True, env="ENABLE_TIER2", description="Enable Tier 2 LLM classification"
    )
    enable_tier3: bool = Field(
        default=True, env="ENABLE_TIER3", description="Enable Tier 3 expert extraction"
    )

    # Tier 1 Configuration
    tier1_min_keywords: int = Field(
        default=1,
        env="TIER1_MIN_KEYWORDS",
        description="Minimum keyword matches for Tier 1 to pass (lowered from 2 to 1 for fail-open recall)",
    )
    tier1_use_llm: bool = Field(
        default=False,
        env="TIER1_USE_LLM",
        description="Use lightweight LLM for Tier 1 instead of keywords",
    )

    # Tier 2 Configuration
    tier2_temperature: float = Field(
        default=0.1, env="TIER2_TEMPERATURE", description="Temperature for Tier 2 LLM"
    )
    tier2_max_tokens: int = Field(
        default=300,
        env="TIER2_MAX_TOKENS",
        description="Max tokens for Tier 2 LLM response (raised from 150 to 300 for the variant-aware prompt, which produces longer reasons that were getting truncated)",
    )
    tier2_confidence_threshold: float = Field(
        default=0.3,
        env="TIER2_CONFIDENCE_THRESHOLD",
        description="Minimum confidence to pass Tier 2",
    )
    tier2_reasoning_effort: Optional[str] = Field(
        default=None,
        env="TIER2_REASONING_EFFORT",
        description=(
            "OpenAI-style reasoning effort for Tier 2 relevance models "
            "(minimal|low|medium|high). None = provider default. Tier 2 is a "
            "high-volume yes/no classifier, so 'minimal' is the cost/latency "
            "saver for gpt-5-family models rather than a quality lever. Ignored "
            "by models without an effort knob."
        ),
    )

    # Tier 3 Configuration
    tier3_temperature: float = Field(
        default=0.0, env="TIER3_TEMPERATURE", description="Temperature for Tier 3 LLM"
    )
    tier3_reasoning_effort: Optional[str] = Field(
        default=None,
        env="TIER3_REASONING_EFFORT",
        description=(
            "OpenAI-style reasoning effort for Tier 3 extraction models "
            "(minimal|low|medium|high). None = provider default. This is the "
            "stage where unique-variant recall and carrier-count accuracy are "
            "decided, so it is the prime candidate for 'medium'/'high' — measure "
            "with scripts/run_recall_suite.py before keeping it on. Ignored by "
            "models without an effort knob."
        ),
    )
    tier3_max_tokens: int = Field(
        default=16000,
        env="TIER3_MAX_TOKENS",
        description="Max tokens for Tier 3 LLM response (increased for large variant tables)",
    )
    tier3_threshold: int = Field(
        default=1,
        env="TIER3_THRESHOLD",
        description="Try next model if first finds fewer variants than this (0 = only use first model)",
    )
    enable_tier3_ensemble_qa: bool = Field(
        default=True,
        env="ENABLE_TIER3_ENSEMBLE_QA",
        description=(
            "Run a compact second-model adjudication pass for high-risk "
            "extractions instead of sending every full paper to every model."
        ),
    )
    tier3_adjudicator_models: Union[str, List[str]] = Field(
        default="anthropic/claude-sonnet-4-6",
        env="TIER3_ADJUDICATOR_MODELS",
        description=(
            "Comma-separated adjudicator models used for high-risk Tier 3 "
            "outputs. These models see compact evidence packets, not full papers."
        ),
    )
    tier3_adjudication_risk_threshold: int = Field(
        default=2,
        env="TIER3_ADJUDICATION_RISK_THRESHOLD",
        description="Minimum extraction risk score required to trigger adjudication.",
    )
    tier3_evidence_packet_max_chars: int = Field(
        default=24000,
        env="TIER3_EVIDENCE_PACKET_MAX_CHARS",
        description="Maximum characters sent to adjudicator evidence packets.",
    )
    tier3_adjudication_max_tokens: int = Field(
        default=12000,
        env="TIER3_ADJUDICATION_MAX_TOKENS",
        description="Maximum output tokens for Tier 3 adjudicator responses.",
    )
    tier3_max_verifier_cards: int = Field(
        default=20,
        env="TIER3_MAX_VERIFIER_CARDS",
        description="Maximum per-variant evidence cards to verify in one extraction.",
    )
    count_guard_policy: str = Field(
        default="off",
        env="COUNT_GUARD_POLICY",
        description=(
            "Policy for the per-PMID count-outlier guard "
            "(pipeline/count_outlier_guard.py) applied to freshly extracted "
            "variants before they are saved (off|flag|clear). 'off' (default) is "
            "a strict no-op. 'flag' only annotates suspected study-wide-N reuse "
            "with count_outlier_flags metadata and leaves raw counts intact. "
            "'clear' also zeros the flagged counts (the raw value is preserved "
            "under the flag). Gold-free; works on new gene-diseases."
        ),
    )
    count_classifier_policy: str = Field(
        default="off",
        env="COUNT_CLASSIFIER_POLICY",
        description=(
            "Policy for the per-variant count classifier "
            "(pipeline/count_classifier.py) applied to freshly extracted "
            "variants before they are saved (off|flag|clear). 'off' (default) is "
            "a strict no-op. 'flag' only annotates counts whose LLM-declared "
            "count_provenance is not per_variant_carrier with "
            "count_classifier_flags metadata and leaves raw counts intact. "
            "'clear' also zeros those counts (raw preserved under the flag). "
            "Gold-free; silently skips variants without count_provenance."
        ),
    )
    strict_cohort_labels: bool = Field(
        default=False,
        env="GVF_STRICT_COHORT_LABELS",
        description=(
            "When True, the deterministic table column mapper "
            "(pipeline/table_router.py) refuses to map an ambiguous case/control "
            "cohort column (header matching only case/disease or "
            "control/healthy/normal) to affected/unaffected on case-control or "
            "assay-style tables, where such columns are cohort totals rather than "
            "per-variant counts. Default False preserves current behavior; this is "
            "a count-accuracy (MAE) lever that should be measured by re-extraction "
            "before being enabled by default. Gold-free."
        ),
    )
    early_debate_models: Union[str, List[str]] = Field(
        default="",
        env="EARLY_DEBATE_MODELS",
        description=(
            "Comma-separated early debate critic models. If unset, defaults to "
            "GPT54_DEPLOYMENT and DEEPSEEK_DEPLOYMENT via Azure."
        ),
    )

    # Table-router (router-first extraction): the LLM classifies which tables
    # contain variant data; a deterministic parser then reads the cells. Falls
    # back to the full-text Tier-3 path when no usable tables are detected.
    enable_table_router: bool = Field(
        default=True,
        env="ENABLE_TABLE_ROUTER",
        description="Try router-first table extraction before sending full text to Tier 3",
    )
    table_router_model: str = Field(
        default="azure_ai/Kimi-K2.6-1",
        env="TABLE_ROUTER_MODEL",
        description="LLM used to classify tables and emit column mappings",
    )
    table_router_max_tokens: int = Field(
        default=8192,
        env="TABLE_ROUTER_MAX_TOKENS",
        description=(
            "Max tokens for the table-router response. Kimi-K2.6 is a reasoning"
            " model that consumes the budget on hidden reasoning tokens before"
            " emitting visible content; 1-2k is too small for any non-trivial"
            " multi-table prompt and produces empty content with"
            " finish_reason='length'. 8k gives reliable JSON output."
        ),
    )
    table_router_reasoning_effort: Optional[str] = Field(
        default=None,
        env="TABLE_ROUTER_REASONING_EFFORT",
        description=(
            "OpenAI-style reasoning effort for the table-router model "
            "(minimal|low|medium|high). None = provider default. Table "
            "classification is a moderate disambiguation task; 'low'/'medium' is "
            "the sensible range to trial. Ignored by models without an effort knob."
        ),
    )

    # Paper Sourcing Configuration
    use_pubmind: bool = Field(
        default=True,
        env="USE_PUBMIND",
        description="Use PubMind as primary literature source",
    )
    use_pubmed: bool = Field(
        default=True,
        env="USE_PUBMED",
        description="Use PubMed API as additional source",
    )
    use_europepmc: bool = Field(
        default=False,
        env="USE_EUROPEPMC",
        description="Use EuropePMC as additional source",
    )
    pubmind_only: bool = Field(
        default=False,
        env="PUBMIND_ONLY",
        description="Use ONLY PubMind (ignore PubMed/EuropePMC)",
    )
    max_papers_per_source: int = Field(
        default=1500,
        env="MAX_PAPERS_PER_SOURCE",
        description="Max papers to fetch per source",
    )

    # Genetic Data Scout Configuration
    scout_enabled: bool = Field(
        default=True,
        env="SCOUT_ENABLED",
        description="Enable Genetic Data Scout to create condensed DATA_ZONES.md files",
    )
    scout_min_relevance: float = Field(
        default=0.3,
        env="SCOUT_MIN_RELEVANCE",
        description="Minimum relevance score (0.0-1.0) for zones to be included",
    )
    scout_max_zones: int = Field(
        default=30,
        env="SCOUT_MAX_ZONES",
        description="Maximum number of data zones to identify per paper",
    )
    scout_use_condensed: bool = Field(
        default=True,
        env="SCOUT_USE_CONDENSED",
        description="Prefer DATA_ZONES.md over FULL_CONTEXT.md for extraction",
    )

    # Extraction Tuning
    extraction_max_chars: int = Field(
        default=60_000,
        env="EXTRACTION_MAX_CHARS",
        description="Max characters sent to LLM in a single extraction prompt",
    )
    scanner_merge_confidence: float = Field(
        default=0.6,
        env="SCANNER_MERGE_CONFIDENCE",
        description="Minimum scanner confidence to merge into LLM results",
    )
    scanner_max_hints: int = Field(
        default=50,
        env="SCANNER_MAX_HINTS",
        description="Maximum scanner hints to include in LLM prompt",
    )
    max_workers: int = Field(
        default=8,
        env="MAX_WORKERS",
        description=(
            "Explicit global override for parallel workers. When set via the"
            " MAX_WORKERS env var it wins over the provider-aware defaults"
            " (anthropic_max_workers / azure_max_workers); otherwise the"
            " get_max_workers() helper picks by MODEL_PROVIDER."
        ),
    )

    # Provider-aware rate limit & concurrency. Anthropic's tier-4 RPM headroom
    # lets us push the pipeline harder than Azure's tighter per-deployment
    # quota. Set via env (ANTHROPIC_RPM / AZURE_RPM / ANTHROPIC_MAX_WORKERS /
    # AZURE_MAX_WORKERS) for runtime tuning. Explicit LLM_REQUESTS_PER_MINUTE
    # / MAX_WORKERS / FILTER_MAX_WORKERS env vars still take precedence at the
    # call sites that read them.
    anthropic_rpm: int = Field(
        default=1000,
        env="ANTHROPIC_RPM",
        description="Default LLM requests-per-minute when MODEL_PROVIDER=anthropic",
    )
    azure_rpm: int = Field(
        default=50,
        env="AZURE_RPM",
        description="Default LLM requests-per-minute when MODEL_PROVIDER=azure",
    )
    deepseek_rpm: int = Field(
        default=4,
        env="DEEPSEEK_RPM",
        description=(
            "DeepSeek deployment RPM cap. DeepSeek-V4-Pro is constrained by a "
            "20k tokens/minute quota, so this stays below the 20 requests/minute "
            "request cap for evidence-card prompts."
        ),
    )
    anthropic_max_workers: int = Field(
        default=20,
        env="ANTHROPIC_MAX_WORKERS",
        description="Default parallel workers (extraction) for Anthropic",
    )
    azure_max_workers: int = Field(
        default=3,
        env="AZURE_MAX_WORKERS",
        description="Default parallel workers (extraction) for Azure",
    )
    # Filter-stage workers default higher than extraction. Tier-2 LLM filter
    # calls are small (single abstract + 200-token JSON response) so Anthropic
    # comfortably handles 20+ concurrent. Extraction calls are 10–100× larger
    # and benefit from MAX_WORKERS=1–4. Keeping filter concurrency independent
    # means the filter stage hits ~200+ RPM even when the user sets
    # MAX_WORKERS=1 for extraction reasons.
    filter_max_workers: int = Field(
        default=20,
        env="FILTER_MAX_WORKERS",
        description=(
            "Default parallel workers for the Tier-2 LLM filter (cli/discover, "
            "filter step of extract). Independent of MAX_WORKERS used for "
            "extraction. Anthropic recommended: 20+."
        ),
    )

    # Figure/Image Extraction Configuration
    extract_figures: bool = Field(
        default=True,
        env="EXTRACT_FIGURES",
        description="Extract images from PDFs during harvesting",
    )
    extract_pedigrees: bool = Field(
        default=True,
        env="EXTRACT_PEDIGREES",
        description="Run pedigree detection and extraction on extracted figures",
    )
    vision_model: str = Field(
        default="azure_ai/gpt-5.3-codex-1",
        env="VISION_MODEL",
        description=(
            "Vision-capable model for pedigree analysis. Azure gpt-5 family"
            " deployments (e.g. gpt-5.3-codex-1) are routed through the"
            " Responses API automatically; other models go through the"
            " standard chat-completions path."
        ),
    )
    vision_reasoning_effort: Optional[str] = Field(
        default=None,
        env="VISION_REASONING_EFFORT",
        description=(
            "OpenAI-style reasoning effort for figure/pedigree vision models "
            "(minimal|low|medium|high). None = provider default. NOTE: the "
            "figure/pedigree call sites bypass BaseLLMCaller and call "
            "litellm.completion() directly, so this is inert until those sites "
            "are wired (see pedigree_extractor.py / figure_*_reader.py)."
        ),
    )

    @field_validator(
        "tier2_reasoning_effort",
        "tier3_reasoning_effort",
        "table_router_reasoning_effort",
        "vision_reasoning_effort",
    )
    @classmethod
    def _validate_reasoning_effort(cls, v: Optional[str]) -> Optional[str]:
        """Normalize reasoning-effort settings; reject unknown levels early."""
        if v is None:
            return None
        v_norm = str(v).strip().lower()
        if v_norm == "":
            return None
        allowed = {"minimal", "low", "medium", "high"}
        if v_norm not in allowed:
            raise ValueError(
                "reasoning_effort must be one of "
                f"{sorted(allowed)} or unset; got {v!r}"
            )
        return v_norm

    @field_validator("count_guard_policy", "count_classifier_policy")
    @classmethod
    def _validate_count_policy(cls, v: str) -> str:
        """Normalize count-hygiene policy settings; reject unknown modes early."""
        v_norm = str(v).strip().lower()
        if v_norm == "":
            return "off"
        allowed = {"off", "flag", "clear"}
        if v_norm not in allowed:
            raise ValueError(
                f"count policy must be one of {sorted(allowed)}; got {v!r}"
            )
        return v_norm

    pedigree_confidence_threshold: float = Field(
        default=0.7,
        env="PEDIGREE_CONFIDENCE_THRESHOLD",
        description="Minimum confidence (0.0-1.0) to classify an image as a pedigree",
    )

    # Tier 3.5: Browser-based HTML fallback for free-after-embargo papers.
    # Drives Playwright to fetch fully-rendered HTML from publishers that block
    # requests-based access; reuses the existing SupplementScraper to parse the
    # resulting DOM. ON by default — every extraction run automatically attempts
    # browser scraping for paywalled papers. Disable with
    # --no-browser-html-fallback or ENABLE_BROWSER_HTML_FALLBACK=false.
    enable_browser_html_fallback: bool = Field(
        default=True,
        env="ENABLE_BROWSER_HTML_FALLBACK",
        description="Enable Tier 3.5 browser-based HTML harvesting (Playwright)",
    )
    browser_html_publisher_allowlist: Union[str, List[str]] = Field(
        default="aha,oxford,wiley,elsevier_open,generic",
        env="BROWSER_HTML_PUBLISHER_ALLOWLIST",
        description="Comma-separated strategy NAMEs to enable for Tier 3.5",
    )
    browser_html_min_embargo_months: Optional[int] = Field(
        default=None,
        env="BROWSER_HTML_MIN_EMBARGO_MONTHS",
        description=(
            "Optional global override for minimum embargo months. If unset"
            " (default), each strategy uses its own EMBARGO_MONTHS. If set,"
            " takes max(strategy_embargo, this value)."
        ),
    )
    browser_html_headless: bool = Field(
        default=True,
        env="BROWSER_HTML_HEADLESS",
        description="Run Tier 3.5 browser headless",
    )
    browser_html_use_profile: bool = Field(
        default=False,
        env="BROWSER_HTML_USE_PROFILE",
        description=(
            "Run Tier 3.5 in a persistent Chrome profile so institutional "
            "SSO/OpenAthens cookies can be reused for paywalled publisher pages."
        ),
    )
    browser_html_profile_path: Optional[str] = Field(
        default=None,
        env="BROWSER_HTML_PROFILE_PATH",
        description=(
            "Persistent Chrome user-data directory for Tier 3.5. Use a dedicated "
            "GVF profile, not the daily Chrome profile that may be locked."
        ),
    )
    browser_html_channel: str = Field(
        default="chrome",
        env="BROWSER_HTML_CHANNEL",
        description=(
            "Browser channel for persistent profile mode. 'chrome' uses the "
            "installed Chrome with normal profile/cookie behavior."
        ),
    )
    browser_html_slow_mo_ms: int = Field(
        default=0,
        env="BROWSER_HTML_SLOW_MO_MS",
        description="Optional Playwright slow_mo delay for Tier 3.5 browser actions.",
    )
    browser_html_max_per_run: int = Field(
        default=0,
        env="BROWSER_HTML_MAX_PER_RUN",
        description=(
            "Hard cap on Tier 3.5 attempts per harvest run. 0 means unlimited; "
            "set a positive value for smoke tests or quota-limited runs."
        ),
    )
    browser_html_per_paper_timeout_s: int = Field(
        default=90,
        env="BROWSER_HTML_PER_PAPER_TIMEOUT_S",
        description="Per-paper timeout for Tier 3.5 attempts (seconds)",
    )

    model_config = SettingsConfigDict(
        env_file=_ENV_PATH,
        env_file_encoding="utf-8",
        extra="ignore",
        populate_by_name=True,
    )

    @field_validator("tier3_models", mode="after")
    @classmethod
    def split_str(cls, v):
        if isinstance(v, str):
            # Handle empty string by returning default
            if not v.strip():
                return ["gpt-4o-mini", "gpt-4o"]
            return [item.strip() for item in v.split(",") if item.strip()]
        return v

    @field_validator("tier3_adjudicator_models", mode="after")
    @classmethod
    def split_adjudicator_models(cls, v):
        if isinstance(v, str):
            return [item.strip() for item in v.split(",") if item.strip()]
        return v

    @field_validator("early_debate_models", mode="after")
    @classmethod
    def split_early_debate_models(cls, v):
        if isinstance(v, str):
            return [item.strip() for item in v.split(",") if item.strip()]
        return v

    @field_validator("browser_html_publisher_allowlist", mode="after")
    @classmethod
    def split_browser_html_allowlist(cls, v):
        if isinstance(v, str):
            if not v.strip():
                return ["aha", "oxford", "wiley", "elsevier_open", "generic"]
            return [item.strip().lower() for item in v.split(",") if item.strip()]
        if isinstance(v, list):
            return [str(item).strip().lower() for item in v if str(item).strip()]
        return v

    _PLACEHOLDER_VALUES: ClassVar[set] = {
        "your_api_key",
        "your_email@example.com",
        "changeme",
        "change_me",
        "your-openai-api-key",
        "your-anthropic-api-key",
    }

    @staticmethod
    def _is_placeholder(value: str) -> bool:
        normalized = value.strip().lower()
        return normalized in Settings._PLACEHOLDER_VALUES or "example.com" in normalized

    @model_validator(mode="after")
    def validate_required_settings(self):
        # NCBI email is always required (PubMed compliance).
        # An LLM provider is required, but any of OPENAI_API_KEY,
        # AZURE_AI_API_KEY, or ANTHROPIC_API_KEY satisfies it — the pipeline
        # uses whichever the configured tier model strings point at and
        # LiteLLM resolves credentials per-call.
        always_required = {"ncbi_email": "NCBI_EMAIL"}
        llm_provider_required = {
            "openai_api_key": "OPENAI_API_KEY",
            "azure_ai_api_key": "AZURE_AI_API_KEY",
            "anthropic_api_key": "ANTHROPIC_API_KEY",
        }

        missing = []
        placeholders = []

        for field_name, env_var in always_required.items():
            value = getattr(self, field_name)
            if value is None or not str(value).strip():
                missing.append(env_var)
                continue
            if self._is_placeholder(str(value)):
                placeholders.append(env_var)

        # At least one LLM provider key must be set.
        llm_provider_present = False
        llm_placeholders = []
        for field_name, env_var in llm_provider_required.items():
            value = getattr(self, field_name)
            if value is None or not str(value).strip():
                continue
            if self._is_placeholder(str(value)):
                llm_placeholders.append(env_var)
                continue
            llm_provider_present = True
        if not llm_provider_present:
            if llm_placeholders:
                placeholders.extend(llm_placeholders)
            else:
                missing.append(" or ".join(llm_provider_required.values()))

        error_parts = []
        if missing:
            error_parts.append(
                "Missing required settings: "
                + ", ".join(missing)
                + ". Set these environment variables or update your .env file."
            )
        if placeholders:
            error_parts.append(
                "The following settings still use placeholder values and must be updated: "
                + ", ".join(placeholders)
                + "."
            )

        if error_parts:
            raise ValueError(" ".join(error_parts))

        return self

    @model_validator(mode="after")
    def validate_source_configuration(self):
        effective_pubmed = self.use_pubmed
        effective_europepmc = self.use_europepmc

        if self.pubmind_only:
            effective_pubmed = False
            effective_europepmc = False

        sources_enabled = [
            ("PubMind", self.use_pubmind),
            ("PubMed", effective_pubmed),
            ("EuropePMC", effective_europepmc),
        ]
        active_sources = [name for name, enabled in sources_enabled if enabled]

        if not active_sources:
            raise ValueError(
                "No literature sources are enabled. By default PUBMIND_ONLY=false enables PubMind, PubMed, and EuropePMC. "
                "Set PUBMIND_ONLY=true to restrict to PubMind only or toggle USE_PUBMED/USE_EUROPEPMC to disable specific sources."
            )

        logger.info(
            "Literature sourcing enabled for %s. Set PUBMIND_ONLY=true or disable USE_PUBMED/USE_EUROPEPMC to limit discovery.",
            ", ".join(active_sources),
        )

        return self

    # ------------------------------------------------------------------
    # Provider-aware model resolution
    # ------------------------------------------------------------------
    # The pipeline asks "which model should this tier use?" via these
    # helpers. They prefer (1) an explicit TIER*_MODEL env var, then
    # (2) the Anthropic default when MODEL_PROVIDER=anthropic, then
    # (3) the original tier field default. This keeps every tier
    # independently configurable while letting --model-provider switch
    # them in bulk.

    def _is_anthropic(self) -> bool:
        return (self.model_provider or "").strip().lower() == "anthropic"

    def _tier_env_set(self, field_name: str) -> bool:
        env_var = _TIER_ENV_VARS.get(field_name)
        if not env_var:
            return False
        value = os.getenv(env_var)
        if value is not None:
            return value.strip() != ""
        if field_name in self.model_fields_set:
            configured = getattr(self, field_name, None)
            if isinstance(configured, str):
                return configured.strip() != ""
            if isinstance(configured, list):
                return bool(configured)
            return configured is not None
        return False

    def get_tier2_model(self) -> str:
        if self._is_anthropic() and not self._tier_env_set("tier2_model"):
            return self.anthropic_tier2_model
        return self.tier2_model

    def get_tier3_models(self) -> List[str]:
        if self._is_anthropic() and not self._tier_env_set("tier3_models"):
            value = self.anthropic_tier3_models
            return [m.strip() for m in value.split(",") if m.strip()]
        models = self.tier3_models
        if isinstance(models, str):
            return [m.strip() for m in models.split(",") if m.strip()]
        return list(models)

    def get_tier3_adjudicator_models(self) -> List[str]:
        models = self.tier3_adjudicator_models
        if isinstance(models, str):
            return [m.strip() for m in models.split(",") if m.strip()]
        return list(models)

    def get_table_router_model(self) -> str:
        if self._is_anthropic() and not self._tier_env_set("table_router_model"):
            return self.anthropic_table_router_model
        return self.table_router_model

    def get_vision_model(self) -> str:
        if self._is_anthropic() and not self._tier_env_set("vision_model"):
            return self.anthropic_vision_model
        return self.vision_model

    @staticmethod
    def _azure_model_string(deployment: str | None) -> str | None:
        if deployment is None:
            return None
        value = deployment.strip()
        if not value:
            return None
        if value.startswith("azure_ai/"):
            return value
        return f"azure_ai/{value}"

    def get_experimental_azure_models(self) -> List[str]:
        """Return opt-in Azure deployments for side-by-side model probes."""
        models = [
            self._azure_model_string(self.azure_deployment_gpt54),
            self._azure_model_string(self.azure_deployment_deepseek),
        ]
        return [model for model in models if model]

    def get_early_debate_models(self) -> List[str]:
        """Return models used for the cheap reciprocal debate tier."""
        if isinstance(self.early_debate_models, str):
            configured = [
                item.strip()
                for item in self.early_debate_models.split(",")
                if item.strip()
            ]
        else:
            configured = list(self.early_debate_models)
        if configured:
            return configured
        defaults = [
            self._azure_model_string(self.azure_deployment_gpt54 or "gpt-5.4"),
            self._azure_model_string(
                self.azure_deployment_deepseek or "DeepSeek-V4-Pro"
            ),
        ]
        return [model for model in defaults if model]

    # ------------------------------------------------------------------
    # Provider-aware rate-limit / concurrency
    # ------------------------------------------------------------------
    # Anthropic's tier-4 RPM budget (200+ for Sonnet/Haiku) is wide enough
    # that the previous global default of 50 RPM bottlenecked the pipeline.
    # Azure deployments share a much tighter quota and need to stay at 3
    # workers / 50 RPM to avoid 429 storms during filter sweeps. Callers
    # should prefer these helpers over reading max_workers directly.

    def get_requests_per_minute(self) -> int:
        """Resolve the active LLM RPM cap.

        Honors LLM_REQUESTS_PER_MINUTE as an explicit override (set in the
        process env, not pydantic — we read it directly so test monkeypatches
        and one-shot CLI invocations work without rebuilding Settings). When
        unset, picks anthropic_rpm vs azure_rpm based on MODEL_PROVIDER.
        """
        raw = os.getenv("LLM_REQUESTS_PER_MINUTE", "").strip()
        if raw:
            try:
                value = int(raw)
                if value > 0:
                    return value
            except ValueError:
                pass  # fall through to provider default
        return self.anthropic_rpm if self._is_anthropic() else self.azure_rpm

    def get_max_workers(self) -> int:
        """Resolve the active parallel-worker cap.

        MAX_WORKERS env var wins (explicit user override). Otherwise picks
        anthropic_max_workers vs azure_max_workers based on MODEL_PROVIDER.
        """
        if os.getenv("MAX_WORKERS", "").strip():
            return self.max_workers
        return (
            self.anthropic_max_workers
            if self._is_anthropic()
            else self.azure_max_workers
        )


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return a cached Settings instance so configuration is loaded once."""
    return Settings()


def reset_settings_cache() -> None:
    """Clear the cached Settings so the next get_settings() re-reads env vars.

    Used by the CLI when --model-provider is supplied: the flag mutates
    MODEL_PROVIDER in the process env and we need the next settings load to
    pick that up rather than returning a stale cached instance.
    """
    get_settings.cache_clear()
