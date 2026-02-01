"""Validated configuration for the Gene Variant Fetcher pipeline."""

import logging
from functools import lru_cache
from pathlib import Path
from typing import ClassVar, List, Optional, Union

from pydantic import Field, model_validator, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_ENV_PATH = Path(__file__).resolve().parent.parent / ".env"
logger = logging.getLogger(__name__)


class Settings(BaseSettings):
    """Centralized application settings loaded from environment variables."""

    # API Keys - all loaded from .env file
    openai_api_key: Optional[str] = Field(default=None, env="OPENAI_API_KEY")
    anthropic_api_key: Optional[str] = Field(default=None, env="ANTHROPIC_API_KEY")
    gemini_api_key: Optional[str] = Field(default=None, env="GEMINI_API_KEY")
    ncbi_email: Optional[str] = Field(default=None, env="NCBI_EMAIL")
    ncbi_api_key: Optional[str] = Field(default=None, env="NCBI_API_KEY")
    elsevier_api_key: Optional[str] = Field(default=None, env="ELSEVIER_API_KEY")
    wiley_api_key: Optional[str] = Field(default=None, env="WILEY_API_KEY")

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
        default=2,
        env="TIER1_MIN_KEYWORDS",
        description="Minimum keyword matches for Tier 1 to pass",
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
        default=150,
        env="TIER2_MAX_TOKENS",
        description="Max tokens for Tier 2 LLM response",
    )
    tier2_confidence_threshold: float = Field(
        default=0.5,
        env="TIER2_CONFIDENCE_THRESHOLD",
        description="Minimum confidence to pass Tier 2",
    )

    # Tier 3 Configuration
    tier3_temperature: float = Field(
        default=0.0, env="TIER3_TEMPERATURE", description="Temperature for Tier 3 LLM"
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
        default=100,
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
        default="gpt-4o",
        env="VISION_MODEL",
        description="Vision-capable model for pedigree analysis (gpt-4o, gemini-1.5-pro, claude-3-5-sonnet)",
    )
    pedigree_confidence_threshold: float = Field(
        default=0.7,
        env="PEDIGREE_CONFIDENCE_THRESHOLD",
        description="Minimum confidence (0.0-1.0) to classify an image as a pedigree",
    )

    model_config = SettingsConfigDict(
        env_file=_ENV_PATH,
        env_file_encoding="utf-8",
        extra="ignore",
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
        required = {
            "openai_api_key": "OPENAI_API_KEY",
            "ncbi_email": "NCBI_EMAIL",
        }

        missing = []
        placeholders = []

        for field_name, env_var in required.items():
            value = getattr(self, field_name)
            if value is None or not str(value).strip():
                missing.append(env_var)
                continue
            if self._is_placeholder(str(value)):
                placeholders.append(env_var)

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


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return a cached Settings instance so configuration is loaded once."""
    return Settings()
