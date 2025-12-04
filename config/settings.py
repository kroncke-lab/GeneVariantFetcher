"""Validated configuration for the Gene Variant Fetcher pipeline."""

from functools import lru_cache
from pathlib import Path
from typing import ClassVar, List, Optional, Union

from pydantic import Field, model_validator, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_ENV_PATH = Path(__file__).resolve().parent.parent / ".env"


class Settings(BaseSettings):
    """Centralized application settings loaded from environment variables."""

    # API Keys
    openai_api_key: Optional[str] = Field(default=None, env="OPENAI_API_KEY")
    anthropic_api_key: Optional[str] = Field(default=None, env="ANTHROPIC_API_KEY")
    ncbi_email: Optional[str] = Field(default=None, env="NCBI_EMAIL")
    ncbi_api_key: Optional[str] = Field(default=None, env="NCBI_API_KEY")

    # Model Configuration
    tier1_model: Optional[str] = Field(default=None, env="TIER1_MODEL", description="Optional LLM for Tier 1 (if using LLM-based Tier 1)")
    tier2_model: str = Field(default="gpt-4o-mini", env="TIER2_MODEL", description="Model for Tier 2 classification (cheap)")
    tier3_models: Union[str, List[str]] = Field(
        default="gpt-4o-mini,gpt-4o",
        env="TIER3_MODELS",
        description="Comma-separated list of models for Tier 3 extraction (e.g., 'gpt-4o-mini,gpt-4o')"
    )

    # Legacy aliases for backward compatibility
    intern_model: str = Field(default="gpt-4o-mini", env="INTERN_MODEL")
    rate_limit_per_minute: int = Field(default=60, env="RATE_LIMIT_PER_MINUTE")

    # Tiered Classification Configuration
    enable_tier1: bool = Field(default=True, env="ENABLE_TIER1", description="Enable Tier 1 keyword/heuristic filtering")
    enable_tier2: bool = Field(default=True, env="ENABLE_TIER2", description="Enable Tier 2 LLM classification")
    enable_tier3: bool = Field(default=True, env="ENABLE_TIER3", description="Enable Tier 3 expert extraction")

    # Tier 1 Configuration
    tier1_min_keywords: int = Field(default=2, env="TIER1_MIN_KEYWORDS", description="Minimum keyword matches for Tier 1 to pass")
    tier1_use_llm: bool = Field(default=False, env="TIER1_USE_LLM", description="Use lightweight LLM for Tier 1 instead of keywords")

    # Tier 2 Configuration
    tier2_temperature: float = Field(default=0.1, env="TIER2_TEMPERATURE", description="Temperature for Tier 2 LLM")
    tier2_max_tokens: int = Field(default=150, env="TIER2_MAX_TOKENS", description="Max tokens for Tier 2 LLM response")
    tier2_confidence_threshold: float = Field(default=0.5, env="TIER2_CONFIDENCE_THRESHOLD", description="Minimum confidence to pass Tier 2")

    # Tier 3 Configuration
    tier3_temperature: float = Field(default=0.0, env="TIER3_TEMPERATURE", description="Temperature for Tier 3 LLM")
    tier3_max_tokens: int = Field(default=8000, env="TIER3_MAX_TOKENS", description="Max tokens for Tier 3 LLM response")
    tier3_threshold: int = Field(default=1, env="TIER3_THRESHOLD", description="Try next model if first finds fewer variants than this (0 = only use first model)")

    # Paper Sourcing Configuration
    use_pubmind: bool = Field(default=True, env="USE_PUBMIND", description="Use PubMind as primary literature source")
    use_pubmed: bool = Field(default=False, env="USE_PUBMED", description="Use PubMed API as additional source")
    use_europepmc: bool = Field(default=False, env="USE_EUROPEPMC", description="Use EuropePMC as additional source")
    pubmind_only: bool = Field(default=True, env="PUBMIND_ONLY", description="Use ONLY PubMind (ignore PubMed/EuropePMC)")
    max_papers_per_source: int = Field(default=100, env="MAX_PAPERS_PER_SOURCE", description="Max papers to fetch per source")

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
            "anthropic_api_key": "ANTHROPIC_API_KEY",
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


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return a cached Settings instance so configuration is loaded once."""
    return Settings()
