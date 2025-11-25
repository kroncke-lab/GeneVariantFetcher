"""Validated configuration for the Gene Variant Fetcher pipeline."""

from functools import lru_cache
from pathlib import Path
from typing import ClassVar

from pydantic import Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_ENV_PATH = Path(__file__).resolve().parent.parent / ".env"


class Settings(BaseSettings):
    """Centralized application settings loaded from environment variables."""

    openai_api_key: str | None = Field(default=None, env="OPENAI_API_KEY")
    anthropic_api_key: str | None = Field(default=None, env="ANTHROPIC_API_KEY")
    ncbi_email: str | None = Field(default=None, env="NCBI_EMAIL")
    ncbi_api_key: str | None = Field(default=None, env="NCBI_API_KEY")

    intern_model: str = Field(default="gpt-4o-mini", env="INTERN_MODEL")
    extractor_model: str = Field(default="gpt-4o", env="EXTRACTOR_MODEL")
    rate_limit_per_minute: int = Field(default=60, env="RATE_LIMIT_PER_MINUTE")

    model_config = SettingsConfigDict(
        env_file=_ENV_PATH,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    _PLACEHOLDER_VALUES: ClassVar[set[str]] = {
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
