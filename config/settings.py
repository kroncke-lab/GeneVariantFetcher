"""
Configuration management for the Gene Variant Fetcher pipeline.
Loads settings from environment variables and .env file using singleton pattern.
"""

import os
from functools import lru_cache
from pathlib import Path
from dotenv import load_dotenv

_ENV_PATH = Path(__file__).resolve().parent.parent / ".env"
load_dotenv(_ENV_PATH)


class Settings:
    """Centralized application settings loaded from environment variables."""

    # API Keys
    openai_api_key: str
    anthropic_api_key: str

    # PubMed/NCBI Configuration
    ncbi_email: str
    ncbi_api_key: str

    # Pipeline Configuration
    intern_model: str
    extractor_model: str
    rate_limit_per_minute: int

    def __init__(self):
        # API Keys
        self.openai_api_key = os.getenv("OPENAI_API_KEY", "")
        self.anthropic_api_key = os.getenv("ANTHROPIC_API_KEY", "")

        # PubMed/NCBI Configuration
        self.ncbi_email = os.getenv("NCBI_EMAIL", "your_email@example.com")
        self.ncbi_api_key = os.getenv("NCBI_API_KEY", "")

        # Pipeline Configuration
        self.intern_model = os.getenv("INTERN_MODEL", "gpt-4o-mini")
        self.extractor_model = os.getenv("EXTRACTOR_MODEL", "gpt-4o")
        self.rate_limit_per_minute = int(os.getenv("RATE_LIMIT_PER_MINUTE", "60"))


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return a cached Settings instance so configuration is loaded once."""
    return Settings()