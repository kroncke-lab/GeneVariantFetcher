"""
Configuration management for the Tiered Biomedical Extraction Pipeline.
Loads settings from environment variables and .env file.
"""

import os
from typing import Optional
from dotenv import load_dotenv

# Load .env file if it exists
load_dotenv()


class Config:
    """Configuration settings loaded from environment variables."""

    # =============================================================================
    # API Keys
    # =============================================================================

    OPENAI_API_KEY: Optional[str] = os.getenv("OPENAI_API_KEY")
    ANTHROPIC_API_KEY: Optional[str] = os.getenv("ANTHROPIC_API_KEY")

    # =============================================================================
    # PubMed/NCBI Configuration
    # =============================================================================

    NCBI_EMAIL: str = os.getenv("NCBI_EMAIL", "your_email@example.com")
    NCBI_API_KEY: Optional[str] = os.getenv("NCBI_API_KEY")

    # =============================================================================
    # Pipeline Configuration
    # =============================================================================

    # Model selection
    INTERN_MODEL: str = os.getenv("INTERN_MODEL", "gpt-4o-mini")
    EXTRACTOR_MODEL: str = os.getenv("EXTRACTOR_MODEL", "gpt-4o")

    # Filter thresholds
    MIN_KEYWORD_MATCHES: int = int(os.getenv("MIN_KEYWORD_MATCHES", "2"))

    # Processing limits
    MAX_PAPERS_DEFAULT: Optional[int] = (
        int(os.getenv("MAX_PAPERS_DEFAULT"))
        if os.getenv("MAX_PAPERS_DEFAULT")
        else None
    )

    # =============================================================================
    # Logging Configuration
    # =============================================================================

    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")
    LOG_FILE: Optional[str] = os.getenv("LOG_FILE")

    # =============================================================================
    # Tier Toggles
    # =============================================================================

    ENABLE_TIER1: bool = os.getenv("ENABLE_TIER1", "true").lower() == "true"
    ENABLE_TIER2: bool = os.getenv("ENABLE_TIER2", "true").lower() == "true"
    ENABLE_TIER3: bool = os.getenv("ENABLE_TIER3", "true").lower() == "true"

    # =============================================================================
    # Advanced Settings
    # =============================================================================

    # Retry configuration
    MAX_RETRIES: int = int(os.getenv("MAX_RETRIES", "3"))
    RETRY_MIN_WAIT: int = int(os.getenv("RETRY_MIN_WAIT", "2"))
    RETRY_MAX_WAIT: int = int(os.getenv("RETRY_MAX_WAIT", "10"))

    # LLM parameters
    INTERN_TEMPERATURE: float = float(os.getenv("INTERN_TEMPERATURE", "0.1"))
    EXTRACTOR_TEMPERATURE: float = float(os.getenv("EXTRACTOR_TEMPERATURE", "0.0"))
    INTERN_MAX_TOKENS: int = int(os.getenv("INTERN_MAX_TOKENS", "150"))
    EXTRACTOR_MAX_TOKENS: int = int(os.getenv("EXTRACTOR_MAX_TOKENS", "4000"))

    @classmethod
    def validate(cls) -> bool:
        """
        Validate that required configuration is present.

        Returns:
            True if configuration is valid, False otherwise.
        """
        errors = []

        # Check for at least one LLM API key
        if not cls.OPENAI_API_KEY and not cls.ANTHROPIC_API_KEY:
            errors.append("No LLM API key found (OPENAI_API_KEY or ANTHROPIC_API_KEY)")

        # Check for NCBI email
        if cls.NCBI_EMAIL == "your_email@example.com":
            errors.append("NCBI_EMAIL not configured (required for PubMed API)")

        if errors:
            print("❌ Configuration Errors:")
            for error in errors:
                print(f"  - {error}")
            return False

        return True

    @classmethod
    def print_config(cls):
        """Print current configuration (hiding sensitive values)."""
        print("\n" + "="*80)
        print("PIPELINE CONFIGURATION")
        print("="*80)

        print("\n🔑 API Keys:")
        print(f"  OPENAI_API_KEY: {'✓ Set' if cls.OPENAI_API_KEY else '✗ Not set'}")
        print(f"  ANTHROPIC_API_KEY: {'✓ Set' if cls.ANTHROPIC_API_KEY else '✗ Not set'}")

        print("\n📧 NCBI Configuration:")
        print(f"  Email: {cls.NCBI_EMAIL}")
        print(f"  API Key: {'✓ Set' if cls.NCBI_API_KEY else '✗ Not set'}")

        print("\n🤖 Models:")
        print(f"  Intern (Tier 2): {cls.INTERN_MODEL}")
        print(f"  Extractor (Tier 3): {cls.EXTRACTOR_MODEL}")

        print("\n⚙️ Pipeline Settings:")
        print(f"  Tier 1 (Keyword): {'Enabled' if cls.ENABLE_TIER1 else 'Disabled'}")
        print(f"  Tier 2 (Intern): {'Enabled' if cls.ENABLE_TIER2 else 'Disabled'}")
        print(f"  Tier 3 (Extractor): {'Enabled' if cls.ENABLE_TIER3 else 'Disabled'}")
        print(f"  Min Keyword Matches: {cls.MIN_KEYWORD_MATCHES}")
        print(f"  Max Papers: {cls.MAX_PAPERS_DEFAULT or 'Unlimited'}")

        print("\n📊 Logging:")
        print(f"  Level: {cls.LOG_LEVEL}")
        print(f"  File: {cls.LOG_FILE or 'Console only'}")

        print("\n" + "="*80 + "\n")


# Create a singleton instance
config = Config()


if __name__ == "__main__":
    # When run directly, print configuration
    config.print_config()

    if not config.validate():
        print("\n⚠️  Please configure the required settings before running the pipeline.")
        print("Copy .env.example to .env and fill in your values.\n")
    else:
        print("✅ Configuration is valid!\n")
