"""Shared pytest fixtures and configuration for GeneVariantFetcher tests.

This module provides:
- Path constants for test fixtures
- Common fixtures used across unit and integration tests
- Test configuration helpers
"""

import os
from pathlib import Path

import pytest
from dotenv import load_dotenv

# Load environment variables for tests
load_dotenv()


def _has_placeholder_ncbi_email() -> bool:
    value = os.getenv("NCBI_EMAIL", "").strip().lower()
    return (
        not value
        or value in {"your_email@example.com", "you@example.com"}
        or "example.com" in value
    )


# =============================================================================
# PATH CONSTANTS
# =============================================================================

# Base paths
TESTS_DIR = Path(__file__).parent
FIXTURES_DIR = TESTS_DIR / "fixtures"

# Fixture subdirectories
PMIDS_DIR = FIXTURES_DIR / "pmids"
HARVEST_DATA_DIR = FIXTURES_DIR / "harvest_data"

# Specific test data directories (paths only — tests skip themselves
# at runtime when these directories aren't present locally)
KCNH2_19716085_DIR = HARVEST_DATA_DIR / "test_kcnh2_19716085"
TEST_HARVEST_DIR = HARVEST_DATA_DIR / "test_harvest"
TEST_PMC_HARVEST_DIR = HARVEST_DATA_DIR / "test_pmc_harvest"

# Test output directory (created during test runs)
TEST_RUN_OUTPUT_DIR = FIXTURES_DIR / "test_run_output"


# =============================================================================
# PYTEST CONFIGURATION
# =============================================================================


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "requires_api: marks tests that require API keys"
    )
    config.addinivalue_line(
        "markers", "requires_network: marks tests that require network access"
    )


def pytest_collection_modifyitems(config, items):
    """Auto-skip tests that require network access unless explicitly selected."""
    if config.getoption("-m", default=None) and "requires_network" in config.getoption(
        "-m"
    ):
        return  # User explicitly selected network tests
    skip_network = pytest.mark.skip(
        reason="requires network access (run with -m 'requires_network')"
    )
    for item in items:
        if "requires_network" in item.keywords:
            item.add_marker(skip_network)


# =============================================================================
# SHARED FIXTURES
# =============================================================================


@pytest.fixture
def fixtures_dir():
    """Return the path to the fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def pmids_dir():
    """Return the path to the pmids fixtures directory."""
    return PMIDS_DIR


@pytest.fixture
def test_pmids_file():
    """Return the path to the test PMIDs file."""
    return PMIDS_DIR / "test_pmids.txt"


@pytest.fixture
def temp_output_dir(tmp_path):
    """Provide a temporary output directory for tests."""
    output_dir = tmp_path / "test_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


# =============================================================================
# SETTINGS FIXTURES
# =============================================================================


@pytest.fixture
def clear_settings_cache():
    """Clear the settings cache before and after a test.

    Use this fixture when testing configuration to ensure
    environment changes take effect.
    """
    from config.settings import get_settings

    get_settings.cache_clear()
    yield
    get_settings.cache_clear()


@pytest.fixture(autouse=True)
def isolated_settings(monkeypatch):
    """Keep offline unit tests independent of a developer's local .env."""
    from config.settings import get_settings

    if _has_placeholder_ncbi_email():
        monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    # A developer .env may set GVF_COOKIE_FILE (loaded above via load_dotenv), which
    # would route every load_chrome_cookies() call through the file path and break
    # tests exercising the Chrome cookie path. Tests that want the file path set it
    # explicitly (see test_cookie_file_loader.py).
    monkeypatch.delenv("GVF_COOKIE_FILE", raising=False)
    get_settings.cache_clear()
    yield
    get_settings.cache_clear()


# =============================================================================
# SKIP CONDITIONS
# =============================================================================


def has_llm_provider_key():
    """Check if at least one supported LLM provider key is available."""
    return any(
        bool(os.getenv(env_var))
        for env_var in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "AZURE_AI_API_KEY")
    )


def has_openai_key():
    """Backward-compatible alias for older tests."""
    return has_llm_provider_key()


def has_ncbi_email():
    """Check if NCBI email is configured."""
    return bool(os.getenv("NCBI_EMAIL"))


requires_llm = pytest.mark.skipif(
    not has_llm_provider_key(),
    reason="no LLM provider key set (ANTHROPIC_API_KEY, OPENAI_API_KEY, or AZURE_AI_API_KEY)",
)
requires_openai = requires_llm

requires_ncbi = pytest.mark.skipif(not has_ncbi_email(), reason="NCBI_EMAIL not set")
