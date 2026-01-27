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

# =============================================================================
# PATH CONSTANTS
# =============================================================================

# Base paths
TESTS_DIR = Path(__file__).parent
FIXTURES_DIR = TESTS_DIR / "fixtures"

# Fixture subdirectories
PMIDS_DIR = FIXTURES_DIR / "pmids"
HARVEST_DATA_DIR = FIXTURES_DIR / "harvest_data"
COMPARISON_DATA_DIR = FIXTURES_DIR / "comparison_data"

# Specific test data directories
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
def harvest_data_dir():
    """Return the path to the harvest data fixtures directory."""
    return HARVEST_DATA_DIR


@pytest.fixture
def comparison_data_dir():
    """Return the path to the comparison data fixtures directory."""
    return COMPARISON_DATA_DIR


@pytest.fixture
def test_pmids_file():
    """Return the path to the test PMIDs file."""
    return PMIDS_DIR / "test_pmids.txt"


@pytest.fixture
def kcnh2_19716085_dir():
    """Return the path to the KCNH2 19716085 test data directory."""
    return KCNH2_19716085_DIR


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


# =============================================================================
# SKIP CONDITIONS
# =============================================================================


def has_openai_key():
    """Check if OpenAI API key is available."""
    return bool(os.getenv("OPENAI_API_KEY"))


def has_ncbi_email():
    """Check if NCBI email is configured."""
    return bool(os.getenv("NCBI_EMAIL"))


requires_openai = pytest.mark.skipif(
    not has_openai_key(), reason="OPENAI_API_KEY not set"
)

requires_ncbi = pytest.mark.skipif(not has_ncbi_email(), reason="NCBI_EMAIL not set")
