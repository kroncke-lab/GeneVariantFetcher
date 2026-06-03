"""
Comprehensive Integration Test Suite for GeneVariantFetcher

Tests all new integrations against gold standard KCNH2 PMIDs with real API calls.

Gold standard PMIDs (verified via API):
    - 24667783: PMC4266740, DOI 10.1038/ejhg.2014.54 (Nature/EJHG)
    - 19841300: PMC3025752, DOI 10.1161/CIRCULATIONAHA.109.863076 (Circulation)
    - 20173333: No PMC, DOI 10.1159/000285512 (Karger)
    - 30036649: PMC6198685, DOI 10.1016/j.ipej.2018.07.007 (Elsevier/IPEJ)

Usage:
    pytest tests/test_all_integrations.py -v
    pytest tests/test_all_integrations.py -v -m "not slow"
"""

import logging
import sys
import time
import traceback
from pathlib import Path

import pytest

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Gold standard KCNH2 test PMIDs (values verified against live APIs)
GOLD_PMIDS = {
    "24667783": {"has_pmc": True, "pmcid": "PMC4266740"},
    "19841300": {"has_pmc": True, "pmcid": "PMC3025752"},
    "20173333": {"has_pmc": False, "pmcid": None},
    "30036649": {"has_pmc": True, "pmcid": "PMC6198685"},
}

# Known DOIs (verified against live APIs)
GOLD_DOIS = {
    "24667783": "10.1038/ejhg.2014.54",
    "19841300": "10.1161/CIRCULATIONAHA.109.863076",
    "20173333": "10.1159/000285512",
    "30036649": "10.1016/j.ipej.2018.07.007",
}


# ============================================================================
# Helpers
# ============================================================================


class _ResultTracker:
    """Track individual test results for summary."""

    def __init__(self):
        self.results = []

    def record(self, name: str, passed: bool, detail: str = ""):
        self.results.append({"name": name, "passed": passed, "detail": detail})

    def summary(self) -> str:
        passed = sum(1 for r in self.results if r["passed"])
        total = len(self.results)
        lines = [
            f"\n{'=' * 70}",
            f"INTEGRATION TEST SUMMARY: {passed}/{total} tests passed",
            f"{'=' * 70}",
        ]
        for r in self.results:
            status = "PASS" if r["passed"] else "FAIL"
            detail = f" - {r['detail']}" if r["detail"] else ""
            lines.append(f"  [{status}] {r['name']}{detail}")
        lines.append(f"{'=' * 70}\n")
        return "\n".join(lines)


# Module-level tracker for standalone summary
_tracker = _ResultTracker()


# ============================================================================
# Europe PMC Tests (via the production query_europepmc entry point)
# ============================================================================


@pytest.mark.requires_network
class TestEuropePMC:
    """Smoke-test the Europe PMC capability the pipeline actually uses.

    The bespoke per-paper Europe PMC client/harvester classes were retired
    (unused by the pipeline, which discovers PMIDs via
    utils.pubmed_utils.query_europepmc). These tests cover that production entry
    point: gene symbol -> set of PMIDs.
    """

    def test_query_returns_pmids(self):
        """A known gene returns a non-empty set of PMID strings."""
        from utils.pubmed_utils import query_europepmc

        pmids = query_europepmc("KCNH2", max_results=20)

        assert isinstance(pmids, set), "query_europepmc should return a set"
        assert len(pmids) > 0, "Should find KCNH2 papers in Europe PMC"
        assert all(p.isdigit() for p in pmids), "Results should be PMID strings"

        _tracker.record(
            "EuropePMC: query_europepmc('KCNH2')",
            True,
            f"found {len(pmids)} PMIDs",
        )

    def test_query_respects_max_results(self):
        """max_results caps the returned set."""
        from utils.pubmed_utils import query_europepmc

        pmids = query_europepmc("KCNH2", max_results=5)

        assert len(pmids) <= 5, "Should respect max_results"

        _tracker.record("EuropePMC: query_europepmc respects max_results", True)

    def test_query_nonsense_gene_is_graceful(self):
        """A bogus gene token returns an empty set, not an error."""
        from utils.pubmed_utils import query_europepmc

        pmids = query_europepmc("ZZZNOTAGENEZZZ", max_results=5)

        assert isinstance(pmids, set), "Should return a set even with no hits"

        _tracker.record("EuropePMC: query_europepmc bogus-gene handling", True)


# ============================================================================
# PMC API Client Tests
# ============================================================================


@pytest.mark.requires_network
class TestPMCAPIClient:
    """Test harvesting.pmc_api.PMCAPIClient."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from harvesting.pmc_api import PMCAPIClient

        self.client = PMCAPIClient()

    def test_pmid_to_pmcid(self):
        """Convert a known PMID to PMCID."""
        pmcid = self.client.pmid_to_pmcid("24667783")

        assert pmcid is not None, "Should find PMCID for 24667783"
        assert pmcid == "PMC4266740", f"Expected PMC4266740, got {pmcid}"

        _tracker.record("PMC API: pmid_to_pmcid(24667783)", True, f"got {pmcid}")

    def test_pmid_to_pmcid_multiple(self):
        """Verify PMID-to-PMCID conversion for multiple gold standard papers."""
        for pmid, info in GOLD_PMIDS.items():
            if info["pmcid"]:
                pmcid = self.client.pmid_to_pmcid(pmid)
                assert (
                    pmcid == info["pmcid"]
                ), f"PMID {pmid}: expected {info['pmcid']}, got {pmcid}"
                time.sleep(0.4)

        _tracker.record("PMC API: pmid_to_pmcid for gold PMIDs", True)

    def test_pmid_without_pmc(self):
        """PMID 20173333 has no PMC entry - should return None."""
        pmcid = self.client.pmid_to_pmcid("20173333")
        assert pmcid is None, f"Expected None for PMID without PMC, got {pmcid}"

        _tracker.record("PMC API: PMID without PMC returns None", True)

    def test_get_doi_from_pmid(self):
        """Fetch DOI for a known PMID."""
        doi = self.client.get_doi_from_pmid("24667783")

        assert doi is not None, "Should find DOI for 24667783"
        assert (
            doi == GOLD_DOIS["24667783"]
        ), f"Expected {GOLD_DOIS['24667783']}, got {doi}"

        _tracker.record("PMC API: get_doi_from_pmid(24667783)", True, f"doi={doi}")

    def test_get_doi_all_gold(self):
        """Verify DOI retrieval for all gold standard PMIDs."""
        for pmid, expected_doi in GOLD_DOIS.items():
            doi = self.client.get_doi_from_pmid(pmid)
            assert doi is not None, f"No DOI for PMID {pmid}"
            # Compare case-insensitively since DOIs can vary in case
            assert (
                doi.lower() == expected_doi.lower()
            ), f"PMID {pmid}: expected {expected_doi}, got {doi}"
            time.sleep(0.4)

        _tracker.record("PMC API: DOI retrieval for all gold PMIDs", True)

    def test_nonexistent_pmid(self):
        """Non-existent PMID should return None."""
        pmcid = self.client.pmid_to_pmcid("99999999999")
        assert pmcid is None, "Should return None for bogus PMID"

        _tracker.record("PMC API: nonexistent PMID handling", True)


# ============================================================================
# Unpaywall API Tests
# ============================================================================


@pytest.mark.requires_network
class TestUnpaywallAPI:
    """Test harvesting.unpaywall_api.UnpaywallClient."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from harvesting.unpaywall_api import UnpaywallClient

        self.client = UnpaywallClient(email="user@example.org")

    def test_find_open_access(self):
        """Look up a known paper by DOI."""
        doi = GOLD_DOIS["24667783"]
        result, error = self.client.find_open_access(doi)

        assert error is None, f"Unpaywall error: {error}"
        assert result is not None, "Should return result for valid DOI"
        assert result["doi"] is not None, "Result should have DOI"
        assert "is_oa" in result, "Result should have is_oa field"
        assert "oa_status" in result, "Result should have oa_status field"

        _tracker.record(
            "Unpaywall: find_open_access(24667783 DOI)",
            True,
            f"is_oa={result['is_oa']}, status={result['oa_status']}",
        )

    def test_find_open_access_all_gold(self):
        """Check Unpaywall for all gold standard DOIs."""
        oa_count = 0
        for pmid, doi in GOLD_DOIS.items():
            result, error = self.client.find_open_access(doi)
            assert (
                result is not None or error is not None
            ), f"Should return either result or error for {doi}"
            if result and result.get("is_oa"):
                oa_count += 1
            time.sleep(0.2)

        _tracker.record(
            "Unpaywall: OA check for all gold DOIs", True, f"{oa_count}/4 are OA"
        )

    def test_doi_with_url_prefix(self):
        """Verify DOI cleaning handles https://doi.org/ prefix."""
        doi = f"https://doi.org/{GOLD_DOIS['24667783']}"
        result, error = self.client.find_open_access(doi)

        assert error is None, f"Should handle doi.org prefix: {error}"
        assert result is not None

        _tracker.record("Unpaywall: DOI with https://doi.org/ prefix", True)

    def test_invalid_doi(self):
        """Invalid DOI should return error, not crash."""
        result, error = self.client.find_open_access("10.9999/nonexistent.fake.doi")

        # Should get None result or an error message, not an exception
        assert (
            result is None or error is not None or (result and not result.get("is_oa"))
        ), "Invalid DOI should not return a valid OA result"

        _tracker.record("Unpaywall: invalid DOI handling", True)

    def test_empty_doi(self):
        """Empty DOI should return error gracefully."""
        result, error = self.client.find_open_access("")

        assert result is None, "Empty DOI should return None"
        assert error is not None, "Empty DOI should have error message"

        _tracker.record("Unpaywall: empty DOI handling", True)


# ============================================================================
# Supplement Scraper Tests
# ============================================================================


@pytest.mark.requires_network
class TestSupplementScraper:
    """Test harvesting.supplement_scraper.SupplementScraper."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from harvesting.supplement_scraper import SupplementScraper

        self.scraper = SupplementScraper()

    def test_import_supplement_scraper(self):
        """Verify SupplementScraper class imports and instantiates."""
        from harvesting.supplement_scraper import SupplementScraper

        scraper = SupplementScraper()
        assert scraper is not None

        _tracker.record("Supplement Scraper: class import", True)

    def test_has_publisher_methods(self):
        """Verify all publisher scrape methods exist."""
        publishers = [
            "scrape_nature_supplements",
            "scrape_elsevier_supplements",
            "scrape_springer_supplements",
            "scrape_oxford_supplements",
            "scrape_wiley_supplements",
            "scrape_karger_supplements",
            "scrape_generic_supplements",
        ]
        missing = [p for p in publishers if not hasattr(self.scraper, p)]
        assert not missing, f"Missing publisher methods: {missing}"

        _tracker.record(
            "Supplement Scraper: all publisher methods exist",
            True,
            f"{len(publishers)} methods verified",
        )

    def test_extract_fulltext(self):
        """Verify extract_fulltext method exists and is callable."""
        assert hasattr(
            self.scraper, "extract_fulltext"
        ), "Should have extract_fulltext method"
        assert callable(self.scraper.extract_fulltext)

        _tracker.record("Supplement Scraper: extract_fulltext method exists", True)


# ============================================================================
# PMID Status Tracker Tests
# ============================================================================


class TestPMIDStatus:
    """Test utils.pmid_status module (local, no network needed)."""

    def test_import(self):
        """Verify pmid_status module imports."""
        from utils.pmid_status import (
            get_pmid_status,
            get_failed_pmids,
            get_stats_summary,
        )

        assert callable(get_pmid_status)
        assert callable(get_failed_pmids)
        assert callable(get_stats_summary)

        _tracker.record("PMID Status: module import", True)

    def test_get_pmid_status_missing(self):
        """Querying a non-existent status file should return None."""
        from utils.pmid_status import get_pmid_status

        result = get_pmid_status("/nonexistent/path", "99999")
        assert result is None, "Should return None for missing status"

        _tracker.record("PMID Status: missing status returns None", True)


# ============================================================================
# Cross-Module Integration Tests
# ============================================================================


@pytest.mark.requires_network
@pytest.mark.slow
class TestCrossModuleIntegration:
    """End-to-end tests combining multiple modules."""

    def test_pmid_to_metadata_to_unpaywall(self):
        """Full pipeline: PMID -> DOI (via PMC API) -> Unpaywall OA check."""
        from harvesting.pmc_api import PMCAPIClient
        from harvesting.unpaywall_api import UnpaywallClient

        pmc = PMCAPIClient()
        unpaywall = UnpaywallClient(email="user@example.org")

        pmid = "24667783"

        # Step 1: Get DOI from PMID
        doi = pmc.get_doi_from_pmid(pmid)
        assert doi is not None, f"Failed to get DOI for PMID {pmid}"
        time.sleep(0.4)

        # Step 2: Check Unpaywall for OA
        result, error = unpaywall.find_open_access(doi)
        assert result is not None, f"Unpaywall failed for DOI {doi}: {error}"

        _tracker.record(
            "Cross-module: PMID->DOI->Unpaywall pipeline",
            True,
            f"doi={doi}, is_oa={result.get('is_oa')}",
        )


# ============================================================================
# Summary fixture - prints results after all tests
# ============================================================================


@pytest.fixture(scope="session", autouse=True)
def print_summary(request):
    """Print test summary after all tests complete."""
    yield
    print(_tracker.summary())


# ============================================================================
# Standalone runner
# ============================================================================


def run_standalone():
    """Run tests outside pytest for quick smoke testing.

    Usage: python tests/test_all_integrations.py
    """
    tracker = _ResultTracker()

    print(f"\n{'=' * 70}")
    print("GeneVariantFetcher Integration Test Suite (standalone)")
    print(f"{'=' * 70}\n")

    # --- Europe PMC (query_europepmc) ---
    print("[1/4] Europe PMC query...")
    try:
        from utils.pubmed_utils import query_europepmc

        pmids = query_europepmc("KCNH2", max_results=10)
        tracker.record(
            "EuropePMC: query_europepmc('KCNH2')",
            isinstance(pmids, set) and len(pmids) > 0,
            f"{len(pmids)} PMIDs",
        )

    except Exception as e:
        tracker.record("EuropePMC: query_europepmc", False, str(e)[:80])
        traceback.print_exc()

    # --- PMC API ---
    print("[2/4] PMC API Client...")
    try:
        from harvesting.pmc_api import PMCAPIClient

        pmc = PMCAPIClient()

        pmcid = pmc.pmid_to_pmcid("24667783")
        tracker.record("PMC API: pmid_to_pmcid", pmcid == "PMC4266740", f"got {pmcid}")

        time.sleep(0.4)
        doi = pmc.get_doi_from_pmid("24667783")
        tracker.record(
            "PMC API: get_doi_from_pmid",
            doi is not None and doi == GOLD_DOIS["24667783"],
            f"doi={doi}",
        )

    except Exception as e:
        tracker.record("PMC API: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- Unpaywall ---
    print("[3/4] Unpaywall API...")
    try:
        from harvesting.unpaywall_api import UnpaywallClient

        uw = UnpaywallClient(email="user@example.org")

        result, error = uw.find_open_access(GOLD_DOIS["24667783"])
        tracker.record(
            "Unpaywall: find_open_access",
            result is not None and error is None,
            f"is_oa={result.get('is_oa')}" if result else f"error={error}",
        )

        time.sleep(0.2)
        result2, error2 = uw.find_open_access("")
        tracker.record(
            "Unpaywall: empty DOI handling",
            result2 is None and error2 is not None,
            "correctly rejected",
        )

    except Exception as e:
        tracker.record("Unpaywall: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- PMID Status ---
    print("[4/4] PMID Status Tracker...")
    try:
        from utils.pmid_status import get_pmid_status

        result = get_pmid_status("/nonexistent", "99999")
        tracker.record("PMID Status: missing status returns None", result is None)

    except Exception as e:
        tracker.record("PMID Status: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- Summary ---
    print(tracker.summary())
    passed = sum(1 for r in tracker.results if r["passed"])
    return passed == len(tracker.results)


if __name__ == "__main__":
    success = run_standalone()
    sys.exit(0 if success else 1)
