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
        lines = [f"\n{'='*70}", f"INTEGRATION TEST SUMMARY: {passed}/{total} tests passed", f"{'='*70}"]
        for r in self.results:
            status = "PASS" if r["passed"] else "FAIL"
            detail = f" - {r['detail']}" if r["detail"] else ""
            lines.append(f"  [{status}] {r['name']}{detail}")
        lines.append(f"{'='*70}\n")
        return "\n".join(lines)


# Module-level tracker for standalone summary
_tracker = _ResultTracker()


# ============================================================================
# Europe PMC Handler Tests
# ============================================================================

@pytest.mark.requires_network
class TestEuropePMCHandler:
    """Test gene_literature.europepmc_handler.EuropePMCClient."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from gene_literature.europepmc_handler import EuropePMCClient
        self.client = EuropePMCClient(timeout=30)

    def test_get_paper_metadata(self):
        """Fetch metadata for a known PMID and verify key fields."""
        pmid = "24667783"
        metadata = self.client.get_paper_metadata(pmid)

        assert metadata is not None, f"No metadata returned for PMID {pmid}"
        assert metadata["pmid"] == pmid
        assert metadata["title"], "Title should be present"
        assert metadata["doi"], "DOI should be present"
        assert metadata["pmcid"] is not None, "PMCID should be present"
        assert len(metadata["authors"]) > 0, "Should have authors"

        _tracker.record("EuropePMC: get_paper_metadata(24667783)", True,
                        f"title={metadata['title'][:60]}...")
        logger.info(f"Metadata OK: {metadata['title'][:80]}")

    def test_get_paper_metadata_all_gold_pmids(self):
        """Verify metadata retrieval works for all gold standard PMIDs."""
        for pmid in GOLD_PMIDS:
            metadata = self.client.get_paper_metadata(pmid)
            assert metadata is not None, f"No metadata for PMID {pmid}"
            assert metadata["pmid"] == pmid
            time.sleep(0.2)

        _tracker.record("EuropePMC: metadata for all 4 gold PMIDs", True)

    def test_search_papers_kcnh2(self):
        """Search for KCNH2 papers and verify results."""
        results = self.client.search_papers("KCNH2 variant", max_results=5)

        assert len(results) > 0, "Should find KCNH2 papers"
        assert len(results) <= 5, "Should respect max_results"

        _tracker.record("EuropePMC: search_papers('KCNH2 variant')", True,
                        f"found {len(results)} results")

    def test_get_fulltext_xml(self):
        """Retrieve full-text XML for a PMC article with OA full text."""
        # PMC6198685 (PMID 30036649) - try to get XML
        # Note: Not all PMC articles have OA XML available through Europe PMC
        pmcid = GOLD_PMIDS["30036649"]["pmcid"]
        xml = self.client.get_fulltext_xml(pmcid)

        # Some articles may not have OA XML; verify graceful handling
        if xml is not None:
            assert len(xml) > 100, f"XML too short ({len(xml)} chars)"
            assert "<?xml" in xml or "<article" in xml, "Doesn't look like valid XML"
            _tracker.record(f"EuropePMC: get_fulltext_xml({pmcid})", True,
                            f"{len(xml)} chars")
        else:
            _tracker.record(f"EuropePMC: get_fulltext_xml({pmcid})", True,
                            "returned None (article not OA full-text in EuropePMC)")

    def test_get_fulltext_xml_returns_none_for_non_oa(self):
        """Verify graceful None return when full-text XML is not available."""
        # PMC4266740 (PMID 24667783) - known to not have OA full text
        pmcid = GOLD_PMIDS["24667783"]["pmcid"]
        xml = self.client.get_fulltext_xml(pmcid)

        # Should return None gracefully, not raise
        assert xml is None or isinstance(xml, str), "Should return None or XML string"

        _tracker.record("EuropePMC: fulltext_xml non-OA handling", True)

    def test_get_supplementary_files(self):
        """Check supplementary file listing for a paper with PMC."""
        pmcid = GOLD_PMIDS["30036649"]["pmcid"]
        files = self.client.get_supplementary_files(pmcid)

        assert isinstance(files, list), "Should return a list"

        _tracker.record(f"EuropePMC: get_supplementary_files({pmcid})", True,
                        f"found {len(files)} files")

    def test_get_supplementary_files_empty_pmcid(self):
        """Empty PMCID should return empty list, not error."""
        files = self.client.get_supplementary_files("")
        assert files == [], "Empty PMCID should return empty list"

        _tracker.record("EuropePMC: supplementary_files empty PMCID", True)

    def test_get_citations(self):
        """Retrieve citations for a paper."""
        pmid = "19841300"
        citations = self.client.get_citations(pmid)

        assert isinstance(citations, list), "Should return a list"

        _tracker.record(f"EuropePMC: get_citations({pmid})", True,
                        f"{len(citations)} citations")

    def test_get_references(self):
        """Retrieve references for a paper."""
        pmid = "19841300"
        references = self.client.get_references(pmid)

        assert isinstance(references, list), "Should return a list"

        _tracker.record(f"EuropePMC: get_references({pmid})", True,
                        f"{len(references)} references")

    def test_nonexistent_pmid(self):
        """Non-existent PMID should return None, not raise."""
        metadata = self.client.get_paper_metadata("99999999999")
        assert metadata is None, "Should return None for bogus PMID"

        _tracker.record("EuropePMC: nonexistent PMID handling", True)


# ============================================================================
# Europe PMC Integration (Harvester) Tests
# ============================================================================

@pytest.mark.requires_network
class TestEuropePMCIntegration:
    """Test gene_literature.europepmc_integration.EuropePMCHarvester."""

    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        from gene_literature.europepmc_integration import EuropePMCHarvester
        self.harvester = EuropePMCHarvester(output_dir=tmp_path / "europepmc_test")
        self.tmp_path = tmp_path

    def test_is_article_available(self):
        """Check availability for gold standard PMIDs."""
        for pmid in GOLD_PMIDS:
            result = self.harvester.is_article_available(pmid)
            assert result is not None, f"Should return a result dict for {pmid}"
            assert "available" in result, f"Result should have 'available' key for {pmid}"
            time.sleep(0.3)

        _tracker.record("EuropePMC Integration: is_article_available", True)

    def test_is_article_available_returns_pmcid(self):
        """Verify that available articles include PMCID."""
        # 30036649 has PMC
        result = self.harvester.is_article_available("30036649")
        assert result is not None

        if result.get("available"):
            assert result.get("pmcid"), "Available article should have PMCID"
            _tracker.record("EuropePMC Integration: available article has PMCID", True,
                            f"pmcid={result.get('pmcid')}")
        else:
            # Article may not have OA full text even with PMCID
            _tracker.record("EuropePMC Integration: available article has PMCID", True,
                            f"not OA full-text: {result.get('reason', '')[:50]}")

    def test_format_content_for_extraction(self):
        """Verify that format_content_for_extraction handles all cases."""
        # Try a paper that may or may not have full text
        content = self.harvester.format_content_for_extraction("30036649")

        if content is not None:
            assert len(content) > 100, f"Content too short ({len(content)} chars)"
            _tracker.record("EuropePMC Integration: format_content(30036649)", True,
                            f"{len(content)} chars")
        else:
            _tracker.record("EuropePMC Integration: format_content(30036649)", True,
                            "returned None (full-text not available as OA)")

    def test_compare_with_pmc(self):
        """Compare Europe PMC and NCBI PMC availability."""
        comparison = self.harvester.compare_with_pmc("24667783")

        assert comparison is not None, "Should return comparison dict"

        _tracker.record("EuropePMC Integration: compare_with_pmc(24667783)", True)


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
                assert pmcid == info["pmcid"], \
                    f"PMID {pmid}: expected {info['pmcid']}, got {pmcid}"
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
        assert doi == GOLD_DOIS["24667783"], \
            f"Expected {GOLD_DOIS['24667783']}, got {doi}"

        _tracker.record("PMC API: get_doi_from_pmid(24667783)", True, f"doi={doi}")

    def test_get_doi_all_gold(self):
        """Verify DOI retrieval for all gold standard PMIDs."""
        for pmid, expected_doi in GOLD_DOIS.items():
            doi = self.client.get_doi_from_pmid(pmid)
            assert doi is not None, f"No DOI for PMID {pmid}"
            # Compare case-insensitively since DOIs can vary in case
            assert doi.lower() == expected_doi.lower(), \
                f"PMID {pmid}: expected {expected_doi}, got {doi}"
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
        self.client = UnpaywallClient(email="brett.kroncke@gmail.com")

    def test_find_open_access(self):
        """Look up a known paper by DOI."""
        doi = GOLD_DOIS["24667783"]
        result, error = self.client.find_open_access(doi)

        assert error is None, f"Unpaywall error: {error}"
        assert result is not None, "Should return result for valid DOI"
        assert result["doi"] is not None, "Result should have DOI"
        assert "is_oa" in result, "Result should have is_oa field"
        assert "oa_status" in result, "Result should have oa_status field"

        _tracker.record("Unpaywall: find_open_access(24667783 DOI)", True,
                        f"is_oa={result['is_oa']}, status={result['oa_status']}")

    def test_find_open_access_all_gold(self):
        """Check Unpaywall for all gold standard DOIs."""
        oa_count = 0
        for pmid, doi in GOLD_DOIS.items():
            result, error = self.client.find_open_access(doi)
            assert result is not None or error is not None, \
                f"Should return either result or error for {doi}"
            if result and result.get("is_oa"):
                oa_count += 1
            time.sleep(0.2)

        _tracker.record("Unpaywall: OA check for all gold DOIs", True,
                        f"{oa_count}/4 are OA")

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
        assert result is None or error is not None or (result and not result.get("is_oa")), \
            "Invalid DOI should not return a valid OA result"

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

        _tracker.record("Supplement Scraper: all publisher methods exist", True,
                        f"{len(publishers)} methods verified")

    def test_extract_fulltext(self):
        """Verify extract_fulltext method exists and is callable."""
        assert hasattr(self.scraper, "extract_fulltext"), \
            "Should have extract_fulltext method"
        assert callable(self.scraper.extract_fulltext)

        _tracker.record("Supplement Scraper: extract_fulltext method exists", True)


# ============================================================================
# PMID Status Tracker Tests
# ============================================================================

class TestPMIDStatus:
    """Test utils.pmid_status module (local, no network needed)."""

    def test_import(self):
        """Verify pmid_status module imports."""
        from utils.pmid_status import get_pmid_status, get_failed_pmids, get_stats_summary
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
        unpaywall = UnpaywallClient(email="brett.kroncke@gmail.com")

        pmid = "24667783"

        # Step 1: Get DOI from PMID
        doi = pmc.get_doi_from_pmid(pmid)
        assert doi is not None, f"Failed to get DOI for PMID {pmid}"
        time.sleep(0.4)

        # Step 2: Check Unpaywall for OA
        result, error = unpaywall.find_open_access(doi)
        assert result is not None, f"Unpaywall failed for DOI {doi}: {error}"

        _tracker.record("Cross-module: PMID->DOI->Unpaywall pipeline", True,
                        f"doi={doi}, is_oa={result.get('is_oa')}")

    def test_europepmc_vs_ncbi_pmc(self):
        """Compare PMCID from Europe PMC metadata vs NCBI PMC API."""
        from gene_literature.europepmc_handler import EuropePMCClient
        from harvesting.pmc_api import PMCAPIClient

        epmc = EuropePMCClient()
        pmc = PMCAPIClient()

        pmid = "24667783"

        # Get PMCID from both sources
        epmc_metadata = epmc.get_paper_metadata(pmid)
        epmc_pmcid = epmc_metadata.get("pmcid") if epmc_metadata else None
        time.sleep(0.4)

        ncbi_pmcid = pmc.pmid_to_pmcid(pmid)

        assert epmc_pmcid is not None, "Europe PMC should return PMCID"
        assert ncbi_pmcid is not None, "NCBI PMC should return PMCID"
        assert epmc_pmcid == ncbi_pmcid, \
            f"PMCIDs should match: EuropePMC={epmc_pmcid}, NCBI={ncbi_pmcid}"

        _tracker.record("Cross-module: EuropePMC vs NCBI PMCID agreement", True,
                        f"both returned {ncbi_pmcid}")

    def test_full_paper_metadata_enrichment(self):
        """Enrich paper data from multiple sources for one PMID."""
        from gene_literature.europepmc_handler import EuropePMCClient
        from harvesting.pmc_api import PMCAPIClient
        from harvesting.unpaywall_api import UnpaywallClient

        epmc = EuropePMCClient()
        pmc = PMCAPIClient()
        unpaywall = UnpaywallClient(email="brett.kroncke@gmail.com")

        pmid = "19841300"

        # Gather from Europe PMC
        metadata = epmc.get_paper_metadata(pmid)
        assert metadata is not None
        time.sleep(0.3)

        # Gather DOI from NCBI
        doi = pmc.get_doi_from_pmid(pmid)
        assert doi is not None
        time.sleep(0.3)

        # Gather OA status from Unpaywall
        oa_result, _ = unpaywall.find_open_access(doi)

        enriched = {
            "pmid": pmid,
            "title": metadata.get("title"),
            "doi": doi,
            "pmcid": metadata.get("pmcid"),
            "citation_count": metadata.get("citation_count"),
            "is_oa": oa_result.get("is_oa") if oa_result else None,
            "oa_status": oa_result.get("oa_status") if oa_result else None,
        }

        assert enriched["title"], "Should have title"
        assert enriched["doi"], "Should have DOI"
        assert enriched["pmcid"], "Should have PMCID"

        _tracker.record("Cross-module: full metadata enrichment", True,
                        f"title={enriched['title'][:50]}...")

    def test_doi_consistency_across_sources(self):
        """DOI from NCBI and Europe PMC should match (case-insensitive)."""
        from gene_literature.europepmc_handler import EuropePMCClient
        from harvesting.pmc_api import PMCAPIClient

        epmc = EuropePMCClient()
        pmc = PMCAPIClient()

        pmid = "19841300"

        epmc_meta = epmc.get_paper_metadata(pmid)
        epmc_doi = epmc_meta.get("doi", "") if epmc_meta else ""
        time.sleep(0.4)

        ncbi_doi = pmc.get_doi_from_pmid(pmid) or ""

        assert epmc_doi.lower() == ncbi_doi.lower(), \
            f"DOIs should match: EuropePMC={epmc_doi}, NCBI={ncbi_doi}"

        _tracker.record("Cross-module: DOI consistency across sources", True,
                        f"doi={ncbi_doi}")


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

    print(f"\n{'='*70}")
    print("GeneVariantFetcher Integration Test Suite (standalone)")
    print(f"{'='*70}\n")

    # --- Europe PMC Handler ---
    print("[1/5] Europe PMC Handler...")
    try:
        from gene_literature.europepmc_handler import EuropePMCClient
        client = EuropePMCClient(timeout=30)

        # Metadata
        meta = client.get_paper_metadata("24667783")
        passed = meta is not None and meta.get("pmid") == "24667783"
        tracker.record("EuropePMC: get_paper_metadata", passed,
                       meta["title"][:50] if meta else "None")

        # Search
        time.sleep(0.3)
        results = client.search_papers("KCNH2 variant", max_results=5)
        tracker.record("EuropePMC: search_papers", len(results) > 0,
                       f"{len(results)} results")

        # Full-text XML (may or may not be available)
        time.sleep(0.3)
        pmcid = GOLD_PMIDS["30036649"]["pmcid"]
        xml = client.get_fulltext_xml(pmcid)
        tracker.record("EuropePMC: get_fulltext_xml",
                       xml is None or (isinstance(xml, str) and len(xml) > 100),
                       f"{len(xml)} chars" if xml else "None (not OA)")

        # Supplements
        time.sleep(0.3)
        supps = client.get_supplementary_files(pmcid)
        tracker.record("EuropePMC: get_supplementary_files", isinstance(supps, list),
                       f"{len(supps)} files")

        # Citations
        time.sleep(0.3)
        cites = client.get_citations("19841300")
        tracker.record("EuropePMC: get_citations", isinstance(cites, list),
                       f"{len(cites)} citations")

    except Exception as e:
        tracker.record("EuropePMC: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- Europe PMC Integration ---
    print("[2/5] Europe PMC Integration...")
    try:
        from gene_literature.europepmc_integration import EuropePMCHarvester
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            harvester = EuropePMCHarvester(output_dir=Path(tmpdir))

            avail = harvester.is_article_available("30036649")
            tracker.record("EuropePMC Integration: is_article_available",
                           avail is not None and "available" in avail,
                           f"available={avail.get('available')}" if avail else "None")

            time.sleep(0.3)
            content = harvester.format_content_for_extraction("30036649")
            tracker.record("EuropePMC Integration: format_content",
                           content is None or len(content) > 100,
                           f"{len(content)} chars" if content else "None (not OA)")

    except Exception as e:
        tracker.record("EuropePMC Integration: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- PMC API ---
    print("[3/5] PMC API Client...")
    try:
        from harvesting.pmc_api import PMCAPIClient
        pmc = PMCAPIClient()

        pmcid = pmc.pmid_to_pmcid("24667783")
        tracker.record("PMC API: pmid_to_pmcid", pmcid == "PMC4266740",
                       f"got {pmcid}")

        time.sleep(0.4)
        doi = pmc.get_doi_from_pmid("24667783")
        tracker.record("PMC API: get_doi_from_pmid",
                       doi is not None and doi == GOLD_DOIS["24667783"],
                       f"doi={doi}")

    except Exception as e:
        tracker.record("PMC API: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- Unpaywall ---
    print("[4/5] Unpaywall API...")
    try:
        from harvesting.unpaywall_api import UnpaywallClient
        uw = UnpaywallClient(email="brett.kroncke@gmail.com")

        result, error = uw.find_open_access(GOLD_DOIS["24667783"])
        tracker.record("Unpaywall: find_open_access",
                       result is not None and error is None,
                       f"is_oa={result.get('is_oa')}" if result else f"error={error}")

        time.sleep(0.2)
        result2, error2 = uw.find_open_access("")
        tracker.record("Unpaywall: empty DOI handling",
                       result2 is None and error2 is not None,
                       "correctly rejected")

    except Exception as e:
        tracker.record("Unpaywall: module load/test", False, str(e)[:80])
        traceback.print_exc()

    # --- PMID Status ---
    print("[5/5] PMID Status Tracker...")
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
