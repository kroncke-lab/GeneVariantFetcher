"""Tests for the unified supplement retrieval system.

Uses real PMIDs with known supplements to validate the fetcher pipeline.
These are integration tests that require network access.

Test PMIDs:
    31983221 - PMC7004454: Circulation OA paper with supplements (PDF + XLSX)
    24667783 - PMC4266740: Non-OA paper (tests graceful fallback)
    35443093 - Elsevier/GIM paper with supplements (mmc format)
"""

import pytest
from dotenv import load_dotenv

load_dotenv()

from gene_literature.supplements import (
    PMCSupplementFetcher,
    ElsevierSupplementFetcher,
    UnifiedSupplementFetcher,
    SupplementFile,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pmc_fetcher():
    return PMCSupplementFetcher(timeout=30)


@pytest.fixture(scope="module")
def elsevier_fetcher():
    return ElsevierSupplementFetcher(timeout=30)


@pytest.fixture(scope="module")
def unified_fetcher():
    return UnifiedSupplementFetcher(timeout=30)


# ---------------------------------------------------------------------------
# SupplementFile dataclass
# ---------------------------------------------------------------------------

class TestSupplementFile:
    def test_to_dict(self):
        sf = SupplementFile(url="https://example.com/mmc1.pdf", name="mmc1.pdf", source="test")
        d = sf.to_dict()
        assert d == {"url": "https://example.com/mmc1.pdf", "name": "mmc1.pdf"}

    def test_from_dict(self):
        d = {"url": "https://example.com/supp.xlsx", "name": "supp.xlsx"}
        sf = SupplementFile.from_dict(d, source="scraper")
        assert sf.url == d["url"]
        assert sf.name == d["name"]
        assert sf.source == "scraper"

    def test_extension(self):
        sf = SupplementFile(url="https://example.com/data.xlsx", name="data.xlsx")
        assert sf.extension == "xlsx"

    def test_normalized_url(self):
        sf = SupplementFile(url="https://example.com/file.pdf#page=2", name="file.pdf")
        assert sf.normalized_url == "https://example.com/file.pdf"

        sf2 = SupplementFile(url="https://example.com/file.pdf/", name="file.pdf")
        assert sf2.normalized_url == "https://example.com/file.pdf"


# ---------------------------------------------------------------------------
# PMC Fetcher (Tier 1)
# ---------------------------------------------------------------------------

@pytest.mark.requires_network
class TestPMCFetcher:
    """Test PMC supplement fetcher with known papers."""

    def test_pmcid_resolution(self, pmc_fetcher):
        """PMID 31983221 should resolve to PMC7004454."""
        pmcid = pmc_fetcher._resolve_pmcid("31983221")
        assert pmcid is not None, "Should resolve PMID to PMCID"
        assert "PMC" in pmcid, f"PMCID should start with PMC, got {pmcid}"

    def test_fetch_known_supplements(self, pmc_fetcher):
        """PMID 31983221 (PMC7004454) is an OA paper with known supplements."""
        results = pmc_fetcher.fetch("31983221")

        assert isinstance(results, list), "Should return a list"
        assert len(results) > 0, (
            "PMID 31983221 should have at least one supplement in PMC"
        )

        # Validate structure
        first = results[0]
        assert isinstance(first, SupplementFile)
        assert first.url.startswith("http"), f"URL should be absolute: {first.url}"
        assert first.name, "Should have a filename"
        assert first.source.startswith("pmc"), f"Source should be pmc_*: {first.source}"

        print(f"\n  PMC fetcher found {len(results)} supplements for PMID 31983221:")
        for s in results:
            print(f"    - {s.name} ({s.source}): {s.url[:80]}...")

    def test_non_oa_paper_graceful_fallback(self, pmc_fetcher):
        """PMID 24667783 is non-OA - should return empty or fall back."""
        results = pmc_fetcher.fetch("24667783")
        # Non-OA papers may have no supplements accessible via PMC
        assert isinstance(results, list), "Should return a list even for non-OA"

    def test_no_supplements_for_non_pmc(self, pmc_fetcher):
        """A paper not in PMC should return an empty list, not raise."""
        results = pmc_fetcher.fetch("99999999")
        assert results == [], "Non-existent PMID should return empty list"

    def test_legacy_dict_compat(self, pmc_fetcher):
        """SupplementFile.to_dict() should produce the legacy format."""
        results = pmc_fetcher.fetch("31983221")
        if results:
            d = results[0].to_dict()
            assert "url" in d
            assert "name" in d
            assert isinstance(d["url"], str)


# ---------------------------------------------------------------------------
# Elsevier Fetcher (Tier 2)
# ---------------------------------------------------------------------------

@pytest.mark.requires_network
class TestElsevierFetcher:
    """Test Elsevier supplement fetcher."""

    def test_available_check(self, elsevier_fetcher):
        """Should report availability based on API key presence."""
        # Just check the property works, value depends on .env
        assert isinstance(elsevier_fetcher.available, bool)

    def test_no_doi_returns_empty(self, elsevier_fetcher):
        """Without a DOI, Elsevier fetcher should return empty."""
        results = elsevier_fetcher.fetch("24667783", doi="")
        assert results == []

    @pytest.mark.skipif(
        not ElsevierSupplementFetcher().available,
        reason="ELSEVIER_API_KEY not configured"
    )
    def test_fetch_with_doi(self, elsevier_fetcher):
        """Test Elsevier fetch with a known Elsevier DOI."""
        # PMID 35443093 - Genetics in Medicine (Elsevier)
        results = elsevier_fetcher.fetch(
            "35443093",
            doi="10.1016/j.gim.2022.03.014"
        )
        # May or may not find supplements depending on API access level
        assert isinstance(results, list)
        if results:
            first = results[0]
            assert isinstance(first, SupplementFile)
            assert first.url.startswith("http")
            print(f"\n  Elsevier fetcher found {len(results)} supplements")


# ---------------------------------------------------------------------------
# Unified Fetcher
# ---------------------------------------------------------------------------

@pytest.mark.requires_network
class TestUnifiedFetcher:
    """Test the unified multi-tier fetcher."""

    def test_fetch_all_with_known_pmc_paper(self, unified_fetcher):
        """Unified fetcher should find supplements for PMID 31983221."""
        results = unified_fetcher.fetch_all("31983221")

        assert isinstance(results, list)
        assert len(results) > 0, "Unified fetcher should find supplements via PMC"

        # Check deduplication: no duplicate URLs
        urls = [s.normalized_url for s in results]
        assert len(urls) == len(set(urls)), "Results should be deduplicated by URL"

        print(f"\n  Unified fetcher found {len(results)} supplements for PMID 31983221:")
        for s in results:
            print(f"    [{s.source}] {s.name}: {s.url[:80]}...")

    def test_fetch_all_returns_empty_for_unknown(self, unified_fetcher):
        """Unknown PMID should return empty, not raise."""
        results = unified_fetcher.fetch_all("99999999")
        assert results == []

    def test_to_legacy_format(self, unified_fetcher):
        """Legacy format conversion should work."""
        results = unified_fetcher.fetch_all("31983221")
        legacy = unified_fetcher.to_legacy_format(results)

        assert isinstance(legacy, list)
        if legacy:
            assert isinstance(legacy[0], dict)
            assert "url" in legacy[0]
            assert "name" in legacy[0]

    def test_fetch_tier1_only(self, unified_fetcher):
        """Tier 1 fetch should only use PMC."""
        results = unified_fetcher.fetch_tier1("31983221")
        assert isinstance(results, list)
        for s in results:
            assert s.source.startswith("pmc"), f"Tier 1 result should be from PMC: {s.source}"

    def test_deduplication_across_tiers(self, unified_fetcher):
        """If the same supplement is found by multiple tiers, keep only one."""
        results = unified_fetcher.fetch_all("31983221")
        seen = set()
        for s in results:
            norm = s.normalized_url
            assert norm not in seen, f"Duplicate URL found: {norm}"
            seen.add(norm)
