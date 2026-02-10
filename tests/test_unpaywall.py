"""Tests for the Unpaywall integration.

These tests hit the live Unpaywall API and require network access.
"""

import pytest

from gene_literature.unpaywall import UnpaywallClient


@pytest.mark.requires_network
class TestUnpaywallClient:
    """Integration tests against the live Unpaywall API."""

    @pytest.fixture(autouse=True)
    def client(self):
        self.client = UnpaywallClient()

    def test_find_oa_url_returns_url(self):
        """Known OA paper (PLOS ONE) should return a URL."""
        doi = "10.1371/journal.pone.0185542"
        url = self.client.find_oa_url(doi)
        assert url is not None, f"Expected OA URL for {doi}"
        assert url.startswith("http")

    def test_find_oa_url_second_doi(self):
        """Oxford Academic OA paper should return a URL."""
        doi = "10.1093/eurheartj/ehv442"
        url = self.client.find_oa_url(doi)
        assert url is not None, f"Expected OA URL for {doi}"
        assert url.startswith("http")

    def test_is_open_access_true(self):
        doi = "10.1093/eurheartj/ehv442"
        assert self.client.is_open_access(doi) is True

    def test_get_oa_locations_non_empty(self):
        doi = "10.1093/eurheartj/ehv442"
        locs = self.client.get_oa_locations(doi)
        assert len(locs) > 0
        assert "url" in locs[0] or "url_for_pdf" in locs[0]

    def test_caching(self):
        """Second call for same DOI should use cache (no extra API hit)."""
        doi = "10.1093/eurheartj/ehv442"
        self.client.find_oa_url(doi)
        assert doi in self.client._cache
        # Second call should not raise or fail
        url = self.client.find_oa_url(doi)
        assert url is not None

    def test_nonexistent_doi(self):
        """Bogus DOI should return None, not raise."""
        url = self.client.find_oa_url("10.9999/this-doi-does-not-exist-xyz")
        assert url is None

    def test_is_open_access_bogus(self):
        assert self.client.is_open_access("10.9999/this-doi-does-not-exist-xyz") is False

    def test_get_oa_locations_bogus(self):
        assert self.client.get_oa_locations("10.9999/this-doi-does-not-exist-xyz") == []
