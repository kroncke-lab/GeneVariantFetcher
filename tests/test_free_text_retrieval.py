"""Tests for detecting free full text availability."""

import importlib.util
import sys
import types
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
PMC_API_PATH = ROOT_DIR / "harvesting" / "pmc_api.py"


class _DummySession:
    def __init__(self):
        self.headers = {}

    def get(self, *args, **kwargs):
        raise RuntimeError("Network calls are disabled in tests")


requests_stub = types.SimpleNamespace(Session=_DummySession)
sys.modules.setdefault("requests", requests_stub)


class _DummyHandle:
    def read(self):
        return b""

    def close(self):
        return None


class _DummyEntrez:
    email = "test@example.com"
    tool = "pmc-test"
    api_key = None

    @staticmethod
    def elink(**kwargs):
        return _DummyHandle()

    @staticmethod
    def efetch(**kwargs):
        return _DummyHandle()

    @staticmethod
    def read(handle):
        return []


bio_module = types.ModuleType("Bio")
bio_module.Entrez = _DummyEntrez
sys.modules.setdefault("Bio", bio_module)
sys.modules.setdefault("Bio.Entrez", _DummyEntrez)

spec = importlib.util.spec_from_file_location("pmc_api", PMC_API_PATH)
pmc_api = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(pmc_api)
PMCAPIClient = pmc_api.PMCAPIClient


def test_free_text_detection(monkeypatch):
    """Ensure free full text detection flags expected PMIDs."""
    client = PMCAPIClient()

    test_cases = {
        "10336646": {
            "doi": "10.1038/10248",
            "linkouts": [
                {"provider": "Open Access", "url": "https://example.org/10336646", "attributes": ["free article"]}
            ],
        },
        "10827225": {
            "doi": "10.1054/jbin.2000.0128",
            "linkouts": [
                {"provider": "Publisher Free", "url": "https://example.org/10827225", "attributes": ["free"]}
            ],
        },
    }

    monkeypatch.setattr(client, "pmid_to_pmcid", lambda pmid: None)
    monkeypatch.setattr(client, "get_doi_from_pmid", lambda pmid: test_cases[pmid]["doi"])
    monkeypatch.setattr(client, "get_pubmed_linkout_urls", lambda pmid: test_cases[pmid]["linkouts"])

    for pmid, expected in test_cases.items():
        pmcid = client.pmid_to_pmcid(pmid)
        assert pmcid is None, f"PMID {pmid} should not have a PMCID"

        doi = client.get_doi_from_pmid(pmid)
        assert doi == expected["doi"], f"Incorrect DOI for {pmid}"

        is_free, free_url = client.is_free_full_text(pmid)
        assert is_free is True, f"PMID {pmid} should be detected as free"
        assert free_url == expected["linkouts"][0]["url"], "Free URL should come from linkouts"


def test_generic_full_text_indicator_is_ignored(monkeypatch):
    """A generic 'full text' indicator should not trigger a free result."""
    client = PMCAPIClient()

    monkeypatch.setattr(client, "pmid_to_pmcid", lambda pmid: None)
    monkeypatch.setattr(client, "get_pubmed_linkout_urls", lambda pmid: [{
        "provider": "Full Text at Publisher",
        "url": "https://example.org/nonfree",
        "attributes": ["full text"],
    }])

    is_free, free_url = client.is_free_full_text("00000000")
    assert is_free is False
    assert free_url is None
