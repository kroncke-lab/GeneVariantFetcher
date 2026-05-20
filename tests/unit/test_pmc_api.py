"""Unit tests for PMC API helper behavior."""

from harvesting import pmc_api
from harvesting.pmc_api import PMCAPIClient


class _DummyHandle:
    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


def test_pmid_to_pmcid_retries_transient_entrez_read(monkeypatch):
    """Transient Entrez read failures should not drop PMC-linked papers."""

    calls = {"read": 0}

    def fake_elink(**kwargs):
        return _DummyHandle()

    def fake_read(handle):
        calls["read"] += 1
        if calls["read"] == 1:
            raise RuntimeError("Read failed: EOF")
        return [{"LinkSetDb": [{"Link": [{"Id": "123456"}]}]}]

    monkeypatch.setattr(PMCAPIClient, "_rate_limit", lambda self: None)
    monkeypatch.setattr(pmc_api.Entrez, "elink", fake_elink)
    monkeypatch.setattr(pmc_api.Entrez, "read", fake_read)
    monkeypatch.setattr(pmc_api.time, "sleep", lambda seconds: None)

    assert PMCAPIClient().pmid_to_pmcid("999") == "PMC123456"
    assert calls["read"] == 2


def test_is_free_full_text_skips_non_article_linkouts(monkeypatch):
    """Database/tool LinkOuts should not become generic free-text candidates."""

    def fake_elink(**kwargs):
        return _DummyHandle()

    def fake_read(handle):
        return [
            {
                "IdUrlList": {
                    "IdUrlSet": [
                        {
                            "ObjUrl": [
                                {
                                    "Attribute": ["free full text"],
                                    "Url": "https://scite.ai/reports/27225049",
                                },
                                {
                                    "Attribute": ["free full text"],
                                    "Url": "https://www.ahajournals.org/doi/10.1161/test",
                                },
                            ]
                        }
                    ]
                }
            }
        ]

    monkeypatch.setattr(PMCAPIClient, "_rate_limit", lambda self: None)
    monkeypatch.setattr(pmc_api.Entrez, "elink", fake_elink)
    monkeypatch.setattr(pmc_api.Entrez, "read", fake_read)

    is_free, url = PMCAPIClient().is_free_full_text("27225049")

    assert is_free is True
    assert url == "https://www.ahajournals.org/doi/10.1161/test"
