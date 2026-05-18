from gene_literature.supplements import PMCSupplementFetcher, UnifiedSupplementFetcher
from gene_literature.supplements.base import SupplementFile


def test_pmc_fetch_uses_supplied_pmcid_without_resolution(monkeypatch):
    fetcher = PMCSupplementFetcher()
    calls = []

    def fail_resolve(_pmid):
        raise AssertionError("PMCID should not be re-resolved when provided")

    def fake_europepmc(pmcid):
        calls.append(pmcid)
        return [
            SupplementFile(
                url="https://example.org/supplement.zip",
                name="supplement.zip",
                source="pmc_europepmc",
                pmcid=pmcid,
            )
        ]

    monkeypatch.setattr(fetcher, "_resolve_pmcid", fail_resolve)
    monkeypatch.setattr(fetcher, "_fetch_europepmc_supplements", fake_europepmc)

    results = fetcher.fetch("123456", pmcid="PMC123456")

    assert calls == ["PMC123456"]
    assert len(results) == 1
    assert results[0].pmcid == "PMC123456"


def test_pmc_fetch_normalizes_numeric_pmcid(monkeypatch):
    fetcher = PMCSupplementFetcher()
    calls = []

    monkeypatch.setattr(
        fetcher,
        "_fetch_europepmc_supplements",
        lambda pmcid: calls.append(pmcid) or [],
    )
    monkeypatch.setattr(fetcher, "_fetch_from_pmc_xml", lambda _pmcid: [])
    monkeypatch.setattr(fetcher, "_fetch_from_oa_service", lambda _pmcid: [])

    assert fetcher.fetch("123456", pmcid="123456") == []
    assert calls == ["PMC123456"]


def test_unified_fetcher_passes_known_pmcid_to_pmc_tier(monkeypatch):
    fetcher = UnifiedSupplementFetcher()
    calls = []

    def fake_pmc_fetch(pmid, doi="", pmcid=None):
        calls.append((pmid, doi, pmcid))
        return []

    monkeypatch.setattr(fetcher.pmc_fetcher, "fetch", fake_pmc_fetch)
    monkeypatch.setattr(fetcher.elsevier_fetcher, "fetch", lambda _pmid, _doi="": [])
    monkeypatch.setattr(fetcher.crossref_fetcher, "fetch", lambda _pmid, _doi="": [])

    assert fetcher.fetch_all("123456", "10.1000/example", pmcid="PMC123456") == []
    assert calls == [("123456", "10.1000/example", "PMC123456")]
