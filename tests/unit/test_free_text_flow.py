"""Tests for initial free-text fallback flow helper."""

from pathlib import Path

from harvesting.free_text_flow import initialize_free_text_access


class _DummyPMCAPI:
    def __init__(self, is_free=False, free_url=None):
        self._is_free = is_free
        self._free_url = free_url

    def is_free_full_text(self, pmid):
        return self._is_free, self._free_url


class _DummyUnpaywall:
    def __init__(self, oa=None, download_ok=False):
        self._oa = oa
        self._download_ok = download_ok

    def find_open_access(self, doi):
        return self._oa, None

    def download_pdf(self, pdf_url, out_path):
        Path(out_path).write_text("pdf", encoding="utf-8")
        return self._download_ok, None


class _DummyConverter:
    def __init__(self, markdown=""):
        self._markdown = markdown

    def pdf_to_markdown(self, path):
        return self._markdown


class _DummyProvider:
    def __init__(self, available=False):
        self.is_available = available

    def is_elsevier_doi(self, doi):
        return False

    def is_springer_doi(self, doi):
        return False

    def is_wiley_doi(self, doi):
        return False


class _DummyDOIResolver:
    def resolve_and_scrape_supplements(self, doi, pmid, scraper):
        return []


class _DummyLogger:
    def warning(self, msg):
        return None


def test_initialize_free_text_access_marks_paywalled_when_unavailable(tmp_path):
    status_calls = []
    paywalled_calls = []

    state = initialize_free_text_access(
        pmid="123",
        doi=None,
        output_dir=tmp_path,
        success_log=tmp_path / "success.csv",
        pmc_api=_DummyPMCAPI(is_free=False),
        unpaywall=_DummyUnpaywall(),
        converter=_DummyConverter(),
        elsevier_api=_DummyProvider(),
        springer_api=_DummyProvider(),
        wiley_api=_DummyProvider(),
        try_elsevier_api=lambda doi, pmid: (None, None),
        try_springer_api=lambda doi, pmid: (None, None),
        try_wiley_api=lambda doi, pmid: (None, None),
        doi_resolver=_DummyDOIResolver(),
        scraper=object(),
        write_pmid_status=lambda pmid, status, details: status_calls.append(
            (pmid, status, details)
        ),
        log_paywalled=lambda pmid, reason, url: paywalled_calls.append(
            (pmid, reason, url)
        ),
        logger=_DummyLogger(),
    )

    assert state.early_result == (False, "No PMCID", None)
    assert paywalled_calls
    assert status_calls and status_calls[0][1] == "paywalled"


def test_initialize_free_text_access_promotes_unpaywall_pdf(tmp_path):
    state = initialize_free_text_access(
        pmid="456",
        doi="10.1000/test",
        output_dir=tmp_path,
        success_log=tmp_path / "success.csv",
        pmc_api=_DummyPMCAPI(is_free=False),
        unpaywall=_DummyUnpaywall(
            oa={"pdf_url": "https://example.org/paper.pdf", "oa_status": "gold"},
            download_ok=True,
        ),
        converter=_DummyConverter(markdown=("a" * 600)),
        elsevier_api=_DummyProvider(),
        springer_api=_DummyProvider(),
        wiley_api=_DummyProvider(),
        try_elsevier_api=lambda doi, pmid: (None, None),
        try_springer_api=lambda doi, pmid: (None, None),
        try_wiley_api=lambda doi, pmid: (None, None),
        doi_resolver=_DummyDOIResolver(),
        scraper=object(),
        write_pmid_status=lambda pmid, status, details: None,
        log_paywalled=lambda pmid, reason, url: None,
        logger=_DummyLogger(),
    )

    assert state.early_result is None
    assert state.is_free
    assert state.free_url == "https://example.org/paper.pdf"
