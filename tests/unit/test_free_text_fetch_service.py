"""Tests for free-text main content fetch service."""

from harvesting.free_text_fetch_service import fetch_main_content_for_free_text


class _Provider:
    is_available = False

    def is_elsevier_doi(self, doi):
        return False

    def is_springer_doi(self, doi):
        return False

    def is_wiley_doi(self, doi):
        return False

    def is_elsevier_url(self, url):
        return False

    def extract_pii_from_url(self, url):
        return None

    def is_wiley_url(self, url):
        return False

    def extract_doi_from_url(self, url):
        return None


class _Session:
    def get(self, *args, **kwargs):
        raise AssertionError("Session.get should not be called in this test")


class _Scraper:
    pass


class _DOIResolver:
    def __init__(self, markdown=None):
        self.markdown = markdown

    def resolve_and_scrape_supplements(self, doi, pmid, scraper):
        return []

    def resolve_and_fetch_fulltext(self, doi, pmid, scraper):
        return self.markdown, "https://doi.org/" + doi, []


def test_returns_error_when_no_doi_and_no_free_url():
    result = fetch_main_content_for_free_text(
        pmid="1",
        doi=None,
        free_url=None,
        suspicious_free_url_domains=set(),
        elsevier_api=_Provider(),
        springer_api=_Provider(),
        wiley_api=_Provider(),
        elsevier_client=object(),
        wiley_client=object(),
        session=_Session(),
        scraper=_Scraper(),
        doi_resolver=_DOIResolver(),
        try_elsevier_api=lambda doi, pmid: (None, None),
        try_wiley_api=lambda doi, pmid: (None, None),
        try_springer_api=lambda doi, pmid: (None, None),
        validate_content_quality=lambda text, url: (True, "ok"),
        log_paywalled=lambda pmid, reason, url: None,
    )
    assert result.early_result == (False, "No DOI or URL for free text", None)


def test_returns_error_for_suspicious_free_url():
    result = fetch_main_content_for_free_text(
        pmid="2",
        doi=None,
        free_url="https://antibodies.cancer.gov/x",
        suspicious_free_url_domains={"antibodies.cancer.gov"},
        elsevier_api=_Provider(),
        springer_api=_Provider(),
        wiley_api=_Provider(),
        elsevier_client=object(),
        wiley_client=object(),
        session=_Session(),
        scraper=_Scraper(),
        doi_resolver=_DOIResolver(),
        try_elsevier_api=lambda doi, pmid: (None, None),
        try_wiley_api=lambda doi, pmid: (None, None),
        try_springer_api=lambda doi, pmid: (None, None),
        validate_content_quality=lambda text, url: (True, "ok"),
        log_paywalled=lambda pmid, reason, url: None,
    )
    assert result.early_result == (False, "Suspicious free URL", None)


def test_uses_doi_resolver_when_content_available():
    result = fetch_main_content_for_free_text(
        pmid="3",
        doi="10.1000/test",
        free_url=None,
        suspicious_free_url_domains=set(),
        elsevier_api=_Provider(),
        springer_api=_Provider(),
        wiley_api=_Provider(),
        elsevier_client=object(),
        wiley_client=object(),
        session=_Session(),
        scraper=_Scraper(),
        doi_resolver=_DOIResolver(markdown="abstract methods results references " * 60),
        try_elsevier_api=lambda doi, pmid: (None, None),
        try_wiley_api=lambda doi, pmid: (None, None),
        try_springer_api=lambda doi, pmid: (None, None),
        validate_content_quality=lambda text, url: (True, "valid"),
        log_paywalled=lambda pmid, reason, url: None,
    )
    assert result.early_result is None
    assert result.main_markdown is not None
