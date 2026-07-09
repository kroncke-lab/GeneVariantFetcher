"""Tests for Google Scholar title-based PDF fallback."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
import requests

from harvesting.html_body_fetcher import HtmlBodyFetcher
from harvesting import scholar_pdf_fallback as spf
from harvesting.scholar_pdf_fallback import ScholarPDFRecovery
from scripts import fetch_paywalled as fp


class _Resp:
    def __init__(
        self,
        *,
        text: str = "",
        content: bytes = b"",
        status_code: int = 200,
        url: str = "https://example.org/",
        headers: dict | None = None,
        payload: dict | None = None,
    ):
        self.text = text
        self.content = content or text.encode("utf-8")
        self.status_code = status_code
        self.url = url
        self.headers = headers or {}
        self._payload = payload or {}

    def json(self):
        return self._payload


@pytest.fixture(autouse=True)
def reset_scholar_state(monkeypatch):
    monkeypatch.setattr(spf, "_LAST_REQUEST", 0.0)
    monkeypatch.setattr(spf, "_SCHOLAR_BLOCKED_UNTIL", 0.0)


def test_scholar_pdf_recovery_fetches_pdf_from_title_result(tmp_path: Path):
    calls: list[str] = []

    class _Session(requests.Session):
        def get(self, url, **kwargs):  # noqa: ANN001
            calls.append(url)
            if "scholar.google.com" in url:
                return _Resp(
                    text="""
                    <html><body>
                      <div class="gs_or_ggsm">
                        <a href="https://lab.example.edu/kcnh2-paper.pdf">[PDF] lab.example.edu</a>
                      </div>
                    </body></html>
                    """,
                    url="https://scholar.google.com/scholar?q=x",
                )
            if url == "https://lab.example.edu/kcnh2-paper.pdf":
                return _Resp(
                    content=b"%PDF-1.4\n" + (b"x" * 50_000),
                    headers={"content-type": "application/pdf"},
                    url=url,
                )
            raise AssertionError(f"unexpected URL {url}")

    converter = MagicMock()
    converter.pdf_to_markdown.return_value = "Methods results discussion " * 400
    client = ScholarPDFRecovery(
        session=_Session(),
        converter=converter,
        email="test@example.org",
        quality_gate=lambda text: (True, "ok") if text else (False, "empty"),
        min_request_interval=0,
    )

    result = client.recover(pmid="1", title="KCNH2 cohort paper")

    assert result.success is True
    assert result.source_url == "https://lab.example.edu/kcnh2-paper.pdf"
    assert converter.pdf_to_markdown.call_count == 1
    assert calls == [
        "https://scholar.google.com/scholar",
        "https://lab.example.edu/kcnh2-paper.pdf",
    ]


def test_scholar_pdf_recovery_treats_captcha_as_clean_miss():
    class _Session(requests.Session):
        def get(self, url, **kwargs):  # noqa: ANN001
            assert "scholar.google.com" in url
            return _Resp(
                text="Our systems have detected unusual traffic from your computer network. CAPTCHA.",
                url="https://scholar.google.com/sorry/index",
            )

    client = ScholarPDFRecovery(
        session=_Session(),
        converter=MagicMock(),
        email="test@example.org",
        quality_gate=lambda text: (True, "ok"),
        min_request_interval=0,
    )

    result = client.recover(pmid="1", title="Blocked paper")

    assert result.success is False
    assert result.error == "Scholar captcha/block page"
    assert result.attempts == [("scholar_search", "miss", "Scholar captcha/block page")]


def test_scholar_pdf_recovery_skips_search_during_block_backoff(monkeypatch):
    calls: list[str] = []

    class _Session(requests.Session):
        def get(self, url, **kwargs):  # noqa: ANN001
            calls.append(url)
            return _Resp(status_code=429, url="https://scholar.google.com/scholar")

    monkeypatch.setattr(spf.time, "time", lambda: 1000.0)
    monkeypatch.setattr(spf, "DEFAULT_SCHOLAR_BLOCK_BACKOFF_SECONDS", 300)
    client = ScholarPDFRecovery(
        session=_Session(),
        converter=MagicMock(),
        email="test@example.org",
        quality_gate=lambda text: (True, "ok"),
        min_request_interval=0,
    )

    first = client.recover(pmid="1", title="Blocked paper")
    second = client.recover(pmid="1", title="Blocked paper")

    assert first.error == "Scholar blocked HTTP 429"
    assert second.error == "Scholar temporarily disabled after block (300s remaining)"
    assert calls == ["https://scholar.google.com/scholar"]


def test_scholar_pdf_recovery_rejects_tiny_pdf_before_conversion():
    class _Session(requests.Session):
        def get(self, url, **kwargs):  # noqa: ANN001
            if "scholar.google.com" in url:
                return _Resp(
                    text='<a href="https://lab.example.edu/error.pdf">[PDF]</a>',
                    url="https://scholar.google.com/scholar?q=x",
                )
            return _Resp(
                content=b"%PDF-1.4 tiny",
                headers={"content-type": "application/pdf"},
                url=url,
            )

    converter = MagicMock()
    client = ScholarPDFRecovery(
        session=_Session(),
        converter=converter,
        email="test@example.org",
        quality_gate=lambda text: (True, "ok") if text else (False, "empty"),
        min_request_interval=0,
    )

    result = client.recover(pmid="1", title="Tiny PDF")

    assert result.success is False
    assert "PDF too small" in result.attempts[-1][2]
    converter.pdf_to_markdown.assert_not_called()


def test_scholar_pdf_recovery_treats_converter_sentinel_as_failure():
    class _Session(requests.Session):
        def get(self, url, **kwargs):  # noqa: ANN001
            if "scholar.google.com" in url:
                return _Resp(
                    text='<a href="https://lab.example.edu/paper.pdf">[PDF]</a>',
                    url="https://scholar.google.com/scholar?q=x",
                )
            return _Resp(
                content=b"%PDF-1.4\n" + (b"x" * 50_000),
                headers={"content-type": "application/pdf"},
                url=url,
            )

    converter = MagicMock()
    converter.pdf_to_markdown.return_value = "[PDF file available at: failed]"
    client = ScholarPDFRecovery(
        session=_Session(),
        converter=converter,
        email="test@example.org",
        quality_gate=lambda text: (True, "ok") if text else (False, "empty"),
        min_request_interval=0,
    )

    result = client.recover(pmid="1", title="Broken PDF")

    assert result.success is False
    assert result.attempts[-1][2].endswith("PDF conversion failed")


def test_scholar_pdf_candidate_parser_blocks_known_bad_hosts():
    html = """
    <a href="https://sci-hub.se/example.pdf">[PDF]</a>
    <a href="https://www.researchgate.net/publication/example.pdf">[PDF]</a>
    <a href="https://lab.example.edu/paper.pdf">[PDF]</a>
    """

    assert spf._extract_pdf_urls(html) == ["https://lab.example.edu/paper.pdf"]


@pytest.mark.parametrize(
    ("doi", "expected_prior_attempt"),
    [
        ("10.1000/example", ("unpaywall_pdf", "miss", "miss")),
        (None, ("unpaywall_html", "skip", "no DOI")),
    ],
)
def test_html_body_fetcher_tries_scholar_after_standard_misses_or_no_doi(
    monkeypatch, doi, expected_prior_attempt
):
    fetcher = HtmlBodyFetcher(
        email="test@example.org",
        session=requests.Session(),
        min_request_interval=0,
    )
    body = "Methods results discussion references " * 300
    fetcher.pmc.get_doi_from_pmid = lambda pmid: doi
    monkeypatch.setattr(fetcher, "_try_europepmc_xml", lambda pmid: (None, "miss"))
    monkeypatch.setattr(fetcher, "_try_ncbi_elink_pmc", lambda pmid: (None, "miss"))
    monkeypatch.setattr(fetcher, "_try_unpaywall_html", lambda doi: (None, "miss"))
    monkeypatch.setattr(fetcher, "_try_unpaywall_pdf", lambda doi: (None, "miss"))
    monkeypatch.setattr(fetcher, "_try_scholar_pdf", lambda pmid: (body, "pdf"))

    result = fetcher.fetch_body("12345678")

    assert result.success is True
    assert result.source == "scholar_pdf"
    assert result.markdown == body
    assert expected_prior_attempt in result.attempts
    assert result.attempts[-1] == ("scholar_pdf", "ok", "pdf")


def test_fetch_paywalled_scholar_fallback_writes_canonical_files(
    monkeypatch, tmp_path: Path
):
    markdown = "Methods results discussion references " * 300
    fake_result = SimpleNamespace(
        success=True,
        markdown=markdown,
        source_url="https://lab.example.edu/paper.pdf",
        title="Recovered title",
        detail="https://lab.example.edu/paper.pdf",
        attempts=[("scholar_pdf", "ok", "accepted")],
    )
    monkeypatch.setattr(fp, "try_scholar_pdf", lambda **kwargs: fake_result)

    row = fp.try_scholar_pdf_fallback(
        pmid="99999999",
        output_dir=tmp_path,
        session=requests.Session(),
        converter=MagicMock(),
    )

    assert row is not None
    assert row["reason"].startswith("Google Scholar PDF (Recovered title):")
    assert row["final_url"] == "https://lab.example.edu/paper.pdf"
    assert (tmp_path / "99999999_FULL_CONTEXT.md").read_text() == markdown
    assert (tmp_path / "99999999" / "FULL_CONTEXT.md").read_text() == markdown


def test_fetch_one_uses_scholar_when_no_doi(monkeypatch, tmp_path: Path):
    scholar_row = {
        "reason": "Google Scholar PDF: https://lab.example.edu/paper.pdf",
        "chars": 9000,
        "final_url": "https://lab.example.edu/paper.pdf",
        "path": str(tmp_path / "1_FULL_CONTEXT.md"),
        "canonical_path": str(tmp_path / "1_FULL_CONTEXT.md"),
        "per_pmid_path": str(tmp_path / "1" / "FULL_CONTEXT.md"),
    }
    monkeypatch.setattr(
        fp, "try_scholar_pdf_fallback", lambda *args, **kwargs: scholar_row
    )
    # No DOI on input AND none resolvable from the PMID -> Scholar is the last
    # resort. Stub the resolver so the path is exercised without network.
    monkeypatch.setattr(fp, "resolve_doi_for_pmid", lambda _pmid: None)
    fetcher = SimpleNamespace(converter=MagicMock(), session=requests.Session())

    row = fp.fetch_one(fetcher, "1", None, tmp_path)

    assert row["outcome"] == "success_via_scholar_pdf"
    assert row["strategy"] == "google_scholar_pdf"
    assert row["final_url"] == "https://lab.example.edu/paper.pdf"
