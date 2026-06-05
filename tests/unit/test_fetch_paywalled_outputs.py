import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from harvesting.paywall_context_enrichment import EnrichmentResult
from scripts import fetch_paywalled


def test_write_outputs_preserves_supplement_only_recovery(tmp_path, monkeypatch):
    result = SimpleNamespace(
        main_markdown=None,
        main_html="<html></html>",
        supp_files=[{"url": "https://example.org/supp.doc", "name": "supp.doc"}],
        figure_paths=[],
        notes=[],
        error=None,
        publisher="generic",
        final_url="https://example.org/article",
    )

    def fake_enrich(**kwargs):
        assert "supplement-only recovery" in kwargs["body_markdown"]
        return EnrichmentResult(
            unified_markdown=kwargs["body_markdown"] + "VARIANT_TABLE",
            supplement_results=[SimpleNamespace(downloaded=True)],
        )

    monkeypatch.setattr(fetch_paywalled, "enrich_paywall_full_context", fake_enrich)

    canonical, per_pmid, body, enrichment = fetch_paywalled.write_outputs(
        "12345",
        result,
        tmp_path,
    )

    assert body is None
    assert enrichment.supplement_count == 1
    assert "VARIANT_TABLE" in canonical.read_text(encoding="utf-8")
    assert per_pmid.read_text(encoding="utf-8") == canonical.read_text(encoding="utf-8")


def test_write_outputs_passes_supplement_download_fallback(tmp_path, monkeypatch):
    result = SimpleNamespace(
        main_markdown="# MAIN TEXT\n\nbody",
        main_html="<html></html>",
        supp_files=[{"url": "https://example.org/supp.xlsx", "name": "supp.xlsx"}],
        figure_paths=[],
        notes=[],
        error=None,
        publisher="generic",
        final_url="https://example.org/article",
    )

    fallback = object()

    def fake_enrich(**kwargs):
        assert kwargs["supplement_download_fallback"] is fallback
        assert kwargs["source_url"] == "https://example.org/article"
        return EnrichmentResult(unified_markdown=kwargs["body_markdown"])

    monkeypatch.setattr(fetch_paywalled, "enrich_paywall_full_context", fake_enrich)

    fetch_paywalled.write_outputs(
        "12345",
        result,
        tmp_path,
        supplement_download_fallback=fallback,
    )


def test_browser_supplement_download_fallback_uses_context_request(tmp_path):
    payload = b"PK\x03\x04xlsx"
    calls = []

    class _Response:
        ok = True
        status = 200
        headers = {
            "content-type": (
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        }

        def body(self):
            return payload

    class _Request:
        def get(self, url, **kwargs):
            calls.append((url, kwargs))
            return _Response()

    class _Context:
        request = _Request()

    class _Page:
        context = _Context()

        def set_default_timeout(self, _timeout_ms):
            return None

        def goto(self, *_args, **_kwargs):
            return None

    class _PageManager:
        def __enter__(self):
            return _Page()

        def __exit__(self, *_exc):
            return False

    class _Pool:
        def page(self):
            return _PageManager()

    fetcher = SimpleNamespace(_pool=_Pool(), _timeout_s=90)
    fallback = fetch_paywalled.make_browser_supplement_download_fallback(fetcher)
    assert fallback is not None

    out = tmp_path / "supplements" / "table-s1.xlsx"
    ok = fallback(
        "https://publisher.example/supp/table-s1.xlsx",
        out,
        "24680",
        "table-s1.xlsx",
        {"source_url": "https://publisher.example/article"},
    )

    assert ok is True
    assert out.read_bytes() == payload
    assert calls[0][0] == "https://publisher.example/supp/table-s1.xlsx"
    assert calls[0][1]["headers"]["Referer"] == "https://publisher.example/article"


def test_elsevier_api_fallback_preserves_recovered_supplement_markdown(
    tmp_path, monkeypatch
):
    class _Client:
        is_available = True

        @staticmethod
        def is_elsevier_doi(_doi):
            return True

        def __init__(self, **_kwargs):
            return None

        def fetch_fulltext(self, **_kwargs):
            return "# MAIN TEXT\n\nAPI body\n", None

    monkeypatch.setattr(fetch_paywalled, "ElsevierAPIClient", _Client)
    monkeypatch.setattr(
        fetch_paywalled,
        "validate_article_content",
        lambda _markdown: (True, "ok"),
    )

    row = fetch_paywalled.try_elsevier_api_fallback(
        "11111",
        "10.1016/j.hrthm.2018.01.014",
        tmp_path,
        supplement_markdown="\n\n# SUPPLEMENTAL FILE 1: mmc1.docx\n\np.Arg1Gln\n",
        supplements_downloaded=1,
    )

    assert row is not None
    assert row["supplements_downloaded"] == 1
    text = Path(row["canonical_path"]).read_text(encoding="utf-8")
    assert "API body" in text
    assert "RECOVERED PUBLISHER SUPPLEMENTS" in text
    assert "p.Arg1Gln" in text


def test_wiley_api_fallback_writes_fulltext_artifacts(tmp_path, monkeypatch):
    class _Client:
        is_available = True

        @staticmethod
        def is_wiley_doi(_doi):
            return True

        def __init__(self, **_kwargs):
            return None

        def fetch_fulltext(self, **_kwargs):
            return "# MAIN TEXT\n\nWiley body with methods and results.\n", None

    monkeypatch.setattr(fetch_paywalled, "WileyAPIClient", _Client)
    monkeypatch.setattr(
        fetch_paywalled,
        "validate_article_content",
        lambda _markdown: (True, "ok"),
    )

    row = fetch_paywalled.try_wiley_api_fallback(
        "22222",
        "10.1111/j.1540-8167.2006.00349.x",
        tmp_path,
        final_url="https://onlinelibrary.wiley.com/doi/full/example",
    )

    assert row is not None
    assert row["outcome"] == "success_via_wiley_api"
    assert row["strategy"] == "wiley_api"
    text = Path(row["canonical_path"]).read_text(encoding="utf-8")
    assert "Wiley body" in text
    meta = tmp_path / "22222" / "result.json"
    assert '"publisher": "wiley_api"' in meta.read_text(encoding="utf-8")


def test_springer_api_fallback_writes_fulltext_artifacts(tmp_path, monkeypatch):
    class _Client:
        is_available = True

        @staticmethod
        def is_springer_doi(_doi):
            return True

        def __init__(self, **_kwargs):
            return None

        def fetch_article(self, _doi):
            return "# MAIN TEXT\n\nSpringer body with methods and results.\n", {}, None

    monkeypatch.setattr(fetch_paywalled, "SpringerAPIClient", _Client)
    monkeypatch.setattr(
        fetch_paywalled,
        "validate_article_content",
        lambda _markdown: (True, "ok"),
    )

    row = fetch_paywalled.try_springer_api_fallback(
        "33333",
        "10.1007/s00125-000-0000-0",
        tmp_path,
    )

    assert row is not None
    assert row["outcome"] == "success_via_springer_api"
    assert row["strategy"] == "springer_api"
    text = Path(row["canonical_path"]).read_text(encoding="utf-8")
    assert "Springer body" in text


def test_fetch_one_uses_publisher_api_when_no_browser_strategy(tmp_path, monkeypatch):
    canonical = tmp_path / "22222_FULL_CONTEXT.md"
    per_pmid = tmp_path / "22222" / "FULL_CONTEXT.md"
    per_pmid.parent.mkdir()
    canonical.write_text("# MAIN TEXT\n\napi body\n", encoding="utf-8")
    per_pmid.write_text(canonical.read_text(encoding="utf-8"), encoding="utf-8")

    fetcher = SimpleNamespace(
        fetch=lambda **_kwargs: (_ for _ in ()).throw(AssertionError("no fetch")),
        converter=object(),
        scraper=object(),
        session=object(),
    )
    scholar_calls = []

    monkeypatch.setattr(fetch_paywalled, "find_strategy", lambda **_kwargs: None)
    monkeypatch.setattr(
        fetch_paywalled,
        "try_publisher_api_fallback",
        lambda *_args, **_kwargs: {
            "strategy": "wiley_api",
            "outcome": "success_via_wiley_api",
            "api_label": "Wiley TDM API",
            "chars": 9000,
            "reason": "ok",
            "path": str(canonical),
            "canonical_path": str(canonical),
            "per_pmid_path": str(per_pmid),
            "supp_files": 0,
            "supplements_downloaded": 0,
            "final_url": "https://api.example/article",
        },
    )
    monkeypatch.setattr(
        fetch_paywalled,
        "try_scholar_pdf_fallback",
        lambda *_args, **_kwargs: scholar_calls.append(True),
    )

    row = fetch_paywalled.fetch_one(
        fetcher,
        "22222",
        "10.1111/j.1540-8167.2006.00349.x",
        tmp_path,
        pmc_session=object(),
    )

    assert row["outcome"] == "success_via_wiley_api"
    assert row["strategy"] == "wiley_api"
    assert row["reason"] == "Wiley TDM API fallback: ok"
    assert scholar_calls == []


def test_fetch_one_counts_generic_publisher_api_default_as_success(
    tmp_path, monkeypatch
):
    canonical = tmp_path / "44444_FULL_CONTEXT.md"
    per_pmid = tmp_path / "44444" / "FULL_CONTEXT.md"
    per_pmid.parent.mkdir()
    canonical.write_text("# MAIN TEXT\n\ngeneric api body\n", encoding="utf-8")
    per_pmid.write_text(canonical.read_text(encoding="utf-8"), encoding="utf-8")

    fetcher = SimpleNamespace(
        fetch=lambda **_kwargs: (_ for _ in ()).throw(AssertionError("no fetch")),
        converter=object(),
        scraper=object(),
        session=object(),
    )

    monkeypatch.setattr(fetch_paywalled, "find_strategy", lambda **_kwargs: None)
    monkeypatch.setattr(
        fetch_paywalled,
        "try_publisher_api_fallback",
        lambda *_args, **_kwargs: {
            "strategy": "publisher_api",
            "api_label": "Publisher API",
            "chars": 9000,
            "reason": "ok",
            "path": str(canonical),
            "canonical_path": str(canonical),
            "per_pmid_path": str(per_pmid),
            "supp_files": 0,
            "supplements_downloaded": 0,
            "final_url": "https://api.example/article",
        },
    )

    row = fetch_paywalled.fetch_one(
        fetcher,
        "44444",
        "10.0000/example",
        tmp_path,
        pmc_session=object(),
    )

    assert row["outcome"] == "success_via_publisher_api"
    assert row["outcome"] in fetch_paywalled.FETCH_SUCCESS_OUTCOMES


def test_ncbi_email_requires_env(monkeypatch):
    monkeypatch.delenv("ENTREZ_EMAIL", raising=False)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)

    with pytest.raises(RuntimeError, match="ENTREZ_EMAIL"):
        fetch_paywalled._ncbi_email()


def test_main_requires_explicit_pmids(monkeypatch, tmp_path):
    monkeypatch.setattr(
        sys,
        "argv",
        ["fetch_paywalled", "--output", str(tmp_path), "--no-cookies"],
    )

    with pytest.raises(SystemExit) as exc:
        fetch_paywalled.main()

    assert exc.value.code == 2


def test_main_wires_persistent_profile_and_auth_url(monkeypatch, tmp_path):
    captured = {}

    class _Pool:
        def __init__(
            self,
            *,
            cookies,
            headless,
            use_chrome_channel,
            persistent_profile_path=None,
        ):
            captured["pool"] = {
                "cookies": cookies,
                "headless": headless,
                "use_chrome_channel": use_chrome_channel,
                "persistent_profile_path": persistent_profile_path,
            }

        def close(self):
            captured["closed"] = True

    class _Fetcher:
        converter = object()
        scraper = object()
        session = object()

        def __init__(self, **kwargs):
            captured["fetcher_pool"] = kwargs["pool"]

    profile_dir = tmp_path / "gvf-browser-profile"
    output_dir = tmp_path / "out"

    monkeypatch.setenv("ENTREZ_EMAIL", "curator@example.org")
    monkeypatch.setattr(fetch_paywalled, "load_chrome_cookies", lambda **_kwargs: [])
    monkeypatch.setattr(fetch_paywalled, "cookie_domain_summary", lambda _cookies: {})
    monkeypatch.setattr(fetch_paywalled, "AuthenticatedBrowserPool", _Pool)
    monkeypatch.setattr(fetch_paywalled, "BrowserHTMLFetcher", _Fetcher)
    monkeypatch.setattr(
        fetch_paywalled, "hydrate_session_with_browser_cookies", lambda *_args: 0
    )
    monkeypatch.setattr(
        fetch_paywalled, "pubmed_resolve_doi", lambda *_args: "10.1002/example"
    )
    monkeypatch.setattr(
        fetch_paywalled,
        "prime_authenticated_browser",
        lambda pool, urls, timeout_s: captured.update(
            {"auth_urls": urls, "auth_timeout_s": timeout_s}
        ),
    )
    monkeypatch.setattr(
        fetch_paywalled,
        "fetch_one",
        lambda *_args, **_kwargs: {
            "pmid": "22222",
            "strategy": "wiley",
            "outcome": "success",
            "chars": 12000,
            "supp_files": 1,
            "supplements_downloaded": 1,
            "final_url": "https://onlinelibrary.wiley.com/doi/full/example",
            "reason": "ok",
        },
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "fetch_paywalled",
            "--pmid",
            "22222",
            "--output",
            str(output_dir),
            "--no-cookies",
            "--no-headless",
            "--browser-profile-dir",
            str(profile_dir),
            "--auth-url",
            "https://onlinelibrary.wiley.com",
            "--auth-timeout-s",
            "45",
        ],
    )

    assert fetch_paywalled.main() == 0

    assert captured["pool"]["headless"] is False
    assert captured["pool"]["persistent_profile_path"] == str(profile_dir.resolve())
    assert captured["auth_urls"] == ["https://onlinelibrary.wiley.com"]
    assert captured["auth_timeout_s"] == 45
    assert captured["closed"] is True


def test_fetch_one_preserves_failed_publisher_supplements_for_pmc_fallback(
    tmp_path, monkeypatch
):
    canonical = tmp_path / "12345_FULL_CONTEXT.md"
    per_pmid = tmp_path / "12345" / "FULL_CONTEXT.md"
    per_pmid.parent.mkdir()
    canonical.write_text("# MAIN TEXT\n\napi body\n", encoding="utf-8")
    per_pmid.write_text(canonical.read_text(encoding="utf-8"), encoding="utf-8")

    result = SimpleNamespace(
        main_markdown="stub",
        main_html="<html></html>",
        supp_files=[
            {"url": "https://publisher.example/supp.docx", "name": "supp.docx"}
        ],
        figure_paths=[],
        notes=[],
        error=None,
        publisher="elsevier_open",
        final_url="https://publisher.example/article",
    )
    fetcher = SimpleNamespace(
        fetch=lambda **_kwargs: result,
        converter=object(),
        scraper=object(),
        session=object(),
    )
    pmc_calls = []

    monkeypatch.setattr(
        fetch_paywalled, "find_strategy", lambda **_kwargs: SimpleNamespace(NAME="test")
    )
    monkeypatch.setattr(
        fetch_paywalled,
        "write_outputs",
        lambda *_args, **_kwargs: (
            canonical,
            per_pmid,
            "stub",
            EnrichmentResult(unified_markdown="stub", supplement_results=[]),
        ),
    )
    monkeypatch.setattr(
        fetch_paywalled,
        "validate_article_content",
        lambda _body: (False, "stub body"),
    )
    monkeypatch.setattr(
        fetch_paywalled,
        "try_publisher_api_fallback",
        lambda *_args, **_kwargs: {
            "strategy": "elsevier_api",
            "outcome": "success_via_elsevier_api",
            "api_label": "Elsevier Article Retrieval API",
            "chars": 9000,
            "reason": "ok",
            "path": str(canonical),
            "canonical_path": str(canonical),
            "per_pmid_path": str(per_pmid),
            "supp_files": 0,
            "supplements_downloaded": 0,
            "final_url": "https://api.example/article",
        },
    )

    def fake_pmc(*_args, **_kwargs):
        pmc_calls.append(True)
        return None

    monkeypatch.setattr(fetch_paywalled, "try_pmc_fallback", fake_pmc)

    row = fetch_paywalled.fetch_one(
        fetcher,
        "12345",
        "10.1016/j.example.2020.01.001",
        tmp_path,
        pmc_session=object(),
    )

    assert row["outcome"] == "success_via_elsevier_api"
    assert row["supp_files"] == 1
    assert pmc_calls == [True]
