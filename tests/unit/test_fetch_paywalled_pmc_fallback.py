"""Tests for ``scripts.fetch_paywalled.try_pmc_fallback``.

The PMC fallback runs when Tier 3.5 leaves us with a paywall stub. It must
not only rescue the body but also pull supplement links from the PMC HTML so
the enricher can download/convert convertible supplements — variant tables
frequently live in those supplement spreadsheets.

These tests mock the network completely; no real PMC traffic.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List
from unittest.mock import MagicMock

import pytest
import requests

from scripts import fetch_paywalled as fp


_PMC_HTML = """
<html><body>
  <article>
    <h1>NIH-deposited paper</h1>
    <p>We report KCNH2 variants in 86 probands. p.Arg176Trp identified.</p>
    <p>Cohort and methods paragraph with plenty of substantive body text so
    the quality gate doesn't reject this article. We include extensive
    detail on family screening, electrocardiographic findings, and
    genotype-phenotype correlations across the entire enrolled cohort. The
    variants span the full coding region of KCNH2 with multiple novel
    missense substitutions described in detail.</p>
    <figure>
      <figcaption><strong>Figure 1.</strong> Pedigree carrying KCNH2 G604S.</figcaption>
    </figure>
    <table><caption>Table 1. Variant counts.</caption><tr><td>x</td></tr></table>
    <a href="/articles/PMC123/bin/supp-table-s1.xlsx">Supplementary Table S1</a>
    <p>%s</p>
  </article>
</body></html>
""" % ("Filler paragraph to clear the 5KB byte floor for the PMC fetch guard. " * 200)


@pytest.fixture
def converter_stub() -> MagicMock:
    c = MagicMock()
    c.excel_to_markdown.return_value = (
        "| variant | n |\n|---|---|\n| p.Ala561Val | 1 |\n\n"
    )
    return c


def _make_session(html: str, supp_bytes: bytes = b"PK\x03\x04xlsx-bytes"):
    """Build a fake requests session covering Europe PMC + PMC HTML + supp."""

    class _Resp:
        def __init__(self, content: bytes, status_code: int = 200, text: str = ""):
            self.content = content
            self.status_code = status_code
            self.text = text or content.decode("utf-8", errors="ignore")
            self.ok = 200 <= status_code < 300

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def iter_content(self, chunk_size: int = 1024):
            yield self.content

        def json(self):
            import json as _json

            return _json.loads(self.text)

        def raise_for_status(self):
            if not self.ok:
                raise RuntimeError(f"status {self.status_code}")

    calls: List[str] = []

    def get(url: str, **kwargs: Any):
        calls.append(url)
        if "ebi.ac.uk/europepmc" in url or "europepmc.org" in url:
            return _Resp(
                b"",
                text='{"resultList":{"result":[{"pmcid":"PMC123"}]}}',
            )
        if "ncbi.nlm.nih.gov/pmc/articles/PMC123" in url and not url.endswith(".xlsx"):
            return _Resp(html.encode("utf-8"), text=html)
        if url.endswith(".xlsx"):
            return _Resp(supp_bytes)
        return _Resp(b"", status_code=404)

    return SimpleNamespace(get=get, _calls=calls)


def test_pmc_fallback_pulls_supplements_and_returns_counts(
    tmp_path: Path, converter_stub: MagicMock
):
    """try_pmc_fallback scrapes supplements and reports counts back."""
    session = _make_session(_PMC_HTML)

    row = fp.try_pmc_fallback(
        pmid="99999999",
        output_dir=tmp_path,
        session=session,
        converter=converter_stub,
    )

    assert row is not None, "expected PMC fallback to succeed on rich HTML"
    assert row["pmcid"] == "PMC123"
    assert row["outcome"] == "success"
    # At least one supplement link was scraped from the PMC HTML.
    assert row["supp_files"] >= 1
    # Enrichment downloaded the xlsx supplement and converted it.
    assert row["supplements_downloaded"] == 1
    converter_stub.excel_to_markdown.assert_called_once()
    # Captions made it through to the row.
    assert row["figure_captions"] >= 1
    assert row["table_captions"] >= 1
    # Unified FULL_CONTEXT.md contains the converted supplement payload.
    full_ctx = Path(row["path"]).read_text(encoding="utf-8")
    assert "p.Ala561Val" in full_ctx
    assert "## FIGURE CAPTIONS" in full_ctx
    # Supplement file was persisted under the per-PMID supplements dir.
    assert (tmp_path / "99999999_supplements" / "supp-table-s1.xlsx").exists()


def test_pmc_fallback_writes_flat_mirror_with_identical_content(
    tmp_path: Path, converter_stub: MagicMock
):
    """PMC fallback writes both per-PMID and flat-mirror FULL_CONTEXT.md.

    The downstream extractor (``cli.extract.find_input_files``) discovers
    rescued papers by globbing ``{output_dir}/{PMID}_FULL_CONTEXT.md``, so the
    flat mirror MUST exist and match the per-PMID copy byte-for-byte.
    """
    session = _make_session(_PMC_HTML)
    pmid = "99999999"

    row = fp.try_pmc_fallback(
        pmid=pmid,
        output_dir=tmp_path,
        session=session,
        converter=converter_stub,
    )

    assert row is not None
    per_pmid = tmp_path / pmid / "FULL_CONTEXT.md"
    flat = tmp_path / f"{pmid}_FULL_CONTEXT.md"
    assert per_pmid.exists(), "per-PMID FULL_CONTEXT.md missing"
    assert flat.exists(), "flat-mirror FULL_CONTEXT.md missing"
    assert per_pmid.read_text(encoding="utf-8") == flat.read_text(encoding="utf-8")
    # row.path should report the canonical (flat) location for discovery.
    assert row["path"] == str(flat)
    assert row["canonical_path"] == str(flat)
    assert row["per_pmid_path"] == str(per_pmid)


def test_pmc_fallback_reuses_passed_scraper(tmp_path: Path, converter_stub: MagicMock):
    """Caller-supplied scraper is used rather than instantiating a new one."""
    session = _make_session(_PMC_HTML)

    scraper = MagicMock()
    scraper.scrape_generic_supplements.return_value = [
        {
            "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC123/bin/supp-table-s1.xlsx",
            "name": "supp-table-s1.xlsx",
        }
    ]

    row = fp.try_pmc_fallback(
        pmid="12345678",
        output_dir=tmp_path,
        session=session,
        converter=converter_stub,
        scraper=scraper,
    )

    assert row is not None
    scraper.scrape_generic_supplements.assert_called_once()
    assert row["supp_files"] == 1


def test_write_outputs_writes_flat_mirror_with_identical_content(tmp_path: Path):
    """``write_outputs`` writes the per-PMID file AND a flat canonical mirror.

    The flat ``{output_dir}/{PMID}_FULL_CONTEXT.md`` is what the extraction
    discovery path globs for; this test pins down that the two files exist
    and have identical content, that the returned tuple shape is
    ``(canonical_path, per_pmid_path, body, enrichment)``, and that
    ``result.json`` records both locations under canonical/per-PMID keys.
    """
    pmid = "12345678"
    body = (
        "# Recovered Article\n\n"
        "This is the main body markdown produced by the DOM extractor.\n"
        "It contains the variant p.Arg176Trp and surrounding cohort context.\n"
    )
    result = SimpleNamespace(
        main_markdown=body,
        main_html="<html><body><article>x</article></body></html>",
        supp_files=[],
        publisher="aha",
        final_url="https://www.ahajournals.org/doi/10.1161/x",
        figure_paths=[],
        notes=[],
        error=None,
    )

    canonical_path, per_pmid_path, body_for_gate, enrichment = fp.write_outputs(
        pmid=pmid,
        result=result,
        output_dir=tmp_path,
        enrich=False,
    )

    flat = tmp_path / f"{pmid}_FULL_CONTEXT.md"
    per_pmid = tmp_path / pmid / "FULL_CONTEXT.md"
    assert flat.exists(), "flat-mirror FULL_CONTEXT.md missing"
    assert per_pmid.exists(), "per-PMID FULL_CONTEXT.md missing"
    assert flat.read_text(encoding="utf-8") == per_pmid.read_text(encoding="utf-8")
    assert flat.read_text(encoding="utf-8") == body
    # Returned tuple advertises the canonical path first.
    assert canonical_path == flat
    assert per_pmid_path == per_pmid
    assert body_for_gate == body
    assert enrichment is None  # enrich=False
    # result.json carries both locations under canonical/per-PMID metadata keys.
    meta = json.loads((tmp_path / pmid / "result.json").read_text(encoding="utf-8"))
    assert meta["canonical_full_context_path"] == str(flat)
    assert meta["per_pmid_full_context_path"] == str(per_pmid)


def test_fetch_one_promotes_pmc_when_publisher_supplements_fail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """A good publisher body is not enough when every supplement failed."""
    pmid = "19841300"
    body = "Full article body with KCNQ1 A300T and enough clinical context."
    canonical = tmp_path / f"{pmid}_FULL_CONTEXT.md"
    per_pmid = tmp_path / pmid / "FULL_CONTEXT.md"

    fetcher = SimpleNamespace(
        converter=MagicMock(),
        scraper=MagicMock(),
        session=SimpleNamespace(),
        fetch=MagicMock(
            return_value=SimpleNamespace(
                main_markdown=body,
                main_html="<html></html>",
                supp_files=[{"url": "https://publisher.example/supp.pdf"}],
                publisher="aha",
                final_url="https://publisher.example/article",
                figure_paths=[],
                notes=[],
                error=None,
            )
        ),
    )

    monkeypatch.setattr(
        fp,
        "find_strategy",
        lambda doi, allowlist=None: SimpleNamespace(NAME="aha"),
    )
    monkeypatch.setattr(fp, "validate_article_content", lambda text: (True, "ok"))

    def fake_write_outputs(*args, **kwargs):
        canonical.parent.mkdir(parents=True, exist_ok=True)
        per_pmid.parent.mkdir(parents=True, exist_ok=True)
        canonical.write_text(body, encoding="utf-8")
        per_pmid.write_text(body, encoding="utf-8")
        enrichment = SimpleNamespace(
            unified_markdown=body,
            supplement_count=0,
            figure_caption_count=0,
            table_caption_count=0,
        )
        return canonical, per_pmid, body, enrichment

    monkeypatch.setattr(fp, "write_outputs", fake_write_outputs)
    monkeypatch.setattr(
        fp,
        "try_pmc_fallback",
        lambda *args, **kwargs: {
            "pmid": pmid,
            "pmcid": "PMC3025752",
            "outcome": "success",
            "reason": "ok",
            "chars": 75000,
            "final_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3025752/",
            "path": str(tmp_path / "pmc_FULL_CONTEXT.md"),
            "canonical_path": str(tmp_path / "pmc_FULL_CONTEXT.md"),
            "per_pmid_path": str(tmp_path / pmid / "PMC_FULL_CONTEXT.md"),
            "supp_files": 2,
            "figure_captions": 5,
            "table_captions": 3,
            "supplements_downloaded": 1,
        },
    )

    row = fp.fetch_one(
        fetcher=fetcher,
        pmid=pmid,
        doi="10.1161/CIRCULATIONAHA.109.863076",
        output_dir=tmp_path,
        pmc_session=SimpleNamespace(),
    )

    assert row["outcome"] == "success_via_pmc"
    assert row["strategy"] == "aha+pmc_fallback"
    assert row["pmcid"] == "PMC3025752"
    assert row["supplements_downloaded"] == 1


def test_fetch_one_promotes_elsevier_api_when_browser_returns_stub(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    pmid = "15840476"
    stub_body = "log in, subscribe or purchase for full access"
    api_body = (
        "# MAIN TEXT\n\n"
        "### Results\n\n"
        "Full article body with KCNH2 variant tables and enough clinical context.\n\n"
        "Table 3 Summary of KCNH2 LQT2-associated variants\n"
        "| No. | Exon | Nucleotide | Variant | Location | No. of patients |\n"
        "| --- | --- | --- | --- | --- | --- |\n"
        "| 27 | 6 | 1264 G>A | A422T | S1 | 1 |\n"
    )

    fetcher = SimpleNamespace(
        converter=MagicMock(),
        scraper=MagicMock(),
        session=SimpleNamespace(),
        fetch=MagicMock(
            return_value=SimpleNamespace(
                main_markdown=stub_body,
                main_html="<html></html>",
                supp_files=[],
                publisher="elsevier_open",
                final_url="https://www.heartrhythmjournal.com/article/x/abstract",
                figure_paths=[],
                notes=[],
                error=None,
            )
        ),
    )

    class FakeElsevierClient:
        def __init__(self, *args, **kwargs):
            self.is_available = True

        @staticmethod
        def is_elsevier_doi(doi: str) -> bool:
            return doi.startswith("10.1016/")

        def fetch_fulltext(self, doi=None, pii=None, url=None):
            return api_body, None

    monkeypatch.setattr(
        fp,
        "find_strategy",
        lambda doi, allowlist=None: SimpleNamespace(NAME="elsevier_open"),
    )
    monkeypatch.setattr(fp, "ElsevierAPIClient", FakeElsevierClient)
    monkeypatch.setattr(
        fp,
        "validate_article_content",
        lambda text: (
            (False, "paywall") if "log in, subscribe" in text else (True, "ok")
        ),
    )

    row = fp.fetch_one(
        fetcher=fetcher,
        pmid=pmid,
        doi="10.1016/j.hrthm.2005.01.020",
        output_dir=tmp_path,
        pmc_session=SimpleNamespace(),
    )

    assert row["outcome"] == "success_via_elsevier_api"
    assert row["strategy"] == "elsevier_open+elsevier_api"
    assert row["chars"] == len(api_body)
    assert row["outcome"] in fp.FETCH_SUCCESS_OUTCOMES
    assert "A422T" in (tmp_path / f"{pmid}_FULL_CONTEXT.md").read_text(encoding="utf-8")


def test_pmc_fallback_returns_none_when_no_pmcid(
    tmp_path: Path, converter_stub: MagicMock
):
    """No PMC deposit => fallback is silently skipped."""

    class _Resp:
        status_code = 200
        ok = True
        text = '{"resultList":{"result":[]}}'
        content = b""

        def json(self) -> Dict[str, Any]:
            return {"resultList": {"result": []}}

        def raise_for_status(self):
            return None

    session = SimpleNamespace(get=lambda *a, **k: _Resp())

    row = fp.try_pmc_fallback(
        pmid="0",
        output_dir=tmp_path,
        session=session,
        converter=converter_stub,
    )
    assert row is None


def test_hydrate_session_with_browser_cookies_authenticates_download_hosts():
    """Browser cookies must also be available to requests supplement fetches."""
    session = requests.Session()

    count = fp.hydrate_session_with_browser_cookies(
        session,
        [
            {
                "name": "PublisherSession",
                "value": "abc123",
                "domain": ".neurology.org",
                "path": "/",
            },
            {
                "name": "HostCookie",
                "value": "from-url",
                "url": "https://www.heartrhythmjournal.com/article/example",
                "path": "/",
            },
        ],
    )

    assert count == 2

    neurology_request = session.prepare_request(
        requests.Request(
            "GET",
            "https://www.neurology.org/doi/suppl/10.1212/x/suppl_file/johnson.pdf",
        )
    )
    assert "PublisherSession=abc123" in neurology_request.headers["Cookie"]

    heart_request = session.prepare_request(
        requests.Request(
            "GET",
            "https://www.heartrhythmjournal.com/pb/assets/raw/foo.doc",
        )
    )
    assert "HostCookie=from-url" in heart_request.headers["Cookie"]


def test_hydrate_session_with_browser_cookies_skips_unusable_records():
    """Malformed cookie records should not poison the requests cookie jar."""
    session = requests.Session()

    count = fp.hydrate_session_with_browser_cookies(
        session,
        [
            {"value": "missing-name", "domain": ".example.org"},
            {"name": "missing-value", "domain": ".example.org"},
            {"name": "missing-domain", "value": "x"},
            {"name": "ok", "value": "", "domain": ".example.org"},
        ],
    )

    assert count == 1
    prepared = session.prepare_request(
        requests.Request("GET", "https://www.example.org/supplement.pdf")
    )
    assert "ok=" in prepared.headers["Cookie"]
