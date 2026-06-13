"""Tests for ``harvesting.paywall_context_enrichment``.

These cover the high-leverage paywall-recovery improvement that landed with
the v2.1.x recall push: when ``scripts/fetch_paywalled.py`` rescues a
paywalled paper, the rescued FULL_CONTEXT.md now also carries the article's
figure / table captions and any supplement files we could download.
Variants frequently live in those regions and the previous body-only
recovery silently dropped them.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict
from unittest.mock import MagicMock

import pytest

from harvesting.figure_extractor import (
    CaptionExtractionResult,
    FigureCaption,
    SupplementDescription,
    TableCaption,
)
from harvesting.paywall_context_enrichment import (
    EnrichmentResult,
    enrich_paywall_full_context,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


SAMPLE_HTML = """
<html>
  <body>
    <article>
      <h1>Family-screening study</h1>
      <p>We screened 86 probands and identified KCNH2 variants.</p>
      <figure id="fig1">
        <img src="https://example.org/fig1.png"/>
        <figcaption>
          <strong>Figure 1.</strong> Pedigree of family carrying KCNH2 G604S.
          Filled symbols are affected; open symbols are unaffected.
        </figcaption>
      </figure>
      <table id="t1">
        <caption>Table 1. Variant carrier counts.</caption>
        <tr><th>Variant</th><th>Carriers</th></tr>
        <tr><td>p.Arg176Trp</td><td>4</td></tr>
      </table>
      <section id="supplementary-materials">
        <p><strong>Supplementary Table S1.</strong>
           Extended KCNH2 variant list — 86 missense variants.
        </p>
      </section>
    </article>
  </body>
</html>
"""


@pytest.fixture
def converter() -> MagicMock:
    """A converter stub that records what it would convert.

    We don't want the test depending on markitdown/pdfplumber/etc.; the
    SupplementProcessingService dispatches on the file extension, so we
    just need the matching method to return canned markdown.
    """
    c = MagicMock()
    c.excel_to_markdown.return_value = (
        "| variant | n |\n|---|---|\n| p.Ala561Val | 1 |\n\n"
    )
    c.html_supplement_to_markdown.return_value = "Extended variant list HTML.\n\n"
    c.pdf_to_markdown.return_value = "PDF body text.\n\n"
    return c


@pytest.fixture
def fake_session() -> SimpleNamespace:
    """A session double whose ``get`` writes canned bytes to disk.

    Mirrors ``requests.Session.get(stream=True)``: returns a context manager
    yielding an object that exposes ``status_code`` and ``iter_content``.
    """
    written: Dict[str, bytes] = {}

    class _Resp:
        def __init__(self, content: bytes, status_code: int = 200):
            self.content = content
            self.status_code = status_code

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def iter_content(self, chunk_size: int = 1024):
            yield self.content

    def get(url: str, stream: bool = False, timeout: int = 60):
        # Pretend any URL containing ".xlsx" returns a small xlsx payload.
        if url.endswith(".xlsx"):
            payload = b"PK\x03\x04dummy-xlsx-content"
        elif url.endswith(".html"):
            payload = b"<html><body><table><tr><td>p.Ala561Val</td></tr></table></body></html>"
        else:
            return _Resp(b"", status_code=404)
        written[url] = payload
        return _Resp(payload)

    return SimpleNamespace(get=get, _written=written)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_enrich_appends_captions_to_body(tmp_path: Path, converter: MagicMock):
    """Caption block lands after body when HTML has captions."""
    body = "# MAIN TEXT\n\nWe screened 86 probands. p.Arg176Trp identified.\n"

    result = enrich_paywall_full_context(
        body_markdown=body,
        html=SAMPLE_HTML,
        supp_files=None,
        pmid="99999999",
        output_dir=tmp_path,
        converter=converter,
        session=None,
        download_supplements=False,
    )

    assert isinstance(result, EnrichmentResult)
    md = result.unified_markdown
    assert md.startswith(body)
    assert "## FIGURE CAPTIONS" in md
    assert "Pedigree of family carrying KCNH2 G604S" in md
    assert "## TABLE CAPTIONS" in md
    assert "Variant carrier counts" in md
    assert result.figure_caption_count >= 1
    assert result.table_caption_count >= 1
    # Body content is preserved verbatim — no rewriting of p.Arg176Trp.
    assert "p.Arg176Trp" in md
    # Captions sidecar is written next to the figures dir.
    expected_index = tmp_path / "99999999_figures" / "captions_index.json"
    assert expected_index.exists()


def test_enrich_with_no_html_still_returns_body(tmp_path: Path, converter: MagicMock):
    """Empty HTML must not crash and must preserve the rescued body intact."""
    body = "# MAIN TEXT\n\nBody only.\n"

    result = enrich_paywall_full_context(
        body_markdown=body,
        html=None,
        supp_files=None,
        pmid="00000001",
        output_dir=tmp_path,
        converter=converter,
        session=None,
    )

    assert result.unified_markdown == body
    assert result.figure_caption_count == 0
    assert result.table_caption_count == 0
    assert result.supplement_count == 0
    # No captions => no sidecar.
    assert (
        not (tmp_path / "00000001_figures").exists()
        or not (tmp_path / "00000001_figures" / "captions_index.json").exists()
    )


def test_enrich_downloads_and_converts_supplements(
    tmp_path: Path, converter: MagicMock, fake_session: SimpleNamespace
):
    """Convertible supplement links are fetched and appended after captions."""
    body = "# MAIN TEXT\n\nMain body.\n"
    supp_files = [
        {
            "url": "https://publisher.example/suppl/table-s1.xlsx",
            "name": "table-s1.xlsx",
        },
        # No URL => skipped silently.
        {"url": "", "name": "broken.pdf"},
        # Unknown extension => skipped (we don't speak .mp4).
        {
            "url": "https://publisher.example/suppl/movie.mp4",
            "name": "movie.mp4",
        },
    ]

    result = enrich_paywall_full_context(
        body_markdown=body,
        html=SAMPLE_HTML,
        supp_files=supp_files,
        pmid="12345678",
        output_dir=tmp_path,
        converter=converter,
        session=fake_session,
    )

    # Excel converter got called for the xlsx supplement.
    converter.excel_to_markdown.assert_called_once()
    assert "p.Ala561Val" in result.unified_markdown
    assert "# SUPPLEMENTAL FILE 1: table-s1.xlsx" in result.unified_markdown
    assert result.supplement_count == 1
    # Captions still appear, ordered before supplement markdown.
    assert result.unified_markdown.index(
        "## FIGURE CAPTIONS"
    ) < result.unified_markdown.index("# SUPPLEMENTAL FILE 1")
    # Downloaded file lives in the per-PMID supplements directory.
    assert (tmp_path / "12345678_supplements" / "table-s1.xlsx").exists()
    # Skipped entries didn't leak through.
    assert not (tmp_path / "12345678_supplements" / "movie.mp4").exists()


def test_enrich_uses_supplement_download_fallback_after_requests_failure(
    tmp_path: Path, converter: MagicMock
):
    body = "# MAIN TEXT\n\nMain body.\n"

    class _Resp:
        status_code = 403

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def iter_content(self, chunk_size: int = 1024):
            yield b""

    class _FailingSession:
        def get(self, *_args, **_kwargs):
            return _Resp()

    fallback_calls = []

    def fallback(url: str, path: Path, pmid: str, filename: str, supp: Dict[str, Any]):
        fallback_calls.append((url, pmid, filename, supp.get("source_url")))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"PK\x03\x04fallback-xlsx")
        return True

    result = enrich_paywall_full_context(
        body_markdown=body,
        html="",
        supp_files=[
            {
                "url": "https://publisher.example/cms/attachment/mmc1.xlsx",
                "name": "mmc1.xlsx",
            }
        ],
        pmid="24680",
        output_dir=tmp_path,
        converter=converter,
        session=_FailingSession(),
        source_url="https://publisher.example/article",
        supplement_download_fallback=fallback,
    )

    assert fallback_calls == [
        (
            "https://publisher.example/cms/attachment/mmc1.xlsx",
            "24680",
            "mmc1.xlsx",
            "https://publisher.example/article",
        )
    ]
    assert result.supplement_count == 1
    assert "p.Ala561Val" in result.unified_markdown
    assert (tmp_path / "24680_supplements" / "mmc1.xlsx").exists()


def test_enrich_downloads_caption_table_images_for_injected_vision(
    tmp_path: Path, converter: MagicMock
):
    """Image-only PMC table blocks are downloaded for the optional vision pass."""
    html = """
    <html><body>
      <article>
        <section class="tw xbox" id="t2-6124">
          <div class="caption p">
            <p><strong>Table 2</strong> Mutations in patients with LQT2.</p>
          </div>
          <p class="img-box">
            <img class="graphic" src="/pmc/blobs/2677528/T2-6124.jpg"/>
          </p>
        </section>
      </article>
    </body></html>
    """

    class _Resp:
        status_code = 200

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def iter_content(self, chunk_size: int = 1024):
            yield b"fake-jpeg-bytes"

    calls: list[str] = []

    class _Session:
        def get(self, url, stream=False, timeout=60):
            calls.append(url)
            return _Resp()

    captured_paths: list[Path] = []

    def _extractor(paths: list[Path]) -> str:
        captured_paths.extend(paths)
        return "\n\n## FIGURE IMAGE TEXT\n\nR176W A193fsX T613M\n"

    result = enrich_paywall_full_context(
        body_markdown="# MAIN TEXT\n\nBody.\n",
        html=html,
        supp_files=None,
        pmid="19038855",
        output_dir=tmp_path,
        converter=converter,
        session=_Session(),
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2677528/",
        download_supplements=False,
        image_text_extractor=_extractor,
    )

    assert calls == ["https://pmc.ncbi.nlm.nih.gov/pmc/blobs/2677528/T2-6124.jpg"]
    assert captured_paths
    assert captured_paths[0].exists()
    assert captured_paths[0].parent == tmp_path / "19038855_figures" / "html_images"
    assert "R176W A193fsX T613M" in result.unified_markdown


def test_oversized_supplement_is_dropped(
    tmp_path: Path, converter: MagicMock, monkeypatch: pytest.MonkeyPatch
):
    """Supplement downloads above the size cap should never land on disk."""

    class _BigResp:
        status_code = 200

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def iter_content(self, chunk_size: int = 1024):
            # 4 chunks of 8 bytes each = 32 bytes > tiny limit below.
            for _ in range(4):
                yield b"x" * 8

    class _Session:
        def get(self, url, stream=False, timeout=60):
            return _BigResp()

    body = "# MAIN TEXT\n\nBody.\n"
    result = enrich_paywall_full_context(
        body_markdown=body,
        html=None,
        supp_files=[
            {"url": "https://example.org/big.pdf", "name": "big.pdf"},
        ],
        pmid="55555",
        output_dir=tmp_path,
        converter=converter,
        session=_Session(),
        supplement_size_limit_bytes=16,
    )

    assert result.supplement_count == 0
    # The partially-written file must have been cleaned up.
    assert not (tmp_path / "55555_supplements" / "big.pdf").exists()
    # Body is still preserved.
    assert result.unified_markdown.startswith(body)


def test_extra_captions_merge_with_html_captions(tmp_path: Path, converter: MagicMock):
    """JATS-side captions can be merged in alongside HTML-extracted ones."""
    extra = CaptionExtractionResult(
        figures=[
            FigureCaption(
                label="Figure 2",
                title="Additional figure from JATS",
                text="Captured from XML.",
            )
        ],
        tables=[],
        supplements=[
            SupplementDescription(
                label="Supplementary Table S2",
                title="JATS-supplied supplement",
                text="Describes 12 additional variants.",
            )
        ],
    )

    result = enrich_paywall_full_context(
        body_markdown="# MAIN TEXT\n\nBody.\n",
        html=SAMPLE_HTML,
        supp_files=None,
        pmid="42",
        output_dir=tmp_path,
        converter=converter,
        session=None,
        extra_captions=extra,
        download_supplements=False,
    )

    md = result.unified_markdown
    # Both HTML-derived and JATS-derived figures are present.
    assert "Pedigree of family carrying KCNH2 G604S" in md
    assert "Additional figure from JATS" in md
    # Supplement description from JATS shows up under SUPPLEMENT DESCRIPTIONS.
    assert "## SUPPLEMENT DESCRIPTIONS" in md
    assert "JATS-supplied supplement" in md


def test_process_supplement_files_receives_extract_figures_and_figures_dir(
    tmp_path: Path, converter: MagicMock, monkeypatch: pytest.MonkeyPatch
):
    """process_supplement_files must be called with extract_figures=True and the
    correct figures_dir so PDF/ZIP supplement images are extracted alongside captions."""
    import harvesting.paywall_context_enrichment as _mod

    calls: list = []

    def _fake_process_supplement_files(**kwargs):
        calls.append(kwargs)
        from harvesting.supplement_processing_service import SupplementProcessingResult

        return SupplementProcessingResult(
            supplement_markdown="",
            downloaded_count=1,
            total_figures_extracted=0,
            file_results=[],
        )

    monkeypatch.setattr(
        _mod, "process_supplement_files", _fake_process_supplement_files
    )

    enrich_paywall_full_context(
        body_markdown="# BODY\n",
        html=None,
        supp_files=[{"url": "https://example.org/s1.pdf", "name": "s1.pdf"}],
        pmid="77777777",
        output_dir=tmp_path,
        converter=converter,
        session=MagicMock(
            **{
                "get.return_value.__enter__": lambda s: s,
                "get.return_value.__exit__": lambda s, *a: False,
                "get.return_value.status_code": 200,
                "get.return_value.iter_content.return_value": [b"data"],
            }
        ),
    )

    assert len(calls) == 1
    assert calls[0]["extract_figures"] is True
    assert calls[0]["figures_dir"] == tmp_path / "77777777_figures"


def test_caption_block_skipped_when_no_captions(tmp_path: Path, converter: MagicMock):
    """A caption-less page must not introduce empty FIGURE CAPTIONS headers."""
    html = "<html><body><p>Body-only article with no figures.</p></body></html>"
    result = enrich_paywall_full_context(
        body_markdown="# MAIN TEXT\n\nBody.\n",
        html=html,
        supp_files=None,
        pmid="1",
        output_dir=tmp_path,
        converter=converter,
        session=None,
        download_supplements=False,
    )
    assert "## FIGURE CAPTIONS" not in result.unified_markdown
    assert "## TABLE CAPTIONS" not in result.unified_markdown
