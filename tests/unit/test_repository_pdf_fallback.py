"""Recovery must find alternate copies without accepting a different paper."""

import json
from types import SimpleNamespace
from unittest.mock import Mock

import fitz
import pytest
import requests

from harvesting.repository_pdf_fallback import (
    RepositoryPDFRecovery,
    RepositoryPDFResult,
    article_identity,
    write_repository_source,
)
from pipeline.source_quality import (
    is_reusable_fulltext_source,
    is_usable_fulltext_source,
)

TITLE = "Clinical characterization of inherited cardiac channel variants"
DOI = "10.1234/cardiac.456"
PDF_URL = "https://repository.example.edu/record/123/document"


def pdf_bytes(title=TITLE, doi=DOI):
    doc = fitz.open()
    for section in ["Abstract", "Methods", "Results", "Discussion"]:
        page = doc.new_page(width=900, height=1100)
        lines = [title, f"doi: {doi}", section]
        for i in range(15):
            lines.append(
                f"Observation {i}: Study participants underwent clinical examination and "
                "electrocardiographic characterization. The investigators recorded patient "
                "genotypes, symptoms and uncertain findings separately for each family. "
                "Variant carriers were assessed during follow-up using the stated methods."
            )
        assert (
            page.insert_textbox(
                fitz.Rect(30, 30, 870, 1070), "\n\n".join(lines), fontsize=11
            )
            >= 0
        )
    body = doc.tobytes()
    doc.close()
    return body


class Response:
    def __init__(self, url, *, data=None, body=b"", status=200, headers=None):
        self.url, self.data, self.content = url, data, body
        self.status_code = status
        self.headers = headers or {}
        self.is_redirect = status in (301, 302, 307, 308)

    def json(self):
        return self.data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, _):
        yield self.content

    def close(self):
        pass


class Session:
    def __init__(self, responses):
        self.responses = responses
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        assert url in self.responses, f"unexpected network request: {url}"
        return self.responses[url]


def recovery(tmp_path, urls, responses, **extra):
    oa = {"doi": DOI, "title": TITLE, "oa_locations": urls}
    session = Session(responses)
    client = RepositoryPDFRecovery(
        email="test@example.org",
        session=session,
        unpaywall=SimpleNamespace(find_open_access=lambda doi: (oa, None)),
    )
    # Unrelated providers have no hits in these isolated download tests.
    client._openalex_candidates = lambda r: []
    client._hal_candidates = lambda r: []
    result = client.recover(
        pmid="12345678", doi=DOI, title=TITLE, output_dir=tmp_path, **extra
    )
    return result, session


def test_failed_best_copy_continues_to_repository_pdf_and_retains_source(tmp_path):
    bad = "https://publisher.example.org/article.pdf"
    result, session = recovery(
        tmp_path,
        [
            {"url_for_pdf": bad, "host_type": "repository"},
            {"url_for_pdf": PDF_URL, "host_type": "repository"},
        ],
        {bad: Response(bad, status=403), PDF_URL: Response(PDF_URL, body=pdf_bytes())},
    )
    assert result.success
    assert [x[0] for x in session.calls] == [bad, PDF_URL]
    assert result.pdf_sha256
    assert "<!-- PDF page 2 -->" in result.markdown
    assert result.candidate["url"] == PDF_URL
    path, content = write_repository_source(result, tmp_path, gene="TEST")
    assert path.read_text() == content
    assert is_usable_fulltext_source(path)
    assert not is_reusable_fulltext_source(path)
    artifact = json.loads((tmp_path / "12345678_artifacts.json").read_text())
    assert artifact["main_text"]["pdf_sha256"] == result.pdf_sha256
    assert artifact["supplement_surface_status"] == "unavailable"


def test_metadata_only_repository_landing_page_can_supply_public_pdf(tmp_path):
    landing = "https://repository.example.edu/record/123"
    client = RepositoryPDFRecovery(email="test@example.org")
    candidates = client._locations(
        [
            {
                "landing_page_url": landing,
                "pdf_url": None,
                "is_oa": False,
                "source": {"type": "repository"},
            }
        ],
        "openalex",
        "https://api.openalex.org/works/example",
    )
    session = Session(
        {
            landing: Response(
                landing,
                body=b'<html><meta name="citation_pdf_url" content="/record/123/document"></html>',
            ),
            PDF_URL: Response(PDF_URL, body=pdf_bytes()),
        }
    )
    client.session = session
    client._unpaywall_candidates = lambda r: []
    client._openalex_candidates = lambda r: candidates
    client._hal_candidates = lambda r: []
    result = client.recover(pmid="12345678", doi=DOI, title=TITLE, output_dir=tmp_path)
    assert result.success
    assert result.candidate["discovered_on"] == landing
    assert result.source_url == PDF_URL


@pytest.mark.parametrize(
    "body",
    [
        b"<html>Sign in to access this PDF</html>",
        pdf_bytes(title="A different article about the cardiac literature"),
        pdf_bytes(doi="10.1234/different.999"),
    ],
)
def test_wrong_payload_or_wrong_article_is_rejected_without_context_write(
    tmp_path, body
):
    result, _ = recovery(
        tmp_path, [{"url_for_pdf": PDF_URL}], {PDF_URL: Response(PDF_URL, body=body)}
    )
    assert not result.success
    assert not list(tmp_path.glob("*_FULL_CONTEXT.md"))
    assert json.loads((tmp_path / "12345678_repository_recovery.json").read_text())[
        "attempts"
    ]


def test_citing_paper_and_hal_cover_are_not_identity_witnesses():
    cover = f"To cite this version\n{TITLE}\ndoi: {DOI}\nHAL Id: hal-123"
    other = (
        f"Another cardiac study\nAbstract\nClinical material\nReferences\n{TITLE} {DOI}"
    )
    assert not article_identity([cover, other], TITLE, DOI)[0]
    assert article_identity(
        [cover, f"{TITLE}\ndoi: {DOI}\nAbstract\nClinical material"], TITLE, DOI
    )[0]


def test_background_in_title_does_not_end_header():
    title = "Genetic background of inherited cardiac variants in families"
    assert article_identity([f"{title}\nAbstract\nThe article"], title, DOI)[0]


@pytest.mark.parametrize(
    "cover",
    [
        f"Amsterdam UMC\n{TITLE}\nDOI: {DOI}\nDocument Version\nLink to publication",
        f"This is the peer reviewed version of the following article:\n{TITLE}\n{DOI}",
    ],
)
def test_repository_citation_cover_does_not_validate_wrong_attachment(cover):
    assert not article_identity(
        [cover, "Another article\nAbstract: a different study"], TITLE, DOI
    )[0]
    assert article_identity(
        [cover, f"{TITLE}\nAbstract: Study participants"], TITLE, DOI
    )[0]


def test_title_quoted_after_inline_abstract_heading_is_not_article_header():
    assert not article_identity(
        [f"Another paper\nAbstract: A review of {TITLE}"], TITLE, DOI
    )[0]


def test_exact_doi_metadata_can_omit_specific_subtitle():
    client = RepositoryPDFRecovery(email="test@example.org")
    result = RepositoryPDFResult(
        doi=DOI, title=TITLE + ": a comprehensive clinical study"
    )
    assert client._bound({"doi": DOI, "title": TITLE}, result)
    assert article_identity(
        [f"{TITLE}\nAbstract\nClinical observations"], result.title, DOI
    )[0]
    assert not client._bound({"doi": "10.1234/different", "title": TITLE}, result)
    short = RepositoryPDFResult(doi=DOI, title="Genetics: a detailed clinical study")
    assert not client._bound({"doi": DOI, "title": "Genetics"}, short)


def test_conflicting_index_doi_is_rejected_before_download(tmp_path):
    client = RepositoryPDFRecovery(
        email="test@example.org",
        session=Session({}),
        unpaywall=SimpleNamespace(
            find_open_access=lambda doi: (
                {
                    "doi": "10.1234/wrong",
                    "title": TITLE,
                    "oa_locations": [{"url_for_pdf": PDF_URL}],
                },
                None,
            )
        ),
    )
    client._openalex_candidates = lambda r: []
    client._hal_candidates = lambda r: []
    result = client.recover(pmid="12345678", doi=DOI, title=TITLE, output_dir=tmp_path)
    assert not result.success
    assert not client.session.calls


def test_missing_doi_resolved_by_exact_pmid_metadata(tmp_path):
    epmc = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    session = Session(
        {
            epmc: Response(
                epmc,
                data={
                    "resultList": {
                        "result": [
                            {
                                "id": "12345678",
                                "source": "MED",
                                "doi": DOI,
                                "title": TITLE,
                            }
                        ]
                    }
                },
            ),
            PDF_URL: Response(PDF_URL, body=pdf_bytes()),
        }
    )
    client = RepositoryPDFRecovery(
        email="test@example.org",
        session=session,
        unpaywall=SimpleNamespace(
            find_open_access=lambda doi: (
                {
                    "doi": DOI,
                    "title": TITLE,
                    "oa_locations": [{"url_for_pdf": PDF_URL}],
                },
                None,
            )
        ),
    )
    result = client.recover(pmid="12345678", output_dir=tmp_path)
    assert result.success and result.doi == DOI


def test_redirect_to_local_address_is_not_requested(tmp_path):
    result, session = recovery(
        tmp_path,
        [{"url_for_pdf": PDF_URL}],
        {
            PDF_URL: Response(
                PDF_URL,
                status=302,
                headers={"Location": "http://127.0.0.1/private.pdf"},
            )
        },
    )
    assert not result.success
    assert len(session.calls) == 1


def test_normal_harvester_uses_fallback_after_publisher_failure_with_browser_disabled(
    tmp_path, monkeypatch
):
    from harvesting.orchestrator import PMCHarvester

    harvester = PMCHarvester.__new__(PMCHarvester)
    harvester.browser_html = None
    harvester._download_free_text_pmid = Mock(
        return_value=(False, "publisher HTTP 403", None)
    )
    harvester._download_repository_pmid = Mock(
        return_value=(True, str(tmp_path / "context.md"), "recovered PDF page text")
    )
    assert (
        harvester._download_with_tier35_fallback("12345678", DOI)[2]
        == "recovered PDF page text"
    )
    harvester._download_repository_pmid.assert_called_once_with("12345678", DOI)
    harvester._download_repository_pmid.reset_mock()
    harvester._download_free_text_pmid.return_value = (
        True,
        "publisher.md",
        "good publisher body",
    )
    assert (
        harvester._download_with_tier35_fallback("12345678", DOI)[2]
        == "good publisher body"
    )
    harvester._download_repository_pmid.assert_not_called()


def test_paywall_entrypoint_uses_same_repository_fallback(tmp_path, monkeypatch):
    from scripts import fetch_paywalled as fp

    monkeypatch.setattr(fp, "find_strategy", lambda **kw: None)
    monkeypatch.setattr(fp, "try_publisher_api_fallback", lambda *a, **kw: None)
    repo = Mock(
        return_value={
            "outcome": "success_via_repository_body_only",
            "chars": 10000,
            "reason": "verified",
        }
    )
    monkeypatch.setattr(fp, "try_repository_pdf_fallback", repo)
    scholar = Mock(
        side_effect=AssertionError("must not continue after repository success")
    )
    monkeypatch.setattr(fp, "try_scholar_pdf_fallback", scholar)
    row = fp.fetch_one(
        SimpleNamespace(converter=None),
        "12345678",
        DOI,
        tmp_path,
        pmc_session=object(),
        title=TITLE,
    )
    assert row["outcome"] == "success_via_repository_body_only"
    assert repo.call_args.kwargs["title"] == TITLE


def test_invalid_pmc_body_tries_repository_and_carries_supplement_text(tmp_path):
    from harvesting.orchestrator import PMCHarvester
    from harvesting.figure_extractor import CaptionExtractionResult

    harvester = PMCHarvester.__new__(PMCHarvester)
    harvester.output_dir, harvester.gene_symbol = tmp_path, "TEST"
    harvester.pmc_api = SimpleNamespace(
        get_doi_from_pmid=lambda pmid: DOI,
        pmid_to_pmcid=lambda pmid: "PMC123",
        get_fulltext_xml=lambda pmcid: "<article><body>Access denied</body></article>",
    )
    harvester.converter = SimpleNamespace(xml_to_markdown=lambda xml: "Access denied")
    harvester._extract_pmc_figures = Mock(return_value=([], CaptionExtractionResult()))
    harvester._markdown_missing_body = lambda body: False
    harvester._process_supplements = Mock(
        return_value=("Supplement observations", 1, [])
    )
    harvester._record_supplement_results = Mock()
    harvester._build_unified_content = Mock(return_value="Access denied")
    harvester._log_paywalled = Mock()
    harvester._download_repository_pmid = Mock(
        return_value=(True, "repository.md", "Body and supplement observations")
    )
    harvester._build_abstract_only_fallback = Mock(
        side_effect=AssertionError("repository succeeded")
    )
    assert harvester.download_pmid("12345678")[0]
    harvester._download_repository_pmid.assert_called_once_with(
        "12345678", DOI, supplement_markdown="Supplement observations"
    )


def test_persisting_recovered_body_preserves_existing_supplement_metadata(tmp_path):
    artifact = tmp_path / "12345678_artifacts.json"
    artifact.write_text(
        json.dumps(
            {
                "supplements": [{"path": "observations.xlsx"}],
                "main_text": {"figure_captions": 2, "table_captions": 3},
            }
        )
    )
    result = RepositoryPDFResult(
        success=True, pmid="12345678", doi=DOI, title=TITLE, markdown="Body"
    )
    _, text = write_repository_source(
        result, tmp_path, supplement_markdown="\n## Supplement\nObservations"
    )
    assert "Observations" in text
    assert json.loads(artifact.read_text())["supplements"] == [
        {"path": "observations.xlsx"}
    ]
    assert json.loads(artifact.read_text())["main_text"]["figure_captions"] == 2
    assert json.loads(artifact.read_text())["main_text"]["table_captions"] == 3


def test_header_doi_with_glued_publisher_word_matches():
    assert article_identity(
        [f"{TITLE}\ndoi:{DOI}Received\nAbstract\nStudy"], TITLE, DOI
    )[0]


def test_enabled_browser_gets_completion_chance_before_repository(monkeypatch):
    from harvesting.orchestrator import PMCHarvester
    import harvesting.browser_html as browser

    monkeypatch.setattr(browser, "get_pub_date_from_pmid", lambda *args: None)
    harvester = PMCHarvester.__new__(PMCHarvester)
    harvester.pmc_api = Mock()
    harvester._download_free_text_pmid = Mock(
        return_value=(False, "publisher failed", None)
    )
    harvester.browser_html = Mock()
    harvester.browser_html.is_enabled.return_value = True
    harvester.browser_html.fetch.return_value.is_usable.return_value = True
    harvester._finalize_browser_html = Mock(
        return_value=(True, "publisher.md", "Body plus supplements")
    )
    harvester._download_repository_pmid = Mock(
        side_effect=AssertionError("browser succeeded")
    )
    assert (
        harvester._download_with_tier35_fallback("12345678", DOI)[2]
        == "Body plus supplements"
    )
    harvester._download_repository_pmid.assert_not_called()


@pytest.mark.parametrize(
    "headers", [{}, {"Location": "https://repository.example.edu/broken"}]
)
def test_malformed_or_looping_redirect_keeps_next_candidate(tmp_path, headers):
    broken = "https://repository.example.edu/broken"
    result, session = recovery(
        tmp_path,
        [
            {"url_for_pdf": broken, "host_type": "repository"},
            {"url_for_pdf": PDF_URL, "host_type": "repository"},
        ],
        {
            broken: Response(broken, status=302, headers=headers),
            PDF_URL: Response(PDF_URL, body=pdf_bytes()),
        },
    )
    assert result.success
    assert len(session.calls) == 2
