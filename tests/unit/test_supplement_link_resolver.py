"""Tests for harvesting.supplement_link_resolver.

The resolver exists because recorded supplement links used to be write-only:
rendered into FULL_CONTEXT.md as ``_link_:`` lines and read by nothing. These
tests pin the two halves — recovering links from what is on disk, and resolving
them into the ``supp_files`` dicts the existing download path consumes.
"""

from __future__ import annotations

import textwrap

import pytest

from harvesting.supplement_link_resolver import (
    ZIP_MAGIC,
    SupplementLink,
    doi_from_href,
    europepmc_archive_job,
    links_from_artifacts,
    looks_like_supplement_file,
    parse_links_from_markdown,
    pmcid_from_artifacts,
    resolve_link,
    resolve_links,
    to_supplement_jobs,
)

# The exact shape figure_extractor._format_supplement_desc_md writes.
MARKDOWN_FIXTURE = textwrap.dedent("""\
    ## SUPPLEMENT DESCRIPTIONS

    ### Supplementary Appendix

    **Supplementary Appendix (nejm_imboden_2744sa1.pdf)**

    Download 32.38 KB

    _link_: /doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm_imboden_2744sa1.pdf

    ### Supplementary Table S2

    Extended list of 90 missense variants.

    _link_: mmc1.xls
    """)


# ---------------------------------------------------------------------------
# Recovering links from disk
# ---------------------------------------------------------------------------


def test_parse_links_from_markdown_recovers_hrefs_with_labels():
    links = parse_links_from_markdown(MARKDOWN_FIXTURE)
    assert [link.href for link in links] == [
        "/doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm_imboden_2744sa1.pdf",
        "mmc1.xls",
    ]
    assert links[0].label == "Supplementary Appendix"
    assert links[1].label == "Supplementary Table S2"
    assert links[0].source == "full_context_markdown"


def test_parse_links_from_markdown_dedupes_repeated_hrefs():
    text = "### A\n\n_link_: mmc1.xls\n\n### B\n\n_link_: mmc1.xls\n"
    assert len(parse_links_from_markdown(text)) == 1


def test_parse_links_from_markdown_handles_no_links():
    assert parse_links_from_markdown("") == []
    assert parse_links_from_markdown("## FIGURE CAPTIONS\n\n### Figure 1\n") == []


def test_link_line_must_be_its_own_line():
    """An inline mention must not be mistaken for a recorded link."""
    assert parse_links_from_markdown("see _link_: mmc1.xls in the text") == []


def test_links_from_artifacts_reads_structured_entries():
    manifest = {
        "pmcid": "PMC4535901",
        "supplement_links": [
            {"label": "Table S1", "href": "mmc1.xls", "source": "pmc"},
            {"label": "no href", "href": ""},
        ],
    }
    links = links_from_artifacts(manifest)
    assert [link.href for link in links] == ["mmc1.xls"]
    assert links[0].source == "pmc"
    assert pmcid_from_artifacts(manifest) == "PMC4535901"


def test_links_from_artifacts_tolerates_older_manifests():
    assert links_from_artifacts({"pmid": "123"}) == []
    assert links_from_artifacts({}) == []
    assert pmcid_from_artifacts({"pmcid": None}) is None


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "href",
    [
        "mmc1.xls",
        "1471-2407-11-265-S3.XLS",  # PMC serves uppercase extensions
        "NIHMS475442-supplement-1.pdf",
        "/doi/suppl/10.1056/NEJMoa042786/suppl_file/x.pdf",
        "https://example.org/bin/table%20s1.xlsx",  # percent-encoded
    ],
)
def test_looks_like_supplement_file_accepts_documents(href):
    assert looks_like_supplement_file(href)


@pytest.mark.parametrize(
    "href",
    ["", "#supplementary", "figs/fig1.jpg", "/articles/PMC123/", "icon.svg"],
)
def test_looks_like_supplement_file_rejects_non_documents(href):
    assert not looks_like_supplement_file(href)


def test_extension_is_judged_on_the_path_not_the_query():
    assert not looks_like_supplement_file("https://x.org/view?file=a.pdf")


def test_doi_from_href_extracts_an_embedded_doi():
    href = "/doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm_imboden_2744sa1.pdf"
    assert doi_from_href(href) == "10.1056/NEJMoa042786"
    assert doi_from_href("mmc1.xls") is None


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def test_bare_filename_resolves_against_pmc_bin_for_a_pmc_paper():
    job = resolve_link(SupplementLink(href="mmc1.xls"), pmcid="PMC4535901")
    assert job.is_fetchable()
    assert job.url == "https://pmc.ncbi.nlm.nih.gov/articles/PMC4535901/bin/mmc1.xls"
    assert job.name == "mmc1.xls"
    # base_url is what lets download_supplement expand the PMC variant list.
    assert job.base_url == "https://pmc.ncbi.nlm.nih.gov/articles/PMC4535901/"
    assert "/bin/" in (job.original_url or "")


def test_root_relative_pmc_href_is_rebuilt_from_the_filename():
    """The recorded /articles/instance/... path is discarded for the canonical one."""
    job = resolve_link(
        SupplementLink(href="/articles/instance/4535901/bin/mmc2.pdf"),
        pmcid="PMC4535901",
    )
    assert job.url == "https://pmc.ncbi.nlm.nih.gov/articles/PMC4535901/bin/mmc2.pdf"


def test_absolute_href_is_used_as_recorded():
    url = "https://ars.els-cdn.com/content/image/1-s2.0-x-mmc1.xlsx"
    job = resolve_link(SupplementLink(href=url), pmcid="PMC1")
    assert job.url == url


def test_page_url_resolves_a_relative_href_without_a_pmcid():
    job = resolve_link(
        SupplementLink(href="/doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm.pdf"),
        page_url="https://www.nejm.org/doi/full/10.1056/NEJMoa042786",
    )
    assert job.url == (
        "https://www.nejm.org/doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm.pdf"
    )


def test_doi_host_resolver_recovers_a_stripped_host():
    calls = []

    def resolver(doi):
        calls.append(doi)
        return "www.nejm.org"

    job = resolve_link(
        SupplementLink(
            href="/doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm_2744sa1.pdf"
        ),
        doi_host_resolver=resolver,
    )
    assert calls == ["10.1056/NEJMoa042786"]
    assert job.url == (
        "https://www.nejm.org/doi/suppl/10.1056/NEJMoa042786"
        "/suppl_file/nejm_2744sa1.pdf"
    )


def test_unresolvable_link_reports_instead_of_guessing_a_host():
    """A wrong host silently downloads an error page — record the failure."""
    job = resolve_link(SupplementLink(href="/suppl/table_s1.xlsx"))
    assert not job.is_fetchable()
    assert job.url == ""
    assert "no PMCID" in (job.unresolved_reason or "")
    assert job.name == "table_s1.xlsx"


def test_doi_resolver_is_not_consulted_when_a_pmcid_is_available():
    """PMC is free and deterministic, so it must win over a network lookup."""

    def resolver(doi):  # pragma: no cover - must not be reached
        raise AssertionError("DOI resolution attempted despite a PMCID")

    job = resolve_link(
        SupplementLink(href="/doi/suppl/10.1056/NEJMoa042786/suppl_file/x.pdf"),
        pmcid="PMC7",
        doi_host_resolver=resolver,
    )
    assert job.url == "https://pmc.ncbi.nlm.nih.gov/articles/PMC7/bin/x.pdf"


# ---------------------------------------------------------------------------
# Batch resolution
# ---------------------------------------------------------------------------


def test_resolve_links_skips_non_documents_and_splits_by_fetchability():
    links = parse_links_from_markdown(MARKDOWN_FIXTURE) + [
        SupplementLink(href="figs/fig1.jpg"),
        SupplementLink(href="#supplementary-section"),
    ]
    report = resolve_links(links, pmid="17192539", pmcid=None)
    assert len(report.jobs) == 2  # the two PDFs/XLS only
    assert report.fetchable == []  # no PMCID, no page URL
    assert len(report.unresolved) == 2


def test_resolve_links_dedupes_hrefs_that_resolve_to_one_url():
    links = [
        SupplementLink(href="mmc1.xls"),
        SupplementLink(href="/articles/instance/99/bin/mmc1.xls"),
    ]
    report = resolve_links(links, pmcid="PMC99")
    assert len(report.jobs) == 1


def test_to_supplement_jobs_emits_the_download_path_shape():
    jobs = to_supplement_jobs(
        parse_links_from_markdown(MARKDOWN_FIXTURE),
        pmid="26484152",
        pmcid="PMC4535901",
    )
    assert len(jobs) == 2
    entry = jobs[0]
    # These are the exact keys process_supplement_files / download_supplement read.
    assert set(entry) >= {"url", "name", "base_url", "original_url"}
    assert entry["name"] == "nejm_imboden_2744sa1.pdf"
    assert entry["description"] == "Supplementary Appendix"


def test_to_supplement_jobs_omits_unresolved_links():
    assert to_supplement_jobs([SupplementLink(href="/suppl/x.xlsx")]) == []


# ---------------------------------------------------------------------------
# Europe PMC archive route
#
# NCBI's per-file bin/ URLs 403 an unattended client on every variant the
# scraper generates (verified live 2026-07-24), so the archive endpoint is the
# route that actually works for open-access papers.
# ---------------------------------------------------------------------------


def test_europepmc_archive_job_targets_the_supplementary_files_endpoint():
    job = europepmc_archive_job("PMC2072960")
    assert job.is_fetchable()
    assert job.url == (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2072960/supplementaryFiles"
    )
    assert job.name == "PMC2072960_supplements.zip"
    assert job.source == "pmc_europepmc_archive"


def test_archive_job_carries_no_base_url_so_pmc_variants_are_not_applied():
    """The endpoint is not a PMC bin/ URL; variant rewriting would corrupt it."""
    entry = europepmc_archive_job("PMC1").to_dict()
    assert "base_url" not in entry
    assert entry["name"].endswith(".zip")


def test_zip_magic_distinguishes_an_archive_from_an_error_bean():
    error_bean = b'<?xml version="1.0"?><errorBean><errCode>0</errCode></errorBean>'
    assert not error_bean.startswith(ZIP_MAGIC)
    assert b"PK\x03\x04somezipbody".startswith(ZIP_MAGIC)
