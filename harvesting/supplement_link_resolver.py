"""Turn recorded supplement links into fetchable download jobs.

``figure_extractor`` records every supplement link it finds in a paper's markup
("Supplementary Appendix … Download 32.38 KB"), but for years those links were
*write-only*: rendered into ``FULL_CONTEXT.md`` as ``_link_: <href>`` and read
by nothing. The real download path discovers supplements through the PMC and
publisher APIs instead, so on any paper those APIs missed, the links we already
had sat unused. Measured on the corpus 2026-07-24: **470 papers carry
file-shaped supplement links and have an empty supplements directory**, 430 of
them with a PMCID on record — 1,467 of 1,550 links resolvable against PMC for
free.

This module closes that loop. It does two things and deliberately no more:

1. Recover links from what is already on disk — parse ``_link_:`` lines back
   out of harvested markdown (:func:`parse_links_from_markdown`), or read the
   structured entries from an artifacts log (:func:`links_from_artifacts`).
2. Resolve each link to the ``{url, name, base_url, original_url}`` dicts that
   :meth:`harvesting.orchestrator.PubMedFullTextFetcher.download_supplement`
   and :func:`harvesting.supplement_processing_service.process_supplement_files`
   already consume (:func:`to_supplement_jobs`).

Downloading, PMC URL-variant retries, conversion, and figure extraction stay
where they are. Nothing here fetches anything.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, List, Optional
from urllib.parse import quote, unquote, urljoin, urlparse

# Extensions worth downloading as a body-text supplement. Anything else (images,
# icons, .html landing pages) is either handled by the figure path or is chrome.
SUPPLEMENT_EXTENSIONS: tuple[str, ...] = (
    ".xlsx",
    ".xls",
    ".csv",
    ".tsv",
    ".zip",
    ".pdf",
    ".doc",
    ".docx",
    ".txt",
    ".rtf",
    ".xml",
)

# ``_link_: <href>`` is the line ``figure_extractor._format_supplement_desc_md``
# writes. ``### <label>`` immediately precedes it in the same block.
_LINK_LINE_RE = re.compile(r"^_link_:[ \t]*(\S+)[ \t]*$", re.MULTILINE)
_HEADING_RE = re.compile(r"^###[ \t]+(.+?)[ \t]*$", re.MULTILINE)

# A DOI embedded in a publisher path, e.g.
# /doi/suppl/10.1056/NEJMoa042786/suppl_file/nejm_imboden_2744sa1.pdf
_DOI_IN_PATH_RE = re.compile(r"(10\.\d{4,9}/[^/\s]+)")

_PMCID_RE = re.compile(r"(PMC\d+)", re.IGNORECASE)

EUROPEPMC_REST = "https://www.ebi.ac.uk/europepmc/webservices/rest"

# A ZIP always starts "PK\x03\x04". Europe PMC answers a non-OA article with an
# XML errorBean under a 200, so the status code alone does not tell you whether
# you got an archive.
ZIP_MAGIC = b"PK\x03\x04"


@dataclass
class SupplementLink:
    """One supplement link as recorded by the caption extractor."""

    href: str
    label: str = ""
    title: str = ""
    media_type: Optional[str] = None
    source: str = ""  # where we recovered the link from, for provenance


@dataclass
class SupplementJob:
    """A resolved, fetchable supplement.

    ``to_dict`` produces exactly the shape ``process_supplement_files`` expects,
    including the ``base_url`` / ``original_url`` pair that lets
    ``download_supplement`` generate PMC URL variants.
    """

    url: str
    name: str
    base_url: Optional[str] = None
    original_url: Optional[str] = None
    description: str = ""
    source: str = ""
    unresolved_reason: Optional[str] = None
    doi: Optional[str] = None

    def is_fetchable(self) -> bool:
        return self.unresolved_reason is None and bool(self.url)

    def to_dict(self) -> dict:
        entry = {
            "url": self.url,
            "name": self.name,
            "description": self.description,
            "source": self.source,
        }
        if self.base_url:
            entry["base_url"] = self.base_url
        if self.original_url:
            entry["original_url"] = self.original_url
        return entry


@dataclass
class ResolutionReport:
    """What resolving one paper's links produced."""

    pmid: str = ""
    pmcid: Optional[str] = None
    jobs: List[SupplementJob] = field(default_factory=list)

    @property
    def fetchable(self) -> List[SupplementJob]:
        return [j for j in self.jobs if j.is_fetchable()]

    @property
    def unresolved(self) -> List[SupplementJob]:
        return [j for j in self.jobs if not j.is_fetchable()]


# ---------------------------------------------------------------------------
# Recovering links from what is already on disk
# ---------------------------------------------------------------------------


def looks_like_supplement_file(href: str) -> bool:
    """True when *href* names a downloadable supplement document.

    Judged on the URL path only — query strings routinely carry unrelated
    extensions — and case-insensitively, since PMC serves ``.PDF`` and ``.XLS``
    as often as lowercase.
    """
    if not href:
        return False
    path = urlparse(href.strip()).path or href.strip()
    return unquote(path).lower().endswith(SUPPLEMENT_EXTENSIONS)


def parse_links_from_markdown(
    text: str, *, source: str = "full_context_markdown"
) -> List[SupplementLink]:
    """Recover supplement links from harvested markdown.

    Reads back the ``_link_:`` lines that ``render_captions_markdown`` wrote,
    pairing each with the nearest preceding ``###`` heading as its label. This
    is the only way to reach the links on the thousands of papers harvested
    before they were recorded structurally.
    """
    if not text:
        return []

    headings = [(m.start(), m.group(1).strip()) for m in _HEADING_RE.finditer(text)]
    out: List[SupplementLink] = []
    seen: set[str] = set()

    for match in _LINK_LINE_RE.finditer(text):
        href = match.group(1).strip()
        if not href or href in seen:
            continue
        seen.add(href)
        label = ""
        for pos, heading in headings:
            if pos < match.start():
                label = heading
            else:
                break
        out.append(SupplementLink(href=href, label=label, source=source))

    return out


def links_from_artifacts(manifest: dict) -> List[SupplementLink]:
    """Read structured supplement links out of an artifacts manifest.

    Prefers the ``supplement_links`` list written by newer harvests; falls back
    to nothing when a manifest predates it (use
    :func:`parse_links_from_markdown` on the markdown instead).
    """
    entries = (manifest or {}).get("supplement_links") or []
    out: List[SupplementLink] = []
    for entry in entries:
        # Tolerate a bare href string, and never crash a paper's recovery on a
        # malformed manifest entry.
        if isinstance(entry, str):
            entry = {"href": entry}
        if not isinstance(entry, dict):
            continue
        href = entry.get("href") or ""
        if not href:
            continue
        out.append(
            SupplementLink(
                href=href,
                label=entry.get("label") or "",
                title=entry.get("title") or "",
                media_type=entry.get("media_type"),
                source=entry.get("source") or "artifacts",
            )
        )
    return out


def pmcid_from_artifacts(manifest: dict) -> Optional[str]:
    """Normalise the PMCID out of an artifacts manifest, or None."""
    raw = (manifest or {}).get("pmcid") or ""
    match = _PMCID_RE.search(str(raw))
    return match.group(1).upper() if match else None


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def doi_from_href(href: str) -> Optional[str]:
    """Extract a DOI embedded in a supplement path, if there is one.

    Publisher supplement URLs commonly carry the article DOI
    (``/doi/suppl/10.1056/NEJMoa042786/suppl_file/x.pdf``), which is enough to
    recover the host that a stripped relative href lost.
    """
    if not href:
        return None
    match = _DOI_IN_PATH_RE.search(href)
    if not match:
        return None
    return match.group(1).rstrip(".,;)")


def encode_path_segment(filename: str) -> str:
    """Percent-encode a decoded filename for use as a single URL path segment.

    Filenames are decoded on the way in so they are correct on disk, but a
    decoded name cannot be pasted back into a URL: ``Supplementary Table 1.xlsx``
    would emit a raw space, and a name holding an encoded separator would turn
    into two path segments and address a different object.
    """
    return quote(filename, safe="")


def pmc_supplement_url(pmcid: str, filename: str) -> str:
    """The canonical PMC ``bin/`` URL for a supplement filename.

    One URL only — ``download_supplement`` expands it into the full variant
    list (instance paths, legacy domain, Europe PMC) via ``base_url``, so
    duplicating that logic here would be a second source of truth.

    Note that as of 2026-07-24 NCBI answers these with **HTTP 403** for an
    unattended client, on every variant. Prefer
    :func:`europepmc_archive_job` for any open-access paper and keep this as
    the per-file fallback.
    """
    return (
        f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/bin/"
        f"{encode_path_segment(filename)}"
    )


def pmc_article_url(pmcid: str) -> str:
    return f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"


def europepmc_archive_url(pmcid: str) -> str:
    """Europe PMC's one-request archive of *every* supplement for an article."""
    return f"{EUROPEPMC_REST}/{pmcid}/supplementaryFiles"


def europepmc_archive_job(pmcid: str, *, description: str = "") -> SupplementJob:
    """A single ZIP job covering all of an OA article's supplements.

    This is the route that actually works. Verified live 2026-07-24: NCBI's
    per-file ``bin/`` URLs return 403 for an unattended client on every variant
    the scraper generates, while this endpoint returns ``200 application/zip``
    with the full set (e.g. PMC2072960 → 172 KB holding all five ``.doc``
    supplements NCBI refused). One request per paper instead of one per file,
    no credentials, and the existing ``.zip`` conversion path already unpacks
    and folds the contents.

    Non-open-access articles get a clean ``application/xml`` ``errorBean``
    instead of a ZIP, so callers must verify the payload is really an archive
    before keeping it.
    """
    return SupplementJob(
        url=europepmc_archive_url(pmcid),
        name=f"{pmcid}_supplements.zip",
        description=description or "Europe PMC supplementary files archive",
        source="pmc_europepmc_archive",
    )


def resolve_link(
    link: SupplementLink,
    *,
    pmcid: Optional[str] = None,
    page_url: Optional[str] = None,
    doi_host_resolver: Optional[Callable[[str], Optional[str]]] = None,
) -> SupplementJob:
    """Resolve one recorded link into a fetchable job.

    Resolution order, cheapest and most certain first:

    1. Already absolute — use as recorded.
    2. A PMCID is on record — rebuild the PMC ``bin/`` URL from the filename.
       This covers the large majority of the corpus backlog and costs nothing.
    3. A page URL is known — ``urljoin`` against it (what should have happened
       at capture time; see ``figure_extractor._resolve_href``).
    4. A DOI is embedded in the path and *doi_host_resolver* is supplied —
       resolve the DOI to a host and rebuild. Requires network, so the caller
       opts in by passing the resolver.

    Otherwise the job comes back with ``unresolved_reason`` set rather than a
    guessed URL: a wrong host silently downloads a publisher error page, which
    is worse than a recorded failure.
    """
    href = (link.href or "").strip()
    description = (link.label or link.title or "").strip()
    source = link.source or "supplement_link"

    if not href:
        return SupplementJob(
            url="",
            name="",
            description=description,
            source=source,
            unresolved_reason="empty href",
        )

    # Basename first, decode second: unquoting up front would let an encoded
    # separator ("a%2Fb.csv") split into two segments and silently address a
    # different file.
    filename = unquote(Path(urlparse(href).path or href).name)
    doi = doi_from_href(href)

    def job(
        url: str, base_url: Optional[str], original: Optional[str]
    ) -> SupplementJob:
        return SupplementJob(
            url=url,
            name=filename or unquote(Path(urlparse(url).path).name),
            base_url=base_url,
            original_url=original,
            description=description,
            source=source,
            doi=doi,
        )

    if urlparse(href).scheme in {"http", "https"}:
        # base_url drives PMC URL-variant generation downstream, so it is only
        # meaningful when it describes the same host the href points at. A PMC
        # href gets its own article base; a foreign CDN href gets none.
        href_host = (urlparse(href).hostname or "").lower()
        same_host = href_host == (urlparse(page_url).hostname or "").lower()
        if pmcid and "ncbi.nlm.nih.gov" in href_host:
            return job(href, pmc_article_url(pmcid), href)
        return job(href, page_url if page_url and same_host else None, href)

    if pmcid:
        # PMC serves every supplement out of the article's bin/ directory, so
        # the filename alone is sufficient — which is why a bare filename is
        # fully recoverable for any paper with a PMCID.
        if not filename:
            return SupplementJob(
                url="",
                name="",
                description=description,
                source=source,
                unresolved_reason="no filename in href",
                doi=doi,
            )
        resolved = pmc_supplement_url(pmcid, filename)
        return job(resolved, pmc_article_url(pmcid), resolved)

    if page_url:
        return job(urljoin(page_url, href), page_url, href)

    if doi and doi_host_resolver is not None:
        host = doi_host_resolver(doi)
        if host:
            base = host if "://" in host else f"https://{host}"
            return job(urljoin(base, href), base, href)

    return SupplementJob(
        url="",
        name=filename,
        description=description,
        source=source,
        doi=doi,
        unresolved_reason=(
            "relative href with no PMCID, page URL, or resolvable DOI host"
        ),
    )


def resolve_links(
    links: Iterable[SupplementLink],
    *,
    pmid: str = "",
    pmcid: Optional[str] = None,
    page_url: Optional[str] = None,
    doi_host_resolver: Optional[Callable[[str], Optional[str]]] = None,
    supplements_only: bool = True,
) -> ResolutionReport:
    """Resolve a paper's recorded links, deduplicating by resolved URL."""
    report = ResolutionReport(pmid=pmid, pmcid=pmcid)
    seen: set[str] = set()

    for link in links:
        if supplements_only and not looks_like_supplement_file(link.href):
            continue
        resolved = resolve_link(
            link,
            pmcid=pmcid,
            page_url=page_url,
            doi_host_resolver=doi_host_resolver,
        )
        key = resolved.url or f"unresolved:{link.href}"
        if key in seen:
            continue
        seen.add(key)
        report.jobs.append(resolved)

    return report


def to_supplement_jobs(
    links: Iterable[SupplementLink],
    *,
    pmid: str = "",
    pmcid: Optional[str] = None,
    page_url: Optional[str] = None,
    doi_host_resolver: Optional[Callable[[str], Optional[str]]] = None,
) -> List[dict]:
    """Convenience wrapper: resolved, fetchable ``supp_files`` dicts only."""
    report = resolve_links(
        links,
        pmid=pmid,
        pmcid=pmcid,
        page_url=page_url,
        doi_host_resolver=doi_host_resolver,
    )
    return [job.to_dict() for job in report.fetchable]
