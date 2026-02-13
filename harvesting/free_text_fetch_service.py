"""Service helpers for main-content retrieval in free-text flow."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Callable, Optional, Tuple
from urllib.parse import urlparse

import requests


@dataclass
class FreeTextContentResult:
    """Result state for free-text main-content retrieval."""

    main_markdown: Optional[str] = None
    final_url: Optional[str] = None
    supp_files: list = field(default_factory=list)
    used_elsevier_api: bool = False
    used_wiley_api: bool = False
    early_result: Optional[Tuple[bool, str, Optional[str]]] = None


def fetch_main_content_for_free_text(
    *,
    pmid: str,
    doi: Optional[str],
    free_url: Optional[str],
    suspicious_free_url_domains: set[str],
    elsevier_api: Any,
    springer_api: Any,
    wiley_api: Any,
    elsevier_client: Any,
    wiley_client: Any,
    session: Any,
    scraper: Any,
    doi_resolver: Any,
    try_elsevier_api: Callable[[str, str], tuple[Optional[str], Optional[str]]],
    try_wiley_api: Callable[[str, str], tuple[Optional[str], Optional[str]]],
    try_springer_api: Callable[[str, str], tuple[Optional[str], Optional[str]]],
    validate_content_quality: Callable[[str, Optional[str]], tuple[bool, str]],
    log_paywalled: Callable[[str, str, str], None],
) -> FreeTextContentResult:
    """Retrieve main markdown + supplement list for a PMID free-text flow."""
    result = FreeTextContentResult()
    api_insufficient_content = False

    if doi and elsevier_api.is_elsevier_doi(doi):
        elsevier_markdown, elsevier_error = try_elsevier_api(doi, pmid)
        if elsevier_markdown:
            result.main_markdown = elsevier_markdown
            result.final_url = f"https://doi.org/{doi}"
            result.used_elsevier_api = True
            result.supp_files = doi_resolver.resolve_and_scrape_supplements(
                doi, pmid, scraper
            )
        elif elsevier_error and "insufficient content" in elsevier_error.lower():
            api_insufficient_content = True
            print("  → Elsevier API returned abstract only, falling back to web scraping...")

    if not result.main_markdown and doi and wiley_api.is_wiley_doi(doi):
        wiley_markdown, wiley_error = try_wiley_api(doi, pmid)
        if wiley_markdown:
            result.main_markdown = wiley_markdown
            result.final_url = f"https://doi.org/{doi}"
            result.used_wiley_api = True
            result.supp_files = doi_resolver.resolve_and_scrape_supplements(
                doi, pmid, scraper
            )
        elif wiley_error and "insufficient content" in wiley_error.lower():
            api_insufficient_content = True
            print("  → Wiley API returned abstract only, falling back to web scraping...")

    if not result.main_markdown and doi and springer_api.is_springer_doi(doi):
        springer_markdown, springer_error = try_springer_api(doi, pmid)
        if springer_markdown:
            result.main_markdown = springer_markdown
            result.final_url = f"https://doi.org/{doi}"
            result.supp_files = doi_resolver.resolve_and_scrape_supplements(
                doi, pmid, scraper
            )
        elif springer_error and "not openaccess" in springer_error.lower():
            print("  → Springer API: article not OpenAccess, falling back to web scraping...")

    if not result.main_markdown and doi:
        result.main_markdown, result.final_url, result.supp_files = (
            doi_resolver.resolve_and_fetch_fulltext(doi, pmid, scraper)
        )
        if result.main_markdown:
            is_valid, reason = validate_content_quality(
                result.main_markdown, f"https://doi.org/{doi}"
            )
            if not is_valid:
                print(f"  ⚠ DOI content validation failed: {reason}")
                print("  ❌ Skipping - does not contain valid full article text")
                log_paywalled(
                    pmid,
                    f"DOI content validation failed: {reason}",
                    f"https://doi.org/{doi}",
                )
                result.main_markdown = None

    if not result.main_markdown and free_url:
        parsed_url = urlparse(free_url)
        if parsed_url.netloc in suspicious_free_url_domains:
            print(
                f"  - Skipping suspicious free URL on {parsed_url.netloc} (likely non-article content)"
            )
            log_paywalled(pmid, "Suspicious free URL skipped", free_url)
            result.early_result = (False, "Suspicious free URL", None)
            return result

        print(f"  - No DOI, attempting to fetch from free URL: {free_url}")

        if elsevier_api.is_available and elsevier_api.is_elsevier_url(free_url):
            pii = elsevier_api.extract_pii_from_url(free_url)
            if pii:
                print(f"  Trying Elsevier API (with circuit breaker) for PII: {pii}")
                try:
                    elsevier_markdown, _ = elsevier_client.fetch_fulltext(pii=pii)
                    if elsevier_markdown:
                        result.main_markdown = elsevier_markdown
                        result.final_url = free_url
                        result.used_elsevier_api = True
                        print(
                            f"  ✓ Full text retrieved via Elsevier API ({len(result.main_markdown)} characters)"
                        )
                except Exception as exc:
                    if "circuit breaker" in str(exc).lower():
                        print(
                            f"  ⚠ Circuit breaker protection activated for Elsevier: {exc}"
                        )
                    else:
                        print(f"  - Elsevier API failed: {exc}")
                    try:
                        response = session.get(free_url, allow_redirects=True, timeout=10)
                        response.raise_for_status()
                        result.supp_files = scraper.scrape_elsevier_supplements(
                            response.text, response.url
                        )
                    except Exception:
                        result.supp_files = []

        if (
            not result.main_markdown
            and wiley_api.is_available
            and wiley_api.is_wiley_url(free_url)
        ):
            extracted_doi = wiley_api.extract_doi_from_url(free_url)
            if extracted_doi:
                print(
                    f"  Trying Wiley API (with circuit breaker) for DOI: {extracted_doi}"
                )
                try:
                    wiley_markdown, _ = wiley_client.fetch_fulltext(doi=extracted_doi)
                    if wiley_markdown:
                        result.main_markdown = wiley_markdown
                        result.final_url = free_url
                        result.used_wiley_api = True
                        print(
                            f"  ✓ Full text retrieved via Wiley API ({len(result.main_markdown)} characters)"
                        )
                except Exception as exc:
                    if "circuit breaker" in str(exc).lower():
                        print(
                            f"  ⚠ Circuit breaker protection activated for Wiley: {exc}"
                        )
                    else:
                        print(f"  - Wiley API failed: {exc}")
                    try:
                        response = session.get(free_url, allow_redirects=True, timeout=10)
                        response.raise_for_status()
                        result.supp_files = scraper.scrape_generic_supplements(
                            response.text, response.url
                        )
                    except Exception:
                        result.supp_files = []

        if not result.main_markdown:
            try:
                response = session.get(free_url, allow_redirects=True, timeout=10)
                response.raise_for_status()
                result.final_url = response.url
                print("  ✓ Retrieved free full text page")

                domain = urlparse(result.final_url).netloc
                if "linkinghub.elsevier.com" in domain:
                    try:
                        pii_match = re.search(r"/pii/([^/?]+)", result.final_url)
                        if pii_match:
                            pii = pii_match.group(1)
                            sciencedirect_url = (
                                f"https://www.sciencedirect.com/science/article/pii/{pii}"
                            )
                            print(
                                f"  → Attempting to access ScienceDirect page: {sciencedirect_url}"
                            )
                            redirect_response = session.get(
                                sciencedirect_url, allow_redirects=True, timeout=30
                            )
                            redirect_response.raise_for_status()
                            result.final_url = redirect_response.url
                            response = redirect_response
                            domain = urlparse(result.final_url).netloc
                    except Exception as exc:
                        print(f"  - Could not follow redirect from linkinghub: {exc}")

                html_content = response.text
                result.main_markdown, _ = scraper.extract_fulltext(
                    html_content, result.final_url
                )
                if result.main_markdown:
                    is_valid, reason = validate_content_quality(
                        result.main_markdown, result.final_url
                    )
                    if is_valid:
                        print(
                            f"  ✓ Extracted full text ({len(result.main_markdown)} characters)"
                        )
                    else:
                        print(f"  ⚠ Content validation failed: {reason}")
                        print(
                            "  ❌ Skipping this URL as it doesn't contain valid article content"
                        )
                        log_paywalled(
                            pmid,
                            f"Content validation failed: {reason}",
                            free_url,
                        )
                        result.main_markdown = None
                else:
                    print("  ❌ Could not extract full text from page")

                domain = urlparse(result.final_url).netloc
                if "nature.com" in domain:
                    result.supp_files = scraper.scrape_nature_supplements(
                        html_content, result.final_url
                    )
                elif any(
                    d in domain
                    for d in ["gimjournal.org", "sciencedirect.com", "elsevier.com"]
                ):
                    result.supp_files = scraper.scrape_elsevier_supplements(
                        html_content, result.final_url
                    )
                else:
                    result.supp_files = scraper.scrape_generic_supplements(
                        html_content, result.final_url
                    )
            except requests.exceptions.RequestException as exc:
                print(f"  ❌ Failed to fetch free full text from {free_url}: {exc}")
                log_paywalled(pmid, f"Free full text fetch failed: {exc}", free_url)
                result.early_result = (False, "Free text fetch failed", None)
                return result
    elif not result.main_markdown:
        print("  ❌ No DOI or free URL available to fetch full text")
        print("     (PubMed metadata indicates free access but provides no usable link)")
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        log_paywalled(
            pmid,
            "Free full text indicated but no DOI or URL in PubMed metadata",
            pubmed_url,
        )
        result.early_result = (False, "No DOI or URL for free text", None)
        return result

    if not result.main_markdown:
        print("  ❌ Could not retrieve full text from publisher")
        fallback_url = result.final_url or (
            f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        )
        log_paywalled(pmid, "Free full text extraction failed", fallback_url)
        result.early_result = (False, "Free text extraction failed", None)
        return result

    print(
        f"  ✓ Full text retrieved from publisher ({len(result.main_markdown)} characters)"
    )
    _ = api_insufficient_content  # Maintained for behavioral parity/debug trace.
    return result
