"""Helpers for the initial free-text access fallback flow."""

from __future__ import annotations

import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

from .free_text_output_service import (
    publisher_api_fallback_source,
    write_free_text_output,
)
from .publisher_strategy import PublisherAttempt, build_publisher_attempt_plan


@dataclass
class FreeTextInitState:
    """State returned from initial free-text fallback checks."""

    is_free: bool
    free_url: Optional[str]
    early_result: Optional[Tuple[bool, str, Optional[str]]] = None


def _build_supplement_markdown(
    supp_files: list,
    converter: Any,
    logger: Any,
    pmid: str,
) -> str:
    """Build supplement markdown for mixed supplement metadata/file inputs."""
    supplement_markdown = ""
    for supp_file in supp_files:
        try:
            if isinstance(supp_file, dict):
                supp_name = supp_file.get("name", "supplement")
                supp_url = supp_file.get("url", "")
                supplement_markdown += (
                    f"\n\n## Supplement: {supp_name}\n\n[Available at: {supp_url}]\n"
                )
                continue

            supp_path = Path(supp_file)
            if supp_path.suffix.lower() == ".pdf":
                supp_md = converter.pdf_to_markdown(str(supp_path))
                if supp_md:
                    supplement_markdown += (
                        f"\n\n## Supplement: {supp_path.name}\n\n{supp_md}"
                    )
        except (TypeError, ValueError) as exc:
            logger.warning(f"Skipping supplement for PMID {pmid}: {exc}")
            continue
    return supplement_markdown


def initialize_free_text_access(
    *,
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    success_log: Path,
    pmc_api: Any,
    converter: Any,
    elsevier_api: Any,
    springer_api: Any,
    wiley_api: Any,
    try_elsevier_api: Callable[[str, str], Tuple[Optional[str], Optional[str]]],
    try_springer_api: Callable[[str, str], Tuple[Optional[str], Optional[str]]],
    try_wiley_api: Callable[[str, str], Tuple[Optional[str], Optional[str]]],
    doi_resolver: Any,
    scraper: Any,
    write_pmid_status: Callable[[str, str, Dict[str, Any]], None],
    log_paywalled: Callable[[str, str, str], None],
    logger: Any,
) -> FreeTextInitState:
    """Run initial checks before full publisher free-text flow."""
    is_free, free_url = pmc_api.is_free_full_text(pmid)

    if is_free:
        return FreeTextInitState(is_free=True, free_url=free_url)

    if doi:
        # OA PDFs are handled by the shared repository fallback after publisher
        # failure, regardless of PubMed's free flag. Never turn converted PDF
        # text back into a URL for the HTML scraper to fetch a second time.
        print("  - Trying publisher APIs based on DOI prefix...")
        provider_plan = build_publisher_attempt_plan(
            doi,
            [
                *(
                    [
                        PublisherAttempt(
                            name="Elsevier",
                            should_try=elsevier_api.is_elsevier_doi(doi),
                            try_fetch=try_elsevier_api,
                        )
                    ]
                    if elsevier_api.is_available
                    else []
                ),
                *(
                    [
                        PublisherAttempt(
                            name="Springer",
                            should_try=springer_api.is_springer_doi(doi),
                            try_fetch=try_springer_api,
                        )
                    ]
                    if springer_api.is_available
                    else []
                ),
                *(
                    [
                        PublisherAttempt(
                            name="Wiley",
                            should_try=wiley_api.is_wiley_doi(doi),
                            try_fetch=try_wiley_api,
                        )
                    ]
                    if wiley_api.is_available
                    else []
                ),
            ],
        )

        main_markdown = None
        for provider in provider_plan:
            if main_markdown:
                break
            try:
                md, _ = provider.try_fetch(doi, pmid)
            except Exception as exc:
                logger.warning(
                    "PMID %s: %s API fallback failed: %s",
                    pmid,
                    provider.name,
                    exc,
                )
                continue
            if md:
                main_markdown = md
                if provider.should_try:
                    print(f"  ✓ Retrieved via {provider.name} API fallback")
                else:
                    print(
                        f"  ✓ Retrieved via {provider.name} API (unmatched DOI prefix)"
                    )
                break

        if main_markdown:
            supp_files = doi_resolver.resolve_and_scrape_supplements(doi, pmid, scraper)
            supplement_markdown = _build_supplement_markdown(
                supp_files=supp_files, converter=converter, logger=logger, pmid=pmid
            )
            output_file, unified_content = write_free_text_output(
                output_dir=output_dir,
                success_log=success_log,
                pmid=pmid,
                main_markdown=main_markdown,
                supplement_markdown=supplement_markdown,
                downloaded_count=len(supp_files),
                source=publisher_api_fallback_source(),
                write_pmid_status=write_pmid_status,
                download_label="Downloaded via publisher API",
                log_paywalled=log_paywalled,
            )
            # Handle validation failure (output_file is None)
            if output_file is None:
                return FreeTextInitState(
                    is_free=False,
                    free_url=None,
                    early_result=(False, unified_content, None),
                )
            return FreeTextInitState(
                is_free=False,
                free_url=None,
                early_result=(True, str(output_file), unified_content),
            )

    print("  - Publisher APIs unavailable; repository/browser fallbacks remain")
    pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    log_paywalled(
        pmid,
        "No PMCID found, not free full text, publisher APIs failed",
        pubmed_url,
    )
    write_pmid_status(
        pmid,
        "paywalled",
        {
            "download_timestamp": datetime.datetime.now().isoformat(),
            "failure_reason": "No PMCID found, not free full text, publisher APIs failed",
            "source": "none",
        },
    )
    return FreeTextInitState(
        is_free=False, free_url=None, early_result=(False, "No PMCID", None)
    )
