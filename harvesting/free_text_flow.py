"""Helpers for the initial free-text access fallback flow."""

from __future__ import annotations

import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

from .persistence import append_success_entry
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
                    supplement_markdown += f"\n\n## Supplement: {supp_path.name}\n\n{supp_md}"
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
    unpaywall: Any,
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

    print("  - No PMCID and not marked as free, trying Unpaywall...")

    if doi:
        unpaywall_result, unpaywall_error = unpaywall.find_open_access(doi)
        if unpaywall_result and unpaywall_result.get("pdf_url"):
            print(
                f"  ✓ Unpaywall found OA version: {unpaywall_result.get('oa_status')}"
            )
            pdf_url = unpaywall_result["pdf_url"]
            pdf_path = output_dir / f"{pmid}_unpaywall.pdf"
            success, dl_error = unpaywall.download_pdf(pdf_url, str(pdf_path))
            if success:
                main_markdown = converter.pdf_to_markdown(str(pdf_path))
                if main_markdown and len(main_markdown) > 500:
                    print(f"  ✓ Retrieved via Unpaywall ({len(main_markdown)} chars)")
                    return FreeTextInitState(is_free=True, free_url=pdf_url)
                print("  - Unpaywall PDF conversion failed or content too short")
            else:
                print(f"  - Unpaywall download failed: {dl_error}")
        elif unpaywall_result and unpaywall_result.get("landing_page"):
            print("  - Unpaywall found landing page but no direct PDF")
        else:
            print(f"  - Unpaywall: {unpaywall_error or 'No OA version found'}")

    if doi:
        print("  - Unpaywall failed, trying publisher APIs based on DOI prefix...")
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
            except Exception:
                continue
            if md:
                main_markdown = md
                if provider.should_try:
                    print(f"  ✓ Retrieved via {provider.name} API fallback")
                else:
                    print(f"  ✓ Retrieved via {provider.name} API (unmatched DOI prefix)")
                break

        if main_markdown:
            supp_files = doi_resolver.resolve_and_scrape_supplements(doi, pmid, scraper)
            supplement_markdown = _build_supplement_markdown(
                supp_files=supp_files, converter=converter, logger=logger, pmid=pmid
            )
            unified_content = main_markdown + supplement_markdown
            output_file = output_dir / f"{pmid}_FULL_CONTEXT.md"
            output_file.write_text(unified_content, encoding="utf-8")
            print(
                f"  ✅ Downloaded via publisher API: {output_file.name} ({len(supp_files)} supplements)"
            )

            append_success_entry(success_log, pmid, "publisher-api", len(supp_files))
            write_pmid_status(
                pmid,
                "extracted",
                {
                    "download_timestamp": datetime.datetime.now().isoformat(),
                    "variant_count": 0,
                    "source": "publisher-api-fallback",
                },
            )
            return FreeTextInitState(
                is_free=False, free_url=None, early_result=(True, str(output_file), unified_content)
            )

    print("  ❌ No PMCID and not available via any method (likely paywalled)")
    pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    log_paywalled(
        pmid,
        "No PMCID found, not free full text, Unpaywall failed, publisher APIs failed",
        pubmed_url,
    )
    write_pmid_status(
        pmid,
        "paywalled",
        {
            "download_timestamp": datetime.datetime.now().isoformat(),
            "failure_reason": "No PMCID found, not free full text, all methods failed",
            "source": "none",
        },
    )
    return FreeTextInitState(is_free=False, free_url=None, early_result=(False, "No PMCID", None))
