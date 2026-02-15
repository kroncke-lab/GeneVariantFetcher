"""Shared output writing and source metadata helpers for free-text retrieval."""

from __future__ import annotations

import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

from .content_validation import validate_content_quality
from .persistence import append_success_entry


@dataclass(frozen=True)
class FreeTextOutputSource:
    """Source metadata used for free-text success logging and status files."""

    success_marker: str
    status_source: str
    source_tag: str = ""


def source_from_free_text_flags(
    *, used_elsevier_api: bool, used_wiley_api: bool
) -> FreeTextOutputSource:
    """Build source metadata from API usage flags in free-text flow."""
    if used_elsevier_api:
        return FreeTextOutputSource(
            success_marker="ELSEVIER_API",
            status_source="elsevier-api",
            source_tag="[via Elsevier API]",
        )
    if used_wiley_api:
        return FreeTextOutputSource(
            success_marker="WILEY_API",
            status_source="wiley-api",
            source_tag="[via Wiley API]",
        )
    return FreeTextOutputSource(
        success_marker="PUBLISHER_FREE",
        status_source="publisher-free",
        source_tag="[from publisher]",
    )


def publisher_api_fallback_source() -> FreeTextOutputSource:
    """Source metadata for pre-web-scrape publisher API fallback success."""
    return FreeTextOutputSource(
        success_marker="publisher-api",
        status_source="publisher-api-fallback",
    )


def write_free_text_output(
    *,
    output_dir: Path,
    success_log: Path,
    pmid: str,
    main_markdown: str,
    supplement_markdown: str,
    downloaded_count: int,
    source: FreeTextOutputSource,
    write_pmid_status: Optional[Callable[[str, str, Dict[str, object]], None]] = None,
    download_label: str = "Downloaded",
    log_paywalled: Optional[Callable[[str, str, str], None]] = None,
) -> Tuple[Path, str]:
    """Write unified markdown output and append success/status entries.
    
    Returns (output_file, unified_content) on success.
    Returns (None, error_message) if content validation fails.
    """
    unified_content = main_markdown + supplement_markdown
    
    # Validate content quality before writing (catches binary/garbage content)
    is_valid, validation_reason = validate_content_quality(unified_content)
    if not is_valid:
        print(f"  ❌ Content validation failed: {validation_reason}")
        if log_paywalled:
            log_paywalled(
                pmid,
                f"Content validation failed: {validation_reason}",
                f"source: {source.status_source}",
            )
        return None, f"Content validation failed: {validation_reason}"
    
    output_file = output_dir / f"{pmid}_FULL_CONTEXT.md"
    output_file.write_text(unified_content, encoding="utf-8")

    source_suffix = f" {source.source_tag}" if source.source_tag else ""
    print(
        f"  ✅ {download_label}: {output_file.name} ({downloaded_count} supplements){source_suffix}"
    )

    append_success_entry(success_log, pmid, source.success_marker, downloaded_count)

    if write_pmid_status is not None:
        write_pmid_status(
            pmid,
            "extracted",
            {
                "download_timestamp": datetime.datetime.now().isoformat(),
                "variant_count": 0,
                "source": source.status_source,
            },
        )

    return output_file, unified_content
