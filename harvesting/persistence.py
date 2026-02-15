"""Persistence helpers for harvesting logs and PMID status files."""

from __future__ import annotations

import csv
import datetime
import json
from pathlib import Path
from typing import Any


def initialize_harvest_logs(paywalled_log: Path, success_log: Path) -> None:
    """Create/reset harvesting CSV logs with headers."""
    with open(paywalled_log, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "PMID",
                "Reason",
                "URL",
                "Classification",
                "Abstract_Carriers",
                "Affected_Count",
                "Unaffected_Count",
                "Variants_Mentioned",
                "Extraction_Confidence",
                "More_In_Fulltext_Probability",
                "Priority_Score",
                "Notes",
            ]
        )

    with open(success_log, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PMID", "PMCID", "Supplements_Downloaded"])


def append_paywalled_entry(
    paywalled_log: Path,
    pmid: str,
    reason: str,
    url: str,
    classification: str = "",
) -> None:
    """Append one paywalled/missing entry row.

    Args:
        paywalled_log: Path to CSV log file.
        pmid: PubMed ID.
        reason: Why the paper couldn't be downloaded.
        url: URL attempted.
        classification: One of PAYWALLED, CAPTCHA_BLOCKED,
            INSTITUTIONAL_ACCESS, SUPPLEMENT_ONLY, API_LIMIT, or empty.
    """
    with open(paywalled_log, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                pmid,
                reason,
                url,
                classification,
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ]
        )


def append_success_entry(
    success_log: Path, pmid: str, pmcid: str, supplements_downloaded: int
) -> None:
    """Append one successful download row."""
    with open(success_log, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([pmid, pmcid, supplements_downloaded])


def write_pmid_status_file(
    output_dir: Path,
    pmid: str,
    status: str,
    details: dict[str, Any] | None = None,
) -> Path:
    """Write status JSON for one PMID and return file path."""
    details = details or {}
    status_dir = output_dir / "pmid_status"
    status_dir.mkdir(exist_ok=True)
    status_file = status_dir / f"{pmid}.json"
    now = datetime.datetime.now().isoformat()

    data = {
        "pmid": pmid,
        "status": status,
        "download_timestamp": details.get("download_timestamp", now),
        "extract_timestamp": details.get("extract_timestamp", now),
        "variant_count": details.get("variant_count"),
        "failure_reason": details.get("failure_reason"),
        "source": details.get("source"),
    }

    with open(status_file, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    return status_file
