"""Abstract harvesting utilities.

This module fetches PubMed metadata and abstracts for a list of PMIDs and
persists them to JSON files for downstream processing. Each output file follows
this schema:

```
{
  "metadata": {
    "pmid": str,
    "title": str | null,
    "authors": [str],
    "journal": str | null,
    "year": str | null
  },
  "abstract": str | null
}
```
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from utils.pubmed_utils import batch_fetch_metadata, fetch_paper_abstract

logger = logging.getLogger(__name__)


def _normalize_authors(raw_authors: Optional[List[Any]]) -> List[str]:
    """Convert Entrez author list into plain strings."""

    authors: List[str] = []

    if not raw_authors:
        return authors

    for author in raw_authors:
        if isinstance(author, dict):
            name = author.get("Name")
            if not name:
                name = " ".join(
                    part for part in [author.get("ForeName"), author.get("LastName")] if part
                )
            if name:
                authors.append(str(name))
        else:
            authors.append(str(author))

    return authors


def _extract_year(metadata: Dict[str, str]) -> Optional[str]:
    """Derive publication year from common PubMed metadata fields."""

    for field in ("PubDate", "EPubDate", "SortPubDate"):
        date_value = metadata.get(field)
        if not date_value:
            continue

        year = str(date_value).strip().split(" ")[0]
        if year and year[:4].isdigit():
            return year[:4]

    return None


def fetch_and_save_abstracts(
    pmids: List[str], output_dir: str, email: Optional[str] = None
) -> Dict[str, Path]:
    """Fetch metadata and abstracts for PMIDs and persist them to JSON files."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    saved_files: Dict[str, Path] = {}

    all_metadata = batch_fetch_metadata(pmids, email=email)
    logger.info("Fetched metadata for %d PMIDs.", len(all_metadata))

    for pmid in pmids:
        logger.info("Processing PMID %s", pmid)

        metadata = all_metadata.get(pmid, {})
        abstract = fetch_paper_abstract(pmid, email=email)

        record = {
            "metadata": {
                "pmid": pmid,
                "title": metadata.get("Title"),
                "authors": _normalize_authors(metadata.get("AuthorList")),
                "journal": metadata.get("FullJournalName") or metadata.get("Source"),
                "year": _extract_year(metadata),
            },
            "abstract": abstract,
        }

        output_file = output_path / f"{pmid}.json"
        with output_file.open("w", encoding="utf-8") as f:
            json.dump(record, f, indent=2)

        saved_files[pmid] = output_file
        logger.info("Saved abstract JSON to %s", output_file)

    return saved_files
