"""Unified supplement fetcher - orchestrates all tiers.

Tries sources in priority order and deduplicates results:
    Tier 1: PMC (free, via Europe PMC + NCBI eUtils)
    Tier 2: Publisher APIs (Elsevier, Wiley, Springer)
    Tier 3: Web scraping (existing harvesting/supplement_scraper.py)

Usage:
    fetcher = UnifiedSupplementFetcher()
    supplements = fetcher.fetch_all("24667783", "10.1093/hmg/ddu...")
    for supp in supplements:
        fetcher.download_file(supp.url, output_dir)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

from .base import SupplementFetcher, SupplementFile
from .pmc_fetcher import PMCSupplementFetcher
from .elsevier_fetcher import ElsevierSupplementFetcher

logger = logging.getLogger(__name__)


class UnifiedSupplementFetcher:
    """Orchestrates supplement retrieval across all tiers.

    Deduplicates by normalized URL so the same supplement isn't
    returned multiple times from different sources.
    """

    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self.pmc_fetcher = PMCSupplementFetcher(timeout=timeout)
        self.elsevier_fetcher = ElsevierSupplementFetcher(timeout=timeout)

        # Ordered list of fetchers to try
        self._fetchers: List[SupplementFetcher] = [
            self.pmc_fetcher,
            self.elsevier_fetcher,
        ]

    def fetch_all(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch supplements from all available tiers, deduplicated.

        Args:
            pmid: PubMed ID
            doi: Digital Object Identifier (optional)

        Returns:
            Deduplicated list of SupplementFile objects from all sources
        """
        all_supplements: List[SupplementFile] = []
        seen_urls: set = set()

        for fetcher in self._fetchers:
            fetcher_name = fetcher.__class__.__name__
            try:
                results = fetcher.fetch(pmid, doi)
                new_count = 0
                for supp in results:
                    norm_url = supp.normalized_url
                    if norm_url not in seen_urls:
                        seen_urls.add(norm_url)
                        supp.pmid = pmid
                        all_supplements.append(supp)
                        new_count += 1

                if results:
                    logger.info(
                        f"{fetcher_name}: {len(results)} found, "
                        f"{new_count} new (after dedup)"
                    )
            except Exception as e:
                logger.warning(f"{fetcher_name} failed for PMID {pmid}: {e}")
                continue

        logger.info(
            f"Unified fetch for PMID {pmid}: {len(all_supplements)} total supplements"
        )
        return all_supplements

    def fetch_tier1(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch from PMC only (free tier)."""
        return self.pmc_fetcher.fetch(pmid, doi)

    def fetch_tier2(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch from publisher APIs only."""
        results: List[SupplementFile] = []
        seen_urls: set = set()

        for fetcher in self._fetchers:
            if isinstance(fetcher, PMCSupplementFetcher):
                continue  # Skip tier 1
            try:
                for supp in fetcher.fetch(pmid, doi):
                    if supp.normalized_url not in seen_urls:
                        seen_urls.add(supp.normalized_url)
                        supp.pmid = pmid
                        results.append(supp)
            except Exception as e:
                logger.warning(f"{fetcher.__class__.__name__} failed: {e}")

        return results

    def to_legacy_format(
        self, supplements: List[SupplementFile]
    ) -> List[Dict[str, str]]:
        """Convert SupplementFile list to legacy dict format.

        For compatibility with existing harvesting/supplement_scraper.py
        consumers that expect List[Dict] with 'url' and 'name' keys.
        """
        return [s.to_dict() for s in supplements]

    def download_all(
        self,
        supplements: List[SupplementFile],
        output_dir: Path,
    ) -> List[Path]:
        """Download all supplement files to a directory.

        Args:
            supplements: List of SupplementFile objects to download
            output_dir: Target directory

        Returns:
            List of successfully downloaded file paths
        """
        output_dir = Path(output_dir)
        downloaded = []

        for supp in supplements:
            path = self.pmc_fetcher.download_file(supp.url, output_dir)
            if path:
                downloaded.append(path)

        logger.info(f"Downloaded {len(downloaded)}/{len(supplements)} supplements")
        return downloaded
