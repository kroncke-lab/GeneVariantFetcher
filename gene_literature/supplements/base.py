"""Abstract base class for supplement fetchers.

Defines the interface that all tier-specific fetchers must implement,
plus shared utilities for downloading and normalizing supplement files.
"""

import logging
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urlparse

import requests

from utils.retry_utils import api_retry

logger = logging.getLogger(__name__)


@dataclass
class SupplementFile:
    """Metadata for a single supplementary file.

    Compatible with the existing harvesting system's Dict[str, str] format
    (with 'url' and 'name' keys) but adds richer metadata.
    """

    url: str
    name: str
    source: str = ""       # e.g. "pmc", "elsevier", "scraper"
    pmid: str = ""
    pmcid: str = ""
    mime_type: str = ""
    size_bytes: int = 0
    description: str = ""

    def to_dict(self) -> Dict[str, str]:
        """Convert to legacy dict format for compatibility with existing scrapers."""
        return {"url": self.url, "name": self.name}

    @classmethod
    def from_dict(cls, d: Dict[str, str], source: str = "scraper") -> "SupplementFile":
        """Create from legacy dict format."""
        return cls(url=d["url"], name=d.get("name", ""), source=source)

    @property
    def extension(self) -> str:
        """File extension (lowercase, without dot)."""
        path = urlparse(self.url).path
        suffix = Path(path).suffix.lower()
        return suffix.lstrip(".")

    @property
    def normalized_url(self) -> str:
        """URL normalized for deduplication (no fragment, no trailing slash)."""
        return self.url.split("#")[0].rstrip("/")


class SupplementFetcher(ABC):
    """Abstract base for tier-specific supplement fetchers.

    Subclasses must implement:
        fetch(pmid, doi) -> List[SupplementFile]
        list_supplements(pmid, doi) -> List[SupplementFile]

    Provides shared utilities:
        download_file(url, dest) -> Path
        _clean_filename(name) -> str
    """

    VALID_EXTENSIONS = {
        "pdf", "docx", "doc", "xlsx", "xls", "csv", "tsv",
        "zip", "rar", "gz", "tar", "txt", "xml", "html",
        "pptx", "ppt", "json", "rtf", "tiff", "tif", "png",
        "jpg", "jpeg", "svg", "eps", "mp4", "avi", "mov",
    }

    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "GeneVariantFetcher/1.0 (Bot for clinical variant extraction)",
            "Accept": "*/*",
        })

    @abstractmethod
    def fetch(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch supplement file metadata for a paper.

        Args:
            pmid: PubMed ID
            doi: Digital Object Identifier (optional but improves coverage)

        Returns:
            List of SupplementFile objects found for this paper
        """

    @abstractmethod
    def list_supplements(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """List available supplements without downloading.

        Alias for fetch() in most implementations; some fetchers may
        return lighter metadata here (no size checks, etc.).
        """

    @api_retry
    def download_file(self, url: str, dest: Path) -> Optional[Path]:
        """Download a supplement file to disk.

        Args:
            url: URL of the supplement file
            dest: Destination directory (file will be named from URL)

        Returns:
            Path to downloaded file, or None on failure
        """
        dest = Path(dest)
        dest.mkdir(parents=True, exist_ok=True)

        try:
            response = self.session.get(url, timeout=self.timeout, stream=True)
            response.raise_for_status()

            # Determine filename
            filename = self._filename_from_response(response, url)
            filepath = dest / filename

            with open(filepath, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            logger.info(f"Downloaded {filepath.name} ({filepath.stat().st_size:,} bytes)")
            return filepath

        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            return None

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _clean_filename(name: str) -> str:
        """Sanitize a filename for safe filesystem use."""
        cleaned = re.sub(r"[^\w\-_.\s]", "", name)
        cleaned = re.sub(r"\s+", "_", cleaned.strip())
        return cleaned or "supplement"

    @staticmethod
    def _filename_from_response(response: requests.Response, url: str) -> str:
        """Extract filename from Content-Disposition header or URL."""
        cd = response.headers.get("Content-Disposition", "")
        if cd:
            match = re.search(r'filename="?([^";\n]+)"?', cd)
            if match:
                return SupplementFetcher._clean_filename(match.group(1))

        # Fall back to URL path
        path = urlparse(url).path
        name = Path(path).name
        return SupplementFetcher._clean_filename(name) if name else "supplement"

    @staticmethod
    def _has_valid_extension(url: str) -> bool:
        """Check if a URL points to a file with a recognized extension."""
        path = urlparse(url).path
        suffix = Path(path).suffix.lower().lstrip(".")
        return suffix in SupplementFetcher.VALID_EXTENSIONS
