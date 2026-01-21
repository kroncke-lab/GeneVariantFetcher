"""NCBI PMC supplementary material harvesting utilities.

This module prefers the official BioC supplementary REST API and falls back to
E-utilities JATS XML when BioC data are unavailable. It provides a small
orchestration layer for downloading supplements and extracting best-effort
text content.
"""

from __future__ import annotations

import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import urljoin
import logging
import requests
import xml.etree.ElementTree as ET


@dataclass
class Supplement:
    """Container for supplementary material metadata and extracted text."""

    pmcid: str
    source: str  # "bioc" or "jats"
    index: Optional[int]
    label: Optional[str]
    href: Optional[str]
    filename: Optional[str]
    mime_type: Optional[str]
    text: Optional[str]


class NCBISettings:
    """Central configuration for NCBI requests and rate limiting."""

    def __init__(
        self,
        user_agent: Optional[str] = None,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
        rate_limit_seconds: float = 0.34,
    ) -> None:
        self.user_agent = user_agent or os.getenv(
            "NCBI_USER_AGENT",
            "GeneVariantFetcher pmc_supp_harvest (contact: dev@example.com)",
        )
        self.email = email or os.getenv("NCBI_EMAIL")
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.rate_limit_seconds = rate_limit_seconds
        self.bioc_base = (
            "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/supplmat.cgi"
        )
        self.efetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        self.elink_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
        self.pmc_root = "https://pmc.ncbi.nlm.nih.gov/"


class _RateLimiter:
    def __init__(self, min_interval: float) -> None:
        self.min_interval = min_interval
        self._last_call = 0.0

    def wait(self) -> None:
        now = time.monotonic()
        delta = now - self._last_call
        if delta < self.min_interval:
            time.sleep(self.min_interval - delta)
        self._last_call = time.monotonic()


class PMCSupplementClient:
    """Client for BioC and JATS-based supplementary retrieval."""

    def __init__(self, settings: Optional[NCBISettings] = None) -> None:
        self.settings = settings or NCBISettings()
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": self.settings.user_agent})
        self._rate_limiter = _RateLimiter(self.settings.rate_limit_seconds)
        self.logger = logging.getLogger(__name__)

    # --------------- HTTP helpers ---------------
    def _get(self, url: str, **kwargs) -> requests.Response:
        self._rate_limiter.wait()
        params = kwargs.pop("params", {})
        if self.settings.api_key and "api_key" not in params:
            params["api_key"] = self.settings.api_key
        if self.settings.email and "email" not in params:
            params["email"] = self.settings.email
        response = self.session.get(url, params=params, timeout=30, **kwargs)
        response.raise_for_status()
        return response

    # --------------- PMID/PMCID helpers ---------------
    def pmid_to_pmcid(self, pmid: str) -> str:
        """Resolve a PMID to a PMCID via E-utilities elink.

        Raises:
            ValueError: when no PMCID mapping is found.
        """
        url = self.settings.elink_base
        params = {
            "dbfrom": "pubmed",
            "db": "pmc",
            "id": pmid,
            "cmd": "neighbor",
            "retmode": "xml",
        }
        resp = self._get(url, params=params)
        root = ET.fromstring(resp.text)
        for link_set in root.findall("LinkSet/LinkSetDb"):
            # Newer ELink responses use <DbTo> instead of LinkSetDbName
            name = link_set.findtext("LinkSetDbName") or link_set.findtext("DbTo")
            if name and name.lower() == "pmc":
                for id_elem in link_set.findall("Link/Id"):
                    pmcid_num = id_elem.text
                    if pmcid_num:
                        return f"PMC{pmcid_num.strip()}"
        raise ValueError(f"No PMCID mapping found for PMID {pmid}")

    # --------------- BioC retrieval ---------------
    def fetch_supplements_bioc(self, pmcid: str) -> List[Supplement]:
        """Retrieve supplementary materials using the BioC REST API.

        This method first enumerates supplements with the `/list` endpoint and
        then downloads combined content via `/all`. It returns Supplement objects
        with extracted text when available.
        """
        list_url = f"{self.settings.bioc_base}/bioc_xml/{pmcid}/list"
        all_url = f"{self.settings.bioc_base}/bioc_xml/{pmcid}/all"

        list_resp = self._get(list_url)
        list_root = ET.fromstring(list_resp.text)
        items: List[Dict[str, Optional[str]]] = []

        for idx, doc in enumerate(list_root.findall(".//document"), start=1):
            info = {
                child.get("key"): (child.text or "") for child in doc.findall("infon")
            }
            items.append(
                {
                    "index": idx,
                    "label": info.get("label") or doc.findtext("id"),
                    "href": info.get("href"),
                }
            )

        all_resp = self._get(all_url)
        all_root = ET.fromstring(all_resp.text)
        supplements: List[Supplement] = []

        # Each <document> typically contains a supplementary item
        documents = all_root.findall(".//document")
        for idx, doc in enumerate(documents, start=1):
            info = {
                child.get("key"): (child.text or "") for child in doc.findall("infon")
            }
            label = info.get("label") or doc.findtext("id")
            href = info.get("href")
            text_chunks: List[str] = []
            for passage in doc.findall("passage"):
                text = passage.findtext("text")
                if text:
                    text_chunks.append(text.strip())
            if not text_chunks:
                doc_text = doc.findtext("passage/text") or doc.findtext("text")
                if doc_text:
                    text_chunks.append(doc_text.strip())
            text_content = "\n\n".join(chunk for chunk in text_chunks if chunk)
            supplements.append(
                Supplement(
                    pmcid=pmcid,
                    source="bioc",
                    index=idx,
                    label=label,
                    href=href,
                    filename=_filename_from_href(href),
                    mime_type=None,
                    text=text_content or None,
                )
            )

        # If list had more metadata than documents, merge href/label where possible
        if items and supplements and len(items) == len(supplements):
            for item, supp in zip(items, supplements):
                supp.label = supp.label or item.get("label")
                supp.href = supp.href or item.get("href")
                if not supp.filename:
                    supp.filename = _filename_from_href(supp.href)

        if not supplements:
            raise ValueError("BioC API returned no supplementary documents")
        return supplements

    # --------------- JATS retrieval ---------------
    def fetch_supplements_jats(self, pmcid: str) -> List[Supplement]:
        """Fetch supplementary file metadata from JATS XML via efetch."""
        params = {"db": "pmc", "id": pmcid, "rettype": "xml"}
        resp = self._get(self.settings.efetch_base, params=params)
        root = ET.fromstring(resp.text)
        supplements: List[Supplement] = []
        for idx, supp in enumerate(root.findall(".//supplementary-material"), start=1):
            href = supp.get("{http://www.w3.org/1999/xlink}href")
            label = supp.findtext("label") or supp.get("id")
            mime = supp.get("mimetype")
            resolved_href = urljoin(self.settings.pmc_root, href) if href else None
            supplements.append(
                Supplement(
                    pmcid=pmcid,
                    source="jats",
                    index=idx,
                    label=label,
                    href=resolved_href,
                    filename=_filename_from_href(resolved_href),
                    mime_type=mime,
                    text=None,
                )
            )
        return supplements

    # --------------- Download and extraction ---------------
    def download_and_extract_text(self, supp: Supplement, out_dir: Path) -> Supplement:
        """Download a supplementary file and attempt to extract text."""
        out_dir.mkdir(parents=True, exist_ok=True)
        if not supp.href:
            supp.text = None
            return supp

        filename = supp.filename or f"supplement_{supp.index or 0}"
        target_path = out_dir / filename

        try:
            resp = self._get(supp.href, stream=True)
        except Exception as exc:  # noqa: BLE001
            self.logger.warning("Failed to download %s: %s", supp.href, exc)
            supp.text = None
            return supp

        mime = resp.headers.get("Content-Type")
        supp.mime_type = mime

        with open(target_path, "wb") as handle:
            for chunk in resp.iter_content(chunk_size=8192):
                handle.write(chunk)

        supp.text = _extract_text_by_type(target_path, mime_type=mime)
        return supp

    # --------------- High-level orchestration ---------------
    def get_all_supplement_text(self, pmcid: str, out_dir: Path) -> List[Supplement]:
        """Retrieve all supplements for a PMCID, preferring BioC over JATS."""
        try:
            supplements = self.fetch_supplements_bioc(pmcid)
            return supplements
        except Exception as exc:  # noqa: BLE001
            self.logger.info(
                "BioC retrieval failed for %s: %s; falling back to JATS", pmcid, exc
            )

        supplements = self.fetch_supplements_jats(pmcid)
        downloaded: List[Supplement] = []
        for supp in supplements:
            downloaded.append(self.download_and_extract_text(supp, out_dir))
        return downloaded

    def harvest_supplements(
        self, ids: Iterable[str], out_dir: Path
    ) -> Dict[str, List[Supplement]]:
        """Harvest supplements for PMCIDs/PMIDs.

        PMIDs are resolved to PMCIDs via elink; if no mapping exists a helpful
        ValueError is raised.
        """
        results: Dict[str, List[Supplement]] = {}
        for raw_id in ids:
            pmcid = raw_id
            if not raw_id.upper().startswith("PMC"):
                pmcid = self.pmid_to_pmcid(raw_id)
            supplements = self.get_all_supplement_text(pmcid, out_dir)
            results[raw_id] = supplements
        return results


def _filename_from_href(href: Optional[str]) -> Optional[str]:
    if not href:
        return None
    return Path(href).name or None


def _extract_text_by_type(path: Path, mime_type: Optional[str]) -> Optional[str]:
    """Extract best-effort text based on MIME type or file extension."""
    ext = path.suffix.lower()
    if mime_type and "text" in mime_type:
        return path.read_text(errors="ignore")
    if ext in {".txt", ".xml", ".html", ".htm"}:
        return path.read_text(errors="ignore")
    if ext in {".doc", ".docx"}:
        return _convert_doc_to_text(path)
    # Unsupported type; return None but keep metadata
    return None


def _convert_doc_to_text(path: Path) -> Optional[str]:
    """Convert a Word document to plain text using LibreOffice if available."""
    try:
        output_dir = path.parent
        subprocess.run(
            [
                "soffice",
                "--headless",
                "--convert-to",
                "txt:Text",
                "--outdir",
                str(output_dir),
                str(path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        txt_path = output_dir / f"{path.stem}.txt"
        if txt_path.exists():
            return txt_path.read_text(errors="ignore")
    except FileNotFoundError:
        logging.getLogger(__name__).warning(
            "LibreOffice not available for doc conversion"
        )
    except subprocess.CalledProcessError as exc:  # noqa: BLE001
        logging.getLogger(__name__).warning("LibreOffice conversion failed: %s", exc)
    return None


if __name__ == "__main__":
    import argparse
    import sys

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Harvest PMC supplementary materials")
    parser.add_argument("ids", nargs="+", help="List of PMCIDs or PMIDs")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("supplements"),
        help="Directory to store downloaded supplements",
    )
    args = parser.parse_args()

    client = PMCSupplementClient()
    try:
        results = client.harvest_supplements(args.ids, args.out)
    except ValueError as exc:
        sys.exit(str(exc))

    for raw_id, supplements in results.items():
        extracted = sum(1 for s in supplements if s.text)
        print(
            f"{raw_id}: {len(supplements)} supplements found; {extracted} with extractable text"
        )
