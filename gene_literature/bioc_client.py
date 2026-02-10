"""BioC client for fetching structured text from NCBI.

Provides access to:
- PubTator3 API for NER annotations (genes, mutations, diseases, etc.)
- PMC Open Access BioC service for full-text structured sections
"""

import logging
import re
import time
from typing import Any, Dict, List, Optional
from xml.etree import ElementTree as ET

import requests

logger = logging.getLogger(__name__)

# Section classification mapping
_SECTION_PATTERNS = {
    "TITLE": re.compile(r"^(title|article.title)$", re.IGNORECASE),
    "ABSTRACT": re.compile(r"^(abstract|abs|summary)$", re.IGNORECASE),
    "INTRODUCTION": re.compile(
        r"(intro|background|significance)", re.IGNORECASE
    ),
    "RESULTS": re.compile(
        r"(result|finding|observation|outcome|data|phenotyp|genotype.phenotype|"
        r"functional|electrophysiol|patch.clamp|voltage.clamp|characteriz)",
        re.IGNORECASE,
    ),
    "METHODS": re.compile(
        r"(method|material|experiment|procedure|protocol|patient|cohort|"
        r"subject|sample|study.design|statistical|genotyp)",
        re.IGNORECASE,
    ),
    "DISCUSSION": re.compile(
        r"(discuss|conclusion|limitation|perspective|implication)",
        re.IGNORECASE,
    ),
    "TABLE": re.compile(r"(table|tab\b|supplementa)", re.IGNORECASE),
    "FIGURE": re.compile(r"(figure|fig\b)", re.IGNORECASE),
    "REFERENCES": re.compile(r"(reference|bibliograph|cited)", re.IGNORECASE),
}

# Sections most likely to contain variant data
_VARIANT_RICH_SECTIONS = frozenset(
    {"TITLE", "ABSTRACT", "RESULTS", "METHODS", "TABLE"}
)


class BioCClient:
    """Client for NCBI BioC APIs (PubTator3 and PMC Open Access).

    Fetches structured, annotated text from PubMed articles and
    PMC full-text papers via BioC XML endpoints.
    """

    PUBTATOR3_URL = (
        "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/"
        "publications/export/biocxml"
    )
    PMCOA_URL = (
        "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi"
    )

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        max_retries: int = 3,
        min_request_interval: float = 0.34,
    ) -> None:
        self.timeout = timeout
        self.max_retries = max_retries
        self._last_request_time: float = 0
        self._min_request_interval = min_request_interval
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": (
                    "GeneVariantFetcher/1.0 "
                    "(https://github.com/kroncke-lab; "
                    "mailto:brett.kroncke@vanderbilt.edu)"
                ),
            }
        )

    # ------------------------------------------------------------------
    # Rate limiting
    # ------------------------------------------------------------------

    def _rate_limit(self) -> None:
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    # ------------------------------------------------------------------
    # HTTP helper
    # ------------------------------------------------------------------

    def _request(self, url: str, params: Optional[dict] = None) -> str:
        """GET *url* with retries and rate limiting. Returns response text."""
        last_err: Optional[Exception] = None
        for attempt in range(self.max_retries):
            self._rate_limit()
            try:
                resp = self.session.get(
                    url, params=params, timeout=self.timeout
                )
                resp.raise_for_status()
                return resp.text
            except requests.RequestException as exc:
                last_err = exc
                status = getattr(exc.response, "status_code", None)
                if status and 400 <= status < 500 and status != 429:
                    logger.warning(
                        "BioC request failed (HTTP %s, non-retryable): %s",
                        status,
                        url,
                    )
                    raise
                delay = min(2 ** attempt, 30)
                logger.info(
                    "BioC request attempt %d/%d failed (%s), "
                    "retrying in %.1fs",
                    attempt + 1,
                    self.max_retries,
                    exc,
                    delay,
                )
                time.sleep(delay)
        raise last_err  # type: ignore[misc]

    # ------------------------------------------------------------------
    # PubTator3
    # ------------------------------------------------------------------

    def get_pubtator_annotations(self, pmid: str) -> Dict[str, Any]:
        """Fetch PubTator3 annotations for a PubMed article.

        Returns a dict with keys:
            title       – article title (str)
            abstract    – abstract text (str)
            annotations – list of annotation dicts, each containing:
                type, text, start, end, identifiers
        """
        pmid = str(pmid).strip()
        xml_text = self._request(
            self.PUBTATOR3_URL, params={"pmids": pmid}
        )
        return self._parse_pubtator_xml(xml_text)

    def _parse_pubtator_xml(self, xml_text: str) -> Dict[str, Any]:
        root = ET.fromstring(xml_text)
        result: Dict[str, Any] = {
            "title": "",
            "abstract": "",
            "annotations": [],
        }

        for document in root.iter("document"):
            for passage in document.iter("passage"):
                infon_type = self._get_infon(passage, "type")
                text = (passage.findtext("text") or "").strip()

                if infon_type and "title" in infon_type.lower():
                    result["title"] = text
                elif infon_type and "abstract" in infon_type.lower():
                    result["abstract"] = text

                offset_el = passage.find("offset")
                passage_offset = (
                    int(offset_el.text) if offset_el is not None and offset_el.text else 0
                )

                for ann in passage.iter("annotation"):
                    ann_dict = self._parse_annotation(ann, passage_offset)
                    if ann_dict:
                        result["annotations"].append(ann_dict)

        return result

    @staticmethod
    def _get_infon(element: ET.Element, key: str) -> Optional[str]:
        for infon in element.iter("infon"):
            if infon.get("key") == key:
                return infon.text
        return None

    def _parse_annotation(
        self, ann: ET.Element, passage_offset: int
    ) -> Optional[Dict[str, Any]]:
        ann_type = self._get_infon(ann, "type")
        text = (ann.findtext("text") or "").strip()
        if not text:
            return None

        identifiers: List[str] = []
        identifier = self._get_infon(ann, "identifier")
        if identifier:
            identifiers = [
                i.strip() for i in identifier.split(";") if i.strip()
            ]

        loc = ann.find("location")
        start = end = 0
        if loc is not None:
            start = int(loc.get("offset", 0))
            end = start + int(loc.get("length", len(text)))

        return {
            "type": ann_type,
            "text": text,
            "start": start,
            "end": end,
            "identifiers": identifiers,
        }

    # ------------------------------------------------------------------
    # PMC Open Access full text
    # ------------------------------------------------------------------

    def get_pmc_fulltext(self, pmcid: str) -> Dict[str, Any]:
        """Fetch full-text sections from PMC Open Access BioC service.

        *pmcid* should be like ``PMC1234567`` or ``1234567``.

        Returns a dict with keys:
            sections – list of dicts with ``name``, ``type``, ``text``
            raw_text – concatenated full text
        """
        pmcid = str(pmcid).strip()
        if not pmcid.upper().startswith("PMC"):
            pmcid = f"PMC{pmcid}"

        url = f"{self.PMCOA_URL}/BioC_xml/{pmcid}/unicode"
        xml_text = self._request(url)
        return self._parse_pmc_xml(xml_text)

    def _parse_pmc_xml(self, xml_text: str) -> Dict[str, Any]:
        root = ET.fromstring(xml_text)
        sections: List[Dict[str, str]] = []
        text_parts: List[str] = []

        for document in root.iter("document"):
            for passage in document.iter("passage"):
                section_type = self._get_infon(passage, "section_type") or ""
                section_name = self._get_infon(passage, "type") or section_type
                text = (passage.findtext("text") or "").strip()
                if not text:
                    continue

                canonical = self.classify_section(
                    section_name or section_type
                )
                sections.append(
                    {
                        "name": section_name,
                        "type": canonical,
                        "text": text,
                    }
                )
                text_parts.append(text)

        return {
            "sections": sections,
            "raw_text": "\n\n".join(text_parts),
        }

    # ------------------------------------------------------------------
    # Section classification
    # ------------------------------------------------------------------

    @staticmethod
    def classify_section(name: str) -> str:
        """Map a free-form section name to a canonical type.

        Returns one of: TITLE, ABSTRACT, INTRODUCTION, METHODS, RESULTS,
        DISCUSSION, TABLE, FIGURE, REFERENCES, or OTHER.
        """
        if not name:
            return "OTHER"
        for canonical, pattern in _SECTION_PATTERNS.items():
            if pattern.search(name):
                return canonical
        return "OTHER"

    # ------------------------------------------------------------------
    # Variant-rich text
    # ------------------------------------------------------------------

    def get_variant_rich_text(self, pmid: str) -> Dict[str, Any]:
        """Return text from sections most likely to contain variant data.

        Combines PubTator3 annotations with section filtering.  If PMC
        full text is unavailable, falls back to title + abstract only.

        Returns a dict with keys:
            pmid        – the queried PMID
            title       – article title
            abstract    – abstract text
            sections    – list of variant-rich section dicts
            annotations – PubTator3 annotations (mutations only)
            text        – concatenated variant-rich text
        """
        pmid = str(pmid).strip()

        # Always fetch PubTator3 annotations
        pubtator = self.get_pubtator_annotations(pmid)

        mutation_annotations = [
            a
            for a in pubtator["annotations"]
            if a["type"]
            and ("mut" in a["type"].lower() or "variant" in a["type"].lower())
        ]

        result: Dict[str, Any] = {
            "pmid": pmid,
            "title": pubtator["title"],
            "abstract": pubtator["abstract"],
            "sections": [],
            "annotations": mutation_annotations,
            "text": "",
        }

        text_parts = [pubtator["title"], pubtator["abstract"]]

        # Try to get PMC full text
        pmc_data = self._try_pmc_fulltext(pmid)
        if pmc_data:
            for section in pmc_data["sections"]:
                if section["type"] in _VARIANT_RICH_SECTIONS:
                    result["sections"].append(section)
                    text_parts.append(section["text"])

        result["text"] = "\n\n".join(t for t in text_parts if t)
        return result

    def _try_pmc_fulltext(self, pmid: str) -> Optional[Dict[str, Any]]:
        """Try to fetch PMC full text, returning None on failure."""
        try:
            # PubTator3 doesn't give us the PMCID, so we try the
            # NCBI ID converter
            pmcid = self._pmid_to_pmcid(pmid)
            if not pmcid:
                return None
            return self.get_pmc_fulltext(pmcid)
        except Exception as exc:
            logger.debug(
                "PMC full text unavailable for PMID %s: %s", pmid, exc
            )
            return None

    def _pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """Convert PMID to PMCID via NCBI ID Converter API."""
        self._rate_limit()
        try:
            resp = self.session.get(
                "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
                params={
                    "ids": pmid,
                    "format": "json",
                    "tool": "GeneVariantFetcher",
                    "email": "brett.kroncke@vanderbilt.edu",
                },
                timeout=self.timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            records = data.get("records", [])
            if records and "pmcid" in records[0]:
                return records[0]["pmcid"]
        except Exception as exc:
            logger.debug("PMID->PMCID conversion failed for %s: %s", pmid, exc)
        return None
