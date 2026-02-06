"""
Orchestrator Module

Main PMCHarvester class that coordinates all harvesting operations:
- Converts PMIDs to PMCIDs
- Downloads full-text and supplemental files
- Creates unified markdown files for LLM processing
"""

import os
import re
import time
import csv
import logging
import requests
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlparse

from .pmc_api import PMCAPIClient
from .doi_resolver import DOIResolver
from .supplement_scraper import SupplementScraper
from .format_converters import FormatConverter
from .elsevier_api import ElsevierAPIClient
from .wiley_api import WileyAPIClient
from .springer_api import SpringerAPIClient
from .unpaywall_api import UnpaywallClient, get_unpaywall_client
from .retry_manager import RetryManager, RetryConfig
from .priority_queue import PriorityQueue, Priority, Status as QueueStatus
from .core_api import COREAPIClient, get_core_client
from .supplement_reference_parser import (
    parse_supplement_references,
    extract_supplement_urls_from_text,
    check_supplement_gap,
)
from ..utils.resilience import CircuitBreaker, ResilientAPIClient


# Setup logging for circuit breaker monitoring
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# INPUT VALIDATION
# =============================================================================

# PMID pattern: 1-8 digit numeric string
PMID_PATTERN = re.compile(r"^\d{1,8}$")


class ValidationError(Exception):
    """Raised when input validation fails."""

    pass


def validate_pmid_format(pmid: str) -> bool:
    """
    Validate that a PMID has correct format.

    Args:
        pmid: PMID string to validate

    Returns:
        True if valid, False otherwise
    """
    if not pmid:
        return False
    # Strip whitespace and check format
    return bool(PMID_PATTERN.match(str(pmid).strip()))


def validate_pmid_list(pmids: List[str]) -> Tuple[List[str], List[str]]:
    """
    Validate a list of PMIDs and return valid/invalid lists.

    Args:
        pmids: List of PMID strings

    Returns:
        Tuple of (valid_pmids, invalid_pmids)
    """
    valid = []
    invalid = []

    for pmid in pmids:
        pmid_str = str(pmid).strip()
        if validate_pmid_format(pmid_str):
            valid.append(pmid_str)
        else:
            invalid.append(pmid_str)

    return valid, invalid


def validate_output_directory(output_dir: Path) -> None:
    """
    Validate that output directory is writable.

    Args:
        output_dir: Path to output directory

    Raises:
        ValidationError: If output directory is not writable
    """
    if output_dir.exists():
        if not output_dir.is_dir():
            raise ValidationError(
                f"Output path exists but is not a directory: {output_dir}"
            )
        if not os.access(output_dir, os.W_OK):
            raise ValidationError(f"Output directory is not writable: {output_dir}")
    else:
        # Check if parent directory is writable (so we can create output_dir)
        parent = output_dir.parent
        if parent.exists() and not os.access(parent, os.W_OK):
            raise ValidationError(
                f"Cannot create output directory (parent not writable): {output_dir}"
            )


def validate_harvest_inputs(pmids: List[str], output_dir: Path) -> List[str]:
    """
    Validate all inputs for the harvest operation.

    Args:
        pmids: List of PMIDs to process
        output_dir: Path to output directory

    Returns:
        List of valid PMIDs (after filtering out invalid ones)

    Raises:
        ValidationError: If critical validation fails
    """
    # Check PMID list is not empty
    if not pmids:
        raise ValidationError(
            "PMID list is empty. Provide at least one PMID to process."
        )

    # Validate PMIDs format
    valid_pmids, invalid_pmids = validate_pmid_list(pmids)

    if invalid_pmids:
        logger.warning(
            f"Skipping {len(invalid_pmids)} invalid PMIDs: {invalid_pmids[:5]}"
            + (
                f"... and {len(invalid_pmids) - 5} more"
                if len(invalid_pmids) > 5
                else ""
            )
        )

    if not valid_pmids:
        raise ValidationError(
            f"No valid PMIDs found. All {len(pmids)} PMIDs have invalid format. "
            f"PMIDs must be 1-8 digit numbers (e.g., 12345678)."
        )

    # Validate output directory
    validate_output_directory(output_dir)

    return valid_pmids


# Import scout components (with fallback for import errors)
try:
    from pipeline.data_scout import GeneticDataScout
    from config.settings import get_settings

    SCOUT_AVAILABLE = True
except ImportError:
    SCOUT_AVAILABLE = False
    get_settings = None

# Import pedigree extractor (with fallback)
try:
    from pipeline.pedigree_extractor import PedigreeExtractor

    PEDIGREE_EXTRACTOR_AVAILABLE = True
except ImportError:
    PEDIGREE_EXTRACTOR_AVAILABLE = False
    PedigreeExtractor = None

# Import manifest utilities for tracking download outcomes
from utils.manifest import Manifest, ManifestEntry, Status, Stage

import json
import logging

logger = logging.getLogger(__name__)


class PMCHarvester:
    """Harvests full-text and supplemental materials from PubMed Central."""

    SUSPICIOUS_FREE_URL_DOMAINS = {"antibodies.cancer.gov"}

    def __init__(self, output_dir: str = "pmc_harvest", gene_symbol: str = None):
        """
        Initialize PMC Harvester.

        Args:
            output_dir: Directory to save harvested files
            gene_symbol: Target gene symbol for data scout analysis
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.paywalled_log = self.output_dir / "paywalled_missing.csv"
        self.success_log = self.output_dir / "successful_downloads.csv"
        self.gene_symbol = gene_symbol

        # Initialize session with browser-like headers
        self.session = requests.Session()
        HEADERS = {
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
            "Accept-Language": "en-US,en;q=0.9",
            "Accept-Encoding": "gzip, deflate, br",
            "Referer": "https://pubmed.ncbi.nlm.nih.gov/",
            "DNT": "1",
            "Connection": "keep-alive",
            "Upgrade-Insecure-Requests": "1",
        }
        self.session.headers.update(HEADERS)

        # Initialize component modules
        self.pmc_api = PMCAPIClient(session=self.session)
        self.doi_resolver = DOIResolver(
            session=self.session, paywalled_log=self.paywalled_log
        )
        self.scraper = SupplementScraper()
        self.converter = FormatConverter()

        # Initialize publisher API clients (optional - uses API keys from settings/env if available)
        elsevier_api_key = None
        wiley_api_key = None
        springer_api_key = None
        if get_settings is not None:
            try:
                settings = get_settings()
                elsevier_api_key = settings.elsevier_api_key
                wiley_api_key = settings.wiley_api_key
                springer_api_key = getattr(settings, "springer_api_key", None)
            except Exception:
                pass  # Settings validation may fail if other keys are missing

        # Fall back to environment variables if not in settings
        if springer_api_key is None:
            springer_api_key = os.environ.get("SPRINGER_API_KEY")

        self.elsevier_api = ElsevierAPIClient(
            api_key=elsevier_api_key, session=self.session
        )
        self.wiley_api = WileyAPIClient(api_key=wiley_api_key, session=self.session)
        self.springer_api = SpringerAPIClient(
            api_key=springer_api_key, session=self.session
        )
        
        # Initialize Unpaywall client (free service, uses NCBI_EMAIL)
        ncbi_email = os.environ.get("NCBI_EMAIL", "gvf@example.com")
        self.unpaywall = UnpaywallClient(email=ncbi_email, session=self.session)
        
        # Initialize CORE API client (optional, uses CORE_API_KEY)
        core_api_key = os.environ.get("CORE_API_KEY")
        self.core_api = COREAPIClient(api_key=core_api_key, session=self.session)
        
        # Initialize retry manager for transient failures
        self.retry_manager = RetryManager(
            config=RetryConfig(max_retries=3, base_delay=2.0),
            log_file=self.output_dir / "retry_log.jsonl"
        )
        
        # Initialize circuit breakers for API resilience
        self.elsevier_circuit = CircuitBreaker("elsevier", max_failures=5, reset_timeout=60)
        self.wiley_circuit = CircuitBreaker("wiley", max_failures=5, reset_timeout=60)
        self.springer_circuit = CircuitBreaker("springer", max_failures=5, reset_timeout=60)
        
        # Create resilient API clients with circuit breaker protection
        self.elsevier_client = ResilientAPIClient(self.elsevier_api, self.elsevier_circuit)
        self.wiley_client = ResilientAPIClient(self.wiley_api, self.wiley_circuit)
        self.springer_client = ResilientAPIClient(self.springer_api, self.springer_circuit)
        
        # Initialize priority queue for manual acquisition
        self.priority_queue = PriorityQueue(self.output_dir / "manual_acquisition_queue.json")

        # Initialize log files with extended columns for carrier data
        with open(self.paywalled_log, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "PMID",
                    "Reason",
                    "URL",
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

        with open(self.success_log, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["PMID", "PMCID", "Supplements_Downloaded"])

    def _log_paywalled(self, pmid: str, reason: str, url: str) -> None:
        """
        Log a paper to the paywalled/missing CSV with placeholder columns.

        Args:
            pmid: PubMed ID
            reason: Why the paper couldn't be downloaded
            url: URL attempted or PubMed URL
        """
        with open(self.paywalled_log, "a", newline="") as f:
            writer = csv.writer(f)
            # Write with empty placeholders for carrier columns
            writer.writerow(
                [
                    pmid,
                    reason,
                    url,
                    "",
                    "",
                    "",  # Abstract_Carriers, Affected_Count, Unaffected_Count
                    "",
                    "",  # Variants_Mentioned, Extraction_Confidence
                    "",
                    "",
                    "",  # More_In_Fulltext_Probability, Priority_Score, Notes
                ]
            )
        
        # Also add to priority queue for manual acquisition
        try:
            self.priority_queue.add_paper(
                pmid=pmid,
                failure_reason=reason,
                priority=Priority.MEDIUM,  # Default priority
                failed_sources=[url] if url else [],
            )
        except Exception as e:
            logger.warning(f"Failed to add {pmid} to priority queue: {e}")

    def _write_pmid_status(self, pmid: str, status: str, details: Dict[str, Any] = None) -> None:
        """
        Write the status of a PMID to a JSON file.

        Args:
            pmid: PubMed ID
            status: Status of the PMID (downloaded, extracted, failed, paywalled, no_variants)
            details: Dictionary of details to write to the file
        """
        import json
        import datetime
        status_dir = self.output_dir / "pmid_status"
        status_dir.mkdir(exist_ok=True)
        status_file = status_dir / f"{pmid}.json"

        data = {
            "pmid": pmid,
            "status": status,
            "download_timestamp": details.get("download_timestamp", datetime.datetime.now().isoformat()),
            "extract_timestamp": details.get("extract_timestamp", datetime.datetime.now().isoformat()),
            "variant_count": details.get("variant_count", None),
            "failure_reason": details.get("failure_reason", None),
            "source": details.get("source", None),
        }

        with open(status_file, "w") as f:
            json.dump(data, f, indent=2)


    # Patterns that indicate junk/non-article content
    JUNK_CONTENT_PATTERNS = [
        "we've detected you are running an older version",
        "browsehappy.com",
        "phosphositeplus",
        "uniprot database entry",
        "sign in to access",
        "subscribe to read",
        "please enable javascript",
        "your browser does not support",
        "access denied",
        "403 forbidden",
        "404 not found",
        "page not found",
        "cookies must be enabled",
        "this site requires javascript",
    ]

    # Domains that return non-article content
    JUNK_CONTENT_DOMAINS = {
        "assays.cancer.gov",
        "antibodies.cancer.gov",
        "biocyc.org",
        "glygen.org",
        "malacards.org",
        "lens.org",
        "clinicaltrials.gov",
        "medlineplus.gov",
    }

    def _validate_content_quality(
        self, content: str, source_url: str = None
    ) -> Tuple[bool, str]:
        """
        Validate that extracted content is actual paper text, not junk.

        Args:
            content: The extracted markdown/text content
            source_url: Optional URL the content came from

        Returns:
            Tuple of (is_valid, reason)
        """
        if not content:
            return False, "Empty content"

        content_lower = content.lower()

        # Check for junk patterns
        for pattern in self.JUNK_CONTENT_PATTERNS:
            if pattern in content_lower:
                return False, f"Junk content detected: '{pattern}'"

        # Check source URL domain
        if source_url:
            try:
                domain = urlparse(source_url).netloc.lower()
                for junk_domain in self.JUNK_CONTENT_DOMAINS:
                    if junk_domain in domain:
                        return False, f"Content from non-article domain: {junk_domain}"
            except Exception:
                pass

        # Check minimum length
        if len(content) < 1500:
            return False, f"Content too short ({len(content)} chars)"

        # Check for paper-like structure (at least 2 of these indicators)
        paper_indicators = [
            "abstract",
            "introduction",
            "methods",
            "results",
            "discussion",
            "conclusion",
            "references",
            "materials and methods",
            "patients and methods",
            "study population",
        ]
        indicator_count = sum(1 for ind in paper_indicators if ind in content_lower)
        if indicator_count < 2:
            # Also check for variant-related content as an alternative indicator
            variant_patterns = [
                r"[A-Z]\d{2,4}[A-Z]",  # A561V
                r"p\.[A-Z][a-z]{2}\d+",  # p.Ala561
                r"c\.\d+[ACGT]>[ACGT]",  # c.1234A>G
                r"mutation",
                r"variant",
            ]
            variant_matches = sum(
                1 for pat in variant_patterns if re.search(pat, content, re.IGNORECASE)
            )
            if variant_matches < 2:
                return (
                    False,
                    f"Missing paper structure (only {indicator_count} indicators) and no variant content",
                )

        return True, "Valid content"

    def _extract_pmc_figures(self, pmcid: str, pmid: str) -> int:
        """
        Extract figure images from a PMC article page.

        PMC articles have figures embedded in the HTML, not as separate downloads.
        This method scrapes the article page to find and download figure images.

        Args:
            pmcid: PubMed Central ID
            pmid: PubMed ID (for output directory naming)

        Returns:
            Number of figures extracted
        """
        # Check if figure extraction is enabled
        try:
            if get_settings is not None:
                settings = get_settings()
                if not settings.extract_figures:
                    return 0
        except Exception:
            pass

        if not pmcid:
            return 0

        # Create figures directory
        figures_dir = self.output_dir / f"{pmid}_figures"
        figures_dir.mkdir(exist_ok=True)

        try:
            from bs4 import BeautifulSoup

            # Fetch the PMC article page
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            response = self.session.get(pmc_url, timeout=30)
            response.raise_for_status()

            soup = BeautifulSoup(response.text, "html.parser")

            # Find figure images - PMC uses cdn.ncbi.nlm.nih.gov URLs
            figure_urls = []
            for img in soup.find_all("img"):
                src = img.get("src", "") or img.get("data-src", "")
                if src and "cdn.ncbi.nlm.nih.gov/pmc" in src:
                    # Skip logos and small images
                    if "logo" in src.lower() or "icon" in src.lower():
                        continue
                    # Only keep figure images (usually .jpg or .png)
                    if src.endswith((".jpg", ".jpeg", ".png", ".gif")):
                        figure_urls.append(src)

            if not figure_urls:
                logger.debug(f"No figures found in PMC article {pmcid}")
                return 0

            # Download each figure
            downloaded = 0
            for i, url in enumerate(figure_urls, 1):
                try:
                    # Determine file extension
                    ext = ".jpg"
                    for e in [".png", ".gif", ".jpeg"]:
                        if e in url.lower():
                            ext = e
                            break

                    filename = f"fig_pmc_{i}{ext}"
                    filepath = figures_dir / filename

                    # Download the image
                    img_response = self.session.get(url, timeout=30)
                    img_response.raise_for_status()

                    # Validate it's actually an image (not HTML error page)
                    content = img_response.content
                    if content[:4] in (
                        b"\x89PNG",
                        b"\xff\xd8\xff\xe0",
                        b"\xff\xd8\xff\xe1",
                        b"GIF8",
                    ):
                        filepath.write_bytes(content)
                        downloaded += 1
                        logger.debug(
                            f"Downloaded figure: {filename} ({len(content)} bytes)"
                        )
                    elif content[:4] == b"\xff\xd8\xff\xdb":  # Another JPEG variant
                        filepath.write_bytes(content)
                        downloaded += 1
                    else:
                        logger.debug(f"Skipped non-image content from {url}")

                except Exception as e:
                    logger.warning(f"Failed to download figure from {url}: {e}")
                    continue

            if downloaded > 0:
                print(f"  ✓ Extracted {downloaded} figures from PMC article")
                logger.info(f"Extracted {downloaded} figures from PMC article {pmcid}")

            return downloaded

        except Exception as e:
            logger.warning(f"Failed to extract figures from PMC article {pmcid}: {e}")
            return 0

    def _run_pedigree_extraction(self, pmid: str) -> Optional[str]:
        """
        Run pedigree extraction on figures extracted from this paper.

        Args:
            pmid: PubMed ID

        Returns:
            Markdown summary of pedigree findings, or None if no pedigrees found
        """
        if not PEDIGREE_EXTRACTOR_AVAILABLE:
            return None

        try:
            settings = get_settings()
            if not settings.extract_pedigrees:
                return None
        except Exception:
            return None

        # Look for figures directory
        figures_dir = self.output_dir / f"{pmid}_figures"
        if not figures_dir.exists():
            logger.debug(f"No figures directory for PMID {pmid}")
            return None

        # Check if there are any images
        image_files = list(figures_dir.glob("*.png")) + list(figures_dir.glob("*.jpg"))
        if not image_files:
            return None

        try:
            settings = get_settings()
            extractor = PedigreeExtractor(
                model=settings.vision_model,
                detection_confidence_threshold=settings.pedigree_confidence_threshold,
            )

            print(f"  Analyzing {len(image_files)} figures for pedigrees...")
            results = extractor.process_figures_directory(figures_dir)

            if results:
                # Save pedigree results as JSON
                pedigree_json_path = self.output_dir / f"{pmid}_PEDIGREES.json"
                with open(pedigree_json_path, "w") as f:
                    json.dump(results, f, indent=2, default=str)

                # Generate markdown summary
                summary = extractor.summarize_for_extraction(results)

                # Also save as standalone markdown
                pedigree_md_path = self.output_dir / f"{pmid}_PEDIGREES.md"
                pedigree_md_path.write_text(summary)

                pedigree_count = len(results)
                total_individuals = sum(len(r.get("individuals", [])) for r in results)
                print(
                    f"  ✓ Pedigree analysis: {pedigree_count} pedigree(s), "
                    f"{total_individuals} individuals extracted"
                )
                return summary
            else:
                print(f"  - No pedigrees detected in figures")
                return None

        except Exception as e:
            logger.warning(f"Pedigree extraction failed for PMID {pmid}: {e}")
            print(f"  - Pedigree extraction failed: {e}")
            return None

    def _run_data_scout(self, pmid: str, unified_content: str) -> bool:
        """
        Run the Genetic Data Scout to identify high-value data zones.

        Creates two additional files:
        - {PMID}_DATA_ZONES.json: Zone metadata for debugging/analysis
        - {PMID}_DATA_ZONES.md: Condensed markdown with only high-value zones

        Args:
            pmid: PubMed ID
            unified_content: Full markdown content to analyze

        Returns:
            True if scout ran successfully, False otherwise
        """
        if not SCOUT_AVAILABLE:
            return False

        if not self.gene_symbol:
            print(f"  - Skipping data scout: no gene symbol provided")
            return False

        try:
            settings = get_settings()
            if not settings.scout_enabled:
                return False

            scout = GeneticDataScout(
                gene_symbol=self.gene_symbol,
                min_relevance_score=settings.scout_min_relevance,
                max_zones=settings.scout_max_zones,
            )

            report = scout.scan(unified_content, pmid=pmid)

            # Write zone metadata JSON
            zones_json_path = self.output_dir / f"{pmid}_DATA_ZONES.json"
            zones_json_path.write_text(scout.to_json(report))

            # Write condensed markdown with high-value zones only
            zones_md_path = self.output_dir / f"{pmid}_DATA_ZONES.md"
            zones_md_path.write_text(scout.format_markdown(report, unified_content))

            print(
                f"  ✓ Data Scout: {report.zones_kept}/{report.total_zones_found} zones kept ({report.compression_ratio:.0%} of original)"
            )
            return True

        except Exception as e:
            print(f"  - Data scout failed: {e}")
            return False

    @staticmethod
    def _markdown_missing_body(markdown: str) -> bool:
        """
        Detect whether the converted markdown is missing main body text.

        Returns True when only the abstract is present or the content
        is suspiciously short, signalling that we should try an HTML scrape.
        """
        if not markdown:
            return True
        headings = [line for line in markdown.splitlines() if line.startswith("### ")]
        non_abstract = [h for h in headings if not h.startswith("### Abstract")]
        return (len(non_abstract) == 0) or (len(markdown) < 2000)

    def get_supplemental_files(
        self, pmcid: str, pmid: str, doi: str
    ) -> List[Dict[str, Any]]:
        """
        Orchestrates fetching supplemental files, first via API, then by scraping.

        Args:
            pmcid: PubMed Central ID
            pmid: PubMed ID (for logging)
            doi: Digital Object Identifier (for fallback scraping)

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys,
            plus 'base_url' and 'original_url' for PMC supplements
        """
        pmc_base_url = (
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" if pmcid else None
        )

        # 1. First, try the EuropePMC API
        api_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        try:
            response = self.session.get(api_url, timeout=30)
            response.raise_for_status()

            # Check if response body is not empty
            if not response.text.strip():
                raise ValueError("Empty response from API")

            content_type = response.headers.get("Content-Type", "").lower()

            # Try to parse as JSON first
            if "application/json" in content_type:
                data = response.json()
                # The API can return a success response with an empty list of files
                if data.get("result", {}).get("supplementaryFiles"):
                    print("  ✓ Found supplemental files via EuropePMC API (JSON)")
                    supp_files = data["result"]["supplementaryFiles"]
                    # Add PMC base_url to enable URL variant generation
                    if pmc_base_url:
                        for f in supp_files:
                            f["base_url"] = pmc_base_url
                            f["original_url"] = f.get("url", "")
                    return supp_files
            # Try to parse as XML if JSON fails
            elif "application/xml" in content_type or "text/xml" in content_type:
                try:
                    from xml.etree import ElementTree as ET

                    root = ET.fromstring(response.text)
                    # Look for supplementary file elements in the XML
                    supp_files = []
                    # Common XML paths for supplementary files
                    for elem in root.iter():
                        if (
                            "supplement" in elem.tag.lower()
                            or "supplementary" in elem.tag.lower()
                        ):
                            url = elem.get("url") or elem.get("href") or elem.text
                            name = (
                                elem.get("name")
                                or elem.get("title")
                                or elem.get("filename")
                            )
                            if url and name:
                                file_info = {"url": url, "name": name}
                                if pmc_base_url:
                                    file_info["base_url"] = pmc_base_url
                                    file_info["original_url"] = url
                                supp_files.append(file_info)
                    if supp_files:
                        print(
                            f"  ✓ Found {len(supp_files)} supplemental files via EuropePMC API (XML)"
                        )
                        return supp_files
                    # If XML parsing succeeded but no files found, continue to fallback
                except ET.ParseError as e:
                    raise ValueError(f"Could not parse XML response: {e}")
                # If we get here, XML was parsed but no files found - continue to fallback
            else:
                raise ValueError(
                    f"Expected JSON or XML but got Content-Type: {content_type}"
                )
        except requests.exceptions.RequestException as e:
            print(
                f"  - EuropePMC supplemental files API request failed for {pmcid}: {e}"
            )
        except ValueError as e:
            print(
                f"  - EuropePMC supplemental files API returned invalid response for {pmcid}: {e}"
            )
        except Exception as e:
            print(f"  - EuropePMC supplemental files API failed for {pmcid}: {e}")

        # 2. If API fails, try scraping PMC page directly
        if pmcid:
            try:
                pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                print(f"  - Trying to scrape PMC page directly: {pmc_url}")
                pmc_response = self.session.get(
                    pmc_url, timeout=30, allow_redirects=True
                )
                pmc_response.raise_for_status()

                # Use the final URL after redirects as the base URL
                final_pmc_url = pmc_response.url

                # Use generic scraper for PMC pages
                pmc_files = self.scraper.scrape_generic_supplements(
                    pmc_response.text, final_pmc_url
                )
                if pmc_files:
                    print(
                        f"  ✓ Found {len(pmc_files)} supplemental files via PMC page scraping"
                    )
                    return pmc_files
            except Exception as e:
                print(f"  - PMC page scraping failed: {e}")

        # 3. If PMC scraping fails, fall back to DOI-based scraping
        if doi:
            print(f"  - Falling back to DOI scraping for {doi}")
            return self.doi_resolver.resolve_and_scrape_supplements(
                doi, pmid, self.scraper
            )
        else:
            print("  - API failed and no DOI available. Cannot scrape.")
            pmc_url = (
                f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                if pmcid
                else "N/A"
            )
            self._log_paywalled(pmid, "Supplemental files API failed, no DOI", pmc_url)
            return []

    def download_supplement(
        self,
        url: str,
        output_path: Path,
        pmid: str,
        filename: str,
        base_url: str = None,
        original_url: str = None,
    ) -> bool:
        """
        Download a supplemental file, trying multiple URL variants for PMC supplements.

        Args:
            url: URL of the supplement file
            output_path: Path to save the file
            pmid: PubMed ID (for logging)
            filename: Name of the file (for logging)
            base_url: Base URL of the article page (for PMC URL variant generation)
            original_url: Original URL before normalization (for PMC URL variant generation)

        Returns:
            True if download succeeded, False otherwise
        """
        # Generate URL variants for PMC supplements
        urls_to_try = [url]
        if base_url and (
            "ncbi.nlm.nih.gov/pmc/" in base_url or "pmc.ncbi.nlm.nih.gov/" in base_url
        ):
            # Use the scraper's method to generate URL variants
            source_url = original_url if original_url else url
            variants = self.scraper.get_pmc_supplement_url_variants(
                source_url, base_url
            )
            # Put primary URL first, then add unique variants
            urls_to_try = [url] + [v for v in variants if v != url]

        last_error = None
        for try_url in urls_to_try:
            try:
                # Use the session to download, which includes our headers
                response = self.session.get(
                    try_url, timeout=60, stream=True, allow_redirects=True
                )
                response.raise_for_status()

                # Collect content to validate before writing
                content_chunks = []
                for chunk in response.iter_content(chunk_size=8192):
                    content_chunks.append(chunk)

                content = b"".join(content_chunks)

                # Validate the downloaded content
                ext = output_path.suffix.lower()

                # Check for HTML error pages disguised as other file types
                # HTML pages typically start with <!DOCTYPE, <html, or have these early in the content
                content_start = (
                    content[:4096].lower() if len(content) >= 4096 else content.lower()
                )
                is_html_page = (
                    content_start.startswith(b"<!doctype")
                    or content_start.startswith(b"<html")
                    or b"<!doctype html" in content_start
                    or b"<html" in content_start[:500]
                    or
                    # Some servers return XHTML
                    b"<?xml" in content_start[:100]
                    and b"<html" in content_start[:500]
                )

                # Detect specific access denial patterns in HTML
                access_denied_patterns = [
                    # HTTP error codes
                    b"403 forbidden",
                    b"401 unauthorized",
                    b"access denied",
                    # Paywall messages
                    b"subscription required",
                    b"purchase this article",
                    b"institutional access",
                    b"sign in to access",
                    b"log in to access",
                    b"login required",
                    # Content not available
                    b"pdf not available",
                    b"full text not available",
                    b"content not available",
                    b"article not found",
                    # ScienceDirect specific
                    b"sciencedirect",
                    b"elsevier",
                    b"get access",
                    # Wiley specific
                    b"wiley online library",
                    # Generic paywall indicators
                    b"buy this article",
                    b"rent this article",
                    b"get full access",
                    b"view full text",
                ]
                is_access_denied = is_html_page and any(
                    pattern in content_start for pattern in access_denied_patterns
                )

                # Check for specific publisher error page signatures
                is_publisher_error = is_html_page and (
                    b"error" in content_start[:200]
                    or b"not found" in content_start[:500]
                    or b"unavailable" in content_start[:500]
                )

                # PDF validation: PDFs should start with %PDF
                if ext == ".pdf":
                    if not content.startswith(b"%PDF"):
                        if is_access_denied:
                            last_error = (
                                f"{filename}: ACCESS DENIED (paywall/login required)"
                            )
                            print(
                                f"    ⚠ {filename}: Access denied - likely paywall or login required"
                            )
                            continue  # Try next URL variant
                        elif is_publisher_error:
                            last_error = f"{filename}: Publisher error page received instead of PDF"
                            print(
                                f"    ⚠ {filename}: Received publisher error page instead of PDF"
                            )
                            continue  # Try next URL variant
                        elif is_html_page:
                            last_error = f"{filename}: HTML page received instead of PDF (file may be paywalled)"
                            print(
                                f"    ⚠ {filename}: Received HTML page instead of PDF - file may be paywalled or URL expired"
                            )
                            continue  # Try next URL variant
                        else:
                            # Content doesn't start with %PDF but isn't HTML either
                            # Could be corrupted or wrong content type
                            content_preview = (
                                content[:50].hex() if len(content) > 0 else "(empty)"
                            )
                            last_error = f"{filename}: Invalid PDF (starts with: {content_preview[:20]}...)"
                            print(
                                f"    ⚠ {filename}: Invalid PDF file (wrong magic bytes)"
                            )
                            continue  # Try next URL variant

                # For other file types, check if we got an HTML error page
                elif ext in [".docx", ".xlsx", ".xls", ".doc", ".zip"]:
                    if is_access_denied:
                        last_error = (
                            f"{filename}: ACCESS DENIED - file appears paywalled"
                        )
                        print(
                            f"    ⚠ {filename}: Access denied - file may be paywalled"
                        )
                        continue  # Try next URL variant
                    elif is_html_page:
                        last_error = (
                            f"{filename}: HTML page received instead of {ext} file"
                        )
                        print(
                            f"    ⚠ {filename}: Received HTML page instead of {ext} file - may be paywalled"
                        )
                        continue  # Try next URL variant

                # Write the validated content
                with open(output_path, "wb") as f:
                    f.write(content)

                if try_url != url:
                    print(f"    ✓ Downloaded from alternate URL: {try_url}")
                return True

            except Exception as e:
                last_error = str(e)
                continue  # Try next URL variant

        # All URL variants failed
        print(f"    Error downloading {filename}: {last_error}")
        # Log failed supplemental download
        self._log_paywalled(pmid, f"Supplemental file download failed: {filename}", url)
        return False

    def _process_supplements(self, pmid: str, pmcid: str, doi: str) -> Tuple[str, int]:
        """
        Download supplemental files and convert them to markdown.

        Also extracts figures from PDF supplements when extract_figures is enabled.

        Args:
            pmid: PubMed ID
            pmcid: PubMed Central ID
            doi: Digital Object Identifier for fallback scraping

        Returns:
            Tuple containing the combined supplement markdown and number of downloaded files
        """
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)

        # Check if figure extraction is enabled
        extract_figures = False
        try:
            if get_settings is not None:
                settings = get_settings()
                extract_figures = settings.extract_figures
        except Exception:
            pass

        # Create figures directory if extracting
        figures_dir = None
        if extract_figures:
            figures_dir = self.output_dir / f"{pmid}_figures"
            figures_dir.mkdir(exist_ok=True)

        supp_files = self.get_supplemental_files(pmcid, pmid, doi)
        print(f"  Found {len(supp_files)} supplemental files")

        supplement_markdown = ""
        downloaded_count = 0
        total_figures_extracted = 0

        for idx, supp in enumerate(supp_files, 1):
            url = supp.get("url", "")
            filename = supp.get("name", f"supplement_{idx}")
            # Get PMC-specific URL info for trying multiple variants
            base_url = supp.get("base_url")
            original_url = supp.get("original_url")

            if not url:
                continue

            file_path = supplements_dir / filename
            print(f"    Downloading: {filename}")

            if self.download_supplement(
                url, file_path, pmid, filename, base_url, original_url
            ):
                downloaded_count += 1

                ext = file_path.suffix.lower()
                supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"

                if ext in [".xlsx", ".xls"]:
                    supplement_markdown += self.converter.excel_to_markdown(file_path)
                elif ext == ".docx":
                    supplement_markdown += self.converter.docx_to_markdown(file_path)
                elif ext == ".doc":
                    supplement_markdown += self.converter.doc_to_markdown(file_path)
                elif ext == ".pdf":
                    # Use image extraction if enabled
                    if extract_figures and figures_dir:
                        text, images = self.converter.pdf_to_markdown_with_images(
                            file_path,
                            output_dir=figures_dir,
                        )
                        supplement_markdown += text
                        if images:
                            total_figures_extracted += len(images)
                            logger.info(
                                f"Extracted {len(images)} figures from {filename}"
                            )
                    else:
                        supplement_markdown += self.converter.pdf_to_markdown(file_path)
                elif ext in [".txt", ".csv"]:
                    try:
                        text = file_path.read_text(encoding="utf-8", errors="ignore")
                        supplement_markdown += text + "\n\n"
                    except Exception as e:
                        supplement_markdown += f"[Error reading text file: {e}]\n\n"
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"

            time.sleep(0.5)

        if total_figures_extracted > 0:
            print(
                f"  ✓ Extracted {total_figures_extracted} figures from PDF supplements"
            )

        return supplement_markdown, downloaded_count

    def download_pmid(self, pmid: str) -> Tuple[bool, str, Optional[str]]:
        """
        Download content for a single PMID and run cheap analysis.

        This method handles the download phase:
        - PMID to PMCID conversion
        - Full-text XML download
        - Figure extraction from PMC
        - Supplement download and conversion
        - Create unified markdown file
        - Run data scout (cheap text analysis)

        NO pedigree extraction - that's expensive (GPT-4o vision) and happens
        in run_post_processing() after all downloads complete.

        Args:
            pmid: PubMed ID to process

        Returns:
            Tuple of (success: bool, result: str, unified_content: Optional[str])
            - success: Whether download completed successfully
            - result: Output file path on success, or error message on failure
            - unified_content: The full markdown content (for passing to post-processing)
        """
        print(f"\nProcessing PMID: {pmid}")

        # Get DOI first (needed for supplement fallback and free text retrieval)
        doi = self.pmc_api.get_doi_from_pmid(pmid)
        if doi:
            print(f"  ✓ DOI: {doi}")
        else:
            print("  - No DOI found for this PMID.")

        # Convert PMID to PMCID
        pmcid = self.pmc_api.pmid_to_pmcid(pmid)

        if not pmcid:
            # No PMCID - check if this is a free full text article via publisher
            print(f"  - No PMCID found, checking for free full text via publisher...")
            return self._download_free_text_pmid(pmid, doi)

        print(f"  ✓ PMCID: {pmcid}")

        # Get full-text XML from PMC
        xml_content = self.pmc_api.get_fulltext_xml(pmcid)

        if not xml_content:
            print(f"  ❌ Full-text not available from PMC")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            self._log_paywalled(pmid, "Full-text not available", pmc_url)
            return False, "No full-text", None

        print(f"  ✓ Full-text XML retrieved from PMC")

        # Extract figures from the PMC article page
        self._extract_pmc_figures(pmcid, pmid)

        # Convert main text to markdown
        main_markdown = self.converter.xml_to_markdown(xml_content)

        # Some PMC entries block XML full text (only abstract). Fall back to scraping
        # the PMC HTML page when the markdown looks incomplete.
        if self._markdown_missing_body(main_markdown):
            pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
            try:
                html_response = self.session.get(
                    pmc_url, timeout=30, allow_redirects=True
                )
                html_response.raise_for_status()
                html_markdown = self.converter.pmc_html_to_markdown(html_response.text)
                if html_markdown and len(html_markdown) > len(main_markdown):
                    main_markdown = html_markdown
                    print(
                        "  ✓ Recovered full text from PMC HTML (XML body unavailable)"
                    )
            except Exception as e:
                print(f"  - PMC HTML fallback failed: {e}")

        supplement_markdown, downloaded_count = self._process_supplements(
            pmid, pmcid, doi
        )

        # Create unified markdown file (WITHOUT pedigree extraction)
        unified_content = main_markdown + supplement_markdown

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(unified_content)

        print(f"  ✅ Downloaded: {output_file.name} ({downloaded_count} supplements)")

        # Check for supplement gaps (Task 5: post-download validation)
        gap_result = check_supplement_gap(
            text=main_markdown,
            downloaded_count=downloaded_count,
            extracted_variant_count=0,  # Will be filled in after extraction
        )
        if gap_result["has_gap"]:
            for warning in gap_result["warnings"]:
                print(f"  ⚠️  PMID {pmid}: {warning}")
                logger.warning(f"PMID {pmid}: {warning}")

        # Log success
        with open(self.success_log, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([pmid, pmcid, downloaded_count])

        self._write_pmid_status(pmid, "extracted", {
            "download_timestamp": datetime.datetime.now().isoformat(),
            "variant_count": 0, #FIXME: Add actual variant count
            "source": "pmc"
        })

        return True, str(output_file), unified_content

    def _download_free_text_pmid(
        self, pmid: str, doi: str
    ) -> Tuple[bool, str, Optional[str]]:
        """
        Download content for a PMID that has no PMCID but may have free full text via publisher.

        This handles download + data scout (cheap text analysis).
        NO pedigree extraction - that happens in run_post_processing().

        Args:
            pmid: PubMed ID to process
            doi: DOI for the article (may be None)

        Returns:
            Tuple of (success: bool, result: str, unified_content: Optional[str])
        """
        # Check if article is marked as free full text
        is_free, free_url = self.pmc_api.is_free_full_text(pmid)

        if not is_free:
            # Not marked as free - try Unpaywall as last resort
            print(f"  - No PMCID and not marked as free, trying Unpaywall...")
            if doi:
                unpaywall_result, unpaywall_error = self.unpaywall.find_open_access(doi)
                if unpaywall_result and unpaywall_result.get("pdf_url"):
                    print(f"  ✓ Unpaywall found OA version: {unpaywall_result.get('oa_status')}")
                    # Download the PDF and convert to markdown
                    pdf_url = unpaywall_result["pdf_url"]
                    pdf_path = self.output_dir / f"{pmid}_unpaywall.pdf"
                    success, dl_error = self.unpaywall.download_pdf(pdf_url, str(pdf_path))
                    if success:
                        # Convert PDF to markdown
                        main_markdown = self.converter.pdf_to_markdown(str(pdf_path))
                        if main_markdown and len(main_markdown) > 500:
                            print(f"  ✓ Retrieved via Unpaywall ({len(main_markdown)} chars)")
                            # Continue with the rest of the flow
                            is_free = True
                            free_url = pdf_url
                        else:
                            print(f"  - Unpaywall PDF conversion failed or content too short")
                    else:
                        print(f"  - Unpaywall download failed: {dl_error}")
                elif unpaywall_result and unpaywall_result.get("landing_page"):
                    print(f"  - Unpaywall found landing page but no direct PDF")
                else:
                    print(f"  - Unpaywall: {unpaywall_error or 'No OA version found'}")
            
            if not is_free:
                print(f"  ❌ No PMCID and not available via any method (likely paywalled)")
                pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                self._log_paywalled(pmid, "No PMCID found, not free full text, Unpaywall failed", pubmed_url)
                self._write_pmid_status(pmid, "paywalled", {
                    "download_timestamp": datetime.datetime.now().isoformat(),
                    "failure_reason": "No PMCID found, not free full text, Unpaywall failed",
                    "source": "unpaywall"
                })
                return False, "No PMCID", None

        print(f"  ✓ Article marked as free full text on PubMed")

        # Initialize variables for tracking content source
        main_markdown = None
        final_url = None
        supp_files = []
        used_elsevier_api = False
        used_wiley_api = False

        # Track whether we need to try web scraping as fallback
        api_insufficient_content = False

        # Try Elsevier API first if this is an Elsevier article
        if doi and self.elsevier_api.is_elsevier_doi(doi):
            elsevier_markdown, elsevier_error = self._try_elsevier_api(doi, pmid)
            if elsevier_markdown:
                main_markdown = elsevier_markdown
                final_url = f"https://doi.org/{doi}"
                used_elsevier_api = True
                # Still need to get supplements via scraping
                supp_files = self.doi_resolver.resolve_and_scrape_supplements(
                    doi, pmid, self.scraper
                )
            elif elsevier_error and "insufficient content" in elsevier_error.lower():
                api_insufficient_content = True
                print(
                    f"  → Elsevier API returned abstract only, falling back to web scraping..."
                )

        # Try Wiley API if this is a Wiley article and Elsevier didn't work
        if not main_markdown and doi and self.wiley_api.is_wiley_doi(doi):
            wiley_markdown, wiley_error = self._try_wiley_api(doi, pmid)
            if wiley_markdown:
                main_markdown = wiley_markdown
                final_url = f"https://doi.org/{doi}"
                used_wiley_api = True
                # Still need to get supplements via scraping
                supp_files = self.doi_resolver.resolve_and_scrape_supplements(
                    doi, pmid, self.scraper
                )
            elif wiley_error and "insufficient content" in wiley_error.lower():
                api_insufficient_content = True
                print(
                    f"  → Wiley API returned abstract only, falling back to web scraping..."
                )

        # Try Springer API if this is a Springer/Nature/BMC article
        if not main_markdown and doi and self.springer_api.is_springer_doi(doi):
            springer_markdown, springer_error = self._try_springer_api(doi, pmid)
            if springer_markdown:
                main_markdown = springer_markdown
                final_url = f"https://doi.org/{doi}"
                # Still need to get supplements via scraping
                supp_files = self.doi_resolver.resolve_and_scrape_supplements(
                    doi, pmid, self.scraper
                )
            elif springer_error and "not openaccess" in springer_error.lower():
                print(
                    f"  → Springer API: article not OpenAccess, falling back to web scraping..."
                )

        # Fall back to DOI resolver if publisher APIs didn't work or returned insufficient content
        if not main_markdown and doi:
            # Try to fetch full text and supplements from publisher using DOI
            main_markdown, final_url, supp_files = (
                self.doi_resolver.resolve_and_fetch_fulltext(doi, pmid, self.scraper)
            )

            # Validate content quality (not just length)
            if main_markdown:
                is_valid, reason = self._validate_content_quality(
                    main_markdown, f"https://doi.org/{doi}"
                )
                if not is_valid:
                    print(f"  ⚠ DOI content validation failed: {reason}")
                    print(f"  ❌ Skipping - does not contain valid full article text")
                    self._log_paywalled(
                        pmid,
                        f"DOI content validation failed: {reason}",
                        f"https://doi.org/{doi}",
                    )
                    main_markdown = None  # Clear invalid content

        # If no DOI or DOI-based methods failed, try free URL
        if not main_markdown and free_url:
            # No DOI, but we have a direct URL to the free full text
            parsed_url = urlparse(free_url)
            if parsed_url.netloc in self.SUSPICIOUS_FREE_URL_DOMAINS:
                print(
                    f"  - Skipping suspicious free URL on {parsed_url.netloc} (likely non-article content)"
                )
                self._log_paywalled(pmid, "Suspicious free URL skipped", free_url)
                return False, "Suspicious free URL", None

            print(f"  - No DOI, attempting to fetch from free URL: {free_url}")

            # Check if this is an Elsevier URL and try API first
            if self.elsevier_api.is_available and self.elsevier_api.is_elsevier_url(
                free_url
            ):
                pii = self.elsevier_api.extract_pii_from_url(free_url)
                if pii:
                    print(f"  Trying Elsevier API (with circuit breaker) for PII: {pii}")
                    try:
                        elsevier_markdown, elsevier_error = (
                            self.elsevier_client.fetch_fulltext(pii=pii)
                        )
                        if elsevier_markdown:
                            main_markdown = elsevier_markdown
                            final_url = free_url
                            used_elsevier_api = True
                            print(
                                f"  ✓ Full text retrieved via Elsevier API ({len(main_markdown)} characters)"
                            )
                    except Exception as e:
                        if "circuit breaker" in str(e).lower():
                            print(f"  ⚠ Circuit breaker protection activated for Elsevier: {e}")
                        else:
                            print(f"  - Elsevier API failed: {e}")
                        # Get supplements via web scraping
                        try:
                            response = self.session.get(
                                free_url, allow_redirects=True, timeout=10
                            )
                            response.raise_for_status()
                            supp_files = self.scraper.scrape_elsevier_supplements(
                                response.text, response.url
                            )
                        except Exception:
                            supp_files = []

            # Check if this is a Wiley URL and try API
            if (
                not main_markdown
                and self.wiley_api.is_available
                and self.wiley_api.is_wiley_url(free_url)
            ):
                extracted_doi = self.wiley_api.extract_doi_from_url(free_url)
                if extracted_doi:
                    print(f"  Trying Wiley API (with circuit breaker) for DOI: {extracted_doi}")
                    try:
                        wiley_markdown, wiley_error = self.wiley_client.fetch_fulltext(
                            doi=extracted_doi
                        )
                        if wiley_markdown:
                            main_markdown = wiley_markdown
                            final_url = free_url
                            used_wiley_api = True
                            print(
                                f"  ✓ Full text retrieved via Wiley API ({len(main_markdown)} characters)"
                            )
                    except Exception as e:
                        if "circuit breaker" in str(e).lower():
                            print(f"  ⚠ Circuit breaker protection activated for Wiley: {e}")
                        else:
                            print(f"  - Wiley API failed: {e}")
                        # Get supplements via web scraping
                        try:
                            response = self.session.get(
                                free_url, allow_redirects=True, timeout=10
                            )
                            response.raise_for_status()
                            supp_files = self.scraper.scrape_generic_supplements(
                                response.text, response.url
                            )
                        except Exception:
                            supp_files = []

            # Fall back to web scraping if API didn't work
            if not main_markdown:
                try:
                    response = self.session.get(
                        free_url, allow_redirects=True, timeout=10
                    )
                    response.raise_for_status()
                    final_url = response.url
                    print(f"  ✓ Retrieved free full text page")

                    # Handle Elsevier linkinghub redirects
                    domain = urlparse(final_url).netloc
                    if "linkinghub.elsevier.com" in domain:
                        try:
                            pii_match = re.search(r"/pii/([^/?]+)", final_url)
                            if pii_match:
                                pii = pii_match.group(1)
                                sciencedirect_url = f"https://www.sciencedirect.com/science/article/pii/{pii}"
                                print(
                                    f"  → Attempting to access ScienceDirect page: {sciencedirect_url}"
                                )
                                redirect_response = self.session.get(
                                    sciencedirect_url, allow_redirects=True, timeout=30
                                )
                                redirect_response.raise_for_status()
                                final_url = redirect_response.url
                                response = redirect_response
                                domain = urlparse(final_url).netloc
                        except Exception as e:
                            print(f"  - Could not follow redirect from linkinghub: {e}")

                    # Extract full text and supplements from the page
                    html_content = response.text
                    main_markdown, title = self.scraper.extract_fulltext(
                        html_content, final_url
                    )

                    # Validate content quality (not just length, but actual paper content)
                    if main_markdown:
                        is_valid, reason = self._validate_content_quality(
                            main_markdown, final_url
                        )
                        if is_valid:
                            print(
                                f"  ✓ Extracted full text ({len(main_markdown)} characters)"
                            )
                        else:
                            print(f"  ⚠ Content validation failed: {reason}")
                            print(
                                f"  ❌ Skipping this URL as it doesn't contain valid article content"
                            )
                            self._log_paywalled(
                                pmid,
                                f"Content validation failed: {reason}",
                                free_url,
                            )
                            main_markdown = None  # Clear invalid content
                    else:
                        print(f"  ❌ Could not extract full text from page")

                    # Get supplements from the page
                    domain = urlparse(final_url).netloc

                    # Route to domain-specific scraper
                    if "nature.com" in domain:
                        supp_files = self.scraper.scrape_nature_supplements(
                            html_content, final_url
                        )
                    elif any(
                        d in domain
                        for d in ["gimjournal.org", "sciencedirect.com", "elsevier.com"]
                    ):
                        supp_files = self.scraper.scrape_elsevier_supplements(
                            html_content, final_url
                        )
                    else:
                        supp_files = self.scraper.scrape_generic_supplements(
                            html_content, final_url
                        )

                except requests.exceptions.RequestException as e:
                    print(f"  ❌ Failed to fetch free full text from {free_url}: {e}")
                    self._log_paywalled(
                        pmid, f"Free full text fetch failed: {e}", free_url
                    )
                    return False, "Free text fetch failed", None
        elif not main_markdown:
            # No DOI and no free URL - PubMed indicates "free" but lacks actionable links
            print(f"  ❌ No DOI or free URL available to fetch full text")
            print(
                f"     (PubMed metadata indicates free access but provides no usable link)"
            )
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            self._log_paywalled(
                pmid,
                "Free full text indicated but no DOI or URL in PubMed metadata",
                pubmed_url,
            )
            return False, "No DOI or URL for free text", None

        if not main_markdown:
            print(f"  ❌ Could not retrieve full text from publisher")
            fallback_url = final_url or (
                f"https://doi.org/{doi}"
                if doi
                else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            )
            self._log_paywalled(pmid, "Free full text extraction failed", fallback_url)
            return False, "Free text extraction failed", None

        print(
            f"  ✓ Full text retrieved from publisher ({len(main_markdown)} characters)"
        )

        # Create supplements directory
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)

        # Check if figure extraction is enabled
        extract_figures = False
        try:
            if get_settings is not None:
                settings = get_settings()
                extract_figures = settings.extract_figures
        except Exception:
            pass

        # Create figures directory if extracting
        figures_dir = None
        if extract_figures:
            figures_dir = self.output_dir / f"{pmid}_figures"
            figures_dir.mkdir(exist_ok=True)

        print(f"  Found {len(supp_files)} supplemental files")

        supplement_markdown = ""
        downloaded_count = 0
        total_figures_extracted = 0

        # Download and convert each supplement
        for idx, supp in enumerate(supp_files, 1):
            url = supp.get("url", "")
            filename = supp.get("name", f"supplement_{idx}")

            if not url:
                continue

            file_path = supplements_dir / filename
            print(f"    Downloading: {filename}")

            if self.download_supplement(url, file_path, pmid, filename):
                downloaded_count += 1

                ext = file_path.suffix.lower()

                supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"

                # Convert supplement to markdown based on file type
                if ext in [".xlsx", ".xls"]:
                    supplement_markdown += self.converter.excel_to_markdown(file_path)
                elif ext == ".docx":
                    supplement_markdown += self.converter.docx_to_markdown(file_path)
                elif ext == ".doc":
                    supplement_markdown += self.converter.doc_to_markdown(file_path)
                elif ext == ".pdf":
                    # Use image extraction if enabled
                    if extract_figures and figures_dir:
                        text, images = self.converter.pdf_to_markdown_with_images(
                            file_path,
                            output_dir=figures_dir,
                        )
                        supplement_markdown += text
                        if images:
                            total_figures_extracted += len(images)
                    else:
                        supplement_markdown += self.converter.pdf_to_markdown(file_path)
                elif ext in [".txt", ".csv"]:
                    try:
                        text = file_path.read_text(encoding="utf-8", errors="ignore")
                        supplement_markdown += text + "\n\n"
                    except Exception as e:
                        supplement_markdown += f"[Error reading text file: {e}]\n\n"
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"

            time.sleep(0.5)

        if total_figures_extracted > 0:
            print(
                f"  ✓ Extracted {total_figures_extracted} figures from PDF supplements"
            )

        # Create unified markdown file (WITHOUT pedigree extraction)
        unified_content = main_markdown + supplement_markdown

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(unified_content)

        if used_elsevier_api:
            source_tag = "[via Elsevier API]"
        elif used_wiley_api:
            source_tag = "[via Wiley API]"
        else:
            source_tag = "[from publisher]"
        print(
            f"  ✅ Downloaded: {output_file.name} ({downloaded_count} supplements) {source_tag}"
        )

        # Log success with special marker for publisher-sourced content
        if used_elsevier_api:
            source_marker = "ELSEVIER_API"
        elif used_wiley_api:
            source_marker = "WILEY_API"
        else:
            source_marker = "PUBLISHER_FREE"
        with open(self.success_log, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([pmid, source_marker, downloaded_count])

        return True, str(output_file), unified_content

    def run_post_processing(self, pmid: str, unified_content: str = None) -> bool:
        """
        Run post-processing on a downloaded paper (pedigree extraction only).

        This method should be called AFTER download_pmid() completes for all papers.
        It handles:
        - Pedigree extraction from figures (expensive GPT-4o vision calls)
        - Appending pedigree summary to the unified markdown file

        Note: Data scout runs during the download phase (it's cheap text analysis).

        Args:
            pmid: PubMed ID to post-process
            unified_content: Optional pre-loaded content. If not provided, reads from file.

        Returns:
            True if post-processing completed successfully, False otherwise
        """
        print(f"\n  Post-processing PMID: {pmid}")

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"

        # Load unified content if not provided
        if unified_content is None:
            if not output_file.exists():
                print(f"  ❌ Cannot post-process: {output_file.name} not found")
                return False
            unified_content = output_file.read_text(encoding="utf-8")

        # Run pedigree extraction on figures (if any were extracted)
        # This is the expensive GPT-4o vision step - batched after all downloads complete
        pedigree_summary = self._run_pedigree_extraction(pmid)
        if pedigree_summary:
            # Append pedigree summary to unified content and re-save
            unified_content += pedigree_summary
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(unified_content)
            print(f"  ✓ Updated {output_file.name} with pedigree data")

        return True

    def batch_data_scout(self, pmids: List[str]) -> Tuple[int, int]:
        """
        Run data scout on a batch of downloaded papers.

        Data scout identifies high-value data zones in the text - this is cheap
        text analysis that should run before expensive pedigree extraction.

        Args:
            pmids: List of PubMed IDs to analyze

        Returns:
            Tuple of (successful_count, failed_count)
        """
        if not pmids:
            print("  No PMIDs for data scout")
            return 0, 0

        print(f"\n{'=' * 60}")
        print(f"DATA SCOUT PHASE: {len(pmids)} papers")
        print(f"(Identifying high-value data zones)")
        print(f"{'=' * 60}")

        successful = 0
        failed = 0

        for idx, pmid in enumerate(pmids, 1):
            print(f"\n[{idx}/{len(pmids)}] Scouting PMID {pmid}...", end="")

            try:
                output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
                if not output_file.exists():
                    print(f" ❌ File not found")
                    failed += 1
                    continue

                unified_content = output_file.read_text(encoding="utf-8")
                if self._run_data_scout(pmid, unified_content):
                    successful += 1
                else:
                    # Data scout may return False if disabled or no gene symbol
                    successful += 1  # Still count as processed
            except Exception as e:
                print(f" ❌ Failed: {e}")
                failed += 1

        print(f"\n  Data scout complete: {successful} succeeded, {failed} failed")
        return successful, failed

    def _has_pedigree_indicators(self, pmid: str) -> bool:
        """
        Check if a paper's figure legends suggest pedigree content.

        Scans the FULL_CONTEXT.md for figure-related text containing
        pedigree indicators like "pedigree", "family", "proband", etc.

        Args:
            pmid: PubMed ID to check

        Returns:
            True if pedigree-related figure content is likely present
        """
        # Keywords that suggest pedigree figures
        pedigree_keywords = [
            "pedigree",
            "family tree",
            "proband",
            "affected",
            "carrier",
            "unaffected",
            "heterozygous",
            "homozygous",
            "kindred",
            "index patient",
            "index case",
            "familial",
            "inheritance",
            "autosomal dominant",
            "autosomal recessive",
            "segregat",
            "generation",
            "siblings",
            "offspring",
            "consanguineous",
        ]

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        if not output_file.exists():
            return False

        try:
            content = output_file.read_text(encoding="utf-8").lower()

            # Look for figure legends/captions containing pedigree keywords
            # Common patterns: "Figure 1.", "Fig. 1:", "**Figure 1**", etc.
            import re

            # Find all figure caption regions (figure reference + next ~500 chars)
            figure_pattern = r"(fig\.?(?:ure)?\.?\s*\d+[^a-z])"
            matches = list(re.finditer(figure_pattern, content, re.IGNORECASE))

            for match in matches:
                # Get text around the figure reference (caption area)
                start = match.start()
                end = min(start + 500, len(content))
                caption_region = content[start:end]

                # Check for pedigree keywords in this region
                for keyword in pedigree_keywords:
                    if keyword in caption_region:
                        return True

            # Also check for standalone pedigree mentions that might indicate figures
            if "pedigree" in content and ("figure" in content or "fig." in content):
                return True

        except Exception as e:
            logger.warning(f"Error checking pedigree indicators for {pmid}: {e}")
            return True  # Err on the side of processing if we can't check

        return False

    def batch_pedigree_extraction(self, pmids: List[str]) -> Tuple[int, int]:
        """
        Run pedigree extraction on a batch of downloaded papers.

        This method filters papers in two stages:
        1. Must have figures extracted
        2. Figure legends must suggest pedigree content (keyword scan)

        This dramatically reduces expensive GPT-4o vision API calls.
        Should be called after data scout completes.

        Args:
            pmids: List of PubMed IDs to process

        Returns:
            Tuple of (successful_count, failed_count)
        """
        if not pmids:
            print("  No PMIDs for pedigree extraction")
            return 0, 0

        # Stage 1: Filter to only papers that have figures
        pmids_with_figures = []
        for pmid in pmids:
            figures_dir = self.output_dir / f"{pmid}_figures"
            if figures_dir.exists():
                image_files = list(figures_dir.glob("*.png")) + list(
                    figures_dir.glob("*.jpg")
                )
                if image_files:
                    pmids_with_figures.append(pmid)

        no_figures_count = len(pmids) - len(pmids_with_figures)

        # Stage 2: Filter by figure legend content (pedigree keywords)
        print(
            f"\n  Scanning {len(pmids_with_figures)} papers for pedigree indicators..."
        )
        pmids_likely_pedigree = []
        for pmid in pmids_with_figures:
            if self._has_pedigree_indicators(pmid):
                pmids_likely_pedigree.append(pmid)

        no_indicators_count = len(pmids_with_figures) - len(pmids_likely_pedigree)

        print(f"\n{'=' * 60}")
        print(f"PEDIGREE EXTRACTION PHASE")
        print(f"  {len(pmids_likely_pedigree)} papers likely have pedigrees")
        print(
            f"  (Skipped: {no_figures_count} no figures, {no_indicators_count} no pedigree keywords)"
        )
        print(f"  (GPT-4o vision calls)")
        print(f"{'=' * 60}")

        if not pmids_likely_pedigree:
            print("  No papers appear to have pedigree figures")
            return 0, 0

        successful = 0
        failed = 0

        for idx, pmid in enumerate(pmids_likely_pedigree, 1):
            print(f"\n[{idx}/{len(pmids_likely_pedigree)}]", end="")

            try:
                if self.run_post_processing(pmid):
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"  ❌ Pedigree extraction failed for {pmid}: {e}")
                failed += 1

        print(
            f"\n  Pedigree extraction complete: {successful} succeeded, {failed} failed"
        )
        return successful, failed

    def batch_post_process(self, pmids: List[str]) -> Tuple[int, int]:
        """
        Run full post-processing on a batch of downloaded papers.

        This runs both data scout AND pedigree extraction in sequence.
        For more control, use batch_data_scout() and batch_pedigree_extraction() separately.

        Args:
            pmids: List of PubMed IDs to post-process

        Returns:
            Tuple of (successful_count, failed_count) for pedigree extraction
        """
        if not pmids:
            print("  No PMIDs to post-process")
            return 0, 0

        # Phase 1: Data scout (cheap text analysis)
        self.batch_data_scout(pmids)

        # Phase 2: Pedigree extraction (expensive GPT-4o)
        return self.batch_pedigree_extraction(pmids)

    def process_pmid(self, pmid: str) -> Tuple[bool, str]:
        """
        Process a single PMID: convert to PMCID, download content, create unified markdown.

        For papers without PMCIDs but marked as "Free Full Text" on PubMed, this method
        will attempt to fetch full text directly from the publisher's website.

        This is the original combined method - downloads AND post-processes in one pass.
        Kept for backward compatibility. For batch processing, use download_pmid() +
        batch_post_process() instead.

        Args:
            pmid: PubMed ID to process

        Returns:
            Tuple of (success: bool, result: str) where result is output file path or error message
        """
        # Phase 1: Download
        success, result, unified_content = self.download_pmid(pmid)

        if not success:
            return False, result

        # Phase 2: Post-process (pedigree extraction + data scout)
        self.run_post_processing(pmid, unified_content)

        return True, result

    def _try_elsevier_api(
        self, doi: str, pmid: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to fetch full text via Elsevier API if applicable.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Full text as markdown if successful, None otherwise
            - error_message: Error description if failed, None if successful
        """
        if not self.elsevier_api.is_available:
            return None, "Elsevier API key not configured"

        if not doi or not self.elsevier_api.is_elsevier_doi(doi):
            return None, "Not an Elsevier DOI"

        print(f"  Trying Elsevier API (with circuit breaker) for DOI: {doi}")
        
        try:
            markdown, error = self.elsevier_client.fetch_fulltext(doi=doi)

            if markdown:
                print(
                    f"  ✓ Full text retrieved via Elsevier API ({len(markdown)} characters)"
                )
                return markdown, None
            else:
                print(f"  - Elsevier API: {error}")
                return None, error
        except Exception as e:
            if "circuit breaker" in str(e).lower():
                print(f"  ⚠ Circuit breaker protection activated: {e}")
            return None, str(e)

    def _try_wiley_api(
        self, doi: str, pmid: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to fetch full text via Wiley TDM API if applicable.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Full text as markdown if successful, None otherwise
            - error_message: Error description if failed, None if successful
        """
        if not self.wiley_api.is_available:
            return None, "Wiley API key not configured"

        if not doi or not self.wiley_api.is_wiley_doi(doi):
            return None, "Not a Wiley DOI"

        print(f"  Trying Wiley API (with circuit breaker) for DOI: {doi}")
        
        try:
            markdown, error = self.wiley_client.fetch_fulltext(doi=doi)

            if markdown:
                print(f"  ✓ Full text retrieved via Wiley API ({len(markdown)} characters)")
                return markdown, None
            else:
                print(f"  - Wiley API: {error}")
                return None, error
        except Exception as e:
            if "circuit breaker" in str(e).lower():
                print(f"  ⚠ Circuit breaker protection activated: {e}")
            return None, str(e)

    def _try_springer_api(
        self, doi: str, pmid: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to fetch full text via Springer Nature OpenAccess API if applicable.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Full text as markdown if successful, None otherwise
            - error_message: Error description if failed, None if successful
        """
        if not self.springer_api.is_available:
            return None, "Springer API key not configured"

        if not doi or not self.springer_api.is_springer_doi(doi):
            return None, "Not a Springer/Nature/BMC DOI"

        print(f"  Trying Springer OpenAccess API (with circuit breaker) for DOI: {doi}")
        
        try:
            markdown, metadata, error = self.springer_client.fetch_article(doi)

            if markdown:
                print(
                    f"  ✓ Full text retrieved via Springer API ({len(markdown)} characters)"
                )
                return markdown, None
            elif metadata:
                # Got metadata but no full text - article may not be open access
                print(f"  - Springer API: Full text unavailable (not OpenAccess)")
                return None, "Article not available in OpenAccess"
            else:
                print(f"  - Springer API: {error}")
                return None, error
        except Exception as e:
            if "circuit breaker" in str(e).lower():
                print(f"  ⚠ Circuit breaker protection activated: {e}")
            return None, str(e)

    def _process_free_text_pmid(self, pmid: str, doi: str) -> Tuple[bool, str]:
        """
        Process a PMID that has no PMCID but may have free full text via publisher.

        This method checks if the article is marked as "Free Full Text" on PubMed
        and attempts to fetch the content from the publisher's website.
        For Elsevier articles, it first tries the Elsevier API if an API key is configured.

        This is the original combined method - kept for backward compatibility.
        Internally calls _download_free_text_pmid() + run_post_processing().

        Args:
            pmid: PubMed ID to process
            doi: DOI for the article (may be None)

        Returns:
            Tuple of (success: bool, result: str) where result is output file path or error message
        """
        # Phase 1: Download
        success, result, unified_content = self._download_free_text_pmid(pmid, doi)

        if not success:
            return False, result

        # Phase 2: Post-process (pedigree extraction + data scout)
        self.run_post_processing(pmid, unified_content)

        return True, result

    def harvest(
        self,
        pmids: List[str],
        delay: float = 2.0,
        run_scout: bool = False,
        manifest_path: Optional[str] = None,
    ):
        """
        Harvest full-text and supplements for a list of PMIDs.

        By default, only performs downloads. Post-processing (data scout and pedigree
        extraction) can be enabled via run_scout flag or called separately via
        batch_data_scout() and batch_pedigree_extraction().

        Args:
            pmids: List of PubMed IDs to process
            delay: Delay in seconds between processing each PMID
            run_scout: If True, run data scout and pedigree extraction after downloads.
                       Default False - download only. Call batch_post_process() separately
                       for post-processing.
            manifest_path: Path to write manifest.json. If None, defaults to
                          {output_dir}/manifest.json

        Raises:
            ValidationError: If input validation fails
        """
        # Validate inputs and filter to valid PMIDs only
        valid_pmids = validate_harvest_inputs(pmids, self.output_dir)

        if len(valid_pmids) < len(pmids):
            print(f"⚠ Filtered {len(pmids) - len(valid_pmids)} invalid PMIDs")

        print(f"Starting harvest for {len(valid_pmids)} PMIDs")
        print(f"Output directory: {self.output_dir.absolute()}\n")

        # Use validated PMIDs for processing
        pmids = valid_pmids

        # Initialize manifest for tracking download outcomes
        manifest = Manifest(stage=Stage.DOWNLOAD, gene=self.gene_symbol)
        if manifest_path is None:
            manifest_path = self.output_dir / "manifest.json"
        else:
            manifest_path = Path(manifest_path)

        # ============================================
        # PHASE 1: DOWNLOAD ALL PAPERS
        # ============================================
        print(f"{'=' * 60}")
        print(f"DOWNLOAD PHASE: {len(pmids)} papers")
        print(f"{'=' * 60}")

        downloaded_pmids = []
        download_successful = 0
        download_failed = 0

        for idx, pmid in enumerate(pmids, 1):
            print(f"[{idx}/{len(pmids)}]", end=" ")

            success, result, content = self.download_pmid(pmid)

            if success:
                download_successful += 1
                downloaded_pmids.append(pmid)
                # Track files created for this PMID
                files_created = []
                full_context = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
                if full_context.exists():
                    files_created.append(str(full_context.name))
                # Check for supplements directory
                supp_dir = self.output_dir / f"{pmid}_supplements"
                if supp_dir.exists():
                    files_created.extend(
                        [f.name for f in supp_dir.iterdir() if f.is_file()]
                    )
                # Check for figures directory
                fig_dir = self.output_dir / f"{pmid}_figures"
                if fig_dir.exists():
                    files_created.extend(
                        [f.name for f in fig_dir.iterdir() if f.is_file()]
                    )

                manifest.add_entry(
                    ManifestEntry(
                        pmid=pmid,
                        status=Status.SUCCESS,
                        files_created=files_created,
                    )
                )
            else:
                download_failed += 1
                # Map error messages to Status enum
                error_msg = result.lower() if result else ""
                if "timeout" in error_msg:
                    status = Status.TIMEOUT
                elif "captcha" in error_msg:
                    status = Status.CAPTCHA
                elif any(
                    kw in error_msg
                    for kw in ["paywall", "full-text", "access", "no pmcid"]
                ):
                    status = Status.PAYWALL
                else:
                    status = Status.FAILED

                manifest.add_entry(
                    ManifestEntry(
                        pmid=pmid,
                        status=status,
                        error_message=result,
                    )
                )

            if idx < len(pmids):
                time.sleep(delay)

        # Save manifest after download phase
        manifest.save(manifest_path)
        print(f"\n  📋 Manifest saved: {manifest_path}")

        print(f"\n{'=' * 60}")
        print(f"Download phase complete!")
        print(f"  ✅ Downloaded: {download_successful}")
        print(f"  ❌ Failed: {download_failed}")
        print(f"{'=' * 60}")

        # ============================================
        # PHASE 2: POST-PROCESS ALL DOWNLOADED PAPERS
        # (Only runs if run_scout=True)
        # ============================================
        post_successful, post_failed = 0, 0
        if run_scout and downloaded_pmids:
            post_successful, post_failed = self.batch_post_process(downloaded_pmids)
        elif downloaded_pmids and not run_scout:
            print(f"\n  ℹ️  Skipping post-processing (run_scout=False)")
            print(
                f"  To analyze downloaded papers, call batch_post_process() separately."
            )

        # ============================================
        # FINAL SUMMARY
        # ============================================
        print(f"\n{'=' * 60}")
        print(f"Harvest complete!")
        print(f"  ✅ Successful: {download_successful}")
        print(f"  ❌ Failed: {download_failed}")
        print(f"  📊 Post-processed: {post_successful}/{len(downloaded_pmids)}")
        print(f"  Output directory: {self.output_dir.absolute()}")
        print(f"  Success log: {self.success_log}")
        print(f"  Paywalled log: {self.paywalled_log}")
        print(f"{'=' * 60}")

    # Backward-compatible methods for tests and legacy code
    def pmid_to_pmcid(self, pmid: str):
        """Backward-compatible wrapper for PMC API."""
        return self.pmc_api.pmid_to_pmcid(pmid)

    def get_doi_from_pmid(self, pmid: str):
        """Backward-compatible wrapper for PMC API."""
        return self.pmc_api.get_doi_from_pmid(pmid)

    def get_fulltext_xml(self, pmcid: str):
        """Backward-compatible wrapper for PMC API."""
        return self.pmc_api.get_fulltext_xml(pmcid)

    def _get_supplemental_files_from_doi(self, doi: str, pmid: str):
        """Backward-compatible wrapper for DOI resolver."""
        return self.doi_resolver.resolve_and_scrape_supplements(doi, pmid, self.scraper)
