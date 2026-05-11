"""
Orchestrator Module

Main PMCHarvester class that coordinates all harvesting operations:
- Converts PMIDs to PMCIDs
- Downloads full-text and supplemental files
- Creates unified markdown files for LLM processing
"""

import csv
import datetime
import json
import os
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlparse

import requests
from bs4 import BeautifulSoup
from utils.bootstrap import initialize_runtime
from utils.logging_utils import get_logger

from utils.resilience import CircuitBreaker, ResilientAPIClient

from .core_api import COREAPIClient, get_core_client
from .doi_resolver import DOIResolver
from .elsevier_api import ElsevierAPIClient
from .format_converters import FormatConverter
from .free_text_fetch_service import fetch_main_content_for_free_text
from .free_text_flow import initialize_free_text_access
from .free_text_output_service import (
    source_from_free_text_flags,
    write_free_text_output,
)
from .pmc_api import PMCAPIClient
from .content_validation import validate_content_quality
from .persistence import (
    append_paywalled_entry,
    append_success_entry,
    initialize_harvest_logs,
    write_pmid_status_file,
)
from .supplement_processing_service import process_supplement_files
from .priority_queue import Priority, PriorityQueue
from .priority_queue import Status as QueueStatus
from .retry_manager import RetryConfig, RetryManager
from .springer_api import SpringerAPIClient
from .supplement_reference_parser import (
    check_supplement_gap,
    extract_supplement_urls_from_text,
    parse_supplement_references,
)
from gene_literature.supplements import UnifiedSupplementFetcher
from .artifacts_log import ArtifactsLog
from .figure_extractor import (
    CaptionExtractionResult,
    extract_from_html as extract_captions_from_html,
    extract_from_jats_xml as extract_captions_from_jats_xml,
    merge_results as merge_caption_results,
    render_captions_markdown,
    save_captions_json,
)
from .supplement_scraper import SupplementScraper
from .unpaywall_api import UnpaywallClient, get_unpaywall_client
from .wiley_api import WileyAPIClient

logger = get_logger(__name__)


# =============================================================================
# INPUT VALIDATION
# =============================================================================

# PMID pattern: 1-8 digit numeric string
PMID_PATTERN = re.compile(r"^\d{1,8}$")
ABSTRACT_ONLY_FALLBACK_MARKER = "abstract-only fallback"
THIN_FULL_CONTEXT_BYTES = 6_000
PUBLISHER_SHELL_MARKERS = (
    "eletters should relate to an article recently published",
    "comments and feedback on aha/asa scientific statements",
    "access through your institution",
    "institutional sign in",
    "please enable javascript",
    "enable cookies",
    "purchase access",
    "rent this article",
    "subscribe to this journal",
    "sign in to view",
)


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


def full_context_needs_retry(full_context_path: Path, output_dir: Path) -> bool:
    """Return True when an existing FULL_CONTEXT artifact is not usable full text."""
    if not full_context_path.exists():
        return True

    try:
        file_size = full_context_path.stat().st_size
    except OSError:
        return True

    if file_size == 0:
        return True

    try:
        with open(full_context_path, encoding="utf-8", errors="ignore") as f:
            head = f.read(12_000)
    except OSError:
        return True

    normalized = head.lower()
    if ABSTRACT_ONLY_FALLBACK_MARKER in normalized:
        return True

    if any(marker in normalized for marker in PUBLISHER_SHELL_MARKERS):
        return True

    if file_size >= THIN_FULL_CONTEXT_BYTES:
        return False

    pmid = full_context_path.name.replace("_FULL_CONTEXT.md", "")
    status_file = output_dir / "pmid_status" / f"{pmid}.json"
    if status_file.exists():
        try:
            with open(status_file, encoding="utf-8") as f:
                status_data = json.load(f)
        except (OSError, json.JSONDecodeError):
            status_data = {}
        status = str(status_data.get("status") or "").lower()
        source = str(status_data.get("source") or "").lower()
        if status == "abstract_only" or source == "pubmed_abstract":
            return True

    return False


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
    from config.settings import get_settings
    from pipeline.data_scout import GeneticDataScout

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

from utils.manifest import Manifest, ManifestEntry, Stage, Status


class PMCHarvester:
    """Harvests full-text and supplemental materials from PubMed Central."""

    SUSPICIOUS_FREE_URL_DOMAINS = {"antibodies.cancer.gov"}

    @property
    def _should_extract_figures(self) -> bool:
        """Check if figure extraction is enabled in settings."""
        try:
            if get_settings is not None:
                return get_settings().extract_figures
        except Exception:
            return False
        return False

    def __init__(self, output_dir: str = "pmc_harvest", gene_symbol: str = None):
        """
        Initialize PMC Harvester.

        Args:
            output_dir: Directory to save harvested files
            gene_symbol: Target gene symbol for data scout analysis
        """
        initialize_runtime()

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.paywalled_log = self.output_dir / "paywalled_missing.csv"
        self.success_log = self.output_dir / "successful_downloads.csv"
        self.abstract_only_log = self.output_dir / "abstract_only_fallback.csv"
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
        elsevier_insttoken = None
        wiley_api_key = None
        springer_api_key = None
        if get_settings is not None:
            try:
                settings = get_settings()
                elsevier_api_key = settings.elsevier_api_key
                elsevier_insttoken = settings.elsevier_insttoken
                wiley_api_key = settings.wiley_api_key
                springer_api_key = getattr(settings, "springer_api_key", None)
            except Exception:
                pass  # Settings validation may fail if other keys are missing

        # Fall back to environment variables if not in settings
        if springer_api_key is None:
            springer_api_key = os.environ.get("SPRINGER_API_KEY")

        self.elsevier_api = ElsevierAPIClient(
            api_key=elsevier_api_key,
            insttoken=elsevier_insttoken,
            session=self.session,
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
            log_file=self.output_dir / "retry_log.jsonl",
        )

        # Initialize circuit breakers for API resilience
        self.elsevier_circuit = CircuitBreaker(
            "elsevier", max_failures=5, reset_timeout=60
        )
        self.wiley_circuit = CircuitBreaker("wiley", max_failures=5, reset_timeout=60)
        self.springer_circuit = CircuitBreaker(
            "springer", max_failures=5, reset_timeout=60
        )

        # Create resilient API clients with circuit breaker protection
        self.elsevier_client = ResilientAPIClient(
            self.elsevier_api, self.elsevier_circuit
        )
        self.wiley_client = ResilientAPIClient(self.wiley_api, self.wiley_circuit)
        self.springer_client = ResilientAPIClient(
            self.springer_api, self.springer_circuit
        )

        # Initialize priority queue for manual acquisition
        self.priority_queue = PriorityQueue(
            self.output_dir / "manual_acquisition_queue.json"
        )

        initialize_harvest_logs(self.paywalled_log, self.success_log)

        # Tier 3.5: Browser HTML fallback (lazy — no browser spawned until used).
        # Disabled unless ENABLE_BROWSER_HTML_FALLBACK is true.
        self.browser_html = None
        try:
            from .browser_html import BrowserHTMLFetcher

            settings_obj = get_settings() if get_settings is not None else None
            self.browser_html = BrowserHTMLFetcher(
                scraper=self.scraper,
                converter=self.converter,
                session=self.session,
                output_dir=self.output_dir,
                settings=settings_obj,
                validate_content_quality=validate_content_quality,
            )
        except Exception as e:
            logger.debug(f"Browser HTML fallback not available: {e}")
            self.browser_html = None

    def _log_paywalled(
        self, pmid: str, reason: str, url: str, classification: str = ""
    ) -> None:
        """
        Log a paper to the paywalled/missing CSV.

        Args:
            pmid: PubMed ID
            reason: Why the paper couldn't be downloaded
            url: URL attempted or PubMed URL
            classification: One of PAYWALLED, CAPTCHA_BLOCKED,
                INSTITUTIONAL_ACCESS, SUPPLEMENT_ONLY, API_LIMIT, or empty.
        """
        append_paywalled_entry(
            self.paywalled_log, pmid, reason, url, classification=classification
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

    def _write_pmid_status(
        self, pmid: str, status: str, details: Dict[str, Any] = None
    ) -> None:
        """
        Write the status of a PMID to a JSON file.

        Args:
            pmid: PubMed ID
            status: Status of the PMID (downloaded, extracted, failed, paywalled, no_variants)
            details: Dictionary of details to write to the file
        """
        write_pmid_status_file(self.output_dir, pmid, status, details)

    def _extract_pmc_figures(
        self,
        pmcid: str,
        pmid: str,
        artifacts: Optional["ArtifactsLog"] = None,
    ) -> Tuple[int, "CaptionExtractionResult"]:
        """
        Extract figure images and captions from a PMC article page.

        PMC articles embed figures in the HTML rather than ship them as
        separate downloads. We:
          1. Fetch the article HTML once.
          2. Pull every ``<figure>`` block via ``figure_extractor`` so we
             have caption + image_url pairs.
          3. Download each image (existing behaviour) and persist a
             ``{PMID}_figures/captions.json`` mapping filenames to
             captions for downstream auditability.

        Args:
            pmcid: PubMed Central ID
            pmid: PubMed ID (for output directory naming)
            artifacts: Optional artifacts log to populate.

        Returns:
            Tuple of (downloaded_count, caption_extraction_result). The
            caption result is empty when nothing usable was found.
        """
        empty_captions = CaptionExtractionResult()
        if not self._should_extract_figures:
            return 0, empty_captions

        if not pmcid:
            return 0, empty_captions

        # Create figures directory
        figures_dir = self.output_dir / f"{pmid}_figures"
        figures_dir.mkdir(exist_ok=True)

        try:
            # Fetch the PMC article page
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            response = self.session.get(pmc_url, timeout=30)
            response.raise_for_status()
            html_text = response.text
        except Exception as e:
            logger.warning(f"Failed to fetch PMC article HTML for {pmcid}: {e}")
            return 0, empty_captions

        # Extract structured caption metadata from the page.
        try:
            captions_result = extract_captions_from_html(html_text)
        except Exception as e:
            logger.warning(f"Caption extraction failed for {pmcid}: {e}")
            captions_result = CaptionExtractionResult()

        # Map image_url -> caption (best effort) for filename pairing below.
        caption_by_url: Dict[str, Any] = {}
        for fig in captions_result.figures:
            if fig.image_url:
                caption_by_url[fig.image_url] = fig

        try:
            soup = BeautifulSoup(html_text, "html.parser")

            figure_urls: List[str] = []
            for img in soup.find_all("img"):
                src = img.get("src", "") or img.get("data-src", "")
                if src and "cdn.ncbi.nlm.nih.gov/pmc" in src:
                    if "logo" in src.lower() or "icon" in src.lower():
                        continue
                    if src.endswith((".jpg", ".jpeg", ".png", ".gif")):
                        figure_urls.append(src)

            if not figure_urls:
                logger.debug(f"No figures found in PMC article {pmcid}")
                return 0, captions_result

            downloaded = 0
            captions_index: List[Dict[str, Any]] = []
            ordered_captions = list(captions_result.figures)

            for i, url in enumerate(figure_urls, 1):
                try:
                    ext = ".jpg"
                    for e in [".png", ".gif", ".jpeg"]:
                        if e in url.lower():
                            ext = e
                            break

                    filename = f"fig_pmc_{i}{ext}"
                    filepath = figures_dir / filename

                    img_response = self.session.get(url, timeout=30)
                    img_response.raise_for_status()

                    content = img_response.content
                    valid = (
                        content[:4]
                        in (
                            b"\x89PNG",
                            b"\xff\xd8\xff\xe0",
                            b"\xff\xd8\xff\xe1",
                            b"GIF8",
                        )
                        or content[:4] == b"\xff\xd8\xff\xdb"
                    )

                    if not valid:
                        logger.debug(f"Skipped non-image content from {url}")
                        continue

                    filepath.write_bytes(content)
                    downloaded += 1
                    logger.debug(
                        f"Downloaded figure: {filename} ({len(content)} bytes)"
                    )

                    # Attach caption: prefer URL-keyed match; otherwise
                    # use the i-th caption as a positional fallback.
                    cap = caption_by_url.get(url)
                    if cap is None and i - 1 < len(ordered_captions):
                        cap = ordered_captions[i - 1]
                    cap_label = cap.label if cap else None
                    cap_title = cap.title if cap else None
                    cap_text = cap.text if cap else None
                    fig_id = cap.figure_id if cap else None

                    captions_index.append(
                        {
                            "filename": filename,
                            "label": cap_label,
                            "title": cap_title,
                            "text": cap_text,
                            "figure_id": fig_id,
                            "source_url": url,
                        }
                    )

                    if artifacts is not None:
                        artifacts.record_figure(
                            ArtifactsLog.figure_artifact_from_path(
                                filepath,
                                source="pmc_html",
                                source_url=url,
                                caption_label=cap_label,
                                caption_title=cap_title,
                                caption_text=cap_text,
                                figure_id=fig_id,
                            )
                        )

                except Exception as e:
                    logger.warning(f"Failed to download figure from {url}: {e}")
                    continue

            # Persist captions sidecar — gives a quick audit answer for
            # "did we capture caption N for figure N?".
            if captions_index:
                try:
                    captions_path = figures_dir / "captions.json"
                    import json as _json

                    captions_path.write_text(
                        _json.dumps(captions_index, indent=2),
                        encoding="utf-8",
                    )
                except Exception as e:
                    logger.warning(f"Failed to write captions.json for {pmid}: {e}")

            if downloaded > 0:
                print(f"  ✓ Extracted {downloaded} figures from PMC article")
                logger.info(f"Extracted {downloaded} figures from PMC article {pmcid}")

            return downloaded, captions_result
        except Exception as e:
            logger.warning(f"Failed to extract figures from PMC article {pmcid}: {e}")
            return 0, captions_result

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
                print("  - No pedigrees detected in figures")
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
            print("  - Skipping data scout: no gene symbol provided")
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
        all_supplements = []

        # 1+2. Unified fetcher: PMC (Europe PMC + XML + NCBI OA) + Elsevier API
        try:
            unified = UnifiedSupplementFetcher(timeout=30)
            supplements = unified.fetch_all(pmid, doi or "")
            if supplements:
                legacy = unified.to_legacy_format(supplements)
                # Add PMC base_url for URL variant generation if available
                if pmcid:
                    pmc_base_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                    for f in legacy:
                        f["base_url"] = pmc_base_url
                        f["original_url"] = f.get("url", "")
                logger.info(
                    f"Unified fetcher found {len(legacy)} supplements for PMID {pmid}"
                )
                print(f"  ✓ Found {len(legacy)} supplemental files via unified fetcher")
                all_supplements.extend(legacy)
        except Exception as e:
            logger.warning(f"Unified supplement fetcher failed for {pmid}: {e}")

        # 3. Fallback: PMC HTML scraping for non-OA PMC articles
        if not all_supplements and pmcid:
            pmc_html_supps = self.scrape_pmc_html_supplements(pmcid)
            if pmc_html_supps:
                logger.info(
                    "PMC HTML scraper found %s supplements for PMID %s (%s)",
                    len(pmc_html_supps),
                    pmid,
                    pmcid,
                )
                print(
                    f"  ✓ Found {len(pmc_html_supps)} supplemental files via PMC HTML"
                )
                all_supplements.extend(pmc_html_supps)

        # 4. Fallback: DOI-based web scraping (existing logic)
        if not all_supplements and doi:
            print(f"  - Falling back to DOI scraping for {doi}")
            return self.doi_resolver.resolve_and_scrape_supplements(
                doi, pmid, self.scraper
            )
        elif not all_supplements:
            print("  - API failed and no DOI available. Cannot scrape.")
            pmc_url = (
                f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                if pmcid
                else "N/A"
            )
            self._log_paywalled(
                pmid,
                "Supplemental files API failed, no DOI",
                pmc_url,
                classification="SUPPLEMENT_ONLY",
            )

        return all_supplements

    def scrape_pmc_html_supplements(self, pmcid: str) -> List[Dict[str, Any]]:
        """
        Scrape supplement links from PMC HTML when OA APIs return no files.

        This targets public-access (non-OA) PMC pages where supplement links are
        visible in HTML but not exposed by OA-specific supplement APIs.
        """
        if not pmcid:
            return []

        base_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
        found_files: List[Dict[str, Any]] = []
        seen: set[str] = set()

        try:
            # NCBI guidance: keep requests at or below 3 req/s.
            time.sleep(0.35)
            response = self.session.get(base_url, timeout=30)
            response.raise_for_status()
        except Exception as exc:
            logger.warning("PMC HTML supplement fetch failed for %s: %s", pmcid, exc)
            return []

        soup = BeautifulSoup(response.text, "html.parser")

        # Prefer supplementary sections, then fall back to all links.
        candidate_links = []
        selectors = [
            "div#mc-supplementary-materials a[href]",
            "div.supplementary-material a[href]",
            "div.suppl-data a[href]",
            'section[id*="supplement"] a[href]',
        ]
        for selector in selectors:
            candidate_links.extend(soup.select(selector))
        if not candidate_links:
            candidate_links = soup.select("a[href]")

        file_exts = (".xlsx", ".xls", ".csv", ".tsv", ".zip", ".pdf", ".doc", ".docx")
        for link in candidate_links:
            href = (link.get("href") or "").strip()
            if not href:
                continue

            link_text = link.get_text(" ", strip=True).lower()
            href_l = href.lower()
            is_suppish = any(
                token in (link_text + " " + href_l)
                for token in (
                    "supplement",
                    "supplementary",
                    "table s",
                    "appendix",
                    "/bin/",
                )
            )
            if not is_suppish and not href_l.endswith(file_exts):
                continue

            normalized = self.scraper._normalize_pmc_url(href, base_url)
            filename = Path(urlparse(normalized).path).name
            if not filename:
                continue
            if filename.lower().startswith("pmc") and "." not in filename:
                continue
            if normalized in seen:
                continue

            seen.add(normalized)
            found_files.append(
                {
                    "url": normalized,
                    "name": filename,
                    "base_url": base_url,
                    "original_url": href,
                }
            )

        return found_files

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
        self._log_paywalled(
            pmid,
            f"Supplemental file download failed: {filename}",
            url,
            classification="SUPPLEMENT_ONLY",
        )
        return False

    def _process_supplements(
        self, pmid: str, pmcid: str, doi: str
    ) -> Tuple[str, int, list]:
        """
        Download supplemental files and convert them to markdown.

        Also extracts figures from PDF supplements when extract_figures is enabled.

        Args:
            pmid: PubMed ID
            pmcid: PubMed Central ID
            doi: Digital Object Identifier for fallback scraping

        Returns:
            Tuple of ``(supplement_markdown, downloaded_count, file_results)``.
            ``file_results`` is a list of ``SupplementFileResult`` rows used to
            populate the per-PMID artifacts audit log.
        """
        extract_figures = self._should_extract_figures

        # Create figures directory if extracting
        figures_dir = None
        if extract_figures:
            figures_dir = self.output_dir / f"{pmid}_figures"
            figures_dir.mkdir(exist_ok=True)

        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supp_files = self.get_supplemental_files(pmcid, pmid, doi)
        print(f"  Found {len(supp_files)} supplemental files")
        if not supp_files:
            self._log_paywalled(
                pmid,
                "SUPPLEMENT_NOT_FOUND: checked API + HTML/DOI fallbacks",
                f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/" if pmcid else "",
                classification="SUPPLEMENT_ONLY",
            )
            return "", 0, []

        result = process_supplement_files(
            supp_files=supp_files,
            supplements_dir=supplements_dir,
            pmid=pmid,
            converter=self.converter,
            download_callback=lambda url,
            file_path,
            pmid,
            filename,
            supp: self.download_supplement(
                url,
                file_path,
                pmid,
                filename,
                supp.get("base_url"),
                supp.get("original_url"),
            ),
            extract_figures=extract_figures,
            figures_dir=figures_dir,
            logger=logger,
            sleep_seconds=0.5,
        )

        if result.total_figures_extracted > 0:
            print(
                f"  ✓ Extracted {result.total_figures_extracted} figures from PDF supplements"
            )

        return (
            result.supplement_markdown,
            result.downloaded_count,
            result.file_results,
        )

    # ------------------------------------------------------------------
    # FULL_CONTEXT assembly + artifacts logging helpers
    # ------------------------------------------------------------------

    def _record_supplement_results(
        self,
        artifacts: "ArtifactsLog",
        file_results: list,
    ) -> None:
        """Push per-supplement-file metadata into the artifacts log."""
        for fr in file_results:
            artifacts.record_supplement_dict(
                filename=fr.filename,
                path=fr.path,
                url=fr.url,
                size_bytes=fr.size_bytes,
                source=fr.source,
                description=fr.description or None,
                converted=fr.converted_chars > 0,
                converted_chars=fr.converted_chars,
                figures_extracted=fr.figures_extracted,
                nested_files=list(fr.nested_files),
                error=fr.error,
            )

    def _record_html_figure_artifacts(
        self,
        artifacts: "ArtifactsLog",
        captions: "CaptionExtractionResult",
        source_label: str,
    ) -> None:
        """Record figure captions extracted from HTML when no images were downloaded.

        These artifact rows have no on-disk filename — they are caption-only
        records useful for audit purposes (we know the paper had figure N
        with this caption, even if we didn't store the image).
        """
        for fig in captions.figures:
            artifacts.record_figure_dict(
                filename=fig.figure_id or fig.label,
                path="",
                size_bytes=0,
                source_url=fig.image_url,
                source=source_label,
                caption_label=fig.label,
                caption_title=fig.title,
                caption_text=fig.text,
                figure_id=fig.figure_id,
            )

    def _build_unified_content(
        self,
        main_markdown: str,
        captions: "CaptionExtractionResult",
        supplement_markdown: str,
    ) -> str:
        """Concatenate main text + caption block + supplement markdown.

        Captions land between main text and supplements so the variant
        scanner / LLM sees them in document order.
        """
        captions_md = render_captions_markdown(captions)
        return (main_markdown or "") + (captions_md or "") + (supplement_markdown or "")

    def _log_abstract_only(self, pmid: str, reason: str, output_path: str) -> None:
        """Append a row to the abstract-only fallback log (creates header if missing)."""
        log_path = self.abstract_only_log
        write_header = not log_path.exists() or log_path.stat().st_size == 0
        with open(log_path, "a", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow(["PMID", "Reason", "Output_Path", "Timestamp"])
            writer.writerow(
                [
                    pmid,
                    reason,
                    output_path,
                    datetime.datetime.now().isoformat(),
                ]
            )

    def _build_abstract_only_fallback(
        self, pmid: str, reason: str
    ) -> Tuple[bool, str, Optional[str]]:
        """Last-resort fallback: build a minimal FULL_CONTEXT.md from the abstract.

        Called only after all full-text harvest paths (PMC XML, PMC HTML,
        publisher APIs, Tier 3.5 browser fallback) have failed. The resulting
        file is clearly marked as "ABSTRACT-ONLY FALLBACK" in its header so
        downstream extraction can distinguish it from full-text successes.
        Logs to ``abstract_only_fallback.csv`` and writes a pmid_status of
        ``abstract_only`` for audit. Returns (False, reason, None) when even
        the abstract is unavailable.
        """
        from utils.pubmed_utils import batch_fetch_metadata, fetch_paper_abstract

        email = os.environ.get("NCBI_EMAIL") or "gvf@example.com"
        try:
            meta_dict = batch_fetch_metadata([pmid], email=email)
        except Exception as e:
            logger.warning(
                f"Abstract-only fallback metadata fetch failed for {pmid}: {e}"
            )
            meta_dict = {}
        metadata = meta_dict.get(pmid, {}) or {}

        try:
            abstract = fetch_paper_abstract(pmid, email=email)
        except Exception as e:
            logger.warning(
                f"Abstract-only fallback abstract fetch failed for {pmid}: {e}"
            )
            abstract = None

        if not abstract and not metadata:
            return (
                False,
                "Abstract-only fallback: no abstract or metadata available",
                None,
            )

        title = (metadata.get("Title") or "").strip()
        journal = (
            metadata.get("FullJournalName") or metadata.get("Source") or ""
        ).strip()
        pubdate = (metadata.get("PubDate") or "").strip()

        authors: List[str] = []
        for author in metadata.get("AuthorList") or []:
            if isinstance(author, dict):
                name = author.get("Name") or " ".join(
                    p for p in [author.get("ForeName"), author.get("LastName")] if p
                )
                if name:
                    authors.append(str(name))
            elif author:
                authors.append(str(author))

        mesh_terms: List[str] = []
        for m in metadata.get("MeshHeadingList") or []:
            if isinstance(m, dict):
                desc = m.get("DescriptorName") or m.get("Descriptor")
                if desc:
                    mesh_terms.append(str(desc))
            elif m:
                mesh_terms.append(str(m))

        parts = [
            "# ABSTRACT-ONLY FALLBACK",
            "",
            f"> **WARNING:** Full text could not be retrieved for PMID {pmid}.",
            f"> Reason: {reason}",
            ">",
            "> This document contains only the PubMed abstract and metadata.",
            "> Variants extracted here are lower-confidence than from full-text",
            "> papers. Treat as preliminary and seek the full text when possible.",
            "",
        ]
        if title:
            parts += ["## Title", title, ""]
        if authors:
            parts += ["## Authors", ", ".join(authors), ""]
        if journal or pubdate:
            citation = " ".join(
                p for p in [journal, f"({pubdate})" if pubdate else ""] if p
            ).strip()
            parts += ["## Citation", citation, ""]
        if mesh_terms:
            parts += ["## MeSH Terms", "; ".join(mesh_terms), ""]
        if abstract:
            parts += ["## Abstract", abstract.strip(), ""]

        content = "\n".join(parts)

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(content)

        self._log_abstract_only(pmid, reason, str(output_file))

        self._write_pmid_status(
            pmid,
            "abstract_only",
            {
                "download_timestamp": datetime.datetime.now().isoformat(),
                "failure_reason": reason,
                "source": "pubmed_abstract",
            },
        )

        print(
            f"  ⚠️  Abstract-only fallback: wrote {output_file.name} "
            f"(no full text — {reason})"
        )
        logger.info(f"PMID {pmid}: abstract-only fallback used (reason: {reason})")
        return True, str(output_file), content

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
            print("  - No PMCID found, checking for free full text via publisher...")
            success, result, content = self._download_with_tier35_fallback(pmid, doi)
            if success:
                return success, result, content
            return self._build_abstract_only_fallback(
                pmid,
                reason=f"No PMCID and free-text/Tier 3.5 fallback failed: {result}",
            )

        print(f"  ✓ PMCID: {pmcid}")

        # Get full-text XML from PMC
        xml_content = self.pmc_api.get_fulltext_xml(pmcid)

        if not xml_content:
            print("  - Full-text XML not available from PMC, trying publisher APIs...")
            # Don't give up — try publisher APIs via DOI before marking paywalled
            if doi:
                success, result, content = self._download_with_tier35_fallback(
                    pmid, doi
                )
                if success:
                    return success, result, content
                return self._build_abstract_only_fallback(
                    pmid,
                    reason=(
                        f"PMC XML unavailable and publisher/Tier 3.5 fallback "
                        f"failed: {result}"
                    ),
                )
            else:
                print("  ❌ No DOI available for publisher API fallback")
                pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                self._log_paywalled(
                    pmid,
                    "Full-text not available from PMC and no DOI",
                    pmc_url,
                    classification="PAYWALLED",
                )
                return self._build_abstract_only_fallback(
                    pmid,
                    reason="Full-text not available from PMC and no DOI",
                )

        print("  ✓ Full-text XML retrieved from PMC")

        # Per-PMID artifact log built up across the rest of this method.
        artifacts = ArtifactsLog(
            pmid=pmid,
            output_dir=self.output_dir,
            pmcid=pmcid,
            doi=doi,
            gene_symbol=self.gene_symbol,
        )

        # Extract figures from the PMC article page (also yields caption metadata).
        _, html_captions = self._extract_pmc_figures(pmcid, pmid, artifacts=artifacts)

        # Convert main text to markdown
        main_markdown = self.converter.xml_to_markdown(xml_content)
        main_text_source = "pmc_xml"

        # Pull caption-like content (figures, table-wraps, supplementary-material
        # descriptions) from the JATS XML — these are silently dropped by
        # xml_to_markdown which only walks <sec>/<p>.
        try:
            xml_captions = extract_captions_from_jats_xml(xml_content)
        except Exception as e:
            logger.warning(f"JATS caption extraction failed for {pmid}: {e}")
            xml_captions = CaptionExtractionResult()

        captions = merge_caption_results(xml_captions, html_captions)

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
                    main_text_source = "pmc_html"
                    print(
                        "  ✓ Recovered full text from PMC HTML (XML body unavailable)"
                    )
                # Even if XML markdown was longer, the HTML may carry
                # additional <figure> blocks that the JATS extractor missed.
                try:
                    fallback_captions = extract_captions_from_html(html_response.text)
                    captions = merge_caption_results(captions, fallback_captions)
                except Exception as exc:
                    logger.warning(
                        f"PMC HTML caption extraction failed for {pmid}: {exc}"
                    )
            except Exception as e:
                print(f"  - PMC HTML fallback failed: {e}")

        (
            supplement_markdown,
            downloaded_count,
            supp_results,
        ) = self._process_supplements(pmid, pmcid, doi)
        self._record_supplement_results(artifacts, supp_results)

        # Persist a sidecar of the full caption set next to the figures dir.
        try:
            captions_path = self.output_dir / f"{pmid}_figures" / "captions_index.json"
            if captions_path.parent.exists():
                save_captions_json(captions, captions_path)
        except Exception as e:
            logger.warning(f"Failed to write captions_index.json for {pmid}: {e}")

        # Create unified markdown file (WITHOUT pedigree extraction).
        # Captions block lives between main text and supplements so it is
        # contiguous with the article body for the scanner / LLM.
        unified_content = self._build_unified_content(
            main_markdown=main_markdown,
            captions=captions,
            supplement_markdown=supplement_markdown,
        )

        artifacts.record_main_text(
            source=main_text_source,
            chars=len(main_markdown or ""),
            figure_captions=len(captions.figures),
            table_captions=len(captions.tables),
            supplement_descriptions=len(captions.supplements),
        )

        # Validate content quality before writing (catches binary/garbage content)
        is_valid, validation_reason = validate_content_quality(unified_content)
        if not is_valid:
            print(f"  ❌ Content validation failed: {validation_reason}")
            self._log_paywalled(
                pmid,
                f"Content validation failed: {validation_reason}",
                f"PMCID: {pmcid}",
                classification="CONTENT_INVALID",
            )
            artifacts.add_note(f"content_validation_failed: {validation_reason}")
            artifacts.save()
            return self._build_abstract_only_fallback(
                pmid,
                reason=f"Content validation failed: {validation_reason}",
            )

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(unified_content)
        artifacts.save()

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
        append_success_entry(self.success_log, pmid, pmcid, downloaded_count)

        self._write_pmid_status(
            pmid,
            "extracted",
            {
                "download_timestamp": datetime.datetime.now().isoformat(),
                "variant_count": 0,  # Placeholder; extraction updates this later
                "source": "pmc",
            },
        )

        return True, str(output_file), unified_content

    def _download_with_tier35_fallback(
        self, pmid: str, doi: str
    ) -> Tuple[bool, str, Optional[str]]:
        """Run free-text download; if it fails and Tier 3.5 is enabled, retry via browser HTML.

        Tier 3.5 (browser-based HTML) is a strict additive fallback. The
        existing free-text path runs unchanged; only on failure do we attempt
        the browser-driven path. Tier 3.5 successes/failures are logged to
        a separate ``browser_html_log.csv`` so existing paywall metrics are
        unaffected.
        """
        success, result, content = self._download_free_text_pmid(pmid, doi)
        if success:
            return success, result, content

        if self.browser_html is None or not self.browser_html.is_enabled():
            return success, result, content
        if not doi:
            return success, result, content

        try:
            from .browser_html import get_pub_date_from_pmid

            pub_date = get_pub_date_from_pmid(pmid, self.pmc_api)
        except Exception:
            pub_date = None

        print("  → Tier 3.5: trying browser-based HTML fallback...")
        fr = self.browser_html.fetch(pmid=pmid, doi=doi, pub_date=pub_date)
        if fr is None or not fr.is_usable():
            return success, result, content

        return self._finalize_browser_html(pmid, doi, fr)

    def _finalize_browser_html(
        self, pmid: str, doi: str, fr
    ) -> Tuple[bool, str, Optional[str]]:
        """Persist a successful Tier 3.5 fetch using existing free-text plumbing."""
        from .browser_html import FetchResult  # for type clarity only

        supp_files = list(fr.supp_files or [])
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        extract_figures = self._should_extract_figures
        figures_dir = None
        if extract_figures:
            figures_dir = self.output_dir / f"{pmid}_figures"
            figures_dir.mkdir(exist_ok=True)

        # Per-PMID artifact log for this Tier 3.5 fetch.
        artifacts = ArtifactsLog(
            pmid=pmid,
            output_dir=self.output_dir,
            pmcid=None,
            doi=doi,
            gene_symbol=self.gene_symbol,
        )

        # Pull captions from the rendered HTML when the strategy preserved it.
        captions = CaptionExtractionResult()
        raw_html = getattr(fr, "main_html", None)
        if raw_html:
            try:
                captions = extract_captions_from_html(raw_html)
            except Exception as exc:
                logger.warning(f"Tier 3.5 caption extraction failed for {pmid}: {exc}")
                captions = CaptionExtractionResult()
        captions_md = render_captions_markdown(captions)

        # Record any figure files Tier 3.5 already downloaded (best-effort
        # caption pairing — list order is preserved by the strategies).
        figure_paths = list(getattr(fr, "figure_paths", []) or [])
        ordered_caps = list(captions.figures)
        for i, fp in enumerate(figure_paths):
            cap = ordered_caps[i] if i < len(ordered_caps) else None
            artifacts.record_figure(
                ArtifactsLog.figure_artifact_from_path(
                    Path(fp),
                    source=f"browser_html:{fr.publisher or 'generic'}",
                    caption_label=cap.label if cap else None,
                    caption_title=cap.title if cap else None,
                    caption_text=cap.text if cap else None,
                    figure_id=cap.figure_id if cap else None,
                )
            )
        # Caption-only rows for figures we didn't (or couldn't) save.
        if len(ordered_caps) > len(figure_paths):
            for cap in ordered_caps[len(figure_paths) :]:
                artifacts.record_figure_dict(
                    filename=cap.figure_id or cap.label,
                    path="",
                    size_bytes=0,
                    source_url=cap.image_url,
                    source=f"browser_html:{fr.publisher or 'generic'}",
                    caption_label=cap.label,
                    caption_title=cap.title,
                    caption_text=cap.text,
                    figure_id=cap.figure_id,
                )

        if supp_files:
            supp_result = process_supplement_files(
                supp_files=supp_files,
                supplements_dir=supplements_dir,
                pmid=pmid,
                converter=self.converter,
                download_callback=lambda url,
                file_path,
                pmid,
                filename,
                supp: self.download_supplement(
                    url,
                    file_path,
                    pmid,
                    filename,
                    supp.get("base_url"),
                    supp.get("original_url"),
                ),
                extract_figures=extract_figures,
                figures_dir=figures_dir,
                logger=logger,
                sleep_seconds=0.5,
            )
            supplement_markdown = supp_result.supplement_markdown
            downloaded_count = supp_result.downloaded_count
            self._record_supplement_results(artifacts, supp_result.file_results)
        else:
            supplement_markdown = ""
            downloaded_count = 0

        from .free_text_output_service import FreeTextOutputSource

        source = FreeTextOutputSource(
            success_marker="BROWSER_HTML",
            status_source=f"browser-html-{fr.publisher or 'generic'}",
            source_tag=f"[via browser HTML: {fr.publisher or 'generic'}]",
        )

        artifacts.record_main_text(
            source=source.status_source,
            chars=len(fr.main_markdown or ""),
            figure_captions=len(captions.figures),
            table_captions=len(captions.tables),
            supplement_descriptions=len(captions.supplements),
        )

        output_file, unified_content = write_free_text_output(
            output_dir=self.output_dir,
            success_log=self.success_log,
            pmid=pmid,
            main_markdown=fr.main_markdown,
            supplement_markdown=supplement_markdown,
            downloaded_count=downloaded_count,
            source=source,
            write_pmid_status=self._write_pmid_status,
            log_paywalled=None,  # don't pollute the paywall log on Tier 3.5
            captions_markdown=captions_md,
        )
        if output_file is None:
            artifacts.add_note(f"content_validation_failed: {unified_content}")
            artifacts.save()
            return False, unified_content, None

        artifacts.save()
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
        init_state = initialize_free_text_access(
            pmid=pmid,
            doi=doi,
            output_dir=self.output_dir,
            success_log=self.success_log,
            pmc_api=self.pmc_api,
            unpaywall=self.unpaywall,
            converter=self.converter,
            elsevier_api=self.elsevier_api,
            springer_api=self.springer_api,
            wiley_api=self.wiley_api,
            try_elsevier_api=self._try_elsevier_api,
            try_springer_api=self._try_springer_api,
            try_wiley_api=self._try_wiley_api,
            doi_resolver=self.doi_resolver,
            scraper=self.scraper,
            write_pmid_status=self._write_pmid_status,
            log_paywalled=self._log_paywalled,
            logger=logger,
        )

        if init_state.early_result is not None:
            return init_state.early_result

        is_free = init_state.is_free
        free_url = init_state.free_url

        print("  ✓ Article marked as free full text on PubMed")

        content_result = fetch_main_content_for_free_text(
            pmid=pmid,
            doi=doi,
            free_url=free_url,
            suspicious_free_url_domains=self.SUSPICIOUS_FREE_URL_DOMAINS,
            elsevier_api=self.elsevier_api,
            springer_api=self.springer_api,
            wiley_api=self.wiley_api,
            elsevier_client=self.elsevier_client,
            wiley_client=self.wiley_client,
            session=self.session,
            scraper=self.scraper,
            doi_resolver=self.doi_resolver,
            try_elsevier_api=self._try_elsevier_api,
            try_wiley_api=self._try_wiley_api,
            try_springer_api=self._try_springer_api,
            validate_content_quality=validate_content_quality,
            log_paywalled=self._log_paywalled,
        )
        if content_result.early_result is not None:
            return content_result.early_result

        main_markdown = content_result.main_markdown
        final_url = content_result.final_url
        supp_files = content_result.supp_files
        used_elsevier_api = content_result.used_elsevier_api
        used_wiley_api = content_result.used_wiley_api

        # Per-PMID artifact log for the free-text path.
        artifacts = ArtifactsLog(
            pmid=pmid,
            output_dir=self.output_dir,
            pmcid=None,
            doi=doi,
            gene_symbol=self.gene_symbol,
        )

        # Pull captions from raw HTML if the scraper preserved it (the API
        # branches return markdown directly so we have nothing to parse).
        captions = CaptionExtractionResult()
        if getattr(content_result, "raw_html", None):
            try:
                captions = extract_captions_from_html(content_result.raw_html)
            except Exception as exc:
                logger.warning(f"Free-text caption extraction failed for {pmid}: {exc}")
                captions = CaptionExtractionResult()
        captions_md = render_captions_markdown(captions)

        # Even when no images were downloaded, log caption-only artifact rows
        # so we have an audit trail for what the paper contained.
        if captions.figures:
            self._record_html_figure_artifacts(
                artifacts, captions, source_label="publisher_html"
            )

        supplements_dir = self.output_dir / f"{pmid}_supplements"

        extract_figures = self._should_extract_figures

        # Create figures directory if extracting
        figures_dir = None
        if extract_figures:
            figures_dir = self.output_dir / f"{pmid}_figures"
            figures_dir.mkdir(exist_ok=True)

        print(f"  Found {len(supp_files)} supplemental files")
        if not supp_files:
            self._log_paywalled(
                pmid,
                "SUPPLEMENT_NOT_FOUND: checked publisher/API fallbacks",
                final_url or free_url or f"https://doi.org/{doi}" if doi else "",
                classification="SUPPLEMENT_ONLY",
            )
            supplement_markdown = ""
            downloaded_count = 0
            source = source_from_free_text_flags(
                used_elsevier_api=used_elsevier_api,
                used_wiley_api=used_wiley_api,
            )
            artifacts.record_main_text(
                source=source.status_source,
                chars=len(main_markdown or ""),
                figure_captions=len(captions.figures),
                table_captions=len(captions.tables),
                supplement_descriptions=len(captions.supplements),
            )
            output_file, unified_content = write_free_text_output(
                output_dir=self.output_dir,
                success_log=self.success_log,
                pmid=pmid,
                main_markdown=main_markdown,
                supplement_markdown=supplement_markdown,
                downloaded_count=downloaded_count,
                source=source,
                log_paywalled=self._log_paywalled,
                captions_markdown=captions_md,
            )
            # Handle validation failure
            if output_file is None:
                artifacts.add_note(f"content_validation_failed: {unified_content}")
                artifacts.save()
                return (
                    False,
                    unified_content,
                    None,
                )  # unified_content contains error message
            artifacts.save()
            return True, str(output_file), unified_content

        free_text_supp_result = process_supplement_files(
            supp_files=supp_files,
            supplements_dir=supplements_dir,
            pmid=pmid,
            converter=self.converter,
            download_callback=lambda url,
            file_path,
            pmid,
            filename,
            supp: self.download_supplement(url, file_path, pmid, filename),
            extract_figures=extract_figures,
            figures_dir=figures_dir,
            logger=logger,
            sleep_seconds=0.5,
        )

        if free_text_supp_result.total_figures_extracted > 0:
            print(
                f"  ✓ Extracted {free_text_supp_result.total_figures_extracted} figures from PDF supplements"
            )

        supplement_markdown = free_text_supp_result.supplement_markdown
        downloaded_count = free_text_supp_result.downloaded_count
        self._record_supplement_results(artifacts, free_text_supp_result.file_results)
        source = source_from_free_text_flags(
            used_elsevier_api=used_elsevier_api,
            used_wiley_api=used_wiley_api,
        )
        artifacts.record_main_text(
            source=source.status_source,
            chars=len(main_markdown or ""),
            figure_captions=len(captions.figures),
            table_captions=len(captions.tables),
            supplement_descriptions=len(captions.supplements),
        )
        output_file, unified_content = write_free_text_output(
            output_dir=self.output_dir,
            success_log=self.success_log,
            pmid=pmid,
            main_markdown=main_markdown,
            supplement_markdown=supplement_markdown,
            downloaded_count=downloaded_count,
            source=source,
            log_paywalled=self._log_paywalled,
            captions_markdown=captions_md,
        )
        artifacts.save()

        # Handle validation failure
        if output_file is None:
            return (
                False,
                unified_content,
                None,
            )  # unified_content contains error message
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
        print("(Identifying high-value data zones)")
        print(f"{'=' * 60}")

        successful = 0
        failed = 0

        for idx, pmid in enumerate(pmids, 1):
            print(f"\n[{idx}/{len(pmids)}] Scouting PMID {pmid}...", end="")

            try:
                output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
                if not output_file.exists():
                    print(" ❌ File not found")
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
        print("PEDIGREE EXTRACTION PHASE")
        print(f"  {len(pmids_likely_pedigree)} papers likely have pedigrees")
        print(
            f"  (Skipped: {no_figures_count} no figures, {no_indicators_count} no pedigree keywords)"
        )
        print("  (GPT-4o vision calls)")
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
                print(
                    f"  ✓ Full text retrieved via Wiley API ({len(markdown)} characters)"
                )
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
                print("  - Springer API: Full text unavailable (not OpenAccess)")
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
        # PHASE 1: DOWNLOAD ALL PAPERS (with resume support)
        # ============================================
        # Skip PMIDs that already have usable FULL_CONTEXT.md (enables resume
        # after crash), but retry abstract-only fallbacks and publisher shells.
        requested_pmids = {str(p) for p in pmids}
        already_downloaded = set()
        retryable_existing = set()
        for f in self.output_dir.glob("*_FULL_CONTEXT.md"):
            existing_pmid = f.name.replace("_FULL_CONTEXT.md", "")
            if existing_pmid not in requested_pmids:
                continue
            if full_context_needs_retry(f, self.output_dir):
                retryable_existing.add(existing_pmid)
            else:
                already_downloaded.add(existing_pmid)

        remaining_pmids = [p for p in pmids if str(p) not in already_downloaded]

        print(f"{'=' * 60}")
        print(
            f"DOWNLOAD PHASE: {len(pmids)} total, {len(already_downloaded)} already done, {len(remaining_pmids)} to download"
        )
        if retryable_existing:
            print(
                f"Retrying {len(retryable_existing)} existing abstract-only/thin FULL_CONTEXT artifacts"
            )
        print(f"{'=' * 60}")

        downloaded_pmids = list(already_downloaded)
        download_successful = len(already_downloaded)
        download_failed = 0

        for idx, pmid in enumerate(remaining_pmids, 1):
            print(f"[{idx}/{len(remaining_pmids)}]", end=" ")

            try:
                success, result, content = self.download_pmid(pmid)
            except Exception as e:
                logger.error(f"Unexpected error downloading PMID {pmid}: {e}")
                success = False
                result = str(e)
                content = None

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

            if idx < len(remaining_pmids):
                time.sleep(delay)

        # Save manifest after download phase
        manifest.save(manifest_path)
        print(f"\n  📋 Manifest saved: {manifest_path}")

        print(f"\n{'=' * 60}")
        print("Download phase complete!")
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
            print("\n  ℹ️  Skipping post-processing (run_scout=False)")
            print(
                "  To analyze downloaded papers, call batch_post_process() separately."
            )

        # Close Tier 3.5 browser pool if it was opened during this run.
        if self.browser_html is not None:
            try:
                self.browser_html.close()
            except Exception as e:
                logger.debug(f"BrowserHTMLFetcher close failed: {e}")

        # ============================================
        # FINAL SUMMARY
        # ============================================
        print(f"\n{'=' * 60}")
        print("Harvest complete!")
        print(f"  ✅ Successful: {download_successful}")
        print(f"  ❌ Failed: {download_failed}")
        print(f"  📊 Post-processed: {post_successful}/{len(downloaded_pmids)}")
        print(f"  Output directory: {self.output_dir.absolute()}")
        print(f"  Success log: {self.success_log}")
        print(f"  Paywalled log: {self.paywalled_log}")
        if self.browser_html is not None and self.browser_html.attempts_made > 0:
            print(
                f"  Tier 3.5 attempts: {self.browser_html.attempts_made}"
                f" (log: {self.browser_html.log_path})"
            )
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
