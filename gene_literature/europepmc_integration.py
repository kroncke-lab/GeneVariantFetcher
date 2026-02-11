"""
Europe PMC Harvesting Integration

Integrates Europe PMC full-text API with the existing PMCHarvester workflow.
Provides an alternative source for papers when regular PMC doesn't have full-text.
"""

import logging
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .europepmc_handler import EuropePMCClient
from harvesting.format_converters import FormatConverter

logger = logging.getLogger(__name__)


class EuropePMCHarvester:
    """
    Europe PMC-based harvesting orchestrator.

    This class provides an alternative pathway for obtaining full-text articles
    when traditional PubMed Central methods fail. It leverages Europe PMC's
    broader coverage of 6.5+ million open access articles.
    """

    def __init__(self, output_dir: Path = None):
        """
        Initialize Europe PMC harvester.

        Args:
            output_dir: Directory to save downloaded files
        """
        self.output_dir = (
            Path(output_dir) if output_dir else Path("europepmc_downloads")
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.client = EuropePMCClient()
        self.converter = FormatConverter()

    def is_article_available(self, pmid: str) -> Dict:
        """
        Determine if a paper is available through Europe PMC.

        Args:
            pmid: PubMed ID to check

        Returns:
            Dictionary with availability status and metadata
        """
        metadata = self.client.get_paper_metadata(pmid)

        if not metadata:
            return {
                "available": False,
                "pmid": pmid,
                "reason": "No metadata found in Europe PMC",
            }

        # Check for required components
        has_pmcid = bool(metadata.get("pmcid"))
        has_fulltext = metadata.get("fulltext_frequency", "N") == "Y"
        is_oa = metadata.get("is_open_access", False)

        if not has_pmcid:
            return {
                "available": False,
                "pmid": pmid,
                "reason": "No PMCID found - full-text not available",
            }

        if not has_fulltext:
            return {
                "available": False,
                "pmid": pmid,
                "reason": "Only abstract/metadata available (no full-text)",
                "pmcid": metadata.get("pmcid"),
            }

        return {
            "available": True,
            "pmid": pmid,
            "pmcid": metadata.get("pmcid"),
            "is_open_access": is_oa,
            "title": metadata.get("title", ""),
            "links": metadata.get("links", {}),
        }

    def download_full_text(self, pmid: str, variants: List[str] = None) -> Dict:
        """
        Download full-text article and supplementary materials from Europe PMC.

        Args:
            pmid: PubMed ID to download
            variants: List of gene variants being searched (for filtering/validation)

        Returns:
            Dictionary with download results
        """
        logger.info(f"Starting Europe PMC download for PMID: {pmid}")

        # Check availability first
        availability = self.is_article_available(pmid)
        if not availability["available"]:
            logger.warning(f"PMID {pmid} not available: {availability['reason']}")
            return {
                "success": False,
                "pmid": pmid,
                "error": availability["reason"],
                "files": [],
            }

        # Create output subfolder
        paper_dir = self.output_dir / f"PMID_{pmid}"
        paper_dir.mkdir(exist_ok=True)

        # Download complete package
        try:
            download_result = self.client.download_full_paper(pmid, paper_dir)

            if not download_result["success"]:
                error_msg = "; ".join(download_result.get("errors", ["Unknown error"]))
                return {
                    "success": False,
                    "pmid": pmid,
                    "error": error_msg,
                    "files": download_result["files_downloaded"],
                }

            # Create unified markdown if we have XML
            xml_files = [
                f for f in download_result["files_downloaded"] if f.endswith(".xml")
            ]
            markdown_file = None

            if xml_files:
                xml_path = Path(xml_files[0])
                try:
                    with open(xml_path, "r", encoding="utf-8") as f:
                        xml_content = f.read()

                    # Convert XML to markdown
                    markdown_content = self.converter.xml_to_markdown(xml_content)

                    # Create unified markdown file
                    markdown_file = paper_dir / f"{pmid}_europepmc_unified.md"
                    with open(markdown_file, "w", encoding="utf-8") as f:
                        f.write(markdown_content)

                    logger.info(f"Created unified markdown: {markdown_file}")

                except Exception as e:
                    logger.warning(f"Failed to create unified markdown: {e}")

            # Return comprehensive result
            return {
                "success": True,
                "pmid": pmid,
                "pmcid": availability["pmcid"],
                "title": availability["title"],
                "files": download_result["files_downloaded"],
                "markdown_file": str(markdown_file) if markdown_file else None,
                "is_open_access": availability["is_open_access"],
                "citation_count": int(
                    download_result["metadata"].get("citation_count", 0)
                ),
            }

        except Exception as e:
            logger.error(f"Unexpected error during download for PMID {pmid}: {e}")
            return {
                "success": False,
                "pmid": pmid,
                "error": f"Unexpected error: {str(e)}",
            }

    def download_batch(self, pmids: List[str], output_dir: Path = None) -> List[Dict]:
        """
        Download multiple papers from Europe PMC in batch.

        Args:
            pmids: List of PubMed IDs to download
            output_dir: Optional custom output directory

        Returns:
            List of download results for each PMID
        """
        if output_dir:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)

        results = []

        logger.info(f"Starting batch download of {len(pmids)} papers from Europe PMC")

        for i, pmid in enumerate(pmids, 1):
            logger.info(f"Processing {i}/{len(pmids)}: PMID {pmid}")

            try:
                result = self.download_full_text(pmid)
                results.append(result)

                # Log summary of result
                if result["success"]:
                    files_count = len(result["files"])
                    logger.info(f"✓ Successfully downloaded → {files_count} files")
                else:
                    logger.warning(f"✗ Failed → {result['error']}")

            except Exception as e:
                error_result = {
                    "success": False,
                    "pmid": pmid,
                    "error": f"Batch processing failed: {str(e)}",
                }
                results.append(error_result)
                logger.error(f"Critical error processing PMID {pmid}: {e}")

        # Summary report
        successful = [r for r in results if r["success"]]
        failed = [r for r in results if not r["success"]]

        logger.info(f"Batch download complete:")
        logger.info(f"  ✅ Successful: {len(successful)}/{len(pmids)}")
        logger.info(f"  ❌ Failed: {len(failed)}/{len(pmids)}")

        if failed:
            logger.info(f"Failed PMIDs: {[r['pmid'] for r in failed]}")

        return results

    def search_and_download(
        self, gene: str, max_results: int = 50, output_dir: Path = None
    ) -> List[Dict]:
        """
        Search for papers by gene symbol and download matching articles.

        Args:
            gene: Gene symbol to search for
            max_results: Maximum number of papers to return
            output_dir: Optional custom output directory

        Returns:
            List of download results
        """
        if output_dir:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)

        # Search for papers with the gene
        search_results = self.client.search_papers(
            f"{gene}[Title/Abstract]", max_results
        )

        if not search_results:
            logger.warning(f"No papers found for gene: {gene}")
            return []

        # Extract PMIDs from search results
        pmids = [
            str(paper.get("pmid", "")) for paper in search_results if paper.get("pmid")
        ]

        logger.info(f"Found {len(pmids)} papers for gene {gene}")

        # Download the papers in batch
        return self.download_batch(pmids, output_dir)

    def get_unavailable_papers(self, pmids: List[str]) -> List[Dict]:
        """
        Identify which papers from a list are not available in Europe PMC.

        Args:
            pmids: List of PubMed IDs to check

        Returns:
            List of unavailable papers with reasons
        """
        unavailable = []

        for pmid in pmids:
            availability = self.is_article_available(pmid)
            if not availability["available"]:
                unavailable.append({"pmid": pmid, "reason": availability["reason"]})

        return unavailable

    def compare_with_pmc(self, pmid: str) -> Dict:
        """
        Compare availability between PubMed Central and Europe PMC.

        Args:
            pmid: PubMed ID to compare

        Returns:
            Dictionary with comparison results
        """
        from harvesting.pmc_api import PMCAPIClient

        pmc_client = PMCAPIClient()
        europepmc_result = self.is_article_available(pmid)

        # Check PMC availability
        pmcid = pmc_client.pmid_to_pmcid(pmid)
        xml_content = pmc_client.get_fulltext_xml(pmcid) if pmcid else None

        pmc_available = xml_content and len(xml_content) > 1000

        return {
            "pmid": pmid,
            "pmc_available": pmc_available,
            "pmcid": pmcid,
            "europepmc_available": europepmc_result["available"],
            "europepmc_pmcid": europepmc_result.get("pmcid"),
            "europepmc_is_oa": europepmc_result.get("is_open_access"),
            "coverage_improvement": europepmc_result["available"] and not pmc_available,
            "both_available": pmc_available and europepmc_result["available"],
        }

    def format_content_for_extraction(self, pmid: str) -> Optional[str]:
        """
        Prepare content from Europe PMC for variant extraction.

        Args:
            pmid: PubMed ID

        Returns:
            Formatted markdown content ready for extraction, or None if unavailable
        """
        availability = self.is_article_available(pmid)
        if not availability["available"]:
            return None

        # Check if we have the unified markdown file
        paper_dir = self.output_dir / f"PMID_{pmid}"
        markdown_file = paper_dir / f"{pmid}_europepmc_unified.md"

        if markdown_file.exists():
            with open(markdown_file, "r", encoding="utf-8") as f:
                content = f.read()
        else:
            # Reconstruct from XML if needed
            xml_files = list(paper_dir.glob(f"{pmid}_europepmc_full.xml"))
            if xml_files:
                with open(xml_files[0], "r", encoding="utf-8") as f:
                    xml_content = f.read()
                content = self.converter.xml_to_markdown(xml_content)
            else:
                return None

        # Add metadata header
        header = f"""# Europe PMC Extraction: PM${pmid}
**Title**: {availability.get("title", "N/A")}
**PMC ID**: {availability.get("pmcid", "N/A")}
**Open Access**: {"Yes" if availability.get("is_open_access") else "No"}
**Source**: Europe PMC RESTful API

---

"""

        return header + content
