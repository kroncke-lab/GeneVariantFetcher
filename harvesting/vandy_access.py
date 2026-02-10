"""
Vanderbilt Institutional Access Manager

Centralized configuration and management for publisher API access
and proxy authentication for paper acquisition pipeline.

Handles:
1. Direct API access (Elsevier, Wiley, Springer)
2. Institutional proxy fallback
3. Access status tracking
4. Batch processing orchestration

Usage:
    from harvesting.vandy_access import InstitutionalAccessManager
    manager = InstitutionalAccessManager()
    results = manager.download_papers(pmids, trackers=True)
"""

import os
import json
import time
import hashlib
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from datetime import datetime
from dataclasses import dataclass
import logging

from harvesting.elsevier_api import ElsevierAPIClient
from harvesting.wiley_api import WileyAPIClient
from harvesting.springer_api import SpringerAPIClient

# Configure logging for this module
logger = logging.getLogger(__name__)

@dataclass
class DownloadResult:
    pmid: str
    provider: str
    status: str
    method: str
    size: int
    error: Optional[str] = None
    filename: Optional[str] = None
    processing_time: float = 0.0

class InstitutionalAccessManager:
    """Main orchestrator for institutional paper access."""
    
    def __init__(self, output_dir: Optional[str] = None):
        """Initialize with institutional access configuration."""
        self.output_dir = Path(output_dir or "/mnt/temp2/kronckbm/gvf_output/institutional_downloads")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize API clients
        self.elsevier = ElsevierAPIClient(api_key=os.getenv("ELSEVIER_API_KEY"))
        self.wiley = WileyAPIClient(api_key=os.getenv("WILEY_API_KEY"))
        self.springer = SpringerAPIClient(api_key=os.getenv("SPRINGER_API_KEY"))
        
        # Vanderbilt proxy configuration
        self.proxy_prefix = os.getenv("VANDERBILT_PROXY_PREFIX", 
                                     "http://proxy.library.vanderbilt.edu/login?url=")
        self.proxy_enabled = os.getenv("VANDERBILT_PROXY_ENABLED", "true").lower() == "true"
        
        # Rate limiting
        self.rate_limits = {
            'elsevier': 0.2,    # 5 req/sec max
            'wiley': 0.5,      # 2 req/sec max  
            'springer': 0.5,   # 2 req/sec max
        }
    
    def get_access_summary(self) -> Dict[str, str]:
        """Get current access status for all providers."""
        return {
            'elsevier': 'VALID' if self.elsevier.is_available else 'NEED_KEY',
            'wiley': 'VALID' if self.wiley.is_available else 'NEED_KEY',
            'springer': 'VALID' if self.springer.is_available else 'NEED_KEY',
            'proxy_enabled': 'ENABLED' if self.proxy_enabled else 'DISABLED'
        }
    
    def validate_access(self) -> Dict[str, Dict]:
        """Validate connectivity to all institutional resources."""
        results = {}
        
        # Test Elsevier API
        try:
            if self.elsevier.is_available:
                # Test with a known DOI 
                xml_content, error = self.elsevier.get_fulltext_by_doi(
                    "10.1016/j.humet.2022.102093"  # Sample Heart Rhythm DOI
                )
                results['elsevier'] = {
                    'status': 'VALID' if xml_content else 'API_LIMIT_REACHED',
                    'key_valid': True,
                    'error': error
                }
            else:
                results['elsevier'] = {'status': 'NO_KEY', 'key_valid': False}
        except Exception as e:
            results['elsevier'] = {'status': 'ERROR', 'error': str(e)}
        
        # Test Wiley API
        try:
            if self.wiley.is_available:
                content, error = self.wiley.get_fulltext_by_doi(
                    "10.1002/humu.25678"  # Sample AJHG DOI
                )
                results['wiley'] = {
                    'status': 'VALID' if content else 'API_LIMIT_REACHED',
                    'key_valid': True,
                    'error': error
                }
            else:
                results['wiley'] = {'status': 'NO_KEY', 'key_valid': False}
        except Exception as e:
            results['wiley'] = {'status': 'ERROR', 'error': str(e)}
        
        # Test Springer API
        try:
            if self.springer.is_available:
                content, error = self.springer.get_fulltext_by_doi(
                    "10.1007/s00439-022-02409-7"  # Sample J Med Genet DOI
                )
                results['springer'] = {
                    'status': 'VALID' if content else 'API_LIMIT_REACHED',
                    'key_valid': True,
                    'error': error
                }
            else:
                results['springer'] = {'status': 'NO_KEY', 'key_valid': False}
        except Exception as e:
            results['springer'] = {'status': 'ERROR', 'error': str(e)}
        
        return results
    
    def download_papers(self, pmids: List[str], 
                       test_count: int = 10,
                       enable_tracker: bool = True) -> List[DownloadResult]:
        """Download papers with institutional access."""
        
        # Create working directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        work_dir = self.output_dir / f"{timestamp}_institutional"
        work_dir.mkdir(exist_ok=True)
        
        # Initialize tracker
        if enable_tracker:
            tracker_path = work_dir / "download_tracker.md"
            self._initialize_tracker(tracker_path, pmids)
        
        # Limit to test count if specified
        target_pmids = pmids[:test_count] if test_count else pmids
        
        results = []
        for i, pmid in enumerate(target_pmids, 1):
            print("[%d/%d] Processing PMID: %s" % (i, len(target_pmids), pmid))
            
            result = self._attempt_download(
                pmid, work_dir, tracker_path if enable_tracker else None
            )
            results.append(result)
            
            # Brief pause between papers
            time.sleep(1)
        
        # Generate final report
        if enable_tracker:
            self._generate_summary_report(work_dir, results)
        
        return results
    
    def _attempt_download(self, pmid: str, work_dir: Path, 
                         tracker_path: Optional[Path]) -> DownloadResult:
        """Attempt download through multiple channels."""
        
        start_time = time.time()
        
        # Get DOI and publisher info (would typically be resolved from PubMed)
        doi, publisher = self._resolve_paper_info(pmid)
        
        if not doi:
            result = DownloadResult(
                pmid=pmid, provider="unknown", 
                status="FAIL", method="resolution_failed",
                size=0, error="Could not resolve DOI from PMID"
            )
            self._log_result(result, tracker_path)
            return result
        
        # Try providers in order of preference
        providers = [
            ('elsevier', self._download_elsevier, self.elsevier),
            ('wiley', self._download_wiley, self.wiley), 
            ('springer', self._download_springer, self.springer)
        ]
        
        for provider_name, handler, client in providers:
            if not client.is_available:
                continue
                
            try:
                print(f"  Trying {provider_name}...")
                content, filename = handler(doi, work_dir)
                if content:
                    processing_time = time.time() - start_time
                    result = DownloadResult(
                        pmid=pmid, provider=provider_name,
                        status="SUCCESS", method="API",
                        size=len(content), filename=filename,
                        processing_time=processing_time
                    )
                    self._log_result(result, tracker_path)
                    return result
                    
            except Exception as e:
                print(f"  {provider_name} failed: {e}")
                continue
        
        # Try proxy fallback
        if self.proxy_enabled:
            content, filename = self._download_via_proxy(pmid, work_dir)
            if content:
                processing_time = time.time() - start_time
                result = DownloadResult(
                    pmid=pmid, provider="proxy",
                    status="SUCCESS", method="PROXY_FALLBACK",
                    size=len(content), filename=filename,
                    processing_time=processing_time
                )
                self._log_result(result, tracker_path)
                return result
        
        # All methods failed
        processing_time = time.time() - start_time
        result = DownloadResult(
            pmid=pmid, provider="multiple",
            status="FAIL", method="all_exhausted",
            size=0, error="All access methods failed (API + proxy)",
            processing_time=processing_time
        )
        self._log_result(result, tracker_path)
        return result
    
    def _resolve_paper_info(self, pmid: str) -> Tuple[Optional[str], Optional[str]]:
        """Resolve PMID to DOI and determine publisher."""
        # This is a placeholder - in production you'd query PubMed
        # For now, we'll simulate resolution based on PMID pattern
        
        # Simulated mapping for testing
        pmid_to_publisher = {
            "37254019": ("10.1016/j.heart.2023.123456", "elsevier"),
            "36912273": ("10.1002/humu.25678", "wiley"),
            "36880547": ("10.1007/s00439-022-02409-7", "springer"),
            # Add more mappings for testing
        }
        
        return pmid_to_publisher.get(pmid, (None, None))
    
    def _download_elsevier(self, doi: str, work_dir: Path) -> Tuple[Optional[bytes], Optional[str]]:
        """Download from Elsevier API."""
        xml_content, error = self.elsevier.get_fulltext_by_doi(doi)
        if xml_content:
            filename = f"{doi.replace('/', '_')}_elsevier.xml"
            filepath = work_dir / filename
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(xml_content)
            return xml_content.encode(), filename
        return None, None
    
    def _download_wiley(self, doi: str, work_dir: Path) -> Tuple[Optional[bytes], Optional[str]]:
        """Download from Wiley API."""
        content, error = self.wiley.get_fulltext_by_doi(doi)
        if content:
            filename = f"{doi.replace('/', '_')}_wiley.html"
            filepath = work_dir / filename
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            return content.encode(), filename
        return None, None
    
    def _download_springer(self, doi: str, work_dir: Path) -> Tuple[Optional[bytes], Optional[str]]:
        """Download from Springer API."""
        content, error = self.springer.get_fulltext_by_doi(doi)
        if content:
            filename = f"{doi.replace('/', '_')}_springer.xml"
            filepath = work_dir / filename
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            return content.encode(), filename
        return None, None
    
    def _download_via_proxy(self, pmid: str, work_dir: Path) -> Tuple[Optional[bytes], Optional[str]]:
        """Download via institutional proxy."""
        # This would involve browser automation for authenticated sessions
        # Placeholder for proxy-based downloads
        return None, None
    
    def _initialize_tracker(self, tracker_path: Path, pmids: List[str]):
        """Initialize download tracker."""
        with open(tracker_path, 'w') as f:
            f.write(f"# Institutional Paper Acquisition Tracker\n")
            f.write(f"**Started:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Target Count:** {len(pmids)} papers\n\n")
            
            f.write("## API Status\n")
            for provider, status in self.get_access_summary().items():
                f.write(f"- **{provider.title()}**: {status}\n")
            f.write("\n")
            
            f.write("## Progress Log\n")
            f.write("| Time | PMID | Provider | Status | Size | Method | Error\n")
            f.write("|------|------|----------|--------|------|--------|------\n")
    
    def _log_result(self, result: DownloadResult, tracker_path: Optional[Path]):
        """Log result to tracker."""
        if not tracker_path:
            return
            
        with open(tracker_path, 'a') as f:
            f.write(f"| {datetime.now().strftime('%H:%M:%S')}")
            f.write(f"| {result.pmid}")
            f.write(f"| {result.provider}")
            f.write(f"| {result.status}")
            f.write(f"| {result.size}")
            f.write(f"| {result.method}")
            f.write(f"| {'None' if not result.error else result.error[:50]}\n")
    
    def _generate_summary_report(self, work_dir: Path, results: List[DownloadResult]):
        """Generate final summary report."""
        success_count = sum(1 for r in results if r.status == "SUCCESS")
        failure_count = len(results) - success_count
        
        # Detailed analysis
        providers = {}
        methods = {}
        for r in results:
            providers[r.provider] = providers.get(r.provider, 0) + (1 if r.status == "SUCCESS" else 0)
            methods[r.method] = methods.get(r.method, 0) + 1
        
        report_path = work_dir / "final_report.md"
        with open(report_path, 'w') as f:
            f.write(f"# Institutional Access Summary\n\n")
            f.write(f"**Completed:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write(f"## Results Summary\n")
            f.write(f"- **Total Papers:** {len(results)}\n")
            f.write(f"- **Successful:** {success_count} ({success_count/len(results)*100:.1f}%)\n")
            f.write(f"- **Failed:** {failure_count} ({failure_count/len(results)*100:.1f}%)\n\n")
            
            f.write("## Provider Analysis\n")
            for provider, count in providers.items():
                f.write(f"- **{provider.title()}**: {count} successes\n")
            
            f.write("\n## Method Distribution\n")
            for method, count in methods.items():
                f.write(f"- **{method}**: {count} attempts\n")
            
            f.write("\n## Recommendations\n")
            if success_count >= len(results) * 0.9:
                f.write("✅ **Ready for full 224-paper queue**\n")
            elif success_count >= len(results) * 0.7:
                f.write("⚠️ **Consider manual intervention for failed papers**\n")
            else:
                f.write("❌ **Address API issues before scaling**\n")

if __name__ == "__main__":
    # For testing
    manager = InstitutionalAccessManager()
    
    # Validate access
    print("Validating access...")
    access_status = manager.validate_access()
    for provider, status in access_status.items():
        print(f"{provider}: {status}")
    
    # Test with sample data
    sample_pmids = ["37254019", "36912273", "36880547", "36712345", "36654321"]
    print(f"\nTesting with {len(sample_pmids)} sample papers...")
    results = manager.download_papers(sample_pmids, test_count=5)
    
    # Display summary
    successful = sum(1 for r in results if r.status == "SUCCESS")
    print(f"\nTest Results: {successful}/{len(results)} successful")
    print(f"Output saved to: {manager.output_dir}")