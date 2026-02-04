#!/usr/bin/env python3
"""
Bulk Download Script for GVF

Downloads all papers for a gene using available APIs:
1. PMC (free, no key needed)
2. Elsevier API (requires ELSEVIER_API_KEY)
3. Wiley API (requires WILEY_TDM_TOKEN)

Tracks progress and failures for later browser fallback.
"""

import os
import sys
import json
import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Tuple

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # dotenv not available, assume env vars already set

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from harvesting.orchestrator import PMCHarvester
from harvesting.elsevier_api import ElsevierAPIClient
from harvesting.wiley_api import WileyAPIClient
from harvesting.pmc_api import PMCAPIClient

# Environment already loaded at top of file

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_pmids(pmid_file: Path) -> List[str]:
    """Load PMIDs from a text file (one per line)."""
    pmids = []
    with open(pmid_file, 'r') as f:
        for line in f:
            pmid = line.strip()
            if pmid and pmid.isdigit():
                pmids.append(pmid)
    return pmids


def get_already_downloaded(output_dir: Path) -> Set[str]:
    """Get set of PMIDs that already have full-text downloaded."""
    downloaded = set()
    
    # Check for FULL_CONTEXT.md files
    for f in output_dir.glob("*_FULL_CONTEXT.md"):
        # Extract PMID from filename like "12345678_FULL_CONTEXT.md"
        parts = f.stem.split('_')
        if parts and parts[-2].isdigit():
            downloaded.add(parts[-2])
        elif len(parts) >= 3 and parts[-3].isdigit():
            # Handle GENE_PMID_12345678_FULL_CONTEXT.md
            downloaded.add(parts[-3])
    
    # Also check the success log
    success_log = output_dir / "successful_downloads.csv"
    if success_log.exists():
        import csv
        with open(success_log, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'PMID' in row:
                    downloaded.add(str(row['PMID']))
    
    return downloaded


def get_doi_for_pmid(pmid: str, pmc_client: PMCAPIClient) -> Tuple[str, str]:
    """Get DOI and PMCID for a PMID."""
    try:
        result = pmc_client.pmid_to_ids(pmid)
        return result.get('doi', ''), result.get('pmcid', '')
    except Exception as e:
        logger.debug(f"Could not get IDs for PMID {pmid}: {e}")
        return '', ''


def bulk_download(
    pmid_file: Path,
    output_dir: Path,
    gene_symbol: str = "KCNH2",
    max_concurrent: int = 1,
    delay: float = 0.5,
    resume: bool = True,
) -> Dict:
    """
    Bulk download papers using all available APIs.
    
    Args:
        pmid_file: Path to file with PMIDs (one per line)
        output_dir: Directory to save downloaded content
        gene_symbol: Gene symbol for naming
        max_concurrent: Max concurrent downloads (keep low for rate limits)
        delay: Delay between requests in seconds
        resume: Skip already downloaded papers
        
    Returns:
        Dict with download statistics
    """
    
    # Initialize API clients
    elsevier_key = os.getenv('ELSEVIER_API_KEY')
    wiley_key = os.getenv('WILEY_API_KEY') or os.getenv('WILEY_TDM_TOKEN')
    
    elsevier = ElsevierAPIClient(api_key=elsevier_key) if elsevier_key else None
    wiley = WileyAPIClient(api_key=wiley_key) if wiley_key else None
    pmc = PMCAPIClient()
    
    logger.info("=" * 70)
    logger.info(f"BULK DOWNLOAD FOR {gene_symbol}")
    logger.info("=" * 70)
    logger.info(f"Elsevier API: {'✓ Available' if elsevier and elsevier.is_available else '✗ Not configured'}")
    logger.info(f"Wiley API: {'✓ Available' if wiley and wiley.is_available else '✗ Not configured'}")
    logger.info(f"PMC API: ✓ Always available")
    logger.info("=" * 70)
    
    # Load PMIDs
    all_pmids = load_pmids(pmid_file)
    logger.info(f"Loaded {len(all_pmids)} PMIDs from {pmid_file}")
    
    # Get already downloaded
    if resume:
        already_done = get_already_downloaded(output_dir)
        logger.info(f"Already downloaded: {len(already_done)} papers")
        pmids_to_download = [p for p in all_pmids if p not in already_done]
    else:
        pmids_to_download = all_pmids
        
    logger.info(f"Papers to download: {len(pmids_to_download)}")
    
    if not pmids_to_download:
        logger.info("Nothing to download!")
        return {"total": len(all_pmids), "already_done": len(already_done), "downloaded": 0}
    
    # Create output dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Track results
    stats = {
        "total": len(all_pmids),
        "attempted": 0,
        "pmc_success": 0,
        "elsevier_success": 0,
        "wiley_success": 0,
        "failed": 0,
        "errors": [],
    }
    
    failures = []  # Track failed PMIDs for browser fallback
    
    # Progress tracking
    start_time = time.time()
    
    for i, pmid in enumerate(pmids_to_download):
        stats["attempted"] += 1
        
        # Progress update every 10 papers
        if (i + 1) % 10 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(pmids_to_download) - i - 1) / rate if rate > 0 else 0
            logger.info(f"Progress: {i+1}/{len(pmids_to_download)} "
                       f"({100*(i+1)/len(pmids_to_download):.1f}%) "
                       f"- ETA: {remaining/60:.1f} min")
        
        success = False
        source = None
        
        try:
            # Get DOI and PMCID
            doi, pmcid = get_doi_for_pmid(pmid, pmc)
            
            # Try 1: PMC (if has PMCID)
            if pmcid:
                try:
                    harvester = PMCHarvester(output_dir)
                    harvester.harvest([pmid], delay=delay)
                    
                    # Check if file was created
                    expected_file = output_dir / f"{gene_symbol}_PMID_{pmid}_FULL_CONTEXT.md"
                    alt_file = output_dir / f"{pmid}_FULL_CONTEXT.md"
                    
                    if expected_file.exists() or alt_file.exists():
                        success = True
                        source = "PMC"
                        stats["pmc_success"] += 1
                except Exception as e:
                    logger.debug(f"PMC failed for {pmid}: {e}")
            
            # Try 2: Elsevier API (if DOI matches)
            if not success and doi and elsevier and elsevier.is_available:
                if elsevier.is_elsevier_doi(doi):
                    try:
                        content = elsevier.get_article_fulltext(doi)
                        if content:
                            out_file = output_dir / f"{gene_symbol}_PMID_{pmid}_FULL_CONTEXT.md"
                            out_file.write_text(content, encoding='utf-8')
                            success = True
                            source = "Elsevier"
                            stats["elsevier_success"] += 1
                    except Exception as e:
                        logger.debug(f"Elsevier API failed for {pmid}: {e}")
            
            # Try 3: Wiley API (if DOI matches)
            if not success and doi and wiley and wiley.is_available:
                if wiley.is_wiley_doi(doi):
                    try:
                        content = wiley.get_article_fulltext(doi)
                        if content:
                            out_file = output_dir / f"{gene_symbol}_PMID_{pmid}_FULL_CONTEXT.md"
                            out_file.write_text(content, encoding='utf-8')
                            success = True
                            source = "Wiley"
                            stats["wiley_success"] += 1
                    except Exception as e:
                        logger.debug(f"Wiley API failed for {pmid}: {e}")
            
            if not success:
                stats["failed"] += 1
                failures.append({
                    "pmid": pmid,
                    "doi": doi or "",
                    "pmcid": pmcid or "",
                    "reason": "No API access available"
                })
            else:
                logger.debug(f"✓ Downloaded PMID {pmid} via {source}")
                
        except Exception as e:
            stats["failed"] += 1
            stats["errors"].append(f"{pmid}: {str(e)}")
            failures.append({
                "pmid": pmid,
                "doi": "",
                "pmcid": "",
                "reason": str(e)
            })
        
        # Rate limiting
        time.sleep(delay)
    
    # Save failures for browser fallback
    if failures:
        failures_file = output_dir / "download_failures.json"
        with open(failures_file, 'w') as f:
            json.dump(failures, f, indent=2)
        logger.info(f"Saved {len(failures)} failures to {failures_file}")
    
    # Save stats
    stats_file = output_dir / "bulk_download_stats.json"
    stats["duration_seconds"] = time.time() - start_time
    stats["timestamp"] = datetime.now().isoformat()
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Summary
    logger.info("=" * 70)
    logger.info("DOWNLOAD COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Total PMIDs: {stats['total']}")
    logger.info(f"Attempted: {stats['attempted']}")
    logger.info(f"PMC success: {stats['pmc_success']}")
    logger.info(f"Elsevier API success: {stats['elsevier_success']}")
    logger.info(f"Wiley API success: {stats['wiley_success']}")
    logger.info(f"Failed (need browser): {stats['failed']}")
    logger.info(f"Duration: {stats['duration_seconds']/60:.1f} minutes")
    logger.info("=" * 70)
    
    return stats


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Bulk download papers for GVF")
    parser.add_argument("--pmids", required=True, help="Path to PMID file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--gene", default="KCNH2", help="Gene symbol")
    parser.add_argument("--delay", type=float, default=0.5, help="Delay between requests")
    parser.add_argument("--no-resume", action="store_true", help="Don't skip already downloaded")
    
    args = parser.parse_args()
    
    bulk_download(
        pmid_file=Path(args.pmids),
        output_dir=Path(args.output),
        gene_symbol=args.gene,
        delay=args.delay,
        resume=not args.no_resume,
    )
