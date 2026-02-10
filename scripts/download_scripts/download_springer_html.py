#!/usr/bin/env python3
"""
Download Springer papers via HTML scraping with institutional access.

The Springer OpenAccess API only returns content for Open Access papers.
For paywalled papers, we use DOI resolution + HTML scraping, which works
with Vanderbilt institutional access.

Usage:
    python download_springer_html.py [--all | --dois DOI1 DOI2 ...]
"""

import os
import sys
import time
import re
import argparse
from pathlib import Path
from datetime import datetime

# Add GVF to path
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher')
os.chdir('/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher')

# Load environment
from dotenv import load_dotenv
load_dotenv()

import requests
from harvesting.supplement_scraper import SupplementScraper

# Configuration
OUTPUT_DIR = Path("/mnt/temp2/kronckbm/gvf_output/papers")
LOG_FILE = Path("/mnt/temp2/kronckbm/gvf_output/springer_html_download_log.txt")
RATE_LIMIT = 0.5  # seconds between requests

# Known Springer DOIs that failed with OpenAccess API
SPRINGER_DOIS = [
    "10.1007/s00109-003-0504-1",
    "10.1007/s00115-006-2118-7",
    "10.1007/s00246-008-9377-y",
    "10.1007/s00246-009-9417-2",
    "10.1007/s00246-018-1951-3",
    "10.1007/s00380-015-0693-x",
    "10.1007/s00392-003-0961-0",
    "10.1007/s00414-011-0572-7",
    "10.1007/s00414-013-0853-4",
    "10.1007/s00439-012-1156-4",
    "10.1007/s11596-011-0670-2",
]


def log_message(msg):
    """Log message to both stdout and log file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, 'a') as f:
        f.write(line + "\n")


def download_via_html(doi, session, scraper):
    """
    Download full text via DOI resolution and HTML scraping.
    
    This works with institutional access (Vanderbilt proxy).
    
    Returns:
        tuple: (success: bool, content_or_error: str, char_count: int)
    """
    time.sleep(RATE_LIMIT)
    
    url = f"https://doi.org/{doi}"
    
    try:
        response = session.get(url, allow_redirects=True, timeout=30)
        response.raise_for_status()
        
        final_url = response.url
        
        # Check if we got SpringerLink
        if 'springer.com' not in final_url and 'nature.com' not in final_url:
            return False, f"Redirected to non-Springer URL: {final_url}", 0
        
        # Extract full text
        markdown, title = scraper.extract_fulltext(response.text, final_url)
        
        if markdown and len(markdown) > 1000:
            return True, markdown, len(markdown)
        else:
            # Check if we're behind a paywall (no content extracted)
            if 'Access provided by' in response.text:
                # We have access but extraction failed
                return False, "Extraction failed despite having access", 0
            else:
                return False, "No access or content too short", 0
                
    except requests.exceptions.HTTPError as e:
        return False, f"HTTP error: {e.response.status_code}", 0
    except requests.exceptions.RequestException as e:
        return False, f"Request error: {str(e)}", 0
    except Exception as e:
        return False, f"Error: {str(e)}", 0


def download_supplements(doi, session, scraper):
    """
    Download supplementary files for a Springer paper.
    
    Returns:
        list: List of supplement info dicts
    """
    time.sleep(RATE_LIMIT)
    
    url = f"https://doi.org/{doi}"
    
    try:
        response = session.get(url, allow_redirects=True, timeout=30)
        response.raise_for_status()
        
        # Get supplements
        supplements = scraper.scrape_springer_supplements(response.text, response.url)
        return supplements
        
    except Exception as e:
        return []


def main():
    parser = argparse.ArgumentParser(description="Download Springer papers via HTML scraping")
    parser.add_argument('--all', action='store_true', help='Download all known failed Springer papers')
    parser.add_argument('--dois', nargs='+', help='Specific DOIs to download')
    parser.add_argument('--dry-run', action='store_true', help='Check access without saving')
    args = parser.parse_args()
    
    if args.dois:
        dois = args.dois
    elif args.all:
        dois = SPRINGER_DOIS
    else:
        # Default: download all
        dois = SPRINGER_DOIS
    
    log_message("="*60)
    log_message(f"Springer HTML Download - {len(dois)} papers")
    log_message("="*60)
    
    # Create session with browser-like headers
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9",
        "Accept-Encoding": "gzip, deflate, br",
        "Referer": "https://pubmed.ncbi.nlm.nih.gov/",
        "DNT": "1",
        "Connection": "keep-alive",
    })
    
    scraper = SupplementScraper()
    
    success_count = 0
    failed_count = 0
    results = []
    
    for i, doi in enumerate(dois, 1):
        log_message(f"\n[{i}/{len(dois)}] Processing DOI: {doi}")
        
        success, content_or_error, char_count = download_via_html(doi, session, scraper)
        
        if success:
            if args.dry_run:
                log_message(f"  ✓ Would download: {char_count} chars")
                results.append({'doi': doi, 'status': 'would_succeed', 'chars': char_count})
            else:
                # Save the markdown content
                safe_doi = doi.replace('/', '_').replace(':', '_')
                
                # Delete old empty JATS file if exists
                old_jats = OUTPUT_DIR / f"DOI{safe_doi}__springer_jats.xml"
                if old_jats.exists():
                    old_jats.unlink()
                    log_message(f"  Removed old empty JATS file")
                
                # Save markdown
                md_file = OUTPUT_DIR / f"DOI{safe_doi}__springer_html.md"
                with open(md_file, 'w', encoding='utf-8') as f:
                    f.write(f"# Source: SpringerLink HTML via institutional access\n")
                    f.write(f"# DOI: {doi}\n")
                    f.write(f"# Downloaded: {datetime.now().isoformat()}\n\n")
                    f.write(content_or_error)
                
                log_message(f"  ✓ Downloaded: {md_file.name} ({char_count} chars)")
                results.append({'doi': doi, 'status': 'success', 'chars': char_count, 'file': str(md_file)})
                
                # Also get supplements
                supplements = download_supplements(doi, session, scraper)
                if supplements:
                    log_message(f"  ✓ Found {len(supplements)} supplements")
                    # Save supplement info
                    supp_file = OUTPUT_DIR / f"DOI{safe_doi}__springer_supplements.txt"
                    with open(supp_file, 'w') as f:
                        for supp in supplements:
                            f.write(f"{supp.get('name', 'unknown')}: {supp.get('url', 'no-url')}\n")
                
            success_count += 1
        else:
            log_message(f"  ✗ Failed: {content_or_error}")
            results.append({'doi': doi, 'status': 'failed', 'error': content_or_error})
            failed_count += 1
    
    # Summary
    log_message("\n" + "="*60)
    log_message("SUMMARY")
    log_message("="*60)
    log_message(f"Total processed: {len(dois)}")
    log_message(f"Success: {success_count}")
    log_message(f"Failed: {failed_count}")
    
    if failed_count > 0:
        log_message("\nFailed DOIs:")
        for r in results:
            if r['status'] == 'failed':
                log_message(f"  - {r['doi']}: {r['error']}")
    
    return success_count > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
