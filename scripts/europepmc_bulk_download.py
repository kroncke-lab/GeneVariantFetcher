#!/usr/bin/env python3
"""
Bulk download from Europe PMC for gold standard KCNH2 PMIDs.

Loads 262 PMIDs from the Excel baseline file and uses EuropePMCClient to:
- Check if each paper exists in Europe PMC
- Download full-text XML if available (via PMCID)
- Download supplementary files if available
- Track success rate

Usage:
    python scripts/europepmc_bulk_download.py
"""

import sys
import os
import json
import time
import logging
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Add user site-packages for xlrd
sys.path.insert(0, '/home/kronckbm/.local/lib/python3.9/site-packages')

import xlrd
from gene_literature.europepmc_handler import EuropePMCClient

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
    ]
)
logger = logging.getLogger(__name__)

# Configuration
EXCEL_FILE = "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
OUTPUT_DIR = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2/europepmc_downloads")
PMID_COLUMN = 24  # Column Y (0-indexed)
RATE_LIMIT_DELAY = 0.5  # seconds between requests to be polite


def load_pmids_from_excel(excel_path: str) -> list:
    """Extract unique PMIDs from the gold standard Excel file."""
    wb = xlrd.open_workbook(excel_path)
    sheet = wb.sheet_by_index(0)

    pmids = set()
    for r in range(1, sheet.nrows):
        val = sheet.cell_value(r, PMID_COLUMN)
        if val:
            if isinstance(val, float):
                val = str(int(val))
            else:
                val = str(val).strip()
            for part in val.replace(';', ',').split(','):
                part = part.strip()
                if part and part.isdigit() and len(part) >= 6:
                    pmids.add(part)

    return sorted(pmids)


def main():
    # Load PMIDs
    logger.info(f"Loading PMIDs from {EXCEL_FILE}")
    pmids = load_pmids_from_excel(EXCEL_FILE)
    logger.info(f"Found {len(pmids)} unique PMIDs")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load progress file to support resuming
    progress_file = OUTPUT_DIR / "download_progress.json"
    if progress_file.exists():
        with open(progress_file, 'r') as f:
            progress = json.load(f)
        logger.info(f"Resuming from previous run: {len(progress['completed'])} already done")
    else:
        progress = {
            "completed": [],
            "results": {}
        }

    # Initialize client
    client = EuropePMCClient(timeout=30)

    # Tracking
    stats = {
        "total": len(pmids),
        "processed": 0,
        "found_in_europepmc": 0,
        "has_pmcid": 0,
        "fulltext_xml_downloaded": 0,
        "supplements_found": 0,
        "supplements_downloaded": 0,
        "metadata_only": 0,
        "not_found": 0,
        "errors": 0,
    }

    for i, pmid in enumerate(pmids, 1):
        if pmid in progress["completed"]:
            # Recount stats from saved results
            r = progress["results"].get(pmid, {})
            stats["processed"] += 1
            if r.get("found"):
                stats["found_in_europepmc"] += 1
            if r.get("pmcid"):
                stats["has_pmcid"] += 1
            if r.get("fulltext_xml"):
                stats["fulltext_xml_downloaded"] += 1
            if r.get("supplements_count", 0) > 0:
                stats["supplements_found"] += 1
                stats["supplements_downloaded"] += r.get("supplements_downloaded", 0)
            if r.get("found") and not r.get("pmcid"):
                stats["metadata_only"] += 1
            if not r.get("found"):
                stats["not_found"] += 1
            continue

        logger.info(f"[{i}/{len(pmids)}] Processing PMID {pmid}...")
        result = {
            "pmid": pmid,
            "found": False,
            "pmcid": None,
            "fulltext_xml": False,
            "supplements_count": 0,
            "supplements_downloaded": 0,
            "files": [],
            "errors": [],
        }

        try:
            # Step 1: Get metadata
            metadata = client.get_paper_metadata(pmid)
            time.sleep(RATE_LIMIT_DELAY)

            if not metadata:
                result["errors"].append("Not found in Europe PMC")
                stats["not_found"] += 1
                logger.warning(f"  PMID {pmid}: Not found in Europe PMC")
            else:
                result["found"] = True
                stats["found_in_europepmc"] += 1

                # Save metadata
                pmid_dir = OUTPUT_DIR / pmid
                pmid_dir.mkdir(exist_ok=True)

                meta_file = pmid_dir / f"{pmid}_metadata.json"
                with open(meta_file, 'w') as f:
                    json.dump(metadata, f, indent=2)
                result["files"].append(str(meta_file))

                pmcid = metadata.get("pmcid")
                result["pmcid"] = pmcid

                if pmcid:
                    stats["has_pmcid"] += 1
                    logger.info(f"  PMCID: {pmcid}")

                    # Step 2: Download full-text XML
                    xml_content = client.get_fulltext_xml(pmcid)
                    time.sleep(RATE_LIMIT_DELAY)

                    if xml_content:
                        xml_file = pmid_dir / f"{pmid}_fulltext.xml"
                        with open(xml_file, 'w', encoding='utf-8') as f:
                            f.write(xml_content)
                        result["fulltext_xml"] = True
                        result["files"].append(str(xml_file))
                        stats["fulltext_xml_downloaded"] += 1
                        logger.info(f"  Full-text XML: Downloaded ({len(xml_content)} chars)")
                    else:
                        logger.info(f"  Full-text XML: Not available")

                    # Step 3: Get supplementary files
                    supplements = client.get_supplementary_files(pmcid)
                    time.sleep(RATE_LIMIT_DELAY)

                    result["supplements_count"] = len(supplements)
                    if supplements:
                        stats["supplements_found"] += 1
                        logger.info(f"  Supplements: {len(supplements)} found")

                        supp_dir = pmid_dir / "supplements"
                        supp_dir.mkdir(exist_ok=True)

                        for idx, supp in enumerate(supplements, 1):
                            supp_url = supp.get('downloadUrl', '') or supp.get('url', '')
                            supp_name = supp.get('filename', '') or supp.get('name', f"supplement_{idx}")

                            if not supp_url:
                                continue

                            try:
                                resp = client.session.get(supp_url, timeout=30)
                                resp.raise_for_status()
                                supp_file = supp_dir / supp_name
                                with open(supp_file, 'wb') as f:
                                    f.write(resp.content)
                                result["supplements_downloaded"] += 1
                                stats["supplements_downloaded"] += 1
                                result["files"].append(str(supp_file))
                                time.sleep(RATE_LIMIT_DELAY)
                            except Exception as e:
                                result["errors"].append(f"Supplement {supp_name}: {str(e)}")
                    else:
                        logger.info(f"  Supplements: None found")
                else:
                    stats["metadata_only"] += 1
                    logger.info(f"  No PMCID - metadata only (not in PMC)")

        except Exception as e:
            stats["errors"] += 1
            result["errors"].append(str(e))
            logger.error(f"  Error processing PMID {pmid}: {e}")

        stats["processed"] += 1
        progress["completed"].append(pmid)
        progress["results"][pmid] = result

        # Save progress every 10 papers
        if stats["processed"] % 10 == 0:
            with open(progress_file, 'w') as f:
                json.dump(progress, f, indent=2)
            _print_stats(stats)

    # Final save
    with open(progress_file, 'w') as f:
        json.dump(progress, f, indent=2)

    # Save final summary
    summary_file = OUTPUT_DIR / "download_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(stats, f, indent=2)

    logger.info("\n" + "=" * 60)
    logger.info("FINAL RESULTS")
    logger.info("=" * 60)
    _print_stats(stats)


def _print_stats(stats):
    total = stats["total"]
    logger.info(f"  Processed:           {stats['processed']}/{total}")
    logger.info(f"  Found in Europe PMC: {stats['found_in_europepmc']}/{stats['processed']} ({100*stats['found_in_europepmc']/max(1,stats['processed']):.1f}%)")
    logger.info(f"  Has PMCID:           {stats['has_pmcid']}/{stats['processed']} ({100*stats['has_pmcid']/max(1,stats['processed']):.1f}%)")
    logger.info(f"  Full-text XML:       {stats['fulltext_xml_downloaded']}/{stats['processed']} ({100*stats['fulltext_xml_downloaded']/max(1,stats['processed']):.1f}%)")
    logger.info(f"  With supplements:    {stats['supplements_found']}/{stats['processed']}")
    logger.info(f"  Suppl files saved:   {stats['supplements_downloaded']}")
    logger.info(f"  Metadata only:       {stats['metadata_only']}")
    logger.info(f"  Not found:           {stats['not_found']}")
    logger.info(f"  Errors:              {stats['errors']}")


if __name__ == "__main__":
    main()
