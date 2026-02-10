#!/usr/bin/env python3
"""
Download supplementary files from Europe PMC for papers that have PMCIDs.

The Europe PMC supplementary files endpoint returns a ZIP file (not JSON),
so we download and extract it for each PMCID.

Usage:
    python scripts/europepmc_download_supplements.py
"""

import sys
import os
import json
import time
import io
import zipfile
import logging
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)

OUTPUT_DIR = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2/europepmc_downloads")
RATE_LIMIT_DELAY = 0.5
BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"


def download_supplements_zip(session, pmcid: str, output_dir: Path) -> dict:
    """Download supplementary files ZIP from Europe PMC and extract."""
    result = {"pmcid": pmcid, "files": [], "error": None}

    url = f"{BASE_URL}/{pmcid}/supplementaryFiles"

    try:
        resp = session.get(url, timeout=60)

        if resp.status_code == 404:
            result["error"] = "No supplements available (404)"
            return result

        resp.raise_for_status()

        content_type = resp.headers.get('content-type', '')

        if 'zip' in content_type or resp.content[:2] == b'PK':
            # It's a ZIP file
            try:
                zf = zipfile.ZipFile(io.BytesIO(resp.content))
                supp_dir = output_dir / "supplements"
                supp_dir.mkdir(exist_ok=True)

                for name in zf.namelist():
                    # Skip directory entries
                    if name.endswith('/'):
                        continue
                    # Extract file
                    data = zf.read(name)
                    # Use just the filename, not the path
                    basename = Path(name).name
                    out_file = supp_dir / basename
                    with open(out_file, 'wb') as f:
                        f.write(data)
                    result["files"].append(str(out_file))

                logger.info(f"  Extracted {len(result['files'])} supplement files from ZIP")
            except zipfile.BadZipFile:
                result["error"] = "Invalid ZIP file received"
        elif 'xml' in content_type or 'html' in content_type:
            # Some endpoints return XML/HTML error pages
            if '<error>' in resp.text.lower() or '404' in resp.text[:200]:
                result["error"] = "No supplements (XML error response)"
            else:
                # Save as-is
                supp_dir = output_dir / "supplements"
                supp_dir.mkdir(exist_ok=True)
                out_file = supp_dir / f"{pmcid}_supplements.xml"
                with open(out_file, 'w') as f:
                    f.write(resp.text)
                result["files"].append(str(out_file))
        else:
            result["error"] = f"Unexpected content-type: {content_type}"

    except requests.exceptions.HTTPError as e:
        result["error"] = str(e)
    except Exception as e:
        result["error"] = str(e)

    return result


def main():
    # Load progress from initial download
    progress_file = OUTPUT_DIR / "download_progress.json"
    with open(progress_file) as f:
        progress = json.load(f)

    # Find all PMIDs with PMCIDs
    pmcid_papers = []
    for pmid, result in progress["results"].items():
        pmcid = result.get("pmcid")
        if pmcid:
            pmcid_papers.append((pmid, pmcid))

    logger.info(f"Found {len(pmcid_papers)} papers with PMCIDs to check for supplements")

    session = requests.Session()
    session.headers.update({
        'User-Agent': 'GeneVariantFetcher/1.0 (Bot for clinical variant extraction)',
        'Accept': '*/*',
    })

    stats = {
        "total": len(pmcid_papers),
        "with_supplements": 0,
        "total_files": 0,
        "no_supplements": 0,
        "errors": 0,
    }

    supplement_results = {}

    for i, (pmid, pmcid) in enumerate(pmcid_papers, 1):
        logger.info(f"[{i}/{len(pmcid_papers)}] PMID {pmid} ({pmcid})")

        pmid_dir = OUTPUT_DIR / pmid
        pmid_dir.mkdir(exist_ok=True)

        result = download_supplements_zip(session, pmcid, pmid_dir)
        supplement_results[pmid] = result

        if result["files"]:
            stats["with_supplements"] += 1
            stats["total_files"] += len(result["files"])
            logger.info(f"  -> {len(result['files'])} files downloaded")
        elif result["error"]:
            if "404" in str(result["error"]) or "No supplements" in str(result["error"]):
                stats["no_supplements"] += 1
                logger.info(f"  -> No supplements available")
            else:
                stats["errors"] += 1
                logger.warning(f"  -> Error: {result['error']}")
        else:
            stats["no_supplements"] += 1

        time.sleep(RATE_LIMIT_DELAY)

    # Save supplement results
    supp_summary = OUTPUT_DIR / "supplement_download_results.json"
    with open(supp_summary, 'w') as f:
        json.dump({"stats": stats, "results": supplement_results}, f, indent=2)

    logger.info("\n" + "=" * 60)
    logger.info("SUPPLEMENT DOWNLOAD RESULTS")
    logger.info("=" * 60)
    logger.info(f"  Papers checked:      {stats['total']}")
    logger.info(f"  With supplements:    {stats['with_supplements']}")
    logger.info(f"  Total files:         {stats['total_files']}")
    logger.info(f"  No supplements:      {stats['no_supplements']}")
    logger.info(f"  Errors:              {stats['errors']}")


if __name__ == "__main__":
    main()
