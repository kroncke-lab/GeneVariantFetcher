#!/usr/bin/env python3
"""
Institutional paper downloader - simplified GVF integration
Explicit Vanderbilt institutional access via GVF APIs
"""

import os
import json
import time
import hashlib
from datetime import datetime
from pathlib import Path

# Import GVF APIs directly
from genevariantfetcher.provider_api.elsevier_api import ElsevierAPI
from genevariantfetcher.utils.file_utils import save_article_pdf


def main():
    output_dir = Path(
        "/mnt/temp2/kronckbm/gvf_output/institutional_downloads/2026-02-07"
    )
    target_pmids_file = output_dir / "target_pmids.txt"
    tracker_path = output_dir / "paper_acquisition_tracker.md"

    # Initialize Elsevier API for institutional access
    elsevier_api = ElsevierAPI()

    # Load target PMIDs
    with open(target_pmids_file) as f:
        pmids = [line.strip() for line in f if line.strip()]

    print(f"Starting institutional download for {len(pmids)} Elsevier papers")

    success_count = 0
    failure_count = 0
    downloaded_hashes = set()

    # Update tracker for start
    with open(tracker_path, "a") as f:
        f.write(f"\n### Starting Downloads @ {datetime.now().strftime('%H:%M:%S')}\n")

    for idx, pmid in enumerate(pmids, 1):
        print(f"[{idx}/{len(pmids)}] Processing PMID: {pmid}")

        try:
            # Attempt Elsevier download
            content = elsevier_api.download_article(pmid)
            if content:
                # Check for duplicates
                content_hash = hashlib.md5(content).hexdigest()
                if content_hash in downloaded_hashes:
                    print(f"  Duplicate detected for PMID {pmid}")
                    with open(tracker_path, "a") as f:
                        f.write(
                            f"| {datetime.now().strftime('%H:%M')} | Elsevier | {pmid} | DUPLICATE | Skipped |"
                        )
                    continue

                downloaded_hashes.add(content_hash)

                # Save PDF
                filename = f"{pmid}_elsevier.pdf"
                filepath = output_dir / filename
                with open(filepath, "wb") as f:
                    f.write(content)

                # Save metadata
                metadata = {
                    "pmid": pmid,
                    "provider": "Elsevier",
                    "method": "institutional_access",
                    "timestamp": datetime.now().isoformat(),
                    "file_hash": content_hash,
                    "file_size": len(content),
                }

                metadata_path = filepath.with_suffix(".json")
                with open(metadata_path, "w") as f:
                    json.dump(metadata, f, indent=2)

                success_count += 1
                print(f"  ✓ Downloaded ({len(content)} bytes)")

                # Update tracker
                with open(tracker_path, "a") as f:
                    f.write(
                        f"| {datetime.now().strftime('%H:%M')} | Elsevier | {pmid} | SUCCESS | {len(content)} bytes downloaded |\n"
                    )

                # Checkpoint every 25 papers
                if success_count % 25 == 0:
                    with open(tracker_path, "a") as f:
                        f.write(
                            f"| {datetime.now().strftime('%H:%M')} | SYSTEM | CHECKPOINT | PROGRESS | {success_count}/40 successful |\n"
                        )
            else:
                print(f"  ✗ Failed - no content retrieved")
                with open(tracker_path, "a") as f:
                    f.write(
                        f"| {datetime.now().strftime('%H:%M')} | Elsevier | {pmid} | FAILED | No content from API |\n"
                    )
                failure_count += 1

        except Exception as e:
            print(f"  ✗ Error: {e}")
            with open(tracker_path, "a") as f:
                f.write(
                    f"| {datetime.now().strftime('%H:%M')} | Elsevier | {pmid} | ERROR | {str(e)} |\n"
                )
            failure_count += 1

        # Brief pause
        time.sleep(1)

    # Final summary
    with open(tracker_path, "a") as f:
        f.write(f"\n### Download Complete @ {datetime.now().strftime('%H:%M:%S')}\n")
        f.write(f"**Total Downloads**: {success_count}\n")
        f.write(f"**Failures**: {failure_count}\n")

    print(f"\nDownload complete!")
    print(f"Success: {success_count}")
    print(f"Failures: {failure_count}")


if __name__ == "__main__":
    main()
