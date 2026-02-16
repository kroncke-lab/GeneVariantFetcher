#!/usr/bin/env python3
"""
Run full GVF download batch for KCNH2.

This script runs the download-only phase (no extraction) for all PMIDs.
"""

import os
import sys
from datetime import datetime
from pathlib import Path

# Ensure GVF is importable
sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv
load_dotenv()

from cli import automated_variant_extraction_workflow

def main():
    gene = "KCNH2"
    email = os.environ.get("NCBI_EMAIL", "brett.kroncke@gmail.com")
    
    # Create timestamped output dir
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"/mnt/temp2/kronckbm/gvf_output/{gene}/{timestamp}"
    
    print(f"=" * 60)
    print(f"GVF Full Batch Download")
    print(f"=" * 60)
    print(f"Gene: {gene}")
    print(f"Output: {output_dir}")
    print(f"Started: {datetime.now().isoformat()}")
    print(f"=" * 60)
    
    # Run the workflow with download-only settings
    automated_variant_extraction_workflow(
        gene_symbol=gene,
        email=email,
        output_dir=output_dir,
        max_pmids=10000,  # High limit to get all
        max_papers_to_download=10000,  # Download all available
        tier_threshold=0,  # Don't do extraction (download only)
        use_clinical_triage=True,  # Use clinical triage filter
        scout_first=True,  # Run data scout
    )
    
    print(f"\n{'=' * 60}")
    print(f"Completed: {datetime.now().isoformat()}")
    print(f"Output: {output_dir}")
    print(f"=" * 60)

if __name__ == "__main__":
    main()
