#!/usr/bin/env python3
"""
Phase 1: Download Verification Summary
Run after batch download to validate results and provide final report
"""

import json
import csv
from pathlib import Path
from datetime import datetime

def summarize_download_results():
    output_dir = Path("/mnt/temp2/kronckbm/gvf_output/download_batch_20260206")
    summary_path = Path("/mnt/temp2/kronckbm/gvf_output/download_log_20260206.md")
    
    print("=" * 80)
    print("GVF PHASE 1: DOWNLOAD VERIFICATION")
    print("=" * 80)
    
    # Check for existence of key files
    manifest_file = output_dir / "manifest.json"
    success_log = output_dir / "successful_downloads.csv"
    paywall_log = output_dir / "paywalled_missing.csv"
    
    success_pmids = []
    failed_pmids = []
    
    # Parse success log
    if success_log.exists():
        with open(success_log) as f:
            reader = csv.DictReader(f)
            for row in reader:
                success_pmids.append({
                    'pmid': row['PMID'],
                    'pmcid': row['PMCID'],
                    'supplements': int(row['Supplements_Downloaded'])
                })
    
    # Parse paywall log
    if paywall_log.exists():
        with open(paywall_log) as f:
            reader = csv.DictReader(f)
            for row in reader:
                failed_pmids.append({
                    'pmid': row['PMID'],
                    'reason': row['Reason'],
                    'url': row['URL']
                })
    
    # Check downloaded files
    full_context_files = list(output_dir.glob("*_FULL_CONTEXT.md"))
    figure_dirs = list(output_dir.glob("*_figures"))
    supplement_dirs = list(output_dir.glob("*_supplements"))
    
    # File sizes for validation
    paper_details = []
    for f in full_context_files:
        pmid = f.stem.replace('_FULL_CONTEXT', '')
        size = f.stat().st_size
        has_figures = any(pmid in str(d) for d in figure_dirs)
        has_supplements = any(pmid in str(d) for d in supplement_dirs)
        
        paper_details.append({
            'pmid': pmid,
            'size': size,
            'size_kb': round(size/1024, 1),
            'has_figures': has_figures,
            'has_supplements': has_supplements
        })
    
    # Summary statistics
    total_processed = len(success_pmids) + len(failed_pmids)
    success_rate = len(success_pmids) / total_processed * 100 if total_processed > 0 else 0
    
    print(f"Total PMIDs processed: {total_processed}")
    print(f"Successfully downloaded: {len(success_pmids)}")
    print(f"Failed downloads: {len(failed_pmids)}")
    print(f"Success rate: {success_rate:.1f}%")
    print()
    
    print("📋 SUCCESSFULLY DOWNLOADED PAPERS:")
    for paper in paper_details:
        print(f"  - PMID {paper['pmid']}: {paper['size_kb']}KB (Figures: {paper['has_figures']}, Supplements: {paper['has_supplements']})")
    
    print()
    print("❌ FAILED DOWNLOADS:")
    for fail in failed_pmids[:5]:  # Show first 5
        print(f"  - PMID {fail['pmid']}: {fail['reason'][:80]}...")
    
    # Save detailed log
    with open(summary_path, "a") as f:
        f.write("\n" + "=" * 50 + "\n")
        f.write("DOWNLOAD VERIFICATION REPORT\n")
        f.write(datetime.now().strftime("Generated: %Y-%m-%d %H:%M:%S\n"))
        f.write("=" * 50 + "\n\n")
        
        f.write("## Downloaded Paper Details:\n")
        for paper in paper_details:
            f.write(f"- **PMID {paper['pmid']}**: {paper['size_kb']}KB, "
                   f"Figures: {paper['has_figures']}, Supplements: {paper['has_supplements']}\n")
        
        f.write(f"\n## Summary:\n")
        f.write(f"- **Total Processed**: {total_processed}\n")
        f.write(f"- **Successfully Downloaded**: {len(success_pmids)}\n")
        f.write(f"- **Failed Downloads**: {len(failed_pmids)}\n")
        f.write(f"- **Success Rate**: {success_rate:.1f}%\n")
    
    print()
    print(f"✅ Verification complete. Results appended to {summary_path}")
    return len(success_pmids), len(failed_pmids)

if __name__ == "__main__":
    success, failed = summarize_download_results()
    print(f"\nFinal validation: {success} successful, {failed} failed")