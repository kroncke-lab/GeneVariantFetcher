#!/usr/bin/env python3
"""Fast GVF extraction and recall calculation."""
import os
import json
import glob
import re
import pandas as pd
import pdfplumber
from concurrent.futures import ProcessPoolExecutor, as_completed

# Configuration
DOWNLOAD_DIR = "/mnt/temp2/kronckbm/gvf_output/verified_downloads_20260208"
GOLD_STANDARD = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
OUTPUT_DIR = "/mnt/temp2/kronckbm/gvf_output/extraction_results_20260208"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def extract_variants_from_text(text):
    """Extract variants from text."""
    variants = set()
    patterns = [
        r'\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})([ACDEFGHIKLMNPQRSTVWY*X])\b',
        r'p\.([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})([ACDEFGHIKLMNPQRSTVWY*X])\b',
        r'\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})(fs[X\d]*)\b',
        r'\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})(del)\b',
    ]
    
    for pattern in patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE):
            aa1 = match.group(1).upper()
            pos = match.group(2)
            aa2 = match.group(3).upper()
            
            aa2 = aa2.replace('*', 'X')
            if 'FS' in aa2.upper():
                aa2 = 'FSX'
            if 'DEL' in aa2.upper():
                aa2 = 'DEL'
            
            try:
                pos_int = int(pos)
                if 1 <= pos_int <= 1200:
                    variants.add(f"{aa1}{pos}{aa2}")
            except ValueError:
                pass
    return variants

def process_pdf(pdf_path):
    """Process a single PDF."""
    try:
        text = ""
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages[:50]:  # Limit pages for speed
                page_text = page.extract_text()
                if page_text:
                    text += page_text + "\n"
        return extract_variants_from_text(text)
    except Exception as e:
        return set()

# Load gold standard
print("Loading gold standard...")
df_gold = pd.read_excel(GOLD_STANDARD)
gold_variants = set()
for v in df_gold['Variant'].dropna():
    v_norm = str(v).upper().strip()
    v_norm = v_norm.replace('*', 'X')
    v_norm = re.sub(r'FSX?\d*$', 'FSX', v_norm)
    gold_variants.add(v_norm)

print(f"Gold standard variants: {len(gold_variants)}")

# Process PDFs
pdfs = glob.glob(os.path.join(DOWNLOAD_DIR, "*.pdf"))
print(f"Processing {len(pdfs)} PDFs...")

all_extracted = set()
processed = 0
errors = 0

for i, pdf_path in enumerate(pdfs):
    if (i + 1) % 20 == 0:
        print(f"  Progress: {i+1}/{len(pdfs)}")
    try:
        variants = process_pdf(pdf_path)
        all_extracted.update(variants)
        processed += 1
    except Exception as e:
        errors += 1

# Calculate recall
matches = gold_variants & all_extracted
novel = all_extracted - gold_variants
missed = gold_variants - all_extracted

print(f"\n{'='*60}")
print("EXTRACTION RESULTS")
print(f"{'='*60}")
print(f"PDFs processed: {processed}/{len(pdfs)} (errors: {errors})")
print(f"Gold standard variants: {len(gold_variants)}")
print(f"Extracted variants: {len(all_extracted)}")
print(f"Matched: {len(matches)}")
print(f"Novel: {len(novel)}")
print(f"Missed: {len(missed)}")
print(f"\nVARIANT EXTRACTION RECALL: {len(matches)}/{len(gold_variants)} = {100*len(matches)/len(gold_variants):.1f}%")
print(f"{'='*60}")

# Save results
results = {
    "pdfs_processed": processed,
    "gold_variants": len(gold_variants),
    "extracted_variants": len(all_extracted),
    "matched_variants": len(matches),
    "recall_percent": round(100*len(matches)/len(gold_variants), 1),
    "matched_list": sorted(list(matches))[:50],
    "sample_novel": sorted(list(novel))[:30],
    "sample_missed": sorted(list(missed))[:30]
}

with open(os.path.join(OUTPUT_DIR, "extraction_recall_results.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nSample matched: {sorted(list(matches))[:15]}")
print(f"Sample missed: {sorted(list(missed))[:15]}")
