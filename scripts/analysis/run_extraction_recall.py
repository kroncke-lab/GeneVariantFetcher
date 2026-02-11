#!/usr/bin/env python3
"""Run GVF extraction on verified downloads and calculate recall against gold standard."""

import os
import sys
import json
import glob
import re
from pathlib import Path

# Add project to path
sys.path.insert(0, "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher")

import pandas as pd

# Configuration
DOWNLOAD_DIR = "/mnt/temp2/kronckbm/gvf_output/verified_downloads_20260208"
GOLD_STANDARD = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
OUTPUT_DIR = "/mnt/temp2/kronckbm/gvf_output/extraction_results_20260208"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load gold standard
df_gold = pd.read_excel(GOLD_STANDARD)
gold_variants = set(df_gold["Variant"].dropna().str.strip().str.upper().unique())

print(f"Gold Standard: {len(gold_variants)} unique variants")
print(f"Sample gold variants: {list(gold_variants)[:10]}")


# Simple variant extraction from text using regex patterns
def extract_variants_from_text(text, gene="KCNH2"):
    """Extract KCNH2 variants from text using regex patterns."""
    variants = set()

    # Amino acid substitution patterns (e.g., G601S, p.G601S, p.Gly601Ser)
    patterns = [
        # Short notation: G601S, A561V
        r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})([ACDEFGHIKLMNPQRSTVWY*X])\b",
        # p. notation: p.G601S, p.Gly601Ser
        r"p\.([ACDEFGHIKLMNPQRSTVWY](?:[a-z]{2})?)(\d{2,4})([ACDEFGHIKLMNPQRSTVWY*X](?:[a-z]{2})?)\b",
        # Frameshift: R1005fs, R1005fsX
        r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})(fs[X\d]*)\b",
        # Deletions: Y475del, A79Del
        r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})(del)\b",
        # Nonsense: R863X, R744*
        r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})([X*])\b",
    ]

    for pattern in patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE):
            if len(match.groups()) >= 3:
                # Normalize to short form
                aa1 = match.group(1).upper()
                pos = match.group(2)
                aa2 = match.group(3).upper()

                # Convert 3-letter codes to 1-letter
                aa_map = {
                    "ALA": "A",
                    "CYS": "C",
                    "ASP": "D",
                    "GLU": "E",
                    "PHE": "F",
                    "GLY": "G",
                    "HIS": "H",
                    "ILE": "I",
                    "LYS": "K",
                    "LEU": "L",
                    "MET": "M",
                    "ASN": "N",
                    "PRO": "P",
                    "GLN": "Q",
                    "ARG": "R",
                    "SER": "S",
                    "THR": "T",
                    "VAL": "V",
                    "TRP": "W",
                    "TYR": "Y",
                }

                if len(aa1) == 3:
                    aa1 = aa_map.get(aa1.upper(), aa1[0])
                if len(aa2) == 3:
                    aa2 = aa_map.get(aa2.upper(), aa2[0])

                # Normalize special cases
                aa2 = aa2.replace("*", "X")
                if "FS" in aa2.upper():
                    aa2 = "fsX" if not aa2.endswith("X") else aa2
                if "DEL" in aa2.upper():
                    aa2 = "Del"

                variant = f"{aa1}{pos}{aa2}"

                # Filter to reasonable KCNH2 positions (1-1159)
                try:
                    pos_int = int(pos)
                    if 1 <= pos_int <= 1200:
                        variants.add(variant.upper())
                except ValueError:
                    pass

    return variants


# Process PDFs
try:
    import pdfplumber

    HAS_PDFPLUMBER = True
except ImportError:
    HAS_PDFPLUMBER = False
    print("WARNING: pdfplumber not available, using basic extraction")


def extract_text_from_pdf(pdf_path):
    """Extract text from PDF."""
    text = ""
    if HAS_PDFPLUMBER:
        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page in pdf.pages:
                    page_text = page.extract_text()
                    if page_text:
                        text += page_text + "\n"
        except Exception as e:
            print(f"  Error with pdfplumber: {e}")
    return text


# Process all PDFs
print(f"\nProcessing PDFs from: {DOWNLOAD_DIR}")
pdfs = glob.glob(os.path.join(DOWNLOAD_DIR, "*.pdf"))
print(f"Found {len(pdfs)} PDFs\n")

all_extracted = set()
pmid_variants = {}
extraction_errors = []

for i, pdf_path in enumerate(pdfs):
    basename = os.path.basename(pdf_path)
    pmid = basename.replace(".pdf", "").split("_")[0]

    if (i + 1) % 20 == 0:
        print(f"Processing {i + 1}/{len(pdfs)}: {basename}")

    text = extract_text_from_pdf(pdf_path)
    if not text:
        extraction_errors.append(pmid)
        continue

    variants = extract_variants_from_text(text)
    pmid_variants[pmid] = list(variants)
    all_extracted.update(variants)

# Calculate recall
print(f"\n{'=' * 70}")
print("EXTRACTION RESULTS")
print(f"{'=' * 70}")


# Normalize gold variants for comparison
def normalize_variant(v):
    """Normalize variant for comparison."""
    v = str(v).upper().strip()
    # Remove 'P.' prefix
    if v.startswith("P."):
        v = v[2:]
    # Replace * with X
    v = v.replace("*", "X")
    # Handle fsX variations
    v = re.sub(r"FSX?\d*$", "FSX", v)
    return v


gold_normalized = {normalize_variant(v) for v in gold_variants}
extracted_normalized = {normalize_variant(v) for v in all_extracted}

# Find matches
matches = gold_normalized & extracted_normalized
novel = extracted_normalized - gold_normalized
missed = gold_normalized - extracted_normalized

print(f"\nTotal PDFs processed:           {len(pdfs)}")
print(f"PDFs with extraction errors:    {len(extraction_errors)}")
print(f"\nGold standard variants:         {len(gold_normalized)}")
print(f"Extracted variants (unique):    {len(extracted_normalized)}")
print(f"Matched variants:               {len(matches)}")
print(f"Novel (not in gold):            {len(novel)}")
print(f"Missed (in gold, not extracted): {len(missed)}")

print(f"\n{'=' * 70}")
print(
    f"VARIANT EXTRACTION RECALL:      {len(matches)}/{len(gold_normalized)} = {100 * len(matches) / len(gold_normalized):.1f}%"
)
print(f"{'=' * 70}")

# Sample matched variants
print(f"\nSample matched variants (first 20):")
for v in sorted(matches)[:20]:
    print(f"  {v}")

# Save results
results = {
    "pdfs_processed": len(pdfs),
    "extraction_errors": len(extraction_errors),
    "gold_variants": len(gold_normalized),
    "extracted_variants": len(extracted_normalized),
    "matched_variants": len(matches),
    "novel_variants": len(novel),
    "missed_variants": len(missed),
    "recall_percent": round(100 * len(matches) / len(gold_normalized), 1),
    "matched_list": sorted(matches),
    "novel_list": sorted(list(novel)[:100]),  # First 100
    "missed_list": sorted(list(missed)[:100]),  # First 100
    "pmid_variants": pmid_variants,
}

with open(os.path.join(OUTPUT_DIR, "extraction_recall_results.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {OUTPUT_DIR}/extraction_recall_results.json")
