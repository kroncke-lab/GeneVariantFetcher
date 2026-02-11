#!/usr/bin/env python3
"""
Refined Metadata Detection - Distinguish real PDFs from XML stubs
"""

import os
import re
import subprocess
import sys
from pathlib import Path


def detect_real_pdf(filepath):
    """Distinguish real PDFs from XML stubs or metadata files"""
    if not filepath.exists():
        return {"status": "missing"}

    result = {
        "status": "stub",
        "file_size": filepath.stat().st_size,
        "type": "unknown",
        "pages": 0,
    }

    try:
        # Read first few bytes to detect file type
        with open(filepath, "rb") as f:
            header = f.read(100).decode("utf-8", errors="ignore")

        # Check for PDF signature %PDF
        is_pdf_signature = b"%PDF" in open(filepath, "rb").read(10)

        # Check for XML content signifying stubs
        xml_indicators = [
            "xml version",
            "pmc-articleset",
            "article-meta",
            "journal-meta",
            "<article",
            "<doi>",
            "<pmid>",
            "<pmc-id",
        ]

        is_xml_stub = any(indicator in header.lower() for indicator in xml_indicators)

        # Check actual PDF content using pdfinfo
        try:
            pdf_result = subprocess.run(
                ["pdfinfo", str(filepath)], capture_output=True, text=True, timeout=10
            )

            if pdf_result.returncode == 0:
                # Real PDF - extract info
                pages_match = re.search(r"Pages:\s*(\d+)", pdf_result.stdout)
                title_match = re.search(r"Title:\s*(.+)", pdf_result.stdout)

                result.update(
                    {
                        "status": "real_pdf",
                        "type": "pdf",
                        "pages": int(pages_match.group(1)) if pages_match else 0,
                        "title": title_match.group(1).strip()
                        if title_match
                        else "No title",
                        "info": pdf_result.stdout,
                    }
                )
            else:
                # Not a valid PDF
                if is_xml_stub:
                    result["type"] = "xml_stub"

        except (FileNotFoundError, subprocess.TimeoutExpired):
            if is_xml_stub:
                result["status"] = "xml_stub"
            elif is_pdf_signature:
                result["status"] = "corrupted_pdf"
            elif filepath.suffix.lower() == ".pdf":
                result["status"] = "wrong_extension"
            else:
                result["status"] = "unknown"

    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)

    return result


def scan_all_downloads(base_path="/mnt/temp2/kronckbm/gvf_output"):
    """Scan all downloads and categorize"""
    base = Path(base_path)

    categories = {
        "real_pdfs": [],
        "xml_stubs": [],
        "missing": [],
        "corrupted": [],
        "other": [],
    }

    pdf_files = list(base.rglob("*.pdf"))
    print(f"Found {len(pdf_files)} .pdf files for analysis")

    for pdf_file in pdf_files:
        detection = detect_real_pdf(pdf_file)

        entry = {
            "path": str(pdf_file),
            "stem_pmids": re.findall(r"(\d+)", str(pdf_file.stem)),
            **detection,
        }

        # Determine category
        if detection["status"] == "real_pdf":
            categories["real_pdfs"].append(entry)
        elif detection["status"] == "xml_stub":
            categories["xml_stubs"].append(entry)
        elif detection["status"] == "missing":
            categories["missing"].append(entry)
        elif detection["status"] == "corrupted_pdf":
            categories["corrupted"].append(entry)
        else:
            categories["other"].append(entry)

    return categories


def analyze_venue_sources(base_path="/mnt/temp2/kronckbm/gvf_output"):
    """Analyze which sources are giving real vs stub content"""
    from collections import defaultdict

    categories = scan_all_downloads(base_path)

    venue_analysis = defaultdict(lambda: {"real": 0, "stubs": 0, "other": 0})

    # Analyze real PDFs
    for pdf_info in categories["real_pdfs"]:
        path = Path(pdf_info["path"])
        # Try to determine venue based on path structure
        venue = "unknown"
        if "pmc" in str(path).lower() or "pmc" in pdf_info.get("title", ""):
            venue = "PMC"
        elif "springer" in str(path).lower():
            venue = "Springer"
        elif "wiley" in str(path).lower():
            venue = "Wiley"
        elif "elsevier" in str(path).lower():
            venue = "Elsevier"

        venue_analysis[venue]["real"] += 1

    # Analyze XML stubs
    for stub_info in categories["xml_stubs"]:
        path = Path(stub_info["path"])
        venue = "unknown"
        if "pmc" in str(path).lower():
            venue = "PMC"
        elif "base" in str(path).lower():
            venue = "BASE"
        else:
            venue = "generic"

        venue_analysis[venue]["stubs"] += 1

    return dict(venue_analysis), categories


def main():
    """Run comprehensive analysis"""
    print("=" * 80)
    print("GVF CONTENT ANALYSIS - REAL PDFS VS XML STUBS")
    print("=" * 80)

    # Run analysis
    venue_analysis, categories = analyze_venue_sources()

    print("\nCategorization Results:")
    print(f"Real PDFs: {len(categories['real_pdfs'])}")
    print(f"XML Stubs: {len(categories['xml_stubs'])}")
    print(
        f"Missing/Corrupted: {len(categories['missing']) + len(categories['corrupted'])}"
    )
    print(f"Other types: {len(categories['other'])}")

    print("\nVenue Analysis:")
    for venue, counts in venue_analysis.items():
        total = sum(counts.values())
        real_rate = (counts["real"] / total * 100) if total > 0 else 0
        print(
            f"{venue}: {counts['real']} real, {counts['stubs']} stubs ({real_rate:.1f}% success)"
        )

    # Show examples
    if categories["real_pdfs"]:
        print("\nFirst 3 real PDFs:")
        for pdf in categories["real_pdfs"][:3]:
            print(
                f"  {Path(pdf['path']).name}: {pdf['pages']} pages, {pdf['file_size']} bytes"
            )

    if categories["xml_stubs"]:
        print("\nFirst 3 XML stubs:")
        for stub in categories["xml_stubs"][:3]:
            print(f"  {Path(stub['path']).name}: {stub['file_size']} bytes")

    return categories


if __name__ == "__main__":
    main()
