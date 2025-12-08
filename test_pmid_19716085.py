#!/usr/bin/env python3
"""Test harvesting and extraction for PMID 19716085 (KCNH2)

This script tests:
1. Harvesting full-text and supplements
2. Triage system to determine if paper has clinical data
3. Patient extraction to get variant and individual data
"""

import sys
import os
import json
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from harvesting import PMCHarvester
from harvesting.pmc_supp_harvest import PMCSupplementClient
from pipeline.filters import ClinicalDataTriageFilter
from pipeline.extraction import ExpertExtractor
from utils.models import Paper
from dotenv import load_dotenv

# Load environment variables
load_dotenv()


def extract_title_and_abstract_from_content(content: str) -> tuple:
    """Extract title and abstract from the harvested markdown content."""
    title = None
    abstract = None

    lines = content.splitlines()

    # Look for the title (usually the second ## heading after # MAIN TEXT)
    in_main_text = False
    for i, line in enumerate(lines):
        if line.strip() == "# MAIN TEXT":
            in_main_text = True
            continue
        if in_main_text and line.startswith("## "):
            title = line[3:].strip()
            break

    # Look for abstract section
    in_abstract = False
    abstract_lines = []
    for line in lines:
        if line.strip() == "### Abstract":
            in_abstract = True
            continue
        if in_abstract:
            if line.startswith("#"):
                break
            abstract_lines.append(line)

    abstract = " ".join(abstract_lines).strip()

    return title, abstract


def main():
    pmid = "19716085"
    gene = "KCNH2"
    output_dir = Path("test_kcnh2_19716085")

    print(f"Testing PMID: {pmid}")
    print(f"Gene: {gene}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)

    # Check for API key
    api_key = os.getenv('AI_INTEGRATIONS_OPENAI_API_KEY') or os.getenv('OPENAI_API_KEY')
    if not api_key:
        print("\n⚠️  Warning: No OpenAI API key found!")
        print("   Set AI_INTEGRATIONS_OPENAI_API_KEY or OPENAI_API_KEY for triage/extraction.")
        print("   Continuing with harvesting only...\n")
        run_extraction = False
    else:
        run_extraction = True

    # ========================================
    # STEP 1: Harvest full-text and supplements
    # ========================================
    print("\n" + "=" * 60)
    print("STEP 1: HARVESTING FULL-TEXT AND SUPPLEMENTS")
    print("=" * 60)

    supplement_client = PMCSupplementClient()
    try:
        pmcid = supplement_client.pmid_to_pmcid(pmid)
        print(f"Resolved PMID {pmid} to PMCID {pmcid}")
    except ValueError as exc:
        print(f"\n✗ Failed to resolve PMCID for PMID {pmid}: {exc}")
        return

    supp_output_dir = output_dir / f"{pmid}_supplements_api"
    print(f"\nAttempting BioC/JATS supplement harvest to: {supp_output_dir}")
    supplements = supplement_client.get_all_supplement_text(pmcid, supp_output_dir)
    extracted_count = sum(1 for s in supplements if s.text)
    print(
        f"✓ Retrieved {len(supplements)} supplements via API pipeline; "
        f"{extracted_count} with extractable text"
    )
    for supp in supplements:
        label = supp.label or f"Supplement {supp.index}"
        print(f"  - {label}: {supp.filename or supp.href} (source={supp.source})")

    harvester = PMCHarvester(output_dir=str(output_dir))
    success, result = harvester.process_pmid(pmid)

    if not success:
        print(f"\n✗ FAILED: {result}")
        return

    print(f"\n✓ SUCCESS: Created {result}")

    # Read content
    with open(result, 'r') as f:
        content = f.read()

    print(f"\nContent length: {len(content):,} characters")
    print(f"Line count: {len(content.splitlines()):,} lines")

    # Check for supplements
    if "## Supplementary" in content or "supplement" in content.lower():
        print("✓ Contains supplementary material references")

    # Show section headers
    print("\nSection headers found:")
    for line in content.splitlines():
        if line.startswith('#'):
            print(f"  {line}")

    # Check for supplement files
    supp_dir = output_dir / f"{pmid}_supplements"
    if supp_dir.exists():
        supp_files = list(supp_dir.iterdir())
        print(f"\n✓ Supplement directory exists with {len(supp_files)} files:")
        for f in supp_files:
            print(f"    - {f.name}")

    if not run_extraction:
        print("\n⚠️  Skipping triage and extraction (no API key).")
        return

    # Extract title and abstract from content
    title, abstract = extract_title_and_abstract_from_content(content)

    if not title:
        title = "Spectrum and prevalence of mutations from the first 2,500 consecutive unrelated patients referred for the FAMILION® long QT syndrome genetic test"
    if not abstract:
        # Fallback abstract from the paper
        abstract = content[:2500]

    print(f"\nExtracted title: {title[:80]}...")
    print(f"Abstract length: {len(abstract)} chars")

    # ========================================
    # STEP 2: Triage - Does this paper have clinical data?
    # ========================================
    print("\n" + "=" * 60)
    print("STEP 2: TRIAGE - CHECKING FOR CLINICAL DATA")
    print("=" * 60)

    triage_filter = ClinicalDataTriageFilter()

    triage_result = triage_filter.triage(
        title=title,
        abstract=abstract,
        gene=gene,
        pmid=pmid
    )

    print(f"\nTriage Result:")
    print(f"  Decision: {triage_result['decision']}")
    print(f"  Confidence: {triage_result['confidence']:.2f}")
    print(f"  Reason: {triage_result['reason']}")

    if triage_result['decision'] == 'DROP':
        print("\n⚠️  Paper dropped by triage - no original clinical data detected.")
        print("    (This may be a review, animal study, or purely methodological paper)")
        return

    print("\n✓ Paper passed triage - proceeding to extraction...")

    # ========================================
    # STEP 3: Extract patient/variant information
    # ========================================
    print("\n" + "=" * 60)
    print("STEP 3: EXTRACTING PATIENT AND VARIANT INFORMATION")
    print("=" * 60)

    # Create Paper object for extraction
    paper = Paper(
        pmid=pmid,
        title=title,
        abstract=abstract,
        full_text=content,
        gene_symbol=gene
    )

    extractor = ExpertExtractor()
    extraction_result = extractor.extract(paper)

    if not extraction_result.success:
        print(f"\n✗ Extraction failed: {extraction_result.error}")
        return

    print(f"\n✓ Extraction successful!")
    print(f"  Model used: {extraction_result.model_used}")

    # Display extraction results
    data = extraction_result.extracted_data

    print("\n" + "-" * 50)
    print("EXTRACTION SUMMARY:")
    print("-" * 50)

    if 'paper_metadata' in data:
        meta = data['paper_metadata']
        print(f"\nPaper: {meta.get('title', 'N/A')[:60]}...")
        print(f"Summary: {meta.get('extraction_summary', 'N/A')}")

    if 'extraction_metadata' in data:
        emeta = data['extraction_metadata']
        print(f"\nTotal variants found: {emeta.get('total_variants_found', 0)}")
        print(f"Confidence: {emeta.get('extraction_confidence', 'N/A')}")
        if emeta.get('challenges'):
            print(f"Challenges: {', '.join(emeta['challenges'])}")

    # Show variants
    variants = data.get('variants', [])
    print(f"\n" + "-" * 50)
    print(f"VARIANTS EXTRACTED: {len(variants)}")
    print("-" * 50)

    for i, variant in enumerate(variants[:10], 1):  # Show first 10
        print(f"\n  Variant {i}:")
        print(f"    Gene: {variant.get('gene_symbol', 'N/A')}")
        print(f"    cDNA: {variant.get('cdna_notation', 'N/A')}")
        print(f"    Protein: {variant.get('protein_notation', 'N/A')}")
        print(f"    Significance: {variant.get('clinical_significance', 'N/A')}")

        # Penetrance data
        pdata = variant.get('penetrance_data', {})
        if pdata:
            print(f"    Penetrance:")
            print(f"      Total carriers: {pdata.get('total_carriers_observed', 'N/A')}")
            print(f"      Affected: {pdata.get('affected_count', 'N/A')}")
            print(f"      Unaffected: {pdata.get('unaffected_count', 'N/A')}")

        # Individual records
        individuals = variant.get('individual_records', [])
        if individuals:
            print(f"    Individual records: {len(individuals)}")
            for ind in individuals[:3]:  # Show first 3
                print(f"      - {ind.get('individual_id', '?')}: "
                      f"{ind.get('affected_status', '?')}, "
                      f"age={ind.get('age_at_evaluation', '?')}, "
                      f"sex={ind.get('sex', '?')}")
            if len(individuals) > 3:
                print(f"      ... and {len(individuals) - 3} more")

    if len(variants) > 10:
        print(f"\n  ... and {len(variants) - 10} more variants")

    # Save full extraction to JSON
    output_json = output_dir / f"{pmid}_extraction.json"
    with open(output_json, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\n✓ Full extraction saved to: {output_json}")

    # Summary stats
    print("\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"  PMID: {pmid}")
    print(f"  Gene: {gene}")
    print(f"  Triage: {triage_result['decision']} (confidence: {triage_result['confidence']:.2f})")
    print(f"  Variants extracted: {len(variants)}")

    total_individuals = sum(len(v.get('individual_records', [])) for v in variants)
    print(f"  Individual patient records: {total_individuals}")


if __name__ == "__main__":
    main()
