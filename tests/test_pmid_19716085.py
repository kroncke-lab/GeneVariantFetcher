#!/usr/bin/env python3
"""Test harvesting and extraction for PMID 19716085 (KCNH2)

This module tests:
1. Harvesting full-text and supplements
2. Triage system to determine if paper has clinical data
3. Patient extraction to get variant and individual data

Can be run with pytest or directly as a script.
"""

import os
import json
from pathlib import Path

import pytest

from harvesting import PMCHarvester
from harvesting.pmc_supp_harvest import PMCSupplementClient
from pipeline.filters import ClinicalDataTriageFilter
from pipeline.extraction import ExpertExtractor
from utils.models import Paper
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Test constants
PMID = "19716085"
GENE = "KCNH2"
OUTPUT_DIR = Path("test_kcnh2_19716085")
FULL_TEXT_PATH = OUTPUT_DIR / f"{PMID}_FULL_CONTEXT.md"
DATA_ZONES_PATH = OUTPUT_DIR / f"{PMID}_DATA_ZONES.md"

# Expected title for this paper
EXPECTED_TITLE = "Spectrum and prevalence of mutations from the first 2,500 consecutive unrelated patients referred for the FAMILION®long QT syndrome genetic test"


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


@pytest.fixture(scope="module")
def harvested_content():
    """Fixture to harvest content once per module."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    harvester = PMCHarvester(output_dir=str(OUTPUT_DIR), gene_symbol=GENE)

    # Reuse cached harvest if available to avoid network calls
    if FULL_TEXT_PATH.exists():
        content = FULL_TEXT_PATH.read_text()
        return content, str(FULL_TEXT_PATH), harvester

    # Harvest the full text
    success, result = harvester.process_pmid(PMID)
    assert success, f"Harvesting failed: {result}"

    content = Path(result).read_text()
    return content, result, harvester


class TestPMID19716085Harvesting:
    """Test class for harvesting PMID 19716085."""

    def test_harvest_full_text(self, harvested_content):
        """Test that full-text can be harvested successfully."""
        content, result_path, _ = harvested_content

        assert content, "Content should not be empty"
        assert len(content) > 10000, "Full text should be substantial"
        assert Path(result_path).exists(), "Result file should exist"

    def test_content_has_expected_sections(self, harvested_content):
        """Test that harvested content contains expected sections."""
        content, _, _ = harvested_content

        assert "# MAIN TEXT" in content, "Should have main text section"
        assert "Abstract" in content, "Should have abstract"
        assert "Methods" in content or "METHODS" in content, "Should have methods"
        assert "Results" in content or "RESULTS" in content, "Should have results"

    def test_content_has_supplementary_material(self, harvested_content):
        """Test that supplementary material is referenced."""
        content, _, _ = harvested_content

        has_supplements = (
            "supplement" in content.lower()
            or "Supplementary" in content
            or "SUPPLEMENTAL" in content
        )
        assert has_supplements, "Should reference supplementary material"

    def test_data_zones_created(self, harvested_content):
        """Test that DATA_ZONES.md is created for condensed extraction."""
        content, _, harvester = harvested_content

        # Ensure data zones exist (create if needed)
        if not DATA_ZONES_PATH.exists():
            harvester._run_data_scout(PMID, content)

        assert DATA_ZONES_PATH.exists(), "DATA_ZONES.md should be created"

        zones_content = DATA_ZONES_PATH.read_text()
        assert len(zones_content) > 0, "DATA_ZONES.md should have content"
        # Data zones should be smaller than full content
        assert len(zones_content) <= len(content), "Data zones should be condensed"

    def test_title_extraction(self, harvested_content):
        """Test that title can be extracted from content."""
        content, _, _ = harvested_content

        title, _ = extract_title_and_abstract_from_content(content)
        assert title is not None, "Title should be extractable"
        assert (
            "FAMILION" in title or "long QT" in title.lower()
        ), f"Title should be about LQTS genetic testing, got: {title}"


class TestPMID19716085SupplementAPI:
    """Test supplement harvesting via BioC/JATS API."""

    def test_pmid_to_pmcid_resolution(self):
        """Test PMID to PMCID resolution."""
        client = PMCSupplementClient()
        try:
            pmcid = client.pmid_to_pmcid(PMID)
            assert pmcid.startswith("PMC"), f"Should return PMCID, got: {pmcid}"
        except ValueError as e:
            pytest.skip(f"PMCID resolution unavailable: {e}")

    def test_supplement_retrieval(self):
        """Test that supplement API can be called (may return empty for some papers)."""
        client = PMCSupplementClient()
        try:
            pmcid = client.pmid_to_pmcid(PMID)
            supp_dir = OUTPUT_DIR / f"{PMID}_supplements_api"
            supplements = client.get_all_supplement_text(pmcid, supp_dir)

            # API call succeeded - supplements may be empty for some papers
            assert isinstance(supplements, list), "Should return a list"

            # If supplements were found, verify their structure
            if supplements:
                for supp in supplements:
                    assert supp.pmcid == pmcid
                    assert supp.source in ("bioc", "jats")
        except ValueError as e:
            pytest.skip(f"Supplement retrieval unavailable: {e}")


@pytest.mark.skipif(not os.getenv("OPENAI_API_KEY"), reason="OPENAI_API_KEY not set")
class TestPMID19716085Triage:
    """Test triage for PMID 19716085."""

    def test_triage_keeps_paper(self, harvested_content):
        """Test that triage correctly identifies this paper as having clinical data."""
        content, _, _ = harvested_content

        title, abstract = extract_title_and_abstract_from_content(content)
        if not title:
            title = EXPECTED_TITLE
        if not abstract:
            abstract = content[:2500]

        triage_filter = ClinicalDataTriageFilter()
        result = triage_filter.triage(
            title=title, abstract=abstract, gene=GENE, pmid=PMID
        )

        assert result["decision"] == "KEEP", (
            f"Paper should be kept for extraction, but got: {result['decision']}. "
            f"Reason: {result.get('reason', 'N/A')}"
        )
        assert result["confidence"] >= 0.5, "Confidence should be reasonably high"


@pytest.mark.skipif(not os.getenv("OPENAI_API_KEY"), reason="OPENAI_API_KEY not set")
class TestPMID19716085Extraction:
    """Test extraction for PMID 19716085."""

    def test_extraction_succeeds(self, harvested_content):
        """Test that extraction completes successfully."""
        content, _, _ = harvested_content

        title, abstract = extract_title_and_abstract_from_content(content)
        if not title:
            title = EXPECTED_TITLE
        if not abstract:
            abstract = content[:2500]

        paper = Paper(
            pmid=PMID,
            title=title,
            abstract=abstract,
            full_text=content,
            gene_symbol=GENE,
        )

        extractor = ExpertExtractor(fulltext_dir=str(OUTPUT_DIR))
        result = extractor.extract(paper)

        assert result.success, f"Extraction should succeed: {result.error}"
        assert result.extracted_data is not None, "Should have extracted data"

    def test_extraction_finds_variants(self, harvested_content):
        """Test that extraction finds KCNH2 variants."""
        content, _, _ = harvested_content

        title, abstract = extract_title_and_abstract_from_content(content)
        if not title:
            title = EXPECTED_TITLE
        if not abstract:
            abstract = content[:2500]

        paper = Paper(
            pmid=PMID,
            title=title,
            abstract=abstract,
            full_text=content,
            gene_symbol=GENE,
        )

        extractor = ExpertExtractor(fulltext_dir=str(OUTPUT_DIR))
        result = extractor.extract(paper)

        assert result.success, f"Extraction should succeed: {result.error}"

        data = result.extracted_data
        variants = data.get("variants", [])

        # This paper should have many KCNH2 variants
        assert (
            len(variants) > 50
        ), f"Should find many variants (this paper has ~199), found: {len(variants)}"

        # Check that variants have expected fields
        for v in variants[:5]:
            assert "gene_symbol" in v or "cdna_notation" in v or "protein_notation" in v

    def test_extraction_output_saved(self, harvested_content):
        """Test that extraction output can be saved to JSON."""
        content, _, _ = harvested_content

        title, abstract = extract_title_and_abstract_from_content(content)
        if not title:
            title = EXPECTED_TITLE
        if not abstract:
            abstract = content[:2500]

        paper = Paper(
            pmid=PMID,
            title=title,
            abstract=abstract,
            full_text=content,
            gene_symbol=GENE,
        )

        extractor = ExpertExtractor(fulltext_dir=str(OUTPUT_DIR))
        result = extractor.extract(paper)

        assert result.success

        # Save to JSON
        output_json = OUTPUT_DIR / f"{PMID}_extraction.json"
        with open(output_json, "w") as f:
            json.dump(result.extracted_data, f, indent=2)

        assert output_json.exists(), "Extraction JSON should be saved"

        # Verify it can be loaded back
        with open(output_json) as f:
            loaded = json.load(f)
        assert "variants" in loaded


def main():
    """Run the test as a script with detailed output."""
    print(f"Testing PMID: {PMID}")
    print(f"Gene: {GENE}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("-" * 50)

    # Check for API key
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("\n  Warning: OPENAI_API_KEY not found!")
        print("   Set OPENAI_API_KEY in your .env file for triage/extraction.")
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

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    harvester = PMCHarvester(output_dir=str(OUTPUT_DIR), gene_symbol=GENE)

    # Reuse cached harvest if available to avoid network in restricted runs
    if FULL_TEXT_PATH.exists():
        print(f"Using cached full text: {FULL_TEXT_PATH}")
        success, result = True, str(FULL_TEXT_PATH)
    else:
        # First, use the main harvester which handles both PMC and publisher fallback
        success, result = harvester.process_pmid(PMID)

        # If paper has a PMCID, also try the BioC/JATS supplement harvest
        supplement_client = PMCSupplementClient()
        try:
            pmcid = supplement_client.pmid_to_pmcid(PMID)
            print(f"\nResolved PMID {PMID} to PMCID {pmcid}")

            supp_output_dir = OUTPUT_DIR / f"{PMID}_supplements_api"
            print(f"Attempting BioC/JATS supplement harvest to: {supp_output_dir}")
            supplements = supplement_client.get_all_supplement_text(
                pmcid, supp_output_dir
            )
            extracted_count = sum(1 for s in supplements if s.text)
            print(
                f"  Retrieved {len(supplements)} supplements via API pipeline; "
                f"{extracted_count} with extractable text"
            )
            for supp in supplements:
                label = supp.label or f"Supplement {supp.index}"
                print(
                    f"  - {label}: {supp.filename or supp.href} (source={supp.source})"
                )
        except Exception as exc:
            print(f"\n  Note: Could not resolve PMCID for PMID {PMID}: {exc}")
            print("  (This is expected if the paper is not in PubMed Central)")

        if not success:
            print(f"\n  FAILED: {result}")
            return 1

    print(f"\n  SUCCESS: Created {result}")

    # Read content
    with open(result, "r") as f:
        content = f.read()

    # Ensure condensed {PMID}_DATA_ZONES.md is created (even when using cached full text)
    if harvester.gene_symbol and not DATA_ZONES_PATH.exists():
        print(f"\nRunning data scout to create condensed zones for {PMID}...")
        harvester._run_data_scout(PMID, content)
        if DATA_ZONES_PATH.exists():
            print(f"  Created {PMID}_DATA_ZONES.md: {DATA_ZONES_PATH}")
        else:
            print(f"  Data scout did not create {PMID}_DATA_ZONES.md")

    print(f"\n{PMID}_FULL_CONTEXT.md stats:")
    print(f"  Content length: {len(content):,} characters")
    print(f"  Line count: {len(content.splitlines()):,} lines")

    # Show {PMID}_DATA_ZONES.md stats if it exists
    if DATA_ZONES_PATH.exists():
        zones_content = DATA_ZONES_PATH.read_text()
        compression = (len(zones_content) / len(content) * 100) if content else 0
        print(f"\n{PMID}_DATA_ZONES.md stats (used for extraction):")
        print(f"  Content length: {len(zones_content):,} characters")
        print(f"  Line count: {len(zones_content.splitlines()):,} lines")
        print(f"  Compression: {compression:.1f}% of original")

    # Check for supplements
    if "## Supplementary" in content or "supplement" in content.lower():
        print("  Contains supplementary material references")

    # Show section headers
    print("\nSection headers found:")
    for line in content.splitlines():
        if line.startswith("#"):
            print(f"  {line}")

    # Check for supplement files
    supp_dir = OUTPUT_DIR / f"{PMID}_supplements"
    if supp_dir.exists():
        supp_files = list(supp_dir.iterdir())
        print(f"\n  Supplement directory exists with {len(supp_files)} files:")
        for f in supp_files:
            print(f"    - {f.name}")

    if not run_extraction:
        print("\n  Skipping triage and extraction (no API key).")
        return 0

    # Extract title and abstract from content
    title, abstract = extract_title_and_abstract_from_content(content)

    if not title:
        title = EXPECTED_TITLE
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
        title=title, abstract=abstract, gene=GENE, pmid=PMID
    )

    print(f"\nTriage Result:")
    print(f"  Decision: {triage_result['decision']}")
    print(f"  Confidence: {triage_result['confidence']:.2f}")
    print(f"  Reason: {triage_result['reason']}")

    if triage_result["decision"] == "DROP":
        print("\n  Paper dropped by triage - no original clinical data detected.")
        print(
            "    (This may be a review, animal study, or purely methodological paper)"
        )
        return 1

    print("\n  Paper passed triage - proceeding to extraction...")

    # ========================================
    # STEP 3: Extract patient/variant information
    # ========================================
    print("\n" + "=" * 60)
    print("STEP 3: EXTRACTING PATIENT AND VARIANT INFORMATION")
    print("=" * 60)

    # Create Paper object for extraction
    paper = Paper(
        pmid=PMID, title=title, abstract=abstract, full_text=content, gene_symbol=GENE
    )

    extractor = ExpertExtractor(fulltext_dir=str(OUTPUT_DIR))
    extraction_result = extractor.extract(paper)

    if not extraction_result.success:
        print(f"\n  Extraction failed: {extraction_result.error}")
        return 1

    print(f"\n  Extraction successful!")
    print(f"  Model used: {extraction_result.model_used}")

    # Display extraction results
    data = extraction_result.extracted_data

    print("\n" + "-" * 50)
    print("EXTRACTION SUMMARY:")
    print("-" * 50)

    if "paper_metadata" in data:
        meta = data["paper_metadata"]
        print(f"\nPaper: {meta.get('title', 'N/A')[:60]}...")
        print(f"Summary: {meta.get('extraction_summary', 'N/A')}")

    if "extraction_metadata" in data:
        emeta = data["extraction_metadata"]
        print(f"\nTotal variants found: {emeta.get('total_variants_found', 0)}")
        print(f"Confidence: {emeta.get('extraction_confidence', 'N/A')}")
        if emeta.get("challenges"):
            print(f"Challenges: {', '.join(emeta['challenges'])}")

    # Show variants
    variants = data.get("variants", [])
    print(f"\n" + "-" * 50)
    print(f"VARIANTS EXTRACTED: {len(variants)}")
    print("-" * 50)

    for i, variant in enumerate(variants[:10], 1):  # Show first 10
        print(f"\n  Variant {i}:")
        print(f"    Gene: {variant.get('gene_symbol', 'N/A')}")
        print(f"    cDNA: {variant.get('cdna_notation', 'N/A')}")
        print(f"    Protein: {variant.get('protein_notation', 'N/A')}")
        print(f"    Significance: {variant.get('clinical_significance', 'N/A')}")

        # Patient count
        patients = variant.get("patients", {})
        if patients and patients.get("count"):
            print(
                f"    Patients: {patients.get('count')} ({patients.get('phenotype', 'N/A')})"
            )

        # Penetrance data
        pdata = variant.get("penetrance_data", {})
        if pdata:
            print(f"    Penetrance:")
            print(
                f"      Total carriers: {pdata.get('total_carriers_observed', 'N/A')}"
            )
            print(f"      Affected: {pdata.get('affected_count', 'N/A')}")
            print(f"      Unaffected: {pdata.get('unaffected_count', 'N/A')}")

        # Individual records
        individuals = variant.get("individual_records", [])
        if individuals:
            print(f"    Individual records: {len(individuals)}")
            for ind in individuals[:3]:  # Show first 3
                print(
                    f"      - {ind.get('individual_id', '?')}: "
                    f"{ind.get('affected_status', '?')}, "
                    f"age={ind.get('age_at_evaluation', '?')}, "
                    f"sex={ind.get('sex', '?')}"
                )
            if len(individuals) > 3:
                print(f"      ... and {len(individuals) - 3} more")

    if len(variants) > 10:
        print(f"\n  ... and {len(variants) - 10} more variants")

    # Save full extraction to JSON
    output_json = OUTPUT_DIR / f"{PMID}_extraction.json"
    with open(output_json, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\n  Full extraction saved to: {output_json}")

    # Summary stats
    print("\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"  PMID: {PMID}")
    print(f"  Gene: {GENE}")
    print(
        f"  Triage: {triage_result['decision']} (confidence: {triage_result['confidence']:.2f})"
    )
    print(f"  Variants extracted: {len(variants)}")

    total_individuals = sum(len(v.get("individual_records", [])) for v in variants)
    print(f"  Individual patient records: {total_individuals}")

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
