#!/usr/bin/env python3
"""
Test abstract extraction with implicit count recognition.

These tests verify that the updated LLM prompts correctly recognize
implicit carrier counts from natural language in abstracts.

Run with: python -m pytest tests/test_abstract_extraction.py -v
Or for live API tests: python tests/test_abstract_extraction.py
"""

import pytest
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Sample abstracts for testing implicit count recognition
SAMPLE_ABSTRACTS = {
    # PMID 10790218 - the case that prompted this fix
    "10790218": {
        "title": "A novel nonsense mutation in HERG causes long QT syndrome",
        "abstract": """We identified a novel nonsense mutation (Y667X) in HERG
        (KCNH2) in a patient with long QT syndrome. The single base pair
        substitution (C to A) at nucleotide position 2001 results in a
        premature stop codon. Interestingly, a healthy individual in the
        same family also carried this mutation, demonstrating incomplete
        penetrance of HERG mutations.""",
        "gene": "KCNH2",
        "expected_variants": ["Y667X"],
        "expected_affected": 1,  # "a patient with long QT syndrome"
        "expected_unaffected": 1,  # "a healthy individual"
        "expected_total": 2,
    },
    # Case report with single affected patient
    "case_report_single": {
        "title": "Novel KCNH2 mutation in a patient with LQTS",
        "abstract": """We report a case of a 35-year-old female patient
        presenting with recurrent syncope and prolonged QT interval.
        Genetic testing revealed a novel heterozygous mutation p.Arg534Cys
        in KCNH2. The proband was successfully treated with beta-blockers.""",
        "gene": "KCNH2",
        "expected_variants": ["p.Arg534Cys", "R534C"],
        "expected_affected": 1,  # proband = 1 affected
        "expected_unaffected": 0,
        "expected_total": 1,
    },
    # Family study with multiple carriers
    "family_study": {
        "title": "Incomplete penetrance of KCNH2 G572R mutation",
        "abstract": """In this family study, we identified the G572R mutation
        in KCNH2 in four family members. Two individuals presented with
        long QT syndrome and documented cardiac events. The other two
        carriers remained asymptomatic despite carrying the same mutation,
        suggesting incomplete penetrance.""",
        "gene": "KCNH2",
        "expected_variants": ["G572R"],
        "expected_affected": 2,  # "Two individuals presented with LQTS"
        "expected_unaffected": 2,  # "other two carriers remained asymptomatic"
        "expected_total": 4,
    },
}


class TestPromptContent:
    """Test that the prompts contain the necessary guidance."""

    def test_abstract_prompt_has_implicit_count_guidance(self):
        """Verify AbstractCarrierExtractor prompt includes implicit count rules."""
        from pipeline.abstract_extraction import ABSTRACT_CARRIER_EXTRACTION_PROMPT

        # Check for key phrases that guide implicit count recognition
        assert "a patient" in ABSTRACT_CARRIER_EXTRACTION_PROMPT.lower()
        assert "a healthy individual" in ABSTRACT_CARRIER_EXTRACTION_PROMPT.lower()
        assert "1 carrier" in ABSTRACT_CARRIER_EXTRACTION_PROMPT.lower()
        assert "implicit" in ABSTRACT_CARRIER_EXTRACTION_PROMPT.lower()

    def test_extraction_prompt_has_implicit_count_guidance(self):
        """Verify EXTRACTION_PROMPT includes implicit count rules."""
        from pipeline.prompts import EXTRACTION_PROMPT

        assert "recognizing implicit counts" in EXTRACTION_PROMPT.lower()
        assert "a patient" in EXTRACTION_PROMPT.lower()
        assert "1 carrier" in EXTRACTION_PROMPT.lower()

    def test_extraction_prompt_has_abstract_only_section(self):
        """Verify EXTRACTION_PROMPT has abstract-only extraction guidance."""
        from pipeline.prompts import EXTRACTION_PROMPT

        assert "abstract-only extraction" in EXTRACTION_PROMPT.lower()
        assert "do not skip extraction" in EXTRACTION_PROMPT.lower()


class TestAbstractCarrierExtractor:
    """Test AbstractCarrierExtractor with sample abstracts."""

    @pytest.fixture
    def extractor(self):
        """Create an AbstractCarrierExtractor instance."""
        # Skip if no API key
        if not os.getenv("OPENAI_API_KEY"):
            pytest.skip("No OpenAI API key available")

        from pipeline.abstract_extraction import AbstractCarrierExtractor
        return AbstractCarrierExtractor()

    @pytest.mark.parametrize("pmid", SAMPLE_ABSTRACTS.keys())
    def test_implicit_count_extraction(self, extractor, pmid):
        """Test that implicit counts are correctly extracted from abstracts."""
        from utils.models import Paper

        sample = SAMPLE_ABSTRACTS[pmid]
        paper = Paper(
            pmid=pmid,
            title=sample["title"],
            abstract=sample["abstract"],
            gene_symbol=sample["gene"],
        )

        result = extractor.extract(paper, sample["gene"])

        # Verify extraction succeeded
        assert result.success, f"Extraction failed for {pmid}: {result.error}"

        # Verify carrier counts were extracted
        assert result.has_carrier_counts, (
            f"No carrier counts extracted for {pmid}. "
            f"Expected affected={sample['expected_affected']}, "
            f"unaffected={sample['expected_unaffected']}. "
            f"Notes: {result.extraction_notes}"
        )

        # Verify affected count (allow some flexibility)
        if sample["expected_affected"] > 0:
            assert result.affected_count is not None and result.affected_count > 0, (
                f"Expected affected_count > 0 for {pmid}, got {result.affected_count}"
            )

        # Verify unaffected count
        if sample["expected_unaffected"] > 0:
            assert result.unaffected_count is not None and result.unaffected_count > 0, (
                f"Expected unaffected_count > 0 for {pmid}, got {result.unaffected_count}"
            )

        # Verify variants were found
        if sample["expected_variants"]:
            assert result.variants_mentioned, (
                f"No variants extracted for {pmid}, expected {sample['expected_variants']}"
            )

        print(f"\n✓ {pmid}: affected={result.affected_count}, "
              f"unaffected={result.unaffected_count}, "
              f"variants={result.variants_mentioned}")


def run_live_tests():
    """Run live API tests (for manual testing)."""
    print("=" * 60)
    print("LIVE API TESTS - Abstract Extraction with Implicit Counts")
    print("=" * 60)

    # Check for API key
    if not os.getenv("OPENAI_API_KEY"):
        print("\n❌ No OpenAI API key found. Set OPENAI_API_KEY in .env to run tests.")
        return

    from pipeline.abstract_extraction import AbstractCarrierExtractor
    from utils.models import Paper

    extractor = AbstractCarrierExtractor()

    for pmid, sample in SAMPLE_ABSTRACTS.items():
        print(f"\n{'='*40}")
        print(f"Testing PMID: {pmid}")
        print(f"Title: {sample['title']}")
        print(f"Expected: affected={sample['expected_affected']}, unaffected={sample['expected_unaffected']}")
        print(f"{'='*40}")

        paper = Paper(
            pmid=pmid,
            title=sample["title"],
            abstract=sample["abstract"],
            gene_symbol=sample["gene"],
        )

        result = extractor.extract(paper, sample["gene"])

        if result.success:
            print(f"✓ Extraction successful")
            print(f"  - Total carriers: {result.total_carriers}")
            print(f"  - Affected: {result.affected_count}")
            print(f"  - Unaffected: {result.unaffected_count}")
            print(f"  - Variants: {result.variants_mentioned}")
            print(f"  - Confidence: {result.confidence:.2f}")
            print(f"  - Notes: {result.extraction_notes}")

            # Check if extraction matches expectations
            if result.affected_count == sample["expected_affected"] and \
               result.unaffected_count == sample["expected_unaffected"]:
                print("  ✅ PASS - Counts match expectations!")
            else:
                print("  ⚠️ PARTIAL - Counts differ from expectations")
        else:
            print(f"❌ Extraction failed: {result.error}")

    print("\n" + "=" * 60)
    print("Tests complete!")
    print("=" * 60)


if __name__ == "__main__":
    run_live_tests()
