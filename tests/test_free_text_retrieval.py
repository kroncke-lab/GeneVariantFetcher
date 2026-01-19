#!/usr/bin/env python3
"""
Test script for free full-text retrieval from publisher pages.

Tests the new functionality that fetches full text and supplements from
publisher websites for PubMed articles that are marked as "Free Full Text"
but do not have PMCIDs.

Example PMIDs without PMCID but with free full text:
- 10336646: Free article on Nature Genetics
- 10827225: Free article
"""

import os
import sys
import tempfile

# Ensure NCBI email is set
# NCBI_EMAIL should be set in .env file

from harvesting.orchestrator import PMCHarvester
from harvesting.pmc_api import PMCAPIClient


def test_free_text_detection():
    """Test detection of free full text articles."""
    print("=" * 60)
    print("Testing free full text detection")
    print("=" * 60)

    client = PMCAPIClient()

    # Test PMIDs that the user mentioned as free articles without PMCID
    test_pmids = ["10336646", "10827225"]

    for pmid in test_pmids:
        print(f"\nTesting PMID: {pmid}")

        # Check for PMCID
        pmcid = client.pmid_to_pmcid(pmid)
        print(f"  PMCID: {pmcid or 'None'}")

        # Check for DOI
        doi = client.get_doi_from_pmid(pmid)
        print(f"  DOI: {doi or 'None'}")

        # Check for free full text status
        is_free, free_url = client.is_free_full_text(pmid)
        print(f"  Is Free Full Text: {is_free}")
        print(f"  Free URL: {free_url or 'None'}")

        # Get all LinkOut URLs for debugging
        linkouts = client.get_pubmed_linkout_urls(pmid)
        if linkouts:
            print(f"  LinkOut URLs ({len(linkouts)} found):")
            for link in linkouts[:3]:  # Show first 3
                print(f"    - {link['provider']}: {link['url'][:60]}...")
                if link.get('attributes'):
                    print(f"      Attributes: {link['attributes']}")


def test_full_text_retrieval():
    """Test full text retrieval from publisher pages."""
    print("\n" + "=" * 60)
    print("Testing full text retrieval from publisher")
    print("=" * 60)

    # Use a temporary directory for output
    with tempfile.TemporaryDirectory() as tmpdir:
        harvester = PMCHarvester(output_dir=tmpdir)

        # Test with the example PMIDs
        test_pmids = ["10336646", "10827225"]

        for pmid in test_pmids:
            print(f"\n{'=' * 40}")
            print(f"Processing PMID: {pmid}")
            print("=" * 40)

            success, result = harvester.process_pmid(pmid)

            if success:
                print(f"\n✅ SUCCESS: {result}")
                # Show first 500 characters of output
                with open(result, 'r') as f:
                    content = f.read()
                    print(f"\nOutput preview ({len(content)} total characters):")
                    print("-" * 40)
                    print(content[:500])
                    print("..." if len(content) > 500 else "")
            else:
                print(f"\n❌ FAILED: {result}")


def main():
    """Run all tests."""
    print("=" * 60)
    print("FREE FULL TEXT RETRIEVAL TEST")
    print("=" * 60)

    # Check if required environment variable is set
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        print("\nWARNING: NCBI_EMAIL not set in .env file.")

    # Run tests
    test_free_text_detection()
    test_full_text_retrieval()

    print("\n" + "=" * 60)
    print("TEST COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
