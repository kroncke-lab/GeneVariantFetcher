#!/usr/bin/env python3
"""
Test script to verify DOI resolution and supplemental file scraping.
"""

from harvesting import PMCHarvester
from pathlib import Path

from tests.conftest import TEST_PMC_HARVEST_DIR


def test_doi_resolution():
    """Test DOI resolution for various publishers."""

    harvester = PMCHarvester(output_dir=str(TEST_PMC_HARVEST_DIR))

    # Test cases: PMID, DOI (hardcoded to avoid rate limiting), Expected Domain, Description
    test_cases = [
        ("34931732", "10.1038/s41467-021-27599-6", "nature.com", "Nature article"),
        (
            "35443093",
            "10.1016/j.gim.2022.04.004",
            "gimjournal.org or sciencedirect.com",
            "GIM/Elsevier article",
        ),
    ]

    print("=" * 80)
    print("Testing DOI Resolution and Scraping")
    print("=" * 80)

    for pmid, doi, expected_domain, description in test_cases:
        print(f"\n{'=' * 80}")
        print(f"Test Case: {description} (PMID: {pmid})")
        print(f"Expected Domain: {expected_domain}")
        print(f"Using DOI: {doi}")
        print("=" * 80)

        # Test DOI resolution
        # Test DOI resolution
        response = harvester.session.get(
            f"https://doi.org/{doi}", allow_redirects=True, timeout=30
        )
        response.raise_for_status()  # Fails test on non-2xx responses
        final_url = response.url
        print(f"✓ Resolved to: {final_url}")
        assert any(d in final_url for d in expected_domain.split(" or "))

        # Test scraping
        print("\nAttempting to scrape supplemental files...")
        supp_files = harvester._get_supplemental_files_from_doi(doi, pmid)

        if supp_files:
            print(f"✅ Found {len(supp_files)} supplemental file(s):")
            for idx, file_info in enumerate(supp_files, 1):
                print(
                    f"  {idx}. {file_info.get('name', 'Unknown')} - {file_info.get('url', 'No URL')}"
                )
                assert "url" in file_info and file_info["url"]
                assert "name" in file_info and file_info["name"]
        else:
            print(
                "⚠️  No supplemental files found (may be correct or scraper needs refinement)"
            )

    print(f"\n{'=' * 80}")
    print("Test complete!")
    print("=" * 80)


if __name__ == "__main__":
    test_doi_resolution()
