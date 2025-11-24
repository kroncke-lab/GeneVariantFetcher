import pytest
from harvest_pmc_fulltext import PMCHarvester

# A pytest fixture to create a single harvester instance for all tests.
# The 'scope="module"' part means this will only run once for this test file.
@pytest.fixture(scope="module")
def harvester():
    """Provides a PMCHarvester instance for testing."""
    # Use a separate directory for test outputs
    return PMCHarvester(output_dir="test_harvest")

def test_nature_doi_resolution_and_scraping(harvester):
    """
    Tests the full pipeline for a Nature article (PMID: 34931732).
    This article is expected to fail the initial API lookup and fall back 
    to the DOI scraping mechanism, which should then succeed.
    """
    pmid = '34931732'
    doi = harvester.get_doi_from_pmid(pmid)
    pmcid = harvester.pmid_to_pmcid(pmid)

    assert doi is not None, "Test setup failed: Could not retrieve DOI for Nature PMID"
    assert pmcid is not None, "Test setup failed: Could not retrieve PMCID for Nature PMID"

    # We call the top-level function to test the complete "waterfall" logic
    supp_files = harvester.get_supplemental_files(pmcid, pmid, doi)

    # Assertions
    assert isinstance(supp_files, list), "The function should always return a list."
    assert len(supp_files) > 0, "The Nature scraper should have found at least one supplemental file."
    
    # Check the structure of the returned data
    first_file = supp_files[0]
    assert 'url' in first_file, "Each file dictionary must have a 'url' key."
    assert 'name' in first_file, "Each file dictionary must have a 'name' key."
    assert first_file['url'].startswith('http'), "The file URL should be a full, valid URL."
    
    print(f"\n✅ Nature test passed. Found {len(supp_files)} files. First file: {first_file['name']}")

def test_gim_doi_resolution_and_scraping(harvester):
    """
    Tests the full pipeline for a Genetics in Medicine (Elsevier) article (PMID: 35443093).
    This article is also expected to fall back to the DOI scraper.
    """
    pmid = '35443093'
    doi = harvester.get_doi_from_pmid(pmid)
    pmcid = harvester.pmid_to_pmcid(pmid)

    assert doi is not None, "Test setup failed: Could not retrieve DOI for GIM PMID"
    assert pmcid is not None, "Test setup failed: Could not retrieve PMCID for GIM PMID"

    supp_files = harvester.get_supplemental_files(pmcid, pmid, doi)

    assert isinstance(supp_files, list)
    assert len(supp_files) > 0, "The Elsevier/GIM scraper should have found at least one supplemental file."
    
    first_file = supp_files[0]
    assert 'url' in first_file
    assert 'name' in first_file
    assert first_file['url'].startswith('http')

    print(f"✅ GIM/Elsevier test passed. Found {len(supp_files)} files. First file: {first_file['name']}")

# To run these tests:
# 1. Make sure you have pytest installed: pip install pytest
# 2. Run pytest from your terminal in the project root directory: pytest -v
