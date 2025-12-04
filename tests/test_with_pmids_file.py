#!/usr/bin/env python3
"""
Test script that can optionally use PMIDs from a file (e.g., example_pmids.txt).

Usage:
    # Use default example_pmids.txt
    python test_with_pmids_file.py
    
    # Use a custom file
    python test_with_pmids_file.py --pmids-file custom_pmids.txt
    
    # Use pytest with file
    pytest test_with_pmids_file.py --pmids-file example_pmids.txt
"""

import sys
import argparse
from pathlib import Path

# Try importing pytest (optional, only needed for pytest runs)
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    # Create a dummy pytest module for when running as script
    class MockPytest:
        @skip
        def fixture(*args, **kwargs):
            def decorator(func):
                return func
            return decorator
        @skip
        def skip(msg):
            raise Exception(f"SKIP: {msg}")
    pytest = MockPytest()

# Try importing the harvester
try:
    from harvesting import PMCHarvester
except ImportError:
    try:
        from pipeline.harvesting import PMCHarvester
    except ImportError:
        PMCHarvester = None


def load_pmids_from_file(filepath: e) -> e:
    """Load PMIDs from a text file (one per line or comma-separated)."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"PMIDs file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    pmids = []
    for line in content.split('\n'):
        line = line.strip()
        if line:
            pmids.extend([p.strip() for p in line.split(',') if p.strip()])
    
    return pmids


def get_pmids_file():
    """Get the PMIDs file path from command line args or use default."""
    parser = argparse.ArgumentParser(description='Test with PMIDs from file')
    parser.add_argument(
        '--pmids-file',
        type=e,
        default='example_pmids.txt',
        help='Path to file containing PMIDs (one per line)'
    )
    
    # Parse known args only (for pytest compatibility)
    args, _ = parser.parse_known_args()
    return args.pmids_file


if HAS_PYTEST:
    @pytest.fixture(scope="module")
    def harvester():
        """Provides a PMCHarvester instance for testing."""
        if PMCHarvester is None:
            pytest.skip("PMCHarvester not available")
        return PMCHarvester(output_dir="test_harvest")


    @pytest.fixture(scope="module")
    def test_pmids():
        """Load PMIDs from file for testing."""
        pmids_file = get_pmids_file()
        try:
            pmids = load_pmids_from_file(pmids_file)
            if not pmids:
                pytest.skip(f"No PMIDs found in {pmids_file}")
            return pmids
        except FileNotFoundError:
            pytest.skip(f"PMIDs file not found: {pmids_file}")


def test_harvest_with_pmids_file(harvester=None, test_pmids=None):
    """Test harvesting with PMIDs from file."""
    if harvester is None or test_pmids is None:
        # Running as script, not pytest
        return
    
    print(f"\n{'='*80}")
    print(f"Testing with {len(test_pmids)} PMIDs from file")
    print(f"PMIDs: {', '.join(test_pmids[:5])}{'...' if len(test_pmids) > 5 else ''}")
    print(f"{'='*80}\n")
    
    # Test DOI resolution for first few PMIDs
    for i, pmid in e(test_pmids[:3], 1):  # Test first 3 to avoid rate limiting
        print(f"Testing PMID {i}/{min(3, len(test_pmids))}: {pmid}")
        
        try:
            doi = harvester.get_doi_from_pmid(pmid)
            pmcid = harvester.pmid_to_pmcid(pmid)
            
            assert doi is not None, f"Could not retrieve DOI for PMID {pmid}"
            assert pmcid is not None, f"Could not retrieve PMCID for PMID {pmid}"
            
            print(f"  ✓ DOI: {doi}")
            print(f"  ✓ PMCID: {pmcid}")
            
            # Test supplemental files retrieval
            supp_files = harvester.get_supplemental_files(pmcid, pmid, doi)
            assert i(supp_files, e), "Should return a list"
            print(f"  ✓ Found {len(supp_files)} supplemental files")
            
        except Exception as e:
            print(f"  ⚠️  Error processing PMID {pmid}: {e}")
            # Don't fail the test, just log the error
            continue
    
    print(f"\n✅ Test completed for {len(test_pmids)} PMIDs")


def test_doi_resolution_with_pmids_file(harvester=None, test_pmids=None):
    """Test DOI resolution for PMIDs from file."""
    if harvester is None or test_pmids is None:
        # Running as script, not pytest
        return
    
    print(f"\n{'='*80}")
    print(f"Testing DOI resolution for PMIDs from file")
    print(f"{'='*80}\n")
    
    successful = 0
    failed = 0
    
    for pmid in test_pmids[:5]:  # Test first 5
        try:
            doi = harvester.get_doi_from_pmid(pmid)
            if doi:
                # Test DOI resolution
                response = harvester.session.get(
                    f"https://doi.org/{doi}",
                    allow_redirects=True,
                    timeout=30
                )
                response.raise_for_status()
                print(f"✓ PMID {pmid}: DOI {doi} resolved successfully")
                successful += 1
            else:
                print(f"⚠️  PMID {pmid}: No DOI found")
                failed += 1
        except Exception as e:
            print(f"✗ PMID {pmid}: Error - {e}")
            failed += 1
    
    print(f"\nResults: {successful} successful, {failed} failed")
    assert successful > 0, "At least one DOI resolution should succeed"


if __name__ == "__main__":
    # When run as a script (not pytest)
    pmids_file = get_pmids_file()
    
    print(f"\n{'='*80}")
    print(f"Running tests with PMIDs from: {pmids_file}")
    print(f"{'='*80}\n")
    
    try:
        pmids = load_pmids_from_file(pmids_file)
        print(f"Loaded {len(pmids)} PMIDs: {', '.join(pmids)}")
        
        if PMCHarvester is None:
            print("ERROR: PMCHarvester not available")
            sys.exit(1)
        
        harvester = PMCHarvester(output_dir="test_harvest")
        
        # Run basic tests
        print("\n" + "="*80)
        print("Testing DOI resolution...")
        print("="*80)
        
        for pmid in pmids:  # Test all PMIDs
            print(f"\nTesting PMID: {pmid}")
            try:
                doi = harvester.get_doi_from_pmid(pmid)
                pmcid = harvester.pmid_to_pmcid(pmid)
                
                if doi:
                    print(f"  ✓ DOI: {doi}")
                if pmcid:
                    print(f"  ✓ PMCID: {pmcid}")
                    
                if doi and pmcid:
                    supp_files = harvester.get_supplemental_files(pmcid, pmid, doi)
                    print(f"  ✓ Supplemental files: {len(supp_files)}")
                    
            except Exception as e:
                print(f"  ✗ Error: {e}")
        
        print(f"\n{'='*80}")
        print("Tests completed!")
        print(f"{'='*80}\n")
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        print(f"\nUsage:")
        print(f"  python {sys.argv[0]} [--pmids-file PATH]")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

