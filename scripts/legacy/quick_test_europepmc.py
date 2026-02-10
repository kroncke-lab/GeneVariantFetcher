#!/usr/bin/env python3
"""
Quick Europe PMC Test

Simple script to quickly test Europe PMC integration.
"""

import sys
from pathlib import Path
import json

# Add project root
sys.path.insert(0, str(Path(__file__).parent))

from gene_literature.europepmc_handler import EuropePMCClient

def quick_test():
    """Run a quick functionality test."""
    print("🧪 Running quick Europe PMC test...")
    
    client = EuropePMCClient()
    
    # Test basic connectivity
    print("1. Testing API connectivity...")
    search_results = client.search_papers('KCNH2', max_results=3)
    
    if search_results:
        print(f"   ✅ Found {len(search_results)} results for 'KCNH2'")
        
        # Test metadata retrieval
        pmid = search_results[0].get('pmid')
        if pmid:
            print(f"2. Testing metadata for PMID {pmid}...")
            metadata = client.get_paper_metadata(str(pmid))
            
            if metadata:
                print(f"   ✅ Retrieved metadata for '{metadata.get('title', 'N/A')[:50]}...'")
                print(f"   📄 PMC ID: {metadata.get('pmcid', 'Not found')}")
                print(f"   🔓 Open Access: {'Yes' if metadata.get('is_open_access') else 'No'}")
            else:
                print("   ❌ Failed to retrieve metadata")
        else:
            print("   ❌ No PMID found in search results")
    else:
        print("   ❌ Search failed - no results found")
        return False
    
    # Test individual volunteer-test PMIDs
    test_pmids = ['30036649', '18808722', '9544837']
    print(f"\n3. Testing golden test PMIDs...")
    
    for pmid in test_pmids:
        metadata = client.get_paper_metadata(pmid)
        if metadata and metadata.get('pmcid'):
            print(f"   ✅ PMID {pmid}: Available (PMC ID: {metadata['pmcid']})")
        else:
            reason = metadata.get('reason', 'No PMC ID') if metadata else 'Not found'
            print(f"   ⚠️  PMID {pmid}: {reason}")
    
    print("\n🎉 Quick test complete!")
    return True

if __name__ == "__main__":
    try:
        success = quick_test()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"❌ Test failed: {e}")
        sys.exit(1)