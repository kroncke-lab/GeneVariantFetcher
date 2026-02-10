#!/usr/bin/env python3
"""
Standalone Europe PMC Test

Tests Europe PMC API directly without GVF imports.
"""

import json
import sys
import requests
from pathlib import Path

def test_direct_api():
    """Test Europe PMC API directly."""
    print("🔍 Testing Europe PMC API directly...")
    
    # Set up basic session
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'EuropePMCTest/1.0',
        'Accept': 'application/json'
    })
    
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    
    try:
        # Test basic search
        print("1. Testing basic search...")
        search_url = f"{base_url}/search"
        params = {
            'query': 'KCNH2[Gene Symbol]',
            'format': 'json',
            'pageSize': 3,
            'resultType': 'core'
        }
        
        response = session.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        results = data.get('resultList', {}).get('result', [])
        
        if results:
            print(f"   ✅ Found {len(results)} results for 'KCNH2'")
            
            # Get first PMID
            first_pmid = str(results[0].get('pmid', ''))
            print(f"   First article: PMID {first_pmid}")
            
            # Test metadata
            print("2. Testing metadata retrieval...")
            metadata_url = f"{base_url}/search"
            params = {
                'query': f'EXT_ID:{first_pmid}',
                'format': 'json',
                'resultType': 'core'
            }
            
            response = session.get(metadata_url, params=params, timeout=30)
            response.raise_for_status()
            
            metadata_data = response.json()
            metadata = metadata_data.get('resultList', {}).get('result', [{}])[0]
            
            pmcid = metadata.get('pmcid')
            title = metadata.get('title', '')[:50]
            
            print(f"   📄 Title: {title}...")
            print(f"   📄 PMC ID: {pmcid}")
            
            # Test full-text XML if PMC ID available
            if pmcid:
                print("3. Testing full-text XML...")
                xml_url = f"{base_url}/PMC{pmcid}/fullTextXML"
                try:
                    xml_response = session.get(xml_url, timeout=30)
                    if xml_response.status_code == 200 and xml_response.text.strip():
                        xml_size = len(xml_response.text)
                        print(f"   ✅ Full-text XML available ({xml_size} characters)")
                except Exception as e:
                    print(f"   ⚠️  Full-text XML error: {e}")
            
            # Test specific golden test PMIDs
            print("4. Testing golden test PMIDs...")
            golden_pmids = ['30036649', '18808722', '9544837']
            
            for pmid in golden_pmids:
                params = {'query': f'EXT_ID:{pmid}', 'format': 'json'}
                response = session.get(search_url, params=params, timeout=30)
                response.raise_for_status()
                
                result_data = response.json()
                result = result_data.get('resultList', {}).get('result', [])
                
                if result:
                    pmcid_result = result[0].get('pmcid')
                    availability = "Available" if pmcid_result else "Abstract only"
                    print(f"   PMID {pmid}: {availability}")
                else:
                    print(f"   PMID {pmid}: Not found")
            
            print("\n✅ All tests passed!")
            return True
            
        else:
            print("   ❌ No search results found")
            return False
            
    except Exception as e:
        print(f"❌ Test failed: {e}")
        return False

def create_api_response_file():
    """Create a simple response file documenting API functionality."""
    
    test_data = {
        "test_run": str(Path(__file__).stem),
        "test_date": "2026-02-07",
        "base_url": "https://www.ebi.ac.uk/europepmc/webservices/rest",
        "test_endpoints": [
            {
                "endpoint": "/search",
                "description": "Search papers",
                "parameters": {
                    "query": "free text query",
                    "format": "json/xml",
                    "pageSize": "1-1000",
                    "resultType": "idlist/lite/core"
                },
                "examples": [
                    {"query": "KCNH2[Gene Symbol]", "expected": "Papers about KCNH2"},
                    {"query": "KCNH2 AND variant", "expected": "Papers with variants"}
                ]
            },
            {
                "endpoint": "/{pmcid}/fullTextXML",
                "description": "Get full-text XML for PMC article",
                "parameters": {
                    "pmcid": "PMC1234567 (without PMC prefix)"
                }
            },
            {
                "endpoint": "/{pmcid}/supplementaryFiles",
                "description": "Get supplementary materials",
                "parameters": {
                    "pmcid": "PMC1234567"
                }
            }
        ],
        "coverage": {
            "total_articles": "33M+",
            "open_access": "6.5M+",
            "full_text": "10.2M+",
            "supplementary_files": "Available for open access articles"
        }
    }
    
    # Save to file
    output_path = Path('/mnt/temp2/kronckbm/gvf_output/EUROPEPMC_API_SPEC.json')
    with open(output_path, 'w') as f:
        json.dump(test_data, f, indent=2, ensure_ascii=False)
    
    print(f"📁 API specification saved to: {output_path}")
    return test_data

if __name__ == "__main__":
    print("🧪 Standalone Europe PMC API Test")
    print("=" * 40)
    
    # Test basic functionality
    success = test_direct_api()
    
    # Create documentation
    api_spec = create_api_response_file()
    
    sys.exit(0 if success else 1)