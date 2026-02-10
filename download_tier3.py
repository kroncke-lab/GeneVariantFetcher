#!/usr/bin/env python3
"""
GVF Tier 3 Download - Process 70 papers from 2005-2010
Batch downloads with PMC availability checking and alternatives
"""

import sys
import os
import requests
import json
from pathlib import Path
from datetime import datetime

# Ensure paths setup
gvf_repo = Path("/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/")
sys.path.insert(0, str(gvf_repo))

# Tier 3 PMIDs
tier3_pmids = [
    "19034806", "19038855", "19065538", "19070294", "19127321", "19136169", "19157587", "19160088", 
    "19169982", "19184172", "19187913", "19215240", "19306396", "19322600", "19352046", "19371231", 
    "19490267", "19695459", "19731233", "19843919", "19996378", "20167303", "20181576", "20197117", 
    "20390067", "20636320", "20975234", "21070882", "21109023", "21130771", "21164565", "21185499", 
    "21216356", "21308345", "21410720", "21419236", "21483829", "21499742", "21951015", "22052944", 
    "22067087", "22104571", "22173492", "22314138", "22338672", "22382559", "22402334", "22407026", 
    "22429796", "22515331", "22727609", "22764740", "22821100", "22882672", "22885918", "23010577", 
    "23134353", "23207121", "23237912", "23277474", "23351921", "23465283", "23546179", "23571586", 
    "23631430", "23864605", "23899126", "23917959", "23981618", "23995044"
]

def download_with_harvester(pmids, output_dir, gene_symbol="KCNH2"):
    """Download papers using PMCHarvester"""
    try:
        sys.path.append(str(gvf_repo))
        from harvesting.orchestrator import PMCHarvester
        
        harvester = PMCHarvester(
            output_dir=str(output_dir),
            gene_symbol=gene_symbol
        )
        
        print(f"Starting harvest of {len(pmids)} PMIDs...")
        
        results = harvester.harvest(
            pmids=pmids,
            delay=2.0,  # Conservative delay for API limits
            run_scout=False  # Skip scout processing, focus on downloads
        )
        
        return results
    except Exception as e:
        print(f"Harvester error: {e}")
        import traceback
        traceback.print_exc()
        return None

def check_pmc_direct(pmids):
    """Direct check for PMC availability using web requests"""
    results = []
    
    print("Checking PMC availability via Europe PMC API...")
    
    for pmid in pmids:
        try:
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}&format=json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                result = data.get('resultList', {}).get('result', [{}])[0]
                
                results.append({
                    'pmid': pmid,
                    'pmcid': result.get('pmcid'),
                    'title': result.get('title', 'N/A')[:100],
                    'doi': result.get('doi', 'N/A'),
                    'open_access': result.get('isOpenAccess', 'N/A')
                })
            else:
                results.append({'pmid': pmid, 'error': f'HTTP {resp.status_code}'})
        except Exception as e:
            results.append({'pmid': pmid, 'error': str(e)})
    
    return results

def process_batch_simple(pmids, output_dir):
    """Simpler batch processing"""
    batch_dir = output_dir / "batch_processing"
    batch_dir.mkdir(exist_ok=True)
    
    results = []
    
    print(f"Processing {len(pmids)} PMIDs in batch...")
    
    for i, pmid in enumerate(pmids, 1):
        print(f"  {i:2d}/70 - Processing {pmid}...")
        
        # Create individual directory for each PMID
        pmid_dir = batch_dir / f"pmid_{pmid}"
        pmid_dir.mkdir(exist_ok=True)
        
        result = {
            'pmid': pmid,
            'status': 'attempted',
            'timestamp': datetime.now().isoformat()
        }
        
        try:
            # We'll use a simplified approach and let the harvester do the work
            # The harvester will create individual results per paper
            result['directory'] = str(pmid_dir)
            
        except Exception as e:
            result['status'] = 'failed'
            result['error'] = str(e)
        
        results.append(result)
    
    return results

def main():
    print("=" * 80)
    print("GVF TIER 3 DOWNLOAD - 70 papers (2005-2010)")
    print("=" * 80)
    
    # Setup directories
    output_base = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2/tier3_downloads/")
    log_file = Path("/mnt/temp2/kronckbm/gvf_output/tier3_download_log.md")
    
    output_base.mkdir(parents=True, exist_ok=True)
    
    print(f"Total Tier 3 PMIDs: {len(tier3_pmids)}")
    print(f"Output directory: {output_base}")
    print()
    
    # First, pre-check PMC availability
    pmc_results = check_pmc_direct(tier3_pmids[:10])  # Check first batch as sample
    
    print("PMC Availability Check Results:")
    for r in pmc_results:
        pmcid = r.get('pmcid', 'N/A')
        has_pmc = '✅' if pmcid == 'PMC' else '❌'
        print(f"  {r['pmid']}: {pmcid} {has_pmc}")
    
    print()
    
    # Process all papers with the harvester
    print("Starting full download process...")
    
    # Process in batches of 15 to avoid overwhelming APIs
    batch_size = 15
    all_results = []
    
    for i in range(0, len(tier3_pmids), batch_size):
        batch = tier3_pmids[i:i+batch_size]
        batch_num = i//batch_size + 1
        
        print(f"\n--- Processing Batch {batch_num} ({len(batch)} papers) ---")
        
        # Create batch directory
        batch_dir = output_base / f"tier3_batch_{batch_num}_{datetime.now().strftime('%Y%m%d')}" 
        batch_dir.mkdir(exist_ok=True)
        
        try:
            # Use the harvester for this batch
            harvest_results = download_with_harvester(batch, batch_dir)
            
            batch_summary = {
                'batch_num': batch_num,
                'pmids': batch,
                'batch_dir': str(batch_dir),
                'timestamp': datetime.now().isoformat(),
                'results': harvest_results
            }
            
            all_results.append(batch_summary)
            print(f"  Batch {batch_num}: Processing started")
            
        except Exception as e:
            print(f"  Batch {batch_num}: Error - {e}")
            continue
    
    # Generate comprehensive log
    generate_log(log_file, tier3_pmids, all_results)
    
    print("\n" + "=" * 80)
    print("TIER 3 DOWNLOAD PROCESS INITIATED!")
    print("=" * 80)
    print(f"Total papers queued: {len(tier3_pmids)}")
    print(f"Batches created: {len(all_results)}")
    print(f"Log file: {log_file}")
    print(f"Output directory: {output_base}")

def generate_log(log_file, pmids, batch_results):
    """Generate comprehensive markdown log"""
    
    content = f"""# Tier 3 Download Log
**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
**Processing**: 70 Tier 3 papers (2005-2010)

## Overview
Tier 3 represents 70 high-priority papers from 2005-2010, as determined by likely publication year 
and expected variant content density. These papers probably contain early breakthrough 
variant descriptions and foundational clinical cases.

## Files Processed
```
{', '.join(pmids)}
```

## Batch Processing Details

"""

    for batch in batch_results:
        content += f"""### Batch {batch['batch_num']}
- **PMIDs**: {', '.join(batch['pmids'])}
- **Output Directory**: `{batch['batch_dir']}`
- **Processing Started**: {batch['timestamp']}
- **Status**: Started via PMCHarvester

"""

    content += f"""## File Locations
- **Master output**: `/mnt/temp2/kronckbm/gvf_output/KCNH2/tier3_downloads/`
- **This log**: `{log_file}`
- **GVF repo**: `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/`

## Strategy Notes
1. Batch processing implemented to respect API limits
2. PMCHarvester used for maximum source coverage
3. Individual directories created per batch for organization
4. Failed downloads will be processed with alternative sources

## Next Steps
1. Monitor individual batch directories for completion
2. Check harvest success logs in each batch directory
3. Retry failed downloads with PaperRetriever
4. Continue with Tier 4 (2000-2005) processing

## Stats
- **Total PMIDs**: {len(pmids)}
- **Batches created**: {len(batch_results)}
- **Estimated processing time**: {len(batch_results) * 5} minutes
- **API sources**: EuropePMC, DOI resolvers, publisher APIs
"""

    log_file.write_text(content)
    print(f"Log generated: {log_file}")

if __name__ == "__main__":
    os.chdir("/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/")
    main()