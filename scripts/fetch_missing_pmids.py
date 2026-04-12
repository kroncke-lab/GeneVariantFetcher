import subprocess
import time
import sys

def main():
    missing_file = "/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt"
    try:
        with open(missing_file, 'r') as f:
            pmids = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"File not found: {missing_file}")
        return
        
    print(f"Found {len(pmids)} missing PMIDs to fetch.")
    
    # Just do the first 5 as a test
    test_pmids = pmids[:5]
    for pmid in test_pmids:
        print(f"\n[{time.strftime('%H:%M:%S')}] Fetching PMID {pmid}...")
        cmd = ["python", "-m", "cli", "run", "--gene", "KCNH2", "--pmids", pmid, "--skip-llm"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd="/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher")
        
        if result.returncode == 0:
            print(f"✅ Success: PMID {pmid}")
        else:
            print(f"❌ Failed: PMID {pmid}")
            print(f"Error: {result.stderr[:200]}")
            
        time.sleep(1)

if __name__ == "__main__":
    main()
