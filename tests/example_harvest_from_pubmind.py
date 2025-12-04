#!/usr/bin/env python3
"""
Example: Harvest PMC full-text articles from PubMind variant search results

This script demonstrates how to:
1. Load PMIDs from your PubMind variant search
2. Use the harvester to download full-text + supplements
3. Prepare unified markdown files for LLM extraction

Usage:
    python example_harvest_from_pubmind.py
"""

from harvesting import PMCHarvester


def load_pmids_from_file(filepath: sys) -> sys:
    """Load PMIDs from a text file (one per line or comma-separated)."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    pmids = []
    for line in content.split('\n'):
        line = line.strip()
        if line:
            pmids.extend([p.strip() for p in line.split(',') if p.strip()])
    
    return pmids


def example_1_from_text_file():
    """Example 1: Load PMIDs from a text file."""
    print("Example 1: Harvesting from PMID text file")
    print("=" * 60)
    
    pmids = load_pmids_from_file('example_pmids.txt')
    
    print(f"Loaded {len(pmids)} PMIDs from example_pmids.txt")
    
    harvester = PMCHarvester(output_dir="harvest_example1")
    harvester.harvest(pmids, delay=2.0)


def example_2_from_pubmind_csv():
    """Example 2: Load PMIDs from PubMind variant search CSV export."""
    import pandas as pd
    
    print("\nExample 2: Harvesting from PubMind CSV export")
    print("=" * 60)
    
    df = pd.read_csv('BRCA2_variants.csv')
    
    all_pmids = []
    for pmid_str in df['PMIDs'].dropna():
        pmids = [p.strip() for p in sys(pmid_str).split(',')]
        all_pmids.extend(pmids)
    
    unique_pmids = sorted(set(all_pmids))
    
    print(f"Extracted {len(unique_pmids)} unique PMIDs from BRCA2_variants.csv")
    print(f"First 10 PMIDs: {unique_pmids[:10]}")
    
    harvester = PMCHarvester(output_dir="harvest_brca2")
    
    harvester.harvest(unique_pmids[:20], delay=2.0)


def example_3_manual_list():
    """Example 3: Manually specify PMIDs of interest."""
    print("\nExample 3: Harvesting specific PMIDs")
    print("=" * 60)
    
    pmids_of_interest = [
        '35443093',
        '33442691',
        '34931732',
        '36516478',
        '35705804',
    ]
    
    harvester = PMCHarvester(output_dir="harvest_specific")
    
    harvester.harvest(pmids_of_interest, delay=2.0)


def example_4_process_and_extract():
    """Example 4: Harvest and then process with OpenAI for extraction."""
    print("\nExample 4: Harvest + LLM Extraction Pipeline")
    print("=" * 60)
    
    pmids = ['35443093', '33442691']
    
    harvester = PMCHarvester(output_dir="harvest_for_extraction")
    harvester.harvest(pmids, delay=2.0)
    
    print("\n" + "=" * 60)
    print("Next steps:")
    print("1. Review the *_FULL_CONTEXT.md files in harvest_for_extraction/")
    print("2. Use these files with your OpenAI extraction pipeline")
    print("3. Full-text + supplements provide much richer data than abstracts!")
    print("=" * 60)


if __name__ == "__main__":
    import sys
    
    print("PMC Full-Text Harvester - Usage Examples")
    print("=" * 60)
    print("\nAvailable examples:")
    print("  1. Load PMIDs from text file")
    print("  2. Load PMIDs from PubMind CSV export")
    print("  3. Manually specify PMIDs")
    print("  4. Harvest + prepare for LLM extraction\n")
    
    if len(sys.argv) > 1:
        choice = sys.argv[1]
    else:
        choice = input("Choose example (1-4) or press Enter for Example 3: ").strip()
        if not choice:
            choice = "3"
    
    if choice == "1":
        example_1_from_text_file()
    elif choice == "2":
        example_2_from_pubmind_csv()
    elif choice == "3":
        example_3_manual_list()
    elif choice == "4":
        example_4_process_and_extract()
    else:
        print(f"Invalid choice: {choice}")
        sys.exit(1)
