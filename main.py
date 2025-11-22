"""
GeneVariantFetcher - Main Entry Point

This is a collection of command-line tools for extracting genetic variant data.
For usage instructions, see README.md or run the example scripts:

- python example_automated_workflow.py GENE --email your@email.com
- python example_usage.py
- python example_triage.py
- python harvest_pmc_fulltext.py
"""


def main():
    print("=" * 80)
    print("GeneVariantFetcher - Genetic Variant Data Extraction Tools")
    print("=" * 80)
    print("\nThis package provides several tools for extracting variant data from PubMed.")
    print("\nQuick Start:")
    print("  1. Automated workflow: python example_automated_workflow.py GENE --email your@email.com")
    print("  2. Extraction pipeline: python example_usage.py")
    print("  3. Clinical triage: python example_triage.py")
    print("  4. PMC harvester: python harvest_pmc_fulltext.py")
    print("\nFor detailed documentation, see README.md")
    print("=" * 80)


if __name__ == "__main__":
    main()
