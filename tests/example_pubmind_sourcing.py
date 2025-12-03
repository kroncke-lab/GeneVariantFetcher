"""Example: build a PubMind-derived PMID list for harvesting or extraction."""

from pipeline.sourcing import PaperSourcer


def run_pubmind_only(query: str):
    """Demonstrates sourcing from only PubMind."""
    print(f"Searching PubMind for query: {query}")
    # Initialize the sourcer
    sourcer = PaperSourcer()
    # Fetch papers using only PubMind by disabling other sources
    pmids = sourcer.fetch_papers(
        gene_symbol=query,
        use_pubmed=False,
        use_europepmc=False,
        use_pubmind=True,
        max_results_per_source=200,
    )
    print(f"Found {len(pmids)} PMIDs from PubMind")
    print(pmids[:20])


def run_combined_sources(query: str):
    """Demonstrates sourcing from all available providers."""
    print(f"Aggregating PMIDs for query: {query}")
    sourcer = PaperSourcer()
    pmids = sourcer.fetch_papers(
        gene_symbol=query,
        use_pubmind=True,
        use_pubmed=True,
        use_europepmc=True,
        max_results_per_source=200,
    )
    print(f"Total deduplicated PMIDs across all sources: {len(pmids)}")
    print(pmids[:20])


if __name__ == "__main__":
    target = "BRCA2"
    run_pubmind_only(target)
    print("\n" + "=" * 60 + "\n")
    run_combined_sources(target)
