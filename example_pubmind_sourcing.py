"""Example: build a PubMind-derived PMID list for harvesting or extraction."""

from pipeline.sourcing import PaperSourcer
# PubMindSourcer - TODO: Update if this class still exists


def run_pubmind_only(query: str):
    print(f"Searching PubMind for query: {query}")
    client = PubMindSourcer()
    pmids = client.search(query)
    print(f"Found {len(pmids)} PMIDs from PubMind")
    print(pmids[:20])


def run_combined_sources(query: str):
    print(f"Aggregating PMIDs for query: {query}")
    sourcer = PaperSourcer()
    pmids = sourcer.fetch_papers(
        gene_symbol=query,
        use_pubmind=True,
        pubmind_query=query,
        max_results_per_source=200,
    )
    print(f"Total deduplicated PMIDs across PubMed, EuropePMC, PubMind: {len(pmids)}")
    print(pmids[:20])


if __name__ == "__main__":
    target = "BRCA2"
    run_pubmind_only(target)
    print("\n" + "=" * 60 + "\n")
    run_combined_sources(target)
