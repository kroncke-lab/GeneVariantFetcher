#!/usr/bin/env python3
"""
Check publishers for non-PMC PMIDs by looking up their DOIs.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Bio import Entrez
from pathlib import Path

Entrez.email = os.getenv("NCBI_EMAIL", "brett.kroncke@vumc.org")
Entrez.api_key = os.getenv("NCBI_API_KEY")

RATE_LIMIT = 0.11 if Entrez.api_key else 0.34

# Publisher detection by DOI prefix
PUBLISHER_PREFIXES = {
    "10.1002/": "Wiley",
    "10.1111/": "Wiley-Blackwell",
    "10.1016/": "Elsevier",
    "10.1161/": "AHA (Circulation/JAHA)",
    "10.1093/": "Oxford Academic",
    "10.1038/": "Nature",
    "10.1007/": "Springer",
    "10.1056/": "NEJM",
    "10.1001/": "JAMA",
    "10.1136/": "BMJ",
    "10.1371/": "PLOS",
    "10.1172/": "JCI",
    "10.1073/": "PNAS",
    "10.1126/": "Science",
    "10.1053/": "Elsevier (Mosby)",
    "10.1097/": "Lippincott Williams",
    "10.1152/": "APS Journals",
    "10.1034/": "Wiley-Blackwell (old)",
    "10.1046/": "Wiley-Blackwell (old)",
}


def get_doi(pmid: str) -> str | None:
    """Get DOI for a PMID."""
    try:
        time.sleep(RATE_LIMIT)
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        pubmed_article = records["PubmedArticle"][0]
        article = pubmed_article["MedlineCitation"]["Article"]

        # Check ELocationID
        for item in article.get("ELocationID", []):
            if item.attributes.get("EIdType") == "doi":
                return str(item)

        # Check ArticleIdList in PubmedData
        pubmed_data = pubmed_article.get("PubmedData", {})
        for identifier in pubmed_data.get("ArticleIdList", []):
            if identifier.attributes.get("IdType") == "doi":
                return str(identifier)

        return None
    except Exception as e:
        print(f"  Error for {pmid}: {e}", file=sys.stderr)
        return None


def identify_publisher(doi: str) -> str:
    """Identify publisher from DOI prefix."""
    if not doi:
        return "Unknown (no DOI)"
    for prefix, publisher in PUBLISHER_PREFIXES.items():
        if doi.startswith(prefix):
            return publisher
    return f"Unknown ({doi[:15]}...)"


def main():
    # Read non-PMC PMIDs
    input_file = Path(__file__).parent.parent / "tests/fixtures/pmids/non_pmc_pmids.txt"

    pmids = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            pmid = line.split()[0]
            pmids.append(pmid)

    print(f"Looking up DOIs for {len(pmids)} non-PMC PMIDs...\n")

    results = []
    publisher_counts = {}

    for i, pmid in enumerate(pmids, 1):
        doi = get_doi(pmid)
        publisher = identify_publisher(doi)
        print(f"[{i}/{len(pmids)}] {pmid}: {doi or 'NO DOI'} -> {publisher}")
        results.append((pmid, doi, publisher))
        publisher_counts[publisher] = publisher_counts.get(publisher, 0) + 1

    print("\n" + "=" * 60)
    print("PUBLISHER BREAKDOWN")
    print("=" * 60)

    for publisher, count in sorted(publisher_counts.items(), key=lambda x: -x[1]):
        print(f"  {publisher}: {count}")

    # Group by download strategy
    print("\n" + "=" * 60)
    print("DOWNLOAD STRATEGIES")
    print("=" * 60)

    wiley = [(p, d) for p, d, pub in results if "Wiley" in pub]
    elsevier = [(p, d) for p, d, pub in results if "Elsevier" in pub]
    aha = [(p, d) for p, d, pub in results if "AHA" in pub]
    other = [
        (p, d, pub)
        for p, d, pub in results
        if not any(x in pub for x in ["Wiley", "Elsevier", "AHA"])
    ]

    if wiley:
        print(f"\n🔧 WILEY API ({len(wiley)} papers):")
        for pmid, doi in wiley:
            print(f"   {pmid}: {doi}")

    if elsevier:
        print(f"\n🔧 ELSEVIER API ({len(elsevier)} papers):")
        for pmid, doi in elsevier:
            print(f"   {pmid}: {doi}")

    if aha:
        print(f"\n🔧 AHA JOURNALS ({len(aha)} papers) - often free after embargo:")
        for pmid, doi in aha:
            print(f"   {pmid}: {doi}")

    if other:
        print(f"\n📥 OTHER - may need manual or browser fetch ({len(other)} papers):")
        for pmid, doi, pub in other:
            print(f"   {pmid}: {doi or 'NO DOI'} ({pub})")

    # Write detailed CSV for fetch_manager
    output_csv = (
        Path(__file__).parent.parent / "tests/fixtures/pmids/non_pmc_papers.csv"
    )
    with open(output_csv, "w") as f:
        f.write("PMID,DOI,Publisher,URL,Status\n")
        for pmid, doi, publisher in results:
            url = (
                f"https://doi.org/{doi}"
                if doi
                else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            )
            f.write(f"{pmid},{doi or ''},{publisher},{url},\n")

    print(f"\nWrote detailed CSV to: {output_csv}")


if __name__ == "__main__":
    main()
