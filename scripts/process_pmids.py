#!/usr/bin/env python3
"""
Process paywalled PMIDs for institutional downloads via Vanderbilt proxy
Target: Elsevier (ScienceDirect) → Wiley → Springer
"""

import os
import sys

import pandas as pd


def main():
    # Path setup
    gvf_repo = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher"
    output_dir = "/mnt/temp2/kronckbm/gvf_output/institutional_downloads/2026-02-07"

    # Read paywalled/missing CSV
    csv_path = f"{gvf_repo}/output/KCNH2/KCNH2/20260205_161223/pmc_fulltext/paywalled_missing.csv"

    if not os.path.exists(csv_path):
        print("ERROR: paywalled_missing.csv not found")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    pmids = [
        str(p).strip()
        for p in df["PMID"].unique()
        if str(p).strip() != "" and str(p) != "nan"
    ]

    print(f"Found {len(pmids)} unique PMIDs in paywalled/missing list")

    # Classify by provider
    elsevier_pmids = []
    wiley_pmids = []
    springer_pmids = []
    unknown_pmids = []

    for pmid in pmids:
        # Check URL and DOI patterns for provider classification
        pmid_rows = df[df["PMID"].astype(str) == str(pmid)]

        is_elsevier = False
        is_wiley = False
        is_springer = False

        for _, row in pmid_rows.iterrows():
            url = str(row.get("URL", "")).lower()
            doi = str(row.get("Reason", "")).lower()  # DOI often in Reason field

            if any(
                indicator in url or doi
                for indicator in [
                    "elsevier",
                    "sciencedirect",
                    "10.1016",
                    "10.1093/oxfordjournals",
                    "linkinghub.elsevier",
                ]
            ):
                is_elsevier = True
                break
            elif any(
                indicator in url or doi
                for indicator in ["wiley", "10.1002", "10.1111", "onlinelibrary.wiley"]
            ):
                is_wiley = True
                break
            elif any(
                indicator in url or doi
                for indicator in ["springer", "10.1007", "10.1365", "link.springer"]
            ):
                is_springer = True
                break

        if is_elsevier:
            elsevier_pmids.append(pmid)
        elif is_wiley:
            wiley_pmids.append(pmid)
        elif is_springer:
            springer_pmids.append(pmid)
        else:
            unknown_pmids.append(pmid)

    # Provider counts
    print(f"Elsevier candidates: {len(set(elsevier_pmids))}")
    print(f"Wiley candidates: {len(set(wiley_pmids))}")
    print(f"Springer candidates: {len(set(springer_pmids))}")
    print(f"Unknown providers: {len(set(unknown_pmids))}")

    # Create prioritized download list
    download_order = (
        list(set(elsevier_pmids))
        + list(set(wiley_pmids))
        + list(set(springer_pmids))
        + list(set(unknown_pmids))
    )

    # Save target PMIDs
    target_pmids_path = f"{output_dir}/target_pmids.txt"
    with open(target_pmids_path, "w") as f:
        f.write("\n".join(download_order))

    # Save provider mapping
    provider_mapping = {
        "elsevier": list(set(elsevier_pmids)),
        "wiley": list(set(wiley_pmids)),
        "springer": list(set(springer_pmids)),
        "unknown": list(set(unknown_pmids)),
    }

    import json

    with open(f"{output_dir}/provider_mapping.json", "w") as f:
        json.dump(provider_mapping, f, indent=2)

    print(f"Total target PMIDs: {len(download_order)}")
    print(f"Target PMIDs written to {target_pmids_path}")

    return download_order, provider_mapping


if __name__ == "__main__":
    download_order, provider_mapping = main()
