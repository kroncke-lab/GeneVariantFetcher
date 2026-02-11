#!/usr/bin/env python3
"""
SP Tier 3 - Direct API Download
Process 70 papers manually using web APIs
"""

import json
import os
import time
import urllib.parse
from datetime import datetime
from pathlib import Path

import requests

tier3_pmids = [
    "19034806",
    "19038855",
    "19065538",
    "19070294",
    "19127321",
    "19136169",
    "19157587",
    "19160088",
    "19169982",
    "19184172",
    "19187913",
    "19215240",
    "19306396",
    "19322600",
    "19352046",
    "19371231",
    "19490267",
    "19695459",
    "19731233",
    "19843919",
    "19996378",
    "20167303",
    "20181576",
    "20197117",
    "20390067",
    "20636320",
    "20975234",
    "21070882",
    "21109023",
    "21130771",
    "21164565",
    "21185499",
    "21216356",
    "21308345",
    "21410720",
    "21419236",
    "21483829",
    "21499742",
    "21951015",
    "22052944",
    "22067087",
    "22104571",
    "22173492",
    "22314138",
    "22338672",
    "22382559",
    "22402334",
    "22407026",
    "22429796",
    "22515331",
    "22727609",
    "22764740",
    "22821100",
    "22882672",
    "22885918",
    "23010577",
    "23134353",
    "23207121",
    "23237912",
    "23277474",
    "23351921",
    "23465283",
    "23546179",
    "23571586",
    "23631430",
    "23864605",
    "23899126",
    "23917959",
    "23981618",
    "23995044",
]


def download_pmc_xml(pmcid, output_path):
    """Download full text XML from PMC"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "pmc", "id": pmcid.replace("PMC", ""), "retmode": "xml"}

    try:
        response = requests.get(base_url, params=params, timeout=30)
        if response.status_code == 200 and "<?xml" in response.text:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(response.text)
            return True, len(response.text)
        return False, response.status_code
    except Exception as e:
        return False, str(e)


def check_pmc_availability(pmid):
    """Check if PMID has PMC availability"""
    api_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}&format=json"

    try:
        response = requests.get(api_url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            result = data.get("resultList", {}).get("result", [{}])[0]

            return {
                "pmid": pmid,
                "pmcid": result.get("pmcid"),
                "doi": result.get("doi"),
                "title": result.get("title", "N/A")[:150],
                "has_pmc": bool(result.get("pmcid")),
                "open_access": result.get("isOpenAccess") == "Y",
                "full_source": "PMC" if result.get("pmcid") else None,
            }
    except Exception as e:
        return {"pmid": pmid, "error": str(e)}

    return {"pmid": pmid, "pmcid": None, "has_pmc": False}


def download_elsevier(doi, api_key, output_dir):
    """Download from Elsevier API"""
    if not doi or not doi.startswith("10.1016"):
        return False, "Not Elsevier DOI"

    url = "https://api.elsevier.com/content/article/doi/"
    headers = {"X-ELS-APIKey": api_key, "Accept": "text/xml"}

    try:
        resp = requests.get(url + doi, headers=headers, timeout=30)
        if resp.status_code == 200:
            filename = os.path.join(
                output_dir, f"elsevier_{doi.replace('/', '_').replace('.', '_')}.xml"
            )
            with open(filename, "w", encoding="utf-8") as f:
                f.write(resp.text)
            return True, filename
        return False, resp.status_code
    except Exception as e:
        return False, str(e)


def main():
    print("=== GVF Tier 3 Direct Download ===")
    print(f"Processing {len(tier3_pmids)} papers")

    # Setup
    output_base = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2/tier3_downloads/")
    output_base.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    batch_dir = output_base / f"direct_api_{timestamp}"
    batch_dir.mkdir(exist_ok=True)

    results = []

    print("\n1. Checking PMC availability...")

    # Process in smaller chunks
    chunk_size = 10
    for i in range(0, len(tier3_pmids), chunk_size):
        chunk = tier3_pmids[i : i + chunk_size]
        chunk_num = i // chunk_size + 1

        print(f"\n--- Chunk {chunk_num} (PMIDs {i + 1}-{i + len(chunk)}) ---")

        for pmid in chunk:
            print(f"  Checking PMID {pmid}...", end=" ")

            # Check PMC
            pmc_info = check_pmc_availability(pmid)

            if pmc_info.get("has_pmc"):
                pmcid = pmc_info["pmcid"]
                filename = f"{pmid}_{pmcid}.xml"
                success, size = download_pmc_xml(pmcid, batch_dir / filename)

                if success:
                    print(f"✅ PMC XML: {size // 1000}KB")
                    results.append(
                        {
                            "pmid": pmid,
                            "source": "PMC",
                            "status": "success",
                            "url": pmc_info["pmcid"],
                            "size": size,
                        }
                    )
                else:
                    print(f"❌ PMC failed: {size}")
                    results.append(
                        {"pmid": pmid, "status": "pmc_failed", "error": str(size)}
                    )

            else:
                # Track for Elsevier/DOI processing
                doi = pmc_info.get("doi")
                if doi and doi.startswith("10.1016"):
                    results.append(
                        {
                            "pmid": pmid,
                            "source": "Elsevier_needed",
                            "doi": doi,
                            "title": pmc_info.get("title", ""),
                        }
                    )
                    print(f"🔗 Elsevier: {doi}")
                else:
                    results.append(
                        {
                            "pmid": pmid,
                            "source": "needs_processing",
                            "doi": doi,
                            "title": pmc_info.get("title", ""),
                        }
                    )
                    print("❓ Manual needed")

        # Small delay
        time.sleep(2)

    # Generate log
    log_data = {
        "timestamp": datetime.now().isoformat(),
        "total_pmids": len(tier3_pmids),
        "results": results,
        "summary": {
            "successful": len([r for r in results if r["status"] == "success"]),
            "needs_elsevier": len(
                [r for r in results if r.get("source") == "Elsevier_needed"]
            ),
            "needs_manual": len(
                [r for r in results if r.get("source") == "needs_processing"]
            ),
        },
    }

    log_file = batch_dir / f"tier3_download_results_{timestamp}.json"
    with open(log_file, "w") as f:
        json.dump(log_data, f, indent=2)

    print("\n=== Results ===")
    summary = log_data["summary"]
    print(f"PMC downloads: {summary['successful']}")
    print(f"Elsevier needed: {summary['needs_elsevier']}")
    print(f"Manual processing: {summary['needs_manual']}")
    print(f"Results: {batch_dir}")

    return summary["successful"] + summary["needs_elsevier"] > 0


if __name__ == "__main__":
    main()
