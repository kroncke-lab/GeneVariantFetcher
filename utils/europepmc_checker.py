#!/usr/bin/env python3
"""
Europe PMC Checker for Failed Papers
Checks Europe PMC for full text availability for specific PMIDs
"""

import requests
import json
import os
from datetime import datetime


def check_europe_pmc(pmid):
    """Check Europe PMC for a given PMID"""
    api_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}&format=json"
    full_text_url = f"https://europepmc.org/article/MED/{pmid}"

    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()

        data = response.json()
        result = data.get("resultList", {}).get("result", [{}])[0]

        return {
            "pmid": pmid,
            "api_url": api_url,
            "full_text_url": full_text_url,
            "title": result.get("title", "N/A"),
            "doi": result.get("doi", "N/A"),
            "pmcid": result.get("pmcid", "N/A"),
            "has_pmc_abstract": result.get("hasAbstract", "N/A"),
            "has_text_mined_terms": result.get("hasTextMinedTerms", "N/A"),
            "has_refs": result.get("hasReferences", "N/A"),
            "is_open_access": result.get("isOpenAccess", "N/A"),
            "fulltext_html_url": f"https://europepmc.org/articles/{result.get('pmcid', '')}"
            if result.get("pmcid")
            else "N/A",
            "fulltext_pdf_url": f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={result.get('pmcid', '')}&blobtype=pdf"
            if result.get("pmcid")
            else "N/A",
            "success": True,
            "raw_api_result": result,
        }

    except Exception as e:
        return {
            "pmid": pmid,
            "api_url": api_url,
            "full_text_url": full_text_url,
            "error": str(e),
            "success": False,
        }


def download_file(url, filename):
    """Download file from URL"""
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()

        with open(filename, "wb") as f:
            f.write(response.content)

        return True, len(response.content)
    except Exception as e:
        return False, str(e)


def main():
    failed_pmids = [
        15840476,
        10973849,
        26496715,
        11854117,
        14661677,
        29650123,
        16922724,
        23631430,
    ]

    output_file = "./output/europepmc_attempt.md"
    download_dir = "./output/europepmc_downloads"

    # Create download directory
    os.makedirs(download_dir, exist_ok=True)

    # Initialize results
    results = []

    print("Checking Europe PMC for failed papers...")

    # Check each PMID
    for pmid in failed_pmids:
        print(f"Processing PMID {pmid}...")
        result = check_europe_pmc(pmid)

        # Attempt to download files if open access
        if result["success"] and result.get("is_open_access") == "Y":
            print(f"  Open access found for PMID {pmid}, downloading...")

            # Download HTML full text
            if result["fulltext_html_url"] and result["fulltext_html_url"] != "N/A":
                filename = f"{download_dir}/{pmid}_europepmc_full.html"
                success, size = download_file(result["fulltext_html_url"], filename)
                if success:
                    result["downloaded_html"] = filename
                    result["html_size"] = f"{size} bytes"
                else:
                    result["html_download_error"] = size

            # Download PDF full text
            if result["fulltext_pdf_url"] and result["fulltext_pdf_url"] != "N/A":
                filename = f"{download_dir}/{pmid}_europepmc_full.pdf"
                success, size = download_file(result["fulltext_pdf_url"], filename)
                if success:
                    result["downloaded_pdf"] = filename
                    result["pdf_size"] = f"{size} bytes"
                else:
                    result["pdf_download_error"] = size

        results.append(result)

    # Generate markdown report
    markdown_content = f"""# Europe PMC Attempt Report
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Summary
- **Total failed PMIDs checked**: {len(failed_pmids)}
- **Successfully retrieved metadata**: {len([r for r in results if r["success"]])}
- **Open access papers found**: {len([r for r in results if r.get("is_open_access") == "Y"])}
- **Files downloaded**: {len([r for r in results if r.get("downloaded_pdf") or r.get("downloaded_html")])}

## Details by PMID

"""

    for result in results:
        markdown_content += f"""### PMID: {result["pmid"]}
- **API Status**: {"✅ Success" if result["success"] else "❌ Failed"}
- **Title**: {result.get("title", "N/A")}
- **DOI**: {result.get("doi", "N/A")}
- **PMC ID**: {result.get("pmcid", "N/A")}
- **Open Access**: {result.get("is_open_access", "N/A")}
- **Europe PMC Link**: [View article]({result.get("full_text_url", "#")})
"""

        if result.get("is_open_access") == "Y":
            markdown_content += "- **Download Status**: Available for download\n"
            if result.get("downloaded_html"):
                markdown_content += (
                    f"  - HTML: {result['downloaded_html']} ({result['html_size']})\n"
                )
            if result.get("downloaded_pdf"):
                markdown_content += (
                    f"  - PDF: {result['downloaded_pdf']} ({result['pdf_size']})\n"
                )
            if result.get("html_download_error"):
                markdown_content += f"  - HTML Error: {result['html_download_error']}\n"
            if result.get("pdf_download_error"):
                markdown_content += f"  - PDF Error: {result['pdf_download_error']}\n"
        else:
            markdown_content += "- **Download Status**: Not open access/Timed out\n"

        if result.get("error"):
            markdown_content += f"- **Error**: {result['error']}\n"

        markdown_content += "\n"

    # Failed downloads section
    failed_access = [
        r for r in results if not (r.get("is_open_access") == "Y") and r["success"]
    ]
    if failed_access:
        markdown_content += "## PMIDs Without Open Access\n\n"
        for result in failed_access:
            if result["success"]:  # Successfully retrieved metadata but not open access
                markdown_content += f"- **PMID {result['pmid']}**: {result.get('title', 'Title unavailable')}\n"
                if result.get("pmcid") and result.get("pmcid") != "N/A":
                    markdown_content += (
                        f"  - Found PMC ID: {result['pmcid']} (wasn't found before)\n"
                    )
        markdown_content += "\n"

    # Write report
    with open(output_file, "w") as f:
        f.write(markdown_content)

    print(f"Report saved to: {output_file}")
    print("\nFinal Results Summary:")
    print(f"- Total PMIDs processed: {len(results)}")
    print(
        f"- With PMC IDs: {len([r for r in results if r.get('pmcid') and r.get('pmcid') != 'N/A'])}"
    )
    print(
        f"- Open access downloads: {len([r for r in results if r.get('downloaded_pdf') or r.get('downloaded_html')])}"
    )


if __name__ == "__main__":
    main()
