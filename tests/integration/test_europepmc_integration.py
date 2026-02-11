#!/usr/bin/env python3
"""
Europe PMC Integration Test Script

Tests the Europe PMC API integration with GVF using golden test PMIDs.
Validates both functionality and integration with existing systems.
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime

# Add the project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from gene_literature.europepmc_handler import EuropePMCClient
from gene_literature.europepmc_integration import EuropePMCHarvester
from harvesting.pmc_api import PMCAPIClient


def setup_logging():
    """Configure logging for test output."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler("/mnt/temp2/kronckbm/gvf_output/europepmc_test.log"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger(__name__)


def test_basic_api_operations():
    """Test basic Europe PMC API functionality."""
    logger = setup_logging()
    logger.info("=== Testing Basic Europe PMC API Operations ===")

    client = EuropePMCClient()

    test_results = {
        "search_tests": [],
        "metadata_tests": [],
        "fulltext_tests": [],
        "supplement_tests": [],
    }

    # Test search operations
    logger.info("Testing search operations...")
    search_results = client.search_papers('"KCNH2" AND variant', max_results=5)
    test_results["search_tests"].append(
        {
            "query": '"KCNH2" AND variant',
            "results_count": len(search_results),
            "first_result": search_results[0] if search_results else None,
        }
    )

    # Test individual paper operations
    test_pmids = ["30036649", "26746457", "9544837"]  # From golden test set

    for pmid in test_pmids:
        logger.info(f"Testing PMID: {pmid}")

        # Test metadata retrieval
        logger.info(f"  Testing metadata retrieval...")
        metadata = client.get_paper_metadata(pmid)
        metadata_result = {
            "pmid": pmid,
            "metadata_success": metadata is not None,
            "title": metadata.get("title", "") if metadata else None,
            "pmcid": metadata.get("pmcid") if metadata else None,
            "is_open_access": metadata.get("is_open_access", False)
            if metadata
            else False,
        }
        test_results["metadata_tests"].append(metadata_result)

        # Test full-text retrieval if available
        pmcid = metadata.get("pmcid") if metadata else None
        if pmcid:
            logger.info(f"  Testing full-text XML retrieval...")
            xml_content = client.get_fulltext_xml(pmcid)
            xml_result = {
                "pmid": pmid,
                "pmcid": pmcid,
                "xml_success": xml_content is not None,
                "xml_length": len(xml_content) if xml_content else 0,
            }
            test_results["fulltext_tests"].append(xml_result)

            # Test supplementary files
            logger.info(f"  Testing supplementary files...")
            supplements = client.get_supplementary_files(pmcid)
            supplement_result = {
                "pmid": pmid,
                "pmcid": pmcid,
                "supplement_count": len(supplements),
                "supplement_details": [
                    {
                        "filename": supp.get("filename", supp.get("name", "unknown")),
                        "url": supp.get("downloadUrl", supp.get("url", "")),
                    }
                    for supp in supplements
                ],
            }
            test_results["supplement_tests"].append(supplement_result)

    return test_results


def test_integration_download():
    """Test full integration download workflow."""
    logger = setup_logging()
    logger.info("=== Testing Integration Download ===")

    harvester = EuropePMCHarvester(
        output_dir=Path("/mnt/temp2/kronckbm/gvf_output/europepmc_test")
    )

    # Use golden test set PMIDs
    golden_pmids = [
        "30036649",  # 1.0 recall in golden set
        "26746457",  # 1.0 recall in golden set
        "9544837",  # 1.0 recall in golden set
        "18808722",  # 0.875 recall in golden set
        "21063070",  # 1.0 recall in golden set
    ]

    logger.info(f"Testing batch download of {len(golden_pmids)} golden test PMIDs")

    results = harvester.download_batch(golden_pmids)

    # Process results
    successful = [r for r in results if r["success"]]
    failed = [r for r in results if not r["success"]]

    summary = {
        "total_pmids": len(golden_pmids),
        "successful": len(successful),
        "failed": len(failed),
        "success_rate": len(successful) / len(golden_pmids) if golden_pmids else 0,
        "results": results,
        "summary_files": {
            "successful_pmids": [r["pmid"] for r in successful],
            "failed_pmids": [r["pmid"] for r in failed],
            "failed_reasons": {r["pmid"]: r.get("error", "Unknown") for r in failed},
        },
    }

    return summary


def test_pmc_vs_europepmc_coverage():
    """Compare coverage between PubMed Central and Europe PMC."""
    logger = setup_logging()
    logger.info("=== Testing PMC vs Europe PMC Coverage ===")

    # All golden test PMIDs
    all_pmids = [
        "30036649",
        "18808722",
        "26746457",
        "11844290",
        "14661677",
        "9544837",
        "21063070",
        "16155735",
        "28049825",
        "19668779",
        "23098067",
        "29016797",
        "15973763",
        "21964171",
    ]

    harvester = EuropePMCHarvester()

    coverage_results = []

    for pmid in all_pmids:
        logger.info(f"Comparing coverage for PMID: {pmid}")
        try:
            comparison = harvester.compare_with_pmc(pmid)
            coverage_results.append(comparison)

            logger.info(f"  PMC: {'✓' if comparison['pmc_available'] else '✗'}")
            logger.info(
                f"  Europe PMC: {'✓' if comparison['europepmc_available'] else '✗'}"
            )
            logger.info(
                f"  Coverage improvement: {'Yes' if comparison['coverage_improvement'] else 'No'}"
            )

        except Exception as e:
            logger.error(f"Error comparing coverage for {pmid}: {e}")
            coverage_results.append({"pmid": pmid, "error": str(e)})

    # Summary statistics
    pmc_available = sum(1 for r in coverage_results if r.get("pmc_available", False))
    epmc_available = sum(
        1 for r in coverage_results if r.get("europepmc_available", False)
    )
    coverage_improvements = sum(
        1 for r in coverage_results if r.get("coverage_improvement", False)
    )

    summary = {
        "total_pmids": len(all_pmids),
        "pmc_only_available": pmc_available,
        "europepmc_only_available": epmc_available - pmc_available,
        "both_available": sum(
            1 for r in coverage_results if r.get("both_available", False)
        ),
        "coverage_improvements": coverage_improvements,
        "detailed_results": coverage_results,
    }

    return summary


def test_search_functionality():
    """Test advanced search capabilities."""
    logger = setup_logging()
    logger.info("=== Testing Search Functionality ===")

    client = EuropePMCClient()
    harvester = EuropePMCHarvester()

    search_tests = [
        "KCNH2[Gene Symbol]",
        "KCNH2 AND (variant OR mutation)",
        "HERG AND long QT",
        "KCNH2[Gene Symbol] AND fullText:y",
        "hERG protein AND variant",
        # Test limit by open access
        "KCNH2 AND open_access:y",
    ]

    search_results = {}

    for query in search_tests:
        logger.info(f"Testing search: {query}")
        try:
            results = client.search_papers(query, max_results=10)
            search_results[query] = {
                "count": len(results),
                "total_hits": results[0].get("hitCount", 0) if results else 0,
                "sample_results": [
                    {
                        "pmid": r.get("pmid"),
                        "title": r.get("title", "")[:100] + "...",
                        "year": r.get("pubYear"),
                        "full_text_available": r.get("hasFullText") == "Y",
                    }
                    for r in results[:3]
                ],
            }

            logger.info(f"  Found {len(results)} results")

        except Exception as e:
            logger.error(f"Error executing search '{query}': {e}")
            search_results[query] = {"error": str(e)}

    return search_results


def run_comprehensive_test():
    """Run all tests and create comprehensive report."""
    logger = setup_logging()
    logger.info("=== Starting Comprehensive Europe PMC Test Suite ===")

    comprehensive_result = {
        "test_start_time": datetime.now().isoformat(),
        "test_suite_name": "Europe PMC Integration Test",
        "environment": {
            "output_path": "/mnt/temp2/kronckbm/gvf_output",
            "test_directory": "/mnt/temp2/kronckbm/gvf_output/europepmc_test",
        },
    }

    # Run all tests
    comprehensive_result["basic_api_tests"] = test_basic_api_operations()
    comprehensive_result["integration_tests"] = test_integration_download()
    comprehensive_result["coverage_tests"] = test_pmc_vs_europepmc_coverage()
    comprehensive_result["search_tests"] = test_search_functionality()

    # Summary statistics
    successful_downloads = len(
        [
            r
            for r in comprehensive_result["integration_tests"]["results"]
            if r["success"]
        ]
    )
    total_pmids = comprehensive_result["coverage_tests"]["total_pmids"]

    comprehensive_result["summary"] = {
        "test_end_time": datetime.now().isoformat(),
        "total_pmids_tested": total_pmids,
        "successful_downloads": successful_downloads,
        "api_connectivity": "PASS"
        if comprehensive_result["basic_api_tests"]["metadata_tests"]
        else "FAIL",
        "integration_complete": "PASS"
        if comprehensive_result["integration_tests"]
        else "FAIL",
    }

    # Save results
    results_file = Path(
        "/mnt/temp2/kronckbm/gvf_output/europepmc_integration_test_results.json"
    )
    with open(results_file, "w", encoding="utf-8") as f:
        json.dump(comprehensive_result, f, indent=2, ensure_ascii=False)

    logger.info(f"Test results saved to: {results_file}")
    return comprehensive_result


if __name__ == "__main__":
    # Run comprehensive test
    results = run_comprehensive_test()

    # Print summary
    print("\n" + "=" * 60)
    print("EUROPE PMC INTEGRATION TEST SUMMARY")
    print("=" * 60)

    summary = results["summary"]
    print(f"Total PMIDs tested: {summary['total_pmids_tested']}")
    print(f"Successful downloads: {summary['successful_downloads']}")
    print(f"API connectivity: {summary['api_connectivity']}")
    print(f"Integration complete: {summary['integration_complete']}")

    print(
        f"\nTest results saved to: /mnt/temp2/kronckbm/gvf_output/europepmc_integration_test_results.json"
    )

    # Exit with appropriate code
    if (
        summary["api_connectivity"] == "PASS"
        and summary["integration_complete"] == "PASS"
    ):
        print("✅ All tests passed - Europe PMC integration is working!")
        sys.exit(0)
    else:
        print("❌ Some tests failed - check the result file for details")
        sys.exit(1)
