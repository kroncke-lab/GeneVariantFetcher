#!/usr/bin/env python3
"""
Test script for improved variant normalizer.
Tests against real KCNH2 data from missing_in_sqlite.csv
"""

import csv
import os
import sys

# Add project root to path
sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from utils.variant_normalizer import VariantNormalizer


def load_test_variants(csv_file):
    """Load variants from CSV file for testing."""
    variants = []
    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["excel_variant_raw"] and row["missing_in_sqlite"] == "True":
                variants.append(
                    {
                        "raw": row["excel_variant_raw"],
                        "expected_normalized": row["excel_variant_norm"],
                        "gene": "KCNH2",
                    }
                )
    return variants


def test_normalization_coverage():
    """Test normalization coverage on real data."""
    normalizer = VariantNormalizer("KCNH2")

    # Get path to test data
    test_file = "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/missing_in_sqlite.csv"

    if not os.path.exists(test_file):
        print("Test file not found, using synthetic test cases...")
        test_cases = [
            {"raw": "A561V", "expected": "p.Ala561Val"},
            {"raw": "A561T", "expected": "p.Ala561Thr"},
            {"raw": "A1017fsX", "expected": "p.Ala1017fs*"},
            {"raw": "R752W", "expected": "p.Arg752Trp"},
            {"raw": "D864sp", "expected": "p.Asp864*"},
            {"raw": "G184Del", "expected": "p.Gly184del"},
            {"raw": "G189Ins", "expected": "p.Gly189ins"},
            {"raw": "p.Ala561Val", "expected": "p.Ala561Val"},  # Already correct
            {"raw": "p.A561V", "expected": "p.Ala561Val"},  # Needs expansion
        ]
    else:
        test_cases = load_test_variants(test_file)[:50]  # First 50 for testing

    # Track statistics
    stats = {"total": len(test_cases), "normalized": 0, "failed": 0, "details": {}}

    print("Testing Variant Normalizer Improvements")
    print("=" * 50)

    for test_case in test_cases:
        raw = test_case["raw"]
        expected = test_case["expected_normalized"]
        normalized = normalizer.normalize_protein(raw)

        success = normalized is not None
        is_correct = normalized == expected

        stats["details"][raw] = {
            "raw": raw,
            "normalized": normalized,
            "expected": expected if "expected" in test_case else "N/A",
            "success": success,
            "correct": is_correct,
        }

        if success:
            stats["normalized"] += 1
            if not is_correct:
                stats["failed"] += 1

    # Print results
    print(f"Total test cases: {stats['total']}")
    print(
        f"Successfully normalized: {stats['normalized']} ({100 * stats['normalized'] / stats['total']:.1f}%)"
    )
    print(f"Failed: {stats['failed']}")
    print()

    print("Sample results:")
    for i, (raw, details) in enumerate(stats["details"].items()):
        if i >= 10:  # Show first 10
            print("...")
            break
        status = "✓" if details["success"] else "✗"
        print(
            f"{status} {details['raw']} → {details['normalized']} (expected: {details['expected']})"
        )

    # Print failed cases for debugging
    failed_cases = [
        d for d in stats["details"].values() if not d["success"] or not d["correct"]
    ]
    if failed_cases:
        print(f"\nFailed cases ({len(failed_cases)}):")
        for case in failed_cases[:5]:  # Show first 5
            print(
                f"  {case['raw']} → {case['normalized']} (expected: {case['expected']})"
            )

    return stats


if __name__ == "__main__":
    results = test_normalization_coverage()
    print(
        f"\nOverall improvement: {100 * results['normalized'] / results['total']:.1f}% normalization rate"
    )
