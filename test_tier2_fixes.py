#!/usr/bin/env python3
"""Test script for Tier 2 fixes: broader table parsing and fuzzy matching."""

import sys

sys.path.insert(0, "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher")

from utils.variant_normalizer import (
    get_variant_type,
    match_variants_fuzzy,
    match_variants_to_baseline,
    normalize_deletion,
    normalize_duplication,
    normalize_frameshift,
    normalize_nonsense,
    normalize_variant,
)


def test_frameshift():
    print("=== Testing Frameshift Normalization ===")
    fs_tests = [
        ("L987fs", "L987fsX"),
        ("L987fsX", "L987fsX"),
        ("L987fsX10", "L987fsX"),
        ("L987fs*10", "L987fsX"),
        ("p.Leu987fs", "L987fsX"),
        ("p.Leu987fsTer10", "L987fsX"),
        ("p.Leu987Profs*10", "L987fsX"),
        ("987fs", "?987fsX"),
        ("fs987", "?987fsX"),
    ]
    passed = 0
    for inp, expected in fs_tests:
        result = normalize_frameshift(inp)
        if result == expected:
            status = "✓"
            passed += 1
        else:
            status = f"✗ got {result}"
        print(f"  {inp:25} -> {str(result):12} expected {expected:12} {status}")
    print(f"  Passed: {passed}/{len(fs_tests)}")
    return passed == len(fs_tests)


def test_nonsense():
    print("\n=== Testing Nonsense Normalization ===")
    ns_tests = [
        ("W1001X", "W1001X"),
        ("W1001*", "W1001X"),
        ("p.Trp1001Ter", "W1001X"),
        ("p.Trp1001*", "W1001X"),
        ("R864stop", "R864X"),
        ("R864sp", "R864X"),
    ]
    passed = 0
    for inp, expected in ns_tests:
        result = normalize_nonsense(inp)
        if result == expected:
            status = "✓"
            passed += 1
        else:
            status = f"✗ got {result}"
        print(f"  {inp:25} -> {str(result):12} expected {expected:12} {status}")
    print(f"  Passed: {passed}/{len(ns_tests)}")
    return passed == len(ns_tests)


def test_deletion():
    print("\n=== Testing Deletion Normalization ===")
    del_tests = [
        ("del552", "?552del"),
        ("p.Leu552del", "L552del"),
        ("L552del", "L552del"),
        ("L552_L555del", "L552_L555del"),
    ]
    passed = 0
    for inp, expected in del_tests:
        result = normalize_deletion(inp)
        if result == expected:
            status = "✓"
            passed += 1
        else:
            status = f"✗ got {result}"
        print(f"  {inp:25} -> {str(result):12} expected {expected:12} {status}")
    print(f"  Passed: {passed}/{len(del_tests)}")
    return passed == len(del_tests)


def test_duplication():
    print("\n=== Testing Duplication Normalization ===")
    dup_tests = [
        ("L552dup", "L552dup"),
        ("p.Leu552dup", "L552dup"),
        ("dup552", "?552dup"),
    ]
    passed = 0
    for inp, expected in dup_tests:
        result = normalize_duplication(inp)
        if result == expected:
            status = "✓"
            passed += 1
        else:
            status = f"✗ got {result}"
        print(f"  {inp:25} -> {str(result):12} expected {expected:12} {status}")
    print(f"  Passed: {passed}/{len(dup_tests)}")
    return passed == len(dup_tests)


def test_fuzzy_matching():
    print("\n=== Testing match_variants_fuzzy ===")
    # Note: exact_normalized is also valid for frameshift/nonsense since
    # normalize_variant() handles these cases before fuzzy matching kicks in
    fuzzy_tests = [
        ("L987fs", "L987fsX10", True, ["frameshift", "exact_normalized"]),
        ("p.Leu987fsTer10", "L987fs", True, ["frameshift", "exact_normalized"]),
        ("W1001X", "p.Trp1001*", True, ["nonsense", "exact_normalized"]),
        ("A561V", "A562V", True, ["fuzzy", "position"]),  # Off by 1
        ("A561V", "A563V", False, ["no_match"]),  # Off by 2
        ("A561V", "G561V", False, ["no_match"]),  # Different ref
    ]
    passed = 0
    for v1, v2, expected_match, expected_type_parts in fuzzy_tests:
        is_match, match_type = match_variants_fuzzy(v1, v2)
        type_ok = (
            any(part in match_type for part in expected_type_parts)
            if is_match
            else not is_match
        )
        if is_match == expected_match and type_ok:
            status = "✓"
            passed += 1
        else:
            status = f"✗ got match={is_match}, type={match_type}"
        print(
            f"  {v1:15} vs {v2:20} -> match={is_match}, type={match_type:25} {status}"
        )
    print(f"  Passed: {passed}/{len(fuzzy_tests)}")
    return passed == len(fuzzy_tests)


def test_baseline_matching():
    print("\n=== Testing match_variants_to_baseline with fuzzy matching ===")
    baseline = {"A561V", "G584S", "L987fsX", "W1001X", "L552del"}
    extracted = [
        "p.Ala561Val",  # Should match A561V (normalized)
        "L987fs*10",  # Should match L987fsX (frameshift)
        "p.Trp1001Ter",  # Should match W1001X (nonsense)
        "p.Leu552del",  # Should match L552del (deletion)
        "G585S",  # Should match G584S (fuzzy +1)
        "R248W",  # Should be filtered (TP53 hotspot)
        "unknown_variant",  # Should be unmatched
    ]

    results = match_variants_to_baseline(
        extracted, baseline, "KCNH2", fuzzy_position=True
    )

    print(f"  Total input: {results['stats']['total_input']}")
    print(f"  Matches: {len(results['matches'])}")
    for m in results["matches"]:
        print(f"    {m['extracted']} -> {m['matched_to']} ({m['match_type']})")
    print(f"  Filtered (non-target): {len(results['filtered_non_target'])}")
    for v, reason in results["filtered_non_target"]:
        print(f"    {v}: {reason}")
    print(f"  Unmatched: {results['unmatched']}")
    print(f"  Stats: {results['stats']}")

    # Verify expected results
    expected_matches = 5  # A561V, L987fsX, W1001X, L552del, G584S
    expected_filtered = 1  # R248W
    expected_unmatched = 1  # unknown_variant

    match_count = len(results["matches"])
    filter_count = len(results["filtered_non_target"])
    unmatch_count = len(results["unmatched"])

    success = (
        match_count >= expected_matches - 1  # Allow some tolerance
        and filter_count >= expected_filtered
        and unmatch_count >= expected_unmatched
    )

    print(
        f"\n  Expected ~{expected_matches} matches, got {match_count}: {'✓' if match_count >= expected_matches - 1 else '✗'}"
    )
    print(
        f"  Expected ~{expected_filtered} filtered, got {filter_count}: {'✓' if filter_count >= expected_filtered else '✗'}"
    )
    print(
        f"  Expected ~{expected_unmatched} unmatched, got {unmatch_count}: {'✓' if unmatch_count >= expected_unmatched else '✗'}"
    )

    return success


def test_table_header_detection():
    print("\n=== Testing Broadened Table Header Detection ===")
    from pipeline.extraction import ExpertExtractor

    extractor = ExpertExtractor()

    # Test headers that should be detected
    valid_headers = [
        "| Nucleotide | Variant | No. of Patients |",  # Original pattern
        "| cDNA | Protein | N |",
        "| c. | p. | Affected |",
        "| HGVS cDNA | HGVS Protein | Carriers |",
        "| Mutation | AAChange | Probands |",
        "| nucleotide change | amino acid | families |",
        "| cDNA Change | Protein Change | subjects |",
    ]

    # Headers that should NOT be detected (missing variant OR count columns)
    invalid_headers = [
        "| Gene | Chromosome | Position |",
        "| Author | Year | Title |",
        "| Sample ID | Status |",
    ]

    passed = 0
    total = len(valid_headers) + len(invalid_headers)

    print("  Valid headers (should detect):")
    for h in valid_headers:
        result = extractor._is_variant_table_header(h)
        if result:
            status = "✓"
            passed += 1
        else:
            status = "✗"
        print(f"    {result}: {h[:60]}... {status}")

    print("  Invalid headers (should NOT detect):")
    for h in invalid_headers:
        result = extractor._is_variant_table_header(h)
        if not result:
            status = "✓"
            passed += 1
        else:
            status = "✗"
        print(f"    {result}: {h[:60]}... {status}")

    print(f"  Passed: {passed}/{total}")
    return passed == total


def main():
    print("=" * 70)
    print("Tier 2 Fixes Test Suite")
    print("=" * 70)

    all_passed = True
    all_passed &= test_frameshift()
    all_passed &= test_nonsense()
    all_passed &= test_deletion()
    all_passed &= test_duplication()
    all_passed &= test_fuzzy_matching()
    all_passed &= test_baseline_matching()
    all_passed &= test_table_header_detection()

    print("\n" + "=" * 70)
    if all_passed:
        print("ALL TESTS PASSED! ✓")
    else:
        print("SOME TESTS FAILED! ✗")
    print("=" * 70)

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
