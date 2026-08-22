#!/usr/bin/env python3
"""
Tests for variant normalization functions.

Includes test cases migrated from test_tier2_fixes.py (root) covering:
- Frameshift, nonsense, deletion, duplication normalization
- Fuzzy variant matching
- Baseline matching with normalization
"""

import pytest

from utils.variant_normalizer import (
    VariantNormalizer,
    match_variants_fuzzy,
    match_variants_to_baseline,
    normalize_deletion,
    normalize_duplication,
    normalize_frameshift,
    normalize_nonsense,
    normalize_variant,
    preferred_alias_protein,
    protein_substitution_frameshift_alias,
    variants_match,
)


# ---------------------------------------------------------------------------
# Frameshift normalization
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw, expected",
    [
        ("L987fs", "L987fsX"),
        ("L987fsX", "L987fsX"),
        ("L987fsX10", "L987fsX"),
        ("L987fs*10", "L987fsX"),
        ("p.Leu987fs", "L987fsX"),
        ("p.Leu987fsTer10", "L987fsX"),
        ("p.Leu987Profs*10", "L987fsX"),
        ("p.M815Wfs*10", "M815fsX"),
        ("p.S1722YfsTer4", "S1722fsX"),
        ("987fs", "?987fsX"),
        ("fs987", "?987fsX"),
    ],
)
def test_normalize_frameshift(raw, expected):
    assert normalize_frameshift(raw) == expected


def test_normalize_variant_unifies_one_and_three_letter_frameshifts():
    assert normalize_variant("p.M815Wfs*10", "BRCA2") == normalize_variant(
        "p.Met815Trpfs*10", "BRCA2"
    )


def test_normalize_variant_unifies_synonymous_equals_and_repeated_residue():
    assert normalize_variant("p.Asn961=", "BRCA1") == "N961N"
    assert normalize_variant("N961=", "BRCA1") == "N961N"
    assert normalize_variant("p.Asn961Asn", "BRCA1") == "N961N"


@pytest.mark.parametrize(
    "short, frameshift",
    [
        ("V340G", "p.Val340Glyfs*6"),
        ("P.LYS339F", "p.Lys339PhefsTer8"),
        ("V340G", "p.(Val340Glyfs*6)"),
    ],
)
def test_explicit_frameshift_alias_requires_matching_ref_position_and_alt(
    short, frameshift
):
    assert protein_substitution_frameshift_alias(short, frameshift)
    assert protein_substitution_frameshift_alias(frameshift, short)
    assert preferred_alias_protein(short, frameshift) == frameshift
    assert preferred_alias_protein(frameshift, short) == frameshift


@pytest.mark.parametrize(
    "left, right",
    [
        ("V340A", "p.Val340Glyfs*6"),
        ("V341G", "p.Val340Glyfs*6"),
        ("V340G", "p.Val340fs*6"),
    ],
)
def test_frameshift_alias_never_matches_alt_mismatch_or_codon_only(left, right):
    assert not protein_substitution_frameshift_alias(left, right)


# ---------------------------------------------------------------------------
# Nonsense normalization
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw, expected",
    [
        ("W1001X", "W1001X"),
        ("W1001*", "W1001X"),
        ("p.Trp1001Ter", "W1001X"),
        ("p.Trp1001*", "W1001X"),
        ("R864stop", "R864X"),
        ("R864sp", "R864X"),
    ],
)
def test_normalize_nonsense(raw, expected):
    assert normalize_nonsense(raw) == expected


# ---------------------------------------------------------------------------
# Deletion normalization
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw, expected",
    [
        ("del552", "?552del"),
        ("p.Leu552del", "L552del"),
        ("L552del", "L552del"),
        ("L552_L555del", "L552_L555del"),
    ],
)
def test_normalize_deletion(raw, expected):
    assert normalize_deletion(raw) == expected


# ---------------------------------------------------------------------------
# Duplication normalization
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw, expected",
    [
        ("L552dup", "L552dup"),
        ("p.Leu552dup", "L552dup"),
        ("dup552", "?552dup"),
    ],
)
def test_normalize_duplication(raw, expected):
    assert normalize_duplication(raw) == expected


# ---------------------------------------------------------------------------
# Fuzzy matching
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "v1, v2, expected_match",
    [
        ("L987fs", "L987fsX10", True),
        ("p.Leu987fsTer10", "L987fs", True),
        ("W1001X", "p.Trp1001*", True),
        ("A561V", "A562V", True),  # Off by 1
        ("A561V", "A563V", False),  # Off by 2 -> no match
        ("A561V", "G561V", False),  # Different ref AA
    ],
)
def test_match_variants_fuzzy(v1, v2, expected_match):
    is_match, _match_type = match_variants_fuzzy(v1, v2)
    assert is_match == expected_match


# ---------------------------------------------------------------------------
# Baseline matching
# ---------------------------------------------------------------------------
def test_match_variants_to_baseline():
    baseline = {"A561V", "G584S", "L987fsX", "W1001X", "L552del"}
    extracted = [
        "p.Ala561Val",  # -> A561V (normalized)
        "L987fs*10",  # -> L987fsX (frameshift)
        "p.Trp1001Ter",  # -> W1001X (nonsense)
        "p.Leu552del",  # -> L552del (deletion)
        "G585S",  # -> G584S (fuzzy +1)
        "R248W",  # TP53 hotspot -> filtered
        "unknown_variant",  # unmatched
    ]

    results = match_variants_to_baseline(
        extracted, baseline, "KCNH2", fuzzy_position=True
    )

    assert len(results["matches"]) >= 4  # At least the exact-normalized ones
    assert len(results["filtered_non_target"]) >= 1
    assert len(results["unmatched"]) >= 1


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("A561V*", "A561V"),
        ("p.Arg231Cys*", "R231C"),
        ("M291T*", "M291T"),
        ("R864*", "R864X"),
        ("R176W(het)", "R176W"),
        ("A561V/WT", "A561V"),
        ("G584S(hom)", "G584S"),
        ("R534C/+", "R534C"),
        ("a561v", "A561V"),
        ("G184Del", "G184del"),
    ],
)
def test_normalize_variant_migrated_demo_cases(raw, expected):
    """Keep the executable assertions formerly hidden in the demo runner."""

    assert normalize_variant(raw) == expected


def test_variant_normalizer_get_all_forms_migrated_demo_case():
    forms = VariantNormalizer("KCNH2").get_all_forms("p.Arg534Cys")

    assert forms["single"] == "R534C"
    assert forms["three"] == "p.Arg534Cys"


@pytest.mark.parametrize(
    "variant, expected_non_target",
    [
        ("R248W", True),
        ("G12D", True),
        ("V600E", True),
        ("P2006A", True),
        ("A561V", False),
    ],
)
def test_non_target_detection_migrated_demo_cases(variant, expected_non_target):
    is_non_target, _reason = VariantNormalizer("KCNH2").is_non_target_variant(variant)

    assert is_non_target is expected_non_target


@pytest.mark.parametrize(
    "left, right, expected",
    [
        ("p.Ala561Val", "A561V", True),
        ("p.Arg534Cys", "R534C", True),
        ("A561V", "G584S", False),
    ],
)
def test_variants_match_migrated_demo_cases(left, right, expected):
    assert variants_match(left, right) is expected
