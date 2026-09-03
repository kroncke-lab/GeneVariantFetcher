#!/usr/bin/env python3
"""
Tests for variant normalization functions.

Includes test cases migrated from test_tier2_fixes.py (root) covering:
- Frameshift, nonsense, deletion, duplication normalization
- Fuzzy variant matching
- Baseline matching with normalization
"""

import pytest

import utils.variant_normalizer as variant_normalizer_module
from utils.variant_normalizer import (
    VariantNormalizer,
    _preprocess_cdna_token,
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


def test_gold_free_alias_policy_bypasses_even_a_warm_file_cache(tmp_path, monkeypatch):
    (tmp_path / "kcnh2_variant_aliases.json").write_text(
        '{"aliases": {"C.1A>G": "A1V"}}'
    )
    monkeypatch.setattr(variant_normalizer_module, "_DATA_DIR", tmp_path)
    monkeypatch.setattr(variant_normalizer_module, "_alias_cache", {})
    monkeypatch.delenv("GVF_DISABLE_GOLD_DERIVED_ALIASES", raising=False)

    assert variant_normalizer_module._load_gene_aliases("KCNH2") == {"C.1A>G": "A1V"}
    monkeypatch.setenv("GVF_DISABLE_GOLD_DERIVED_ALIASES", "1")
    assert variant_normalizer_module._load_gene_aliases("KCNH2") == {}


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


@pytest.mark.parametrize(
    "raw",
    ["p.R360_Q361dupRQ", "R360_Q361dupRQ", "p.Arg360_Gln361dupArgGln"],
    ids=["single_letter_prefixed", "single_letter_bare", "three_letter"],
)
def test_range_duplication_drops_a_run_that_restates_the_range(raw):
    """A run whose length equals the coordinate span only restates it, so the
    spellings must converge or the verbatim one stays disjoint from every
    matcher form."""

    assert normalize_duplication(raw) == "R360_Q361dup"
    assert normalize_variant(raw, "KCNQ1") == normalize_variant("R360_Q361dup", "KCNQ1")


@pytest.mark.parametrize(
    "raw",
    ["p.R360_Q361dupQKQR", "R360_Q361dupQKQR", "R360_Q361DUPQKQR"],
    ids=["prefixed", "bare", "gold_spelling"],
)
def test_range_duplication_keeps_a_run_longer_than_the_range(raw):
    """QKQR is four residues across a two-residue range: a distinct insertion
    allele, not a restatement. Collapsing it pooled separate identities, their
    carrier counts and their dedup keys under one key."""

    assert normalize_duplication(raw) == "R360_Q361dupQKQR"
    assert normalize_variant(raw, "KCNQ1") != normalize_variant("R360_Q361dup", "KCNQ1")


def test_unparsed_protein_fall_through_drops_the_presentation_prefix():
    assert normalize_variant("p.Gln1507_Pro1509del", "SCN5A") == "GLN1507_PRO1509DEL"


# ---------------------------------------------------------------------------
# cDNA character tolerance
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw",
    [
        "c.[999–424_1338 + 81del]",
        "c.999−424_1338+81del",
        "c.[999-424_1338+81del]+[=]",
    ],
    ids=["bracketed_en_dash_spaces", "minus_sign", "bracket_plus_reference_allele"],
)
def test_cdna_typography_is_normalized_before_grammar(raw):
    assert VariantNormalizer("SCN5A").normalize_cdna(raw) == "c.999-424_1338+81del"


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        # Literature writes the indel range separator as a hyphen. Adjacent
        # coordinates are unambiguous: an insertion sits between two flanking
        # neighbouring nucleotides, so this is presentation, not identity.
        ("c.2550-2551insTG", "c.2550_2551insTG"),
        ("c.2550-2551del", "c.2550_2551del"),
        ("c.2550-2551dup", "c.2550_2551dup"),
        # Intronic acceptor offsets keep the minus: the right operand is an
        # offset, not the neighbouring coordinate.
        ("c.1234-1del", "c.1234-1del"),
        ("c.1234-12dup", "c.1234-12dup"),
        # A wider span cannot be told apart from an offset, so it fails closed.
        ("c.100-200del", "c.100-200del"),
        # Operands carrying their own sign are never touched.
        ("c.999-424_1338+81del", "c.999-424_1338+81del"),
    ],
    ids=[
        "adjacent_ins_is_a_range",
        "adjacent_del_is_a_range",
        "adjacent_dup_is_a_range",
        "acceptor_offset_del_kept",
        "acceptor_offset_dup_kept",
        "wide_span_fails_closed",
        "offset_range_untouched",
    ],
)
def test_hyphen_indel_range_separator(raw, expected):
    assert _preprocess_cdna_token(raw) == expected


def test_hyphen_range_repair_does_not_touch_substitutions():
    """``c.2550-2A>G`` is a splice-acceptor substitution, not a range."""
    assert _preprocess_cdna_token("c.2550-2A>G") == "c.2550-2A>G"
    assert _preprocess_cdna_token("c.2550-2551A>G") == "c.2550-2551A>G"


def test_hyphen_range_repair_makes_the_two_spellings_one_variant():
    """The whole point: the two spellings must resolve to one identity."""
    assert variants_match("c.2550-2551insTG", "c.2550_2551insTG", "SCN5A")
    assert not variants_match("c.1234-1del", "c.1234_1235del", "SCN5A")


def test_multi_allele_brackets_are_left_alone():
    """``c.[X];[Y]`` carries real two-allele semantics; do not unwrap it."""
    raw = "c.[1234A>G];[2000C>T]"

    assert VariantNormalizer("SCN5A").normalize_cdna(raw) == raw


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
