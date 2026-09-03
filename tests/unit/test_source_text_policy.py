"""The canonical source-text policy, and what it must refuse to change.

Every case here is character-level. Allele semantics (what a hyphen *means* in
a given HGVS token) are deliberately not this layer's job, and the "must not
fold" cases pin that boundary.

Fixtures are synthetic or anonymous cell payloads. No gene, PMID or expected
record from any dataset appears in production code paths exercised here.
"""

from __future__ import annotations

import pytest

from utils.source_text import (
    PRESERVED_CHARACTERS,
    normalize_header_text,
    normalize_notation_token,
    normalize_source_text,
    split_glued_notation,
)

# Raw cell payloads observed in converted publisher sources, each of which one
# route accepted and another dropped before the policy was shared.
DIVERGENT_CELLS = [
    ("c.1471C >T", "c.1471C>T"),  # thin space
    ("c.722 T>C", "c.722T>C"),  # no-break space
    ("c.2550‑2551insTG", "c.2550-2551insTG"),  # non-breaking hyphen
    ("p.Arg491Trp​", "p.Arg491Trp"),  # zero-width space
    ("c.354T＞G", "c.354T>G"),  # fullwidth greater-than
    ("c.1682C→T", "c.1682C>T"),  # arrow
    ("c.999–424_1338 + 81del", "c.999-424_1338+81del"),  # en dash + spaces
]


@pytest.mark.parametrize("raw,expected", DIVERGENT_CELLS)
def test_divergent_cells_reach_one_canonical_token(raw, expected):
    assert normalize_notation_token(raw) == expected


@pytest.mark.parametrize("raw,expected", DIVERGENT_CELLS)
def test_every_route_agrees_on_the_same_cell(raw, expected):
    """The four routes must not disagree about what a cell says.

    This is the regression that motivated the module: the router returned
    ``None`` for cells the normalizer parsed, so one allele was a row on one
    route and absent on another.
    """
    from pipeline.extraction import ExpertExtractor
    from pipeline.table_router import _normalize_cdna, _normalize_protein

    router = _normalize_cdna(raw) or _normalize_protein(raw)
    assert router == expected, f"table_router disagreed on {raw!r}"

    cleaned = ExpertExtractor._clean_table_cell(None, raw)
    assert cleaned == expected, f"extraction._clean_table_cell disagreed on {raw!r}"


def test_normalization_is_idempotent():
    for raw, _ in DIVERGENT_CELLS:
        once = normalize_source_text(raw)
        assert normalize_source_text(once) == once


# --- What the policy must NOT do -------------------------------------------


@pytest.mark.parametrize(
    "token",
    [
        "c.-15C>T",  # 5' UTR coordinate, not a range
        "c.*70A>G",  # 3' UTR coordinate
        "c.123-2A>G",  # intronic acceptor offset
        "c.123+1G>A",  # intronic donor offset
        "c.169-198_273+820del",  # both endpoints carry their own offset
    ],
)
def test_offset_coordinates_survive_unchanged(token):
    """A hyphen that is an offset must never become a range separator.

    Collapsing the Unicode dash family to ASCII is presentation. Deciding a
    hyphen means ``_`` is notation semantics and stays behind the adjacency
    guard in the normalizer, so these tokens must round-trip untouched.
    """
    assert normalize_notation_token(token) == token


@pytest.mark.parametrize("char,reason", sorted(PRESERVED_CHARACTERS.items()))
def test_meaning_bearing_characters_are_preserved(char, reason):
    """Pinned exclusions: NFKC would fold these and change meaning."""
    payload = f"value{char}suffix"
    assert char in normalize_source_text(payload), reason


def test_two_alleles_glued_by_a_zero_width_char_do_not_fuse():
    """The chimera guard.

    Deleting a zero-width joiner between two alleles produces one token that
    appears in no source and parses as neither. The policy substitutes a space
    so they split, and the token-level entry point refuses the cell outright
    rather than returning a fused string.
    """
    glued = "c.1471C>T​p.Arg491Trp"
    assert split_glued_notation(glued) == ["c.1471C>T", "p.Arg491Trp"]
    assert normalize_notation_token(glued) is None


def test_two_alleles_separated_by_whitespace_do_not_fuse():
    assert normalize_notation_token("c.1471C>T p.Arg491Trp") is None


def test_a_single_token_wrapped_by_the_converter_still_joins():
    """One allele broken across a line is still one allele."""
    assert normalize_notation_token("c.1471C\n>T") == "c.1471C>T"


def test_null_markers_become_none():
    for marker in ("", "  ", "-", "N/A", "na", "none", "nd", "."):
        assert normalize_notation_token(marker) is None


def test_header_matching_ignores_presentational_spacing():
    assert normalize_header_text("No. of patients") == normalize_header_text(
        "No. of patients"
    )
    assert normalize_header_text("Unaffected") == "unaffected"


# --- Regressions found by adversarial review -------------------------------


def test_two_unprefixed_residue_labels_glued_by_zero_width_do_not_fuse():
    """Prefix-based splitting alone was not enough.

    ``normalize_notation_token`` originally refused only when two ``c.``/``p.``
    prefixes were present, so two bare residue labels glued by a zero-width
    joiner had zero prefixes, fell through, and were fused into a token that
    appears in no source.
    """
    glued = "Arg491Trp​Gln403Ter"
    assert split_glued_notation(glued) == ["Arg491Trp", "Gln403Ter"]
    assert normalize_notation_token(glued) is None


def test_one_letter_residue_labels_glued_by_zero_width_do_not_fuse():
    assert normalize_notation_token("R491W​Q403*") is None


def test_the_normalizer_does_not_re_fuse_what_the_policy_split():
    """Every route must agree; the scanner split while the normalizer fused.

    ``_preprocess_cdna_token`` applied the shared policy and then deleted all
    whitespace, undoing the deliberate zero-width split and producing a
    chimeric token the scanner would never emit.
    """
    from utils.variant_normalizer import VariantNormalizer

    fused = VariantNormalizer("GENEA").normalize_cdna("c.1471C>T​p.Arg491Trp")
    assert fused is None or "c.1471C>Tp.Arg491Trp" not in fused


def test_a_wrapped_single_token_is_still_joined_by_the_normalizer():
    """The split guard must not break the legitimate wrapped-token case."""
    from utils.variant_normalizer import VariantNormalizer

    assert VariantNormalizer("GENEA").normalize_cdna("c.1471C >T") == "c.1471C>T"


def test_non_adjacent_dash_range_is_not_converted_to_a_range_separator():
    """A dash between non-adjacent coordinates stays a dash.

    Collapsing the Unicode dash family to ASCII must not let a numeric range
    become an intronic offset or vice versa. Only strictly adjacent
    coordinates are repaired to ``_``, and that rule lives in the normalizer.
    """
    from utils.variant_normalizer import VariantNormalizer

    normalizer = VariantNormalizer("GENEA")
    assert normalizer.normalize_cdna("c.100–200del") == "c.100-200del"
    assert normalizer.normalize_cdna("c.2550‑2551insTG") == "c.2550_2551insTG"
