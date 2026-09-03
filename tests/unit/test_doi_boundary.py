"""DOI extraction must stop at a glued-prose boundary.

Four modules each carried their own copy of the same greedy pattern
(``10\\.\\d{4,9}/[^\\s"'<>]+``). That is correct for a DOI in a URL or a
metadata field and wrong for one harvested from *rendered page text*, where a
publisher landing page can print the DOI immediately followed by the next word
with no delimiter. The extracted identifier then cannot resolve, and the
failure surfaces as a fetch-stage error rather than as a parse error.

The repair is a boundary rule, not a publisher rule, so the "must not change"
cases below matter more than the repaired one: a DOI cleaner that trims real
identifiers is worse than one that occasionally over-collects.
"""

from __future__ import annotations

import pytest

from utils.doi import clean_doi, doi_from_text, extract_doi

# Real DOIs spanning the publishers this corpus actually draws from. None may
# be altered by cleaning. Several legitimately end in digits, which is the
# character class the glue rule keys on, so they are the sharp test.
REAL_DOIS = [
    "10.1038/s41586-021-03819-2",
    "10.1016/j.jacc.2019.01.011",
    "10.1016/j.stem.2016.08.019",
    "10.1002/humu.21260",
    "10.1161/CIRCGENETICS.113.000260",
    "10.1093/nar/gkw1121",
    "10.1186/s13073-019-0685-z",
    "10.1186/s12881-018-0608-7",
    "10.1186/1465-9921-12-99",
    "10.1186/1755-8794-1-45",
    "10.1183/13993003.00327-2021",
    "10.5144/0256-4947.2026.111",
    "10.3390/ijms25052734",
]


@pytest.mark.parametrize("doi", REAL_DOIS)
def test_real_dois_are_never_altered(doi):
    assert clean_doi(doi) == doi


@pytest.mark.parametrize("doi", REAL_DOIS)
def test_real_dois_survive_a_resolver_prefix(doi):
    assert clean_doi(f"https://doi.org/{doi}") == doi
    assert clean_doi(f"http://dx.doi.org/{doi}") == doi


@pytest.mark.parametrize("doi", REAL_DOIS)
def test_cleaning_is_idempotent(doi):
    once = clean_doi(doi)
    assert clean_doi(once) == once


def test_prose_glued_onto_a_rendered_doi_is_trimmed():
    """The observed production failure.

    A rendered landing page printed the DOI with the next word attached, and
    the source-recovery stage then failed against an identifier that cannot
    exist.
    """
    rendered = (
        "https://doi.org/10.3390/ijms25052734Submission received: 6 December 2023"
    )
    assert extract_doi(rendered) == "10.3390/ijms25052734"


@pytest.mark.parametrize(
    "rendered,expected",
    [
        ("10.1234/abc123Received", "10.1234/abc123"),
        ("10.1234/abc123PublishedOnline", "10.1234/abc123"),
        ("10.1234/abc123Abstract", "10.1234/abc123"),
        ("https://doi.org/10.1234/xyz999Available online", "10.1234/xyz999"),
    ],
)
def test_other_glued_word_shapes_are_trimmed_without_an_allowlist(rendered, expected):
    """One boundary rule, not a list of metadata words.

    This codebase has already recorded that ignore-lists of this shape do not
    converge, so the rule is structural: capitalised words glued directly onto
    a digit at the end of the token are prose.
    """
    assert extract_doi(rendered) == expected


def test_uppercase_suffix_segments_are_not_mistaken_for_glue():
    """An uppercase run that is part of the identifier must survive."""
    assert clean_doi("10.1161/CIRCEP.111.962019") == "10.1161/CIRCEP.111.962019"
    assert clean_doi("10.1234/ABC") == "10.1234/ABC"


def test_trimming_never_empties_the_suffix():
    """Over-trimming to a bare prefix would be worse than not trimming."""
    cleaned = clean_doi("10.1234/Submission")
    assert cleaned.startswith("10.1234/")
    assert cleaned.split("/", 1)[1]


def test_trailing_punctuation_is_stripped():
    assert extract_doi("see 10.1002/humu.21260.") == "10.1002/humu.21260"
    assert extract_doi("(10.1002/humu.21260)") == "10.1002/humu.21260"


def test_labelled_doi_wins_over_a_reference_list_doi():
    """The citation line names the article; a reference's DOI is not the paper's."""
    text = "References: 10.9999/other.1 ... Citation doi: 10.1234/mine.5"
    assert doi_from_text(text) == "10.1234/mine.5"


def test_empty_and_absent_inputs_are_empty_strings():
    for value in (None, "", "   ", "no identifier here"):
        assert extract_doi(value) == ""


def test_every_module_shares_one_cleaner():
    """The defect was four copies drifting apart, so pin the delegation."""
    from harvesting.orchestrator import _clean_doi as orchestrator_clean
    from scripts.recall_audit.build_acquisition_worklist import (
        _clean_doi as worklist_clean,
    )
    from scripts.recall_audit.source_acquisition_audit import (
        _clean_doi as audit_clean,
    )

    glued = "10.3390/ijms25052734Submission"
    expected = "10.3390/ijms25052734"
    assert orchestrator_clean(glued) == expected
    assert worklist_clean(glued) == expected
    assert audit_clean(glued) == expected
