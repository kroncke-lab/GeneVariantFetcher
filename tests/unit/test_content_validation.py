"""Tests for harvesting content validation helpers."""

from harvesting.content_validation import (
    is_abstract_reference_shell,
    validate_content_quality,
)


def test_rejects_empty_content():
    valid, reason = validate_content_quality("")
    assert not valid
    assert "Empty content" in reason


def test_rejects_known_junk_pattern():
    content = "We've detected you are running an older version of your browser."
    valid, reason = validate_content_quality(content + " " * 2000)
    assert not valid
    assert "Junk content detected" in reason


def test_rejects_known_junk_domain():
    content = "abstract methods results discussion references " * 80
    valid, reason = validate_content_quality(
        content, source_url="https://antibodies.cancer.gov/some-page"
    )
    assert not valid
    assert "non-article domain" in reason


def test_accepts_structured_paper_content():
    content = (
        "abstract introduction methods results discussion references "
        "mutation variant c.1234A>G p.Ala123Val"
    ) * 30
    valid, reason = validate_content_quality(content)
    assert valid
    assert reason == "Valid content"


def test_large_abstract_reference_shell_is_not_full_text():
    content = (
        "# Paper title\n\n## Abstract\n"
        + "Methods and Results: aggregate cohort findings. " * 80
        + "\n\n## References\n"
        + "Reference citation with mutation keyword. " * 100
    )

    assert is_abstract_reference_shell(content) is True
    valid, reason = validate_content_quality(content)
    assert valid is False
    assert "Abstract/reference shell" in reason


def test_case_presentation_between_abstract_and_references_is_body():
    content = (
        "## Abstract\nA short summary.\n"
        "## Case presentation\nA patient carrying p.Arg12Trp was evaluated.\n"
        "## References\nA citation.\n"
    )

    assert is_abstract_reference_shell(content) is False
