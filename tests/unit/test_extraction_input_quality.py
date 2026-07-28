"""Circuit-breaker input-quality gate (`ExpertExtractor._assess_input_quality`).

Pins the 2026-07-27 change to Check 3: markup disqualifies an extraction input
by *density*, not by how many distinct markup regexes appear somewhere in the
document.

The old rule discarded any paper matching >=3 of five patterns anywhere in the
text, which threw away real full text. Two cases from the locked 48-paper
cardiac eval set (benchmarks/codex_paper_eval/runs/20260726_fixed48_production):

* RYR2 PMID 19926015 - 62,988 chars carrying 340 chars of markup (0.54%), 3/5
  patterns present. Discarded whole; gold has 40 variants.
* KCNQ1 PMID 25087618 - 35,474 chars carrying 602 chars of markup (1.70%), 3/5
  patterns present. Discarded whole; gold has 1 variant, which a single-model
  protocol reading the same file extracted correctly.

`html_open` and `html_close` match the two halves of the same tag pair, so
presence-counting made the effective threshold 2, not 3: any document with a few
tags plus one `class=` attribute or one brace-colon fragment tripped it.

The fixtures below are self-checking - each asserts its own measured density -
so a future edit to the pattern list fails with the reason rather than silently
drifting off the band it was built to represent.

Offline: no model calls, no network.
"""

import re
import time
from pathlib import Path

import pytest

from config.constants import (
    MAX_MARKUP_DENSITY,
    MIN_ALPHANUMERIC_RATIO,
    MIN_EXTRACTION_INPUT_SIZE,
)
from pipeline.extraction import ExpertExtractor, _find_data_zones_file

PROSE = (
    "Proband II-3 presented with syncope during exercise and a corrected QT "
    "interval of 480 ms. Targeted sequencing identified the heterozygous "
    "variant c.1021G>A (p.Ala341Val), which segregated with the phenotype in "
    "four of six genotyped relatives. Of the 18 mutation carriers in this "
    "kindred, 11 were symptomatic and 7 remained asymptomatic through the "
    "follow-up period. Resting electrocardiograms were obtained for every "
    "family member and reviewed independently by two cardiologists. "
)


def _body(target_chars: int) -> str:
    """Realistic clinical prose padded to roughly ``target_chars``."""
    repeats = max(1, target_chars // len(PROSE))
    return (PROSE * repeats)[:target_chars]


def _interleave(body: str, snippets: list[str]) -> str:
    """Sprinkle ``snippets`` evenly through ``body`` the way a converted page does."""
    if not snippets:
        return body
    stride = len(body) // (len(snippets) + 1)
    out = []
    for index, snippet in enumerate(snippets):
        out.append(body[index * stride : (index + 1) * stride])
        out.append(snippet)
    out.append(body[len(snippets) * stride :])
    return "".join(out)


def _distinct_pattern_count(text: str) -> int:
    """How many of the five patterns are present - i.e. the pre-2026-07-27 score."""
    return sum(
        1
        for pattern in ExpertExtractor.HTML_GARBAGE_PATTERNS
        if re.search(pattern, text)
    )


@pytest.fixture(scope="module")
def extractor() -> ExpertExtractor:
    return ExpertExtractor(models=["gpt-4"])


@pytest.fixture(scope="module")
def ryr2_shaped_source() -> str:
    """~63 KB of prose carrying ~0.5% markup, after RYR2 PMID 19926015.

    Real profile: html_open x10, html_close x10, class_attr x10 -> 3/5 patterns.
    """
    snippets = []
    for index in range(10):
        snippets.append(f'<div class="section-{index}">')
        snippets.append("</div>")
    return _interleave(_body(62_988), snippets)


@pytest.fixture(scope="module")
def kcnq1_shaped_source() -> str:
    """~35 KB of prose carrying ~1.7% markup, after KCNQ1 PMID 25087618.

    Real profile: html_open x14, html_close x14, json_frag x2 -> 3/5 patterns.
    """
    snippets = []
    for index in range(14):
        snippets.append(f'<span class="small-caps" id="para-{index}">')
        snippets.append("</span>")
    snippets.append('{"figure": "Figure 1 legend"}')
    snippets.append('{"table": "Table 2 legend"}')
    return _interleave(_body(35_474), snippets)


class TestTraceMarkupPasses:
    """Papers that are readable prose with incidental markup must not be discarded."""

    def test_ryr2_shaped_source_passes(self, extractor, ryr2_shaped_source):
        density, content_chars = extractor._markup_profile(ryr2_shaped_source)

        assert 0.003 <= density <= 0.008, f"fixture drifted off ~0.5%: {density:.4%}"
        assert content_chars > 50_000
        assert _distinct_pattern_count(ryr2_shaped_source) >= 3, (
            "fixture must reproduce the >=3-distinct-patterns shape the old rule rejected"
        )

        is_usable, reason = extractor._assess_input_quality(ryr2_shaped_source, "RYR2")
        assert is_usable, f"trace markup must not discard a 63 KB paper: {reason}"

    def test_kcnq1_shaped_source_passes(self, extractor, kcnq1_shaped_source):
        density, content_chars = extractor._markup_profile(kcnq1_shaped_source)

        assert 0.014 <= density <= 0.022, f"fixture drifted off ~1.7%: {density:.4%}"
        assert content_chars > 30_000
        assert _distinct_pattern_count(kcnq1_shaped_source) >= 3

        is_usable, reason = extractor._assess_input_quality(
            kcnq1_shaped_source, "KCNQ1"
        )
        assert is_usable, f"trace markup must not discard a 35 KB paper: {reason}"

    def test_open_close_tag_pair_is_not_two_independent_signals(self, extractor):
        """One tag pair plus one attribute used to reach the old threshold of 3."""
        text = _interleave(_body(20_000), ['<div class="article-body">', "</div>"])

        assert _distinct_pattern_count(text) >= 3
        assert extractor._assess_input_quality(text, "KCNQ1")[0]

    def test_brace_colon_prose_is_not_markup(self, extractor):
        """`json_frag` matches ordinary legend prose; it must not cast a deciding vote."""
        legends = [
            "{Abbreviations: LQTS, long QT syndrome; TdP, torsades de pointes}",
            "{Legend: filled symbols denote mutation carriers}",
        ]
        text = _interleave(_body(20_000), legends + ["<span>", "</span>"])

        json_frag = ExpertExtractor.HTML_GARBAGE_PATTERNS[2]
        assert all(re.search(json_frag, legend) for legend in legends), (
            "fixture must actually trip json_frag for this test to mean anything"
        )
        assert _distinct_pattern_count(text) >= 3
        assert extractor._assess_input_quality(text, "KCNQ1")[0]


class TestMarkupDominatedRejected:
    """The breaker must still fire when markup has crowded the prose out."""

    def test_markup_dominated_document_is_rejected(self, extractor):
        nav = (
            '<div class="nav-item"><span style="display:none">'
            '<span class="icon"></span></span></div>'
        )
        text = nav * 40 + (
            " KCNQ1 variant p.Ala341Val was reported in this cohort of carriers. "
        )

        density, content_chars = extractor._markup_profile(text)
        assert density > 0.30, f"fixture is not markup-dominated: {density:.2%}"
        assert content_chars < MIN_EXTRACTION_INPUT_SIZE
        assert len(text) >= MIN_EXTRACTION_INPUT_SIZE, (
            "must clear Check 1 to reach Check 3"
        )

        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert not is_usable
        assert "markup-dominated" in reason.lower(), reason

    def test_inline_script_body_counts_as_markup(self, extractor):
        """A script blob is markup for its whole length, not just its opening tag."""
        script = "<script>" + ("var pageData={id:1,track:true};" * 60) + "</script>"
        text = script + " KCNQ1 c.1021G>A carrier. "

        density, content_chars = extractor._markup_profile(text)
        assert density > 0.90, (
            f"script body must count against the document: {density:.2%}"
        )
        assert content_chars < MIN_EXTRACTION_INPUT_SIZE

        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert not is_usable
        assert "markup-dominated" in reason.lower(), reason

    def test_high_density_with_real_prose_behind_it_still_passes(self, extractor):
        """Density alone must not reject: the corpus's densest real papers reach ~46%.

        PMID 39162406 is 46% markup and still carries 44 KB of prose.
        """
        nav = '<div class="nav-item"><span style="display:none"></span></div>'
        text = nav * 400 + _body(30_000)

        density, content_chars = extractor._markup_profile(text)
        assert density >= MAX_MARKUP_DENSITY
        assert content_chars >= MIN_EXTRACTION_INPUT_SIZE

        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert is_usable, f"paper with 30 KB of prose behind the markup: {reason}"


class TestOtherChecksStillFire:
    """Check 3 got looser; the checks around it did not."""

    def test_paywall_placeholder_still_rejected(self, extractor):
        text = (
            "[PDF file available at: http://example.com/paper.pdf]\n"
            "Subscription required to view this article.\n"
            "Please login to continue.\n" + _body(2_000)
        )

        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert not is_usable
        assert "failed extraction markers" in reason

    def test_failed_conversion_placeholder_still_rejected(self, extractor):
        text = "[Error converting PDF to text]\n[NO TEXT AVAILABLE]\n" + _body(2_000)

        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert not is_usable
        assert "failed extraction markers" in reason

    def test_low_alphanumeric_ratio_still_rejected(self, extractor):
        text = "!@#$%^&*()_+|:<>?[];',./" * 100

        assert len(text) >= MIN_EXTRACTION_INPUT_SIZE
        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert not is_usable
        assert "alphanumeric" in reason.lower(), reason
        assert f"min={MIN_ALPHANUMERIC_RATIO}" in reason

    def test_short_input_still_rejected(self, extractor):
        is_usable, reason = extractor._assess_input_quality(
            "KCNQ1 p.Ala341Val", "KCNQ1"
        )
        assert not is_usable
        assert "too short" in reason.lower()

    def test_irrelevant_long_text_still_rejected(self, extractor):
        text = (
            "This article discusses regional pasta shapes and the water-to-flour "
            "ratios preferred by home cooks in each province. "
        ) * 40

        is_usable, reason = extractor._assess_input_quality(text, "KCNQ1")
        assert not is_usable
        assert "no variant patterns" in reason.lower()


class TestMarkupProfile:
    """The measurement itself: overlapping matches must not be double-counted."""

    def test_nested_attribute_match_counted_once(self, extractor):
        tag = '<div class="x">'
        text = tag + "a" * 985

        density, content_chars = extractor._markup_profile(text)

        # `class="x"` sits inside the `<div ...>` match; merged, the span is the
        # tag's 15 chars, not 15 + 9.
        assert density == pytest.approx(len(tag) / len(text))
        assert content_chars == 985

    def test_plain_prose_has_zero_density(self, extractor):
        text = _body(5_000)

        density, content_chars = extractor._markup_profile(text)

        assert density == 0.0
        assert content_chars > 4_000

    def test_empty_text(self, extractor):
        assert extractor._markup_profile("") == (0.0, 0)

    def test_closing_tag_case_and_spacing_tolerated(self, extractor):
        text = "<SCRIPT type='text/js'>" + "var a=1;" * 200 + "</Script >" + "z" * 200

        density, content_chars = extractor._markup_profile(text)

        assert density > 0.85, "the whole element is markup, not just its open tag"
        assert content_chars == 200

    def test_unclosed_script_opener_counts_only_its_own_tag(self, extractor):
        """A truncated opener must not swallow the prose that follows it."""
        opener = '<script src="/static/analytics.js">'
        text = opener + _body(20_000)

        density, content_chars = extractor._markup_profile(text)

        assert density == pytest.approx(len(opener) / len(text))
        assert content_chars > 19_000
        assert extractor._assess_input_quality(text, "KCNQ1")[0]

    def test_many_unclosed_openers_scan_linearly(self, extractor):
        """Guards `_markup_block_spans` against a return to regex backtracking.

        A non-greedy `.*?</script>` re-scans the tail once per opener that never
        closes: this input took ~7 s that way and ~8 ms with the linear scan.
        """
        text = ('<script data-x="1" ' + "a" * 80 + ">") * 20_000 + "x" * 100_000

        started = time.perf_counter()
        density, _ = extractor._markup_profile(text)
        elapsed = time.perf_counter() - started

        assert density > 0.9
        assert elapsed < 2.0, f"markup profiling went superlinear: {elapsed:.2f}s"


def test_data_zones_discovery_never_scans_implicit_working_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """A stale file in cwd must not be selected without an explicit run path."""
    stale = tmp_path / "12345678_DATA_ZONES.md"
    stale.write_text("stale cross-run artifact", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    assert _find_data_zones_file("12345678") is None
    assert _find_data_zones_file("12345678", [str(tmp_path)]) == stale
