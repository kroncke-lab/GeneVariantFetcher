"""Data scout hardening: pathological content must not wedge a run.

RYR2 PMID 32508047's FULL_CONTEXT carried a 14.5 MB supplement whose .doc
conversion leaked raw OLE bytes as text. Stray pipe characters in the soup
fabricated ~18K "markdown tables", and the scout's per-candidate bookkeeping
(prefix-length sums, whole-prefix parent-section scans, uncapped caption
walk-backs) went quadratic: the 2026-08-29 reviewer-50 refresh sat 2h45m+ at
100% CPU inside the scan. These tests pin the linear-time rewrite, the input
bounds, and the wall-clock budget net.
"""

import time

import pytest

from pipeline.data_scout import GeneticDataScout
from utils.scout_models import DataZoneReport


def _pipe_soup(target_chars: int) -> str:
    """Binary-garbage-with-pipes shaped like the 32508047 leaked .doc.

    Control-char noise interleaved with pipe rows and |---| separators so the
    markdown-table detector fires constantly, the way the real soup did.
    """
    unit = (
        "\x00\x01g\x00\x00\x0bg\x02\x03\x04\x05\x06\x07\x08\x0e\x0f\x10\x11\n"
        "| g\x00\x00\x08g\x00\x00 | g\x00\x00 |\n"
        "|---|---|---|---|\n"
        "g\x00\x00\x0bg\x00\x00\ng\x00\x00\n"
        "| \x13\x14\x15\x16\x17\x18\x19\x1a | îˇ\x00\x08÷r\x00\x05Òˇí |\n"
        "|---|---|\n"
        "@\x0fÓ\x13Æ\x17©\x1b\x00#\x00\x00ˇˇˇˇˇˇˇˇ\x00#\x00\x00ˇˇˇˇˇˇˇˇ\n"
    )
    return unit * (target_chars // len(unit) + 1)


def test_pipe_soup_scan_is_linear_not_quadratic():
    # Pre-fix, 2 MB of this shape took ~150s and grew quadratically (24s at
    # 1 MB); post-fix it completes in a couple of seconds. The 60s budget
    # turns any reintroduced quadratic into a clear truncation failure
    # instead of a hung suite.
    text = _pipe_soup(2_000_000)
    scout = GeneticDataScout("RYR2", scan_budget_seconds=60.0)
    started = time.monotonic()
    report = scout.scan(text, pmid="soup")
    elapsed = time.monotonic() - started
    assert not report.scan_truncated, (
        f"scan hit its budget on 2 MB of pipe soup ({elapsed:.1f}s): "
        "the quadratic candidate bookkeeping is back"
    )


def test_multi_megabyte_single_line_completes():
    # One physical line of uppercase letters and digits with no hyphen: the
    # unbounded gene-mutation pattern ([A-Z][A-Z0-9]{2,}-...) re-counted the
    # whole run at every start position (O(n^2) — hours at this size). The
    # per-line input cap plus the bounded repeat make it instant.
    line = "\t" + "R2D2C3PQ89" * 300_000  # ~3 MB, tab makes it pseudo-table bait
    scout = GeneticDataScout("RYR2", scan_budget_seconds=60.0)
    started = time.monotonic()
    report = scout.scan(line, pmid="monster-line")
    elapsed = time.monotonic() - started
    assert not report.scan_truncated, (
        f"scan hit its budget on a single 3 MB line ({elapsed:.1f}s)"
    )


def test_budget_exhaustion_returns_truncated_empty_report():
    text = _pipe_soup(300_000)
    scout = GeneticDataScout("RYR2", scan_budget_seconds=1e-6)
    report = scout.scan(text, pmid="over-budget")
    assert report.scan_truncated is True
    assert report.zones == []
    assert report.zones_kept == 0
    assert report.total_chars_original == len(text)
    # The deadline must not leak into a later scan on the same instance.
    scout.scan_budget_seconds = 60.0
    again = scout.scan("## Results\n\nnothing here\n", pmid="fresh")
    assert again.scan_truncated is False


def test_zero_budget_disables_the_deadline():
    scout = GeneticDataScout("RYR2", scan_budget_seconds=0)
    report = scout.scan(_pipe_soup(100_000), pmid="uncapped")
    assert report.scan_truncated is False


def test_normal_paper_zones_survive_the_rewrite():
    # The offsets/window refactor must not move zone boundaries: the
    # condensed markdown has to carry the actual table rows and the case
    # paragraph, with the parent section still attributed.
    filler = (
        "The assay protocol follows the manufacturer instructions without "
        "modification and buffers were prepared fresh for each batch run.\n"
    )
    text = (
        "# MAIN TEXT\n\n"
        "## Results\n\n"
        "Table 1: RYR2 variants in probands\n"
        "| Variant | Carriers | Affected |\n"
        "|---|---|---|\n"
        "| p.Arg420Trp | 5 | 3 |\n"
        "| c.1259G>A | 2 | 2 |\n"
        "\n" + filler + "\n"
        # >100 signal-free chars so the table and patient zones don't merge
        # and each keeps its own metadata.
        "Patient 1 was a 34-year-old male proband diagnosed with CPVT. "
        "Patient 2, an unaffected carrier of p.Arg420Trp in RYR2, was "
        "asymptomatic at evaluation.\n"
    )
    scout = GeneticDataScout("RYR2", min_relevance_score=0.1)
    report = scout.scan(text, pmid="fixture")

    assert report.scan_truncated is False
    table_zones = [z for z in report.zones if z.type == "TABLE"]
    assert table_zones, "markdown table zone lost"
    # The merged table zone inherits its metadata from the earliest candidate
    # (pre-existing behavior), so assert section attribution on the patient
    # TEXT zone — its header sits well inside the new 10K search window.
    assert any(z.type == "TEXT" and z.source_section == "results" for z in report.zones)

    condensed = scout.format_markdown(report, text)
    assert "p.Arg420Trp | 5 | 3" in condensed
    assert "Patient 1" in condensed


def test_line_offsets_match_the_old_prefix_sums():
    lines = "a\n\nbb ccc\n| x |\nlast".split("\n")
    offsets = GeneticDataScout._line_offsets(lines)
    for k in range(len(lines) + 1):
        assert offsets[k] == sum(len(lines[j]) + 1 for j in range(k))


def test_find_parent_section_window_preserves_distance_semantics():
    scout = GeneticDataScout("RYR2")
    text = "## Results\n" + ("x" * 20_000)
    near = len("## Results\n") + 5_000
    far = len("## Results\n") + 15_000
    assert scout._find_parent_section(text, near) == "results"
    assert scout._find_parent_section(text, far) is None


def test_scan_truncated_report_serializes_with_flag():
    report = DataZoneReport(
        gene_symbol="RYR2",
        total_zones_found=0,
        zones_kept=0,
        zones_discarded=0,
        total_chars_original=10,
        total_chars_condensed=0,
        compression_ratio=0.0,
        scan_truncated=True,
    )
    assert '"scan_truncated":true' in report.model_dump_json().replace(" ", "")


def test_steps_runner_skips_data_zones_for_truncated_scan(tmp_path, monkeypatch):
    # The pipeline step must not write a DATA_ZONES file from a truncated
    # report — a partial condensation silently hides the unscanned tail.
    import pipeline.data_scout as data_scout_module
    from pipeline.steps import run_data_scout

    class _BudgetlessScout:
        def __init__(self, *args, **kwargs):
            pass

        def scan(self, text, pmid=None):
            return DataZoneReport(
                pmid=pmid,
                gene_symbol="RYR2",
                total_zones_found=0,
                zones_kept=0,
                zones_discarded=0,
                total_chars_original=len(text),
                total_chars_condensed=0,
                compression_ratio=0.0,
                scan_truncated=True,
            )

    monkeypatch.setattr(data_scout_module, "GeneticDataScout", _BudgetlessScout)
    (tmp_path / "111_FULL_CONTEXT.md").write_text("some text", encoding="utf-8")

    result = run_data_scout(tmp_path, "RYR2")

    assert result.success is True
    assert result.stats["errors"] == 1
    assert not (tmp_path / "111_DATA_ZONES.md").exists()


@pytest.mark.parametrize(
    "pattern_text,should_match",
    [
        ("KCNH2-T895M", True),
        ("SCN5A-G1743R", True),
        ("II-1 and III-2 were probands", True),
        ("age at onset: 45", True),
        ("mean age: 41", True),
        ("1234 G>A", True),
    ],
)
def test_bounded_regexes_still_match_real_notation(pattern_text, should_match):
    scout = GeneticDataScout("KCNH2")
    all_patterns = (
        scout._variant_patterns + scout._individual_patterns + scout._aggregate_patterns
    )
    matched = any(p.search(pattern_text) for p in all_patterns)
    assert matched is should_match
