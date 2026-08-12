"""Coverage for per-gold-row source reachability stratification of a tier gate.

The point of the module under test is to stop a tier gate from reporting an
un-downloadable paper as an extraction failure, without handing out credit the
run did not earn. Both halves of that need pinning: the exclusion must fire when
the source genuinely cannot support a row, and it must NOT fire when the probe
is merely unsure.
"""

from __future__ import annotations

import json

import pytest

from scripts.recall_audit.source_stratified_metrics import (
    ABSTRACT_ONLY_MARKER,
    MIN_FULL_TEXT_CHARS,
)
from scripts.recall_audit.tier_source_reachability import (
    BodyIndex,
    FormIndex,
    GoldRow,
    classify_row,
    is_substitution,
    pair_paper,
    render,
    score_tier,
)


def _body(corpus, gene, pmid, text):
    paper = corpus / gene / pmid
    paper.mkdir(parents=True, exist_ok=True)
    (paper / f"{pmid}_FULL_CONTEXT.md").write_text(text, encoding="utf-8")


def _padded(text):
    """A body long enough to clear the stub threshold."""
    return text + ("\nfiller line to clear the stub threshold." * 400)


def _classify(corpus, gene, pmid, variant):
    row = GoldRow(gene=gene, pmid=pmid, variant=variant, values={})
    return classify_row(row, BodyIndex(corpus), FormIndex())[0]


@pytest.mark.parametrize(
    "variant",
    ["R176W", "p.Arg176Trp", "Arg176Trp", "R176*", "R176X", "p.R176W"],
)
def test_substitutions_are_recognised(variant):
    assert is_substitution(variant)


@pytest.mark.parametrize(
    "variant",
    [
        # The regression this predicate exists for: a three-letter tail that is
        # not an amino acid. A bare length check accepted these and excluded
        # rows the run had actually extracted.
        "P.F1617DEL",
        "P.N1380del",
        "p.Q1507_P1509del",
        "c.526C>T",
        "c.4389_4396delCCTCTTTA",
        "R920fsX",
        "IVS7-2A>G",
        "",
    ],
)
def test_non_substitutions_are_rejected(variant):
    assert not is_substitution(variant)


def test_missing_source_is_body_absent(tmp_path):
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "body_absent"


def test_abstract_only_fallback_is_body_absent(tmp_path):
    _body(tmp_path, "KCNH2", "1", _padded(f"{ABSTRACT_ONLY_MARKER}\nR176W carriers"))
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "body_absent"


def test_stub_shorter_than_a_body_is_body_absent(tmp_path):
    _body(tmp_path, "KCNH2", "1", "R176W")
    assert len("R176W") < MIN_FULL_TEXT_CHARS
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "body_absent"


def test_variant_in_body_is_penalized(tmp_path):
    _body(tmp_path, "KCNH2", "1", _padded("Table 2 lists R176W in 16 carriers."))
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_variant_found_through_a_normalized_form(tmp_path):
    """Gold writes one notation, the paper another. That is not unreachable."""
    _body(tmp_path, "KCNH2", "1", _padded("the p.Arg176Trp proband"))
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_whitespace_in_the_body_does_not_hide_a_variant(tmp_path):
    _body(tmp_path, "KCNH2", "1", _padded("p. Arg176\nTrp"))
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_absent_substitution_is_excluded(tmp_path):
    _body(tmp_path, "KCNH2", "1", _padded("this paper discusses only A561V"))
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "variant_absent_from_body"


def test_absent_indel_is_inconclusive_not_excluded(tmp_path):
    """The conservative half: an unsearchable notation stays penalized.

    A body may write ``c.4389_4396delCCTCTTTA`` where gold writes
    ``P.N1380del``. Absence is not evidence, so the row must not be excluded --
    a weak probe must never manufacture credit.
    """
    _body(tmp_path, "SCN5A", "1", _padded("c.4389_4396delCCTCTTTA in two probands"))
    assert _classify(tmp_path, "SCN5A", "1", "P.N1380del") == "notation_inconclusive"


def test_a_variant_only_in_the_cleaned_file_is_reachable(tmp_path):
    """The extractor may read FULL_CONTEXT; reachability asks about all of disk.

    On Tier 1 this recovered 46 of 246 candidate exclusions whose variant was
    sitting in a sibling file on the same disk.
    """
    paper = tmp_path / "KCNH2" / "1"
    paper.mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(_padded("no variants here"))
    (paper / "1_CLEANED.md").write_text(_padded("Table 2 lists R176W"))
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_a_usable_cleaned_file_does_not_require_full_context(tmp_path):
    paper = tmp_path / "KCNH2" / "1"
    paper.mkdir(parents=True)
    (paper / "1_CLEANED.md").write_text(
        _padded("Table 2 lists R176W"), encoding="utf-8"
    )

    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_a_variant_only_in_a_text_supplement_is_reachable(tmp_path):
    paper = tmp_path / "KCNH2" / "1"
    (paper / "1_supplements").mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(_padded("see supplementary table"))
    (paper / "1_supplements" / "mmc1.csv").write_text("variant,n\nR176W,16\n")
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_a_variant_only_in_reader_pdf_text_is_reachable(tmp_path, monkeypatch):
    paper = tmp_path / "KCNH2" / "1"
    paper.mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(_padded("only A561V in the body"))
    (paper / "supplement.pdf").write_bytes(b"%PDF-placeholder")
    monkeypatch.setattr(
        "scripts.recall_audit.tier_source_reachability.extract_pdf_text",
        lambda paths, max_chars: "PDF table lists p.Arg176Trp in 16 carriers",
    )

    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_pdf_can_recover_a_row_when_body_is_abstract_only(tmp_path, monkeypatch):
    paper = tmp_path / "KCNH2" / "1"
    paper.mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(
        _padded(f"{ABSTRACT_ONLY_MARKER}\nonly A561V"), encoding="utf-8"
    )
    (paper / "article.pdf").write_bytes(b"%PDF-placeholder")
    monkeypatch.setattr(
        "scripts.recall_audit.tier_source_reachability.extract_pdf_text",
        lambda paths, max_chars: "PDF table lists p.Arg176Trp in 16 carriers",
    )

    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "present_in_body"


def test_an_unsearchable_binary_supplement_does_not_grant_reachability(tmp_path):
    paper = tmp_path / "KCNH2" / "1"
    (paper / "1_supplements").mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(_padded("see supplementary table"))
    (paper / "1_supplements" / "mmc1.xlsx").write_bytes(b"PK\x03\x04binary")
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "variant_absent_from_body"


def test_figures_on_disk_block_exclusion(tmp_path):
    """A PNG we cannot grep is still source we acquired.

    The figure lane reads images, so absence from the text is not evidence the
    row was unreachable. Excluding it would credit the run for rows the figure
    reader actually had a shot at -- 18 such rows exist in Tier 1.
    """
    paper = tmp_path / "KCNH2" / "1"
    (paper / "1_figures").mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(_padded("only A561V in the text"))
    (paper / "1_figures" / "fig1.png").write_bytes(b"\x89PNG")
    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "absent_but_figures_present"


def test_figures_block_exclusion_even_without_a_body(tmp_path):
    paper = tmp_path / "KCNH2" / "1" / "figures"
    paper.mkdir(parents=True)
    (paper / "fig1.png").write_bytes(b"\x89PNG")

    assert _classify(tmp_path, "KCNH2", "1", "R176W") == "absent_but_figures_present"


def test_figure_blocked_rows_stay_in_the_scored_population(tmp_path):
    corpus = tmp_path / "corpus"
    paper = corpus / "KCNH2" / "1"
    (paper / "1_figures").mkdir(parents=True)
    (paper / "1_FULL_CONTEXT.md").write_text(_padded("R176W only"))
    (paper / "1_figures" / "fig1.png").write_bytes(b"\x89PNG")
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n", "A561V,1,7,7,0\n"])
    run = _run_dir(
        tmp_path,
        [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
    )

    payload = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)
    assert payload["reachability_classes"]["absent_but_figures_present"] == 1
    assert payload["noted_exclusions"] == []
    assert payload["strata"]["source_reachable"]["gold_rows"] == 2


def test_pairing_is_greedy_and_one_to_one():
    """Two predictions of the same variant cannot both claim the one gold row."""
    gold = [{"variant": "R176W", "carriers": 3}]
    predicted = [
        {"variant": "R176W", "carriers": 3},
        {"variant": "R176W", "carriers": 99},
    ]
    pairs, missed = pair_paper("KCNH2", predicted, gold)
    assert len(pairs) == 1
    assert missed == []
    assert pairs[0][0]["carriers"] == 3


def test_unmatched_gold_rows_are_reported_as_missed():
    gold = [{"variant": "R176W", "carriers": 3}, {"variant": "A561V", "carriers": 1}]
    pairs, missed = pair_paper("KCNH2", [{"variant": "R176W", "carriers": 3}], gold)
    assert len(pairs) == 1
    assert missed == [1]


def _run_dir(tmp_path, variants, published=None):
    run = tmp_path / "run"
    run.mkdir(parents=True, exist_ok=True)
    (run / "predictions.json").write_text(
        json.dumps({"papers": [{"gene": "KCNH2", "pmid": "1", "variants": variants}]}),
        encoding="utf-8",
    )
    (run / "report.json").write_text(
        json.dumps({"overall": published or {}}), encoding="utf-8"
    )
    return run


def _gold_dir(tmp_path, rows):
    gold = tmp_path / "gold"
    gold.mkdir()
    header = "variant,pmid,carriers,affected,unaffected\n"
    (gold / "KCNH2_recall_input.csv").write_text(
        header + "".join(rows), encoding="utf-8"
    )
    return gold


def test_unreachable_rows_leave_the_scored_stratum_and_are_noted(tmp_path):
    corpus = tmp_path / "corpus"
    # The body holds R176W only. A561V is a substitution the body does not
    # contain, so no reader could have produced it.
    _body(corpus, "KCNH2", "1", _padded("Table 2: R176W, 16 carriers, 4 affected"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n", "A561V,1,7,7,0\n"])
    run = _run_dir(
        tmp_path,
        [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
    )

    payload = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)

    assert payload["reachability_classes"] == {
        "present_in_body": 1,
        "variant_absent_from_body": 1,
    }
    all_gold = payload["strata"]["all_gold"]
    reachable = payload["strata"]["source_reachable"]
    excluded = payload["strata"]["excluded_unreachable"]

    assert all_gold["gold_rows"] == 2
    assert reachable["gold_rows"] == 1
    assert excluded["gold_rows"] == 1
    # The whole point: the unreachable row stops depressing count recall.
    assert all_gold["counts"]["carriers"]["recall"] == pytest.approx(0.5)
    assert reachable["counts"]["carriers"]["recall"] == pytest.approx(1.0)
    # ...and it is surfaced rather than silently dropped.
    noted = payload["noted_exclusions"]
    assert [(n["variant"], n["class"]) for n in noted] == [
        ("A561V", "variant_absent_from_body")
    ]
    assert noted[0]["matched_anyway"] is False


def test_inconclusive_notation_still_counts_against_the_score(tmp_path):
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("Table 2: R176W, 16 carriers"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n", "P.N1380del,1,7,7,0\n"])
    run = _run_dir(
        tmp_path,
        [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
    )

    payload = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)

    assert payload["reachability_classes"]["notation_inconclusive"] == 1
    assert payload["noted_exclusions"] == []
    # Both gold rows remain in the scored population, so the missing count for
    # the indel row is still a miss.
    assert payload["strata"]["source_reachable"]["gold_rows"] == 2
    assert payload["strata"]["source_reachable"]["counts"]["carriers"][
        "recall"
    ] == pytest.approx(0.5)


def test_a_paper_with_no_body_at_all_is_fully_excluded(tmp_path):
    """The EZproxy-expired case: nothing on disk, so nothing is charged."""
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n"])
    run = _run_dir(tmp_path, [])

    payload = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)

    assert payload["reachability_classes"] == {"body_absent": 1}
    assert payload["strata"]["source_reachable"]["gold_rows"] == 0
    assert payload["strata"]["excluded_unreachable"]["gold_rows"] == 1
    assert payload["noted_exclusions"][0]["reason"] == "no source on disk"


def test_error_per_gold_row_penalizes_a_null_count(tmp_path):
    """MAE over predictions alone rewards silence; this form does not.

    Two gold rows, one predicted exactly right and one not predicted at all.
    MAE sees a perfect 0.0 because it only averages over predictions;
    error_per_gold_row must not.
    """
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("R176W and A561V both appear here"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n", "A561V,1,7,7,0\n"])
    run = _run_dir(
        tmp_path,
        [
            {"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12},
            {"variant": "A561V", "carriers": None},
        ],
    )

    counts = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)["strata"][
        "source_reachable"
    ]["counts"]["carriers"]

    assert counts["mae"] == pytest.approx(0.0)
    assert counts["recall"] == pytest.approx(0.5)
    assert counts["error_per_gold_row"] == pytest.approx(3.5)
    assert counts["missing_predictions"] == 1

    # And a wrong count is charged over the same denominator.
    run2 = _run_dir(
        tmp_path / "b",
        [
            {"variant": "R176W", "carriers": 20, "affected": 4, "unaffected": 12},
            {"variant": "A561V", "carriers": None},
        ],
    )
    counts2 = score_tier([("KCNH2", "1")], run2, corpus, gold_dir=gold)["strata"][
        "source_reachable"
    ]["counts"]["carriers"]
    assert counts2["error_per_gold_row"] == pytest.approx(5.5)


def test_missing_explicit_zero_still_has_a_one_unit_coverage_loss(tmp_path):
    """A zero-valued key cannot disappear from the coverage-aware metric."""
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("A561V appears here"))
    gold = _gold_dir(tmp_path, ["A561V,1,7,7,0\n"])
    run = _run_dir(tmp_path, [{"variant": "A561V", "unaffected": None}])

    counts = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)["strata"][
        "source_reachable"
    ]["counts"]["unaffected"]

    assert counts["mae"] is None
    assert counts["recall"] == pytest.approx(0.0)
    assert counts["error_per_gold_row"] == pytest.approx(1.0)
    assert counts["missing_predictions"] == 1


def test_pairing_mismatch_is_a_hard_gate(tmp_path):
    """A drifted matcher must withhold, not publish."""
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("R176W here"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n"])
    run = _run_dir(
        tmp_path,
        [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
        published={"tp": 99, "count": {}},
    )

    payload = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)
    assert payload["reproduction"]["pairing_ok"] is False


def test_answer_key_drift_is_reported_without_failing(tmp_path):
    """A gold value adjudicated after the run is drift, not a scoring defect."""
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("R176W here"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n"])
    run = _run_dir(
        tmp_path,
        [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
        published={
            "tp": 1,
            "count": {
                # As-scored the key said 99 carriers; today it says 16.
                "carriers": {"predicted": 1, "gold_asserted": 1, "mae": 83.0},
                "affected": {"predicted": 1, "gold_asserted": 1, "mae": 0.0},
                "unaffected": {"predicted": 1, "gold_asserted": 1, "mae": 0.0},
            },
        },
    )

    rep = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)["reproduction"]
    assert rep["pairing_ok"] is True
    assert [d["metric"] for d in rep["answer_key_drift"]] == ["carriers.mae"]


def test_a_run_that_does_not_cover_the_tier_is_partial_not_drifted(tmp_path):
    """Scoring tier 2 against a tier-1 run is partial coverage, not matcher drift.

    Conflating the two is dangerous in both directions: it raises a false alarm
    about the matcher, and it invites quoting a subset as the tier result.
    """
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("R176W here"))
    _body(corpus, "KCNH2", "2", _padded("A561V here"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n", "A561V,2,7,7,0\n"])
    # The run holds paper 1 only, while the tier asks for 1 and 2.
    run = _run_dir(
        tmp_path,
        [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
        published={"tp": 1, "count": {}},
    )

    payload = score_tier([("KCNH2", "1"), ("KCNH2", "2")], run, corpus, gold_dir=gold)

    rep = payload["reproduction"]
    assert rep["applicable"] is False
    assert rep["pairing_ok"] is None
    assert payload["coverage"]["run_matches_tier"] is False
    assert payload["coverage"]["scored_papers"] == 1
    assert payload["coverage"]["tier_attempts_with_gold"] == 2
    # Reachability still covers the whole tier, so the exclusion list is usable
    # before a run exists for every paper.
    assert sum(payload["reachability_classes"].values()) == 2
    assert "PARTIAL SUBSET" in render("t", payload)


def test_a_run_holding_papers_outside_the_tier_is_partial(tmp_path):
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("R176W here"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n", "A561V,2,7,7,0\n"])
    run = tmp_path / "run"
    run.mkdir()
    (run / "predictions.json").write_text(
        json.dumps(
            {
                "papers": [
                    {"gene": "KCNH2", "pmid": "1", "variants": []},
                    {"gene": "KCNH2", "pmid": "2", "variants": []},
                ]
            }
        ),
        encoding="utf-8",
    )
    (run / "report.json").write_text(json.dumps({"overall": {}}), encoding="utf-8")

    payload = score_tier([("KCNH2", "1")], run, corpus, gold_dir=gold)
    assert payload["coverage"]["run_papers_outside_tier"] == ["KCNH2:2"]
    assert payload["reproduction"]["applicable"] is False


def test_classification_never_reads_a_prediction(tmp_path):
    """Integrity: the same gold row classifies identically under any prediction.

    If the exclusion predicate could see predictions it could be tuned to drop
    whatever scored badly. Same corpus and gold, opposite predictions, identical
    reachability.
    """
    corpus = tmp_path / "corpus"
    _body(corpus, "KCNH2", "1", _padded("only A561V is discussed"))
    gold = _gold_dir(tmp_path, ["R176W,1,16,4,12\n"])

    empty = score_tier(
        [("KCNH2", "1")], _run_dir(tmp_path / "a", []), corpus, gold_dir=gold
    )
    perfect = score_tier(
        [("KCNH2", "1")],
        _run_dir(
            tmp_path / "c",
            [{"variant": "R176W", "carriers": 16, "affected": 4, "unaffected": 12}],
        ),
        corpus,
        gold_dir=gold,
    )
    assert (empty["reachability_classes"] == perfect["reachability_classes"]) and empty[
        "reachability_classes"
    ] == {"variant_absent_from_body": 1}
