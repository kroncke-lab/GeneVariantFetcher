"""Pin the gold-row source-presence classifier.

The sweep decides, blind to predictions, whether each gold variant string is
present in *any* on-disk representation of its paper. Its classes feed the
acquisition-ceiling stratum in the mixed-gold reports, so the order of checks
and the two probe blind spots closed on 2026-09-03 (nonsense ``Ter``/``Stop``
spellings; abstract-length stubs and glyph-garbled PDF text) must not drift.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.recall_audit import gold_source_presence_sweep as sweep


class _NoConverter:
    """Binary supplements are not exercised here; fail loudly if reached."""

    def __getattr__(self, name):  # pragma: no cover - guard
        raise AssertionError(f"converter.{name} should not be called")


def _paper(corpus: Path, gene: str, pmid: str, body: str | None, **files: str) -> Path:
    paper = corpus / gene / pmid
    paper.mkdir(parents=True)
    if body is not None:
        (paper / f"{pmid}_FULL_CONTEXT.md").write_text(body)
    for rel, text in files.items():
        target = paper / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
    return paper


def _index(tmp_path: Path, gene: str, pmid: str) -> sweep.PaperIndex:
    return sweep.build_index(
        gene, pmid, tmp_path / "corpus", tmp_path / "cache", _NoConverter()
    )


def _long_body(core: str) -> str:
    filler = "Long QT syndrome is an inherited arrhythmia disorder. " * 200
    return f"# MAIN TEXT\n\n{filler}\n\n{core}\n\n{filler}"


def test_present_in_body_wins_over_everything(tmp_path):
    _paper(
        tmp_path / "corpus",
        "KCNQ1",
        "1",
        _long_body("We found p.Arg518* in two probands."),
    )
    idx = _index(tmp_path, "KCNQ1", "1")
    forms = sweep.FormIndex().get("KCNQ1", "R518X")
    assert sweep.classify("R518X", forms, idx) == ("present_in_body", "body")


@pytest.mark.parametrize("spelling", ["Arg222Ter", "Arg222Stop", "R222*", "p.Arg222X"])
def test_nonsense_spellings_the_scorer_bridges_are_visible_to_the_probe(
    tmp_path, spelling
):
    _paper(
        tmp_path / "corpus",
        "SCN5A",
        "2",
        _long_body(f"carrying the mutation {spelling}."),
    )
    idx = _index(tmp_path, "SCN5A", "2")
    forms = sweep.FormIndex().get("SCN5A", "R222X")
    assert sweep.classify("R222X", forms, idx)[0] == "present_in_body"


def test_text_supplement_only_is_its_own_class_and_names_the_file(tmp_path):
    _paper(
        tmp_path / "corpus",
        "KCNH2",
        "3",
        _long_body("Mutations are listed in Supplementary Table 1."),
        **{"3_supplements/table_s1.csv": "variant,n\nG628S,4\n"},
    )
    idx = _index(tmp_path, "KCNH2", "3")
    forms = sweep.FormIndex().get("KCNH2", "G628S")
    klass, where = sweep.classify("G628S", forms, idx)
    assert klass == "present_in_supplement_only"
    assert where == "supplement:3_supplements/table_s1.csv"


def test_stub_body_is_an_acquisition_class_not_an_extraction_miss(tmp_path):
    _paper(
        tmp_path / "corpus",
        "SCN5A",
        "4",
        "# MAIN TEXT\n\n### Abstract\n\nShort abstract only.\n",
    )
    idx = _index(tmp_path, "SCN5A", "4")
    assert idx.searchable_chars < sweep.STUB_CHARS
    forms = sweep.FormIndex().get("SCN5A", "E1784K")
    klass, why = sweep.classify("E1784K", forms, idx)
    assert klass == "text_absent_stub_body"
    assert "searchable chars" in why


def test_glyph_garbled_pdf_text_is_flagged_even_when_long(tmp_path):
    garbled = " ".join(f"(cid:{i % 40})" for i in range(2500))
    _paper(
        tmp_path / "corpus",
        "KCNH2",
        "5",
        "# MAIN TEXT\n\nAnnals of Medicine\n\n" + garbled,
    )
    idx = _index(tmp_path, "KCNH2", "5")
    assert idx.body_garbled
    forms = sweep.FormIndex().get("KCNH2", "G572S")
    assert sweep.classify("G572S", forms, idx)[0] == "text_absent_garbled_body"


def test_figures_on_disk_keep_an_absent_row_penalized_as_unknown(tmp_path):
    _paper(
        tmp_path / "corpus",
        "RYR2",
        "6",
        _long_body("Pedigrees are shown in Figure 2."),
        **{"6_figures/fig2.png": "not really a png"},
    )
    idx = _index(tmp_path, "RYR2", "6")
    forms = sweep.FormIndex().get("RYR2", "R420W")
    klass, _ = sweep.classify("R420W", forms, idx)
    assert klass == "text_absent_figures_present"
    assert klass not in sweep.HARD_CEILING
    assert klass in sweep.WIDE_CEILING


def test_indel_notation_is_inconclusive_not_absent(tmp_path):
    _paper(
        tmp_path / "corpus",
        "SCN5A",
        "7",
        _long_body("An in-frame deletion removed residues 1500 to 1502."),
    )
    idx = _index(tmp_path, "SCN5A", "7")
    forms = sweep.FormIndex().get("SCN5A", "p.Lys1500del")
    assert (
        sweep.classify("p.Lys1500del", forms, idx)[0]
        == "text_absent_notation_inconclusive"
    )


def test_absent_substitution_with_a_real_body_is_the_hard_ceiling(tmp_path):
    _paper(
        tmp_path / "corpus",
        "SCN5A",
        "8",
        _long_body("No variant table survived the capture."),
    )
    idx = _index(tmp_path, "SCN5A", "8")
    forms = sweep.FormIndex().get("SCN5A", "A1330T")
    klass, _ = sweep.classify("A1330T", forms, idx)
    assert klass == "text_absent_substitution"
    assert klass in sweep.HARD_CEILING


def test_missing_paper_directory_is_source_absent(tmp_path):
    (tmp_path / "corpus").mkdir()
    idx = _index(tmp_path, "KCNQ1", "9")
    forms = sweep.FormIndex().get("KCNQ1", "A341V")
    assert sweep.classify("A341V", forms, idx) == ("source_absent", "no source on disk")


def test_summary_reports_both_ceilings_and_carrier_rows_behind_them():
    rows = [
        {
            "gene": "KCNQ1",
            "pmid": "1",
            "variant": "A",
            "class": "present_in_body",
            "gold_carriers": 3,
            "tranche": "t1",
            "gold_provenance": "p",
        },
        {
            "gene": "KCNQ1",
            "pmid": "1",
            "variant": "B",
            "class": "text_absent_substitution",
            "gold_carriers": 2,
            "tranche": "t1",
            "gold_provenance": "p",
        },
        {
            "gene": "KCNQ1",
            "pmid": "2",
            "variant": "C",
            "class": "text_absent_figures_present",
            "gold_carriers": 0,
            "tranche": "t1",
            "gold_provenance": "p",
        },
        {
            "gene": "KCNQ1",
            "pmid": "3",
            "variant": "D",
            "class": "source_absent",
            "gold_carriers": 1,
            "tranche": "t2",
            "gold_provenance": "q",
        },
    ]
    papers = {("KCNQ1", "1"): {}, ("KCNQ1", "2"): {}, ("KCNQ1", "3"): {}}
    summary = sweep.summarize(rows, papers)
    assert summary["gold_rows"] == 4
    assert summary["hard_ceiling_rows"] == 2
    assert summary["wide_ceiling_rows"] == 3
    assert summary["max_reachable_recall_if_hard_excluded_pct"] == 50.0
    assert summary["nonzero_carrier_rows"] == 3
    assert summary["nonzero_carrier_rows_behind_hard_ceiling"] == 2
    assert summary["by_tranche"]["t2"] == {"source_absent": 1}
    # The rendered markdown must carry every class column so the table cannot
    # silently drop a stratum.
    md = sweep.render_md(summary, "unit test")
    for klass in sweep.CLASSES:
        assert klass in md
    json.dumps(summary)  # must stay JSON-serialisable for summary.json
