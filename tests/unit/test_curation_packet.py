"""Tests for the cold-start manual-curation packet build + scoring glue."""

import csv
import json

import pytest

from scripts.build_curation_packet import (
    _full_text_papers,
    _protocol_md,
    _read_pmids,
    _split_for_pmid,
    _title,
)
from scripts.score_curation_packet import (
    CurationValidationError,
    convert,
    write_pmid_scope,
    write_recall_input,
)


def _run_dir(tmp_path):
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "111_FULL_CONTEXT.md").write_text("x" * 4000)  # full text
    (ft / "222_FULL_CONTEXT.md").write_text("x" * 1000)  # stub -> excluded
    (ft / "abc_FULL_CONTEXT.md").write_text("x" * 4000)  # non-digit -> excluded
    aj = tmp_path / "abstract_json"
    aj.mkdir()
    (aj / "111.json").write_text(json.dumps({"metadata": {"title": "BRCA2 carriers"}}))
    return tmp_path


# --- build_curation_packet --------------------------------------------------


def test_full_text_papers_filters_stubs_and_nondigits(tmp_path):
    pool = _full_text_papers(_run_dir(tmp_path), min_bytes=3000)
    assert [p[0] for p in pool] == ["111"]
    assert pool[0][2] == 4000


def test_title_reads_abstract_metadata(tmp_path):
    run = _run_dir(tmp_path)
    assert _title(run, "111") == "BRCA2 carriers"
    assert _title(run, "999") == ""  # missing -> empty


def test_protocol_mentions_gene_count_and_schema():
    md = _protocol_md("BRCA2", 50)
    assert "BRCA2" in md and "50 papers" in md
    assert "germline_or_somatic" in md and "NONE" in md  # the key rules
    assert "carrier_count_role" in md and "Family counts are not" in md


def test_exact_roster_reader_and_split_are_stable(tmp_path):
    roster = tmp_path / "pmids.txt"
    roster.write_text("333\n111 # retained\n222\n", encoding="utf-8")
    assert _read_pmids(roster) == ["333", "111", "222"]
    first = {
        pmid: _split_for_pmid(pmid, "BRCA2", 42, 2, _read_pmids(roster))
        for pmid in _read_pmids(roster)
    }
    second = {
        pmid: _split_for_pmid(pmid, "BRCA2", 42, 2, _read_pmids(roster))
        for pmid in reversed(_read_pmids(roster))
    }
    assert first == second
    assert list(first.values()).count("calibration") == 2


def test_exact_roster_reader_rejects_duplicates(tmp_path):
    roster = tmp_path / "pmids.txt"
    roster.write_text("111\n111\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate PMID"):
        _read_pmids(roster)


# --- score_curation_packet conversion ---------------------------------------


def _filled(tmp_path):
    path = tmp_path / "filled.csv"
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "pmid",
                "variant",
                "curation_status",
                "germline_or_somatic",
                "carriers",
                "carrier_count_role",
                "affected",
                "unaffected",
                "evidence_note",
            ]
        )
        w.writerow(
            [
                "111",
                "c.5946delT",
                "complete",
                "germline",
                "14",
                "individual",
                "9",
                "5",
                "Table 2",
            ]
        )
        w.writerow(
            [
                "111",
                "p.Cys1365Tyr",
                "complete",
                "unknown",
                "1",
                "individual",
                "",
                "",
                "Results",
            ]
        )
        w.writerow(
            [
                "222",
                "c.9999A>T",
                "complete",
                "somatic",
                "3",
                "individual",
                "3",
                "",
                "tumor NGS",
            ]
        )
        w.writerow(
            ["333", "NONE", "complete", "", "", "not_reported", "", "", "review"]
        )
        w.writerow(
            [
                "444",
                "NONE",
                "complete",
                "",
                "",
                "not_reported",
                "",
                "",
                "tables checked",
            ]
        )
    return path


def test_convert_excludes_somatic_keeps_germline_and_unknown(tmp_path):
    rows, summary = convert(_filled(tmp_path), include_somatic=False)
    variants = {r["variant"] for r in rows}
    assert variants == {"c.5946delT", "p.Cys1365Tyr"}  # somatic + NONE dropped
    assert summary["curated_pmids"] == 4  # 111,222,333,444 all count as curated
    assert summary["gold_variant_rows"] == 2
    assert summary["no_variant_papers"] == 2
    assert summary["dropped_somatic"] == 1


def test_convert_include_somatic_keeps_it(tmp_path):
    rows, summary = convert(_filled(tmp_path), include_somatic=True)
    assert any(r["variant"] == "c.9999A>T" for r in rows)
    assert summary["dropped_somatic"] == 0


def test_write_recall_input_matches_gold_schema(tmp_path):
    rows, _ = convert(_filled(tmp_path), include_somatic=False)
    path = write_recall_input(rows, "BRCA2", tmp_path / "gold")
    assert path.name == "BRCA2_recall_input.csv"
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == ["variant", "pmid", "carriers", "affected", "unaffected"]
        first = next(reader)
        assert first == ["c.5946delT", "111", "14", "9", "5"]


def test_convert_rejects_unfinished_blank_variant(tmp_path):
    path = _filled(tmp_path)
    text = path.read_text(encoding="utf-8")
    path.write_text(
        text.replace("444,NONE,complete", "444,,complete"), encoding="utf-8"
    )
    with pytest.raises(CurationValidationError, match="blank variant is unfinished"):
        convert(path, include_somatic=False)


def test_convert_rejects_counts_on_explicit_none(tmp_path):
    path = _filled(tmp_path)
    text = path.read_text(encoding="utf-8")
    path.write_text(
        text.replace(
            "333,NONE,complete,,,not_reported,,,review",
            "333,NONE,complete,germline,1,individual,,,review",
        ),
        encoding="utf-8",
    )
    with pytest.raises(CurationValidationError, match="NONE cannot carry"):
        convert(path, include_somatic=False)


def test_convert_rejects_duplicate_variant_row(tmp_path):
    path = _filled(tmp_path)
    with path.open("a", encoding="utf-8") as f:
        f.write("111,c.5946delT,complete,germline,14,individual,9,5,duplicate\n")
    with pytest.raises(CurationValidationError, match="duplicate variant row"):
        convert(path, include_somatic=False)


def test_convert_requires_exact_packet_roster(tmp_path):
    with pytest.raises(CurationValidationError, match="roster mismatch"):
        convert(
            _filled(tmp_path),
            include_somatic=False,
            expected_pmids={"111", "222", "333", "444", "555"},
        )


def test_family_count_is_not_scored_as_individual_mae(tmp_path):
    path = _filled(tmp_path)
    text = path.read_text(encoding="utf-8")
    path.write_text(
        text.replace(
            "p.Cys1365Tyr,complete,unknown,1,individual",
            "p.Cys1365Tyr,complete,unknown,2,family",
        ),
        encoding="utf-8",
    )
    rows, summary = convert(path, include_somatic=False)
    family = next(row for row in rows if row["variant"] == "p.Cys1365Tyr")
    assert family["carriers"] == ""
    assert summary["nonindividual_counts_excluded_from_mae"] == 1


def test_write_pmid_scope_keeps_none_papers(tmp_path):
    path = write_pmid_scope({"333", "111"}, "BRCA2", tmp_path)
    assert path.read_text(encoding="utf-8") == "111\n333\n"
