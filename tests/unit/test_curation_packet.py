"""Tests for the cold-start manual-curation packet build + scoring glue."""

import csv
import json

from scripts.build_curation_packet import _full_text_papers, _protocol_md, _title
from scripts.score_curation_packet import convert, write_recall_input


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


# --- score_curation_packet conversion ---------------------------------------


def _filled(tmp_path):
    path = tmp_path / "filled.csv"
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "pmid",
                "variant",
                "germline_or_somatic",
                "carriers",
                "affected",
                "unaffected",
                "evidence_note",
            ]
        )
        w.writerow(["111", "c.5946delT", "germline", "14", "9", "5", "Table 2"])
        w.writerow(["111", "p.Cys1365Tyr", "unknown", "1", "", "", "Results"])
        w.writerow(["222", "c.9999A>T", "somatic", "3", "3", "", "tumor NGS"])
        w.writerow(["333", "NONE", "", "", "", "", "review"])
        w.writerow(["444", "", "", "", "", "", "blank variant"])
    return path


def test_convert_excludes_somatic_keeps_germline_and_unknown(tmp_path):
    rows, summary = convert(_filled(tmp_path), include_somatic=False)
    variants = {r["variant"] for r in rows}
    assert variants == {"c.5946delT", "p.Cys1365Tyr"}  # somatic + NONE/blank dropped
    assert summary["curated_pmids"] == 4  # 111,222,333,444 all count as curated
    assert summary["gold_variant_rows"] == 2
    assert summary["no_variant_papers"] == 2  # NONE + blank-variant
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
