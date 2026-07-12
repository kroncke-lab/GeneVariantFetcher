"""Offline test for the static provenance/coverage dashboard generator."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from collections import Counter

from cli.dashboard import (
    _delta_compact,
    _delta_html,
    _gene_snapshot,
    _health_band,
    _health_signals,
    _load_paper_final_check,
    _load_prev_snapshot,
    _sanitize_local_paths,
    _worklist_card,
    generate_dashboard,
    md_to_html,
)


def _make_corpus(corpus: Path) -> None:
    (corpus).mkdir(parents=True, exist_ok=True)
    (corpus / "INDEX.csv").write_text(
        "gene,pmid,fulltext,fulltext_bytes,source_sha256,full_text_status,"
        "n_figures,n_supplement_files,n_source_copies,chosen_from\n"
        "TESTGENE,111,corpus/TESTGENE/111/111_FULL_CONTEXT.md,6000,abc,ok,2,1,1,run\n"
        "TESTGENE,222,corpus/TESTGENE/222/222_FULL_CONTEXT.md,300,def,stub,0,0,1,run\n",
        encoding="utf-8",
    )
    d = corpus / "TESTGENE" / "111"
    d.mkdir(parents=True, exist_ok=True)
    (d / "111_FULL_CONTEXT.md").write_text(
        "# Results\n\nThe proband carried p.Arg100Trp and was affected at age 12.\n\n"
        "| Variant | Patient |\n|---|---|\n| p.Gly300Ser | case 3 |\n",
        encoding="utf-8",
    )
    (d / "111_CLEANED.md").write_text(
        "# Cleaned Results\n\nCLEANED SOURCE p.Arg100Trp evidence text.\n\n"
        "File available at: /Users/example/GeneVariantFetcher/results/source.pdf\n",
        encoding="utf-8",
    )
    (d / "111_artifacts.json").write_text(
        json.dumps({"pmid": "111", "doi": "10.1/x", "pmcid": "PMC1"}), encoding="utf-8"
    )


def _make_db(db: Path) -> None:
    con = sqlite3.connect(str(db))
    con.executescript(
        """
        CREATE TABLE papers(pmid TEXT PRIMARY KEY, title TEXT, first_author TEXT,
            journal TEXT, publication_date TEXT, doi TEXT, pmc_id TEXT, gene_symbol TEXT);
        CREATE TABLE variants(variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT, genomic_position TEXT, clinical_significance TEXT);
        CREATE TABLE variant_papers(variant_id INTEGER, pmid TEXT, source_location TEXT,
            additional_notes TEXT, key_quotes TEXT, count_provenance TEXT);
        CREATE TABLE extraction_metadata(pmid TEXT, model_used TEXT, source_file TEXT,
            source_type TEXT, extraction_confidence TEXT, total_variants_found INTEGER,
            notes TEXT, abstract_only INTEGER);
        CREATE TABLE individual_records(variant_id INTEGER, pmid TEXT, individual_id TEXT,
            age_at_evaluation INTEGER, age_at_onset INTEGER, sex TEXT, affected_status TEXT,
            phenotype_details TEXT, evidence_sentence TEXT, ethnicity TEXT, geographic_origin TEXT);
        """
    )
    con.execute(
        "INSERT INTO papers VALUES ('111','A KCNH2 family study','Chen Q','Circulation','1999','10.1/x',NULL,'TESTGENE')"
    )
    # A malformed extraction row whose pmid is "N/A": the slash must NOT abort the
    # build by trying to write papers/N/A.html. It should be silently skipped.
    con.execute(
        "INSERT INTO papers VALUES ('N/A','junk row',NULL,NULL,NULL,NULL,NULL,'TESTGENE')"
    )
    con.execute(
        "INSERT INTO variants VALUES (1,'TESTGENE',NULL,'p.Arg100Trp',NULL,'Pathogenic')"
    )
    con.execute(
        "INSERT INTO variant_papers VALUES (1,'111','Table 1, case 3','proband; LQT2',"
        '\'["carried p.Arg100Trp"]\',\'{"carriers_column_label":"N carriers","carriers_count_type":"per_variant_carrier"}\')'
    )
    source_file = db.parent / "corpus" / "TESTGENE" / "111" / "111_CLEANED.md"
    con.execute(
        "INSERT INTO extraction_metadata VALUES ('111','test-model',?,'fulltext','high',1,'',0)",
        (str(source_file),),
    )
    con.execute(
        "INSERT INTO individual_records VALUES (1,'111','proband',12,10,'F','affected',"
        "'long QT','The proband carried p.Arg100Trp and was affected at age 12.','East Asian','Japan')"
    )
    con.commit()
    con.close()


def test_sanitize_local_paths_redacts_mac_and_windows_paths():
    mac = _sanitize_local_paths("see /Users/brett/GVF/results/source.pdf now")
    assert "/Users/brett" not in mac
    assert "[local path: source.pdf]" in mac

    win_backslash = _sanitize_local_paths("see C:\\Users\\brett\\GVF\\source.pdf now")
    assert "C:\\Users" not in win_backslash
    assert "[local path: source.pdf]" in win_backslash

    win_forwardslash = _sanitize_local_paths("see C:/Users/brett/GVF/source.pdf now")
    assert "C:/Users" not in win_forwardslash
    assert "[local path: source.pdf]" in win_forwardslash


def test_md_to_html_renders_blocks_and_tables():
    html = md_to_html("# Heading\n\npara text\n\n| a | b |\n|---|---|\n| 1 | 2 |\n")
    assert "<h2" in html and "Heading" in html
    assert '<p class="blk">para text</p>' in html
    assert "<table" in html and "<td>1</td>" in html


def test_generate_dashboard_produces_pages_links_and_jump(tmp_path: Path):
    corpus = tmp_path / "corpus"
    _make_corpus(corpus)
    db = tmp_path / "TESTGENE.db"
    _make_db(db)
    out = tmp_path / "dash"

    stats = generate_dashboard(
        out_dir=out,
        corpus_dir=corpus,
        db_map={"TESTGENE": db},
        genes=["TESTGENE"],
        max_papers=0,
        generated="2026-06-04",
        score=False,
    )

    assert stats["genes"] == 1 and stats["paper_pages"] >= 1
    # the malformed "N/A" pmid was skipped, not written as a path
    assert not (out / "TESTGENE" / "papers" / "N").exists()
    assert not (out / "TESTGENE" / "papers" / "N" / "A.html").exists()
    index = (out / "index.html").read_text()
    assert "TESTGENE" in index and "papers in corpus" in index
    assert "TESTGENE/index.html" in index  # per-gene subdir layout
    gene = (out / "TESTGENE" / "index.html").read_text()
    assert "Source coverage" in gene and "Provenance completeness" in gene
    assert "papers/111.html" in gene  # links to the adjudication page

    paper = (out / "TESTGENE" / "papers" / "111.html").read_text()
    assert "pubmed.ncbi.nlm.nih.gov/111/" in paper  # PubMed link
    assert "doi.org/10.1/x" in paper  # DOI link from artifacts.json
    assert "ncbi.nlm.nih.gov/pmc/articles/PMC1/" in paper  # PMC link
    assert "p.Arg100Trp" in paper  # the extracted variant
    assert "Table 1, case 3" in paper  # source_location jump target
    # click-to-highlight wired AND the onclick arg is HTML-escaped (&quot;) so an
    # inner double-quote can't terminate the onclick="..." attribute (the bug that
    # silently broke every clickthrough). Raw jump("..." would be malformed.
    assert "jump(&quot;" in paper
    assert 'onclick="jump("' not in paper
    assert "CLEANED SOURCE p.Arg100Trp evidence text" in paper
    assert "extraction_metadata.source_file" in paper
    assert str(tmp_path) not in paper
    assert "/Users/example" not in paper
    assert "[local path: source.pdf]" in paper
    assert "source_file: 111_CLEANED.md" in paper
    assert "id='src'" in paper  # source pane present
    # new provenance fields
    assert (
        "Chen Q" in paper and "Circulation" in paper
    )  # first author + journal bibliography
    assert "East Asian" in paper and "Japan" in paper  # patient ethnicity / origin
    assert "count basis" in paper and "N carriers" in paper  # count_provenance ("why")

    # variant-centric (website-style) page
    variants = (out / "TESTGENE" / "variants.html").read_text()
    assert "p.Arg100Trp" in variants and "unique variants" in variants
    assert "papers/111.html" in variants  # clickthrough to adjudication


def _add_trust_and_finalcheck(db: Path) -> None:
    """Widen the base test DB with penetrance_data (trust-tiered) and a flagged
    paper_final_check row, so the health / worklist / badge surfaces have data."""
    con = sqlite3.connect(str(db))
    con.executescript(
        """
        CREATE TABLE penetrance_data(penetrance_id INTEGER PRIMARY KEY, variant_id INTEGER,
            pmid TEXT, total_carriers_observed INTEGER, affected_count INTEGER,
            unaffected_count INTEGER, trust_tier TEXT, trust_reasons TEXT, trust_rule_version TEXT);
        CREATE TABLE paper_final_check(pmid TEXT PRIMARY KEY, gene_symbol TEXT, verdict TEXT,
            confidence REAL, n_facts INTEGER, n_flagged INTEGER, flags_json TEXT, summary TEXT,
            checked_at TEXT, completeness_status TEXT, n_missing INTEGER, source_grounded INTEGER);
        CREATE TABLE variants2(x INTEGER);
        """
    )
    con.execute(
        "INSERT INTO penetrance_data VALUES (1,1,'111',10,8,2,'trusted','[]','trust-v1')"
    )
    # a wildly-high carrier count the trust gate soft-quarantined
    con.execute(
        "INSERT INTO penetrance_data VALUES (2,1,'111',999,0,0,'quarantine',"
        "'[\"paper_outlier\"]','trust-v1')"
    )
    con.execute(
        "INSERT INTO paper_final_check VALUES ('111','TESTGENE','flag',0.2,3,1,"
        "'[\"carrier count exceeds cohort\"]','count exceeds cohort size',"
        "'2026-06-04','incomplete',3,1)"
    )
    con.commit()
    con.close()


def test_dashboard_health_worklist_badges_and_delta(tmp_path: Path):
    corpus = tmp_path / "corpus"
    _make_corpus(corpus)
    db = tmp_path / "TESTGENE.db"
    _make_db(db)
    _add_trust_and_finalcheck(db)
    # Two more papers to exercise the extraction split: 333 = usable text on disk
    # but never extracted (throughput lever); 444 = extracted but 0 variants (a
    # genuine miss). These must show up as DISTINCT levers, not one lumped count.
    with (corpus / "INDEX.csv").open("a", encoding="utf-8") as f:
        f.write(
            "TESTGENE,333,corpus/TESTGENE/333/333_FULL_CONTEXT.md,7000,g33,ok,0,0,1,run\n"
            "TESTGENE,444,corpus/TESTGENE/444/444_FULL_CONTEXT.md,7000,g44,ok,0,0,1,run\n"
        )
    con = sqlite3.connect(str(db))
    con.execute(
        "INSERT INTO papers VALUES ('444','Extracted but empty',NULL,NULL,NULL,NULL,NULL,'TESTGENE')"
    )
    con.execute(
        "INSERT INTO extraction_metadata VALUES ('444','test-model','','fulltext','low',0,'',0)"
    )
    # a 'skipped' final-check verdict must render muted, NOT as a green "ok"
    con.execute(
        "INSERT INTO paper_final_check VALUES ('444','TESTGENE','skipped',NULL,0,0,"
        "'[]','no usable source','2026-06-04','',0,0)"
    )
    con.commit()
    con.close()
    out = tmp_path / "dash"

    def build(generated: str):
        return generate_dashboard(
            out_dir=out,
            corpus_dir=corpus,
            db_map={"TESTGENE": db},
            genes=["TESTGENE"],
            max_papers=0,
            generated=generated,
            score=False,
        )

    build("2026-06-04")
    gene = (out / "TESTGENE" / "index.html").read_text()
    # no-gold run-health panel
    assert "Run health" in gene
    assert "trusted / quarantined counts" in gene
    assert "paper_outlier" in gene  # top quarantine reason surfaced
    assert "final-check papers ok / flagged" in gene
    # ranked worklist with linked levers + suggested commands
    assert "What to change next" in gene
    assert "Acquire missing full text" in gene  # PMID 222 is a stub
    assert "Review quarantined counts" in gene
    assert "trust_report.py" in gene  # suggested command
    assert "Adjudicate final-check flags" in gene
    # the extraction split: throughput (333, not yet extracted) is a DIFFERENT
    # lever from a genuine miss (444, extracted but empty)
    assert "Extract usable papers not yet in the DB" in gene
    assert "Investigate extraction misses (extracted, 0 variants)" in gene
    assert "extracted but 0 variants" in gene  # health card line
    # per-paper Check column shows the flag (checker says 3 carriers missed)
    assert "Check" in gene and "flag −3" in gene
    # the 'skipped' verdict renders muted, not as a green pass
    assert "data-v='2'>skipped" in gene
    # first run: delta line explains itself, no chip yet
    assert "First tracked run" in gene

    index = (out / "index.html").read_text()
    assert "trust 1✓ / 1 quar" in index  # overview badge
    assert "check 1 flagged" in index

    # history persisted (gitignored corpus dir), one snapshot so far
    hist = corpus / "dashboard_history" / "TESTGENE.jsonl"
    assert hist.exists() and len([x for x in hist.read_text().splitlines() if x]) == 1

    # second run after the DB gained a unique variant -> a real delta appears
    con = sqlite3.connect(str(db))
    con.execute(
        "INSERT INTO variants VALUES (9,'TESTGENE',NULL,'p.Cys500Tyr',NULL,'Pathogenic')"
    )
    con.execute("INSERT INTO variant_papers VALUES (9,'111','Table 2','','[]',NULL)")
    con.commit()
    con.close()

    build("2026-06-05")
    gene2 = (out / "TESTGENE" / "index.html").read_text()
    assert "Since 2026-06-04" in gene2 and "unique variants +1" in gene2
    index2 = (out / "index.html").read_text()
    assert "uniqV +1 since last run" in index2
    assert len([x for x in hist.read_text().splitlines() if x]) == 2


def _bare_db() -> dict:
    return {
        "trust": {"tiered": False},
        "final_check": {"by_pmid": {}, "counts": Counter()},
        "quarantined": [],
    }


def test_worklist_gold_gap_lever_and_empty_state():
    # gold gene short of 90% -> the headline recall lever renders with the right
    # gap (ceil(100*0.9) - 80 = 10) and the corrected recovery path.
    recall = {"unique_variants": (80, 100, 0.8)}
    per_pmid = {"111": {"gold": 5, "missing": 3}, "222": {"gold": 2, "missing": 0}}
    html = _worklist_card(
        "KCNH2", {"extracted": 1}, _bare_db(), [], [], [], {"111"}, recall, per_pmid
    )
    assert "Close the gap to 90% unique-variant recall (+10)" in html
    assert "scripts/recall_recovery/run_all_layers.py" in html  # corrected path
    assert "papers/111.html" in html  # links the worst-missing paper

    # no-DB gene with usable text on disk -> an "extract these" lever, not the
    # misleading "no open levers" empty state.
    nodb = _worklist_card(
        "KCNH2", {"extracted": 0}, None, [], [], ["333", "444"], set(), None, None
    )
    assert "Extract usable papers not yet in the DB" in nodb

    # genuinely-healthy gene -> the empty-state card
    empty = _worklist_card(
        "KCNH2", {"extracted": 1}, _bare_db(), [], [], [], set(), None, None
    )
    assert "No open levers detected" in empty


def test_delta_html_direction_no_change_and_formatting():
    # a lower_better field improving (5 -> 2) reads green/ok, not red
    prev = {"quarantine": 5, "unique": 10, "generated": "2026-06-01"}
    cur = {"quarantine": 2, "unique": 10}
    html = _delta_html(prev, cur)
    assert "quarantined counts -3" in html and "tag ok" in html

    # nothing tracked changed -> honest message, not "No change" contradicting a
    # snapshot that a None<->number availability flip would still append
    assert "No tracked metric changed" in _delta_html(
        {"unique": 10, "generated": "x"}, {"unique": 10, "quarantine": None}
    )
    # float recall delta formats without scientific notation
    rec = _delta_html(
        {"recall_uv_pct": 80.0, "generated": "x"}, {"recall_uv_pct": 85.5}
    )
    assert "uniqV recall % +5.5" in rec
    # a huge integer delta is not rendered in sci notation
    big = _delta_html(
        {"variant_rows": 0, "generated": "x"}, {"variant_rows": 2_000_000}
    )
    assert "+2000000" in big and "e+0" not in big


def test_gene_snapshot_avoids_none_traps():
    s = {
        "usable": 5,
        "stub": 1,
        "extracted": 3,
        "with_variants": 2,
        "unique": 10,
        "variant_rows": 12,
    }
    # final check ran, ALL ok (no "flag" key in the Counter) -> 0, not None
    db = {
        "trust": {"tiered": True, "total": 4, "quarantine": 1},
        "final_check": {"counts": Counter({"total": 3, "ok": 3})},
    }
    snap = _gene_snapshot(s, db, {"unique_variants": (8, 10, 0.8)})
    assert snap["flagged_papers"] == 0 and snap["missing_carriers"] == 0
    assert snap["quarantine"] == 1 and snap["recall_uv_pct"] == 80.0

    # trust gate ran but tiered 0 facts -> None (must not read as a real 0)
    empty = _gene_snapshot(
        s,
        {"trust": {"tiered": True, "total": 0}, "final_check": {"counts": Counter()}},
        None,
    )
    assert empty["quarantine"] is None and empty["flagged_papers"] is None

    # trust_tier COLUMN present but every row untiered (gate never ran), and a
    # final check that only skipped: both must read as not-measured (None), not 0.
    untiered = _gene_snapshot(
        s,
        {
            "trust": {"tiered": True, "total": 50, "trusted": 0, "quarantine": 0},
            "final_check": {"counts": Counter({"total": 5, "skipped": 5})},
        },
        None,
    )
    assert untiered["quarantine"] is None  # untiered rows != a real 0% quarantine
    assert untiered["flagged_papers"] is None  # all-skipped != a real 0 flags


def test_health_signals_untiered_trust_is_not_healthy():
    # a present trust_tier column whose rows are all untiered must NOT gate green
    db = {
        "trust": {"tiered": True, "total": 50, "trusted": 0, "quarantine": 0},
        "final_check": {"counts": Counter()},
    }
    sig = _health_signals({"extracted": 10}, db, 0, 0, 10)
    assert sig["trust_gated"] is False and sig["quar_rate"] is None


def test_delta_compact_guards_corrupt_values():
    assert _delta_compact({"unique": 3}, {"unique": 5}) == (
        "<span class='tag ok'>uniqV +2 since last run</span>"
    )
    assert _delta_compact({"unique": "bad"}, {"unique": 5}) == ""  # no crash, no chip
    assert _delta_compact({"unique": 5}, {"unique": 5}) == ""


def test_load_prev_snapshot_scans_back_past_corrupt_last_lines(tmp_path: Path):
    hist = tmp_path / "hist"
    hist.mkdir()
    # last well-formed DICT is the first line; the rest are corrupt / non-object.
    # The dict carries a non-ASCII char so the read must be utf-8 (the platform
    # default would UnicodeDecodeError on Windows — a ValueError, not OSError).
    (hist / "G.jsonl").write_text(
        '{"unique": 5, "generated": "dµ"}\n[]\nnot json at all\n', encoding="utf-8"
    )
    assert _load_prev_snapshot(hist, "G") == {"unique": 5, "generated": "dµ"}


def test_load_paper_final_check_coerces_nonstring_columns(tmp_path: Path):
    # SQLite columns can hold unexpected types; text-ish fields must be coerced
    # to str so later slicing / escaping can't TypeError.
    db = tmp_path / "fc.db"
    con = sqlite3.connect(str(db))
    con.row_factory = sqlite3.Row  # load_db sets this; _load_paper_final_check needs it
    con.execute(
        "CREATE TABLE paper_final_check(pmid TEXT, verdict TEXT, n_flagged INTEGER, "
        "n_missing INTEGER, completeness_status, summary)"
    )
    # summary + completeness_status stored as INTEGERs (dynamic typing)
    con.execute("INSERT INTO paper_final_check VALUES ('111','flag',1,2,0,12345)")
    con.commit()
    out = _load_paper_final_check(con)
    con.close()
    rec = out["by_pmid"]["111"]
    assert isinstance(rec["summary"], str) and rec["summary"] == "12345"
    assert isinstance(rec["completeness_status"], str)
    assert out["counts"]["flag"] == 1


def test_health_signals_reviewable_and_denominators():
    db = {
        "trust": {
            "tiered": True,
            "total": 100,
            "quarantine": 8,
            "quarantine_rate": 0.08,
            "trusted": 92,
            "by_reason": {"paper_outlier": 8},
        },
        "final_check": {
            "counts": Counter({"total": 100, "ok": 34, "flag": 6, "skipped": 60})
        },
    }
    sig = _health_signals(
        {"extracted": 50},
        db,
        zero_var_extracted=3,
        usable_unextracted=10,
        usable_extracted=40,
    )
    # flag rate is over reviewable (ok+flag=40), NOT the 100 raw rows
    assert sig["fc_reviewable"] == 40
    assert abs(sig["fc_flag_rate"] - 6 / 40) < 1e-9
    # miss rate is over usable-extracted (40), NOT all extracted (50)
    assert abs(sig["zero_rate"] - 3 / 40) < 1e-9
    assert _health_band(sig)[0] in ("warn", "bad")  # 6/40 = 0.15 -> at least "watch"

    # schema-only trust gate (tiered but 0 facts) is not a healthy 0%
    sig2 = _health_signals(
        {"extracted": 0},
        {"trust": {"tiered": True, "total": 0}, "final_check": {"counts": Counter()}},
        0,
        0,
        0,
    )
    assert sig2["trust_gated"] is False and sig2["quar_rate"] is None
