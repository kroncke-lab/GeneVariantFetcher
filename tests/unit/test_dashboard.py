"""Offline test for the static provenance/coverage dashboard generator."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from cli.dashboard import generate_dashboard, md_to_html


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
