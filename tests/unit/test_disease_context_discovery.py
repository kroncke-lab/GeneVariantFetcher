from gene_literature.disease_context import (
    build_gene_disease_context,
    parse_clinigen_gene_validity_csv,
)
from gene_literature.discovery import build_gene_keyword_queries


def test_pubmed_disease_query_parenthesizes_gene_or_block_and_aliases():
    queries = build_gene_keyword_queries(
        "BRCA1",
        synonyms=["BRCAI", "FANCS"],
        disease="hereditary breast and ovarian cancer",
        disease_terms=["BRCA1-related cancer predisposition", "HBOC"],
    )

    assert queries[0] == (
        '("BRCA1"[Title/Abstract] OR "BRCAI"[Title/Abstract] OR '
        '"FANCS"[Title/Abstract]) AND '
        '("hereditary breast and ovarian cancer"[Title/Abstract] OR '
        '"BRCA1-related cancer predisposition"[Title/Abstract] OR '
        '"HBOC"[Title/Abstract])'
    )
    assert all(query.startswith("(") for query in queries)
    assert all(") AND (" in query for query in queries)


def test_clingen_csv_parser_skips_decorative_rows():
    csv_text = "\n".join(
        [
            '"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""',
            '"FILE CREATED: 2026-06-16","","","","","","","","",""',
            '"WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity","","","","","","","","",""',
            '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++","+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"',
            '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"',
            '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++","+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"',
            '"BRCA1","HGNC:1100","BRCA1-related cancer predisposition","MONDO:0700268","AD","SOP11","Definitive","https://example.test/report","2024-08-29T17:00:00.000Z","Hereditary Cancer Gene Curation Expert Panel"',
        ]
    )

    rows = parse_clinigen_gene_validity_csv(csv_text, "BRCA1")

    assert len(rows) == 1
    assert rows[0].disease_label == "BRCA1-related cancer predisposition"
    assert rows[0].classification == "Definitive"


def test_manual_context_adds_requested_gene_disease_aliases(monkeypatch):
    monkeypatch.setattr(
        "gene_literature.disease_context.fetch_clinigen_gene_validity_curations",
        lambda gene, timeout_s=30: [],
    )

    context = build_gene_disease_context(
        "APOE",
        "Alzheimer disease",
    )

    assert context.disease_terms[:2] == [
        "Alzheimer disease",
        "Alzheimer's disease",
    ]
    assert "LOAD" in context.disease_terms
    assert context.prompt_disease


def test_gene_only_context_includes_all_clingen_phenotypes(monkeypatch):
    from gene_literature.disease_context import ClinGenDiseaseCuration

    curations = [
        ClinGenDiseaseCuration(
            gene_symbol="BRCA1",
            disease_label="BRCA1-related cancer predisposition",
            disease_id="MONDO:0700268",
            mode_of_inheritance="AD",
            classification="Definitive",
            classification_date="2024-08-29",
            gcep="Hereditary Cancer",
            online_report="https://example.test/brca1",
        ),
        ClinGenDiseaseCuration(
            gene_symbol="BRCA1",
            disease_label="Fanconi anemia, complementation group S",
            disease_id="MONDO:0019391",
            mode_of_inheritance="AR",
            classification="Limited",
            classification_date="2024-08-29",
            gcep="Hereditary Cancer",
            online_report="https://example.test/fancs",
        ),
    ]
    monkeypatch.setattr(
        "gene_literature.disease_context.fetch_clinigen_gene_validity_curations",
        lambda gene, timeout_s=30: curations,
    )

    context = build_gene_disease_context("BRCA1", None)

    assert "BRCA1-related cancer predisposition" in context.disease_terms
    assert "Fanconi anemia, complementation group S" in context.disease_terms
    assert len(context.selected_clinigen_curations) == 2


def test_requested_context_records_but_does_not_select_unrelated_clingen_rows(
    monkeypatch,
):
    from gene_literature.disease_context import ClinGenDiseaseCuration

    curations = [
        ClinGenDiseaseCuration(
            gene_symbol="BRCA1",
            disease_label="BRCA1-related cancer predisposition",
            disease_id="MONDO:0700268",
            mode_of_inheritance="AD",
            classification="Definitive",
            classification_date="2024-08-29",
            gcep="Hereditary Cancer",
            online_report="https://example.test/brca1",
        ),
        ClinGenDiseaseCuration(
            gene_symbol="BRCA1",
            disease_label="Fanconi anemia, complementation group S",
            disease_id="MONDO:0019391",
            mode_of_inheritance="AR",
            classification="Limited",
            classification_date="2024-08-29",
            gcep="Hereditary Cancer",
            online_report="https://example.test/fancs",
        ),
    ]
    monkeypatch.setattr(
        "gene_literature.disease_context.fetch_clinigen_gene_validity_curations",
        lambda gene, timeout_s=30: curations,
    )

    context = build_gene_disease_context(
        "BRCA1",
        "hereditary breast and ovarian cancer",
    )

    assert len(context.clinigen_curations) == 2
    assert "Fanconi anemia, complementation group S" not in context.disease_terms
    assert [c.disease_label for c in context.selected_clinigen_curations] == [
        "BRCA1-related cancer predisposition"
    ]


def test_requested_context_can_include_all_clingen_phenotypes(monkeypatch):
    from gene_literature.disease_context import ClinGenDiseaseCuration

    curations = [
        ClinGenDiseaseCuration(
            gene_symbol="MYBPC3",
            disease_label="hypertrophic cardiomyopathy",
            disease_id="MONDO:0005045",
            mode_of_inheritance="AD",
            classification="Definitive",
            classification_date="2024-08-29",
            gcep="Cardiovascular",
            online_report="https://example.test/hcm",
        ),
        ClinGenDiseaseCuration(
            gene_symbol="MYBPC3",
            disease_label="dilated cardiomyopathy",
            disease_id="MONDO:0005021",
            mode_of_inheritance="AD",
            classification="Disputed",
            classification_date="2024-08-29",
            gcep="Cardiovascular",
            online_report="https://example.test/dcm",
        ),
    ]
    monkeypatch.setattr(
        "gene_literature.disease_context.fetch_clinigen_gene_validity_curations",
        lambda gene, timeout_s=30: curations,
    )

    context = build_gene_disease_context(
        "MYBPC3",
        "hypertrophic cardiomyopathy",
        include_all_clinigen_phenotypes=True,
    )

    assert "hypertrophic cardiomyopathy" in context.disease_terms
    assert "dilated cardiomyopathy" in context.disease_terms
    assert len(context.selected_clinigen_curations) == 2


def test_discovery_includes_penetrance_segregation_lane():
    # Criticisms 3 and 6: a dedicated lane must recruit carrier-first
    # segregation / penetrance / cascade studies, not just case series.
    queries = build_gene_keyword_queries("BRCA1")
    lane = [q for q in queries if "penetrance" in q.lower()]
    assert lane, "expected a penetrance/segregation discovery lane"
    joined = lane[0].lower()
    for term in ("segregation", "cascade", "unaffected carrier", "prospective"):
        assert term in joined
