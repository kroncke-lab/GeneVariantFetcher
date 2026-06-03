from scripts.recall_audit import build_acquisition_worklist as worklist


def test_worklist_routes_missing_target_gene_supplement_as_source_gap(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        "# Article\n\nAbstract\nMethods\nResults\n\n"
        "All SCN5A mutations are listed in Supplemental Table 2.\n"
        + ("main text without individual variant identifiers. " * 80),
        encoding="utf-8",
    )
    row = {
        "gene": "SCN5A",
        "pmid": "123",
        "source_status": "recovered_pmc",
        "available_source_status": "recovered_pmc",
        "context_path": str(source),
        "available_context_path": str(source),
        "context_bytes": str(source.stat().st_size),
        "available_context_bytes": str(source.stat().st_size),
        "source_desync": "False",
        "source_unbound": "False",
        "failure_class": "available_source_underextraction",
        "missing_rows": "87",
        "missing_carriers": "109",
        "top_missing_variants": "A1416G;A1698T",
    }

    item = worklist.classify_row(
        row,
        gene="SCN5A",
        missing_distinct_variants=87,
        doi="10.1016/j.hrthm.2018.01.014",
        doi_source="fixture.csv",
    )

    assert item["action"] == "fetch"
    assert item["route"] == "fetch_elsevier_insttoken"
    assert item["likely_missing_variant_supplement"] is True


def test_worklist_does_not_flag_present_supplement_table_body(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        "# Article\n\nAll SCN5A mutations are listed in Supplemental Table 2.\n\n"
        "## Supplemental Table 2\n\n"
        "| Variant | Count |\n| --- | ---: |\n| A1416G | 1 |\n" + ("table text. " * 80),
        encoding="utf-8",
    )
    row = {
        "gene": "SCN5A",
        "pmid": "123",
        "source_status": "recovered_pmc",
        "available_source_status": "recovered_pmc",
        "context_path": str(source),
        "available_context_path": str(source),
        "context_bytes": str(source.stat().st_size),
        "available_context_bytes": str(source.stat().st_size),
        "source_desync": "False",
        "source_unbound": "False",
        "failure_class": "available_source_underextraction",
        "missing_rows": "3",
    }

    item = worklist.classify_row(
        row,
        gene="SCN5A",
        missing_distinct_variants=3,
        doi="",
        doi_source="",
    )

    assert item["route"] == "available_source_underextraction"
    assert item["likely_missing_variant_supplement"] is False
