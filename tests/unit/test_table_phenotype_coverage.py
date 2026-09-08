"""Coverage metadata must never turn a mixed endpoint into a clinical fact."""

import copy

from pipeline.extraction import ExpertExtractor
from pipeline.table_phenotype_coverage import audit_table_phenotype_coverage
from utils.models import ExtractionResult, Paper

SOURCE = """
Table 1: Mutation list
Mutation        Patients
p.Arg12Trp      7

Table 2: Clinical characteristics by mutation at diagnosis
Mutation        Asymptomatic       Syncope
p.Arg12Trp      5                  2

Table 3: Mutations per ECG phenotype at follow-up
Mutation        Negative ECG       Disease
p.Arg12Trp      4                  3
"""


def test_distinct_endpoints_report_candidates_without_mutating_counts():
    data = {
        "variants": [
            {
                "protein_notation": "p.Arg12Trp",
                "penetrance_data": {"total_carriers_observed": 7},
            }
        ]
    }
    before = copy.deepcopy(data)
    audit = audit_table_phenotype_coverage(data, SOURCE)
    assert data == before
    assert audit["supplied_rows"] == {
        "total_carriers_observed": 1,
        "affected_count": 0,
        "unaffected_count": 0,
    }
    assert audit["missing_phenotype_fields"] == ["affected_count", "unaffected_count"]
    assert audit["status"] == "review_candidates"
    assert [c["caption"] for c in audit["clinical_table_candidates"]] == [
        "Table 2: Clinical characteristics by mutation at diagnosis",
        "Table 3: Mutations per ECG phenotype at follow-up",
    ]
    for c in audit["clinical_table_candidates"]:
        assert SOURCE.splitlines()[c["source_line_1based"] - 1] == c["caption"]


def test_no_candidates_is_not_a_claim_of_completeness_and_zero_is_supplied():
    audit = audit_table_phenotype_coverage(
        {
            "variants": [
                {
                    "penetrance_data": {
                        "total_carriers_observed": 4,
                        "affected_count": 0,
                        "unaffected_count": 4,
                    }
                }
            ]
        },
        "unlabelled source",
    )
    assert audit["supplied_rows"]["affected_count"] == 1
    assert audit["missing_phenotype_fields"] == []
    assert audit["status"] == "completeness_not_established"


def test_candidates_are_bounded_and_body_table_mentions_are_not_headings():
    source = "Table 2 shows clinical characteristics.\n" + "\n".join(
        f"Table {i}: Clinical features\n| variant | affected |\n|---|---|\n| A12V | 2 |\n"
        for i in range(1, 31)
    )
    audit = audit_table_phenotype_coverage({"variants": []}, source)
    assert audit["clinical_table_candidate_count"] == 30
    assert len(audit["clinical_table_candidates"]) == 20
    assert audit["candidates_truncated"]


def test_deterministic_success_audits_even_when_model_calls_would_fail(monkeypatch):
    extractor = ExpertExtractor(models=["test"])
    row = {
        "protein_notation": "p.Arg12Trp",
        "penetrance_data": {"total_carriers_observed": 7},
    }
    monkeypatch.setattr(
        extractor,
        "_do_attempt_extraction",
        lambda *args: ExtractionResult(
            pmid="1",
            success=True,
            model_used="deterministic-fixed-width-table-parser",
            extracted_data={"variants": [row]},
        ),
    )

    def no_api(*args, **kwargs):
        raise AssertionError("coverage must never call a model")

    monkeypatch.setattr(extractor, "call_llm_json", no_api)
    monkeypatch.setattr(extractor, "call_llm_json_with_status", no_api)
    result = extractor._attempt_extraction(
        Paper(pmid="1", gene_symbol="SCN5A", full_text=SOURCE), "test", SOURCE
    )
    assert result.success
    assert result.extracted_data["variants"][0]["penetrance_data"] == {
        "total_carriers_observed": 7
    }
    assert (
        result.extracted_data["extraction_metadata"]["table_phenotype_coverage"][
            "status"
        ]
        == "review_candidates"
    )


def test_pdf_form_feed_does_not_shift_source_line_coordinates():
    source = "Page one\n\fTable 2: Clinical phenotypes\n| variant | affected |\n"
    audit = audit_table_phenotype_coverage({"variants": []}, source)
    assert audit["clinical_table_candidates"][0]["source_line_1based"] == 2


def test_crlf_source_and_signal_truncation_are_explicit():
    source = (
        "preamble\r\nTable 1: Clinical phenotypes\r\n"
        "phenotype phenotypic symptomatic symptoms asymptomatic affected unaffected "
        "ECG QTc diagnosis diagnosed cardiac events follow-up\r\n"
    )
    audit = audit_table_phenotype_coverage({"variants": []}, source)
    candidate = audit["clinical_table_candidates"][0]
    assert candidate["source_line_1based"] == 2
    assert len(candidate["clinical_signals"]) == 12
    assert candidate["signals_truncated"]


def test_router_success_uses_original_source_coordinates_not_synthetic_tables(
    monkeypatch,
):
    extractor = ExpertExtractor(models=["test"])
    monkeypatch.setattr(
        extractor,
        "_do_attempt_extraction",
        lambda *args: ExtractionResult(
            pmid="1",
            success=True,
            model_used="router+test",
            extracted_data={"variants": []},
        ),
    )
    monkeypatch.setattr(
        extractor,
        "_augment_pdf_linearized_tables",
        lambda text: text + "\nTable 99: Clinical synthetic reconstruction\n",
    )
    result = extractor._attempt_extraction(
        Paper(pmid="1", full_text="original source"), "test", "prepared text"
    )
    audit = result.extracted_data["extraction_metadata"]["table_phenotype_coverage"]
    assert audit["clinical_table_candidate_count"] == 0
    import hashlib

    assert audit["source_sha256"] == hashlib.sha256(b"original source").hexdigest()
