"""Pin the v2 hardenings drawn from the tranche 01 candidate and its reviews.

Grok 4.6 and Gemini 3.1 Pro each produced a counterexample shape for the
tranche 01 fixes; the scorer pairing error was found in the locked score.
"""

from __future__ import annotations

from benchmarks.codex_paper_eval.run_eval import matches
from harvesting.content_validation import has_substantive_supplement_content
from pipeline.table_router import (
    _is_case_series_people_header,
    _normalize_protein,
    enumerate_markdown_tables,
    parse_routed_table,
)
from pipeline.trust_gate import evaluate_fact
from utils.legacy_notation import promote_source_only_protein_identity

FOLD_BEGIN = "<!-- GVF_FOLDED_SUPPLEMENTS_BEGIN -->"


# 1. A bare case-control "Cases" column is a cohort denominator, not a series.
def test_bare_cases_header_does_not_fire_case_series_rule():
    for label in ("Cases", "Patients", "Probands", "cases (%)"):
        assert not _is_case_series_people_header(label), label
    for label in (
        "Case count",
        "No. of patients",
        "Probands (n)",
        "n patients",
        "Carriers",
    ):
        assert _is_case_series_people_header(label), label


def test_case_control_cases_column_never_becomes_carriers():
    text = """
| rsID | Protein | Cases |
|---|---|---|
| rs1 | p.Arg190Trp | 200 |
"""
    table = enumerate_markdown_tables(text)[0]
    variants = parse_routed_table(table, {"protein": 1, "affected": 2}, "KCNH2")
    # The router either refuses the assay-cued table outright or, if it emits
    # the row, must not turn the cohort "Cases" total into a carrier count.
    for variant in variants:
        assert variant["penetrance_data"]["total_carriers_observed"] is None
        assert variant["count_provenance"]["affected_count_type"] != "case"


# 2. Fold chrome without tables or variants does not rescue a shell.
def _shell() -> str:
    return (
        "# MAIN TEXT\n\n## Title\n\n### Abstract\n\n"
        + "Cohort summary. " * 40
        + "\n\n### References\n\n"
        + "\n".join(f"- Ref {i}." for i in range(60))
    )


def test_fold_marker_with_prose_only_is_not_substantive():
    chrome = (
        _shell()
        + f"\n{FOLD_BEGIN}\n# SUPPLEMENTAL FILE 1: consent.pdf\n\n"
        + (
            "Participants provided written informed consent for the study procedures. "
            * 30
        )
    )
    assert not has_substantive_supplement_content(chrome)


def test_fold_with_variant_tokens_but_no_pipe_table_is_substantive():
    prose = (
        _shell()
        + f"\n{FOLD_BEGIN}\n# SUPPLEMENTAL FILE 1: methods.docx\n\n"
        + ("Variants detected: p.Arg190Trp, p.Gly628Ser, c.1128C>T and A561V. " * 20)
    )
    assert has_substantive_supplement_content(prose)


# 3. Nucleotide-letter payloads (including N) are not protein legacy.
def test_nucleotide_wildcard_payload_is_not_promoted():
    for token in ("1100delN", "1570insA", "1234delGGN"):
        assert (
            promote_source_only_protein_identity({"source_notation": token}, "SCN5A")
            is None
        )
    assert (
        promote_source_only_protein_identity({"source_notation": "1795insD"}, "SCN5A")
        == "1795insD"
    )


# 4. p.(?) and p.(=) are not identities.
def test_unknown_and_silent_protein_markers_are_not_unwrapped():
    assert _normalize_protein("p.(?)") is None
    assert _normalize_protein("p.(=)") is None
    assert _normalize_protein("p.(Arg943Ter)") == "p.Arg943Ter"
    assert _normalize_protein("p.(Y1136*)") == "p.Y1136*"


# 5. Only genotype-frequency/count labels are population counts.
def test_genotype_phenotype_label_is_not_a_population_count():
    assert "population_count" not in evaluate_fact(
        {"carriers": 3},
        provenance={"carriers_column_label": "Genotype-phenotype correlation"},
    )
    assert "population_count" in evaluate_fact(
        {"carriers": 3}, provenance={"carriers_column_label": "Genotype frequency"}
    )


# 6. Conflicting concrete cDNAs disable the codon-position bridges.
def test_same_codon_nonsense_and_deletion_with_different_cdna_do_not_match():
    assert not matches("c.3407_3409delACT", "p.Tyr1136* c.3408C>A", "MYBPC3")
    assert not matches("p.Tyr1136* c.3408C>A", "c.3407_3409delACT", "MYBPC3")
    # The same events still pair, in every spelling the scorer already accepted.
    assert matches("c.3408C>A", "p.Tyr1136* c.3408C>A", "MYBPC3")
    assert matches("c.3407_3409delACT", "p.Tyr1136del c.3407_3409delACT", "MYBPC3")
    assert matches("c.3407_3409delACT", "p.Tyr1136del", "MYBPC3")
    # A cDNA in-frame deletion still bridges to its codon's protein deletion
    # when only one side names the cDNA (c.4850 is codon 1617).
    assert matches("F1617del", "c.4850_4852del", "SCN5A")
    assert matches("c.693delCA", "c.692_693delCA", "KCNH2")
    assert matches("R376H", "c.1127G>A p.Arg376His", "KCNH2")
