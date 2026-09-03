"""Pin the six generic fixes drawn from the tranche 01 baseline root cause.

Each fix is gene-agnostic and was found by tracing a paper-derived miss (or
a counted extra) through the stages the run left on disk. The tests name the
type specimen so the reason survives the code.
"""

from __future__ import annotations

from pathlib import Path

import pipeline.steps as steps
from benchmarks.codex_paper_eval.db_to_predictions import notation
from harvesting.content_validation import (
    has_substantive_supplement_content,
    is_abstract_reference_shell,
    validate_content_quality,
)
from pipeline.source_quality import (
    RICHER_SIBLING_MIN_BYTES,
    is_reusable_fulltext_source,
    is_usable_fulltext_source,
    prefer_richer_sibling,
)
from pipeline.table_router import (
    MarkdownTable,
    _infer_column_mapping_from_headers,
    _normalize_protein,
    enumerate_markdown_tables,
    parse_routed_table,
)
from pipeline.trust_gate import evaluate_fact
from utils.legacy_notation import promote_source_only_protein_identity
from utils.source_layers import infer_source_layer_from_text

FOLD_BEGIN = "<!-- GVF_FOLDED_SUPPLEMENTS_BEGIN -->"


# ---------------------------------------------------------------------------
# 1. A publisher access shell whose folded supplements carry the tables is a
#    usable, reusable source (RYR2 22222782: 235 supplement table rows).
# ---------------------------------------------------------------------------
def _springer_shell(with_supplement: bool) -> str:
    body = (
        "# MAIN TEXT\n\n## Postmortem genetic testing of RYR2\n\n### Abstract\n\n"
        + "Sudden unexplained death cases were screened for RYR2 variants. " * 12
        + "\n\n### Access this article\n\nBuy Now\n\n### References\n\n"
        + "\n".join(f"- Author {i}. Title {i}. Journal 20{i:02d}." for i in range(80))
        + "\n\n### Author information\n\nAuthors and Affiliations\n"
    )
    if with_supplement:
        rows = "\n".join(
            f"| Case {i} | p.Arg{100 + i}Trp | c.{300 + 3 * i}C>T | {i % 3 + 1} |"
            for i in range(40)
        )
        body += (
            f"\n\n{FOLD_BEGIN}\n# SUPPLEMENTAL FILE 1: 414_2011_658_MOESM1_ESM.doc\n\n"
            "| Case | Protein | cDNA | Carriers |\n|---|---|---|---|\n" + rows + "\n"
        )
    return body


def test_shell_with_folded_supplement_tables_is_usable_and_reusable(tmp_path):
    path = tmp_path / "22222782_FULL_CONTEXT.md"
    path.write_text(_springer_shell(with_supplement=True), encoding="utf-8")
    content = path.read_text(encoding="utf-8")
    assert is_abstract_reference_shell(content[:8192])  # the article body is absent
    assert has_substantive_supplement_content(content)
    assert is_usable_fulltext_source(path)
    assert is_reusable_fulltext_source(path)
    valid, _ = validate_content_quality(content)
    assert valid


def test_plain_shell_without_supplements_stays_refused(tmp_path):
    path = tmp_path / "12710526_FULL_CONTEXT.md"
    path.write_text(_springer_shell(with_supplement=False), encoding="utf-8")
    content = path.read_text(encoding="utf-8")
    assert not has_substantive_supplement_content(content)
    assert not is_usable_fulltext_source(path)
    valid, reason = validate_content_quality(content)
    assert not valid and "shell" in reason.lower()


def test_supplement_marker_with_negligible_text_is_not_enough():
    content = (
        _springer_shell(with_supplement=False)
        + f"\n{FOLD_BEGIN}\n# SUPPLEMENTAL FILE 1: x.doc\nshort\n"
    )
    assert not has_substantive_supplement_content(content)


# ---------------------------------------------------------------------------
# 2. A much richer CLEANED sibling is preferred over a worse FULL_CONTEXT, and
#    preprocessing does not overwrite a staged richer rendering (KCNQ1 21131640).
# ---------------------------------------------------------------------------
def _article(chars: int) -> str:
    return (
        "# MAIN TEXT\n\n## Title\n\n## Abstract\n\nSummary.\n\n## Methods\n\n"
        + ("Patients were genotyped for KCNQ1 variants. " * (chars // 45))
        + "\n\n## Results\n\nWe found p.Arg231Cys in three probands.\n"
    )


def test_prefer_richer_sibling_when_cleaned_is_much_larger(tmp_path):
    full = tmp_path / "21131640_FULL_CONTEXT.md"
    cleaned = tmp_path / "21131640_CLEANED.md"
    full.write_text(_article(6_000), encoding="utf-8")
    cleaned.write_text(_article(22_000), encoding="utf-8")
    assert prefer_richer_sibling(full) == cleaned


def test_prefer_richer_sibling_keeps_full_context_when_sizes_are_close(tmp_path):
    full = tmp_path / "1_FULL_CONTEXT.md"
    cleaned = tmp_path / "1_CLEANED.md"
    full.write_text(_article(20_000), encoding="utf-8")
    cleaned.write_text(_article(22_000), encoding="utf-8")
    assert prefer_richer_sibling(full) == full


def test_prefer_richer_sibling_refuses_junk_sibling(tmp_path):
    full = tmp_path / "2_FULL_CONTEXT.md"
    cleaned = tmp_path / "2_CLEANED.md"
    full.write_text(_article(6_000), encoding="utf-8")
    cleaned.write_text("Access denied. " * 2_000, encoding="utf-8")
    assert prefer_richer_sibling(full) == full


def test_corpus_cache_stages_richer_cleaned_sibling_as_full_context(tmp_path):
    corpus = tmp_path / "corpus"
    harvest = tmp_path / "run" / "pmc_fulltext"
    harvest.mkdir(parents=True)
    paper = corpus / "KCNQ1" / "21131640"
    paper.mkdir(parents=True)
    (paper / "21131640_FULL_CONTEXT.md").write_text(_article(6_000), encoding="utf-8")
    (paper / "21131640_CLEANED.md").write_text(_article(22_000), encoding="utf-8")

    recovered = steps._consolidate_from_corpus(["21131640"], harvest, "KCNQ1", corpus)

    assert recovered == {"21131640"}
    staged = (harvest / "21131640_FULL_CONTEXT.md").read_text(encoding="utf-8")
    assert len(staged) > 20_000, "the richer rendering must be what the run reads"


def test_preprocess_keeps_existing_richer_cleaned(tmp_path):
    cleaned = tmp_path / "1_CLEANED.md"
    cleaned.write_text(_article(22_000), encoding="utf-8")
    assert cleaned.stat().st_size >= RICHER_SIBLING_MIN_BYTES
    assert steps._existing_cleaned_is_richer(
        cleaned, "# ABSTRACT-ONLY\n\nshort abstract"
    )
    assert not steps._existing_cleaned_is_richer(cleaned, _article(21_000))
    assert not steps._existing_cleaned_is_richer(tmp_path / "missing_CLEANED.md", "x")


# ---------------------------------------------------------------------------
# 3. Source-only legacy PROTEIN spellings get an identity (SCN5A 16764707).
# ---------------------------------------------------------------------------
def test_promotes_single_letter_frameshift_with_stop_distance():
    variant = {
        "source_notation": "L860fsx89",
        "cdna_notation": None,
        "protein_notation": None,
    }
    assert promote_source_only_protein_identity(variant, "SCN5A") == "L860fsX"
    assert variant["protein_notation"] == "L860fsX"


def test_promotes_residue_position_amino_acid_insertion_as_legacy():
    for token, expected in (
        ("1795insD", "1795insD"),
        ("1570 insI", "1570insI"),
        ("277delS", "277delS"),
    ):
        variant = {"source_notation": token}
        assert promote_source_only_protein_identity(variant, "SCN5A") == expected
        assert variant["legacy_notation"] == expected
        assert not variant.get("cdna_notation")


def test_promotion_leaves_nucleotide_payloads_and_named_rows_alone():
    # A/C/G/T-only payload is ambiguous with a nucleotide indel: not ours.
    variant = {"source_notation": "1570insA"}
    assert promote_source_only_protein_identity(variant, "SCN5A") is None
    assert "legacy_notation" not in variant
    # A row that already has an identity is never rewritten.
    named = {"source_notation": "L860fsx89", "protein_notation": "p.Leu860fs"}
    assert promote_source_only_protein_identity(named, "SCN5A") is None
    assert named["protein_notation"] == "p.Leu860fs"
    # Prose-shaped tokens are not variants.
    assert (
        promote_source_only_protein_identity({"source_notation": "exon 3"}, "SCN5A")
        is None
    )


def test_projection_notation_falls_back_to_legacy_identity():
    assert notation(None, None, legacy="1795insD") == "1795insD"
    assert notation("p.Arg231Cys", None, legacy="1795insD") == "p.Arg231Cys"
    assert notation(None, None) == ""


# ---------------------------------------------------------------------------
# 4. "HGVSp" columns and p.(Xxx###Yyy) predicted-protein values (MYBPC3 33673806).
# ---------------------------------------------------------------------------
PANEL_TABLE = """
| gene_name | gt | classification | hgvsc | hgvsp | Case count |
|---|---|---|---|---|---|
| ACTC1 | HET | pathogenic | c.301G>A | p.(Glu101Lys) | 1 |
| MYBPC3 | HET | likely_pathogenic | c.2528_2536del | p.(Glu843_Arg845del) | 2 |
| MYBPC3 | HET | pathogenic | c.1039G>A | p.(Glu347Lys) | 1 |
| MYBPC3 | HET | pathogenic | c.2827C>T | p.(Arg943Ter) | 3 |
"""


def test_predicted_protein_parentheses_are_normalized():
    assert _normalize_protein("p.(Glu101Lys)") == "p.Glu101Lys"
    assert _normalize_protein("p.(Glu843_Arg845del)") == "p.Glu843_Arg845del"
    assert _normalize_protein("p.Glu101Lys") == "p.Glu101Lys"


def test_hgvsp_header_maps_to_protein_and_case_count_yields_carriers():
    table = enumerate_markdown_tables(PANEL_TABLE)[0]
    mapping = _infer_column_mapping_from_headers(table, target_gene="MYBPC3")
    assert mapping is not None
    assert mapping.get("protein") == 4, mapping
    assert mapping.get("cdna") == 3, mapping
    assert mapping.get("affected") == 5, mapping

    variants = parse_routed_table(table, mapping, "MYBPC3")
    by_protein = {v["protein_notation"]: v for v in variants}
    assert "p.Glu843_Arg845del" in by_protein, sorted(by_protein)
    row = by_protein["p.Glu843_Arg845del"]
    pen = row["penetrance_data"]
    # Closed case series: every counted case is an observed carrier.
    assert pen["affected_count"] == 2
    assert pen["total_carriers_observed"] == 2
    assert pen["unaffected_count"] is None
    prov = row["count_provenance"]
    assert prov["carriers_column_label"] == "Case count"
    assert prov["carriers_count_type"] == "per_variant_carrier"
    assert row["patients"]["column_ref"] == "Case count"
    # Rows of other genes stay out of a MYBPC3 extraction.
    assert "p.Glu101Lys" not in by_protein


def test_case_series_rule_does_not_fire_with_an_unaffected_column():
    text = """
| Variant | Protein | Affected (n) | Unaffected (n) |
|---|---|---|---|
| c.1A>G | p.Met1Val | 2 |  |
"""
    table = enumerate_markdown_tables(text)[0]
    mapping = {"cdna": 0, "protein": 1, "affected": 2, "unaffected": 3}
    variants = parse_routed_table(table, mapping, "KCNQ1")
    pen = variants[0]["penetrance_data"]
    assert pen["affected_count"] == 2
    assert pen["total_carriers_observed"] is None


def test_case_series_rule_ignores_family_and_allele_headers():
    text = """
| Variant | Protein | Families with variant |
|---|---|---|
| c.1A>G | p.Met1Val | 2 |
"""
    table = enumerate_markdown_tables(text)[0]
    mapping = {"cdna": 0, "protein": 1, "affected": 2}
    variants = parse_routed_table(table, mapping, "KCNQ1")
    assert variants[0]["penetrance_data"]["total_carriers_observed"] is None


# ---------------------------------------------------------------------------
# 5. A paper row whose notes mention ClinVar stays a paper row (SCN5A 28469501).
# ---------------------------------------------------------------------------
def test_notes_mentioning_clinvar_do_not_make_a_paper_row_a_linkage_row():
    layer = infer_source_layer_from_text(
        source_location="Results of Molecular Testing, Table 3",
        additional_notes="Previously reported in LQTS (ref 8); ClinVar conflicting/uncertain significance.",
    )
    assert layer == "llm_table"
    assert (
        infer_source_layer_from_text(
            source_location="Discussion", additional_notes="listed in PubTator"
        )
        == "llm_text"
    )


def test_ingest_markers_still_define_linkage_origin():
    assert (
        infer_source_layer_from_text(source_location="ClinVar (PMID citation)")
        == "clinvar"
    )
    assert (
        infer_source_layer_from_text(source_location="PubTator3 (text-mined)")
        == "pubtator"
    )
    assert infer_source_layer_from_text(extraction_source="clinvar") == "clinvar"


# ---------------------------------------------------------------------------
# 6. Genotype-frequency columns are population counts (SCN5A 23414114).
# ---------------------------------------------------------------------------
def test_genotype_frequency_label_is_a_population_count():
    reasons = evaluate_fact(
        {"carriers": 1}, provenance={"carriers_column_label": "Genotype frequency"}
    )
    assert "population_count" in reasons
    clean = evaluate_fact(
        {"carriers": 3}, provenance={"carriers_column_label": "No. of carriers"}
    )
    assert "population_count" not in clean


def _unused(
    _: Path,
) -> None:  # pragma: no cover - keeps Path import honest for tmp_path typing
    return None
