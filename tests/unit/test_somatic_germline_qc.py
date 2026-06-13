"""Tests for the gold-free somatic-vs-germline QC validator.

This guards cold-start CANCER-gene runs (BRCA2, etc.): somatic tumor variants
are not heritable carriers, so they must be flaggable before they pollute
carrier/penetrance counts. Flag-only by default; drop only the clearly-somatic.
"""

import json

import pytest

from pipeline.somatic_germline_qc import (
    AMBIGUOUS,
    GERMLINE,
    SOMATIC,
    UNKNOWN,
    annotate_variants,
    classify_text,
    classify_variant,
)
from scripts.apply_somatic_germline_qc import _paper_context, _pmid


# --- classify_text ----------------------------------------------------------


def test_classify_text_clear_somatic():
    hit = classify_text(
        "Somatic mutation identified in tumor tissue by ctDNA sequencing"
    )
    assert hit.label == SOMATIC
    assert "somatic" in hit.somatic_terms
    assert hit.germline_terms == []


def test_classify_text_clear_germline():
    hit = classify_text(
        "Germline carrier; proband from a familial pedigree, blood sample"
    )
    assert hit.label == GERMLINE
    assert "germline" in hit.germline_terms
    assert hit.somatic_terms == []


def test_classify_text_ambiguous_when_both_present():
    hit = classify_text(
        "Germline carriers were screened; tumor biopsies showed loss of heterozygosity"
    )
    assert hit.label == AMBIGUOUS
    assert hit.somatic_terms and hit.germline_terms


def test_classify_text_unknown_when_neither():
    hit = classify_text("The variant was observed in 12 individuals in Table 2")
    assert hit.label == UNKNOWN
    assert hit.somatic_terms == [] and hit.germline_terms == []


def test_classify_text_handles_none_and_empty():
    assert classify_text(None).label == UNKNOWN
    assert classify_text("").label == UNKNOWN


# --- classify_variant reads all the record's text fields --------------------


def test_classify_variant_reads_individual_record_evidence():
    variant = {
        "source_location": "Table 2",
        "key_quotes": [],
        "individual_records": [
            {"evidence_sentence": "Tumor specimen sequencing revealed the variant."}
        ],
    }
    assert classify_variant(variant).label == SOMATIC


def test_classify_variant_reads_key_quotes_and_paper_context():
    variant = {
        "source_location": "Table 1",
        "key_quotes": ["heterozygous carrier in the family"],
    }
    assert classify_variant(variant).label == GERMLINE
    # paper-level context alone is enough to flag somatic
    bare = {"source_location": "Table 1"}
    assert (
        classify_variant(bare, paper_context="A tumor-only NGS panel study").label
        == SOMATIC
    )


# --- annotate_variants: summary + policies ----------------------------------


def _corpus():
    return [
        {"source_location": "germline carrier proband", "key_quotes": []},
        {"source_location": "somatic tumor biopsy", "key_quotes": []},
        {
            "source_location": "germline pedigree; tumor ctDNA",
            "key_quotes": [],
        },  # ambiguous
        {"source_location": "Table 3, 4 individuals", "key_quotes": []},  # unknown
    ]


def test_annotate_flag_policy_annotates_without_dropping():
    variants = _corpus()
    summary = annotate_variants(variants, policy="flag")
    assert len(variants) == 4  # nothing dropped
    assert all("somatic_germline_flag" in v for v in variants)
    assert summary["counts"] == {GERMLINE: 1, SOMATIC: 1, AMBIGUOUS: 1, UNKNOWN: 1}
    # contamination = somatic + ambiguous = 2/4
    assert summary["somatic_fraction"] == 0.5
    assert summary["dropped_somatic"] == 0


def test_annotate_drop_somatic_removes_only_clear_somatic():
    variants = _corpus()
    summary = annotate_variants(variants, policy="drop_somatic")
    assert summary["dropped_somatic"] == 1
    assert len(variants) == 3  # somatic removed; ambiguous + unknown + germline kept
    labels = {v["somatic_germline_flag"]["label"] for v in variants}
    assert SOMATIC not in labels
    assert AMBIGUOUS in labels and UNKNOWN in labels and GERMLINE in labels


def test_annotate_empty_corpus_is_zero_fraction():
    summary = annotate_variants([], policy="flag")
    assert summary["total"] == 0
    assert summary["somatic_fraction"] == 0.0


def test_annotate_rejects_unknown_policy():
    with pytest.raises(ValueError):
        annotate_variants([], policy="nuke")


# --- apply-script paper-level context (the BRCA2 run carries the abstract
#     in abstract_json/, not in the extraction payload) ----------------------


def test_apply_pmid_reads_paper_metadata():
    assert _pmid({"paper_metadata": {"pmid": "123"}}) == "123"
    assert _pmid({"pmid": "456"}) == "456"
    assert _pmid({}) == ""


def test_apply_paper_context_pulls_abstract_from_source_dir(tmp_path):
    (tmp_path / "abstract_json").mkdir()
    (tmp_path / "abstract_json" / "999.json").write_text(
        json.dumps({"abstract": "Somatic BRCA2 mutations from tumor ctDNA sequencing"})
    )
    payload = {"paper_metadata": {"pmid": "999", "title": "Paper 999"}, "variants": []}
    # with the source dir, the abstract's somatic signal is visible
    assert "somatic" in _paper_context(payload, tmp_path, "999").lower()
    # without it, only the placeholder title -> no signal (the bug this fixes)
    assert "somatic" not in _paper_context(payload, None, "999").lower()
