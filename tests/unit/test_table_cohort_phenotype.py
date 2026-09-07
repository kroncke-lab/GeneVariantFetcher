"""Gold-free tests for the deterministic table-cohort phenotype projection.

Fixtures use invented captions, headers and counts. They pin the contract:
a case-series count column becomes ``affected`` only when the table itself
names the cohort; a control table's count becomes ``unaffected``; family,
relative, allele, population, autopsy, case-control and clinical-column tables
refuse; model-authored rows and existing phenotype values are never touched;
the always-on phenotype guard keeps every stamped value.
"""

from __future__ import annotations

import copy

from pipeline.count_provenance import (
    TABLE_COHORT_PHENOTYPE_SOURCE,
    strip_model_code_owned_phenotype_sources,
)
from pipeline.phenotype_count_guard import apply_phenotype_count_guard
from pipeline.table_cohort_phenotype import (
    TIER_CAPTION,
    TIER_CAPTION_DISEASE,
    TIER_COLUMN,
    TIER_PAPER,
    classify_table_cohort,
    derive_table_cohort_phenotype_counts,
    expand_caption,
)

CASE_SOURCE = """
# MAIN TEXT

Table 2. Summary of putative LQT2-associated mutations in KCNH2

| Region | Nucleotide | Variant | Mutation Type | No. of patients |
|---|---|---|---|---|
| Exon 2 | c.94G>A | p.Ala32Thr | Missense | 3 |
| Exon 7 | c.1682C>T | p.Ala561Val | Missense | 1 |
"""

COMPENDIUM_SOURCE = """
### Table 1

Demographics and mutation yield for 9 Brugada syndrome genetic testing centers

|  | 1 | 2 |
|---|---|---|
| Total | 451 | 365 |

### Table 2

Control variants found in 2,600 reference alleles

| Exon | Nucleotide change | Variant | Number | Status |
|---|---|---|---|---|
| 2 | 52 C>T | R18W | 1 | Rare control |
| 2 | 100 C>T | R34C | 44 | Polymorphism |

### Table 4

Compendium of Brugada syndrome-associated SCN5A mutations

| Region | Nucleotide change | Coding effect | No. of unrelated individuals | Testing center |
|---|---|---|---|---|
| Exon 2 | 3 G>A | M1I* | 1 | 1 |
| Exon 28 | 5350 G>A | E1784K | 14 | 1, 2, 3 |
"""


def table_row(
    *,
    protein="p.Ala32Thr",
    carriers=3,
    affected=None,
    unaffected=None,
    source_table="Table 2. Summary of putative LQT2-associated mutations in KCNH2",
    headers=("Region", "Nucleotide", "Variant", "Mutation Type", "No. of patients"),
    count_label="No. of patients",
    parser="markdown_table",
    provenance_extra=None,
    **overrides,
):
    provenance = {
        "carriers_column_label": count_label,
        "carriers_count_type": "per_variant_carrier",
        "affected_column_label": None,
        "affected_count_type": None,
        "unaffected_column_label": None,
        "unaffected_count_type": None,
    }
    provenance.update(provenance_extra or {})
    row = {
        "gene_symbol": "KCNH2",
        "protein_notation": protein,
        "cdna_notation": None,
        "patients": {
            "count": carriers,
            "phenotype": None,
            "source_ref": source_table,
            "row_ordinal": 1,
            "column_ref": count_label,
            "locator_extra": {"parser": parser},
        },
        "penetrance_data": {
            "total_carriers_observed": carriers,
            "affected_count": affected,
            "unaffected_count": unaffected,
        },
        "source_location": source_table,
        "source_table": source_table,
        "source_table_caption": source_table,
        "source_table_headers": list(headers),
        "source_row": "1",
        "source_column": count_label,
        "locator_extra": {"parser": parser},
        "additional_notes": "Parsed via deterministic table parser",
        "count_provenance": provenance,
    }
    row.update(overrides)
    return row


def extraction(*rows):
    return {"variants": list(rows), "extraction_metadata": {}}


def derive(data, source, **kwargs):
    kwargs.setdefault("enabled", True)
    kwargs.setdefault("allow_paper_tier", False)
    return derive_table_cohort_phenotype_counts(data, source, **kwargs)


def test_case_series_count_column_becomes_affected_and_survives_guard():
    data = extraction(table_row())
    before = copy.deepcopy(data)
    result = derive(data, CASE_SOURCE, gene_symbol="KCNH2")
    assert data == before, "input must not be mutated"
    variant = result["variants"][0]
    assert variant["penetrance_data"]["affected_count"] == 3
    assert variant["penetrance_data"]["unaffected_count"] is None
    provenance = variant["count_provenance"]
    assert provenance["affected_count_type"] == "case"
    assert provenance["affected_source"] == TABLE_COHORT_PHENOTYPE_SOURCE
    assert "No. of patients" in provenance["affected_column_label"]
    derivation = variant["phenotype_derivation"]
    assert derivation["method"] == "derived_from_table_cohort"
    assert derivation["cohort_role"] == "case"
    assert derivation["tier"] == TIER_COLUMN
    facts = [
        f for f in variant["fact_provenance"] if f["fact_type"] == "affected_count"
    ]
    assert facts and facts[0]["fact_value"] == "3"
    meta = result["extraction_metadata"]["table_cohort_phenotype_derivation"]
    assert meta["applied_variant_count"] == 1
    assert meta["outcomes"][0]["status"] == "applied"

    summary = apply_phenotype_count_guard(result["variants"])
    assert summary.cleared == 0
    assert result["variants"][0]["penetrance_data"]["affected_count"] == 3


def test_parser_case_copy_is_stamped_instead_of_cleared():
    """The parser sometimes copies the count onto affected with no provenance.

    Without the stamp the guard clears it as ``copied_carriers_onto_affected``;
    when the table classifies as a case series the copy is the source's claim.
    """
    row = table_row(carriers=5, affected=5)
    unstamped = extraction(copy.deepcopy(row))
    cleared = apply_phenotype_count_guard(unstamped["variants"])
    assert cleared.cleared == 1
    result = derive(extraction(row), CASE_SOURCE)
    variant = result["variants"][0]
    assert (
        result["extraction_metadata"]["table_cohort_phenotype_derivation"][
            "stamped_variant_count"
        ]
        == 1
    )
    assert variant["count_provenance"]["affected_count_type"] == "case"
    assert apply_phenotype_count_guard(result["variants"]).cleared == 0
    assert variant["penetrance_data"]["affected_count"] == 5


def test_bare_table_label_resolves_to_descriptive_caption():
    caption = expand_caption("Table 4", COMPENDIUM_SOURCE)
    assert (
        caption == "Table 4 Compendium of Brugada syndrome-associated SCN5A mutations"
    )
    assert expand_caption("Table 2", COMPENDIUM_SOURCE).endswith(
        "Control variants found in 2,600 reference alleles"
    )
    # "Table 1" must not swallow "Table 10." or a body sentence.
    assert expand_caption("Table 1", COMPENDIUM_SOURCE).startswith(
        "Table 1 Demographics"
    )


def test_wrapped_caption_continues_but_header_fragments_do_not():
    wrapped = (
        "SUPPLEMENTAL TABLES\n\n"
        "Table S1: List of Mutations by Coding Effect, Location, and Frequency in 406\n"
        "LQT3 Patients. The different regions of the channel were defined as\n"
        "the coding sequence involving amino acid residues.\n"
    )
    assert expand_caption(
        "Table S1: List of Mutations by Coding Effect, Location, and Frequency in 406",
        wrapped,
    ) == (
        "Table S1: List of Mutations by Coding Effect, Location, and Frequency "
        "in 406 LQT3 Patients."
    )
    fixed_width = (
        " Table 2. Included SCN5A Mutations and Variants\n"
        " History of Cardiac\n"
        " Nucleotide Aborted Cardiac Events Predicted as\n"
    )
    assert (
        expand_caption("Table 2. Included SCN5A Mutations and Variants", fixed_width)
        == "Table 2. Included SCN5A Mutations and Variants"
    )


def test_bare_label_resolves_inline_caption_before_classification():
    source = (
        "Table 20. Mutations in LQT2 patients\n\n"
        "Table 2 shows the patient counts.\n\n"
        "Table 2. Variants in patients and controls\n\n"
        "| Variant | No. of patients |\n|---|---|\n| p.Ala32Thr | 3 |\n"
    )
    assert expand_caption("Table 2", source) == (
        "Table 2. Variants in patients and controls"
    )
    result = derive(extraction(table_row(source_table="Table 2")), source)
    assert result["variants"][0]["penetrance_data"]["affected_count"] is None
    assert "phenotype_derivation" not in result["variants"][0]


def test_disease_caption_with_people_column_is_a_case_series():
    row = table_row(
        protein="E1784K",
        carriers=14,
        source_table="Table 4",
        headers=(
            "Region",
            "Nucleotide change",
            "Coding effect",
            "No. of unrelated individuals",
            "Testing center",
        ),
        count_label="No. of unrelated individuals",
        gene_symbol="SCN5A",
    )
    result = derive(extraction(row), COMPENDIUM_SOURCE, gene_symbol="SCN5A")
    variant = result["variants"][0]
    assert variant["penetrance_data"]["affected_count"] == 14
    assert variant["phenotype_derivation"]["tier"] == TIER_CAPTION_DISEASE
    assert "Brugada" in variant["phenotype_derivation"]["evidence_quote"]


def test_disease_caption_with_bare_number_column_refuses():
    cohort = classify_table_cohort(
        "Table 3. BrS-associated SCN5A mutations", "Number", ["Variant", "Number"]
    )
    assert cohort.role is None
    assert cohort.reason == "disease_caption_without_people_noun"


def test_caption_naming_disease_patients_with_count_column():
    cohort = classify_table_cohort(
        "Table S1: List of Mutations by Coding Effect, Location, and Frequency "
        "in 406 LQT3 Patients.",
        "COUNT",
        [],
    )
    assert cohort.role == "case"
    assert cohort.tier == TIER_CAPTION


def test_run_disease_phrase_counts_as_disease_evidence():
    cohort = classify_table_cohort(
        "Table 2. BMPR2 variants in the registry of pulmonary hypertension probands",
        "n",
        [],
        run_disease=__import__(
            "pipeline.table_cohort_phenotype", fromlist=["disease_pattern"]
        ).disease_pattern("Pulmonary arterial hypertension"),
    )
    assert cohort.role == "case"


def test_control_people_table_closes_partition_and_guard_keeps_the_zero():
    source = (
        "Table 3. Variants identified in 500 healthy control subjects\n\n"
        "| Variant | Number |\n|---|---|\n| R18W | 4 |\n"
    )
    row = table_row(
        protein="R18W",
        carriers=4,
        source_table="Table 3. Variants identified in 500 healthy control subjects",
        headers=("Variant", "Number"),
        count_label="Number",
    )
    result = derive(extraction(row), source)
    variant = result["variants"][0]
    assert variant["penetrance_data"]["unaffected_count"] == 4
    assert variant["penetrance_data"]["affected_count"] == 0
    provenance = variant["count_provenance"]
    assert provenance["unaffected_count_type"] == "unaffected_control"
    assert provenance["affected_count_type"] == "control"
    assert provenance["affected_source"] == TABLE_COHORT_PHENOTYPE_SOURCE
    assert variant["patients"]["phenotype"] == "unaffected control"
    assert apply_phenotype_count_guard(result["variants"]).cleared == 0
    assert result["variants"][0]["penetrance_data"]["affected_count"] == 0


def test_control_tables_must_pass_caption_header_and_people_count_checks():
    for caption, label, headers in (
        (
            "Table 1. Functional assays in control cells",
            "Number",
            ("Variant", "Number"),
        ),
        (
            "Table 2. Variants in cases and controls",
            "Controls (n)",
            ("Variant", "Cases (n)", "Controls (n)"),
        ),
        (
            "Table 3. Variants in healthy control subjects",
            "Number",
            ("Variant", "Number", "Affected", "Unaffected"),
        ),
        ("Table 4. BrS and control variants", "Number", ("Variant", "Number")),
        (
            "Table 5. Variants in healthy control subjects",
            "Positive tests",
            ("Variant", "Positive tests"),
        ),
        (
            "Table 6. Rare variants",
            "Controls (n)",
            ("Variant", "Cases (n)", "Controls (n)"),
        ),
        (
            "Table 7. Healthy control subjects and their relatives",
            "Number",
            ("Variant", "Number"),
        ),
    ):
        row = table_row(source_table=caption, count_label=label, headers=headers)
        result = derive(extraction(row), caption)
        variant = result["variants"][0]
        assert variant["penetrance_data"] == row["penetrance_data"], caption
        assert "phenotype_derivation" not in variant, caption
        assert apply_phenotype_count_guard(result["variants"]).cleared == 0


def test_model_declared_control_zero_is_still_unsourced():
    row = table_row(
        carriers=4,
        affected=0,
        unaffected=4,
        parser="",
        source_layer="llm_table",
        provenance_extra={"affected_count_type": "control"},
    )
    row["locator_extra"] = {}
    row["patients"]["locator_extra"] = {}
    row["additional_notes"] = ""
    summary = apply_phenotype_count_guard([row])
    assert {a["reason"] for a in summary.annotations} == {"unsourced_zero_affected"}


def test_control_allele_table_is_left_exactly_as_parsed():
    row = table_row(
        protein="R34C",
        carriers=44,
        affected=0,
        unaffected=44,
        source_table="Table 2",
        headers=("Exon", "Nucleotide change", "Variant", "Number", "Status"),
        count_label="Number",
    )
    result = derive(extraction(copy.deepcopy(row)), COMPENDIUM_SOURCE)
    variant = result["variants"][0]
    assert variant["penetrance_data"] == row["penetrance_data"]
    assert "phenotype_derivation" not in variant
    assert variant["count_provenance"].get("unaffected_source") is None
    outcome = result["extraction_metadata"]["table_cohort_phenotype_derivation"][
        "outcomes"
    ][0]
    assert outcome["status"] == "not_derived:control_count_unit_is_alleles"
    # Pre-existing behaviour is untouched: the guard drops the unsourced zero
    # and keeps the parsed unaffected count.
    apply_phenotype_count_guard(result["variants"])
    assert variant["penetrance_data"]["affected_count"] is None
    assert variant["penetrance_data"]["unaffected_count"] == 44


def test_per_person_proband_rows_and_relative_tables_refuse():
    caption = (
        "Table 2 Primary symptom of proband, mutation type, novelty of mutation, "
        "associated asymptomatic/symptomatic RyR2 variant-carrying relatives (n)"
    )
    implicit = classify_table_cohort(
        caption, "implicit one carrier per clinical row", ["Proband", "Mutation"]
    )
    assert implicit.role is None
    assert implicit.reason == "per_person_clinical_row"
    relatives = classify_table_cohort(
        "Table 3. CPVT probands and their relatives", "No. of carriers", []
    )
    assert relatives.role is None
    assert relatives.reason.startswith("count_column_excluded:carrier")
    family = classify_table_cohort(
        "Table 3. LQT1 mutations identified in 40 families", "No. of patients", []
    )
    assert family.role is None
    assert family.reason.startswith("caption_excluded:famil")


def test_pooled_case_and_control_count_label_refuses():
    """A converter can join "BrS (2111) | LQT (2888) | Control (8975)" into one
    label; a count that pools cases with controls is neither class."""
    for label in ("BrS + LQT + Control", "Cases/Controls", "patients and controls (n)"):
        cohort = classify_table_cohort(
            "Supplemental Table 1: Properties of SCN5A nsSNVs", label, []
        )
        assert cohort.role is None, label
        assert cohort.reason == "count_column_mixes_cases_and_controls", label
    plain = classify_table_cohort(
        "Supplemental Table 1: Properties of SCN5A nsSNVs", "Control (8975)", []
    )
    assert plain.role == "control"


def test_case_control_and_clinical_columns_refuse():
    two_arm = classify_table_cohort(
        "Table 2. Rare variants in BrS cases and controls",
        "Cases",
        ["Variant", "Cases", "Controls"],
    )
    assert two_arm.role is None
    assert two_arm.reason == "caption_mixes_cases_and_controls"
    columns_only = classify_table_cohort(
        "Table 2. Rare SCN5A variants",
        "Cases (n)",
        ["Variant", "Cases (n)", "Controls (n)"],
    )
    assert columns_only.role is None
    assert columns_only.reason == "case_and_control_columns_present"
    clinical = classify_table_cohort(
        "Table 2. LQT2 mutations in 30 patients",
        "No. of patients",
        ["Variant", "No. of patients", "Symptoms", "QTc"],
    )
    assert clinical.role is None
    assert clinical.reason.startswith("clinical_column_present")


def test_population_autopsy_and_literature_captions_refuse():
    for caption in (
        "Table 2. Rare variants in 5,000 UK Biobank population participants",
        "Table 3. SCN5A variants in 42 autopsy cases of sudden unexplained death",
        "Table 4. Previously published LQT1 mutations in patients",
        "Table 1. Clinical characteristics of patients with KCNH2 mutations",
        "Table 5. Functional expression of LQT2 mutations in HEK293 cells",
        "Table 4. Gene mutations in Chinese CTEPH patients and PE without PH patients.",
        "Table 2. SCN5A variants in BrS patients compared with AF patients",
    ):
        cohort = classify_table_cohort(caption, "No. of patients", [])
        assert cohort.role is None, caption
        assert cohort.reason.startswith("caption_excluded:"), caption
    allele = classify_table_cohort(
        "Table 2. Variants in 300 LQTS patients", "Allele count", []
    )
    assert allele.role is None
    assert allele.reason.startswith("count_column_excluded:allele")


def test_model_authored_rows_and_existing_phenotypes_are_untouched():
    model_row = {
        "gene_symbol": "KCNH2",
        "protein_notation": "p.Ala561Val",
        "patients": {"count": 6},
        "penetrance_data": {
            "total_carriers_observed": 6,
            "affected_count": None,
            "unaffected_count": None,
        },
        "source_location": "Table 2. Summary of putative LQT2-associated mutations in KCNH2",
        "source_layer": "llm_table",
        "count_provenance": {
            "carriers_column_label": "No. of patients",
            "carriers_count_type": "per_variant_carrier",
        },
    }
    populated = table_row(carriers=7, affected=2, unaffected=5)
    result = derive(extraction(model_row, populated), CASE_SOURCE)
    assert result["variants"][0]["penetrance_data"]["affected_count"] is None
    assert result["variants"][1]["penetrance_data"] == populated["penetrance_data"]
    statuses = [
        o["status"]
        for o in result["extraction_metadata"]["table_cohort_phenotype_derivation"][
            "outcomes"
        ]
    ]
    assert statuses == ["model_authored_row", "phenotype_already_populated"]


def test_non_per_variant_carrier_roles_refuse():
    row = table_row(
        carriers=12,
        count_label="No. of families",
        provenance_extra={"carriers_count_type": "family_count"},
    )
    result = derive(extraction(row), CASE_SOURCE)
    assert result["variants"][0]["penetrance_data"]["affected_count"] is None
    assert (
        result["extraction_metadata"]["table_cohort_phenotype_derivation"]["outcomes"][
            0
        ]["status"]
        == "carrier_type_family_count"
    )


def test_model_stamp_is_scrubbed_before_derivation():
    data = extraction(
        table_row(
            affected=3,
            provenance_extra={
                "affected_count_type": "case",
                "affected_source": TABLE_COHORT_PHENOTYPE_SOURCE,
            },
        )
    )
    assert strip_model_code_owned_phenotype_sources(data) == 1
    assert "affected_source" not in data["variants"][0]["count_provenance"]


def test_disabled_setting_is_a_recorded_noop():
    data = extraction(table_row())
    result = derive_table_cohort_phenotype_counts(data, CASE_SOURCE, enabled=False)
    assert result["variants"][0]["penetrance_data"]["affected_count"] is None
    meta = result["extraction_metadata"]["table_cohort_phenotype_derivation"]
    assert meta == {
        "protocol_version": TABLE_COHORT_PHENOTYPE_SOURCE,
        "attempted": False,
        "applied": False,
        "reason": "disabled_by_setting",
    }


def test_strip_replay_preserves_counts_that_preceded_projection():
    from scripts.replay_table_cohort_phenotype import derived_records, patch_paper

    for row, expected in (
        (table_row(carriers=1, affected=1), {"affected": 1}),
        (table_row(carriers=3, affected=3), {"affected": None}),
        (
            table_row(
                carriers=4,
                affected=0,
                unaffected=4,
                source_table="Table 3. Variants in healthy control subjects",
                count_label="Number",
                headers=("Variant", "Number"),
            ),
            {"affected": None, "unaffected": 4},
        ),
    ):
        result = derive(extraction(row), row["source_table"])
        records, _ = derived_records(
            result, "", gene="KCNH2", disease=None, mode="strip", paper_tier=False
        )
        assert len(records) == 1
        assert records[0]["targets"] == expected
        current = result["variants"][0]["penetrance_data"]
        paper = {
            "variants": [
                {
                    "variant": row["protein_notation"],
                    "affected": current["affected_count"],
                    "unaffected": current["unaffected_count"],
                }
            ]
        }
        patch_paper(paper, records, "strip")
        assert {f: paper["variants"][0][f] for f in expected} == expected


def test_strip_replay_refuses_legacy_stamps_without_previous_counts():
    import pytest
    from scripts.replay_table_cohort_phenotype import derived_records

    result = derive(extraction(table_row()), CASE_SOURCE)
    del result["variants"][0]["phenotype_derivation"][
        "guarded_counts_without_projection"
    ]
    with pytest.raises(ValueError, match="without audited pre-projection"):
        derived_records(
            result, "", gene="KCNH2", disease=None, mode="strip", paper_tier=False
        )


def test_replay_uses_scored_reference_assignment_for_count_audit():
    from scripts.replay_table_cohort_phenotype import scored_gold_rows

    gold = [
        {"variant": "A32T", "affected": 1},
        {"variant": "A32T", "affected": 2},
        {"variant": "A32del", "affected": 3},
    ]
    score = {
        "matched_variants": [
            {"predicted": "p.Ala32Thr", "gold": "A32T"},
            {"predicted": "c.94G>A", "gold": "A32T"},
            {"predicted": "p.Ala32del", "gold": "A32del"},
        ]
    }
    assigned = scored_gold_rows(score, gold)
    assert assigned["p.Ala32Thr"]["affected"] == 1
    assert assigned["c.94G>A"]["affected"] == 2
    assert assigned["p.Ala32del"]["affected"] == 3
    assert "A32T" not in assigned  # an unscored spelling must not reuse gold


def test_paper_ascertainment_tier_is_off_unless_enabled():
    source = (
        "Genotype-Phenotype Correlation of SCN5A Mutation for Probands With Brugada "
        "Syndrome: A Japanese Multicenter Registry\n\n"
        "Table 2. Included SCN5A Mutations and Variants\n\n"
        "| Nucleotide | Coding Effect | n |\n|---|---|---|\n| 163C>T | Q55X | 1 |\n"
    )
    row = table_row(
        protein="Q55X",
        carriers=1,
        source_table="Table 2. Included SCN5A Mutations and Variants",
        headers=("Nucleotide", "Coding Effect", "n"),
        count_label="n",
        gene_symbol="SCN5A",
    )
    title = (
        "Genotype-Phenotype Correlation of SCN5A Mutation for Probands With "
        "Brugada Syndrome: A Japanese Multicenter Registry"
    )
    off = derive(extraction(copy.deepcopy(row)), source, title=title)
    assert off["variants"][0]["penetrance_data"]["affected_count"] is None
    assert (
        off["extraction_metadata"]["table_cohort_phenotype_derivation"]["outcomes"][0][
            "status"
        ]
        == "not_derived:no_cohort_evidence"
    )
    on = derive(
        extraction(copy.deepcopy(row)), source, title=title, allow_paper_tier=True
    )
    variant = on["variants"][0]
    assert variant["penetrance_data"]["affected_count"] == 1
    assert variant["phenotype_derivation"]["tier"] == TIER_PAPER
    assert "Probands With" in variant["phenotype_derivation"]["evidence_quote"]


def test_extraction_success_boundary_applies_projection(monkeypatch):
    """The hook runs for every successful extraction shape, after the audited
    patient-row lane and before the persist-site guard."""
    from pipeline.extraction import ExpertExtractor
    from utils.models import ExtractionResult, Paper

    extractor = ExpertExtractor(models=["gpt-4"])
    row = table_row(carriers=3)
    monkeypatch.setattr(
        extractor,
        "_do_attempt_extraction",
        lambda paper, model, prepared_full_text=None, estimated_variants=None: (
            ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=extraction(row),
                model_used="deterministic-table-parser",
            )
        ),
    )
    paper = Paper(pmid="1", title="Test", full_text=CASE_SOURCE, gene_symbol="KCNH2")
    result = extractor._attempt_extraction(paper, "gpt-4", CASE_SOURCE)
    variant = result.extracted_data["variants"][0]
    assert variant["penetrance_data"]["affected_count"] == 3
    assert (
        variant["count_provenance"]["affected_source"] == TABLE_COHORT_PHENOTYPE_SOURCE
    )
    meta = result.extracted_data["extraction_metadata"]
    assert meta["table_cohort_phenotype_derivation"]["applied_variant_count"] == 1
    assert meta["patient_row_phenotype_derivation"]["attempted"] is True


def test_router_rows_resolve_caption_by_table_id():
    source = (
        "Table 2. Mutations identified in 50 LQT1 patients\n\n"
        "| Variant | No. of subjects |\n|---|---|\n| p.Gly168Arg | 4 |\n"
    )
    row = {
        "gene_symbol": "KCNQ1",
        "protein_notation": "p.Gly168Arg",
        "patients": {
            "count": 4,
            "phenotype": None,
            "source_ref": "Table T1",
            "row_ordinal": 1,
            "column_ref": "No. of subjects",
            "locator_extra": {"parser": "table_router", "table_id": "T1"},
        },
        "penetrance_data": {
            "total_carriers_observed": 4,
            "affected_count": None,
            "unaffected_count": None,
        },
        "source_location": "Table T1, row 1 (router+deterministic)",
        "source_ref": "Table T1",
        "locator_extra": {"parser": "table_router", "table_id": "T1"},
        "count_provenance": {
            "carriers_column_label": "No. of subjects",
            "carriers_count_type": "per_variant_carrier",
        },
    }
    result = derive(extraction(row), source, gene_symbol="KCNQ1")
    variant = result["variants"][0]
    assert variant["penetrance_data"]["affected_count"] == 4
    assert variant["phenotype_derivation"]["source_table"].startswith(
        "Table 2. Mutations identified in 50 LQT1 patients"
    )
