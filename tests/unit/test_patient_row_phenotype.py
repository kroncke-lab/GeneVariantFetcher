from pipeline.count_provenance import strip_model_code_owned_phenotype_sources
from pipeline.patient_row_phenotype import derive_patient_row_phenotype_counts


SOURCE = """
We identified 5 living p.G357S mutation-positive individuals and 2 SCD cases.
Familial genetic analysis in 2 SCD cases showed that all had p.G357S.

**Supplementary Table 3. Basal phenotypic features of living mutation positive subjects**

|  |  |  |  |
| --- | --- | --- | --- |
| **Patient** | **Age** | **Previous symptoms** | **VA in basal test** |
| 1 | 20 | Syncope | No |
| 2 | 21 | No | Yes |
| 3 | 22 | No | No |
| 4 | 23 | No | - |
| 5 | 24 | Dizziness | - |
"""


def extraction(*, total=7, variants=1):
    rows = []
    for index in range(variants):
        rows.append(
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.G357S" if index == 0 else "p.R420W",
                "source_notation": "p.G357S" if index == 0 else "p.R420W",
                "patients": {"count": total},
                "penetrance_data": {
                    "total_carriers_observed": total,
                    "affected_count": None,
                    "unaffected_count": None,
                    "uncertain_count": None,
                },
                "count_provenance": {"carriers_count_type": "per_variant_carrier"},
            }
        )
    return {"variants": rows, "extraction_metadata": {}}


def test_derives_closed_partition_without_treating_missing_as_unaffected():
    result = derive_patient_row_phenotype_counts(extraction(), SOURCE)
    variant = result["variants"][0]
    assert variant["penetrance_data"] == {
        "total_carriers_observed": 7,
        "affected_count": 5,
        "unaffected_count": 1,
        "uncertain_count": 1,
    }
    assert variant["count_provenance"]["affected_count_type"] == (
        "derived_from_patient_rows"
    )
    assert variant["count_provenance"]["unaffected_count_type"] == (
        "derived_from_patient_rows"
    )
    assert variant["count_provenance"]["affected_source"] == "patient_row_phenotype_v2"
    assert (
        variant["count_provenance"]["unaffected_source"] == "patient_row_phenotype_v2"
    )
    assert variant["phenotype_derivation"]["table_affected"] == 3
    assert variant["phenotype_derivation"]["additional_affected"] == 2
    assert variant["phenotype_derivation"]["predicate_tallies"] == {
        "previous_symptoms_positive": 2,
        "va_in_basal_test_positive": 1,
        "any_selected_positive": 3,
        "all_selected_negative": 1,
        "missing_without_positive": 1,
    }


def test_unlinked_additional_carriers_remain_uncertain():
    source = SOURCE.replace(
        "Familial genetic analysis in 2 SCD cases showed that all had p.G357S.",
        "Two additional relatives were not clinically characterized.",
    ).replace("and 2 SCD cases", "and 2 additional relatives")
    result = derive_patient_row_phenotype_counts(extraction(), source)
    pdata = result["variants"][0]["penetrance_data"]
    assert pdata == {
        "total_carriers_observed": 7,
        "affected_count": 3,
        "unaffected_count": 1,
        "uncertain_count": 3,
    }


def test_requires_independent_confirmation_of_complete_table_total():
    source = SOURCE.replace("5 living", "several living")
    original = extraction()
    result = derive_patient_row_phenotype_counts(original, source)
    assert result["variants"][0]["penetrance_data"]["affected_count"] is None
    assert result["extraction_metadata"]["patient_row_phenotype_derivation"] == {
        "protocol_version": "patient_row_phenotype_v2",
        "attempted": True,
        "applied": False,
        "paper_variant_count": 1,
        "parsed_patient_table_count": 1,
        "applied_variant_count": 0,
        "newly_applied_variant_count": 0,
        "outcomes": [
            {
                "variant": "p.G357S",
                "status": "no_eligible_complete_patient_table",
                "eligible_table_count": 0,
            }
        ],
    }


def test_does_not_attribute_one_patient_table_across_multiple_variants():
    original = extraction(variants=2)
    result = derive_patient_row_phenotype_counts(original, SOURCE)
    assert all(
        row["penetrance_data"]["affected_count"] is None for row in result["variants"]
    )
    assert not result["extraction_metadata"]["patient_row_phenotype_derivation"][
        "applied"
    ]


def test_descending_patient_ids_and_generic_sheet_caption_use_single_variant_title():
    source = """
# Genealogy and clinical course of CPVT caused by RYR2 P2328S

A total of 5 mutation carriers were found; 3 subjects were eligible for clinical course analyses.

Sheet: Data

| Patient Id | CPVT syncope | ACA | SCD | Symptomatic at enrollment |
| --- | --- | --- | --- | --- |
| 3 | no | no | no | no |
| 2 | yes | no | no | yes |
| 1 | no | yes | no | yes |
"""
    original = extraction(total=None)
    original["paper_metadata"] = {
        "title": "Genealogy and clinical course of CPVT caused by RYR2 P2328S"
    }
    variant = original["variants"][0]
    variant["protein_notation"] = "p.Pro2328Ser"
    variant["source_notation"] = "RYR2 P2328S or p.(Pro2328Ser)"
    result = derive_patient_row_phenotype_counts(original, source)
    assert result["variants"][0]["penetrance_data"] == {
        "total_carriers_observed": 5,
        "affected_count": 2,
        "unaffected_count": 1,
        "uncertain_count": 2,
    }
    audit = result["variants"][0]["phenotype_derivation"]
    assert audit["table_total"] == 3
    assert audit["additional_uncertain"] == 2
    assert audit["predicate_tallies"] == {
        "cpvt_syncope_positive": 1,
        "aca_positive": 1,
        "symptomatic_at_enrollment_positive": 2,
        "any_selected_positive": 2,
        "all_selected_negative": 1,
        "missing_without_positive": 0,
    }


def test_variant_named_caption_can_bind_one_table_in_multi_variant_paper():
    source = SOURCE.replace(
        "**Supplementary Table 3. Basal phenotypic features of living mutation positive subjects**",
        "**Supplementary Table 3. p.G357S mutation-positive phenotypic features**",
    )
    result = derive_patient_row_phenotype_counts(extraction(variants=2), source)
    first, second = result["variants"]
    assert first["penetrance_data"]["affected_count"] == 5
    assert second["penetrance_data"]["affected_count"] is None
    assert (
        result["extraction_metadata"]["patient_row_phenotype_derivation"][
            "applied_variant_count"
        ]
        == 1
    )


def test_incomplete_model_derivation_is_not_counted_as_protocol_application():
    original = extraction()
    variant = original["variants"][0]
    variant["penetrance_data"]["affected_count"] = 3
    variant["phenotype_derivation"] = {
        "method": "derived_from_patient_rows",
        "complete_table": False,
    }
    result = derive_patient_row_phenotype_counts(original, SOURCE)
    audit = result["extraction_metadata"]["patient_row_phenotype_derivation"]
    assert not audit["applied"]
    assert audit["outcomes"] == [
        {"variant": "p.G357S", "status": "phenotype_already_populated"}
    ]


def test_model_cannot_forge_code_owned_phenotype_sources():
    payload = extraction()
    provenance = payload["variants"][0]["count_provenance"]
    provenance.update(
        {
            "affected_source": "patient_row_phenotype_v2",
            "unaffected_source": "source_bound_phenotype_v1",
            "carriers_source": "count_recovery",
        }
    )

    assert strip_model_code_owned_phenotype_sources(payload) == 2
    assert "affected_source" not in provenance
    assert "unaffected_source" not in provenance
    assert provenance["carriers_source"] == "count_recovery"
