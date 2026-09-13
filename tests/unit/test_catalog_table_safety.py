"""Variant annotations must not become patient observations in the regex lane."""

import pytest

from pipeline.extraction import ExpertExtractor


def parse(headers, rows):
    text = "\n".join(
        "| " + " | ".join(cells) + " |"
        for cells in [headers, ["---"] * len(headers), *rows]
    )
    return ExpertExtractor.__new__(ExpertExtractor)._parse_markdown_table_variants(
        text, "BRCA2"
    )


@pytest.mark.parametrize("condition", ["Cancer syndrome", "Cancer|Other syndrome"])
def test_clinvar_export_is_not_a_patient_table(condition):
    # The source export has no subject column. Its three somatic "clinical"
    # headers previously minted one affected carrier, even with nan cells;
    # an unescaped pipe in Condition(s) was not the cause.
    headers = [
        "Name",
        "Gene(s)",
        "Protein change",
        "Condition(s)",
        "Accession",
        "GRCh37Chromosome",
        "GRCh37Location",
        "GRCh38Chromosome",
        "GRCh38Location",
        "VariationID",
        "AlleleID(s)",
        "dbSNP ID",
        "Canonical SPDI",
        "Variant type",
        "Molecular consequence",
        "Germline classification",
        "Germline date last evaluated",
        "Germline review status",
        "Somatic clinical impact",
        "Somatic clinical impact date last evaluated",
        "Somatic clinical impact review status",
        "Oncogenicity classification",
        "Oncogenicity date last evaluated",
        "Oncogenicity review status",
    ]
    row = [
        "NM_000059.4(BRCA2):c.3128C>T (p.Ala1043Val)",
        "BRCA2",
        "A1043V",
        condition,
        "VCV001002452",
        "13",
        "32911620",
        "13",
        "32337483",
        "1002452",
        "995526",
        "rs1555283059",
        "NC_000013.11:32337482:C:T",
        "single nucleotide variant",
        "missense variant",
        "Conflicting classifications of pathogenicity",
        "2023-03-23",
        "criteria provided, conflicting classifications",
        "nan",
        "nan",
        "nan",
        "nan",
        "nan",
        "nan",
    ]
    assert parse(headers, [row]) == []


@pytest.mark.parametrize(
    "classification",
    [
        "Conflicting classifications of pathogenicity",
        "Likely benign",
        "Pathogenic",
        "Uncertain significance",
    ],
)
def test_annotation_identity_preserves_class_without_inventing_counts(classification):
    [variant] = parse(
        [
            "cDNA",
            "Protein",
            "Germline classification",
            "Somatic clinical impact",
            "VariationID",
        ],
        [["c.3128C>T", "A1043V", classification, "nan", "1002452"]],
    )
    assert variant["protein_notation"] == "A1043V"
    assert variant["clinical_significance"] == classification.lower()
    assert variant["penetrance_data"] == {
        "total_carriers_observed": None,
        "affected_count": None,
        "unaffected_count": None,
    }
    assert variant["patients"]["count"] is None


@pytest.mark.parametrize("metadata", ["nan", "Pathogenic", "normal"])
def test_subject_identity_does_not_turn_classification_into_phenotype(metadata):
    [variant] = parse(
        ["Patient ID", "Protein", "Clinical significance", "Somatic clinical impact"],
        [["P01", "A1043V", "Conflicting classifications of pathogenicity", metadata]],
    )
    assert (
        variant["clinical_significance"]
        == "conflicting classifications of pathogenicity"
    )
    assert variant["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": None,
        "unaffected_count": None,
    }


@pytest.mark.parametrize(
    "phenotype, expected",
    [("Breast cancer", (1, None)), ("unaffected", (0, 1)), ("nan", (None, None))],
)
def test_patient_phenotype_survives_neighboring_annotations(phenotype, expected):
    [variant] = parse(
        [
            "Patient ID",
            "Protein",
            "Phenotype",
            "Clinical significance",
            "gnomAD allele count",
        ],
        [
            [
                "P01",
                "A1043V",
                phenotype,
                "Conflicting classifications of pathogenicity",
                "500",
            ]
        ],
    )
    assert variant["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": expected[0],
        "unaffected_count": expected[1],
    }


def test_explicit_carrier_and_phenotype_counts_survive_annotations():
    [variant] = parse(
        [
            "cDNA",
            "Protein",
            "Carriers",
            "Affected",
            "Unaffected",
            "Germline classification",
            "gnomAD allele count",
        ],
        [
            [
                "c.3128C>T",
                "A1043V",
                "7",
                "2",
                "5",
                "Conflicting classifications of pathogenicity",
                "500",
            ]
        ],
    )
    assert variant["penetrance_data"] == {
        "total_carriers_observed": 7,
        "affected_count": 2,
        "unaffected_count": 5,
    }
    assert (
        variant["clinical_significance"]
        == "conflicting classifications of pathogenicity"
    )


@pytest.mark.parametrize("annotation", ["gnomAD allele count", "REVEL score"])
def test_annotation_without_subject_does_not_infer_one_carrier(annotation):
    [variant] = parse(
        ["cDNA", "Protein", "Clinical features", annotation],
        [["c.3128C>T", "A1043V", "Hereditary cancer", "18"]],
    )
    assert variant["penetrance_data"] == {
        "total_carriers_observed": None,
        "affected_count": None,
        "unaffected_count": None,
    }
    assert variant["clinical_significance"] == "uncertain"
