from pipeline.aggregation import DataAggregator
from pipeline.extraction import _guard_penetrance_claims


def _variant(percentage, facts=None, age_points=None):
    return {
        "gene_symbol": "BMPR2",
        "protein_notation": "p.Arg491Trp",
        "source_notation": "R491W",
        "penetrance_data": {
            "total_carriers_observed": 10,
            "affected_count": 2,
            "penetrance_percentage": percentage,
            "age_dependent_penetrance": age_points or [],
        },
        "fact_provenance": facts or [],
    }


def test_unsourced_penetrance_percentage_is_cleared_and_audited():
    data = {"variants": [_variant(43.9)]}

    guarded = _guard_penetrance_claims(
        data,
        "The paper reports ten BMPR2 carriers and two affected relatives.",
    )

    variant = guarded["variants"][0]
    assert variant["penetrance_data"]["penetrance_percentage"] is None
    assert variant["penetrance_guard_flags"][0]["raw_value"] == 43.9
    assert guarded["extraction_metadata"]["penetrance_claims_guarded"] == 1


def test_explicit_variant_specific_penetrance_quote_is_retained():
    quote = "For BMPR2 R491W, penetrance was 20 percent."
    data = {
        "variants": [
            _variant(
                20,
                facts=[
                    {
                        "fact_type": "penetrance_percentage",
                        "fact_value": "20",
                        "evidence_quote": quote,
                    }
                ],
                age_points=[
                    {
                        "age_range": "by age 50",
                        "penetrance_percentage": 20,
                        "evidence_quote": quote,
                    }
                ],
            )
        ]
    }

    guarded = _guard_penetrance_claims(data, f"Results. {quote} End.")

    pdata = guarded["variants"][0]["penetrance_data"]
    assert pdata["penetrance_percentage"] == 20
    assert len(pdata["age_dependent_penetrance"]) == 1
    assert "penetrance_guard_flags" not in guarded["variants"][0]


def test_aggregation_never_derives_penetrance_from_raw_counts():
    aggregate = DataAggregator().calculate_aggregate_penetrance(
        {
            "individual_records": [],
            "penetrance_data_points": [
                {
                    "_source_pmid": "123",
                    "total_carriers_observed": 10,
                    "affected_count": 10,
                    "unaffected_count": 0,
                    "penetrance_percentage": 100,
                }
            ],
            "source_pmids": ["123"],
        }
    )

    assert aggregate["total_carriers"] == 10
    assert aggregate["penetrance_percentage"] is None
