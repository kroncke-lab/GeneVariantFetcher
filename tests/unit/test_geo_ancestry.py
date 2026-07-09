"""Unit tests for deterministic ethnicity / cohort-origin extraction."""

from utils.geo_ancestry import (
    country_from_affiliation,
    enrich_ethnicity_origin,
    find_ethnicity,
    find_origin,
)


def test_find_ethnicity_prefers_specific_terms():
    assert find_ethnicity("one with Ashkenazi Jewish ancestry") == "Ashkenazi Jewish"
    assert find_ethnicity("an African American proband") == "African American"
    assert find_ethnicity("East Asian cohort") == "East Asian"
    assert find_ethnicity("no ancestry stated") is None


def test_find_origin_countries_and_regions():
    assert find_origin("families from Northern Italy") == "Northern Italy"
    assert find_origin("recruited in the USA") == "United States"
    assert find_origin("recruited in the U.S.") == "United States"
    assert find_origin("a UK biobank sample") == "United Kingdom"
    assert find_origin("a U.K. cohort") == "United Kingdom"
    assert find_origin("Slav peoples of Eastern Europe") == "Eastern Europe"
    assert find_origin("nothing geographic here") is None
    assert find_origin("the variant was reported to us by the family") is None


def test_country_from_affiliation_tail_and_variants():
    assert (
        country_from_affiliation("Dept of Genetics, University of Tokyo, Tokyo, Japan.")
        == "Japan"
    )
    assert (
        country_from_affiliation("Cedars-Sinai Medical Center, Los Angeles, CA, USA")
        == "United States"
    )
    assert (
        country_from_affiliation("Brigham and Women's Hospital, Boston, MA, U.S.")
        == "United States"
    )
    assert (
        country_from_affiliation("Peking Union Medical College, Beijing, P.R. China")
        == "China"
    )
    assert country_from_affiliation("") is None
    assert country_from_affiliation(None) is None


def test_enrich_prefers_stated_then_text_then_author():
    # Stated origin wins over author country.
    assert enrich_ethnicity_origin(
        stated_origin="Japan", author_affiliation="Somewhere, USA"
    ) == (None, "Japan")
    # Text-derived origin when nothing stated.
    eth, origin = enrich_ethnicity_origin(
        text_fields=["Ashkenazi Jewish families from Israel"]
    )
    assert eth == "Ashkenazi Jewish" and origin == "Israel"
    # Author-affiliation last resort is clearly tagged.
    _, origin = enrich_ethnicity_origin(
        text_fields=["breast cancer cohort"],
        author_affiliation="University of Melbourne, Melbourne, Australia",
    )
    assert origin == "Australia (inferred from author affiliation)"
    # Nothing anywhere.
    assert enrich_ethnicity_origin() == (None, None)
