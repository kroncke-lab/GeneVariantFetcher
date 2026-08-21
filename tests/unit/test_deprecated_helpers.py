"""Compatibility-window coverage for legacy importable helper names."""

import pytest

from utils import html_utils, pubmed_utils, variant_normalizer


@pytest.mark.parametrize(
    "call",
    [
        lambda: html_utils.extract_pmids_from_html("<p>PMID: 12345678</p>"),
        lambda: html_utils.extract_dois_from_html("<p>10.1000/example</p>"),
        lambda: html_utils.create_scraping_session(),
        lambda: html_utils.extract_pmids_from_json_results(
            {"results": [{"pmid": "12345678"}]}
        ),
        lambda: variant_normalizer.normalize_protein_variant("p.Ala561Val"),
        lambda: variant_normalizer.normalize_cdna_variant("1682C>T"),
        lambda: variant_normalizer.find_matching_variants(["A561V"], ["A561V"]),
    ],
)
def test_legacy_helpers_emit_deprecation_warning(call):
    with pytest.warns(DeprecationWarning, match="is deprecated"):
        call()


def test_query_pubmed_for_gene_warns_before_delegating(monkeypatch):
    monkeypatch.setattr(
        pubmed_utils,
        "query_pubmed_with_entrez",
        lambda *_args, **_kwargs: ["12345678"],
    )

    with pytest.warns(DeprecationWarning, match="query_pubmed_for_gene"):
        result = pubmed_utils.query_pubmed_for_gene("BRCA1")

    assert result == {"12345678"}


def test_validate_pmid_warns_before_existence_check(monkeypatch):
    monkeypatch.setattr(
        pubmed_utils,
        "fetch_paper_metadata",
        lambda _pmid: {"pmid": "12345678"},
    )

    with pytest.warns(DeprecationWarning, match="validate_pmid"):
        assert pubmed_utils.validate_pmid("12345678") is True
