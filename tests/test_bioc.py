"""Tests for the BioC client (PubTator3 + PMC OA).

Uses real NCBI endpoints — requires network access.
"""

import pytest

from gene_literature.bioc_client import BioCClient

# Well-known KCNH2 / ion-channel papers with PubTator3 coverage
TEST_PMIDS = ["24667783", "19841300", "20173333"]


@pytest.fixture(scope="module")
def client():
    return BioCClient()


# ------------------------------------------------------------------
# Section classifier (offline, fast)
# ------------------------------------------------------------------


class TestClassifySection:
    def test_title(self):
        assert BioCClient.classify_section("TITLE") == "TITLE"
        assert BioCClient.classify_section("article.title") == "TITLE"

    def test_abstract(self):
        assert BioCClient.classify_section("ABSTRACT") == "ABSTRACT"
        assert BioCClient.classify_section("abstract") == "ABSTRACT"

    def test_methods(self):
        assert BioCClient.classify_section("METHODS") == "METHODS"
        assert BioCClient.classify_section("Materials and Methods") == "METHODS"
        assert BioCClient.classify_section("Experimental Procedures") == "METHODS"
        assert BioCClient.classify_section("Genotyping") == "METHODS"
        assert BioCClient.classify_section("Patients and Methods") == "METHODS"

    def test_results(self):
        assert BioCClient.classify_section("RESULTS") == "RESULTS"
        assert BioCClient.classify_section("Results and Discussion") == "RESULTS"
        assert BioCClient.classify_section("Functional Characterization") == "RESULTS"
        assert BioCClient.classify_section("Electrophysiology Results") == "RESULTS"
        assert (
            BioCClient.classify_section("Genotype-Phenotype Correlations") == "RESULTS"
        )

    def test_table(self):
        assert BioCClient.classify_section("TABLE") == "TABLE"
        assert BioCClient.classify_section("Table 1") == "TABLE"
        assert BioCClient.classify_section("Supplementary Table") == "TABLE"

    def test_other(self):
        assert BioCClient.classify_section("Acknowledgements") == "OTHER"
        assert BioCClient.classify_section("") == "OTHER"
        assert BioCClient.classify_section("Competing Interests") == "OTHER"

    def test_references(self):
        assert BioCClient.classify_section("REFERENCES") == "REFERENCES"
        assert BioCClient.classify_section("Bibliography") == "REFERENCES"


# ------------------------------------------------------------------
# PubTator3 integration tests
# ------------------------------------------------------------------


@pytest.mark.requires_network
class TestPubTator3:
    def test_basic_fetch(self, client):
        """PubTator3 returns title, abstract, and annotations."""
        result = client.get_pubtator_annotations("24667783")
        assert result["title"], "Expected non-empty title"
        assert result["abstract"], "Expected non-empty abstract"
        assert isinstance(result["annotations"], list)

    def test_annotations_have_structure(self, client):
        """Each annotation has the expected keys."""
        result = client.get_pubtator_annotations("24667783")
        assert len(result["annotations"]) > 0, "Expected at least one annotation"
        for ann in result["annotations"]:
            assert "type" in ann
            assert "text" in ann
            assert "start" in ann
            assert "end" in ann
            assert "identifiers" in ann

    def test_annotations_contain_gene_or_disease(self, client):
        """PubTator3 returns at least Gene or Disease annotations for LQTS papers.

        Note: PubTator3 does not always produce Mutation annotations for
        all papers.  Gene and Disease are the most reliably annotated types
        for the KCNH2/ion-channel papers used in this test suite.
        """
        result = client.get_pubtator_annotations("24667783")
        ann_types = {a["type"] for a in result["annotations"]}
        assert ann_types & {"Gene", "Disease"}, (
            f"Expected Gene or Disease annotations, got types: {ann_types}"
        )

    @pytest.mark.parametrize("pmid", TEST_PMIDS)
    def test_multiple_pmids(self, client, pmid):
        """All test PMIDs return a non-empty title."""
        result = client.get_pubtator_annotations(pmid)
        assert result["title"], f"PMID {pmid}: expected non-empty title"


# ------------------------------------------------------------------
# Variant-rich text
# ------------------------------------------------------------------


@pytest.mark.requires_network
class TestVariantRichText:
    def test_returns_expected_keys(self, client):
        result = client.get_variant_rich_text("24667783")
        assert result["pmid"] == "24667783"
        assert result["title"]
        assert isinstance(result["annotations"], list)
        assert isinstance(result["text"], str)
        assert len(result["text"]) > 0

    def test_annotations_are_mutations_only(self, client):
        result = client.get_variant_rich_text("24667783")
        for ann in result["annotations"]:
            assert "mut" in ann["type"].lower() or "variant" in ann["type"].lower(), (
                f"Expected only mutation annotations, got {ann['type']}"
            )
