"""Offline tests for the supplement-file value object."""

from gene_literature.supplements import SupplementFile


def test_to_dict():
    supplement = SupplementFile(
        url="https://example.com/mmc1.pdf", name="mmc1.pdf", source="test"
    )

    assert supplement.to_dict() == {
        "url": "https://example.com/mmc1.pdf",
        "name": "mmc1.pdf",
    }


def test_from_dict():
    payload = {"url": "https://example.com/supp.xlsx", "name": "supp.xlsx"}

    supplement = SupplementFile.from_dict(payload, source="scraper")

    assert supplement.url == payload["url"]
    assert supplement.name == payload["name"]
    assert supplement.source == "scraper"


def test_extension():
    supplement = SupplementFile(url="https://example.com/data.xlsx", name="data.xlsx")

    assert supplement.extension == "xlsx"


def test_normalized_url_removes_fragment_and_trailing_slash():
    with_fragment = SupplementFile(
        url="https://example.com/file.pdf#page=2", name="file.pdf"
    )
    with_trailing_slash = SupplementFile(
        url="https://example.com/file.pdf/", name="file.pdf"
    )

    assert with_fragment.normalized_url == "https://example.com/file.pdf"
    assert with_trailing_slash.normalized_url == "https://example.com/file.pdf"
