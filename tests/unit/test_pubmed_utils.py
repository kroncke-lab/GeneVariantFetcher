import pytest
from urllib.error import HTTPError

pytest.importorskip("Bio")
from Bio.Entrez.Parser import ValidationError

from utils import pubmed_utils


class DummyHandle:
    def __init__(self, data=None, text=None):
        self.data = data
        self.text = text

    def read(self):
        return self.text if self.text is not None else ""

    def close(self):
        pass


def test_query_pubmed_with_entrez_retries_on_transient_error(monkeypatch):
    call_count = {"count": 0}

    def fake_esearch(**kwargs):
        call_count["count"] += 1
        if call_count["count"] < 3:
            raise HTTPError(
                url="https://example.com",
                code=500,
                msg="server error",
                hdrs=None,
                fp=None,
            )
        return DummyHandle({"IdList": ["12345"]})

    monkeypatch.setattr(pubmed_utils.Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(pubmed_utils.Entrez, "read", lambda handle: handle.data)

    pmids = pubmed_utils.query_pubmed_with_entrez("test")

    assert pmids == ["12345"]
    assert call_count["count"] == 3


def test_query_pubmed_with_entrez_raises_after_retries(monkeypatch):
    call_count = {"count": 0}

    def failing_esearch(**kwargs):
        call_count["count"] += 1
        raise HTTPError(
            url="https://example.com",
            code=503,
            msg="service unavailable",
            hdrs=None,
            fp=None,
        )

    monkeypatch.setattr(pubmed_utils.Entrez, "esearch", failing_esearch)
    monkeypatch.setattr(pubmed_utils.Entrez, "read", lambda handle: handle.data)

    with pytest.raises(HTTPError):
        pubmed_utils.query_pubmed_with_entrez("test")

    assert call_count["count"] == 3


def test_query_pubmed_with_entrez_handles_parse_errors(monkeypatch):
    call_count = {"count": 0}

    def fake_esearch(**kwargs):
        call_count["count"] += 1
        return DummyHandle({"IdList": ["12345"]})

    def bad_read(handle):
        raise ValidationError("Malformed data")

    monkeypatch.setattr(pubmed_utils.Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(pubmed_utils.Entrez, "read", bad_read)

    pmids = pubmed_utils.query_pubmed_with_entrez("test")

    assert pmids == []
    assert call_count["count"] == 1
