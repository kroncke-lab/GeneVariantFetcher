import sys
import types

import pytest

pydantic_stub = types.SimpleNamespace(
    BaseModel=type("BaseModel", (object,), {}),
    Field=lambda *args, **kwargs: None,
)
sys.modules.setdefault("pydantic", pydantic_stub)
sys.modules.setdefault("litellm", types.SimpleNamespace(completion=lambda *args, **kwargs: None))
sys.modules.setdefault(
    "tenacity",
    types.SimpleNamespace(
        retry=lambda *args, **kwargs: (lambda func: func),
        stop_after_attempt=lambda *args, **kwargs: None,
        wait_exponential=lambda *args, **kwargs: None,
        retry_if_exception_type=lambda *args, **kwargs: None,
    ),
)
sys.modules.setdefault(
    "requests",
    types.SimpleNamespace(
        RequestException=Exception,
        exceptions=types.SimpleNamespace(RequestException=Exception),
        get=lambda *args, **kwargs: None,
        Session=lambda *args, **kwargs: None,
    ),
)
sys.modules.setdefault(
    "bs4",
    types.SimpleNamespace(BeautifulSoup=type("BeautifulSoup", (object,), {})),
)
sys.modules.setdefault("Bio", types.SimpleNamespace(Entrez=type("Entrez", (), {})))
sys.modules.setdefault("dotenv", types.SimpleNamespace(load_dotenv=lambda *args, **kwargs: None))
sys.modules.setdefault(
    "pandas",
    types.SimpleNamespace(
        DataFrame=type("DataFrame", (object,), {}),
        concat=lambda *args, **kwargs: None,
        read_csv=lambda *args, **kwargs: None,
    ),
)

from pipeline.sourcing import PaperSourcer
from utils import pubmed_utils


class DummyHandle:
    def close(self):
        return None


@pytest.fixture
def mock_entrez(monkeypatch):
    recorded = {}

    def fake_esearch(*args, **kwargs):
        recorded["email"] = pubmed_utils.Entrez.email
        return DummyHandle()

    def fake_read(handle):
        return {"IdList": []}

    pubmed_utils.Entrez.email = None
    monkeypatch.setattr(pubmed_utils.Entrez, "esearch", fake_esearch, raising=False)
    monkeypatch.setattr(pubmed_utils.Entrez, "read", fake_read, raising=False)

    return recorded


def test_papersourcer_uses_provided_email(mock_entrez):
    sourcer = PaperSourcer(email="test@example.com")

    sourcer.fetch_papers(
        "GENE",
        use_pubmed=True,
        use_europepmc=False,
        use_pubmind=False,
    )

    assert mock_entrez["email"] == "test@example.com"
