from harvesting.unpaywall_api import UnpaywallClient


class _DummyResponse:
    status_code = 200
    reason = "OK"
    headers = {"content-type": "application/json"}
    text = '{"doi":"10.1000/test","is_oa":false,"oa_status":"closed"}'

    def json(self):
        return {"doi": "10.1000/test", "is_oa": False, "oa_status": "closed"}


class _DummySession:
    def __init__(self):
        self.calls = []

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append(
            {"url": url, "params": params, "timeout": timeout, "headers": headers}
        )
        return _DummyResponse()


def test_find_open_access_uses_api_headers_with_browser_session():
    session = _DummySession()
    client = UnpaywallClient(email="researcher@example.org", session=session)

    result, error = client.find_open_access("10.1000/test")

    assert error is None
    assert result["doi"] == "10.1000/test"
    headers = session.calls[0]["headers"]
    assert headers["Accept"] == "application/json"
    assert headers["Accept-Encoding"] == "identity"
    assert "researcher@example.org" in headers["User-Agent"]
