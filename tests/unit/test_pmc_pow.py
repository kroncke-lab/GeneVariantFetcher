from harvesting.pmc_pow import (
    attach_pmc_pow_cookie,
    is_pmc_pow_challenge,
    solve_pmc_pow_cookie,
)


def test_solves_pmc_pow_cookie_for_download_gate():
    html = """
    <html><body>Preparing to download ...</body>
    <script>
    const POW_CHALLENGE = "unit-test-challenge"
    const POW_DIFFICULTY = "2"
    const POW_COOKIE_NAME = "cloudpmc-viewer-pow"
    </script></html>
    """

    cookie = solve_pmc_pow_cookie(html)

    assert cookie is not None
    assert cookie.startswith("unit-test-challenge,")


def test_declines_pmc_pow_above_difficulty_ceiling():
    html = """
    const POW_CHALLENGE = "unit-test-challenge"
    const POW_DIFFICULTY = "9"
    """

    assert solve_pmc_pow_cookie(html, max_difficulty=5) is None


def test_detects_pmc_pow_challenge_page():
    content = b"""
    <html><body>Preparing to download ...</body>
    <script>
    const POW_CHALLENGE = "abc"
    const POW_COOKIE_NAME = "cloudpmc-viewer-pow"
    </script></html>
    """

    assert is_pmc_pow_challenge(content)


def test_attach_pmc_pow_cookie_sets_session_cookie():
    class _Cookies:
        def __init__(self):
            self.calls = []

        def set(self, name, value, domain=None, path=None):
            self.calls.append((name, value, domain, path))

    class _Session:
        def __init__(self):
            self.cookies = _Cookies()

    html = """
    const POW_CHALLENGE = "abc"
    const POW_DIFFICULTY = "1"
    """
    session = _Session()

    assert attach_pmc_pow_cookie(
        session, html=html, url="https://pmc.ncbi.nlm.nih.gov/articles/x/bin/y.pdf"
    )

    name, value, domain, path = session.cookies.calls[0]
    assert name == "cloudpmc-viewer-pow"
    assert value.startswith("abc,")
    assert domain == "pmc.ncbi.nlm.nih.gov"
    assert path == "/"
