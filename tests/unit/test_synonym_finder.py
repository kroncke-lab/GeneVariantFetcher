"""Gene resolution in SynonymFinder — symbol-collision disambiguation.

Regression: "TTN" resolved to TTR/transthyretin (7276) instead of titin (7273)
because a bare ``[Gene Name]`` esearch ranks the more-studied colliding gene
first. Resolution now tries ``[Preferred Symbol]`` first.
"""

from gene_literature.synonym_finder import SynonymFinder


class _FakeResp:
    def __init__(self, idlist):
        self._idlist = idlist

    def json(self):
        return {"esearchresult": {"idlist": self._idlist}}


def _finder_with(responses):
    """Build a finder whose _request returns idlists keyed by esearch field.

    ``responses`` maps an esearch field name (e.g. "Preferred Symbol") to the
    idlist returned when that ``[field]`` appears in the query term.
    """
    sf = SynonymFinder(email="t@example.org")
    sf._calls = []

    def fake_request(url, params):
        term = params["term"]
        sf._calls.append(term)
        for field, ids in responses.items():
            if f"[{field}]" in term:
                return _FakeResp(ids)
        return _FakeResp([])

    sf._request = fake_request  # type: ignore[method-assign]
    return sf


def test_search_gene_prefers_preferred_symbol():
    # [Preferred Symbol] resolves to titin (7273); [Gene Name] would rank
    # TTR/transthyretin (7276) ahead of titin.
    sf = _finder_with({"Preferred Symbol": ["7273"], "Gene Name": ["7276", "7273"]})
    assert sf._search_gene("TTN") == 7273
    # Preferred Symbol hit short-circuits — the Gene Name fallback is not queried.
    assert any("[Preferred Symbol]" in t for t in sf._calls)
    assert all("[Gene Name]" not in t for t in sf._calls)


def test_search_gene_falls_back_to_gene_name():
    # A full gene name won't match [Preferred Symbol]; fall back to [Gene Name].
    sf = _finder_with({"Preferred Symbol": [], "Gene Name": ["4000"]})
    assert sf._search_gene("lamin A/C") == 4000
    assert any("[Preferred Symbol]" in t for t in sf._calls)
    assert any("[Gene Name]" in t for t in sf._calls)
    # Multi-word terms must be quoted so NCBI treats them as one phrase.
    assert all('"lamin A/C"' in t for t in sf._calls)


def test_search_gene_returns_none_when_no_hit():
    sf = _finder_with({"Preferred Symbol": [], "Gene Name": []})
    assert sf._search_gene("NOTAGENE123") is None
