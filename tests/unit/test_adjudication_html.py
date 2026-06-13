"""Tests for the self-contained HTML adjudication review packet builder."""

import json
import re

from scripts.build_adjudication_html import build_data, render_index


def _record():
    return {
        "pmid": "111",
        "title": "BRCA2 carriers",
        "url": "https://pubmed.ncbi.nlm.nih.gov/111/",
        "fulltext_file": "papers/111.md",
        "abstract": "germline carriers were screened in families",
        "variants": [
            {
                "cdna_notation": "c.1A>T",
                "protein_notation": None,
                "clinical_significance": "pathogenic",
                "penetrance_data": {
                    "total_carriers_observed": 5,
                    "affected_count": 3,
                    "unaffected_count": 2,
                },
                "patients": {"count": None},
                "count_provenance": {
                    "carriers_count_type": "explicit_total",
                    "carriers_column_label": "Carriers",
                },
                "source_location": "Table 2",
                "key_quotes": ["5 carriers, 3 affected"],
                "additional_notes": "from table 2",
            }
        ],
    }


def test_build_data_groups_variants_per_paper_with_provenance():
    (p,) = build_data([_record()])
    assert p["pmid"] == "111"
    assert p["fulltext"] == "papers/111.html"
    assert p["pubmed"].endswith("/111/")
    (v,) = p["variants"]
    assert v["variant"] == "c.1A>T"
    assert v["carriers"] == 5 and v["affected"] == 3 and v["unaffected"] == 2
    assert v["gs"] == "germline"
    assert "explicit_total" in v["prov"]
    assert v["loc"] == "Table 2"
    assert "5 carriers" in v["quote"]


def test_build_data_paper_with_no_variants_is_empty_list():
    rec = {
        "pmid": "222",
        "title": "",
        "url": "u",
        "fulltext_file": "papers/222.html",
        "abstract": "",
        "variants": [],
    }
    (p,) = build_data([rec])
    assert p["variants"] == []


def test_render_index_embeds_parseable_json_and_ui():
    papers = build_data([_record()])
    h = render_index(papers, "BRCA2")
    # placeholders fully substituted
    assert "__DATA__" not in h and "__GENE__" not in h
    # the embedded data block must be valid JSON (a malformed blob breaks the page)
    m = re.search(r"const DATA = (\{.*?\});\nconst KEY", h, re.S)
    assert m, "DATA block not found"
    data = json.loads(m.group(1))
    assert data["gene"] == "BRCA2"
    assert len(data["papers"]) == 1
    # the reviewer UI affordances are present
    assert "Export my answers" in h
    assert "Add a variant the AI missed" in h
