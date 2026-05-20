import sys
from types import SimpleNamespace

import pytest

from harvesting.paywall_context_enrichment import EnrichmentResult
from scripts import fetch_paywalled


def test_write_outputs_preserves_supplement_only_recovery(tmp_path, monkeypatch):
    result = SimpleNamespace(
        main_markdown=None,
        main_html="<html></html>",
        supp_files=[{"url": "https://example.org/supp.doc", "name": "supp.doc"}],
        figure_paths=[],
        notes=[],
        error=None,
        publisher="generic",
        final_url="https://example.org/article",
    )

    def fake_enrich(**kwargs):
        assert "supplement-only recovery" in kwargs["body_markdown"]
        return EnrichmentResult(
            unified_markdown=kwargs["body_markdown"] + "VARIANT_TABLE",
            supplement_results=[SimpleNamespace(downloaded=True)],
        )

    monkeypatch.setattr(fetch_paywalled, "enrich_paywall_full_context", fake_enrich)

    canonical, per_pmid, body, enrichment = fetch_paywalled.write_outputs(
        "12345",
        result,
        tmp_path,
    )

    assert body is None
    assert enrichment.supplement_count == 1
    assert "VARIANT_TABLE" in canonical.read_text(encoding="utf-8")
    assert per_pmid.read_text(encoding="utf-8") == canonical.read_text(encoding="utf-8")


def test_ncbi_email_requires_env(monkeypatch):
    monkeypatch.delenv("ENTREZ_EMAIL", raising=False)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)

    with pytest.raises(RuntimeError, match="ENTREZ_EMAIL"):
        fetch_paywalled._ncbi_email()


def test_main_requires_explicit_pmids(monkeypatch, tmp_path):
    monkeypatch.setattr(
        sys,
        "argv",
        ["fetch_paywalled", "--output", str(tmp_path), "--no-cookies"],
    )

    with pytest.raises(SystemExit) as exc:
        fetch_paywalled.main()

    assert exc.value.code == 2
