"""Regression tests for stub DOI recovery.

Abstract-only / paywall stubs used to reach source recovery with no DOI, so
``fetch_paywalled.fetch_one`` short-circuited on "no DOI" and every publisher /
proxy / PMC route was skipped — even though the DOI is almost always available
(printed in the abstract we already saved, and resolvable from PubMed). These
tests lock in the three fixes that close that gap:

  * the audit fills the DOI from the on-disk abstract text, so the fetch queue
    carries a routable DOI;
  * the harvester's abstract-only fallback extracts/persists the DOI;
  * ``fetch_paywalled`` resolves a missing DOI from the PMID before giving up.
"""

from __future__ import annotations

import csv
from pathlib import Path

from harvesting import orchestrator as orch
from scripts import fetch_paywalled as fp
from scripts.recall_audit import source_acquisition_audit as audit


# --------------------------------------------------------------------------- #
# 1. Audit: recover the DOI from an abstract-only stub's saved text
# --------------------------------------------------------------------------- #

_ABSTRACT_STUB = (
    "# ABSTRACT-ONLY FALLBACK\n\n"
    "> **WARNING:** Full text could not be retrieved for PMID 11332568.\n"
    "> Reason: No PMCID and free-text/Tier 3.5 fallback failed: No PMCID\n\n"
    "## Title\nElectrocardiographic prediction of abnormal genotype.\n\n"
    "## Abstract\n"
    "J Cardiovasc Electrophysiol. 2001 Apr;12(4):455-61. "
    "doi: 10.1046/j.1540-8167.2001.00455.x.\n\n"
    "Congenital long QT syndrome cohort of 101 related family members.\n"
)


def test_doi_from_source_text_prefers_labelled_line(tmp_path: Path):
    stub = tmp_path / "s.md"
    stub.write_text(_ABSTRACT_STUB, encoding="utf-8")
    assert audit._doi_from_source_text(stub) == "10.1046/j.1540-8167.2001.00455.x"


def test_doi_from_source_text_missing_file_is_empty():
    assert audit._doi_from_source_text(None) == ""
    assert audit._doi_from_source_text(Path("/no/such/file.md")) == ""


def test_audit_fills_doi_for_abstract_stub(tmp_path: Path):
    """A DOI-less abstract stub must enter the fetch queue with a DOI + route."""
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    (run_dir / "extractions").mkdir(parents=True)
    harvest_dir.mkdir()
    pmid = "11332568"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        _ABSTRACT_STUB, encoding="utf-8"
    )

    rows, _summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    row = next(r for r in rows if r["pmid"] == pmid)

    assert row["action"] == "fetch"
    assert row["doi"] == "10.1046/j.1540-8167.2001.00455.x"
    # No longer the DOI-less placeholder route that strands the paper.
    assert row["route"] != "doi_lookup_then_fetch"

    # And it is written into fetch_input.csv with the DOI populated.
    fetch_input = tmp_path / "fetch_input.csv"
    audit.write_fetch_input(rows, fetch_input)
    with fetch_input.open(newline="", encoding="utf-8") as handle:
        fetch_rows = list(csv.DictReader(handle))
    assert any(
        r["PMID"] == pmid and r["DOI"] == "10.1046/j.1540-8167.2001.00455.x"
        for r in fetch_rows
    )


# --------------------------------------------------------------------------- #
# 2. Orchestrator: DOI extraction used by the abstract-only fallback
# --------------------------------------------------------------------------- #


def test_orchestrator_doi_from_text_labelled():
    text = "Circulation. 2010. doi: 10.1161/CIRCEP.109.882159. Long QT cohort."
    assert orch._doi_from_text(text) == "10.1161/CIRCEP.109.882159"


def test_orchestrator_doi_from_text_bare_and_url():
    assert orch._doi_from_text(
        "see https://doi.org/10.1016/j.hrthm.2009.10.032 now"
    ) == ("10.1016/j.hrthm.2009.10.032")
    assert orch._doi_from_text("no identifier here") == ""


def test_orchestrator_clean_doi_strips_wrapping():
    assert orch._clean_doi("https://dx.doi.org/10.1002/elps.201000022)") == (
        "10.1002/elps.201000022"
    )


# --------------------------------------------------------------------------- #
# 3. fetch_paywalled: resolve a missing DOI from the PMID before giving up
# --------------------------------------------------------------------------- #


def test_resolve_doi_gated_on_email(monkeypatch):
    fp._DOI_RESOLVE_CACHE.clear()
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    called = {"n": 0}

    def _boom(*_a, **_k):  # must not be called without an email
        called["n"] += 1
        return "10.9999/should-not-happen"

    monkeypatch.setattr(fp, "get_doi_from_pmid", _boom)
    assert fp.resolve_doi_for_pmid("123") is None
    assert called["n"] == 0


def test_resolve_doi_from_pmid_and_cache(monkeypatch):
    fp._DOI_RESOLVE_CACHE.clear()
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    calls = {"n": 0}

    def _lookup(pmid, email=None):
        calls["n"] += 1
        return "10.1093/europace/eup446"

    monkeypatch.setattr(fp, "get_doi_from_pmid", _lookup)
    assert fp.resolve_doi_for_pmid("20123697") == "10.1093/europace/eup446"
    # Cached — a second call does not hit the network again.
    assert fp.resolve_doi_for_pmid("20123697") == "10.1093/europace/eup446"
    assert calls["n"] == 1


def test_resolve_doi_degrades_on_error(monkeypatch):
    fp._DOI_RESOLVE_CACHE.clear()
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")

    def _raise(*_a, **_k):
        raise RuntimeError("network down")

    monkeypatch.setattr(fp, "get_doi_from_pmid", _raise)
    assert fp.resolve_doi_for_pmid("999") is None
