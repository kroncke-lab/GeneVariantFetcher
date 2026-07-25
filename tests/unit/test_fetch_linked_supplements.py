"""Tests for scripts/fetch_linked_supplements.py discovery and routing.

gvf-run's source-recovery step drives this over a *flat* run harvest directory
where every paper shares one parent and is distinguished only by a ``<PMID>_``
prefix. An unscoped glob there reports another paper's files as this one's, so
the PMID scoping is load-bearing, not cosmetic.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "fetch_linked_supplements", REPO / "scripts" / "fetch_linked_supplements.py"
)
fls = importlib.util.module_from_spec(_spec)
sys.modules["fetch_linked_supplements"] = fls
_spec.loader.exec_module(fls)


LINKED_MD = (
    "## SUPPLEMENT DESCRIPTIONS\n\n"
    "### Supplementary Table S1\n\n"
    "_link_: mmc1.xls\n\n"
    "### Supplementary Table S2\n\n"
    "_link_: mmc2.pdf\n"
)


def _corpus_paper(root: Path, gene: str, pmid: str, *, pmcid="PMC1", md=LINKED_MD):
    pdir = root / gene / pmid
    pdir.mkdir(parents=True)
    (pdir / f"{pmid}_FULL_CONTEXT.md").write_text(md)
    (pdir / f"{pmid}_artifacts.json").write_text('{"pmcid": "%s"}' % pmcid)
    return pdir


def _flat_paper(harvest: Path, pmid: str, *, pmcid="PMC1", md=LINKED_MD):
    harvest.mkdir(parents=True, exist_ok=True)
    (harvest / f"{pmid}_FULL_CONTEXT.md").write_text(md)
    (harvest / f"{pmid}_artifacts.json").write_text('{"pmcid": "%s"}' % pmcid)


# ---------------------------------------------------------------------------
# Corpus layout
# ---------------------------------------------------------------------------


def test_corpus_paper_with_links_and_empty_supplements_is_a_target(tmp_path):
    _corpus_paper(tmp_path, "KCNQ1", "111")
    targets = fls.discover_targets(tmp_path)
    assert [(t.gene, t.pmid, t.pmcid) for t in targets] == [("KCNQ1", "111", "PMC1")]
    assert len(targets[0].links) == 2


def test_paper_whose_files_are_all_present_is_not_a_target(tmp_path):
    pdir = _corpus_paper(tmp_path, "KCNQ1", "111")
    supp = pdir / "111_supplements"
    supp.mkdir()
    (supp / "mmc1.xls").write_bytes(b"x")
    (supp / "mmc2.pdf").write_bytes(b"x")
    assert fls.discover_targets(tmp_path) == []


def test_partially_fetched_paper_is_skipped_by_default_and_kept_with_the_flag(tmp_path):
    pdir = _corpus_paper(tmp_path, "KCNQ1", "111")
    supp = pdir / "111_supplements"
    supp.mkdir()
    (supp / "mmc1.xls").write_bytes(b"x")

    assert fls.discover_targets(tmp_path) == []
    partial = fls.discover_targets(tmp_path, include_partial=True)
    assert [t.pmid for t in partial] == ["111"]
    assert [link.href for link in partial[0].links] == ["mmc2.pdf"]


def test_gene_and_pmid_filters_narrow_discovery(tmp_path):
    _corpus_paper(tmp_path, "KCNQ1", "111")
    _corpus_paper(tmp_path, "SCN5A", "222")
    assert [t.pmid for t in fls.discover_targets(tmp_path, genes={"SCN5A"})] == ["222"]
    assert [t.pmid for t in fls.discover_targets(tmp_path, wanted_pmids={"111"})] == [
        "111"
    ]


def test_paper_without_file_shaped_links_is_ignored(tmp_path):
    _corpus_paper(tmp_path, "KCNQ1", "111", md="### Fig 1\n\n_link_: fig1.jpg\n")
    assert fls.discover_targets(tmp_path) == []


# ---------------------------------------------------------------------------
# Flat run harvest layout (what gvf-run passes)
# ---------------------------------------------------------------------------


def test_flat_harvest_dir_discovers_each_paper(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    _flat_paper(harvest, "111", pmcid="PMC1")
    _flat_paper(harvest, "222", pmcid="PMC2")
    targets = fls.discover_targets(tmp_path, harvest_dir=harvest, genes={"KCNQ1"})
    assert sorted(t.pmid for t in targets) == ["111", "222"]
    assert {t.gene for t in targets} == {"KCNQ1"}
    assert all(t.paper_dir == harvest for t in targets)


def test_flat_layout_does_not_credit_one_paper_with_anothers_supplements(tmp_path):
    """The bug an unscoped ``*_supplements`` glob would cause."""
    harvest = tmp_path / "pmc_fulltext"
    _flat_paper(harvest, "111")
    _flat_paper(harvest, "222")
    supp = harvest / "222_supplements"
    supp.mkdir()
    (supp / "mmc1.xls").write_bytes(b"x")
    (supp / "mmc2.pdf").write_bytes(b"x")

    targets = fls.discover_targets(tmp_path, harvest_dir=harvest)
    # 222 is complete; 111 must still be a target with both links outstanding.
    assert [t.pmid for t in targets] == ["111"]
    assert len(targets[0].links) == 2


def test_flat_markdown_scan_is_pmid_scoped(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    _flat_paper(harvest, "111", md="### T\n\n_link_: only_111.xls\n")
    _flat_paper(harvest, "222", md="### T\n\n_link_: only_222.xls\n")
    targets = {t.pmid: t for t in fls.discover_targets(tmp_path, harvest_dir=harvest)}
    assert [link.href for link in targets["111"].links] == ["only_111.xls"]
    assert [link.href for link in targets["222"].links] == ["only_222.xls"]


def test_missing_harvest_dir_is_not_an_error(tmp_path):
    assert fls.discover_targets(tmp_path, harvest_dir=tmp_path / "nope") == []


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------


def test_plan_stays_per_link_so_the_archive_can_genuinely_fall_back(tmp_path):
    """plan() must not substitute the archive job.

    It used to, which made ``report.fetchable`` *be* the archive job — so when
    Europe PMC declined, the "per-file fallback" retried the same archive URL.
    """
    _corpus_paper(tmp_path, "KCNQ1", "111", pmcid="PMC9")
    target = fls.discover_targets(tmp_path)[0]
    urls = [job.url for job in fls.plan(target).fetchable]
    assert urls == [
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9/bin/mmc1.xls",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9/bin/mmc2.pdf",
    ]
    assert not any("europepmc" in u for u in urls)


def test_href_filename_ignores_a_query_string():
    assert fls._href_filename("https://x.org/bin/a.pdf?dl=1") == "a.pdf"
    assert fls._href_filename("/articles/instance/9/bin/B.XLS") == "b.xls"


@pytest.mark.parametrize("body", [b"PK\x03\x04rest", b"<?xml ...errorBean"])
def test_is_zip_accepts_only_a_real_archive(tmp_path, body):
    path = tmp_path / "x.zip"
    path.write_bytes(body)
    assert fls._is_zip(path) == body.startswith(b"PK\x03\x04")


def test_source_override_rows_use_the_recorded_path_not_a_rebuilt_one(tmp_path):
    """In a flat harvest dir there is no <gene>/<pmid>/ path to rebuild from."""
    out = tmp_path / "ov.csv"
    fls._write_source_overrides(
        out,
        [
            {
                "pmid": "111",
                "folded": True,
                "full_context": "/run/pmc_fulltext/111_FULL_CONTEXT.md",
            },
            {"pmid": "222", "folded": False, "full_context": "/run/x.md"},
        ],
    )
    lines = out.read_text().strip().split("\n")
    assert lines[0] == "pmid,source_path"
    assert lines[1] == "111,/run/pmc_fulltext/111_FULL_CONTEXT.md"
    assert len(lines) == 2  # unfolded papers are not offered to refresh_run_db


# ---------------------------------------------------------------------------
# Per-paper failure breaker
# ---------------------------------------------------------------------------


class _FakeHarvester:
    """Records download attempts; every one fails, as a gated article does."""

    def __init__(self, results=None):
        self.attempts = []
        self.results = results or {}

    def download_supplement(self, url, file_path, pmid, filename, base, original):
        self.attempts.append(filename)
        return self.results.get(filename, False)


def _call(cb, name, tmp_path):
    return cb("http://x/" + name, tmp_path / name, "1", name, {})


def test_per_file_route_stops_after_the_budget_of_consecutive_failures(tmp_path):
    h = _FakeHarvester()
    cb = fls._download_callback(h, set(), failure_budget=2)
    for name in ("a.pdf", "b.pdf", "c.pdf", "d.pdf"):
        assert _call(cb, name, tmp_path) is False
    # Only the first two hit the network; the rest are short-circuited.
    assert h.attempts == ["a.pdf", "b.pdf"]


def test_a_success_resets_the_failure_budget(tmp_path):
    h = _FakeHarvester(results={"b.pdf": True})
    cb = fls._download_callback(h, set(), failure_budget=2)
    for name in ("a.pdf", "b.pdf", "c.pdf", "d.pdf"):
        _call(cb, name, tmp_path)
    # b succeeded, so c and d still get their chance.
    assert h.attempts == ["a.pdf", "b.pdf", "c.pdf", "d.pdf"]


def test_the_archive_route_has_no_breaker(tmp_path):
    """One request per paper either way — a budget would only add a failure mode."""
    h = _FakeHarvester()
    cb = fls._download_callback(h, {"z.zip"}, failure_budget=None)
    for name in ("a.zip", "b.zip", "c.zip"):
        _call(cb, name, tmp_path)
    assert h.attempts == ["a.zip", "b.zip", "c.zip"]


def test_a_non_archive_payload_named_zip_is_rejected_and_deleted(tmp_path):
    """Europe PMC returns a 200 XML errorBean for closed-access articles."""
    target = tmp_path / "PMC1_supplements.zip"

    class _Writer(_FakeHarvester):
        def download_supplement(self, url, file_path, pmid, filename, base, original):
            file_path.write_bytes(b'<?xml version="1.0"?><errorBean/>')
            return True

    cb = fls._download_callback(_Writer(), {"PMC1_supplements.zip"})
    assert cb("http://x", target, "1", "PMC1_supplements.zip", {}) is False
    assert not target.exists()


def test_a_real_archive_is_kept(tmp_path):
    target = tmp_path / "PMC1_supplements.zip"

    class _Writer(_FakeHarvester):
        def download_supplement(self, url, file_path, pmid, filename, base, original):
            file_path.write_bytes(b"PK\x03\x04and the rest")
            return True

    cb = fls._download_callback(_Writer(), {"PMC1_supplements.zip"})
    assert cb("http://x", target, "1", "PMC1_supplements.zip", {}) is True
    assert target.exists()
