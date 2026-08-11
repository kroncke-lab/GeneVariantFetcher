"""Integrity checks for the three canonical paper-evaluation tiers."""

import hashlib
import json
from collections import Counter
from pathlib import Path

from utils.pmid_utils import is_valid_pmid


ROOT = Path(__file__).resolve().parents[2]
TIERS = ROOT / "benchmarks" / "evaluation_tiers"
REVIEW = ROOT / "benchmarks" / "curated_extraction_eval"
PAPER_EVAL = ROOT / "benchmarks" / "codex_paper_eval"


def _rows(path):
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        gene, pmid = line.split("\t")
        rows.append((gene, pmid))
    return rows


def _pmids(path):
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]


def test_registry_names_exactly_three_ordered_active_tiers():
    registry = json.loads((TIERS / "registry.json").read_text())
    tiers = registry["tiers"]

    assert registry["count_unit"] == "gene_paper_attempt"
    assert [tier["id"] for tier in tiers] == [
        "gold_50",
        "cardiac_120",
        "reviewer_546",
    ]
    assert [tier["order"] for tier in tiers] == [1, 2, 3]


def test_tier_manifests_match_registry_counts_and_valid_pmids():
    registry = json.loads((TIERS / "registry.json").read_text())

    for tier in registry["tiers"]:
        manifest = TIERS / tier["manifest"]
        rows = _rows(manifest)
        assert hashlib.sha256(manifest.read_bytes()).hexdigest() == tier["sha256"]
        assert len(rows) == tier["attempt_count"]
        assert len(rows) == len(set(rows))
        assert len({pmid for _, pmid in rows}) == tier["unique_pmid_count"]
        assert Counter(gene for gene, _ in rows) == Counter(tier["gene_attempt_counts"])
        assert all(is_valid_pmid(pmid) for _, pmid in rows)


def test_tier1_is_the_collaborator_grounded_48_plus_2_gate():
    active_50 = _rows(TIERS / "tier1_gold_50.tsv")
    source = _rows(PAPER_EVAL / "highcarrier48_plus_brca2_collaborator2_20260811.tsv")

    assert active_50 == source
    assert [pmid for gene, pmid in active_50 if gene == "BRCA2"] == [
        "26833046",
        "26848529",
    ]


def test_tier2_is_top_30_of_each_ranked_cardiac_review_queue():
    actual = _rows(TIERS / "tier2_cardiac_120.tsv")
    expected = []
    for gene in ("KCNH2", "KCNQ1", "RYR2", "SCN5A"):
        expected.extend(
            (gene, pmid)
            for pmid in _pmids(REVIEW / "review_pmids_50" / f"{gene}.txt")[:30]
        )

    assert actual == expected


def test_tier3_tracks_all_eleven_reviewer_workspaces():
    actual = set(_rows(TIERS / "tier3_reviewer_546.tsv"))
    expected = set()

    for gene in ("APOE", "BRCA1", "KCNH2", "KCNQ1", "MYBPC3", "RYR2", "SCN5A"):
        expected.update(
            (gene, pmid) for pmid in _pmids(REVIEW / "review_pmids_50" / f"{gene}.txt")
        )
    expected.update(
        ("BRCA2", pmid)
        for pmid in _pmids(
            REVIEW / "review_pmids_20260811_brca2_provenance" / "BRCA2.txt"
        )
    )

    # The live SCN5A snapshot differs by one paper from the frozen July list.
    expected.remove(("SCN5A", "15687307"))
    expected.add(("SCN5A", "16414944"))
    for gene in ("BMPR2", "LMNA", "TTN"):
        expected.update(
            (gene, pmid)
            for pmid in _pmids(TIERS / "reviewer_pmids_50_20260811" / f"{gene}.txt")
        )

    assert actual == expected


def test_new_reviewer_order_manifests_are_pinned_and_in_full_tier():
    registry = json.loads((TIERS / "registry.json").read_text())
    tier3 = registry["tiers"][2]
    full_tier = set(_rows(TIERS / tier3["manifest"]))
    order_manifests = tier3["review_order_manifests"]

    assert set(order_manifests) == {"BMPR2", "LMNA", "TTN"}
    for gene, metadata in order_manifests.items():
        manifest = TIERS / metadata["manifest"]
        pmids = _pmids(manifest)
        assert hashlib.sha256(manifest.read_bytes()).hexdigest() == metadata["sha256"]
        assert len(pmids) == metadata["attempt_count"] == 50
        assert len(pmids) == len(set(pmids))
        assert all(is_valid_pmid(pmid) for pmid in pmids)
        assert {(gene, pmid) for pmid in pmids} <= full_tier


def test_cardiac_expansion_is_a_subset_of_the_full_reviewer_backlog():
    tier2 = set(_rows(TIERS / "tier2_cardiac_120.tsv"))
    tier3 = set(_rows(TIERS / "tier3_reviewer_546.tsv"))

    assert tier2 < tier3
