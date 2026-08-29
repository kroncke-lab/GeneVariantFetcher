"""Integrity checks for the four canonical paper-evaluation tiers."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from benchmarks.codex_paper_eval.run_eval import (
    DEFAULT_GOLD,
    gold_count_eligible_pmids,
)
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


def test_registry_names_exactly_four_ordered_active_tiers():
    registry = json.loads((TIERS / "registry.json").read_text())
    tiers = registry["tiers"]

    assert registry["count_unit"] == "gene_paper_attempt"
    assert [tier["id"] for tier in tiers] == [
        "gold_50",
        "gold_120",
        "reviewer_545",
        "gold_120b",
    ]
    assert [tier["order"] for tier in tiers] == [1, 2, 3, 4]


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


def test_tier2_is_seeded_30_per_gene_manual_gold_expansion():
    registry = json.loads((TIERS / "registry.json").read_text())
    actual = _rows(TIERS / "tier2_gold_120.tsv")

    assert registry["tiers"][1]["selection_seed"] == 2026081301
    # KCNH2 is 28 after quarantining two non-genetics PMIDs: 10086972 on
    # 2026-08-15 and 14642689 on 2026-08-21. The other genes stay at 30.
    assert Counter(gene for gene, _ in actual) == Counter(
        {"KCNH2": 28, "KCNQ1": 30, "RYR2": 30, "SCN5A": 30}
    )
    assert ("KCNH2", "10086972") not in actual
    assert ("KCNH2", "14642689") not in actual
    for gene in ("KCNH2", "KCNQ1", "RYR2", "SCN5A"):
        eligible = gold_count_eligible_pmids(DEFAULT_GOLD, gene)
        assert {pmid for row_gene, pmid in actual if row_gene == gene} <= eligible


def test_tier3_tracks_all_eleven_reviewer_workspaces():
    actual = set(_rows(TIERS / "tier3_reviewer_545.tsv"))
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


def test_gold_expansion_is_separate_from_the_full_reviewer_backlog():
    tier2 = set(_rows(TIERS / "tier2_gold_120.tsv"))
    tier3 = set(_rows(TIERS / "tier3_reviewer_545.tsv"))

    assert ("BRCA2", "19944633") not in tier3

    assert not tier2 <= tier3


def test_tier4_is_a_seeded_replication_tranche_of_the_same_gold_rule():
    registry = json.loads((TIERS / "registry.json").read_text())
    tier = next(t for t in registry["tiers"] if t["id"] == "gold_120b")
    actual = _rows(TIERS / tier["manifest"])

    assert tier["selection_seed"] == 2026082501
    assert tier["role"] == "scored_gold_replication"
    assert Counter(gene for gene, _ in actual) == Counter(
        {"BRCA2": 5, "KCNH2": 30, "KCNQ1": 30, "RYR2": 30, "SCN5A": 30}
    )
    # Same eligibility helper as tier 2, so the two tranches stay comparable.
    for gene in ("KCNH2", "KCNQ1", "RYR2", "SCN5A", "BRCA2"):
        eligible = gold_count_eligible_pmids(DEFAULT_GOLD, gene)
        assert {pmid for row_gene, pmid in actual if row_gene == gene} <= eligible
    for quarantined in (
        ("KCNH2", "10086972"),
        ("KCNH2", "14642689"),
        ("BRCA2", "19944633"),
    ):
        assert quarantined not in actual


def test_tier4_shares_no_article_with_any_prior_scored_or_staged_surface():
    """Replication only means something on text nothing has been tuned against.

    Attempt-level disjointness is not enough: a multi-gene paper already scored
    under KCNQ1 has had its tables and supplements optimized against, and
    BRCA1/BRCA2 output that differs only in the gene column is a recorded
    failure here. So the check is at PMID level, across every prior surface.
    """
    tier4 = {pmid for _, pmid in _rows(TIERS / "tier4_gold_120b.tsv")}

    prior = set()
    for manifest in (
        "tier1_gold_50.tsv",
        "tier2_gold_120.tsv",
        "tier3_reviewer_545.tsv",
    ):
        prior.update(pmid for _, pmid in _rows(TIERS / manifest))
    gold150 = REVIEW / "gold150_preregistered_20260824"
    for template in gold150.glob("*/*/curation_template.csv"):
        with template.open(newline="", encoding="utf-8-sig") as handle:
            prior.update(
                (row.get("pmid") or "").strip()
                for row in csv.DictReader(handle)
                if (row.get("pmid") or "").strip()
            )

    assert tier4 & prior == set()


def test_tier4_answer_key_is_frozen_and_covers_every_drawn_paper():
    provenance = json.loads((TIERS / "tier4_gold_120b_selection.json").read_text())
    rows = _rows(TIERS / "tier4_gold_120b.tsv")

    assert provenance["gold_value_blinded"] is True
    assert provenance["exclusion_level"] == "pmid"
    assert {(p["gene"], p["pmid"]) for p in provenance["papers"]} == set(rows)
    # Every drawn paper is bound to the exact source bytes it was drawn on.
    assert all(len(p["source_sha256"]) == 64 for p in provenance["papers"])

    key_root = TIERS / "gold_120b_answer_key"
    for entry in provenance["answer_key"]:
        key = key_root / entry["file"]
        assert hashlib.sha256(key.read_bytes()).hexdigest() == entry["sha256"]
        with key.open(newline="", encoding="utf-8-sig") as handle:
            key_pmids = {(r.get("pmid") or "").strip() for r in csv.DictReader(handle)}
        assert key_pmids == {pmid for gene, pmid in rows if gene == entry["gene"]}


def test_corrected_brca2_50_extends_the_safe_45_without_restoring_exclusions():
    active_45 = _pmids(REVIEW / "review_pmids_20260811_brca2_provenance" / "BRCA2.txt")
    corrected_50 = _pmids(REVIEW / "review_pmids_50_20260821" / "BRCA2.txt")

    assert corrected_50[:45] == active_45
    assert corrected_50[45:] == [
        "26183948",
        "25923920",
        "20380699",
        "25884701",
        "26843898",
    ]
    assert len(corrected_50) == len(set(corrected_50)) == 50
    assert all(is_valid_pmid(pmid) for pmid in corrected_50)

    excluded = {"15365993", "18489799", "19944633", "22655046", "25802882"}
    assert excluded.isdisjoint(corrected_50)
