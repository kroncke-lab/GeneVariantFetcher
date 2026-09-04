"""Integrity checks for the Gate-2-sized mixed-gold continuation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from benchmarks.evaluation_tiers.build_mixed_tranches import sha256_file


REPO = Path(__file__).parents[2]
OLD = REPO / "benchmarks" / "evaluation_tiers" / "mixed_gold"
SUITE = REPO / "benchmarks" / "evaluation_tiers" / "mixed_gold_continuation_120"


def manifest_rows(path: Path) -> list[tuple[str, str]]:
    return [
        tuple(line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]


def test_continuation_is_gate2_sized_and_excludes_consumed_articles():
    registry = json.loads((SUITE / "registry.json").read_text())
    opened = set()
    for name in ("tranche_01.tsv", "tranche_02.tsv"):
        opened.update(pmid for _, pmid in manifest_rows(OLD / name))

    assigned: list[tuple[str, str]] = []
    assigned_pmids: dict[str, str] = {}
    for tier in registry["tiers"]:
        rows = manifest_rows(SUITE / tier["manifest"])
        assert len(rows) == tier["attempt_count"]
        assert tier["attempt_count"] in {120, 121}
        assert sha256_file(SUITE / tier["manifest"]) == tier["sha256"]
        for gene, pmid in rows:
            assert pmid not in opened
            assert assigned_pmids.setdefault(pmid, tier["id"]) == tier["id"]
            assigned.append((gene, pmid))

    with (SUITE / "inventory.tsv").open(newline="") as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    included = {
        (row["gene"], row["pmid"]) for row in inventory if row["status"] == "included"
    }
    assert set(assigned) == included
    assert len(assigned) == 1324
    assert registry["inventory"]["previously_consumed_attempts"] == 98
    assert registry["inventory"]["previously_consumed_pmids"] == 90
    assert registry["evaluation_design"]["consume_order"][0] == (
        "mixed_gold_cont120_01"
    )


def test_old_unopened_tiers_have_append_only_abandonment_record():
    events = [
        json.loads(line)
        for line in (OLD / "abandonment_log.jsonl").read_text().splitlines()
        if line.strip()
    ]
    event = events[-1]
    assert event["preserved_consumed_tiers"] == ["mixed_gold_01", "mixed_gold_02"]
    assert event["tier_ids"] == [f"mixed_gold_{number:02d}" for number in range(3, 30)]
    assert event["replacement_registry_sha256"] == sha256_file(SUITE / "registry.json")
