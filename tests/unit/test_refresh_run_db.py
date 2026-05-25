import json

from scripts.refresh_run_db import select_replay_candidates


def test_selects_stale_abstract_only_when_fulltext_exists(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "11111111"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\nKCNH2 c.1601G>A R534Q\n" * 40,
        encoding="utf-8",
    )
    (extraction_dir / f"KCNH2_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [],
                "extraction_metadata": {
                    "abstract_only": True,
                    "notes": "Abstract-only extraction; full text not available",
                },
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="KCNH2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=5,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert "stale_abstract_only" in candidates[0].reasons


def test_does_not_select_stale_abstract_only_when_only_fallback_exists(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "33333333"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n"
        "> **WARNING:** Full text could not be retrieved for PMID 33333333.\n"
        "> This document contains only the PubMed abstract and metadata.\n",
        encoding="utf-8",
    )
    (extraction_dir / f"KCNH2_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [],
                "extraction_metadata": {
                    "abstract_only": True,
                    "notes": "Abstract-only extraction; full text not available",
                },
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="KCNH2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=5,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
    )

    assert candidates == []


def test_selects_deterministic_vertical_table_lift(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "22222222"
    rows = "\n".join(f"RYR2\n{1258 + idx}c>t\nR{420 + idx}W\nNT" for idx in range(6))
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "Table S1. Subject Clinical and Genetic Characteristics\n" + rows,
        encoding="utf-8",
    )
    (extraction_dir / f"RYR2_PMID_{pmid}.json").write_text(
        json.dumps(
            {"variants": [], "extraction_metadata": {"total_variants_found": 0}}
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="RYR2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=5,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert candidates[0].deterministic_variants == 6
    assert any(r.startswith("deterministic_parser_lift") for r in candidates[0].reasons)
