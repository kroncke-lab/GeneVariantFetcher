import csv
import json
from pathlib import Path

from scripts.recall_audit import source_acquisition_audit as audit


def _write_article(path: Path, gene: str = "KCNH2") -> None:
    body = "# Article\n\nAbstract\nIntroduction\nMethods\nResults\n\n" + (
        f"{gene} p.Arg1His cohort table with carriers and probands.\n" * 220
    )
    path.write_text(body, encoding="utf-8")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_infer_publisher_covers_legacy_elsevier_prefixes():
    for doi in (
        "10.1016/example",
        "10.1053/example",
        "10.1067/example",
        "10.1054/example",
        "10.1006/example",
    ):
        assert audit.infer_publisher(doi) == "elsevier"


def test_audit_uses_full_context_when_data_zones_is_stub(tmp_path: Path):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    harvest_dir.mkdir(parents=True)
    extraction_dir.mkdir()

    pmid = "111"
    data_zones = harvest_dir / f"{pmid}_DATA_ZONES.md"
    data_zones.write_text(
        "# DATA ZONES\n\nNo high-value data zones identified.\n",
        encoding="utf-8",
    )
    full_context = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    _write_article(full_context)

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)

    assert len(rows) == 1
    row = rows[0]
    assert row["action"] == "refresh_replay"
    assert row["route"] == "missing_extraction_for_usable_source"
    assert row["source_path"] == str(full_context)
    assert row["available_context_path"] == str(full_context)
    assert row["source_status"] == "recovered_pmc"
    assert summary["pmid_coverage"]["usable_fulltext_current"]["pmids"] == 1
    assert summary["pmid_coverage"]["selected_for_source_refresh"]["pmids"] == 1


def test_audit_fetch_queue_reports_selected_pmids_without_gold(tmp_path: Path):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    (run_dir / "extractions").mkdir(parents=True)
    harvest_dir.mkdir()
    paywalled = harvest_dir / "paywalled_missing.csv"
    paywalled.write_text(
        "PMID,DOI,Reason\n222,10.1016/j.hrthm.2020.01.001,paywall\n",
        encoding="utf-8",
    )

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    fetch_input = tmp_path / "fetch_input.csv"
    audit.write_fetch_input(rows, fetch_input)

    assert rows == [
        {
            "gene": "KCNH2",
            "pmid": "222",
            "action": "fetch",
            "route": "fetch_elsevier_insttoken",
            "priority_score": 100,
            "ev_score": 0.0,
            "ev_est_carriers": 0.0,
            "ev_est_variants": 0.0,
            "ev_p_relevant": 0.0,
            "ev_population_flag": 0,
            "ev_note": "",
            "source_status": "not_attempted",
            "source_path": "",
            "source_bytes": 0,
            "available_context_path": "",
            "available_context_bytes": 0,
            "extraction_path": "",
            "extraction_variants": 0,
            "extraction_abstract_only": False,
            "source_unbound": False,
            "source_sha_mismatch": False,
            "missing_variant_supplement": False,
            "zero_variant_qc": False,
            "single_carrier_qc": False,
            "doi": "10.1016/j.hrthm.2020.01.001",
            "publisher_hint": "elsevier",
            "notes": "no usable run-local full text",
        }
    ]
    assert summary["pmid_coverage"]["selected_for_fetch"] == {
        "pmids": 1,
        "total_pmids": 1,
        "coverage": 1.0,
    }
    assert _read_csv(fetch_input) == [
        {
            "PMID": "222",
            "DOI": "10.1016/j.hrthm.2020.01.001",
            "route": "fetch_elsevier_insttoken",
        }
    ]


def test_fetch_queue_ordered_by_acquisition_ev(tmp_path: Path):
    """The un-downloaded fetch tail is ranked by predicted per-paper yield.

    A high-yield cohort abstract must be fetched before a single-case abstract,
    and every fetch row must carry a populated ev_score column.
    """
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    abstract_dir = run_dir / "abstract_json"
    (run_dir / "extractions").mkdir(parents=True)
    harvest_dir.mkdir()
    abstract_dir.mkdir()

    (abstract_dir / "111.json").write_text(
        json.dumps(
            {
                "metadata": {"pmid": "111", "title": "KCNH2 variants in a LQT2 cohort"},
                "abstract": (
                    "We genotyped 200 probands carrying p.Arg190Gln, c.1000A>G, "
                    "p.Gly628Ser and p.Ala561Val; 150 affected mutation-positive "
                    "carriers across 40 families were phenotyped."
                ),
            }
        ),
        encoding="utf-8",
    )
    (abstract_dir / "222.json").write_text(
        json.dumps(
            {
                "metadata": {"pmid": "222", "title": "A case report"},
                "abstract": "We report a single patient with a novel finding.",
            }
        ),
        encoding="utf-8",
    )
    (harvest_dir / "paywalled_missing.csv").write_text(
        "PMID,DOI,Reason\n"
        "111,10.1016/j.hrthm.2020.01.001,paywall\n"
        "222,10.1016/j.hrthm.2020.02.002,paywall\n",
        encoding="utf-8",
    )

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    by_pmid = {row["pmid"]: row for row in rows}
    assert by_pmid["111"]["ev_score"] > by_pmid["222"]["ev_score"] > 0
    assert by_pmid["111"]["ev_est_carriers"] == 200.0

    fetch_input = tmp_path / "fetch_input.csv"
    audit.write_fetch_input(rows, fetch_input)
    assert [r["PMID"] for r in _read_csv(fetch_input)] == ["111", "222"]

    top = summary["top_fetch_targets_by_ev"]
    assert top and top[0]["pmid"] == "111"
    assert summary["acquisition_ev_available"] is True


def test_audit_routes_known_blocked_publishers_out_of_fetch_queue(tmp_path: Path):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    (run_dir / "extractions").mkdir(parents=True)
    harvest_dir.mkdir()
    (harvest_dir / "paywalled_missing.csv").write_text(
        "PMID,DOI,Reason\n444,10.1089/gtmb.2019.001,paywall\n",
        encoding="utf-8",
    )

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    fetch_input = tmp_path / "fetch_input.csv"
    audit.write_fetch_input(rows, fetch_input)

    assert rows[0]["action"] == "manual_or_blocked"
    assert rows[0]["route"] == "blocked_liebert_sage_cloudflare"
    assert summary["pmid_coverage"]["selected_for_fetch"]["pmids"] == 0
    assert summary["pmid_coverage"]["selected_for_manual_or_blocked"]["pmids"] == 1
    assert _read_csv(fetch_input) == []


def test_audit_emits_source_override_for_zero_variant_fulltext(tmp_path: Path):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    harvest_dir.mkdir(parents=True)
    extraction_dir.mkdir()

    pmid = "333"
    full_context = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    _write_article(full_context)
    extraction = extraction_dir / f"KCNH2_PMID_{pmid}.json"
    extraction.write_text(
        json.dumps(
            {
                "variants": [],
                "extraction_metadata": {
                    "pmid": pmid,
                    "source_file": str(full_context),
                },
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "source_completeness.json").write_text(
        json.dumps({"zero_variant_pmids": [pmid]}),
        encoding="utf-8",
    )

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    override = tmp_path / "source_override.csv"
    audit.write_source_override(rows, override)

    assert len(rows) == 1
    row = rows[0]
    assert row["action"] == "refresh_replay"
    assert row["route"] == "zero_variant_fulltext_qc"
    assert row["zero_variant_qc"] is True
    assert row["missing_variant_supplement"] is False
    assert summary["pmid_coverage"]["zero_variant_usable_fulltext"]["pmids"] == 1
    assert _read_csv(override) == [
        {
            "gene": "KCNH2",
            "pmid": pmid,
            "action": "refresh_replay",
            "route": "zero_variant_fulltext_qc",
            "available_context_path": str(full_context),
            "available_context_bytes": str(full_context.stat().st_size),
            "notes": "usable full text produced zero variants",
        }
    ]


def test_audit_ignores_raw_discovery_pmids_by_default(tmp_path: Path):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    harvest_dir.mkdir(parents=True)
    extraction_dir.mkdir()
    (run_dir / "KCNH2_pmids.txt").write_text("999\n", encoding="utf-8")

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    discovery_rows, discovery_summary = audit.build_audit(
        gene="KCNH2",
        run_dir=run_dir,
        include_discovery_pmids=True,
    )

    assert rows == []
    assert summary["pmids"] == 0
    assert summary["include_discovery_pmids"] is False
    assert [row["pmid"] for row in discovery_rows] == ["999"]
    assert discovery_rows[0]["action"] == "fetch"
    assert discovery_summary["pmids"] == 1
    assert discovery_summary["include_discovery_pmids"] is True


def test_audit_fetches_missing_target_gene_variant_supplement(tmp_path: Path):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    harvest_dir.mkdir(parents=True)
    extraction_dir.mkdir()

    pmid = "555"
    full_context = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    full_context.write_text(
        "# Article\n\nAbstract\nMethods\nResults\n\n"
        "The probands underwent KCNH2 genetic testing. "
        "All KCNH2 mutations are listed in Supplemental Table 2.\n\n"
        "Table 1. Cohort demographics\n\n| Group | N |\n| --- | --- |\n| A | 1 |\n"
        + ("main text without variant identifiers. " * 220),
        encoding="utf-8",
    )
    extraction = extraction_dir / f"KCNH2_PMID_{pmid}.json"
    extraction.write_text(
        json.dumps(
            {
                "variants": [{"protein_notation": "p.Arg1His"}],
                "extraction_metadata": {
                    "pmid": pmid,
                    "source_file": str(full_context),
                    "source_sha256": audit._sha256(full_context),
                },
            }
        ),
        encoding="utf-8",
    )
    result_dir = harvest_dir / pmid
    result_dir.mkdir()
    (result_dir / "result.json").write_text(
        json.dumps({"pmid": pmid, "doi": "10.1016/j.hrthm.2020.01.001"}),
        encoding="utf-8",
    )

    rows, summary = audit.build_audit(gene="KCNH2", run_dir=run_dir)
    fetch_input = tmp_path / "fetch_input.csv"
    supplement_input = tmp_path / "supplement_input.csv"
    source_override = tmp_path / "source_override.csv"
    audit.write_fetch_input(rows, fetch_input)
    audit.write_supplement_input(rows, supplement_input)
    audit.write_source_override(rows, source_override)

    assert len(rows) == 1
    row = rows[0]
    assert row["action"] == "fetch_supplement_only"
    assert row["route"] == "fetch_elsevier_supplements_only"
    assert row["source_status"] == "recovered_pmc"
    assert row["missing_variant_supplement"] is True
    assert (
        row["notes"]
        == "usable main text references a missing target-gene variant supplement"
    )
    assert summary["pmid_coverage"]["missing_variant_supplement"]["pmids"] == 1
    assert _read_csv(fetch_input) == []
    assert _read_csv(supplement_input) == [
        {
            "PMID": "555",
            "DOI": "10.1016/j.hrthm.2020.01.001",
            "route": "fetch_elsevier_supplements_only",
        }
    ]
    assert _read_csv(source_override) == []

    wiley = audit.classify_pmid(
        gene="KCNH2",
        pmid=pmid,
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        doi="10.1002/humu.21126",
        zero_variant_pmids=set(),
        single_carrier_pmids=set(),
    )
    assert wiley["action"] == "fetch"
    assert wiley["route"] == "fetch_wiley_tdm"
