"""Bounded end-to-end regression tests for the recall stack.

These tests don't run the LLM-bound extraction step (that costs minutes per
PMID and 1000s of dollars at scale). Instead they exercise the *scoring
pipeline + recovery driver wiring* against a tiny in-memory fixture so a
break in the scorer, driver, or recall suite gets caught in seconds, before
a full multi-hour run.

The expensive end-to-end path is covered separately by the
``requires_network`` integration tests in tests/integration/.
"""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
RECALL_RECOVERY_DIR = REPO_ROOT / "scripts" / "recall_recovery"
EXTRACT_FIGURES_SCRIPT = REPO_ROOT / "scripts" / "extract_figure_variants.py"
RUN_RECALL_SUITE = REPO_ROOT / "scripts" / "run_recall_suite.py"


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _build_fixture_db(db_path: Path) -> None:
    """Create a tiny SQLite DB with the schema run_recall_suite expects."""
    con = sqlite3.connect(str(db_path))
    cur = con.cursor()
    # Minimal schema (mirrors harvesting/migrate_to_sqlite.py)
    cur.executescript(
        """
        CREATE TABLE papers (
            pmid TEXT PRIMARY KEY,
            title TEXT,
            journal TEXT,
            publication_date TEXT,
            doi TEXT,
            pmc_id TEXT,
            gene_symbol TEXT,
            extraction_summary TEXT,
            extraction_timestamp TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        );
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL,
            cdna_notation TEXT,
            protein_notation TEXT,
            genomic_position TEXT,
            clinical_significance TEXT,
            evidence_level TEXT,
            UNIQUE(gene_symbol, cdna_notation, protein_notation, genomic_position)
        );
        CREATE TABLE variant_papers (
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            source_location TEXT,
            additional_notes TEXT,
            key_quotes TEXT,
            PRIMARY KEY (variant_id, pmid)
        );
        CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            total_carriers_observed INTEGER,
            affected_count INTEGER,
            unaffected_count INTEGER,
            uncertain_count INTEGER,
            penetrance_percentage REAL
        );
        CREATE TABLE individual_records (
            record_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            individual_id TEXT,
            affected_status TEXT,
            age INTEGER,
            sex TEXT,
            additional_notes TEXT,
            confidence REAL,
            extraction_source TEXT,
            family_id TEXT
        );
        """
    )
    # Seed: 2 papers, 3 variants, 2 of which match the gold below.
    cur.execute(
        "INSERT INTO papers (pmid, gene_symbol, extraction_summary) VALUES (?, ?, ?)",
        ("11111111", "TESTGENE", "fixture paper 1"),
    )
    cur.execute(
        "INSERT INTO papers (pmid, gene_symbol, extraction_summary) VALUES (?, ?, ?)",
        ("22222222", "TESTGENE", "fixture paper 2"),
    )
    # R100W is in both DB and gold for pmid 11111111 — match
    cur.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES (?, ?, ?)",
        ("TESTGENE", None, "p.Arg100Trp"),
    )
    cur.execute(
        "INSERT INTO variant_papers (variant_id, pmid, source_location) VALUES (?, ?, ?)",
        (1, "11111111", "fixture"),
    )
    cur.execute(
        "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, affected_count, unaffected_count) VALUES (?, ?, ?, ?, ?)",
        (1, "11111111", 5, 3, 2),
    )
    # L200P is in DB for pmid 22222222 but gold has it for a different pmid — won't match
    cur.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES (?, ?, ?)",
        ("TESTGENE", None, "p.Leu200Pro"),
    )
    cur.execute(
        "INSERT INTO variant_papers (variant_id, pmid, source_location) VALUES (?, ?, ?)",
        (2, "22222222", "fixture"),
    )
    # G300S is in DB and gold for pmid 22222222 — match
    cur.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES (?, ?, ?)",
        ("TESTGENE", None, "p.Gly300Ser"),
    )
    cur.execute(
        "INSERT INTO variant_papers (variant_id, pmid, source_location) VALUES (?, ?, ?)",
        (3, "22222222", "fixture"),
    )
    cur.execute(
        "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, affected_count) VALUES (?, ?, ?, ?)",
        (3, "22222222", 4, 4),
    )
    con.commit()
    con.close()


def _build_fixture_gold(gold_dir: Path) -> Path:
    """Create a gold-standard package containing one gene CSV + a manifest."""
    normalized = gold_dir / "normalized"
    normalized.mkdir(parents=True, exist_ok=True)
    csv_path = normalized / "TESTGENE_recall_input.csv"
    csv_path.write_text(
        "variant,pmid,carriers,affected,unaffected\n"
        # In DB + matches:
        "R100W,11111111,5,3,2\n"
        "G300S,22222222,4,4,0\n"
        # NOT in DB — fixture-side misses (should show up as missing_in_sqlite):
        "Y50N,11111111,1,1,0\n"
        "T400I,33333333,2,1,1\n"
    )
    # Bare-minimum manifest so run_recall_suite.discover_gold_inputs is happy.
    (gold_dir / "manifest.json").write_text(json.dumps({"genes": {"TESTGENE": {}}}))
    return csv_path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_scoring_pipeline_round_trip(tmp_path: Path):
    """The recall suite scores a known fixture DB to the expected numbers.

    Catches breaks in cli.compare_variants, run_recall_suite output shape,
    schema assumptions, or normalization regressions.
    """
    db_path = tmp_path / "TESTGENE.db"
    gold_dir = tmp_path / "gold"
    out_dir = tmp_path / "score_out"
    _build_fixture_db(db_path)
    _build_fixture_gold(gold_dir)

    result = subprocess.run(
        [
            sys.executable,
            str(RUN_RECALL_SUITE),
            "--score",
            "--genes",
            "TESTGENE",
            "--db",
            f"TESTGENE={db_path}",
            "--gold-dir",
            str(gold_dir),
            "--outdir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"recall suite failed: {result.stderr[-400:]}"
    summary_file = out_dir / "summary.json"
    assert summary_file.exists(), "summary.json was not written"
    data = json.loads(summary_file.read_text())
    agg = data["aggregate_recall"]

    # 2 of the 3 unique gold variants match (R100W, G300S match; Y50N + T400I miss).
    # 2 of the 3 distinct gold pmids match (11111111 + 22222222 match; 33333333 misses).
    assert agg["pmids"]["matched"] == 2
    assert agg["pmids"]["gold"] == 3
    # variant_rows: 4 gold rows, 2 matched (R100W + G300S)
    assert agg["variant_rows"]["matched"] == 2
    assert agg["variant_rows"]["gold"] == 4
    # patients = sum carriers across matched rows = 5 + 4 = 9, gold total = 5+4+1+2 = 12
    assert agg["patients"]["matched"] == 9
    assert agg["patients"]["gold"] == 12


def test_run_all_layers_skip_everything(tmp_path: Path):
    """Driver smoke test: with every layer skipped, it should just emit the
    baseline score.

    Catches breaks in the driver's argparse, subprocess wiring, score-parse,
    and progression file emission.
    """
    db_path = tmp_path / "TESTGENE.db"
    gold_dir = tmp_path / "gold"
    out_dir = tmp_path / "driver_out"
    _build_fixture_db(db_path)
    gold_csv = _build_fixture_gold(gold_dir)

    result = subprocess.run(
        [
            sys.executable,
            str(RECALL_RECOVERY_DIR / "run_all_layers.py"),
            "--gene",
            "TESTGENE",
            "--db",
            str(db_path),
            "--gold",
            str(gold_csv),
            "--outdir",
            str(out_dir),
            "--skip",
            "clinvar",
            "--skip",
            "pubtator",
            "--skip",
            "figures",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"driver failed: {result.stderr[-400:]}"

    progression = out_dir / "progression.json"
    assert progression.exists()
    payload = json.loads(progression.read_text())
    rows = payload["progression"]
    assert len(rows) == 1, f"expected only baseline row, got {rows}"
    baseline = rows[0]
    assert baseline["layer"] == "0_baseline"
    # 2 of 3 PMIDs matched → 66.67% (rounded)
    assert baseline["pmids"] == pytest.approx(66.67, abs=0.5)
    # 2 of 4 variant_rows matched → 50.0%
    assert baseline["variant_rows"] == pytest.approx(50.0, abs=0.5)


def test_extract_figure_variants_help():
    """The figure-reader CLI parses --help cleanly. Catches argparse regressions."""
    result = subprocess.run(
        [sys.executable, str(EXTRACT_FIGURES_SCRIPT), "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--auto-pmids" in result.stdout, "auto-pmids flag missing from CLI"
    assert "--force" in result.stdout, "force flag missing from CLI"


def test_run_all_layers_help():
    """The driver CLI parses --help cleanly."""
    result = subprocess.run(
        [sys.executable, str(RECALL_RECOVERY_DIR / "run_all_layers.py"), "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--with-v12" in result.stdout, "with-v12 flag missing from CLI"
    assert "--gene" in result.stdout


# ---------------------------------------------------------------------------
# gvf-run wiring tests
# ---------------------------------------------------------------------------


def test_gvf_run_help():
    """The gvf-run command parses --help cleanly."""
    result = subprocess.run(
        [sys.executable, "-m", "cli", "gvf-run", "--help"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert result.returncode == 0
    assert "--pmid-file" in result.stdout
    assert "--resume-dir" in result.stdout
    assert "--skip" in result.stdout
    assert "source-qc" in result.stdout
    assert "--source-recovery" in result.stdout


def test_gvf_run_helpers_resolve_paths(tmp_path: Path):
    """gvf-run's auto-detect helpers behave correctly on a synthetic layout."""
    from cli.gvf_run import _find_db, _find_gold, _find_v12_db

    # _find_db picks the gene-named DB when present
    run_dir = tmp_path / "RUN"
    run_dir.mkdir()
    (run_dir / "FOO.db").write_bytes(b"")
    assert _find_db(run_dir, "FOO") == run_dir / "FOO.db"

    # _find_db falls back to latest *.db when no gene-named DB
    (run_dir / "FOO.db").unlink()
    (run_dir / "other.db").write_bytes(b"")
    assert _find_db(run_dir, "FOO") == run_dir / "other.db"

    # _find_db returns None when no DB at all
    (run_dir / "other.db").unlink()
    assert _find_db(run_dir, "FOO") is None

    # _find_v12_db refuses non-KCNH2 genes
    assert _find_v12_db("KCNQ1") is None
    assert _find_v12_db("RYR2") is None
