"""Unit tests for the figure-reader precision gate.

These tests are hermetic: no vision/LLM call, no network, no real run dirs.
We feed ``ingest_cached_variants`` a synthetic list of variant dicts against a
throwaway on-disk sqlite (sqlite cannot share an in-memory DB across the two
connections the function opens, so we use ``tmp_path``).

Covered:
  * validate (default): out-of-range protein position is dropped.
  * validate (default): malformed notation is dropped.
  * validate (default): a valid in-range variant passes.
  * off: everything passes (the historical raw behavior).
  * the gate function itself + the dropped-reason breakdown.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from scripts.extract_figure_variants import (  # noqa: E402
    FIGURE_VARIANT_GATE_ENV,
    _figure_variant_passes_gate,
    ingest_cached_variants,
)

GENE = "KCNH2"  # protein length 1159 in PROTEIN_LENGTHS


def _make_db(path: Path) -> Path:
    """Create the minimal schema ingest_cached_variants touches."""
    con = sqlite3.connect(str(path))
    con.executescript(
        """
        CREATE TABLE papers (
            pmid TEXT PRIMARY KEY,
            gene_symbol TEXT,
            extraction_summary TEXT
        );
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL,
            cdna_notation TEXT,
            protein_notation TEXT,
            genomic_position TEXT,
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
        """
    )
    con.commit()
    con.close()
    return path


def _inserted_variants(path: Path):
    """Return (cdna, protein, source_location) rows actually persisted."""
    con = sqlite3.connect(str(path))
    try:
        return con.execute(
            """SELECT v.cdna_notation, v.protein_notation, vp.source_location
               FROM variant_papers vp
               JOIN variants v ON v.variant_id = vp.variant_id
               ORDER BY v.variant_id"""
        ).fetchall()
    finally:
        con.close()


# Synthetic figure-reader output: one good, one out-of-range, one malformed.
SYNTHETIC = [
    {"cdna": "c.1682C>T", "protein": "p.Ala561Val", "context": "good in-range"},
    {"cdna": None, "protein": "p.Arg9999Cys", "context": "position out of range"},
    {"cdna": None, "protein": "notavariant", "context": "malformed"},
]


# ---------------------------------------------------------------------------
# Gate function (pure, no DB)
# ---------------------------------------------------------------------------
def test_gate_validate_passes_in_range():
    passed, reason = _figure_variant_passes_gate(
        GENE, "c.1682C>T", "p.Ala561Val", mode="validate"
    )
    assert passed is True
    assert reason == "ok"


def test_gate_validate_drops_out_of_range_position():
    passed, reason = _figure_variant_passes_gate(
        GENE, None, "p.Arg9999Cys", mode="validate"
    )
    assert passed is False
    # is_non_target_variant catches position>length first.
    assert reason in {"non_target", "position_out_of_range"}


def test_gate_validate_drops_malformed_notation():
    passed, reason = _figure_variant_passes_gate(
        GENE, "definitely not cdna", "notavariant", mode="validate"
    )
    assert passed is False
    assert reason == "malformed_notation"


def test_gate_off_passes_everything():
    for cdna, protein in [
        ("c.1682C>T", "p.Ala561Val"),
        (None, "p.Arg9999Cys"),
        (None, "notavariant"),
    ]:
        passed, reason = _figure_variant_passes_gate(GENE, cdna, protein, mode="off")
        assert passed is True, (cdna, protein)
        assert reason == "gate_off"


def test_gate_unknown_gene_does_not_overreject():
    # A gene with no known length must not be filtered for position bounds.
    passed, reason = _figure_variant_passes_gate(
        "ZZZ9", None, "p.Arg9999Cys", mode="validate"
    )
    assert passed is True
    assert reason == "ok"


def test_gate_corroborate_behaves_like_validate_for_now():
    # Stub: corroborate currently == validate (passes valid, drops bad).
    ok, _ = _figure_variant_passes_gate(
        GENE, "c.1682C>T", "p.Ala561Val", mode="corroborate"
    )
    assert ok is True
    bad, reason = _figure_variant_passes_gate(
        GENE, None, "p.Arg9999Cys", mode="corroborate"
    )
    assert bad is False
    assert reason in {"non_target", "position_out_of_range"}


# ---------------------------------------------------------------------------
# End-to-end ingest against a temp sqlite
# ---------------------------------------------------------------------------
def test_ingest_validate_default_filters_bad_variants(tmp_path, monkeypatch):
    # Default behavior (no env set) is ``validate``.
    monkeypatch.delenv(FIGURE_VARIANT_GATE_ENV, raising=False)
    db = _make_db(tmp_path / "validate_default.db")

    added = ingest_cached_variants(
        pmid="111", gene=GENE, distinct=list(SYNTHETIC), db_path=db
    )

    assert added == 1, "only the valid in-range variant should be inserted"
    rows = _inserted_variants(db)
    assert len(rows) == 1
    cdna, protein, source = rows[0]
    assert cdna == "c.1682C>T"
    assert protein == "p.Ala561Val"
    # Gated-in variants stay traceable as figure-reader output.
    assert source == "figure-reader (cached)"


def test_ingest_explicit_validate_matches_default(tmp_path, monkeypatch):
    monkeypatch.setenv(FIGURE_VARIANT_GATE_ENV, "validate")
    db = _make_db(tmp_path / "validate_explicit.db")
    added = ingest_cached_variants(
        pmid="111", gene=GENE, distinct=list(SYNTHETIC), db_path=db
    )
    assert added == 1


def test_ingest_gate_off_inserts_everything(tmp_path, monkeypatch):
    monkeypatch.setenv(FIGURE_VARIANT_GATE_ENV, "off")
    db = _make_db(tmp_path / "gate_off.db")

    added = ingest_cached_variants(
        pmid="111", gene=GENE, distinct=list(SYNTHETIC), db_path=db
    )

    # All three have a non-empty cdna OR protein, so off == raw behavior.
    assert added == 3
    rows = _inserted_variants(db)
    assert len(rows) == 3
    proteins = {protein for _cdna, protein, _src in rows}
    assert {"p.Ala561Val", "p.Arg9999Cys", "notavariant"} == proteins


def test_ingest_unknown_gate_mode_falls_back_to_validate(tmp_path, monkeypatch):
    monkeypatch.setenv(FIGURE_VARIANT_GATE_ENV, "bogusmode")
    db = _make_db(tmp_path / "bogus.db")
    added = ingest_cached_variants(
        pmid="111", gene=GENE, distinct=list(SYNTHETIC), db_path=db
    )
    assert added == 1, "unknown mode must fall back to the safe validate default"


def test_ingest_empty_notation_dropped(tmp_path, monkeypatch):
    monkeypatch.delenv(FIGURE_VARIANT_GATE_ENV, raising=False)
    db = _make_db(tmp_path / "empty.db")
    added = ingest_cached_variants(
        pmid="111",
        gene=GENE,
        distinct=[{"cdna": "", "protein": "  ", "context": "blank"}],
        db_path=db,
    )
    assert added == 0
    assert _inserted_variants(db) == []


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
