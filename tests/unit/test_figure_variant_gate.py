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

import json
import sqlite3
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from scripts.extract_figure_variants import (  # noqa: E402
    FIGURE_VARIANT_GATE_ENV,
    _discover_pmids_with_figures,
    _figure_variant_passes_gate,
    ingest_cached_variants,
)

GENE = "KCNH2"  # protein length 1159 in PROTEIN_LENGTHS


def test_discovers_nested_corpus_figures(tmp_path: Path):
    fig_dir = tmp_path / "12345" / "12345_figures"
    fig_dir.mkdir(parents=True)
    (fig_dir / "figure.png").write_bytes(b"\x00")

    assert _discover_pmids_with_figures(tmp_path) == ["12345"]


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
        CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            total_carriers_observed INTEGER,
            affected_count INTEGER,
            unaffected_count INTEGER
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


def test_ingest_enriches_an_existing_variant_paper_link(tmp_path, monkeypatch):
    monkeypatch.delenv(FIGURE_VARIANT_GATE_ENV, raising=False)
    db = _make_db(tmp_path / "existing_link.db")
    con = sqlite3.connect(db)
    con.execute(
        "INSERT INTO variants(gene_symbol, cdna_notation, protein_notation) VALUES(?,?,?)",
        (GENE, "c.1682C>T", "p.Ala561Val"),
    )
    variant_id = con.execute("SELECT variant_id FROM variants").fetchone()[0]
    con.execute(
        """INSERT INTO variant_papers(
               variant_id, pmid, source_location, additional_notes, key_quotes
           ) VALUES(?,?,?,?,?)""",
        (variant_id, "111", "Table 2", "text extraction note", "[]"),
    )
    con.commit()
    con.close()

    changed = ingest_cached_variants(
        pmid="111",
        gene=GENE,
        distinct=[
            {
                "cdna": "c.1682C>T",
                "protein": "p.Ala561Val",
                "carriers": 13,
                "affected": 5,
                "context": "pedigree carrier counts",
            }
        ],
        db_path=db,
    )

    con = sqlite3.connect(db)
    notes, layer, provenance = con.execute(
        "SELECT additional_notes, source_layer, count_provenance FROM variant_papers"
    ).fetchone()
    con.close()
    assert changed == 1
    assert json.loads(notes)["prior_notes"] == "text extraction note"
    assert json.loads(notes)["carriers"] == 13
    assert layer == "figure"
    assert json.loads(provenance)["carriers_count_type"] == "per_variant_carriers"


def test_ingest_reuses_a_compatible_partial_variant_identity(tmp_path, monkeypatch):
    monkeypatch.delenv(FIGURE_VARIANT_GATE_ENV, raising=False)
    db = _make_db(tmp_path / "partial_identity.db")
    con = sqlite3.connect(db)
    con.execute(
        "INSERT INTO variants(gene_symbol, protein_notation) VALUES(?,?)",
        (GENE, "p.Ala561Val"),
    )
    con.commit()
    con.close()

    ingest_cached_variants(
        pmid="111",
        gene=GENE,
        distinct=[
            {
                "cdna": "c.1682C>T",
                "protein": "p.Ala561Val",
                "context": "pedigree carrier counts",
            }
        ],
        db_path=db,
    )

    con = sqlite3.connect(db)
    rows = con.execute(
        "SELECT variant_id, cdna_notation, protein_notation FROM variants"
    ).fetchall()
    links = con.execute("SELECT variant_id, pmid FROM variant_papers").fetchall()
    con.close()
    assert rows == [(1, "c.1682C>T", "p.Ala561Val")]
    assert links == [(1, "111")]


def test_ingest_preserves_existing_figure_counts_while_filling_nulls(
    tmp_path, monkeypatch
):
    monkeypatch.delenv(FIGURE_VARIANT_GATE_ENV, raising=False)
    db = _make_db(tmp_path / "fill_null_notes.db")
    con = sqlite3.connect(db)
    con.execute(
        "INSERT INTO variants(gene_symbol, cdna_notation, protein_notation) VALUES(?,?,?)",
        (GENE, "c.1682C>T", "p.Ala561Val"),
    )
    variant_id = con.execute("SELECT variant_id FROM variants").fetchone()[0]
    con.execute(
        """INSERT INTO variant_papers(
               variant_id, pmid, source_location, additional_notes, key_quotes
           ) VALUES(?,?,?,?,?)""",
        (
            variant_id,
            "111",
            "prior figure",
            json.dumps({"carriers": 8, "affected": None, "context": "pedigree"}),
            "[]",
        ),
    )
    con.commit()
    con.close()

    ingest_cached_variants(
        pmid="111",
        gene=GENE,
        distinct=[
            {
                "cdna": "c.1682C>T",
                "protein": "p.Ala561Val",
                "carriers": 3,
                "affected": 1,
                "context": "pedigree reread",
            }
        ],
        db_path=db,
    )

    con = sqlite3.connect(db)
    notes = json.loads(
        con.execute("SELECT additional_notes FROM variant_papers").fetchone()[0]
    )
    con.close()
    assert notes["carriers"] == 8
    assert notes["affected"] == 1
    assert notes["context"] == "pedigree reread"


def test_ingest_does_not_relabel_an_existing_text_count_as_figure(
    tmp_path, monkeypatch
):
    monkeypatch.delenv(FIGURE_VARIANT_GATE_ENV, raising=False)
    db = _make_db(tmp_path / "provenance_scope.db")
    con = sqlite3.connect(db)
    con.execute(
        "INSERT INTO variants(gene_symbol, cdna_notation, protein_notation) VALUES(?,?,?)",
        (GENE, "c.1682C>T", "p.Ala561Val"),
    )
    variant_id = con.execute("SELECT variant_id FROM variants").fetchone()[0]
    con.execute(
        """INSERT INTO variant_papers(
               variant_id, pmid, source_location, additional_notes, key_quotes
           ) VALUES(?,?,?,?,?)""",
        (variant_id, "111", "Table 2", "text extraction", "[]"),
    )
    con.execute(
        """INSERT INTO penetrance_data(
               variant_id, pmid, total_carriers_observed, affected_count
           ) VALUES(?,?,?,?)""",
        (variant_id, "111", None, 5),
    )
    con.commit()
    con.close()

    ingest_cached_variants(
        pmid="111",
        gene=GENE,
        distinct=[
            {
                "cdna": "c.1682C>T",
                "protein": "p.Ala561Val",
                "carriers": 13,
                "affected": 3,
                "context": "pedigree reread",
            }
        ],
        db_path=db,
    )

    con = sqlite3.connect(db)
    provenance = json.loads(
        con.execute("SELECT count_provenance FROM variant_papers").fetchone()[0]
    )
    con.close()
    assert provenance["carriers_source"] == "figure"
    assert "affected_source" not in provenance


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
