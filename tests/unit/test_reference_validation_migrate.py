"""Tests for the B1 reference-validation gate wired into DB migration.

Covers ``harvesting.migrate_to_sqlite._apply_reference_validation`` (the helper
that runs over a variant list right before insertion) and the end-to-end
``migrate_extraction_file`` path under ``REFERENCE_VALIDATION_POLICY=drop``, plus
the new ``Settings.reference_validation_policy`` knob.

The helper-level tests exercise the real committed KCNH2 reference cache
(``data/reference_sequences/KCNH2.fasta`` == NP_000229.1), so the verdicts are
ground-truth: position 561 is Ala (A561V is OK), 628 is Gly (G628S is OK), 178 is
Ser and 341 is Val (so A178T / A341E are reference-residue mismatches). These are
stable: the manifest pins the accession and the loader cross-checks its length.
"""

import json

import pytest

import harvesting.migrate_to_sqlite as migrate
from harvesting.migrate_to_sqlite import (
    _apply_reference_validation,
    create_database_schema,
    migrate_extraction_file,
)
from pipeline.reference_validation import clear_cache


@pytest.fixture(autouse=True)
def _fresh_ref_cache():
    """Drop the in-process reference-sequence cache around each test."""
    clear_cache()
    yield
    clear_cache()


def _kcnh2(protein):
    return {"gene_symbol": "KCNH2", "protein_notation": protein}


def test_off_is_strict_noop():
    variants = [_kcnh2("A178T"), _kcnh2("A561V")]
    kept, stats = _apply_reference_validation(variants, policy="off")
    assert kept is variants  # same object, untouched
    assert stats == {"ok": 0, "unknown": 0, "reject": 0, "flagged": 0, "dropped": 0}
    assert all("reference_validation" not in v for v in kept)


def test_flag_annotates_mismatch_but_keeps_all():
    variants = [_kcnh2("A178T"), _kcnh2("A561V"), _kcnh2("c.1681G>A")]
    kept, stats = _apply_reference_validation(variants, policy="flag")
    assert len(kept) == 3  # nothing dropped in flag mode
    assert stats["ok"] == 1  # A561V
    assert stats["unknown"] == 1  # cDNA-only -> no residue to check
    assert stats["reject"] == 1 and stats["flagged"] == 1 and stats["dropped"] == 0
    flagged = [v for v in kept if "reference_validation" in v]
    assert len(flagged) == 1
    rv = flagged[0]["reference_validation"]
    assert rv["status"] == "reject"
    assert rv["reason"] == "reference_residue_mismatch"
    assert rv["expected"] == "S" and rv["observed"] == "A" and rv["policy"] == "flag"


def test_drop_removes_only_mismatches():
    variants = [
        _kcnh2("A178T"),  # pos 178 is S -> reject
        _kcnh2("A341E"),  # pos 341 is V -> reject
        _kcnh2("A561V"),  # pos 561 is A -> ok
        _kcnh2("G628S"),  # pos 628 is G -> ok
        _kcnh2("c.1681G>A"),  # cDNA-only -> unknown -> pass
    ]
    kept, stats = _apply_reference_validation(variants, policy="drop")
    kept_proteins = {v["protein_notation"] for v in kept}
    assert kept_proteins == {"A561V", "G628S", "c.1681G>A"}
    assert stats["reject"] == 2 and stats["dropped"] == 2 and stats["flagged"] == 2
    assert stats["ok"] == 2 and stats["unknown"] == 1


def test_drop_spares_low_confidence_rejects():
    """drop omits only high-confidence rejects; fs + N-terminal rejects are kept."""
    variants = [
        _kcnh2("A178T"),  # simple sub, pos 178 -> high-confidence -> DROPPED
        _kcnh2("p.Gln1122fs*147"),  # frameshift reject (exp P) -> kept
        _kcnh2("p.Ala8Val"),  # N-terminal reject (pos 8) -> kept
        _kcnh2("A561V"),  # ok
    ]
    kept, stats = _apply_reference_validation(variants, policy="drop")
    kept_proteins = {v["protein_notation"] for v in kept}
    assert kept_proteins == {"p.Gln1122fs*147", "p.Ala8Val", "A561V"}
    assert stats["reject"] == 3  # A178T + fs + A8V all mismatch the reference
    assert stats["flagged"] == 3  # all three rejects are annotated
    assert stats["dropped"] == 1  # but only the high-confidence one is removed
    # The spared rejects are annotated high_confidence=False, dropped=False.
    spared = [v for v in kept if "reference_validation" in v]
    assert {v["protein_notation"] for v in spared} == {
        "p.Gln1122fs*147",
        "p.Ala8Val",
    }
    for v in spared:
        rv = v["reference_validation"]
        assert rv["status"] == "reject"
        assert rv["high_confidence"] is False and rv["dropped"] is False


def test_uncached_gene_always_passes():
    # No FASTA for this gene -> every verdict is 'unknown' -> nothing dropped.
    variants = [
        {"gene_symbol": "ZZ_NOT_A_GENE", "protein_notation": "A178T"},
        {"gene_symbol": "ZZ_NOT_A_GENE", "protein_notation": "R500W"},
    ]
    kept, stats = _apply_reference_validation(variants, policy="drop")
    assert len(kept) == 2 and stats["dropped"] == 0
    assert stats["unknown"] == 2 and stats["reject"] == 0


def test_migrate_extraction_file_drop_excludes_mismatch(tmp_path, monkeypatch):
    """End-to-end: drop policy keeps the mismatch out of the variants table."""
    monkeypatch.setattr(migrate, "_resolve_reference_validation_policy", lambda: "drop")
    json_file = tmp_path / "PMID_12345_KCNH2_extraction.json"
    json_file.write_text(
        json.dumps(
            {
                "paper_metadata": {"pmid": "12345", "title": "KCNH2 variants"},
                "variants": [
                    {"gene_symbol": "KCNH2", "protein_notation": "A561V"},  # ok
                    {"gene_symbol": "KCNH2", "protein_notation": "A178T"},  # reject
                ],
            }
        ),
        encoding="utf-8",
    )

    conn = create_database_schema(str(tmp_path / "variants.db"))
    cursor = conn.cursor()
    ok, msg = migrate_extraction_file(cursor, json_file)
    conn.commit()
    assert ok, msg

    rows = cursor.execute("SELECT protein_notation FROM variants").fetchall()
    proteins = {r[0] for r in rows}
    assert len(rows) == 1, f"expected only the OK variant, got {proteins}"
    assert any("561" in (p or "") for p in proteins)
    assert not any("178" in (p or "") for p in proteins)
    conn.close()


def test_migrate_extraction_file_off_keeps_all(tmp_path, monkeypatch):
    """Default (off) policy inserts both variants, including the mismatch."""
    monkeypatch.setattr(migrate, "_resolve_reference_validation_policy", lambda: "off")
    json_file = tmp_path / "PMID_12345_KCNH2_extraction.json"
    json_file.write_text(
        json.dumps(
            {
                "paper_metadata": {"pmid": "12345", "title": "KCNH2 variants"},
                "variants": [
                    {"gene_symbol": "KCNH2", "protein_notation": "A561V"},
                    {"gene_symbol": "KCNH2", "protein_notation": "A178T"},
                ],
            }
        ),
        encoding="utf-8",
    )
    conn = create_database_schema(str(tmp_path / "variants.db"))
    cursor = conn.cursor()
    ok, _ = migrate_extraction_file(cursor, json_file)
    conn.commit()
    assert ok
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


def _baseline_env(monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "user@example.org")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-anthropic")


def test_settings_reference_policy_default_off(monkeypatch):
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.delenv("REFERENCE_VALIDATION_POLICY", raising=False)
    get_settings.cache_clear()
    try:
        assert get_settings().reference_validation_policy == "off"
    finally:
        get_settings.cache_clear()


def test_settings_reference_policy_accepts_flag_and_drop(monkeypatch):
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    for raw, want in (("FLAG", "flag"), ("drop", "drop"), (" Drop ", "drop")):
        monkeypatch.setenv("REFERENCE_VALIDATION_POLICY", raw)
        get_settings.cache_clear()
        try:
            assert get_settings().reference_validation_policy == want
        finally:
            get_settings.cache_clear()


def test_settings_reference_policy_rejects_unknown(monkeypatch):
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.setenv("REFERENCE_VALIDATION_POLICY", "nonsense")
    get_settings.cache_clear()
    try:
        with pytest.raises(ValueError):
            get_settings()
    finally:
        get_settings.cache_clear()
