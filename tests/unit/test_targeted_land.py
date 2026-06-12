"""Tests for the fast targeted-land candidate scan.

The scan must pick only papers that look over-counted-but-unpaired (where a
deterministic re-parse adds cDNA+protein pairing) and skip already-well-paired or
small extractions, using intrinsic signals only — no gold/PMID/gene literals.
"""

import json

import scripts.targeted_land as T


def _write(run_dir, corpus, gene, pmid, variants):
    (run_dir / "extractions").mkdir(parents=True, exist_ok=True)
    (run_dir / "extractions" / f"{gene}_PMID_{pmid}.json").write_text(
        json.dumps({"variants": variants}), encoding="utf-8"
    )
    src_dir = corpus / gene / pmid
    src_dir.mkdir(parents=True, exist_ok=True)
    (src_dir / f"{pmid}_FULL_CONTEXT.md").write_text("source", encoding="utf-8")


def test_scan_selects_overcounted_cdna_only_skips_others(tmp_path, monkeypatch):
    gene = "KCNQ1"
    run_dir = tmp_path / "run"
    corpus = tmp_path / "corpus"
    # 111: 20 cDNA-only rows, 0 paired -> suspicious (re-parse would add pairing).
    _write(
        run_dir,
        corpus,
        gene,
        "111",
        [{"cdna_notation": f"c.{100 + i}A>G"} for i in range(20)],
    )
    # 222: 20 fully-paired rows -> already good, skipped by the pre-filter.
    _write(
        run_dir,
        corpus,
        gene,
        "222",
        [
            {"cdna_notation": f"c.{100 + i}A>G", "protein_notation": f"p.X{i}Y"}
            for i in range(20)
        ],
    )
    # 333: tiny extraction -> below the suspicious threshold, skipped.
    _write(run_dir, corpus, gene, "333", [{"cdna_notation": "c.1A>G"}])

    # Deterministic re-parse returns 6 paired variants for any suspicious source.
    monkeypatch.setattr(
        T.RDB,
        "deterministic_variant_list",
        lambda _ext, _src, _gene: [
            {"cdna_notation": f"c.{100 + i}A>G", "protein_notation": f"p.X{i}Y"}
            for i in range(6)
        ],
    )

    candidates = T.scan_quality_lift_candidates(gene, run_dir, corpus)
    assert candidates == ["111"]


def test_scan_skips_when_reparse_adds_no_pairing(tmp_path, monkeypatch):
    gene = "KCNQ1"
    run_dir = tmp_path / "run"
    corpus = tmp_path / "corpus"
    _write(
        run_dir,
        corpus,
        gene,
        "111",
        [{"cdna_notation": f"c.{100 + i}A>G"} for i in range(20)],
    )
    # Deterministic re-parse finds only cDNA-only variants too -> no pairing lift.
    monkeypatch.setattr(
        T.RDB,
        "deterministic_variant_list",
        lambda _ext, _src, _gene: [
            {"cdna_notation": f"c.{100 + i}A>G"} for i in range(8)
        ],
    )
    assert T.scan_quality_lift_candidates(gene, run_dir, corpus) == []
