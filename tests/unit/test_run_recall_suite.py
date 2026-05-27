"""Tests for the canonical multi-gene recall runner."""

import sqlite3

from scripts.run_recall_suite import (
    build_gene_results,
    combine_mae,
    combine_recall,
    discover_gold_inputs,
)


def _write_recall_csv(path):
    path.parent.mkdir(parents=True)
    path.write_text(
        "variant,pmid,carriers,affected,unaffected\n"
        "p.Arg1His,111,2,1,1\n"
        "p.Val2Met,222,1,1,0\n",
        encoding="utf-8",
    )


def _write_gvf_db(path):
    path.parent.mkdir(parents=True)
    conn = sqlite3.connect(path)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE variants (
                variant_id INTEGER PRIMARY KEY,
                gene_symbol TEXT,
                cdna_notation TEXT,
                protein_notation TEXT,
                genomic_position TEXT
            )
            """
        )
        cur.execute(
            """
            CREATE TABLE variant_papers (
                variant_id INTEGER,
                pmid TEXT
            )
            """
        )
        cur.execute(
            """
            CREATE TABLE penetrance_data (
                variant_id INTEGER,
                pmid TEXT,
                total_carriers_observed INTEGER,
                affected_count INTEGER,
                unaffected_count INTEGER,
                uncertain_count INTEGER
            )
            """
        )
        cur.execute(
            """
            INSERT INTO variants
            (variant_id, gene_symbol, protein_notation, cdna_notation)
            VALUES (1, 'KCNQ1', 'p.Arg1His', 'c.2G>A')
            """
        )
        cur.execute("INSERT INTO variant_papers VALUES (1, '111')")
        cur.execute("INSERT INTO penetrance_data VALUES (1, '111', 2, 1, 1, 0)")
        conn.commit()
    finally:
        conn.close()


def test_discover_gold_inputs(tmp_path):
    gold_dir = tmp_path / "gold"
    _write_recall_csv(gold_dir / "normalized" / "KCNQ1_recall_input.csv")

    discovered = discover_gold_inputs(gold_dir)

    assert discovered == {"KCNQ1": gold_dir / "normalized" / "KCNQ1_recall_input.csv"}


def test_build_gene_results_scores_available_and_marks_missing(tmp_path):
    gold_dir = tmp_path / "gold"
    results_dir = tmp_path / "results"
    outdir = tmp_path / "metrics"
    gold_path = gold_dir / "normalized" / "KCNQ1_recall_input.csv"
    db_path = results_dir / "KCNQ1" / "run" / "KCNQ1_v2.db"
    _write_recall_csv(gold_path)
    _write_gvf_db(db_path)

    gene_results = build_gene_results(
        genes=["KCNQ1", "SCN5A", "RYR2"],
        gold_inputs={"KCNQ1": gold_path, "SCN5A": gold_path},
        db_overrides={},
        results_dir=results_dir,
        outdir=outdir,
        manifest={"genes": {"RYR2": {"note": "no per-PMID table"}}},
        variant_match_mode="fuzzy",
        fuzzy_threshold=0.80,
    )

    by_gene = {item["gene"]: item for item in gene_results}
    assert by_gene["KCNQ1"]["status"] == "scored"
    assert by_gene["SCN5A"]["status"] == "missing_db"
    assert by_gene["RYR2"]["status"] == "missing_gold"

    recall = by_gene["KCNQ1"]["summary"]["recall"]
    assert recall["variant_rows"]["matched"] == 1
    assert recall["variant_rows"]["gold"] == 2
    assert recall["patients"]["matched"] == 2
    assert recall["patients"]["gold"] == 3

    aggregate = combine_recall([by_gene["KCNQ1"]])
    assert aggregate["variant_rows"]["recall"] == 0.5
    assert aggregate["affected"]["recall"] == 0.5


def test_combine_mae_aggregates_across_genes():
    """Aggregate MAE sums |err| and matched-N across genes, recomputes ratio."""
    scored = [
        {
            "summary": {
                "mae": {
                    "carriers": {"sum_abs_error": 10, "n_matched": 5, "mae": 2.0},
                    "affected": {"sum_abs_error": 3, "n_matched": 5, "mae": 0.6},
                    "unaffected": {"sum_abs_error": 0, "n_matched": 5, "mae": 0.0},
                }
            }
        },
        {
            "summary": {
                "mae": {
                    "carriers": {"sum_abs_error": 30, "n_matched": 15, "mae": 2.0},
                    "affected": {"sum_abs_error": 9, "n_matched": 15, "mae": 0.6},
                    "unaffected": {"sum_abs_error": 5, "n_matched": 15, "mae": 0.33},
                }
            }
        },
    ]

    agg = combine_mae(scored)

    assert agg["carriers"]["sum_abs_error"] == 40
    assert agg["carriers"]["n_matched"] == 20
    assert agg["carriers"]["mae"] == 2.0  # 40/20
    assert agg["unaffected"]["sum_abs_error"] == 5
    assert agg["unaffected"]["n_matched"] == 20
    assert agg["unaffected"]["mae"] == 0.25  # 5/20


def test_combine_mae_handles_missing_mae_blocks_safely():
    """Genes without MAE blocks (older summaries) must not crash combine_mae."""
    scored = [
        {"summary": {}},  # no mae block at all
        {"summary": {"mae": {"carriers": {"sum_abs_error": 4, "n_matched": 2}}}},
    ]
    agg = combine_mae(scored)
    assert agg["carriers"]["sum_abs_error"] == 4
    assert agg["carriers"]["n_matched"] == 2
    assert agg["carriers"]["mae"] == 2.0
    # Fields absent from any gene still produce a zero/n=0 block, not a KeyError
    assert agg["affected"]["n_matched"] == 0
    assert agg["affected"]["mae"] is None
