"""Per-file SAVEPOINT atomicity for directory migration.

Regression guard: a file that fails partway through migration (e.g. a
constraint error raised *after* its paper metadata was already inserted) must
leave no rows behind, and must not stop later files from committing. Before the
savepoint fix, the partial rows rode along on the directory-wide commit.
"""

import json
import sqlite3

import harvesting.migrate_to_sqlite as m
from harvesting.migrate_to_sqlite import (
    create_database_schema,
    migrate_extraction_directory,
)


def _extraction(pmid: str, protein: str = "p.Val254Met") -> dict:
    return {
        "paper_metadata": {"pmid": pmid, "title": "t", "journal": "j"},
        "variants": [
            {
                "gene_symbol": "KCNQ1",
                "protein_notation": protein,
                "penetrance_data": {"total_carriers_observed": 25},
                "patients": {"count": 25, "phenotype": "LQTS"},
            }
        ],
        "extraction_metadata": {"total_variants_found": 1},
    }


def _write(tmp_path, pmid: str) -> None:
    (tmp_path / f"KCNQ1_PMID_{pmid}.json").write_text(json.dumps(_extraction(pmid)))


def test_failed_file_leaves_no_partial_rows(tmp_path, monkeypatch):
    # Bad file sorts first so we also prove a later good file still commits
    # after the failed file's savepoint is rolled back and released.
    _write(tmp_path, "1000001")  # will fail mid-insert (after paper metadata)
    _write(tmp_path, "1000002")  # good

    conn = create_database_schema(str(tmp_path / "t.db"))

    real_insert_variant = m.insert_variant_data

    def flaky(cursor, pmid, variant, **kwargs):
        if pmid == "1000001":
            raise sqlite3.IntegrityError("simulated mid-file failure")
        return real_insert_variant(cursor, pmid, variant, **kwargs)

    monkeypatch.setattr(m, "insert_variant_data", flaky)

    stats = migrate_extraction_directory(conn, tmp_path)

    assert stats["successful"] == 1
    assert stats["failed"] == 1

    pmids = {row[0] for row in conn.execute("SELECT pmid FROM papers")}
    assert "1000002" in pmids, "later good file must still commit"
    assert "1000001" not in pmids, "failed file's paper row must be rolled back"
    orphan_penetrance = conn.execute(
        "SELECT COUNT(*) FROM penetrance_data WHERE pmid = '1000001'"
    ).fetchone()[0]
    assert orphan_penetrance == 0, "failed file must leave no child rows"
    conn.close()


def test_all_good_files_commit(tmp_path):
    _write(tmp_path, "2000001")
    _write(tmp_path, "2000002")

    conn = create_database_schema(str(tmp_path / "t.db"))
    stats = migrate_extraction_directory(conn, tmp_path)

    assert stats["successful"] == 2
    assert stats["failed"] == 0
    pmids = {row[0] for row in conn.execute("SELECT pmid FROM papers")}
    assert {"2000001", "2000002"} <= pmids
    conn.close()
