import sqlite3
import types

from pipeline.carrier_guard import apply_carrier_guard
from pipeline.full_coverage import run_walk_to_taper


def test_carrier_guard_quarantines_and_nulls_implausible_counts(tmp_path):
    db = tmp_path / "gene.db"
    con = sqlite3.connect(db)
    con.execute(
        """CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER,
            pmid TEXT,
            total_carriers_observed INTEGER,
            affected_count INTEGER,
            unaffected_count INTEGER,
            uncertain_count INTEGER,
            penetrance_percentage REAL
        )"""
    )
    con.execute(
        """INSERT INTO penetrance_data (
            variant_id, pmid, total_carriers_observed, affected_count,
            unaffected_count, uncertain_count, penetrance_percentage
        ) VALUES (1, '123', 565000, 10, 2, 1, 75.0)"""
    )
    con.commit()
    con.close()

    result = apply_carrier_guard(db, threshold=100_000)

    assert result["guarded"] == 1
    con = sqlite3.connect(db)
    row = con.execute(
        """SELECT total_carriers_observed, affected_count, unaffected_count,
                  uncertain_count, penetrance_percentage
           FROM penetrance_data"""
    ).fetchone()
    assert row == (None, None, None, None, None)
    audit = con.execute(
        """SELECT total_carriers_observed, threshold, reason
           FROM carrier_guard_quarantine"""
    ).fetchone()
    assert audit[0] == 565000
    assert audit[1] == 100_000
    assert "plausible per-variant ceiling" in audit[2]
    con.close()


def test_walk_to_taper_passes_offset_and_workers(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    priority_dir = run_dir / "extraction_priority"
    priority_dir.mkdir(parents=True)
    (priority_dir / "priority_candidates.tsv").write_text(
        "pmid\tscore\n1\t1\n2\t1\n3\t1\n",
        encoding="utf-8",
    )
    db = run_dir / "GENE.db"
    db.write_text("")
    commands = []

    def fake_run(cmd, **kwargs):
        commands.append((cmd, kwargs))
        return types.SimpleNamespace(returncode=0)

    monkeypatch.setattr("pipeline.full_coverage._distinct_variants", lambda db: 10)
    monkeypatch.setattr("pipeline.full_coverage.subprocess.run", fake_run)

    result = run_walk_to_taper(
        "GENE",
        run_dir,
        model="azure_ai/gpt-5.4",
        max_workers=7,
        step=1,
        start_offset=2,
        min_new_variants=99,
        patience=1,
    )

    assert result["walked"] is True
    assert result["start_offset"] == 2
    assert len(commands) == 1
    cmd = commands[0][0]
    assert "--priority-offset" in cmd
    assert cmd[cmd.index("--priority-offset") + 1] == "2"
    assert "--max-workers" in cmd
    assert cmd[cmd.index("--max-workers") + 1] == "7"
    assert (priority_dir / "full_coverage_walk_summary.json").exists()
