"""Coverage for the deterministic count-repair pass.

Each rule was adopted on measured evidence from the locked Tier-1 baseline, so
the tests pin the behaviour that evidence justifies — including the cases where
a rule must decline to act, which is what keeps it from manufacturing counts.
"""

from __future__ import annotations

import json
import sqlite3

import pytest

from pipeline.count_repair import (
    adopt_figure_counts,
    apply_count_repair,
    figure_context_is_patient_level,
    figure_counts,
    repair_counts,
)


# ----------------------------------------------------------- figure counts


def test_figure_counts_reads_the_serialized_blob():
    notes = json.dumps(
        {
            "carriers": 13,
            "affected": None,
            "unaffected": None,
            "context": "Table 2 carrier counts per variant",
        }
    )
    assert figure_counts(notes) == {"carriers": 13}


def test_a_survival_curve_is_not_a_patient_count():
    """KCNQ1 18713323: a Kaplan-Meier at-risk table yielded 13/8/25/19.

    Those are people remaining at a timepoint, not carriers of the variant.
    """
    notes = json.dumps(
        {"carriers": 13, "context": "Kaplan-Meier plot and patients-at-risk table"}
    )
    assert figure_counts(notes) == {}


@pytest.mark.parametrize(
    "context",
    [
        "Kaplan-Meier survival",
        "patients at risk",
        "cumulative event-free survival",
        "y-axis label",
        "hazard ratio plot",
    ],
)
def test_curve_contexts_are_rejected(context):
    assert not figure_context_is_patient_level(context)


def test_a_missing_context_fails_closed():
    """No description is not evidence of a patient count."""
    assert not figure_context_is_patient_level(None)
    assert not figure_context_is_patient_level("   ")
    assert figure_counts(json.dumps({"carriers": 13})) == {}


def test_figure_counts_tolerates_junk():
    for notes in (None, "", "not json", "[1,2]", json.dumps("scalar")):
        assert figure_counts(notes) == {}


def test_figure_counts_rejects_non_integers():
    notes = json.dumps(
        {"carriers": "13", "affected": 2.5, "unaffected": True, "context": "pedigree"}
    )
    assert figure_counts(notes) == {}


def test_figure_counts_rejects_negatives():
    assert figure_counts(json.dumps({"carriers": -1, "context": "pedigree"})) == {}


def test_adopt_only_applies_to_figure_rows():
    notes = json.dumps({"carriers": 13, "context": "pedigree"})
    assert adopt_figure_counts({"carriers": None}, notes, "llm_text") == {}
    assert adopt_figure_counts({"carriers": None}, notes, "figure") == {"carriers": 13}


def test_adopt_never_overwrites_a_stored_count():
    notes = json.dumps({"carriers": 13, "context": "pedigree"})
    assert adopt_figure_counts({"carriers": 8}, notes, "figure") == {}


def test_adopt_fills_zero_which_is_a_real_count():
    notes = json.dumps({"carriers": 5, "affected": 0, "context": "pedigree"})
    got = adopt_figure_counts({"carriers": None, "affected": None}, notes, "figure")
    assert got == {"carriers": 5, "affected": 0}


# --------------------------------- the rule that was measured and REJECTED


def test_no_arithmetic_derivation_of_unaffected():
    """`unaffected = carriers - affected` must not exist.

    It scored 59/61 on cardiac gold, but every one of those 61 derivations wrote
    0 -- all were carriers == affected rows, 45 with a single carrier, none with
    a cohort of 10+. The arithmetic was never exercised; the rule was a "default
    unaffected to 0", which pipeline/prompts.py forbids verbatim. On a biobank
    cohort it would reclassify never-phenotyped people as confirmed unaffected,
    which pipeline/claim_verifier.py refuses by name.
    """
    import pipeline.count_repair as mod

    assert not hasattr(mod, "derive_unaffected")
    assert "derive_unaffected" not in mod.RULES
    changes = repair_counts(
        {"carriers": 10, "affected": 4, "unaffected": None}, None, "llm_table"
    )
    assert changes == {}


def test_rule_selection_happens_before_evaluation():
    notes = json.dumps(
        {"carriers": 12, "affected": 0, "unaffected": 12, "context": "pedigree"}
    )
    counts = {"carriers": None, "affected": None, "unaffected": None}

    adopted = repair_counts(
        counts,
        notes,
        "figure",
        rules={"adopt_figure_counts"},
    )
    assert adopted == {
        "carriers": (12, "adopt_figure_counts"),
        "affected": (0, "adopt_figure_counts"),
        "unaffected": (12, "adopt_figure_counts"),
    }
    assert repair_counts(counts, notes, "figure", rules=set()) == {}


def test_repair_reports_no_change_when_nothing_applies():
    assert (
        repair_counts({"carriers": 5, "affected": 2, "unaffected": 3}, None, "x") == {}
    )


# ------------------------------------------------------------ database


def _db(tmp_path, rows):
    path = tmp_path / "run.db"
    con = sqlite3.connect(path)
    con.executescript(
        """CREATE TABLE penetrance_data(
             penetrance_id INTEGER PRIMARY KEY, variant_id INTEGER, pmid TEXT,
             total_carriers_observed INTEGER, affected_count INTEGER,
             unaffected_count INTEGER, trust_tier TEXT, trust_reasons TEXT);
           CREATE TABLE variant_papers(
             variant_id INTEGER, pmid TEXT, additional_notes TEXT,
             source_layer TEXT);"""
    )
    for i, r in enumerate(rows, start=1):
        con.execute(
            "INSERT INTO penetrance_data VALUES(?,?,?,?,?,?,?,?)",
            (
                i,
                i,
                r.get("pmid", "1"),
                r.get("carriers"),
                r.get("affected"),
                r.get("unaffected"),
                r.get("tier"),
                r.get("reasons"),
            ),
        )
        con.execute(
            "INSERT INTO variant_papers VALUES(?,?,?,?)",
            (i, r.get("pmid", "1"), r.get("notes"), r.get("layer")),
        )
    con.commit()
    con.close()
    return path


def _counts(path):
    con = sqlite3.connect(path)
    out = list(
        con.execute(
            "SELECT total_carriers_observed, affected_count, unaffected_count "
            "FROM penetrance_data ORDER BY penetrance_id"
        )
    )
    con.close()
    return out


def test_apply_repairs_a_database(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "carriers": None,
                "layer": "figure",
                "notes": json.dumps(
                    {"carriers": 13, "affected": 5, "context": "pedigree"}
                ),
            },
            # Existing raw observations are never destructively cleared.
            {
                "carriers": 124,
                "affected": 0,
                "unaffected": 124,
                "layer": "llm_text",
                "tier": "quarantine",
                "reasons": '["population_count"]',
            },
            # A real control-cohort finding is likewise untouched.
            {
                "carriers": 1,
                "affected": 0,
                "unaffected": 1,
                "layer": "llm_text",
                "tier": "trusted",
                "reasons": "[]",
            },
        ],
    )
    summary = apply_count_repair(db)
    assert summary["rows_changed"] == 1
    assert summary["adopt_figure_counts"] == 2
    assert _counts(db) == [(13, 5, None), (124, 0, 124), (1, 0, 1)]


def test_apply_is_idempotent(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "carriers": 124,
                "affected": 0,
                "unaffected": 124,
                "layer": "llm_text",
                "tier": "quarantine",
            },
            {
                "carriers": None,
                "layer": "figure",
                "notes": json.dumps(
                    {"carriers": 4, "affected": 1, "context": "pedigree"}
                ),
            },
        ],
    )
    apply_count_repair(db)
    first = _counts(db)
    second_summary = apply_count_repair(db)
    assert second_summary["rows_changed"] == 0
    assert _counts(db) == first


def test_dry_run_changes_nothing(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "layer": "figure",
                "notes": json.dumps(
                    {"carriers": 9, "affected": 0, "context": "pedigree"}
                ),
            }
        ],
    )
    before = _counts(db)
    summary = apply_count_repair(db, dry_run=True)
    assert summary["adopt_figure_counts"] == 2
    assert summary["dry_run"] is True
    assert _counts(db) == before
    con = sqlite3.connect(db)
    tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master")}
    con.close()
    assert "count_repair_log" not in tables


def test_every_change_is_logged_reversibly(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "layer": "figure",
                "notes": json.dumps(
                    {"carriers": 13, "affected": 5, "context": "pedigree"}
                ),
            }
        ],
    )
    apply_count_repair(db)
    con = sqlite3.connect(db)
    log = list(
        con.execute(
            "SELECT field, previous_value, new_value, rule FROM count_repair_log "
            "ORDER BY field"
        )
    )
    con.close()
    assert log == [
        ("affected", None, 5, "adopt_figure_counts"),
        ("carriers", None, 13, "adopt_figure_counts"),
    ]


def test_rules_can_be_selected(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "carriers": 124,
                "affected": 0,
                "unaffected": 124,
                "layer": "llm_text",
                "tier": "quarantine",
            }
        ],
    )
    summary = apply_count_repair(db, rules={"adopt_figure_counts"})
    assert summary["rows_changed"] == 0
    assert _counts(db) == [(124, 0, 124)]


def test_missing_penetrance_table_is_not_fatal(tmp_path):
    path = tmp_path / "empty.db"
    sqlite3.connect(path).close()
    assert apply_count_repair(path)["rows_changed"] == 0


def test_a_paper_row_with_no_penetrance_row_gets_one_created(tmp_path):
    """The whole point of the figure rule: there is nothing to UPDATE.

    migrate_to_sqlite only writes a penetrance_data row once some count is
    non-null, so every figure row whose counts sit unread in additional_notes
    has no penetrance row at all. Driving from penetrance_data finds zero of
    them — on the four production databases that mistake made the rule a no-op.
    """
    path = tmp_path / "run.db"
    con = sqlite3.connect(path)
    con.executescript(
        """CREATE TABLE penetrance_data(
             penetrance_id INTEGER PRIMARY KEY, variant_id INTEGER, pmid TEXT,
             total_carriers_observed INTEGER, affected_count INTEGER,
             unaffected_count INTEGER, trust_tier TEXT, trust_reasons TEXT);
           CREATE TABLE variant_papers(
             variant_id INTEGER, pmid TEXT, additional_notes TEXT,
             source_layer TEXT);"""
    )
    con.execute(
        "INSERT INTO variant_papers VALUES(?,?,?,?)",
        (
            1,
            "1",
            json.dumps({"carriers": 13, "affected": 5, "context": "pedigree"}),
            "figure",
        ),
    )
    con.commit()
    con.close()

    summary = apply_count_repair(path)
    assert summary["adopt_figure_counts"] == 2
    assert summary["rows_created"] == 1
    # unaffected stays NULL: no lane asserted it, and nothing derives it.
    assert _counts(path) == [(13, 5, None)]


def test_dry_run_counts_the_rows_it_would_create(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "layer": "figure",
                "notes": json.dumps(
                    {"carriers": 4, "affected": 1, "context": "pedigree"}
                ),
            }
        ],
    )
    # _db always writes a penetrance row, so drop it to model the real shape.
    con = sqlite3.connect(db)
    con.execute("DELETE FROM penetrance_data")
    con.commit()
    con.close()

    summary = apply_count_repair(db, dry_run=True)
    assert summary["rows_created"] == 1
    assert _counts(db) == []


@pytest.mark.parametrize("layer", ["figure", "figure,llm_table", "mixed,figure"])
def test_composite_layer_strings_still_adopt(layer):
    notes = json.dumps({"carriers": 7, "context": "pedigree"})
    assert adopt_figure_counts({"carriers": None}, notes, layer) == {"carriers": 7}


def test_ambiguous_multi_cohort_parent_fails_closed(tmp_path):
    path = tmp_path / "ambiguous.db"
    con = sqlite3.connect(path)
    con.executescript(
        """CREATE TABLE penetrance_data(
             penetrance_id INTEGER PRIMARY KEY, variant_id INTEGER, pmid TEXT,
             total_carriers_observed INTEGER, affected_count INTEGER,
             unaffected_count INTEGER, trust_tier TEXT, trust_reasons TEXT);
           CREATE TABLE variant_papers(
             variant_id INTEGER, pmid TEXT, additional_notes TEXT,
             source_layer TEXT);
           INSERT INTO penetrance_data VALUES(1, 7, '111', NULL, NULL, NULL, 'trusted', '[]');
           INSERT INTO penetrance_data VALUES(2, 7, '111', NULL, NULL, NULL, 'trusted', '[]');
           INSERT INTO variant_papers VALUES(
             7, '111', '{"carriers": 13, "affected": 5, "context": "pedigree"}',
             'figure'
           );"""
    )
    con.commit()
    con.close()

    summary = apply_count_repair(path)

    assert summary["ambiguous_variant_papers"] == 1
    assert summary["rows_changed"] == 0
    assert _counts(path) == [(None, None, None), (None, None, None)]


def test_adopted_figure_fields_land_pending_trust_review(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "carriers": None,
                "layer": "figure",
                "notes": json.dumps(
                    {"carriers": 13, "affected": 5, "context": "pedigree"}
                ),
                "tier": "trusted",
                "reasons": "[]",
            }
        ],
    )
    con = sqlite3.connect(db)
    con.execute("DELETE FROM penetrance_data")
    con.commit()
    con.close()

    apply_count_repair(db)

    con = sqlite3.connect(db)
    tier, reasons = con.execute(
        "SELECT trust_tier, trust_reasons FROM penetrance_data"
    ).fetchone()
    con.close()
    assert tier == "quarantine"
    assert "figure_count_pending_review" in json.loads(reasons)


def test_pending_figure_carriers_use_the_projection_field_name(tmp_path):
    db = _db(
        tmp_path,
        [
            {
                "carriers": None,
                "affected": 5,
                "layer": "figure,llm_text",
                "notes": json.dumps({"carriers": 13, "context": "pedigree"}),
                "tier": "trusted",
                "reasons": "[]",
            }
        ],
    )
    con = sqlite3.connect(db)
    con.execute("ALTER TABLE penetrance_data ADD COLUMN field_trust TEXT")
    con.execute("ALTER TABLE penetrance_data ADD COLUMN trust_sources TEXT")
    con.execute("ALTER TABLE penetrance_data ADD COLUMN trust_rule_version TEXT")
    con.execute(
        "UPDATE penetrance_data SET field_trust=?, trust_sources='[]'",
        (json.dumps({"affected": "trusted"}),),
    )
    con.commit()
    con.close()

    apply_count_repair(db)

    con = sqlite3.connect(db)
    tier, field_trust_raw = con.execute(
        "SELECT trust_tier, field_trust FROM penetrance_data"
    ).fetchone()
    con.close()
    field_trust = json.loads(field_trust_raw)
    assert tier == "quarantine"
    assert field_trust["total_carriers"] == "quarantine"
    assert field_trust["affected"] == "trusted"
    assert "carriers" not in field_trust
