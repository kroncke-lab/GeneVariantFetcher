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
    quarantined,
    refuse_all_unaffected,
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


# ------------------------------------------------------ arithmetic refusal


def test_refuse_clears_the_all_unaffected_contradiction():
    # The measured case: KCNQ1 33141630 T224M, 124/0/124 against gold 124/34/54,
    # which the trust gate had already quarantined as a population_count.
    assert refuse_all_unaffected(
        {"carriers": 124, "affected": 0, "unaffected": 124}, True
    ) == {"affected": None, "unaffected": None}


def test_refuse_requires_independent_corroboration():
    """The arithmetic shape alone must never be enough to delete evidence.

    A control-cohort or benign-variant paper reports exactly this shape as a
    real negative finding. The measured pair separates only on the trust tier:
    SCN5A 26746457 p.Asp1243Asn is 1/0/1, correct, and `trusted`.
    """
    counts = {"carriers": 1, "affected": 0, "unaffected": 1}
    assert refuse_all_unaffected(counts, False) == {}
    assert refuse_all_unaffected(counts, True) != {}


def test_refuse_keeps_the_carrier_total():
    """The carrier count in this pattern is usually exact; only the split is junk."""
    counts = {"carriers": 124, "affected": 0, "unaffected": 124}
    assert "carriers" not in refuse_all_unaffected(counts, True)


def test_refuse_declines_when_the_split_is_informative():
    # A genuine partial split must survive untouched even under quarantine.
    assert (
        refuse_all_unaffected({"carriers": 10, "affected": 3, "unaffected": 7}, True)
        == {}
    )
    assert (
        refuse_all_unaffected({"carriers": 10, "affected": 0, "unaffected": 4}, True)
        == {}
    )
    assert (
        refuse_all_unaffected({"carriers": 10, "affected": 10, "unaffected": 0}, True)
        == {}
    )


def test_refuse_declines_on_a_zero_carrier_row():
    assert (
        refuse_all_unaffected({"carriers": 0, "affected": 0, "unaffected": 0}, True)
        == {}
    )


def test_refuse_declines_on_nulls():
    assert (
        refuse_all_unaffected({"carriers": 5, "affected": None, "unaffected": 5}, True)
        == {}
    )


def test_quarantined_reads_the_tier_not_advisory_reasons():
    assert quarantined("quarantine", None)
    assert quarantined("QUARANTINE", "[]")
    assert not quarantined("trusted", '["population_count"]')
    assert not quarantined(None, None)


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


def test_a_quarantined_partial_split_is_left_alone():
    """Only the exact 100%-non-penetrant shape is refused, even under quarantine."""
    changes = repair_counts(
        {"carriers": 10, "affected": 4, "unaffected": 6},
        None,
        "llm_text",
        trust_tier="quarantine",
    )
    assert changes == {}


# ------------------------------------------------------------- ordering


def test_figure_counts_are_adopted_before_the_refusal_is_judged():
    """The refusal must see the figure reading, not the pre-adoption nulls."""
    notes = json.dumps(
        {"carriers": 12, "affected": 0, "unaffected": 12, "context": "pedigree"}
    )
    changes = repair_counts(
        {"carriers": None, "affected": None, "unaffected": None},
        notes,
        "figure",
        trust_tier="quarantine",
    )
    assert changes["carriers"] == (12, "adopt_figure_counts")
    assert changes["affected"] == (None, "refuse_all_unaffected")
    assert changes["unaffected"] == (None, "refuse_all_unaffected")


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
            # quarantined -> the split is refused, the carrier total survives
            {
                "carriers": 124,
                "affected": 0,
                "unaffected": 124,
                "layer": "llm_text",
                "tier": "quarantine",
                "reasons": '["population_count"]',
            },
            # same shape but trusted -> a real control-cohort finding, untouched
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
    assert summary["rows_changed"] == 2
    assert summary["adopt_figure_counts"] == 2
    assert summary["refuse_all_unaffected"] == 2
    assert _counts(db) == [(13, 5, None), (124, None, None), (1, 0, 1)]


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
                "carriers": 9,
                "affected": 0,
                "unaffected": 9,
                "layer": "llm_text",
                "tier": "quarantine",
            }
        ],
    )
    before = _counts(db)
    summary = apply_count_repair(db, dry_run=True)
    assert summary["refuse_all_unaffected"] == 2
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
                "carriers": 124,
                "affected": 0,
                "unaffected": 124,
                "layer": "llm_text",
                "tier": "quarantine",
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
        ("affected", 0, None, "refuse_all_unaffected"),
        ("unaffected", 124, None, "refuse_all_unaffected"),
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
