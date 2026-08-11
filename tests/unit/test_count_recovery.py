"""Offline tests for pipeline.count_recovery. No network, no model calls."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest

from pipeline.count_recovery import (
    COUNT_RECOVERY_XHIGH_MAX_TOKENS,
    PLAUSIBLE_COUNT_CEILING,
    CountRecoveryError,
    PaperGap,
    VariantGap,
    build_prompt,
    find_count_gaps,
    parse_response,
    quote_supports,
    recover_counts,
    validate_paper_response,
    write_recovered_counts,
)

SOURCE = (
    "Methods. We screened 812 probands for KCNH2 variants.\n"
    "Results. The p.Leu552Ser variant was identified in 44 carriers, of whom 12 "
    "had documented arrhythmic events and 25 remained asymptomatic. "
    "The p.Arg752Trp variant was seen once in this cohort. "
    "The p.Gly628Ser variant is listed in Table 2 without carrier information. "
    "gnomAD reports an allele number of 251384 at this locus.\n"
)


def _make_db(tmp_path, rows):
    """rows: (variant_id, protein, cdna, pmid, carriers, affected, unaffected)."""
    db = tmp_path / "T.db"
    con = sqlite3.connect(db)
    con.executescript(
        """
        CREATE TABLE variants(variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT);
        CREATE TABLE variant_papers(variant_id INTEGER, pmid TEXT,
            source_location TEXT, key_quotes TEXT, source_layer TEXT);
        CREATE TABLE penetrance_data(penetrance_id INTEGER PRIMARY KEY,
            variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
            affected_count INTEGER, unaffected_count INTEGER, uncertain_count INTEGER);
        """
    )
    for vid, prot, cdna, pmid, car, aff, una in rows:
        con.execute("INSERT INTO variants VALUES (?,?,?,?)", (vid, "KCNH2", cdna, prot))
        con.execute(
            "INSERT INTO variant_papers VALUES (?,?,?,?,?)",
            (vid, pmid, "Table 2", "", "llm_table"),
        )
        if (car, aff, una) != (None, None, None) or True:
            con.execute(
                "INSERT INTO penetrance_data (variant_id, pmid, "
                "total_carriers_observed, affected_count, unaffected_count) "
                "VALUES (?,?,?,?,?)",
                (vid, pmid, car, aff, una),
            )
    con.commit()
    con.close()
    return db


def _stub(payload_by_call):
    """Return an llm_caller yielding successive canned response strings."""
    calls = {"n": 0, "prompts": [], "kwargs": []}

    def caller(**kwargs):
        calls["prompts"].append(kwargs["messages"][0]["content"])
        calls["kwargs"].append(kwargs)
        i = calls["n"]
        calls["n"] += 1
        text = payload_by_call[min(i, len(payload_by_call) - 1)]
        return SimpleNamespace(
            choices=[SimpleNamespace(message=SimpleNamespace(content=text))]
        )

    caller.calls = calls
    return caller


class TestFindCountGaps:
    def test_reports_only_null_fields(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", "c.1655T>C", "111", None, None, None),
                (2, "p.Arg752Trp", None, "111", 26, None, None),
            ],
        )
        gaps = find_count_gaps(db, "KCNH2", fields=("carriers",))
        assert len(gaps) == 1
        assert [v.variant_id for v in gaps[0].variants] == [1]

    def test_variant_with_no_penetrance_row_is_a_gap(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        sqlite3.connect(db).execute("DELETE FROM penetrance_data").connection.commit()
        gaps = find_count_gaps(db, "KCNH2")
        assert gaps and gaps[0].variants[0].variant_id == 1

    def test_field_populated_on_any_row_is_not_a_gap(self, tmp_path):
        """A variant reached via two layers, counted on one of them, is complete."""
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        con = sqlite3.connect(db)
        con.execute(
            "INSERT INTO variant_papers VALUES (1,'111','ClinVar','','clinvar')"
        )
        con.execute(
            "INSERT INTO penetrance_data (variant_id,pmid,total_carriers_observed) "
            "VALUES (1,'111',44)"
        )
        con.commit()
        con.close()
        assert find_count_gaps(db, "KCNH2", fields=("carriers",)) == []

    def test_pmid_filter_and_bad_field(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        assert find_count_gaps(db, "KCNH2", pmids=["999"]) == []
        with pytest.raises(ValueError):
            find_count_gaps(db, "KCNH2", fields=("penetrance",))


class TestQuoteGrounding:
    @pytest.mark.parametrize(
        "quote,value",
        [
            ("identified in 44 carriers", 44),
            ("seen once in this cohort", 1),
            ("both siblings carried it", 2),
            ("no carriers were found", 0),
            ("no one carried it", 0),
            ("1,234 patients", 1234),
        ],
    )
    def test_supported(self, quote, value):
        assert quote_supports(quote, value)

    @pytest.mark.parametrize(
        "quote,value",
        [
            ("identified in 44 carriers", 45),
            ("listed in Table 2", 3),
            ("p.Arg752Trp was listed", 752),
            ("none of the subjects carried it", 1),
            ("the notation was reviewed", 0),
            ("eighteen carriers", 8),
            ("one hundred carriers", 1),
            ("", 1),
        ],
    )
    def test_unsupported(self, quote, value):
        assert not quote_supports(quote, value)


class TestValidation:
    def _gap(self):
        return PaperGap(
            "KCNH2",
            "111",
            [
                VariantGap(1, "p.Leu552Ser", ["carriers"]),
                VariantGap(2, "p.Arg752Trp", ["carriers"]),
            ],
        )

    def test_accepts_grounded_count(self):
        res = validate_paper_response(
            self._gap(),
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers",
                }
            ],
            SOURCE,
        )
        assert [(c.variant_id, c.field, c.value) for c in res.accepted] == [
            (1, "carriers", 44)
        ]
        assert res.rejected == []

    def test_accepts_word_number(self):
        """'once' must ground carriers=1 -- strict digit-only grounding drops it."""
        res = validate_paper_response(
            self._gap(),
            [
                {
                    "variant": "p.Arg752Trp",
                    "carriers": 1,
                    "quote": "The p.Arg752Trp variant was seen once in this cohort",
                }
            ],
            SOURCE,
        )
        assert [c.value for c in res.accepted] == [1]

    def test_rejects_quote_not_in_source(self):
        res = validate_paper_response(
            self._gap(),
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "quote": "we observed 44 carriers of this variant",
                }
            ],
            SOURCE,
        )
        assert not res.accepted
        assert "verbatim" in res.rejected[0]["reason"]

    def test_rejects_number_absent_from_quote(self):
        res = validate_paper_response(
            self._gap(),
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 7,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers",
                }
            ],
            SOURCE,
        )
        assert not res.accepted
        assert "number not stated" in res.rejected[0]["reason"]

    def test_rejects_study_wide_count_without_variant_binding(self):
        res = validate_paper_response(
            self._gap(),
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 812,
                    "quote": "Methods. We screened 812 probands for KCNH2 variants.",
                }
            ],
            SOURCE,
        )
        assert not res.accepted
        assert "bind" in res.rejected[0]["reason"]

    def test_rejects_collapsed_multi_variant_table_block(self):
        quote = "453delC* R176W*\n1631delAG* Y569H* G584S*\n1\n4\n1\n1\n10\n8"
        gap = PaperGap("KCNH2", "111", [VariantGap(3, "p.Gly584Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Gly584Ser", "carriers": 1, "quote": quote}],
            quote,
        )
        assert not res.accepted
        assert "bind" in res.rejected[0]["reason"]

    def test_rejects_unlabeled_table_row_with_multiple_integer_values(self):
        quote = "| p.L552S | 73 (31) | 27 |"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 73, "quote": quote}],
            quote,
        )
        assert not res.accepted
        assert "table row" in res.rejected[0]["reason"]

    def test_accepts_single_unambiguous_table_value_with_equivalent_notation(self):
        quote = "| p.L552S | 73 |"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 73, "quote": quote}],
            quote,
        )
        assert [c.value for c in res.accepted] == [73]

    def test_accepts_table_value_mapped_by_header(self):
        quote = "| Variant | Carriers | Affected |\n| p.L552S | 73 | 31 |"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 73, "quote": quote}],
            quote,
        )
        assert [c.value for c in res.accepted] == [73]

    def test_accepts_target_local_count_in_multi_variant_table_row(self):
        quote = (
            "| Mutations by tertile | Y184S | 14 | 0.21 | "
            "W120C | 2 | 0.62 | E160K | 3 | 0.80 |"
        )
        gap = PaperGap("KCNQ1", "111", [VariantGap(1, "Y184S", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "Y184S", "carriers": 14, "quote": quote}],
            quote,
        )
        assert [c.value for c in res.accepted] == [14]

    def test_accepts_uniquely_bound_multi_line_table_group(self):
        quote = (
            "| Genotype-positive G601S families |  |\n"
            "| n | 3 |\n"
            "| Genotype-positive subjects, n (female) | 9 (6) |"
        )
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Gly601Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Gly601Ser", "carriers": 9, "quote": quote}],
            quote,
        )
        assert [c.value for c in res.accepted] == [9]

    def test_rejects_count_from_detached_study_total_table_row(self):
        quote = "| p.Leu552Ser | pathogenic |\n| Study total | 73 |"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 73, "quote": quote}],
            quote,
        )
        assert not res.accepted
        assert "table row" in res.rejected[0]["reason"]

    def test_rejects_different_cdna_change_in_target_protein_codon(self):
        quote = "The c.1202A>T variant was found in 3 carriers."
        gap = PaperGap(
            "KCNH2",
            "111",
            [
                VariantGap(
                    1,
                    "p.Gly401Ser c.1201G>A",
                    ["carriers"],
                    parts=("p.Gly401Ser", "c.1201G>A"),
                )
            ],
        )
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Gly401Ser c.1201G>A",
                    "carriers": 3,
                    "quote": quote,
                }
            ],
            quote,
        )
        assert not res.accepted
        assert "bind" in res.rejected[0]["reason"]

    def test_rejects_position_only_cdna_to_protein_bridge(self):
        quote = "The c.1201G>A variant was found in 3 carriers."
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Gly401Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Gly401Ser", "carriers": 3, "quote": quote}],
            quote,
        )
        assert not res.accepted

    def test_accepts_collective_zero_when_exact_target_is_one_of_named_variants(self):
        quote = (
            "All study individuals were screened for the SCN5A variants "
            "(1062T>C, 354G>C, 287C>T and 1199G>C). However, none of the "
            "patients or controls were positive for these variants."
        )
        gap = PaperGap("SCN5A", "111", [VariantGap(1, "c.1199G>C", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "c.1199G>C", "carriers": 0, "quote": quote}],
            quote,
        )
        assert [c.value for c in res.accepted] == [0]

    def test_rejects_quote_that_binds_target_to_another_gene(self):
        quote = (
            "One participant had both rs199473599 (SCN5A-D1243N) and "
            "rs149955375 (KCNH2-T983I)."
        )
        gap = PaperGap("SCN5A", "111", [VariantGap(1, "T983I", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "T983I", "carriers": 1, "quote": quote}],
            quote,
        )
        assert not res.accepted
        assert "bind" in res.rejected[0]["reason"]

    def test_rejects_single_variant_explicitly_prefixed_by_another_gene(self):
        quote = "The KCNH2-T983I variant was found in 2 carriers."
        gap = PaperGap("SCN5A", "111", [VariantGap(1, "T983I", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "T983I", "carriers": 2, "quote": quote}],
            quote,
        )
        assert not res.accepted
        assert "bind" in res.rejected[0]["reason"]

    def test_rejects_ambiguous_structural_variant_table_row(self):
        quote = (
            r"| SCN5A | p.(Phe1617\_del) | c.4850\_4852del | "
            r"in-frame deletion | 2 | 6 | 3 |"
        )
        gap = PaperGap("SCN5A", "111", [VariantGap(1, "P.PHE1617DEL", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "P.PHE1617DEL", "carriers": 6, "quote": quote}],
            quote,
        )
        assert not res.accepted

    def test_accepts_structural_variant_with_carriers_header(self):
        quote = (
            r"| Gene | Variant | cDNA | Carriers |"
            "\n"
            r"| SCN5A | p.(Phe1617\_del) | c.4850\_4852del | 6 |"
        )
        gap = PaperGap("SCN5A", "111", [VariantGap(1, "P.PHE1617DEL", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "P.PHE1617DEL", "carriers": 6, "quote": quote}],
            quote,
        )
        assert [c.value for c in res.accepted] == [6]

    def test_rejects_family_count_as_carriers(self):
        quote = "The p.Leu552Ser variant was found in 12 families."
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 12, "quote": quote}],
            quote,
        )
        assert not res.accepted

    # ------------------------------------------------------------------
    # Role-gate regression fixtures: the exact quotes the 2026-07-28 audit
    # reproduced as ACCEPTED carrier counts. Prose denominators/cohort totals
    # (P1, both reviews) and explicitly non-carrier table columns (P0-1, the
    # PMID 33013630 annotation-table failure class) must all fail closed.
    # ------------------------------------------------------------------
    ROLE_GATE_PROBES = [
        # (quote, every integer the prior gate accepted as `carriers`)
        ("The p.Leu552Ser variant was found in 12 of 812 subjects.", 12),
        ("The p.Leu552Ser variant was found in 12 of 812 subjects.", 812),
        ("p.Leu552Ser was present in 5/120 patients", 5),
        ("p.Leu552Ser was present in 5/120 patients", 120),
        (
            "Among 44 cases of Long QT, the p.Leu552Ser variant was identified once.",
            44,
        ),
        (
            "Among 44 cases of Long QT, the p.Leu552Ser variant was identified once.",
            1,
        ),
        ("The p.Leu552Ser variant was identified in 3 of 200 probands (1.5%).", 3),
        ("The p.Leu552Ser variant was identified in 3 of 200 probands (1.5%).", 200),
        ("n=7 among 913 individuals carried p.Leu552Ser", 7),
        ("n=7 among 913 individuals carried p.Leu552Ser", 913),
        ("| p.Leu552Ser | allele count | 12 |", 12),
        ("| p.Leu552Ser | gnomAD allele count 37 |", 37),
        ("| Variant | gnomAD AC |\n| p.Leu552Ser | 12 |", 12),
        ("| Variant | Age at diagnosis |\n| p.Leu552Ser | 45 |", 45),
        ("| p.Leu552Ser | families screened | 3 |", 3),
    ]

    @pytest.mark.parametrize("quote,value", ROLE_GATE_PROBES)
    def test_rejects_denominator_and_non_carrier_roles(self, quote, value):
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": value, "quote": quote}],
            quote,
        )
        assert not res.accepted, f"{value} accepted from {quote!r}"
        assert res.rejected and res.rejected[0]["value"] == value

    def test_a_cited_table_number_does_not_block_the_real_count(self):
        """Fail-closed must not become fail-always: 'Table 2' is not a unit."""
        quote = "As shown in Table 2, the p.Leu552Ser variant was found in 12 carriers."
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        accepted = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 12, "quote": quote}],
            quote,
        )
        assert [c.value for c in accepted.accepted] == [12]
        # ...but the table's own number is still not a carrier count.
        refused = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 2, "quote": quote}],
            quote,
        )
        assert not refused.accepted

    def test_measurement_beside_a_number_still_fails_closed(self):
        quote = "The p.Leu552Ser carrier had a QTc of 480 ms at presentation."
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [{"variant": "p.Leu552Ser", "carriers": 480, "quote": quote}],
            quote,
        )
        assert not res.accepted

    def test_model_declared_non_per_variant_role_is_decisive(self):
        """A clean carrier quote is still refused when the model admits the role."""
        quote = "The p.Leu552Ser variant was identified in 44 carriers"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "count_role": "cohort_total",
                    "quote": quote,
                }
            ],
            quote,
        )
        assert not res.accepted
        assert "non-per-variant role" in res.rejected[0]["reason"]

    def test_mismatched_declared_role_is_refused(self):
        quote = "The p.Leu552Ser variant was identified in 44 carriers"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "count_role": "per_variant_affected",
                    "quote": quote,
                }
            ],
            quote,
        )
        assert not res.accepted
        assert "does not match the requested field" in res.rejected[0]["reason"]

    def test_accepted_count_carries_role_and_locator(self):
        quote = "| Variant | Carriers |\n| p.L552S | 73 |"
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 73,
                    "count_role": "per_variant_carriers",
                    "evidence_locator": {"kind": "table_row", "column": "Carriers"},
                    "quote": quote,
                }
            ],
            quote,
        )
        (accepted,) = res.accepted
        assert accepted.count_role == "per_variant_carriers"
        assert accepted.evidence_locator["kind"] == "table_header_column"
        assert accepted.evidence_locator["header"] == "Carriers"
        assert accepted.evidence_locator["model_declared"]["column"] == "Carriers"
        assert accepted.model_declared_role == "per_variant_carriers"

    def test_rejects_variant_not_asked_about(self):
        res = validate_paper_response(
            self._gap(),
            [{"variant": "p.Gly628Ser", "carriers": 3, "quote": SOURCE[:40]}],
            SOURCE,
        )
        assert not res.accepted
        assert res.rejected[0]["reason"] == "not a variant we asked about"

    def test_rejects_population_scale_value(self):
        res = validate_paper_response(
            self._gap(),
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 251384,
                    "quote": "gnomAD reports an allele number of 251384 at this locus",
                }
            ],
            SOURCE,
        )
        assert not res.accepted
        assert str(PLAUSIBLE_COUNT_CEILING) in res.rejected[0]["reason"]

    @pytest.mark.parametrize("bad", [-1, True, 4.5, "44"])
    def test_rejects_non_int(self, bad):
        res = validate_paper_response(
            self._gap(),
            [{"variant": "p.Leu552Ser", "carriers": bad, "quote": SOURCE[:60]}],
            SOURCE,
        )
        assert not res.accepted

    def test_null_is_not_a_rejection(self):
        res = validate_paper_response(
            self._gap(),
            [{"variant": "p.Leu552Ser", "carriers": None, "quote": None}],
            SOURCE,
        )
        assert not res.accepted and not res.rejected

    def test_ignores_field_we_did_not_ask_for(self):
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "affected": 12,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers, "
                    "of whom 12 had documented arrhythmic events",
                }
            ],
            SOURCE,
        )
        assert [c.field for c in res.accepted] == ["carriers"]


class TestParseResponse:
    def test_fenced_and_prose_prefixed(self):
        assert parse_response('```json\n{"variants": [{"variant": "X"}]}\n```')
        assert parse_response('Here you go: {"variants": []} -- done') == []

    @pytest.mark.parametrize(
        "raw", ["", "   ", "no json here", '{"other": 1}', '{"variants": "nope"}']
    )
    def test_rejects_unusable(self, raw):
        with pytest.raises(CountRecoveryError):
            parse_response(raw)


class TestWrite:
    def test_fills_null_and_logs(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers",
                }
            ],
            SOURCE,
        )
        stats = write_recovered_counts(db, [res], model="m", reasoning_effort="high")
        assert stats["counts_written"] == 1
        con = sqlite3.connect(db)
        assert (
            con.execute(
                "SELECT total_carriers_observed FROM penetrance_data WHERE variant_id=1"
            ).fetchone()[0]
            == 44
        )
        log = con.execute(
            "SELECT field, value, quote, model FROM count_recovery_log"
        ).fetchone()
        assert log[0] == "carriers" and log[1] == 44 and log[3] == "m"

    def test_never_overwrites_existing_count(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", 44, None, None)])
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 1,
                    "quote": "The p.Leu552Ser variant was seen once in this cohort",
                }
            ],
            "The p.Leu552Ser variant was seen once in this cohort",
        )
        stats = write_recovered_counts(db, [res], model="m", reasoning_effort=None)
        assert stats["counts_written"] == 0 and stats["already_populated_skipped"] == 1
        con = sqlite3.connect(db)
        assert (
            con.execute(
                "SELECT total_carriers_observed FROM penetrance_data WHERE variant_id=1"
            ).fetchone()[0]
            == 44
        )

    def test_updates_only_one_of_multiple_null_penetrance_rows(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        con = sqlite3.connect(db)
        con.execute(
            "INSERT INTO penetrance_data (variant_id,pmid,total_carriers_observed) "
            "VALUES (1,'111',NULL)"
        )
        con.commit()
        con.close()
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers",
                }
            ],
            SOURCE,
        )
        stats = write_recovered_counts(db, [res], model="m", reasoning_effort=None)
        assert stats["counts_written"] == 1
        con = sqlite3.connect(db)
        values = [
            row[0]
            for row in con.execute(
                "SELECT total_carriers_observed FROM penetrance_data "
                "WHERE variant_id=1 ORDER BY penetrance_id"
            )
        ]
        assert values == [44, None]

    def test_dry_run_writes_nothing(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers",
                }
            ],
            SOURCE,
        )
        write_recovered_counts(
            db, [res], model="m", reasoning_effort=None, dry_run=True
        )
        con = sqlite3.connect(db)
        assert (
            con.execute(
                "SELECT total_carriers_observed FROM penetrance_data WHERE variant_id=1"
            ).fetchone()[0]
            is None
        )
        assert not con.execute(
            "SELECT name FROM sqlite_master WHERE name='count_recovery_log'"
        ).fetchall()

    def test_inserts_row_when_penetrance_row_absent(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        con = sqlite3.connect(db)
        con.execute("DELETE FROM penetrance_data")
        con.commit()
        con.close()
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "quote": "The p.Leu552Ser variant was identified in 44 carriers",
                }
            ],
            SOURCE,
        )
        stats = write_recovered_counts(db, [res], model="m", reasoning_effort=None)
        assert stats["penetrance_rows_inserted"] == 1
        con = sqlite3.connect(db)
        assert (
            con.execute(
                "SELECT total_carriers_observed FROM penetrance_data WHERE variant_id=1"
            ).fetchone()[0]
            == 44
        )


class TestEndToEnd:
    def test_xhigh_gpt56_has_room_for_reasoning_before_json(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        caller = _stub(['{"variants": []}'])

        recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="azure_ai/gpt-5.6-luna",
            reasoning_effort="xhigh",
            dry_run=True,
        )

        assert caller.calls["kwargs"][0]["max_tokens"] == (
            COUNT_RECOVERY_XHIGH_MAX_TOKENS
        )
        assert caller.calls["kwargs"][0]["reasoning_effort"] == "xhigh"

    def test_xhigh_non_gpt56_keeps_compact_output_budget(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        caller = _stub(['{"variants": []}'])

        recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="azure_ai/grok-4.5",
            reasoning_effort="xhigh",
            dry_run=True,
        )

        assert caller.calls["kwargs"][0]["max_tokens"] == 8192

    def test_recovers_and_is_idempotent(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", "c.1655T>C", "111", None, None, None),
                (2, "p.Arg752Trp", None, "111", None, None, None),
                (3, "p.Gly628Ser", None, "111", None, None, None),
            ],
        )
        response = (
            '{"variants": ['
            '{"variant": "p.Leu552Ser c.1655T>C", "carriers": 44, "quote": '
            '"The p.Leu552Ser variant was identified in 44 carriers"},'
            '{"variant": "p.Arg752Trp", "carriers": 1, "quote": '
            '"The p.Arg752Trp variant was seen once in this cohort"},'
            '{"variant": "p.Gly628Ser", "carriers": null, "quote": null}]}'
        )
        caller = _stub([response])
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="stub",
            reasoning_effort="high",
        )
        assert stats["gaps_found"] == 3
        assert stats["counts_accepted"] == 2
        assert stats["counts_written"] == 2
        assert stats["llm_calls"] == 1

        # Second run: the two filled gaps are gone, only the ungrounded one remains.
        caller2 = _stub([response])
        stats2 = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller2,
            model="stub",
            reasoning_effort="high",
        )
        assert stats2["gaps_found"] == 1
        assert stats2["counts_written"] == 0

    def test_missing_source_is_reported_not_silently_dropped(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        caller = _stub(['{"variants": []}'])
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: None,
            llm_caller=caller,
            model="stub",
        )
        assert stats["papers_without_source"] == 1
        assert stats["llm_calls"] == 0 and stats["counts_written"] == 0

    def test_model_failure_on_one_paper_does_not_abort(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", None, "111", None, None, None),
                (2, "p.Arg752Trp", None, "222", None, None, None),
            ],
        )
        good = (
            '{"variants": [{"variant": "p.Arg752Trp", "carriers": 1, "quote": '
            '"The p.Arg752Trp variant was seen once in this cohort"}]}'
        )
        caller = _stub(["this is not json", good])
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="stub",
        )
        assert stats["papers_failed"] == 1
        assert stats["counts_written"] == 1

    def test_batches_respect_max_variants_per_call(self, tmp_path):
        rows = [(i, f"p.Ala{i}Val", None, "111", None, None, None) for i in range(1, 6)]
        db = _make_db(tmp_path, rows)
        caller = _stub(['{"variants": []}'])
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="stub",
            max_variants_per_call=2,
        )
        assert stats["llm_calls"] == 3

    def test_failed_paper_is_counted_once_across_failed_batches(self, tmp_path):
        rows = [(i, f"p.Ala{i}Val", None, "111", None, None, None) for i in range(1, 5)]
        db = _make_db(tmp_path, rows)
        caller = _stub(["not json"])
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="stub",
            max_variants_per_call=1,
        )
        assert stats["papers_failed"] == 1
        assert stats["batch_failures"] == 4

    def test_quote_from_truncated_unseen_tail_is_rejected(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        visible = "Methods without a count. "
        hidden = "The p.Leu552Ser variant was identified in 44 carriers."
        caller = _stub(
            [
                '{"variants": [{"variant": "p.Leu552Ser", "carriers": 44, '
                f'"quote": "{hidden}"}}]}}'
            ]
        )
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: visible + hidden,
            llm_caller=caller,
            model="stub",
            max_source_chars=len(visible),
        )
        assert stats["counts_accepted"] == 0
        assert stats["counts_rejected"] == 1
        assert stats["counts_written"] == 0


class TestPrompt:
    def test_lists_only_requested_variants_and_forbids_new_ones(self):
        p = build_prompt(
            "KCNH2",
            "111",
            SOURCE,
            [VariantGap(1, "p.Leu552Ser", ["carriers"])],
        )
        flat = " ".join(p.lower().split())  # the prompt is hard-wrapped
        assert "p.Leu552Ser" in p
        assert "do not look for new variants" in flat
        assert "verbatim" in flat

    def test_marks_truncation(self):
        p = build_prompt(
            "KCNH2",
            "111",
            "x" * 500,
            [VariantGap(1, "p.Leu552Ser", ["carriers"])],
            max_source_chars=100,
        )
        assert "truncated" in p


class TestDuplicateVariantCollapse:
    """One variant stored under several variant_ids must not be double-counted.

    Production writes the same variant under different notations per layer (e.g.
    `p.Leu552Ser c.1655T>C` from llm_text and `p.L552S` from a PubTator
    linkage). If the extractor counted one of them, filling the other would add
    a second carrier count for one variant in one paper.
    """

    def test_counted_sibling_suppresses_the_gap(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", "c.1655T>C", "111", 42, None, None),
                (341, "p.L552S", None, "111", None, None, None),
            ],
        )
        assert find_count_gaps(db, "KCNH2", fields=("carriers",)) == []

    def test_counted_linkage_sibling_suppresses_paper_derived_gap(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", "c.1655T>C", "111", None, None, None),
                (341, "p.L552S", None, "111", 42, None, None),
            ],
        )
        con = sqlite3.connect(db)
        con.execute(
            "UPDATE variant_papers SET source_layer='clinvar' WHERE variant_id=341"
        )
        con.commit()
        con.close()
        assert find_count_gaps(db, "KCNH2", fields=("carriers",)) == []

    def test_linkage_only_rows_help_dedupe_but_are_never_write_targets(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.L552S", None, "111", None, None, None),
                (341, "p.Leu552Ser", "c.1655T>C", "111", None, None, None),
            ],
        )
        con = sqlite3.connect(db)
        con.execute(
            "UPDATE variant_papers SET source_layer='clinvar' WHERE variant_id=341"
        )
        con.commit()
        con.close()
        gaps = find_count_gaps(db, "KCNH2", fields=("carriers",))
        (variant,) = gaps[0].variants
        assert variant.variant_id == 1
        assert variant.paper_derived
        assert variant.duplicate_ids == [341]

    def test_uncounted_cluster_yields_one_gap_on_the_richest_notation(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", "c.1655T>C", "111", None, None, None),
                (341, "p.L552S", None, "111", None, None, None),
            ],
        )
        gaps = find_count_gaps(db, "KCNH2", fields=("carriers",))
        assert len(gaps) == 1
        (v,) = gaps[0].variants
        assert v.variant_id == 1
        assert v.duplicate_ids == [341]

    def test_protein_only_and_cdna_only_same_codon_are_conservatively_collapsed(
        self, tmp_path
    ):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", None, "111", None, None, None),
                (341, None, "c.1655T>C", "111", None, None, None),
            ],
        )
        gaps = find_count_gaps(db, "KCNH2", fields=("carriers",))
        assert len(gaps) == 1
        assert len(gaps[0].variants) == 1
        variant = gaps[0].variants[0]
        assert {variant.variant_id, *variant.duplicate_ids} == {1, 341}

    def test_counted_cdna_only_alias_suppresses_protein_only_gap(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", None, "111", None, None, None),
                (341, None, "c.1655T>C", "111", 42, None, None),
            ],
        )
        assert find_count_gaps(db, "KCNH2", fields=("carriers",)) == []

    def test_genuinely_different_variants_stay_separate(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", None, "111", None, None, None),
                (2, "p.Arg752Trp", None, "111", None, None, None),
            ],
        )
        gaps = find_count_gaps(db, "KCNH2", fields=("carriers",))
        assert sorted(v.variant_id for v in gaps[0].variants) == [1, 2]

    def test_write_touches_only_the_representative_row(self, tmp_path):
        db = _make_db(
            tmp_path,
            [
                (1, "p.Leu552Ser", "c.1655T>C", "111", None, None, None),
                (341, "p.L552S", None, "111", None, None, None),
            ],
        )
        caller = _stub(
            [
                '{"variants": [{"variant": "p.Leu552Ser c.1655T>C", "carriers": 44, '
                '"quote": "The p.Leu552Ser variant was identified in 44 carriers"}]}'
            ]
        )
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda p: SOURCE,
            llm_caller=caller,
            model="stub",
        )
        assert stats["counts_written"] == 1
        con = sqlite3.connect(db)
        rows = dict(
            con.execute(
                "SELECT variant_id, total_carriers_observed FROM penetrance_data"
            ).fetchall()
        )
        assert rows[1] == 44, "representative row should be filled"
        assert rows[341] is None, "duplicate row must stay NULL (no double count)"


def test_find_count_gaps_restricts_to_requested_gene(tmp_path):
    db = _make_db(
        tmp_path,
        [
            (1, "p.Leu552Ser", None, "111", None, None, None),
            (2, "p.Arg752Trp", None, "222", None, None, None),
        ],
    )
    con = sqlite3.connect(db)
    con.execute("UPDATE variants SET gene_symbol='SCN5A' WHERE variant_id=2")
    con.commit()
    con.close()
    gaps = find_count_gaps(db, "KCNH2")
    assert [(g.pmid, [v.variant_id for v in g.variants]) for g in gaps] == [
        ("111", [1])
    ]


class TestRecoveryProvenanceAndDurability:
    """The trust gate reads role out of count_provenance; recovery must write it."""

    QUOTE = "The p.Leu552Ser variant was identified in 44 carriers"

    def _db_with_provenance(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        con = sqlite3.connect(db)
        try:
            con.execute("ALTER TABLE variant_papers ADD COLUMN count_provenance TEXT")
            con.execute(
                "ALTER TABLE penetrance_data ADD COLUMN trust_tier TEXT "
                "DEFAULT 'trusted'"
            )
            con.execute("ALTER TABLE penetrance_data ADD COLUMN trust_reasons TEXT")
            con.execute("ALTER TABLE penetrance_data ADD COLUMN trust_sources TEXT")
            con.execute(
                "ALTER TABLE penetrance_data ADD COLUMN trust_rule_version TEXT"
            )
            con.commit()
        finally:
            con.close()
        return db

    def _accepted(self):
        gap = PaperGap("KCNH2", "111", [VariantGap(1, "p.Leu552Ser", ["carriers"])])
        res = validate_paper_response(
            gap,
            [
                {
                    "variant": "p.Leu552Ser",
                    "carriers": 44,
                    "count_role": "per_variant_carriers",
                    "quote": self.QUOTE,
                }
            ],
            SOURCE,
        )
        assert res.accepted, res.rejected
        return res

    def test_write_persists_role_and_locator_into_count_provenance(self, tmp_path):
        db = self._db_with_provenance(tmp_path)
        stats = write_recovered_counts(
            db, [self._accepted()], model="m", reasoning_effort="high"
        )
        assert stats["counts_written"] == 1
        assert stats["count_provenance_written"] == 1

        con = sqlite3.connect(db)
        try:
            raw = con.execute(
                "SELECT count_provenance FROM variant_papers WHERE variant_id=1"
            ).fetchone()[0]
            log = con.execute(
                "SELECT count_role, evidence_locator, model_declared_role "
                "FROM count_recovery_log"
            ).fetchone()
        finally:
            con.close()
        provenance = json.loads(raw)
        assert provenance["carriers_count_type"] == "per_variant_carriers"
        assert provenance["carriers_source"] == "count_recovery"
        assert provenance["carriers_recovery"]["evidence_kind"] == "prose_sentence"
        assert provenance["carriers_recovery"]["quote_sha256"]
        assert log[0] == "per_variant_carriers"
        assert json.loads(log[1])["kind"] == "prose_sentence"
        assert log[2] == "per_variant_carriers"

    def test_stale_cohort_total_provenance_is_overwritten_for_this_field(
        self, tmp_path
    ):
        """A correct per-variant fill must not inherit a stale denominator role."""
        db = self._db_with_provenance(tmp_path)
        con = sqlite3.connect(db)
        try:
            con.execute(
                "UPDATE variant_papers SET count_provenance=? WHERE variant_id=1",
                (
                    json.dumps(
                        {
                            "carriers_count_type": "cohort_total",
                            "affected_count_type": "per_variant_affected",
                        }
                    ),
                ),
            )
            con.commit()
        finally:
            con.close()

        write_recovered_counts(db, [self._accepted()], model="m", reasoning_effort=None)

        con = sqlite3.connect(db)
        try:
            provenance = json.loads(
                con.execute(
                    "SELECT count_provenance FROM variant_papers WHERE variant_id=1"
                ).fetchone()[0]
            )
        finally:
            con.close()
        assert provenance["carriers_count_type"] == "per_variant_carriers"
        # Other fields' provenance is merged, not clobbered.
        assert provenance["affected_count_type"] == "per_variant_affected"

    def test_recovered_row_lands_in_quarantine_for_the_trust_gate_to_promote(
        self, tmp_path
    ):
        from pipeline.trust_gate import evaluate_fact

        db = self._db_with_provenance(tmp_path)
        stats = write_recovered_counts(
            db, [self._accepted()], model="m", reasoning_effort=None
        )
        assert stats["landed_in_quarantine"] == 1

        con = sqlite3.connect(db)
        try:
            tier, reasons, sources = con.execute(
                "SELECT trust_tier, trust_reasons, trust_sources FROM penetrance_data"
            ).fetchone()
            provenance = json.loads(
                con.execute(
                    "SELECT count_provenance FROM variant_papers WHERE variant_id=1"
                ).fetchone()[0]
            )
        finally:
            con.close()
        assert tier == "quarantine"
        assert json.loads(reasons) == ["count_recovery_pending_review"]
        assert json.loads(sources) == ["count_recovery"]
        # With the role proven, the gate has no reason to keep it quarantined.
        assert "recovered_count_unverified" not in evaluate_fact(
            {"carriers": 44}, provenance
        )

    def test_dry_run_reports_what_a_real_run_would_actually_write(self, tmp_path):
        """Dry run used to count every accepted value, even onto occupied slots."""
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", 44, None, None)])
        res = self._accepted()

        dry = write_recovered_counts(
            db, [res], model="m", reasoning_effort=None, dry_run=True
        )
        real = write_recovered_counts(db, [res], model="m", reasoning_effort=None)

        assert dry["counts_written"] == real["counts_written"] == 0
        assert (
            dry["already_populated_skipped"] == real["already_populated_skipped"] == 1
        )
        assert dry["dry_run"] is True

    def test_dry_run_does_not_create_the_audit_table_or_touch_the_database(
        self, tmp_path
    ):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        before = db.read_bytes()

        dry = write_recovered_counts(
            db, [self._accepted()], model="m", reasoning_effort=None, dry_run=True
        )

        assert dry["counts_written"] == 1  # it WOULD write this one
        assert db.read_bytes() == before
        con = sqlite3.connect(db)
        try:
            tables = {
                row[0]
                for row in con.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
        finally:
            con.close()
        assert "count_recovery_log" not in tables

    def test_recover_counts_backs_up_before_mutating(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        response = SimpleNamespace(
            choices=[
                SimpleNamespace(
                    message=SimpleNamespace(
                        content=json.dumps(
                            {
                                "variants": [
                                    {
                                        "variant": "p.Leu552Ser",
                                        "carriers": 44,
                                        "count_role": "per_variant_carriers",
                                        "quote": self.QUOTE,
                                    }
                                ]
                            }
                        )
                    )
                )
            ]
        )

        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda pmid: SOURCE,
            llm_caller=lambda **kwargs: response,
            model="m",
        )

        backup = Path(str(db) + ".before_count_recovery.db")
        assert stats["backup_path"] == str(backup)
        assert backup.is_file()
        # The backup predates the write: its slot is still NULL.
        con = sqlite3.connect(backup)
        try:
            assert (
                con.execute(
                    "SELECT total_carriers_observed FROM penetrance_data"
                ).fetchone()[0]
                is None
            )
        finally:
            con.close()
        assert stats["counts_written"] == 1

    def test_stats_are_json_serializable(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        response = SimpleNamespace(
            choices=[
                SimpleNamespace(
                    message=SimpleNamespace(
                        content=json.dumps(
                            {
                                "variants": [
                                    {
                                        "variant": "p.Leu552Ser",
                                        "carriers": 44,
                                        "count_role": "per_variant_carriers",
                                        "quote": self.QUOTE,
                                    }
                                ]
                            }
                        )
                    )
                )
            ]
        )
        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda pmid: SOURCE,
            llm_caller=lambda **kwargs: response,
            model="m",
            dry_run=True,
        )
        stats.pop("result_objects")
        payload = json.loads(json.dumps(stats))  # used to raise on dataclasses
        assert payload["results"][0]["accepted"][0]["count_role"] == (
            "per_variant_carriers"
        )
        assert payload["paper_derived_only"] is True

    def test_total_failure_is_reported_for_the_caller_to_escalate(self, tmp_path):
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])

        def boom(**kwargs):
            raise RuntimeError("provider down")

        stats = recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda pmid: SOURCE,
            llm_caller=boom,
            model="m",
        )

        assert stats["all_batches_failed"] is True
        assert stats["papers_failed"] == 1
        assert stats["counts_accepted"] == 0


class TestPaperDerivedOnly:
    def test_include_linkage_rows_is_reachable(self, tmp_path):
        """`paper_derived_only` was never plumbed, so the documented mode was dead."""
        db = _make_db(tmp_path, [(1, "p.Leu552Ser", None, "111", None, None, None)])
        con = sqlite3.connect(db)
        try:
            con.execute(
                "UPDATE variant_papers SET source_layer='clinvar' WHERE variant_id=1"
            )
            con.commit()
        finally:
            con.close()

        assert find_count_gaps(db, "KCNH2") == []
        (paper,) = find_count_gaps(db, "KCNH2", paper_derived_only=False)
        assert [v.notation for v in paper.variants] == ["p.Leu552Ser"]

        seen: list = []
        recover_counts(
            db,
            "KCNH2",
            source_for_pmid=lambda pmid: seen.append(pmid) or None,
            llm_caller=lambda **kwargs: None,
            model="m",
            paper_derived_only=False,
            dry_run=True,
        )
        assert seen == ["111"]
