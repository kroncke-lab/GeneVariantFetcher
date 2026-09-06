"""Acceptance-policy tests with synthetic evidence, no gold and no API calls."""

from clinical_reader import reconcile


def paper(carriers=None):
    return {
        "gene": "KCNH2",
        "pmid": "111",
        "variants": [
            {
                "variant": "p.Leu552Ser",
                "carriers": carriers,
                "affected": None,
                "unaffected": None,
                "evidence": "Initial paper evidence",
            }
        ],
    }


def proposal(
    value=44,
    basis="explicit",
    quote="The p.Leu552Ser variant was identified in 44 carriers",
):
    return {
        "variants": [
            {
                "variant": "p.L552S",
                "carriers": value,
                "carriers_evidence": {
                    "quote": quote,
                    "basis": basis,
                    "scope": "current_study",
                    "count_role": "per_variant_carriers",
                },
            }
        ]
    }


def test_source_bound_fill_and_baseline_immutable():
    baseline = paper()
    payload = proposal()
    source = payload["variants"][0]["carriers_evidence"]["quote"]
    accepted, raw, audit = reconcile(baseline, payload, source)
    assert accepted[0]["carriers"] == raw[0]["carriers"] == 44
    assert len(audit["accepted"]) == 1
    assert baseline["variants"][0]["carriers"] is None


def test_conflicting_filled_value_is_only_a_contradiction():
    accepted, raw, audit = reconcile(
        paper(40), proposal(), "The p.Leu552Ser variant was identified in 44 carriers"
    )
    assert accepted[0]["carriers"] == raw[0]["carriers"] == 40
    assert len(audit["contradictions"]) == 1
    assert not audit["accepted"]


def test_derived_or_unquoted_proposal_is_raw_only():
    for payload in [
        proposal(basis="individual_rows"),
        proposal(quote="invented quote"),
    ]:
        accepted, raw, audit = reconcile(
            paper(), payload, "The p.Leu552Ser variant was identified in 44 carriers"
        )
        assert accepted[0]["carriers"] is None and raw[0]["carriers"] == 44
        assert audit["rejected"] and not audit["accepted"]


def test_disagreeing_proposals_do_not_choose_a_winner():
    payload = {"variants": proposal(44)["variants"] + proposal(45)["variants"]}
    accepted, raw, audit = reconcile(
        paper(), payload, "The p.Leu552Ser variant was identified in 44 carriers"
    )
    assert accepted[0]["carriers"] is None and raw[0]["carriers"] is None
    assert audit["rejected"][0]["reason"] == "conflicting independent proposals"


def test_rolled_back_fill_does_not_leave_accepted_evidence():
    source = "The p.Leu552Ser variant was identified in 3 affected carriers."
    payload = {
        "variants": [
            {
                "variant": "p.L552S",
                "affected": 3,
                "affected_evidence": {
                    "quote": source,
                    "basis": "explicit",
                    "scope": "current_study",
                    "count_role": "per_variant_affected",
                },
            }
        ]
    }
    accepted, raw, audit = reconcile(paper(2), payload, source)
    assert accepted[0]["affected"] is None
    assert accepted[0]["evidence"] == "Initial paper evidence"
    assert raw[0]["affected"] == 3
    assert any("person totals conflict" in r["reason"] for r in audit["rejected"])


def test_original_row_index_survives_single_target_validation():
    baseline = paper()
    baseline["variants"].insert(
        0,
        {
            "variant": "p.Arg752Trp",
            "carriers": None,
            "affected": None,
            "unaffected": None,
            "evidence": "Different variant",
        },
    )
    source = "The p.Leu552Ser variant was identified in 44 carriers"
    accepted, _, audit = reconcile(baseline, proposal(), source)
    assert accepted[0]["carriers"] is None
    assert accepted[1]["carriers"] == 44
    assert audit["accepted"][0]["variant_id"] == 1


def test_combined_literal_cdna_and_protein_notation_maps_uniquely():
    baseline = paper()
    baseline["variants"][0]["variant"] = "p.Leu552Ser c.1655T>C"
    payload = proposal()
    payload["variants"][0]["variant"] = "c.1655T>C (p.Leu552Ser)"
    accepted, _, audit = reconcile(
        baseline, payload, "The p.Leu552Ser variant was identified in 44 carriers"
    )
    assert accepted[0]["carriers"] == 44
    assert len(audit["accepted"]) == 1


def test_combined_notation_cannot_choose_between_two_retained_variants():
    baseline = paper()
    baseline["variants"].append(
        {
            "variant": "p.Arg752Trp",
            "carriers": None,
            "affected": None,
            "unaffected": None,
            "evidence": "second variant",
        }
    )
    payload = proposal()
    payload["variants"][0]["variant"] = "p.Leu552Ser and p.Arg752Trp"
    accepted, _, audit = reconcile(
        baseline, payload, "The p.Leu552Ser variant was identified in 44 carriers"
    )
    assert all(row["carriers"] is None for row in accepted)
    assert audit["rejected"][0]["matching_rows"] == 2
