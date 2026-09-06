"""Adversarial checks for the opened-source diagnostic's aggregation boundary."""

import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def scorer(monkeypatch, tmp_path):
    directory = (
        Path(__file__).resolve().parents[2] / "docs/evidence/model_followup_20260906"
    )
    monkeypatch.syspath_prepend(str(directory))
    spec = importlib.util.spec_from_file_location(
        "bounded_roster_score", directory / "score.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "HERE", tmp_path)
    return module


def person(ident, endpoint="positive"):
    return dict(
        sources=["L" + ident],
        variants=["V1"],
        family=None,
        person=ident,
        genotype="Y",
        phenotype="disease",
        index_case=None,
        endpoint_status=endpoint,
    )


def grade(scorer, expected, observed, finish="stop"):
    packet = dict(
        name="packet", source_units={r["sources"][0]: r["person"] for r in expected}
    )
    path = scorer.HERE / "response.json"
    path.write_text(
        json.dumps(
            dict(
                status="returned",
                seconds=1,
                accounted_usd=0.01,
                response=dict(
                    choices=[
                        dict(
                            finish_reason=finish,
                            message=dict(
                                content=json.dumps(
                                    dict(
                                        columns=scorer.COLS,
                                        rows=[
                                            [r[c] for c in scorer.COLS]
                                            for r in observed
                                        ],
                                        limitations=[],
                                    )
                                )
                            ),
                        )
                    ]
                ),
            )
        )
    )
    return scorer.grade(
        packet, expected, dict(path="response.json", sha256=scorer.sha(path))
    )


def test_partial_roster_is_not_a_complete_variant_count(scorer):
    out = grade(scorer, [person("1"), person("2")], [person("1")])
    assert out["counts"]["v1"]["carriers"] == 1
    assert out["missing_people"] == 1
    assert out["source_adjudicated_complete_variant_counts"] == {}


def test_duplicate_index_person_never_adds_an_extra_carrier(scorer):
    out = grade(scorer, [person("1")], [person("1"), person("1")])
    assert out["counts"]["v1"]["carriers"] == 1
    assert out["duplicate_people"]
    assert out["source_adjudicated_complete_variant_counts"] == {}


def test_invented_person_cannot_pass_by_matching_a_total(scorer):
    out = grade(scorer, [person("1"), person("2")], [person("1"), person("3")])
    assert out["counts"] == out["reference_counts"]
    assert out["extra_people"] == out["missing_people"] == 1
    assert out["source_adjudicated_complete_variant_counts"] == {}


def test_wrong_endpoint_preserves_carrier_but_fails_evidence_acceptance(scorer):
    out = grade(scorer, [person("1")], [person("1", "negative")])
    assert out["counts"]["v1"]["carriers"] == 1
    assert out["endpoint_errors"]
    assert out["source_adjudicated_complete_variant_counts"] == {}


def test_length_failure_stays_in_count_coverage_denominator(scorer):
    out = grade(scorer, [person("1")], [], finish="length")
    assert out["status"] == "length"
    assert len(out["count_comparison"]) == 3
    assert not any(c["exact"] for c in out["count_comparison"])


def test_correct_empty_negative_control(scorer):
    out = grade(scorer, [], [])
    assert out["correct_empty_roster"] is True
    assert out["counts"] == {}
