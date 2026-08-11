import csv
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
HISTORICAL = (
    ROOT / "benchmarks/count_semantics_eval/runs/20260810_luna_xhigh_56/metrics.json"
)
ACTIVE_DIR = ROOT / "benchmarks/count_semantics_eval/runs/20260811_collaborator_gold_50"


def _totals(rows, field):
    return {
        key: sum(int(row[f"{field}_{key}"]) for row in rows)
        for key in ("assertions", "predicted", "absolute_error")
    }


def test_active_50_projection_recomputes_from_frozen_56_contributions():
    historical = json.loads(HISTORICAL.read_text())
    active = json.loads((ACTIVE_DIR / "metrics.json").read_text())
    with (ACTIVE_DIR / "brca2_paper_contributions.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    kept = [row for row in rows if row["active_collaborator_gold"] == "true"]

    for field in ("carriers", "affected", "unaffected"):
        all_brca2 = _totals(rows, field)
        kept_brca2 = _totals(kept, field)
        for phase in ("before", "after"):
            old = historical["count_fields"][field][phase]
            expected = {
                key: old[key] - all_brca2[key] + kept_brca2[key]
                for key in ("assertions", "predicted", "absolute_error")
            }
            actual = active["count_fields"][field][phase]
            assert {key: actual[key] for key in expected} == expected
            assert actual["mae"] == pytest.approx(
                expected["absolute_error"] / expected["predicted"]
            )
            assert actual["recall"] == pytest.approx(
                expected["predicted"] / expected["assertions"]
            )

    assert active["active_cohort"]["papers"] == 50
    assert active["result"]["brca2_2_carrier_mae"] == pytest.approx(4 / 3)
    assert active["projection"]["predictions_changed"] is False
