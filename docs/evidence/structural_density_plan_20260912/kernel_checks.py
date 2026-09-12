"""Reproduce the proposed kernel preflight without any structure or model jobs.

Run from any directory with Python 3.10+; no third-party packages are required.
Outputs are adjacent to this script. The sine function preserves the historical
formula, replacing rounded 3.14 boundaries with math.pi. Its midpoint is an
absolute weight of 0.5; for small midpoints K(0) is below 1. The normalized
sigmoid has K(0)=1 and K(h)=0.5 exactly (up to floating-point error).
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent
HALF_DISTANCES = (2.0, 3.0, 5.0)
POLYMERS = (
    ("historical_baseline", 3.8, 0.5),
    ("existing_ppa", 5.5, 0.55),
)


def historical_sine(distance: float, midpoint: float) -> float:
    """Historical sin_dist, with pi replacing its rounded 3.14 boundaries."""
    if distance < 0 or midpoint <= 0:
        raise ValueError("Distance must be nonnegative and midpoint positive")
    if distance <= midpoint - math.pi:
        return 1.0
    if distance <= midpoint + math.pi:
        return 0.5 - 0.5 * math.sin((distance - midpoint) / 2)
    return 0.0


def normalized_sigmoid(distance: float, half_distance: float) -> float:
    """Unit weight at zero and half the zero-distance weight at half_distance."""
    if distance < 0 or half_distance <= 0:
        raise ValueError("Distance must be nonnegative and half-distance positive")
    decay = math.exp(-math.log(3) * distance / half_distance)
    return 2 * decay / (1 + decay)


def polymer_distance(separation: int, b: float, nu: float) -> float:
    return b * abs(separation) ** nu


def make_rows() -> list[dict]:
    distances = {
        "zero": 0.0,
        "2A": 2.0,
        "3A": 3.0,
        "adjacent_baseline": 3.8,
        "two_apart_baseline": 3.8 * math.sqrt(2),
        "three_apart_baseline": 3.8 * math.sqrt(3),
    }
    for midpoint in HALF_DISTANCES:
        if midpoint - math.pi >= 0:
            distances[f"sine_a{midpoint:g}_plateau_end"] = midpoint - math.pi
        distances[f"sine_a{midpoint:g}_cutoff"] = midpoint + math.pi
        distances[f"sine_a{midpoint:g}_beyond_cutoff"] = midpoint + math.pi + 1e-9
    positions = [
        ("distance_probe", label, "", "", "", distance)
        for label, distance in distances.items()
    ]
    for label, b, nu in POLYMERS:
        for separation in (*range(11), 20, 30):
            positions.append(
                (
                    "polymer_probe",
                    label,
                    b,
                    nu,
                    separation,
                    polymer_distance(separation, b, nu),
                )
            )
    rows = []
    kernels = (
        ("historical_sine_pi_cleanup", historical_sine),
        ("normalized_sigmoid", normalized_sigmoid),
    )
    for kind, label, b, nu, separation, distance in positions:
        for kernel_name, function in kernels:
            for half_distance in HALF_DISTANCES:
                zero_weight = function(0.0, half_distance)
                value = function(distance, half_distance)
                rows.append(
                    {
                        "probe_type": kind,
                        "probe_label": label,
                        "polymer_b_angstrom": b,
                        "polymer_nu": nu,
                        "sequence_separation": separation,
                        "distance_angstrom": distance,
                        "kernel": kernel_name,
                        "midpoint_or_half_distance_angstrom": half_distance,
                        "weight": value,
                        "weight_at_zero": zero_weight,
                        "relative_to_zero_weight": value / zero_weight,
                    }
                )
    return rows


def check() -> dict:
    parameters = []
    for midpoint in HALF_DISTANCES:
        sine_zero = historical_sine(0.0, midpoint)
        assert math.isclose(historical_sine(midpoint, midpoint), 0.5, abs_tol=1e-12)
        assert historical_sine(midpoint + math.pi, midpoint) == 0.0
        assert historical_sine(midpoint + math.pi + 1e-9, midpoint) == 0.0
        if midpoint >= math.pi:
            assert historical_sine(midpoint - math.pi, midpoint) == 1.0
        assert normalized_sigmoid(0.0, midpoint) == 1.0
        assert math.isclose(normalized_sigmoid(midpoint, midpoint), 0.5, abs_tol=1e-12)
        for function in (historical_sine, normalized_sigmoid):
            weights = [function(index / 100, midpoint) for index in range(3001)]
            assert all(0 <= value <= 1 for value in weights)
            assert all(left >= right for left, right in zip(weights, weights[1:]))
        parameters.append(
            {
                "midpoint_or_half_distance_angstrom": midpoint,
                "sine_zero_weight": sine_zero,
                "sine_midpoint_weight": historical_sine(midpoint, midpoint),
                "sine_midpoint_relative_weight": 0.5 / sine_zero,
                "sine_cutoff_angstrom": midpoint + math.pi,
                "sigmoid_zero_weight": normalized_sigmoid(0, midpoint),
                "sigmoid_half_distance_weight": normalized_sigmoid(midpoint, midpoint),
            }
        )
    support = {}
    for label, b, nu in POLYMERS:
        support[label] = [
            separation
            for separation in range(101)
            if historical_sine(polymer_distance(separation, b, nu), 3.0) > 0
        ]
    assert support["historical_baseline"] == [0, 1, 2]
    assert support["existing_ppa"] == [0, 1]
    return {
        "passed": True,
        "scope": "Numerical preflight only; no structure, density fit, or validation run",
        "source": (
            "https://github.com/kroncke-lab/Bayes_BrS1_Penetrance/"
            "blob/master/func_dist_seq.R"
        ),
        "sine_numerical_change": "Use math.pi in boundaries instead of rounded 3.14",
        "sine_formula": "1 below a-pi; 0.5-0.5*sin((d-a)/2) through a+pi; 0 above",
        "sigmoid_formula": "2/(1+exp(log(3)*d/h))",
        "parameters": parameters,
        "polymer_models": [
            {"label": label, "b_angstrom": b, "nu": nu} for label, b, nu in POLYMERS
        ],
        "nonzero_sine_a3_sequence_separations_tested_0_through_100": support,
        "notes": [
            "The sine midpoint is absolute half weight, not exactly half of K(0) when a<pi.",
            "Sine a=3 K(0) is approximately 0.998747; a=2 K(0) is approximately 0.920735.",
            "N=0 may retain other substitutions at the target residue after variant-only LOO.",
            "No remaining donors means missing density and zero support; do not widen silently.",
            "Positive sigmoid weights at longer distances do not establish useful support.",
        ],
    }


def main() -> None:
    checks = check()
    rows = make_rows()
    with (OUT / "kernel_sensitivity.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    checks["csv_rows"] = len(rows)
    (OUT / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(f"Kernel checks passed; wrote {len(rows)} sensitivity rows")


if __name__ == "__main__":
    main()
