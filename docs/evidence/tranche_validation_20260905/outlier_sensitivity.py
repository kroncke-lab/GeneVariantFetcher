"""Post-hoc count-row influence checks; never change locked scores or gates."""

import json

from summarize import FIELDS, OUT, count_summary, key, load_arm, read


def main():
    arms = {
        (number, arm): load_arm(number, arm)
        for number in ("02", "03")
        for arm in ("baseline", "candidate")
    }
    exposure = read(OUT / "prior_exposure.json")
    novel = {
        key(row)
        for number in ("02", "03")
        for row in exposure[f"tranche_{number}"]["previously_unscored_attempts"]
    }
    result = {
        "classification": "Post-hoc influence sensitivity; descriptive, not an acceptance test",
        "excluded_row": {"gene": "SCN5A", "pmid": "20129283", "gold_variant": "H558R"},
        "reason": "408 is reference-matching but person-versus-allele unit and study attribution remain unresolved.",
        "identity_metrics_and_official_gold_unchanged": True,
        "comparisons": {},
        "input_sha256": {data["run_id"]: data["sha256"] for data in arms.values()},
    }
    for name, numbers, allowed in (
        ("tranche_03", ("03",), None),
        ("pooled", ("02", "03"), None),
        ("pooled_previously_unscored", ("02", "03"), novel),
    ):
        pair = {}
        for arm in ("baseline", "candidate"):
            rows = [
                row
                for number in numbers
                for row in arms[number, arm]["rows"]
                if allowed is None or key(row) in allowed
            ]
            removed = [
                row
                for row in rows
                if key(row) == ("SCN5A", "20129283") and row["gold_variant"] == "H558R"
            ]
            assert len(removed) == 3
            retained = [row for row in rows if row not in removed]
            pair[arm] = {
                "removed_rows": removed,
                "counts": {
                    field: count_summary(
                        [row for row in retained if row["measure"] == field]
                    )
                    for field in FIELDS
                },
                "affected_unaffected_combined": count_summary(
                    [
                        row
                        for row in retained
                        if row["measure"] in ("affected", "unaffected")
                    ]
                ),
            }
        old, new = (
            pair[arm]["affected_unaffected_combined"]["absolute_error_sum"]
            for arm in ("baseline", "candidate")
        )
        pair["combined_au_relative_error_reduction"] = (old - new) / old
        result["comparisons"][name] = pair
    # A separate whole-paper carrier influence check, also post-hoc.
    carriers = {}
    for arm in ("baseline", "candidate"):
        rows = [
            row
            for number in ("02", "03")
            for row in arms[number, arm]["rows"]
            if row["measure"] == "carriers"
            and key(row) not in {("SCN5A", "27566755"), ("SCN5A", "20129283")}
        ]
        carriers[arm] = count_summary(rows)
    old, new = (
        carriers[arm]["absolute_error_sum"] for arm in ("baseline", "candidate")
    )
    carriers["relative_error_reduction"] = (old - new) / old
    result["carrier_sensitivity_excluding_two_whole_papers"] = carriers
    (OUT / "outlier_sensitivity.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        {
            name: data["combined_au_relative_error_reduction"]
            for name, data in result["comparisons"].items()
        }
    )


if __name__ == "__main__":
    main()
