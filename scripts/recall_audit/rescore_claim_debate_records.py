#!/usr/bin/env python3
"""Rescore stored claim-debate records against current gold columns."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

try:
    from scripts.recall_audit.common import (
        canonical_variant,
        parse_int,
        repo_path,
        write_json,
    )
    from scripts.recall_audit.run_claim_verification_pilot import (
        FIELDS,
        GOLD_FIELD_MAP,
        best_gold_row,
        evaluate_field,
        gold_index,
        gold_value,
        gold_value_source,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import canonical_variant, parse_int, repo_path, write_json
    from run_claim_verification_pilot import (
        FIELDS,
        GOLD_FIELD_MAP,
        best_gold_row,
        evaluate_field,
        gold_index,
        gold_value,
        gold_value_source,
    )


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def record_db_values(record: dict[str, Any]) -> dict[str, int | None]:
    values: dict[str, int | None] = {}
    for item in record.get("baseline_field_evaluations", []):
        field = str(item.get("field") or "")
        if field in FIELDS:
            values[field] = parse_int(item.get("db_value"))
    return values


def rescore_items(
    *,
    values: dict[str, int | None],
    gold_row: dict[str, Any],
    gold_value_set: str,
    value_source: str,
    verification: dict[str, Any],
) -> list[dict[str, Any]]:
    return [
        evaluate_field(
            field=field,
            db_value=values.get(field),
            gold_value=gold_value(gold_row, field, gold_value_set=gold_value_set),
            original_gold_value=parse_int(gold_row.get(GOLD_FIELD_MAP[field])),
            gold_value_source=value_source,
            verification=verification,
        )
        for field in FIELDS
    ]


def rescore_record(record: dict[str, Any], *, gold_value_set: str) -> dict[str, Any]:
    rescored = copy.deepcopy(record)
    gene = str(rescored.get("gene") or "").upper()
    pmid = str(rescored.get("pmid") or "")
    variant = str(rescored.get("variant") or "")
    values = record_db_values(rescored)

    gold_rows = gold_index(gene, pmid).get(canonical_variant(variant), [])
    gold_row = best_gold_row(
        gold_rows,
        values,
        gold_value_set=gold_value_set,
    )
    value_source = gold_value_source(gold_row, gold_value_set=gold_value_set)

    rescored["matched_gold_row"] = gold_row
    rescored["gold_value_set"] = gold_value_set
    rescored["gold_value_source"] = value_source
    rescored["baseline_field_evaluations"] = rescore_items(
        values=values,
        gold_row=gold_row,
        gold_value_set=gold_value_set,
        value_source=value_source,
        verification=rescored.get("baseline_verification") or {},
    )
    rescored["debate_field_evaluations"] = rescore_items(
        values=values,
        gold_row=gold_row,
        gold_value_set=gold_value_set,
        value_source=value_source,
        verification=rescored.get("debate") or {},
    )
    return rescored


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-jsonl", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--gold-value-set",
        choices=["v2", "original"],
        default="v2",
        help="Gold value set to use for rescoring stored debate outputs.",
    )
    args = parser.parse_args()

    records = [
        rescore_record(record, gold_value_set=args.gold_value_set)
        for record in load_jsonl(repo_path(args.records_jsonl))
    ]
    out_dir = repo_path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = out_dir / "claim_debate_records.jsonl"
    with jsonl_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    write_json(
        out_dir / "run_summary.json",
        {
            "source_jsonl": repo_path(args.records_jsonl),
            "gold_value_set": args.gold_value_set,
            "records": len(records),
            "jsonl": jsonl_path,
            "llm_calls": 0,
        },
    )
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
