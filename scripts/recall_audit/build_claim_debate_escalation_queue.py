#!/usr/bin/env python3
"""Build a deduplicated escalation queue from claim-debate outputs.

Optionally include direct claim-verification records. Direct records are used to
catch high-risk false consensus cases where models agree on a trusted value but
the evidence pattern is known to be fragile, such as pedigree-derived counts or
aggregate cohort counts.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    from scripts.recall_audit.common import repo_path, write_csv_rows, write_json
except ModuleNotFoundError:  # pragma: no cover
    from common import repo_path, write_csv_rows, write_json

COUNT_FIELDS = ("total_carriers", "affected", "unaffected")
HIGH_RISK_PATTERNS = {
    "pedigree_family": ("family", "families", "pedigree", "sibling", "parent"),
    "segregation_haplotype": ("haplotype", "segregat", "cosegregat"),
    "zygosity_recessive": ("homozygous", "heterozygous", "recessive"),
    "aggregate_counts": (
        "a total of",
        "overall population",
        "all mutations",
        "five mutations",
        "six index",
        "aggregate",
    ),
    "cohort_vs_carrier": (
        "samples",
        "controls",
        "cases",
        "heterozygous",
        "sids",
        "cohort",
    ),
    "table_figure": ("table", "figure", "supplement", "fig."),
}


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def clip(value: Any, limit: int = 240) -> str:
    text = "" if value is None else str(value)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def field_names_with_disagreement(record: dict[str, Any]) -> list[str]:
    field_agreement = (record.get("debate") or {}).get("field_agreement") or {}
    return [
        field
        for field, agreement in field_agreement.items()
        if str(agreement).strip().lower() != "agree"
    ]


def is_escalation_candidate(record: dict[str, Any]) -> bool:
    agreement = str((record.get("debate") or {}).get("agreement") or "").lower()
    return agreement != "agree" or bool(field_names_with_disagreement(record))


def severity(record: dict[str, Any]) -> str:
    agreement = str((record.get("debate") or {}).get("agreement") or "").lower()
    fields = set(field_names_with_disagreement(record))
    if agreement == "disagree" or fields.intersection(COUNT_FIELDS):
        return "high"
    if agreement in {"partial", "abstain"} or fields:
        return "medium"
    return "low"


def high_risk_labels(record: dict[str, Any]) -> list[str]:
    card = record.get("card") or {}
    evidence = " ".join(
        str(value or "")
        for value in (
            card.get("evidence"),
            card.get("source_location"),
            record.get("source_context"),
            record.get("failure_class"),
        )
    ).lower()
    return [
        label
        for label, needles in HIGH_RISK_PATTERNS.items()
        if any(needle in evidence for needle in needles)
    ]


def trusted_field_values(record: dict[str, Any]) -> dict[str, int]:
    values = {}
    for item in record.get("field_evaluations", []):
        field = item.get("field")
        value = item.get("trusted_value")
        if field in COUNT_FIELDS and value is not None:
            values[str(field)] = int(value)
    return values


def count_eval_flags(items: list[dict[str, Any]], key: str) -> int:
    return sum(1 for item in items if item.get(key))


def _queue_item(
    grouped: dict[tuple[str, str, str], dict[str, Any]],
    record: dict[str, Any],
    *,
    escalation_model: str,
) -> dict[str, Any]:
    key = (
        str(record.get("gene") or ""),
        str(record.get("pmid") or ""),
        str(record.get("variant") or ""),
    )
    return grouped.setdefault(
        key,
        {
            "gene": key[0],
            "pmid": key[1],
            "variant": key[2],
            "failure_class": record.get("failure_class") or "",
            "severity": "low",
            "recommended_escalation_model": escalation_model,
            "directions": [],
            "agreement_labels": [],
            "disagreement_fields": set(),
            "baseline_trusted_wrong": 0,
            "debate_trusted_wrong": 0,
            "baseline_flagged": 0,
            "debate_flagged": 0,
            "reasons": [],
            "source_context": record.get("source_context") or "",
        },
    )


def add_high_risk_consensus(
    grouped: dict[tuple[str, str, str], dict[str, Any]],
    records: list[dict[str, Any]],
    *,
    escalation_model: str,
) -> None:
    by_card: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        by_card[
            (
                str(record.get("gene") or ""),
                str(record.get("pmid") or ""),
                str(record.get("variant") or ""),
            )
        ].append(record)

    for _key, card_records in by_card.items():
        labels = sorted({label for r in card_records for label in high_risk_labels(r)})
        if not labels:
            continue
        by_field: dict[str, dict[int, set[str]]] = defaultdict(lambda: defaultdict(set))
        for record in card_records:
            model = str(record.get("model") or "")
            for field, value in trusted_field_values(record).items():
                by_field[field][value].add(model)

        consensus_fields: list[str] = []
        for field, value_to_models in by_field.items():
            if len(value_to_models) != 1:
                continue
            models = next(iter(value_to_models.values()))
            if len(models) >= 2:
                consensus_fields.append(field)
        if not consensus_fields:
            continue

        representative = card_records[0]
        item = _queue_item(grouped, representative, escalation_model=escalation_model)
        item["severity"] = "high"
        item["directions"].append(
            "direct high-risk consensus: "
            + ",".join(sorted({str(r.get("model") or "") for r in card_records}))
        )
        item["agreement_labels"].append("high_risk_consensus")
        item["disagreement_fields"].update(consensus_fields)
        item["baseline_trusted_wrong"] += sum(
            count_eval_flags(r.get("field_evaluations", []), "after_wrong")
            for r in card_records
        )
        item["baseline_flagged"] += sum(
            count_eval_flags(
                r.get("field_evaluations", []), "flagged_not_autopopulated"
            )
            for r in card_records
        )
        item["reasons"].append(
            clip(
                "Trusted consensus on high-risk evidence pattern(s): "
                + ",".join(labels)
                + ". Escalate even when debate agrees."
            )
        )


def build_queue(
    records: list[dict[str, Any]],
    *,
    escalation_model: str,
    verification_records: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], dict[str, Any]] = {}
    severity_rank = {"low": 0, "medium": 1, "high": 2}

    for record in records:
        if not is_escalation_candidate(record):
            continue
        item = _queue_item(grouped, record, escalation_model=escalation_model)
        current_severity = severity(record)
        if severity_rank[current_severity] > severity_rank[item["severity"]]:
            item["severity"] = current_severity
        agreement = (record.get("debate") or {}).get("agreement") or "abstain"
        direction = (
            f"{record.get('baseline_model')} -> "
            f"{record.get('debater_model')}: {agreement}"
        )
        item["directions"].append(direction)
        item["agreement_labels"].append(str(agreement))
        item["disagreement_fields"].update(field_names_with_disagreement(record))
        item["baseline_trusted_wrong"] += count_eval_flags(
            record.get("baseline_field_evaluations", []), "after_wrong"
        )
        item["debate_trusted_wrong"] += count_eval_flags(
            record.get("debate_field_evaluations", []), "after_wrong"
        )
        item["baseline_flagged"] += count_eval_flags(
            record.get("baseline_field_evaluations", []), "flagged_not_autopopulated"
        )
        item["debate_flagged"] += count_eval_flags(
            record.get("debate_field_evaluations", []), "flagged_not_autopopulated"
        )
        reason = (record.get("debate") or {}).get("reason")
        if reason:
            item["reasons"].append(clip(reason))

    if verification_records:
        add_high_risk_consensus(
            grouped,
            verification_records,
            escalation_model=escalation_model,
        )

    rows = []
    for item in grouped.values():
        row = dict(item)
        row["directions"] = " | ".join(row["directions"])
        row["agreement_labels"] = ";".join(sorted(set(row["agreement_labels"])))
        row["disagreement_fields"] = ";".join(sorted(row["disagreement_fields"]))
        row["reasons"] = " | ".join(row["reasons"][:3])
        rows.append(row)
    return sorted(
        rows,
        key=lambda row: (
            -severity_rank[row["severity"]],
            -int(row.get("baseline_trusted_wrong") or 0),
            row["gene"],
            row["pmid"],
            row["variant"],
        ),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--debate-records", action="append", required=True)
    parser.add_argument(
        "--verification-records",
        action="append",
        help=(
            "Optional direct claim-verification JSONL used to add high-risk "
            "trusted-consensus rows."
        ),
    )
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-json")
    parser.add_argument(
        "--escalation-model",
        default="anthropic/claude-opus-4-7",
        help="Model to recommend for the next adjudication tier.",
    )
    args = parser.parse_args()

    records = []
    for path in args.debate_records:
        records.extend(load_jsonl(repo_path(path)))
    verification_records = []
    for path in args.verification_records or []:
        verification_records.extend(load_jsonl(repo_path(path)))
    rows = build_queue(
        records,
        escalation_model=args.escalation_model,
        verification_records=verification_records,
    )
    fieldnames = [
        "severity",
        "recommended_escalation_model",
        "gene",
        "pmid",
        "variant",
        "failure_class",
        "disagreement_fields",
        "agreement_labels",
        "directions",
        "baseline_trusted_wrong",
        "debate_trusted_wrong",
        "baseline_flagged",
        "debate_flagged",
        "source_context",
        "reasons",
    ]
    write_csv_rows(rows, fieldnames, repo_path(args.out_csv))
    if args.out_json:
        write_json(repo_path(args.out_json), rows)
    print(repo_path(args.out_csv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
