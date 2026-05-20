#!/usr/bin/env python3
"""Run small-context claim verification over selected recall-regression papers."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

try:
    from scripts.recall_audit.common import (
        DEFAULT_RESULTS_DIR,
        canonical_variant,
        display_variant,
        find_full_contexts,
        load_gold_rows,
        parse_int,
        query_pipeline_rows,
        read_csv_rows,
        repo_path,
        resolve_gene_db,
        write_csv_rows,
        write_json,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        DEFAULT_RESULTS_DIR,
        canonical_variant,
        display_variant,
        find_full_contexts,
        load_gold_rows,
        parse_int,
        query_pipeline_rows,
        read_csv_rows,
        repo_path,
        resolve_gene_db,
        write_csv_rows,
        write_json,
    )

from config.settings import get_settings
from pipeline.claim_verifier import (
    SUPPORTED_VERDICTS,
    VariantClaimVerifier,
    build_claim_card,
    normalize_verification,
)


FIELDS = ("total_carriers", "affected", "unaffected")
DB_FIELD_MAP = {
    "total_carriers": "total_carriers_observed",
    "affected": "affected_count",
    "unaffected": "unaffected_count",
}
GOLD_FIELD_MAP = {
    "total_carriers": "carriers",
    "affected": "affected",
    "unaffected": "unaffected",
}

DEFAULT_GENE_DISEASES = {
    "KCNH2": "Long QT syndrome type 2, Short QT syndrome",
    "KCNQ1": "Long QT syndrome type 1, Jervell and Lange-Nielsen syndrome",
    "SCN5A": "Long QT syndrome type 3, Brugada syndrome, cardiac conduction disease",
    "RYR2": "Catecholaminergic polymorphic ventricular tachycardia",
}


def row_to_variant(row: dict[str, Any]) -> dict[str, Any]:
    total = parse_int(row.get("total_carriers_observed"))
    return {
        "gene_symbol": row.get("gene_symbol"),
        "cdna_notation": row.get("cdna_notation"),
        "protein_notation": row.get("protein_notation"),
        "genomic_position": row.get("genomic_position"),
        "clinical_significance": row.get("clinical_significance"),
        "patients": {"count": total},
        "penetrance_data": {
            "total_carriers_observed": total,
            "affected_count": parse_int(row.get("affected_count")),
            "unaffected_count": parse_int(row.get("unaffected_count")),
            "uncertain_count": parse_int(row.get("uncertain_count")),
        },
        "source_location": row.get("source_location"),
        "additional_notes": row.get("additional_notes"),
    }


def gold_index(
    gene: str, pmid: str, gold_dir: str | None = None
) -> dict[str, list[dict]]:
    rows = load_gold_rows(gene, pmid, gold_dir)
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[canonical_variant(row.get("variant"))].append(row)
    return dict(grouped)


def trusted_value(
    field: str, db_value: int | None, verification: dict[str, Any]
) -> int | None:
    verdict = (verification.get("field_verdicts") or {}).get(field, "ambiguous")
    if verdict not in SUPPORTED_VERDICTS:
        return None
    corrected = (verification.get("corrected_values") or {}).get(field)
    return parse_int(corrected) if corrected is not None else db_value


def evaluate_field(
    *,
    field: str,
    db_value: int | None,
    gold_value: int | None,
    verification: dict[str, Any],
) -> dict[str, Any]:
    trusted = trusted_value(field, db_value, verification)
    before_wrong = (
        db_value is not None and gold_value is not None and db_value != gold_value
    )
    before_correct = (
        db_value is not None and gold_value is not None and db_value == gold_value
    )
    after_wrong = (
        trusted is not None and gold_value is not None and trusted != gold_value
    )
    after_correct = (
        trusted is not None and gold_value is not None and trusted == gold_value
    )
    return {
        "field": field,
        "db_value": db_value,
        "gold_value": gold_value,
        "trusted_value": trusted,
        "before_wrong": before_wrong,
        "before_correct": before_correct,
        "after_wrong": after_wrong,
        "after_correct": after_correct,
        "flagged_not_autopopulated": db_value is not None and trusted is None,
        "verdict": (verification.get("field_verdicts") or {}).get(field),
    }


def candidate_rows(
    *,
    db_rows: list[dict[str, Any]],
    gold: dict[str, dict],
    max_cards: int,
) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    for row in db_rows:
        key = canonical_variant(display_variant(row))
        if not key or key in grouped:
            continue
        grouped[key] = row

    intersect = [row for key, row in grouped.items() if key in gold]
    extras = [row for key, row in grouped.items() if key not in gold]
    return [*intersect, *extras][:max_cards]


def values_from_row(row: dict[str, Any]) -> dict[str, int | None]:
    return {field: parse_int(row.get(DB_FIELD_MAP[field])) for field in FIELDS}


def best_gold_row(
    gold_rows: list[dict[str, Any]],
    values: dict[str, int | None],
) -> dict[str, Any]:
    if not gold_rows:
        return {}

    def score(row: dict[str, Any]) -> tuple[int, int]:
        exact = 0
        distance = 0
        for field, gold_field in GOLD_FIELD_MAP.items():
            value = values.get(field)
            gold_value = parse_int(row.get(gold_field))
            if value is None or gold_value is None:
                continue
            if value == gold_value:
                exact += 1
            else:
                distance += abs(value - gold_value)
        return exact, -distance

    return max(gold_rows, key=score)


def run_for_case(
    *,
    case: dict[str, str],
    model: str,
    run_dir: Path,
    results_dir: Path,
    disease: str | None,
    max_cards_per_paper: int,
    dry_run: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    gene = case["gene"].upper()
    pmid = case["pmid"]
    contexts = find_full_contexts(results_dir, gene, pmid, global_search=True)
    if not contexts:
        return [], [
            {
                "model": model,
                "gene": gene,
                "pmid": pmid,
                "error": "no_full_context",
            }
        ]
    source_text = contexts[0].read_text(encoding="utf-8", errors="ignore")
    db_rows = query_pipeline_rows(resolve_gene_db(gene, run_dir=run_dir), pmid)
    gold = gold_index(gene, pmid)
    rows = candidate_rows(
        db_rows=db_rows,
        gold=gold,
        max_cards=max_cards_per_paper,
    )
    verifier = None if dry_run else VariantClaimVerifier(model=model)

    records: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for row in rows:
        variant = row_to_variant(row)
        display = display_variant(row)
        card = build_claim_card(
            source_text=source_text,
            gene=gene,
            disease=disease or DEFAULT_GENE_DISEASES.get(gene),
            pmid=pmid,
            title=row.get("title") or case.get("title"),
            variant=variant,
            max_evidence_chars=4_000,
        )
        if card is None:
            continue
        if dry_run:
            verification = {
                "verdict": "not_run",
                "field_verdicts": {
                    "variant": "ambiguous",
                    "total_carriers": "ambiguous",
                    "affected": "ambiguous",
                    "unaffected": "ambiguous",
                },
                "corrected_values": {
                    "total_carriers": None,
                    "affected": None,
                    "unaffected": None,
                },
                "reason": "dry_run",
                "evidence_quote": "",
            }
        else:
            try:
                verification = verifier.verify(card) if verifier else {}
            except Exception as exc:  # noqa: BLE001
                errors.append(
                    {
                        "model": model,
                        "gene": gene,
                        "pmid": pmid,
                        "variant": display,
                        "error": str(exc),
                    }
                )
                continue
            verification = normalize_verification(verification)

        gold_rows = gold.get(canonical_variant(display), [])
        gold_row = best_gold_row(gold_rows, values_from_row(row))
        field_evals = [
            evaluate_field(
                field=field,
                db_value=parse_int(row.get(DB_FIELD_MAP[field])),
                gold_value=parse_int(gold_row.get(GOLD_FIELD_MAP[field])),
                verification=verification,
            )
            for field in FIELDS
        ]
        records.append(
            {
                "model": model,
                "gene": gene,
                "pmid": pmid,
                "failure_class": case.get("failure_class"),
                "variant": display,
                "in_gold": bool(gold_rows),
                "gold_candidate_rows": len(gold_rows),
                "matched_gold_row": gold_row,
                "source_context": str(contexts[0]),
                "card": card.__dict__,
                "verification": verification,
                "field_evaluations": field_evals,
            }
        )
    return records, errors


def summarize(
    records: list[dict[str, Any]], errors: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for record in records:
        key = (record["model"], record["failure_class"] or "")
        grouped[key]["cards"] += 1
        if record["in_gold"]:
            grouped[key]["cards_in_gold"] += 1
        for item in record["field_evaluations"]:
            grouped[key]["fields"] += 1
            if item["before_correct"]:
                grouped[key]["before_correct"] += 1
            if item["before_wrong"]:
                grouped[key]["before_wrong"] += 1
            if item["after_correct"]:
                grouped[key]["after_trusted_correct"] += 1
            if item["after_wrong"]:
                grouped[key]["after_trusted_wrong"] += 1
            if item["flagged_not_autopopulated"]:
                grouped[key]["flagged_not_autopopulated"] += 1
            grouped[key][f"verdict_{item['verdict']}"] += 1

    for error in errors:
        grouped[(error["model"], "ERROR")]["errors"] += 1

    rows = []
    for (model, failure_class), counts in sorted(grouped.items()):
        row: dict[str, Any] = {"model": model, "failure_class": failure_class}
        row.update(counts)
        before_wrong = counts.get("before_wrong", 0)
        after_wrong = counts.get("after_trusted_wrong", 0)
        row["wrong_autopopulated_delta"] = after_wrong - before_wrong
        rows.append(row)
    return rows


def summarize_by_paper(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str], Counter[str]] = defaultdict(Counter)
    for record in records:
        key = (
            record["model"],
            record["gene"],
            record["pmid"],
            record["failure_class"] or "",
        )
        grouped[key]["cards"] += 1
        if record["in_gold"]:
            grouped[key]["cards_in_gold"] += 1
        for item in record["field_evaluations"]:
            grouped[key]["fields"] += 1
            grouped[key]["before_correct"] += bool(item["before_correct"])
            grouped[key]["before_wrong"] += bool(item["before_wrong"])
            grouped[key]["after_trusted_correct"] += bool(item["after_correct"])
            grouped[key]["after_trusted_wrong"] += bool(item["after_wrong"])
            grouped[key]["flagged_not_autopopulated"] += bool(
                item["flagged_not_autopopulated"]
            )

    rows: list[dict[str, Any]] = []
    for (model, gene, pmid, failure_class), counts in sorted(grouped.items()):
        row: dict[str, Any] = {
            "model": model,
            "gene": gene,
            "pmid": pmid,
            "failure_class": failure_class,
        }
        row.update(counts)
        row["wrong_autopopulated_delta"] = counts.get(
            "after_trusted_wrong", 0
        ) - counts.get("before_wrong", 0)
        rows.append(row)
    return rows


def select_cases(
    cases: list[dict[str, str]], *, limit: int, balanced: bool
) -> list[dict[str, str]]:
    if not balanced:
        return cases[:limit]
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for case in cases:
        by_gene[case.get("gene", "").upper()].append(case)
    selected: list[dict[str, str]] = []
    genes = sorted(gene for gene in by_gene if gene)
    while len(selected) < limit and any(by_gene.values()):
        for gene in genes:
            if not by_gene[gene]:
                continue
            selected.append(by_gene[gene].pop(0))
            if len(selected) >= limit:
                break
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases-csv", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--model", action="append", dest="models")
    parser.add_argument("--pmid", action="append", dest="pmids")
    parser.add_argument(
        "--disease",
        help="Optional disease term applied to all claim cards.",
    )
    parser.add_argument("--limit-papers", type=int, default=12)
    parser.add_argument("--max-cards-per-paper", type=int, default=3)
    parser.add_argument("--failure-class", action="append", dest="failure_classes")
    parser.add_argument(
        "--balanced",
        action="store_true",
        help="Select cases round-robin by gene instead of taking the CSV prefix.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cases = read_csv_rows(repo_path(args.cases_csv))
    if args.failure_classes:
        wanted = set(args.failure_classes)
        cases = [case for case in cases if case.get("failure_class") in wanted]
    if args.pmids:
        wanted_pmids = set(args.pmids)
        cases = [case for case in cases if case.get("pmid") in wanted_pmids]
    cases = select_cases(cases, limit=args.limit_papers, balanced=args.balanced)

    models = args.models or get_settings().get_tier3_adjudicator_models()
    out_dir = repo_path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for model in models:
        for case in cases:
            case_records, case_errors = run_for_case(
                case=case,
                model=model,
                run_dir=repo_path(args.run_dir),
                results_dir=repo_path(args.results_dir),
                disease=args.disease,
                max_cards_per_paper=args.max_cards_per_paper,
                dry_run=args.dry_run,
            )
            records.extend(case_records)
            errors.extend(case_errors)
            print(
                f"{model} {case.get('gene')} {case.get('pmid')}: "
                f"{len(case_records)} cards, {len(case_errors)} errors",
                file=sys.stderr,
            )

    jsonl_path = out_dir / "claim_verification_records.jsonl"
    with jsonl_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    write_json(out_dir / "claim_verification_errors.json", errors)

    summary_rows = summarize(records, errors)
    fieldnames = sorted({key for row in summary_rows for key in row})
    write_csv_rows(summary_rows, fieldnames, out_dir / "claim_verification_summary.csv")
    paper_summary_rows = summarize_by_paper(records)
    paper_fieldnames = sorted({key for row in paper_summary_rows for key in row})
    write_csv_rows(
        paper_summary_rows,
        paper_fieldnames,
        out_dir / "claim_verification_by_paper.csv",
    )
    write_json(
        out_dir / "run_summary.json",
        {
            "models": models,
            "cases": len(cases),
            "records": len(records),
            "errors": len(errors),
            "jsonl": jsonl_path,
            "summary_csv": out_dir / "claim_verification_summary.csv",
            "paper_summary_csv": out_dir / "claim_verification_by_paper.csv",
        },
    )
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
