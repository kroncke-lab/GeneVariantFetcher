#!/usr/bin/env python3
"""Run second-model debate over existing claim-verification records."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from config.settings import get_settings
    from scripts.recall_audit.common import repo_path, write_csv_rows, write_json
    from scripts.recall_audit.run_claim_verification_pilot import FIELDS, evaluate_field
except ModuleNotFoundError:  # pragma: no cover
    from config.settings import get_settings
    from common import repo_path, write_csv_rows, write_json
    from run_claim_verification_pilot import FIELDS, evaluate_field

from pipeline.claim_verifier import (
    FIELD_NAMES,
    VariantClaimCard,
    normalize_verification,
)
from utils.llm_utils import BaseLLMCaller, clamp_max_tokens

AGREEMENTS = {"agree", "disagree", "partial", "abstain"}


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def load_queue_keys(path: Path) -> set[tuple[str, str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (
                str(row.get("gene") or ""),
                str(row.get("pmid") or ""),
                str(row.get("variant") or ""),
            )
            for row in csv.DictReader(handle)
        }


def card_from_record(record: dict[str, Any]) -> VariantClaimCard:
    card = record["card"]
    return VariantClaimCard(
        gene=card.get("gene") or record.get("gene") or "",
        disease=card.get("disease"),
        pmid=str(card.get("pmid") or record.get("pmid") or ""),
        title=card.get("title"),
        variant=card.get("variant") or record.get("variant") or "",
        extracted=card.get("extracted") or {},
        evidence=card.get("evidence") or "",
        source_location=card.get("source_location"),
    )


def build_debate_prompt(
    *,
    card: VariantClaimCard,
    baseline_model: str,
    baseline_verification: dict[str, Any],
) -> str:
    baseline = {
        "model": baseline_model,
        "verification": baseline_verification,
    }
    return f"""You are adjudicating another model's biomedical variant/count verification.

You get the same local evidence card plus the baseline model's structured
verdict. Do not defer to the baseline model. First read the evidence yourself,
then decide whether the baseline is supported.

Return your own final field verdicts and corrected count values. Also say
whether you agree with the baseline overall and per field:
- agree: baseline field verdict and value are supported.
- disagree: baseline field verdict or value is wrong.
- partial: baseline is directionally useful but incomplete or overconfident.
- abstain: evidence is insufficient to judge the baseline.

Important rules:
- Do not turn a baseline model's rationale into evidence. Only the evidence
  lines in the claim card can support a value.
- A lower symptom/event subset is not affected_count when the target disease
  patient/proband count is larger.
- If the relevant table/supplement is absent from the evidence, mark
  count fields unsupported or source_missing and leave corrected values null.
- Prefer abstention over autopopulating a count from unlabeled columns,
  study-wide totals, family totals, allele frequencies, or prediction scores.

Claim card:
{card.to_prompt_json()}

Baseline verification:
{json.dumps(baseline, ensure_ascii=False, indent=2)}

Return strict JSON only:
{{
  "agreement": "agree|disagree|partial|abstain",
  "field_agreement": {{
    "variant": "agree|disagree|partial|abstain",
    "total_carriers": "agree|disagree|partial|abstain",
    "affected": "agree|disagree|partial|abstain",
    "unaffected": "agree|disagree|partial|abstain"
  }},
  "verdict": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
  "field_verdicts": {{
    "variant": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "total_carriers": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "affected": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "unaffected": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing"
  }},
  "corrected_values": {{
    "total_carriers": "integer or null",
    "affected": "integer or null",
    "unaffected": "integer or null"
  }},
  "reason": "brief evidence-based explanation of agreement/disagreement",
  "evidence_quote": "short quote or line reference supporting your decision"
}}
"""


def normalize_debate(raw: dict[str, Any], card: VariantClaimCard) -> dict[str, Any]:
    normalized = normalize_verification(raw, card=card)
    agreement = str(raw.get("agreement") or "abstain").strip().lower()
    if agreement not in AGREEMENTS:
        agreement = "abstain"
    raw_field_agreement = raw.get("field_agreement") or {}
    field_agreement = {}
    for field in FIELD_NAMES:
        value = str(raw_field_agreement.get(field) or agreement).strip().lower()
        field_agreement[field] = value if value in AGREEMENTS else "abstain"
    normalized["agreement"] = agreement
    normalized["field_agreement"] = field_agreement
    return normalized


class ClaimDebateVerifier(BaseLLMCaller):
    def __init__(
        self,
        model: str,
        temperature: float = 0.0,
        max_tokens: int = 2500,
    ):
        super().__init__(
            model=model,
            temperature=temperature,
            max_tokens=clamp_max_tokens(model, max_tokens),
        )

    def debate(
        self,
        *,
        card: VariantClaimCard,
        baseline_model: str,
        baseline_verification: dict[str, Any],
    ) -> dict[str, Any]:
        raw = self.call_llm_json(
            build_debate_prompt(
                card=card,
                baseline_model=baseline_model,
                baseline_verification=baseline_verification,
            )
        )
        return normalize_debate(raw, card)


def evaluate_debate(
    record: dict[str, Any],
    debate: dict[str, Any],
) -> list[dict[str, Any]]:
    by_field = {item["field"]: item for item in record.get("field_evaluations", [])}
    items = []
    for field in FIELDS:
        baseline_item = by_field.get(field, {})
        items.append(
            evaluate_field(
                field=field,
                db_value=baseline_item.get("db_value"),
                gold_value=baseline_item.get("gold_value"),
                original_gold_value=baseline_item.get("original_gold_value"),
                gold_value_source=baseline_item.get("gold_value_source") or "original",
                verification=debate,
            )
        )
    return items


def summarize(records: list[dict[str, Any]], errors: list[dict[str, Any]]):
    grouped: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for record in records:
        key = (
            record["baseline_model"],
            record["debater_model"],
            record.get("failure_class") or "",
        )
        counts = grouped[key]
        counts["cards"] += 1
        counts[f"agreement_{record['debate'].get('agreement', 'abstain')}"] += 1
        for item in record.get("baseline_field_evaluations", []):
            if item.get("after_correct"):
                counts["baseline_trusted_correct"] += 1
            if item.get("after_wrong"):
                counts["baseline_trusted_wrong"] += 1
            if item.get("flagged_not_autopopulated"):
                counts["baseline_flagged"] += 1
        for item in record.get("debate_field_evaluations", []):
            counts["fields"] += 1
            if item.get("after_correct"):
                counts["debate_trusted_correct"] += 1
            if item.get("after_wrong"):
                counts["debate_trusted_wrong"] += 1
            if item.get("flagged_not_autopopulated"):
                counts["debate_flagged"] += 1

    for error in errors:
        grouped[
            (
                error.get("baseline_model", ""),
                error.get("debater_model", ""),
                "ERROR",
            )
        ]["errors"] += 1

    rows = []
    for (baseline_model, debater_model, failure_class), counts in sorted(
        grouped.items()
    ):
        row: dict[str, Any] = {
            "baseline_model": baseline_model,
            "debater_model": debater_model,
            "failure_class": failure_class,
        }
        row.update(counts)
        row["trusted_wrong_delta"] = counts.get("debate_trusted_wrong", 0) - counts.get(
            "baseline_trusted_wrong", 0
        )
        rows.append(row)
    return rows


def summarize_by_paper(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str, str], Counter[str]] = defaultdict(Counter)
    for record in records:
        key = (
            record["baseline_model"],
            record["debater_model"],
            record["gene"],
            record["pmid"],
            record.get("failure_class") or "",
        )
        counts = grouped[key]
        counts["cards"] += 1
        counts[f"agreement_{record['debate'].get('agreement', 'abstain')}"] += 1
        for item in record.get("baseline_field_evaluations", []):
            counts["baseline_trusted_correct"] += bool(item.get("after_correct"))
            counts["baseline_trusted_wrong"] += bool(item.get("after_wrong"))
            counts["baseline_flagged"] += bool(item.get("flagged_not_autopopulated"))
        for item in record.get("debate_field_evaluations", []):
            counts["debate_trusted_correct"] += bool(item.get("after_correct"))
            counts["debate_trusted_wrong"] += bool(item.get("after_wrong"))
            counts["debate_flagged"] += bool(item.get("flagged_not_autopopulated"))

    rows = []
    for (baseline_model, debater_model, gene, pmid, failure_class), counts in sorted(
        grouped.items()
    ):
        row: dict[str, Any] = {
            "baseline_model": baseline_model,
            "debater_model": debater_model,
            "gene": gene,
            "pmid": pmid,
            "failure_class": failure_class,
        }
        row.update(counts)
        row["trusted_wrong_delta"] = counts.get("debate_trusted_wrong", 0) - counts.get(
            "baseline_trusted_wrong", 0
        )
        rows.append(row)
    return rows


def select_records(
    records: list[dict[str, Any]],
    *,
    baseline_models: set[str] | None,
    pmids: set[str] | None,
    failure_classes: set[str] | None,
    queue_keys: set[tuple[str, str, str]] | None,
    limit_cards: int,
) -> list[dict[str, Any]]:
    selected = []
    for record in records:
        if baseline_models and record.get("model") not in baseline_models:
            continue
        if pmids and str(record.get("pmid")) not in pmids:
            continue
        if failure_classes and record.get("failure_class") not in failure_classes:
            continue
        if (
            queue_keys
            and (
                str(record.get("gene") or ""),
                str(record.get("pmid") or ""),
                str(record.get("variant") or ""),
            )
            not in queue_keys
        ):
            continue
        selected.append(record)
        if limit_cards and len(selected) >= limit_cards:
            break
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-records", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--debater-model", action="append", dest="debater_models")
    parser.add_argument(
        "--queue-csv",
        help=(
            "Optional escalation queue CSV from build_claim_debate_escalation_queue.py. "
            "When set, only matching gene/PMID/variant cards are debated."
        ),
    )
    parser.add_argument(
        "--final-adjudicator",
        action="store_true",
        help="Use FINAL_ADJUDICATOR_MODELS instead of early debate defaults.",
    )
    parser.add_argument(
        "--final-arbiter",
        action="store_true",
        help="Use FINAL_ARBITER_MODEL as the only debater model.",
    )
    parser.add_argument("--baseline-model", action="append", dest="baseline_models")
    parser.add_argument("--pmid", action="append", dest="pmids")
    parser.add_argument("--failure-class", action="append", dest="failure_classes")
    parser.add_argument("--limit-cards", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.final_adjudicator and args.final_arbiter:
        parser.error("--final-adjudicator and --final-arbiter are mutually exclusive")

    baseline_records = load_jsonl(repo_path(args.baseline_records))
    queue_keys = load_queue_keys(repo_path(args.queue_csv)) if args.queue_csv else None
    selected = select_records(
        baseline_records,
        baseline_models=set(args.baseline_models) if args.baseline_models else None,
        pmids=set(args.pmids) if args.pmids else None,
        failure_classes=set(args.failure_classes) if args.failure_classes else None,
        queue_keys=queue_keys,
        limit_cards=args.limit_cards,
    )

    settings = get_settings()
    if args.debater_models:
        debater_models = args.debater_models
    elif args.final_arbiter:
        debater_models = [settings.get_final_arbiter_model()]
    elif args.final_adjudicator:
        debater_models = settings.get_final_adjudicator_models()
    else:
        debater_models = settings.get_early_debate_models()
    out_dir = repo_path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    debate_records: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for debater_model in debater_models:
        verifier = None if args.dry_run else ClaimDebateVerifier(model=debater_model)
        for record in selected:
            card = card_from_record(record)
            baseline_verification = record.get("verification") or {}
            try:
                debate = (
                    normalize_debate(
                        {
                            "agreement": "abstain",
                            "field_agreement": {},
                            "verdict": "ambiguous",
                            "field_verdicts": {},
                            "corrected_values": {},
                            "reason": "dry_run",
                            "evidence_quote": "",
                        },
                        card,
                    )
                    if args.dry_run
                    else verifier.debate(
                        card=card,
                        baseline_model=record.get("model") or "",
                        baseline_verification=baseline_verification,
                    )
                )
            except Exception as exc:  # noqa: BLE001
                errors.append(
                    {
                        "baseline_model": record.get("model"),
                        "debater_model": debater_model,
                        "gene": record.get("gene"),
                        "pmid": record.get("pmid"),
                        "variant": record.get("variant"),
                        "error": str(exc),
                    }
                )
                continue
            debate_records.append(
                {
                    "baseline_model": record.get("model"),
                    "debater_model": debater_model,
                    "gene": record.get("gene"),
                    "pmid": record.get("pmid"),
                    "failure_class": record.get("failure_class"),
                    "variant": record.get("variant"),
                    "source_context": record.get("source_context"),
                    "card": record.get("card"),
                    "baseline_verification": baseline_verification,
                    "baseline_field_evaluations": record.get("field_evaluations", []),
                    "debate": debate,
                    "debate_field_evaluations": evaluate_debate(record, debate),
                }
            )
            print(
                f"{debater_model} over {record.get('model')} "
                f"{record.get('gene')} {record.get('pmid')} "
                f"{record.get('variant')}: {debate.get('agreement')}",
                file=sys.stderr,
            )

    jsonl_path = out_dir / "claim_debate_records.jsonl"
    with jsonl_path.open("w", encoding="utf-8") as handle:
        for record in debate_records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    write_json(out_dir / "claim_debate_errors.json", errors)

    summary_rows = summarize(debate_records, errors)
    summary_fieldnames = sorted({key for row in summary_rows for key in row})
    write_csv_rows(
        summary_rows, summary_fieldnames, out_dir / "claim_debate_summary.csv"
    )
    paper_rows = summarize_by_paper(debate_records)
    paper_fieldnames = sorted({key for row in paper_rows for key in row})
    write_csv_rows(
        paper_rows,
        paper_fieldnames,
        out_dir / "claim_debate_by_paper.csv",
    )
    write_json(
        out_dir / "run_summary.json",
        {
            "baseline_records": str(repo_path(args.baseline_records)),
            "queue_csv": str(repo_path(args.queue_csv)) if args.queue_csv else None,
            "selected_records": len(selected),
            "debater_models": debater_models,
            "records": len(debate_records),
            "errors": len(errors),
            "jsonl": jsonl_path,
            "summary_csv": out_dir / "claim_debate_summary.csv",
            "paper_summary_csv": out_dir / "claim_debate_by_paper.csv",
        },
    )
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
