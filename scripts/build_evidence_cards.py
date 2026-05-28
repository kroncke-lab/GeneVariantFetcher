#!/usr/bin/env python3
"""Build evidence cards for high-uncertainty count rows in a run (#6).

Walks an extractions directory and produces a JSONL file of EvidenceCards
for downstream verification — either by an LLM verifier (future) or by a
human reviewer. Each card identifies one (PMID × variant × count field) row
that warrants a second look, with the source text excerpt and the triggers
that flagged it.

Modes (dual; same artifact serves both):
  * No-gold (default): triggers fire on outlier-guard flags, classifier
    flags, and ``unknown``/missing provenance with large values.
  * Gold: pass ``--gold <gene_recall_input.csv>`` to additionally emit
    cards for any matched (PMID × variant) row where |gold - extracted| ≥
    ``--gold-error-threshold`` on any count field.

The heuristic verifier is applied to each card so the verdict column gives
a fast first-pass triage. Verdicts: ``reject`` / ``withhold`` (the
heuristic never auto-confirms; that requires positive per-variant
evidence which a downstream verifier must supply).

Typical usage:

    .venv/bin/python scripts/build_evidence_cards.py \\
        --extractions <run>/extractions \\
        --pmc-dir <run>/pmc_fulltext \\
        --output /tmp/evidence_cards.jsonl

With gold (for audit on validation surface only — production paths must
not require gold):

    .venv/bin/python scripts/build_evidence_cards.py \\
        --extractions <run>/extractions \\
        --pmc-dir <run>/pmc_fulltext \\
        --gold gene_variant_fetcher_gold_standard/normalized/KCNQ1_recall_input.csv \\
        --output /tmp/evidence_cards_kcnq1.jsonl
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.evidence_card import (  # noqa: E402
    GOLD_ERROR_THRESHOLD_DEFAULT,
    LARGE_VALUE_THRESHOLD_DEFAULT,
    build_evidence_cards,
)

PMID_FROM_FILENAME = re.compile(r"_PMID_(\d+)\.json$")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--extractions",
        required=True,
        type=Path,
        help="Extractions directory (<GENE>_PMID_<PMID>.json files).",
    )
    p.add_argument(
        "--pmc-dir",
        type=Path,
        default=None,
        help=(
            "pmc_fulltext directory. If supplied, the source excerpt around the "
            "extracted value is included in each card."
        ),
    )
    p.add_argument(
        "--gold",
        type=Path,
        default=None,
        help=(
            "Optional gene_recall_input.csv. Triggers gold-mode (|gold - "
            "extracted| ≥ threshold) in addition to the no-gold triggers."
        ),
    )
    p.add_argument(
        "--gold-error-threshold",
        type=int,
        default=GOLD_ERROR_THRESHOLD_DEFAULT,
        help="Minimum |gold - extracted| for a gold-mode trigger.",
    )
    p.add_argument(
        "--large-value-threshold",
        type=int,
        default=LARGE_VALUE_THRESHOLD_DEFAULT,
        help="Minimum extracted value for the unknown-provenance trigger.",
    )
    p.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Where to write evidence_cards.jsonl.",
    )
    p.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional summary JSON path (totals by trigger and verdict).",
    )
    return p


def _load_gold(gold_path: Path) -> dict[str, dict[str, dict[str, int]]]:
    """Return {pmid: {variant_norm: {carriers/affected/unaffected: int}}}.

    Variant normalization here matches the lookup in build_evidence_cards:
    the LLM's protein_notation is used as the key. The recall_input.csv
    already provides a `variant` column in canonical form.
    """
    gold: dict[str, dict[str, dict[str, int]]] = {}
    with gold_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pmid = str(row.get("pmid") or "").strip()
            variant = (row.get("variant") or "").strip()
            if not pmid or not variant:
                continue
            d = gold.setdefault(pmid, {}).setdefault(variant, {})
            for col, key in (
                ("carriers", "carriers"),
                ("affected", "affected"),
                ("unaffected", "unaffected"),
            ):
                v = row.get(col)
                if v is None or v == "":
                    continue
                try:
                    d[key] = int(float(v))
                except (TypeError, ValueError):
                    continue
    return gold


def _load_source_text(pmc_dir: Path | None, pmid: str) -> str:
    if pmc_dir is None:
        return ""
    # Prefer FULL_CONTEXT.md, then CLEANED.md
    for suffix in ("_FULL_CONTEXT.md", "_CLEANED.md"):
        p = pmc_dir / f"{pmid}{suffix}"
        if p.exists():
            try:
                return p.read_text(encoding="utf-8", errors="replace")
            except Exception:  # noqa: BLE001
                return ""
    return ""


def main() -> int:
    args = build_parser().parse_args()
    extractions = args.extractions.expanduser().resolve()
    if not extractions.is_dir():
        sys.exit(f"extractions dir not found: {extractions}")
    pmc_dir = args.pmc_dir.expanduser().resolve() if args.pmc_dir else None
    gold = _load_gold(args.gold.expanduser().resolve()) if args.gold else None

    args.output.parent.mkdir(parents=True, exist_ok=True)
    totals = {
        "pmids_processed": 0,
        "cards_total": 0,
        "by_trigger": {},
        "by_verdict": {},
    }
    with args.output.open("w", encoding="utf-8") as out:
        for json_path in sorted(extractions.glob("*.json")):
            m = PMID_FROM_FILENAME.search(json_path.name)
            if not m:
                continue
            pmid = m.group(1)
            try:
                payload = json.loads(json_path.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                continue
            variants = payload.get("variants")
            if not isinstance(variants, list):
                continue
            totals["pmids_processed"] += 1
            source_text = _load_source_text(pmc_dir, pmid)
            gold_for_pmid = (gold or {}).get(pmid)
            cards = build_evidence_cards(
                pmid=pmid,
                variants=variants,
                source_text=source_text,
                gold_counts_by_variant=gold_for_pmid if gold is not None else None,
                gold_error_threshold=args.gold_error_threshold,
                large_value_threshold=args.large_value_threshold,
            )
            for card in cards:
                card.run_heuristic_verdict()
                out.write(json.dumps(card.as_dict()) + "\n")
                totals["cards_total"] += 1
                for trig in card.triggers:
                    totals["by_trigger"][trig] = totals["by_trigger"].get(trig, 0) + 1
                if card.verdict:
                    totals["by_verdict"][card.verdict] = (
                        totals["by_verdict"].get(card.verdict, 0) + 1
                    )

    print(json.dumps(totals, indent=2))
    if args.summary_out:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(
            json.dumps(
                {
                    "extractions_dir": str(extractions),
                    "pmc_dir": str(pmc_dir) if pmc_dir else None,
                    "gold": str(args.gold) if args.gold else None,
                    "totals": totals,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
