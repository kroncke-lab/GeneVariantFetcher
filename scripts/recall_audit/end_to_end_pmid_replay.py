#!/usr/bin/env python3
"""Replay extraction for one PMID and diff against DB plus gold rows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from scripts.recall_audit.common import (
        DEFAULT_RESULTS_DIR,
        canonical_variant,
        display_variant,
        find_full_contexts,
        load_gold_rows,
        query_pipeline_rows,
        repo_path,
        resolve_gene_db,
        write_json,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        DEFAULT_RESULTS_DIR,
        canonical_variant,
        display_variant,
        find_full_contexts,
        load_gold_rows,
        query_pipeline_rows,
        repo_path,
        resolve_gene_db,
        write_json,
    )


class RecordingExpertExtractor:
    """Lazy wrapper so importing this script does not import LLM dependencies."""

    def __init__(self, *args, **kwargs):
        from pipeline.extraction import ExpertExtractor

        class _Recorder(ExpertExtractor):
            last_prompt: str | None = None
            last_response: str | None = None

            def call_llm_json_with_status(self, prompt, *call_args, **call_kwargs):
                self.last_prompt = prompt
                data, was_truncated, raw_text = super().call_llm_json_with_status(
                    prompt, *call_args, **call_kwargs
                )
                self.last_response = raw_text
                return data, was_truncated, raw_text

        self.impl = _Recorder(*args, **kwargs)

    def extract(self, paper):
        return self.impl.extract(paper)

    @property
    def last_prompt(self):
        return self.impl.last_prompt

    @property
    def last_response(self):
        return self.impl.last_response


def write_dict_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--pmid", required=True)
    parser.add_argument(
        "--out", required=True, help="Output directory for replay artifacts"
    )
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--run-dir", help="Validation run directory containing dbs/")
    parser.add_argument("--db", help="Explicit SQLite DB path")
    parser.add_argument("--gold-dir", help="Gold-standard normalized directory")
    parser.add_argument("--model", help="Single LiteLLM model to use for replay")
    parser.add_argument("--max-tokens", type=int, help="Override extractor max_tokens")
    args = parser.parse_args()

    gene = args.gene.upper()
    out_dir = repo_path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    contexts = find_full_contexts(repo_path(args.results_dir), gene, args.pmid)
    if not contexts:
        raise SystemExit(f"No FULL_CONTEXT.md found for {gene} PMID {args.pmid}")
    context_path = contexts[0]
    full_text = context_path.read_text(encoding="utf-8", errors="ignore")

    gold_rows = load_gold_rows(gene, args.pmid, args.gold_dir)
    db_rows = query_pipeline_rows(
        resolve_gene_db(gene, args.db, args.run_dir), args.pmid
    )
    write_dict_csv(out_dir / "gold_rows.csv", gold_rows)
    write_dict_csv(out_dir / "db_snapshot_rows.csv", db_rows)

    from utils.models import Paper

    extractor_kwargs = {
        "models": [args.model] if args.model else None,
        "fulltext_dir": str(context_path.parent),
        "tier_threshold": 0,
    }
    if args.max_tokens:
        extractor_kwargs["max_tokens"] = args.max_tokens
    extractor = RecordingExpertExtractor(**extractor_kwargs)

    paper = Paper(
        pmid=str(args.pmid),
        title=f"Replay PMID {args.pmid}",
        full_text=full_text,
        gene_symbol=gene,
    )

    try:
        result = extractor.extract(paper)
    except Exception as exc:
        (out_dir / "error.txt").write_text(
            f"{type(exc).__name__}: {exc}\n", encoding="utf-8"
        )
        raise

    payload = result.model_dump() if hasattr(result, "model_dump") else result.dict()
    write_json(out_dir / "extraction_result.json", payload)
    if extractor.last_prompt is not None:
        (out_dir / "llm_prompt.txt").write_text(extractor.last_prompt, encoding="utf-8")
    if extractor.last_response is not None:
        (out_dir / "llm_response.txt").write_text(
            extractor.last_response, encoding="utf-8"
        )

    new_variants = (
        payload.get("extracted_data", {}).get("variants", [])
        if payload.get("success")
        else []
    )
    new_keys = {canonical_variant(display_variant(row)) for row in new_variants}
    db_keys = {canonical_variant(display_variant(row)) for row in db_rows}
    gold_keys = {canonical_variant(row.get("variant")) for row in gold_rows}
    new_keys.discard("")
    db_keys.discard("")
    gold_keys.discard("")

    diff_rows = [
        {
            "variant": variant,
            "in_replay": variant in new_keys,
            "in_db_snapshot": variant in db_keys,
            "in_gold": variant in gold_keys,
        }
        for variant in sorted(new_keys | db_keys | gold_keys)
    ]
    write_dict_csv(out_dir / "variant_diff.csv", diff_rows)
    write_json(
        out_dir / "diff_summary.json",
        {
            "gene": gene,
            "pmid": str(args.pmid),
            "context_path": context_path,
            "context_bytes": context_path.stat().st_size,
            "replay_variant_count": len(new_keys),
            "db_variant_count": len(db_keys),
            "gold_variant_count": len(gold_keys),
            "replay_only": sorted(new_keys - db_keys - gold_keys),
            "db_only": sorted(db_keys - new_keys),
            "gold_missing_from_replay": sorted(gold_keys - new_keys),
        },
    )
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
