#!/usr/bin/env python
"""Convert production gvf-run SQLite output into codex_paper_eval predictions.json.

The eval harness scores `papers[].variants[]` with {variant, carriers, affected,
unaffected}. Production stores the same facts across three tables, so this is a
straight projection -- no re-interpretation:

  variants.protein_notation + variants.cdna_notation  -> variants[].variant
  penetrance_data.total_carriers_observed             -> carriers
  penetrance_data.affected_count                      -> affected
  penetrance_data.unaffected_count                    -> unaffected
  variant_papers.source_location / key_quotes         -> source_location / evidence

Both notations go into one prediction string because run_eval.matches() expands
each side through variant_candidates(); handing it everything production stored
is the neutral choice (it neither hides nor invents a notation).

--exclude-layers writes the paper-derived-only view (drops ClinVar/PubTator
linkage rows, which are database lookups rather than readings of the paper) so
the comparison against a read-the-paper protocol is legible.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
if str(REPO) not in sys.path:  # run_eval + utils.variant_normalizer live under the repo
    sys.path.insert(0, str(REPO))

TOOL_RATIONALE = (
    "Production gvf-run strategy: deterministic regex/table sweep over the "
    "untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a "
    "60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk "
    "papers, then recovery layers. It reads the same corpus/ markdown the "
    "single-model runs read; --no-source-recovery kept the source identical. "
    "Mapped to the harness 'text' route because no single harness route "
    "describes a multi-source pipeline."
)


def find_db(root: Path, gene: str) -> Path:
    """Newest final <GENE>.db under root/<GENE>/**.

    Excludes the `.before_layers.db` snapshot gvf-run leaves next to the real
    database -- rglob("<GENE>*.db") matches it too, and picking it would
    silently score a pre-recovery-layer run.
    """
    candidates = [
        p
        for p in (root / gene).rglob(f"{gene}*.db")
        if "before_layers" not in p.name and not p.name.endswith(".before_layers.db")
    ]
    if not candidates:
        raise SystemExit(f"no final {gene}*.db found under {root / gene}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def layer_rank(source_layer: str | None) -> int:
    """Lower is more authoritative as the count-bearing record for a variant."""
    layers = {t.strip() for t in (source_layer or "").split(",") if t.strip()}
    if layers & {"llm_table", "llm_text", "regex_table", "regex_text", "figure"}:
        return 0  # read out of the paper
    return 1  # clinvar / pubtator database linkage


def merge_same_variant(rows: list[dict], gene: str) -> list[dict]:
    """Collapse rows that denote the same variant in different notations.

    Production can store one variant several times -- e.g. `p.Leu552Ser
    c.1655T>C` from llm_text with counts plus `p.L552S` from pubtator without.
    Emitting both would charge production a false positive for a variant it got
    right, so rows are merged using the harness's own matches(), and the
    surviving record keeps the counts from the most authoritative layer.
    """
    from benchmarks.codex_paper_eval.run_eval import matches

    kept: list[dict] = []
    for row in sorted(
        rows, key=lambda r: (layer_rank(r["source_layer"]), r["variant"])
    ):
        for existing in kept:
            if matches(row["variant"], existing["variant"], gene):
                for field in ("carriers", "affected", "unaffected"):
                    if existing[field] is None and row[field] is not None:
                        existing[field] = row[field]
                existing["merged_notations"] = existing.get("merged_notations", []) + [
                    row["variant"]
                ]
                break
        else:
            kept.append(row)
    return kept


def notation(protein: str | None, cdna: str | None) -> str:
    parts = [p.strip() for p in (protein, cdna) if p and p.strip()]
    return " ".join(dict.fromkeys(parts))


def rows_for_gene(db: Path, pmids: set[str], exclude_layers: set[str]):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    out: dict[str, list[dict]] = defaultdict(list)
    dropped = defaultdict(int)
    q = """
        SELECT vp.pmid              AS pmid,
               vp.variant_id        AS variant_id,
               vp.source_location   AS source_location,
               vp.key_quotes        AS key_quotes,
               vp.source_layer      AS source_layer,
               v.protein_notation   AS protein_notation,
               v.cdna_notation      AS cdna_notation,
               pd.total_carriers_observed AS carriers,
               pd.affected_count    AS affected,
               pd.unaffected_count  AS unaffected
        FROM variant_papers vp
        JOIN variants v ON v.variant_id = vp.variant_id
        LEFT JOIN penetrance_data pd
               ON pd.variant_id = vp.variant_id AND pd.pmid = vp.pmid
    """
    seen: set[tuple[str, str]] = set()
    for r in con.execute(q):
        pmid = str(r["pmid"])
        if pmid not in pmids:
            continue
        layers = {t.strip() for t in (r["source_layer"] or "").split(",") if t.strip()}
        if exclude_layers and layers and layers <= exclude_layers:
            dropped[pmid] += 1
            continue
        var = notation(r["protein_notation"], r["cdna_notation"])
        if not var:
            dropped[pmid] += 1
            continue
        key = (pmid, var)
        if key in seen:  # identical notation string reached via several layers
            continue
        seen.add(key)
        quote = (r["key_quotes"] or "").strip()
        if quote in ("", "[]", "null"):
            quote = f"no quote captured; source_layer={r['source_layer'] or 'unknown'}"
        out[pmid].append(
            {
                "variant": var,
                "carriers": r["carriers"],
                "affected": r["affected"],
                "unaffected": r["unaffected"],
                "evidence": quote[:2000],
                "source_location": (r["source_location"] or "unspecified")[:500],
                "source_layer": r["source_layer"],
            }
        )
    con.close()
    return out, dropped


def trace_telemetry(trace_root: Path, pmids: set[str]) -> dict[str, dict]:
    """Aggregate exact provider telemetry from one finalized gvf-run trace tree.

    A paper that took only a deterministic extraction path legitimately has
    zero calls/tokens; the presence of the finalized trace index distinguishes
    that from missing telemetry.
    """
    if not (trace_root / "trace_index.jsonl").is_file():
        raise SystemExit(f"missing finalized trace index: {trace_root}")

    totals = {
        pmid: {
            "telemetry_available": True,
            "input_tokens": 0,
            "output_tokens": 0,
            "total_tokens": 0,
            "call_count": 0,
            "elapsed_seconds": 0.0,
            "models": set(),
        }
        for pmid in pmids
    }
    for path in trace_root.rglob("*.json"):
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if record.get("record_type") != "llm_call":
            continue
        context = record.get("context") or {}
        pmid = str(context.get("pmid") or "")
        if pmid not in totals:
            continue
        response = record.get("response") or {}
        usage = response.get("usage") or {}
        item = totals[pmid]
        item["input_tokens"] += int(
            usage.get("prompt_tokens") or usage.get("input_tokens") or 0
        )
        item["output_tokens"] += int(
            usage.get("completion_tokens") or usage.get("output_tokens") or 0
        )
        item["total_tokens"] += int(usage.get("total_tokens") or 0)
        item["elapsed_seconds"] += float(response.get("duration_seconds") or 0.0)
        item["call_count"] += 1
        if context.get("model"):
            item["models"].add(str(context["model"]))

    for item in totals.values():
        item["models"] = sorted(item["models"])
        item["elapsed_seconds"] = round(item["elapsed_seconds"], 6)
    return totals


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, default=HERE)
    ap.add_argument("--production-root", type=Path, default=HERE / "production_runs")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--exclude-layers",
        default="",
        help="comma list, e.g. clinvar,pubtator -- drops rows whose layers are all excluded",
    )
    args = ap.parse_args()

    selection = json.loads((args.run_dir / "selection.json").read_text())
    excl = {t.strip() for t in args.exclude_layers.split(",") if t.strip()}

    wanted: dict[str, set[str]] = defaultdict(set)
    for p in selection["papers"]:
        wanted[p["gene"]].add(str(p["pmid"]))

    per_gene = {}
    per_gene_telemetry = {}
    dbs = {}
    merged_away = 0
    for gene, pmids in sorted(wanted.items()):
        db = find_db(args.production_root, gene)
        dbs[gene] = str(db.relative_to(REPO)) if db.is_relative_to(REPO) else str(db)
        per_gene_telemetry[gene] = trace_telemetry(db.parent / "llm_traces", pmids)
        raw, _ = rows_for_gene(db, pmids, excl)
        collapsed = {}
        for pmid, rows in raw.items():
            kept = merge_same_variant(rows, gene)
            merged_away += len(rows) - len(kept)
            collapsed[pmid] = kept
        per_gene[gene] = collapsed

    papers = []
    for p in selection["papers"]:
        gene, pmid = p["gene"], str(p["pmid"])
        variants = per_gene[gene].get(pmid, [])
        telemetry = per_gene_telemetry[gene][pmid]
        papers.append(
            {
                "gene": gene,
                "pmid": pmid,
                "tool": "text",
                "tool_rationale": TOOL_RATIONALE,
                "source_completeness": "corpus_as_locked",
                "elapsed_seconds": telemetry["elapsed_seconds"],
                "token_usage": {
                    key: value
                    for key, value in telemetry.items()
                    if key != "elapsed_seconds"
                },
                "variants": variants,
            }
        )

    aggregate_usage = {
        "telemetry_available": True,
        "input_tokens": sum(p["token_usage"]["input_tokens"] for p in papers),
        "output_tokens": sum(p["token_usage"]["output_tokens"] for p in papers),
        "total_tokens": sum(p["token_usage"]["total_tokens"] for p in papers),
        "call_count": sum(p["token_usage"]["call_count"] for p in papers),
        "models": sorted(
            {model for p in papers for model in p["token_usage"].get("models", [])}
        ),
        "note": "Summed from finalized gvf-run LLM call traces.",
    }

    predictions = {
        "schema_version": 1,
        "run_id": selection["run_id"],
        "strategy": "production_gvf_run",
        "excluded_source_layers": sorted(excl),
        "source_databases": dbs,
        "started_at": None,
        "extraction_started_at": None,
        "extraction_elapsed_seconds": None,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "token_usage": aggregate_usage,
        "prelock_gold_usage": {
            "read_only_layer_scoring_possible": True,
            "scores_feed_back_into_predictions": False,
            "gold_pmid_enrichment_enabled": False,
            "note": (
                "gvf-run automatically emits read-only layer scorecards when "
                "registered gold is available; this external projection is not "
                "the native lock-before-any-gold-read harness."
            ),
        },
        "papers": papers,
    }
    args.out.write_text(json.dumps(predictions, indent=2, sort_keys=True) + "\n")
    total = sum(len(p["variants"]) for p in papers)
    empty = sum(1 for p in papers if not p["variants"])
    print(
        f"{args.out}: {len(papers)} papers, {total} predicted variants, {empty} empty, "
        f"{merged_away} duplicate-notation rows merged"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
