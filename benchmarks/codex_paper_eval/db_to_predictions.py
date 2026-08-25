#!/usr/bin/env python3
"""Convert production gvf-run SQLite output into codex_paper_eval predictions.json.

The eval harness scores `papers[].variants[]` with {variant, carriers, affected,
unaffected}. Production stores the same facts across three tables, so this is a
straight projection -- no re-interpretation:

  variants protein/cDNA notation, or a structural-only description
                                              -> variants[].variant
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
import hashlib
import re
import json
import sqlite3
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:  # run_eval + utils.variant_normalizer live under the repo
    sys.path.insert(0, str(REPO))

from benchmarks.codex_paper_eval.production_run import (  # noqa: E402
    ProductionRunError,
    resolve_active_gene_run,
)
from utils.variant_normalizer import structural_variant_identity  # noqa: E402

TOOL_RATIONALE = (
    "Production gvf-run calibrated-manifest strategy using the current "
    "deterministic and model-backed extraction, verification, trust, and "
    "recovery layers. Exact code, prompt, model, source, and call-trace "
    "provenance is bound separately by the run and evaluation manifests. "
    "Mapped to the harness 'text' route because no single harness route "
    "describes the multi-stage production pipeline."
)


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def empty_usage() -> dict:
    return {
        "input_tokens": 0,
        "output_tokens": 0,
        "total_tokens": 0,
        "provider_seconds": 0.0,
        "llm_calls": 0,
        "successful_calls": 0,
        "models": {},
    }


def add_usage(target: dict, record: dict) -> None:
    response = record.get("response") or {}
    usage = response.get("usage") or {}
    input_tokens = int(usage.get("prompt_tokens") or usage.get("input_tokens") or 0)
    total_tokens = int(usage.get("total_tokens") or 0)
    output_tokens = max(total_tokens - input_tokens, 0)
    model = str((record.get("context") or {}).get("model") or "unknown")
    target["input_tokens"] += input_tokens
    target["output_tokens"] += output_tokens
    target["total_tokens"] += total_tokens
    target["provider_seconds"] += float(response.get("duration_seconds") or 0)
    target["llm_calls"] += 1
    target["successful_calls"] += int(bool(response.get("success")))
    model_usage = target["models"].setdefault(model, empty_usage())
    # Model buckets do not need a nested copy of themselves.
    model_usage.pop("models", None)
    model_usage["input_tokens"] += input_tokens
    model_usage["output_tokens"] += output_tokens
    model_usage["total_tokens"] += total_tokens
    model_usage["provider_seconds"] += float(response.get("duration_seconds") or 0)
    model_usage["llm_calls"] += 1
    model_usage["successful_calls"] += int(bool(response.get("success")))


def trace_usage(trace_root: Path) -> tuple[dict, dict[str, dict]]:
    total = empty_usage()
    by_pmid: dict[str, dict] = {}
    for path in sorted(trace_root.rglob("*.json")):
        if path.name == "trace_manifest.json":
            continue
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if record.get("record_type") != "llm_call":
            continue
        add_usage(total, record)
        pmid = str((record.get("context") or {}).get("pmid") or "")
        if pmid:
            add_usage(by_pmid.setdefault(pmid, empty_usage()), record)
    return total, by_pmid


def find_db(root: Path, gene: str) -> Path:
    """Resolve the one completed run's declared active DB, never an mtime guess."""
    try:
        _run_dir, db, _status = resolve_active_gene_run(root, gene)
    except ProductionRunError as exc:
        raise SystemExit(str(exc)) from exc
    return db


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
    from benchmarks.codex_paper_eval.run_eval import twin_identical

    def explicit_cdna_tokens(value: str) -> set[str]:
        return {
            match.group(0).casefold()
            for match in re.finditer(
                r"\bc\.\d+(?:[+-]\d+)?(?:_(?:c\.)?\d+(?:[+-]\d+)?)?"
                r"(?:del|dup|ins)[A-Za-z0-9]*",
                value,
                re.IGNORECASE,
            )
        }

    def broad_structural_alias(value: str) -> bool:
        return bool(
            re.search(
                r"\b(?:deletion|duplication|del|dup)\s+of\s+exons?\s*\d+"
                r"|\bexons?\s*\d+\s*(?:deletion|duplication|del|dup)\b",
                value,
                re.IGNORECASE,
            )
        )

    kept: list[dict] = []
    for row in sorted(
        rows, key=lambda r: (layer_rank(r["source_layer"]), r["variant"])
    ):
        for existing in kept:
            row_cdna = explicit_cdna_tokens(row["variant"])
            existing_cdna = explicit_cdna_tokens(existing["variant"])
            if row_cdna and existing_cdna and row_cdna.isdisjoint(existing_cdna):
                # Two exact breakpoints remain distinct even when they share a
                # broad structural alias such as "exon 3 deletion".
                continue
            if (
                bool(row_cdna) != bool(existing_cdna)
                and broad_structural_alias(row["variant"])
                and broad_structural_alias(existing["variant"])
            ):
                # A generic exon event cannot stand in for one exact
                # breakpoint. Preserve both until source adjudication proves
                # that they are the same allele.
                continue
            if twin_identical(row["variant"], existing["variant"], gene):
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


def notation(
    protein: str | None,
    cdna: str | None,
    structural_description: str | None = None,
    variant_class: str | None = None,
) -> str:
    parts = [p.strip() for p in (protein, cdna) if p and p.strip()]
    structural = (structural_description or "").strip()
    if structural and (variant_class or "").strip().lower() in (
        STRUCTURAL_ONLY_VARIANT_CLASSES
    ):
        parts.append(structural)
    return " ".join(dict.fromkeys(parts))


LINKAGE_LAYERS = {"clinvar", "pubtator"}
TRUSTED_UNMATCHED_VF_CLASSES = frozenset(
    {"novel_in_range", "cdna_only_unmatched", "known_isoform_offset"}
)
STRUCTURAL_ONLY_VARIANT_CLASSES = frozenset(
    {
        "large_deletion",
        "large_duplication",
        "cnv",
        "exon_deletion",
        "exon_duplication",
    }
)
TRUSTED_STRUCTURAL_VF_CLASSES = frozenset({"no_notation_suspect"})


def canonical_structural_identity(description: str | None) -> str:
    """Return only structural identities with a canonical exon/gene key."""

    identity = structural_variant_identity(description)
    if re.fullmatch(r"(?:del|dup):(?:exon\d+(?:-\d+)?|wholegene)", identity):
        return identity
    return ""


def origin_layer(source_layer: str | None) -> str:
    raw = (source_layer or "").strip().lower()
    return raw.split(",")[0].strip() if raw else ""


def variant_codons(variant: str | None, gene: str) -> set[int]:
    from benchmarks.codex_paper_eval.run_eval import variant_candidates

    positions: set[int] = set()
    for candidate in variant_candidates(variant or "", gene):
        for match in re.finditer(r"[A-Z](\d{2,5})", candidate.upper()):
            positions.add(int(match.group(1)))
    return positions


def grounded_in_source(variant: str | None, gene: str, text_upper: str) -> bool:
    from benchmarks.codex_paper_eval.run_eval import variant_candidates

    for candidate in variant_candidates(variant or "", gene):
        token = candidate.strip()
        if len(token) < 4 or token.isdigit():
            continue
        if token.upper() in text_upper:
            return True
    return False


def drop_linkage_shadows(
    rows: list[dict], gene: str, source_path: str | None
) -> tuple[list[dict], int]:
    """Exclude DB-linkage rows that are BOTH ungrounded in the paper's locked
    source text AND sit at the exact codon of an independently-extracted row.

    ClinVar submitters cite one landmark PMID for many nearby variants, so the
    citation index dumps same-codon neighbors onto a paper. When the paper's
    own text never mentions the variant and a non-linkage row already occupies
    that codon, the linkage row is an enumeration artifact, not an observation
    from this paper. Measured on the locked gold120 verticalfix predictions:
    55 false positives removed with ZERO true positives lost (the +/-3-residue
    variant of this gate costs a real gold TP and is deliberately not used).

    Rows stay untouched in the production DB — this shapes only the scored
    projection. Papers with no readable source text keep every linkage row:
    on degraded papers the citation index is the only signal, and dropping it
    measurably costs recall (gold120: -53 TP for the grounding-only gate).
    """
    if not rows or not source_path:
        return rows, 0
    path = Path(source_path)
    if not path.is_file():
        return rows, 0
    text_upper = re.sub(r"\s+", " ", path.read_text(errors="replace").upper())
    if not text_upper.strip():
        return rows, 0
    anchor_codons: set[int] = set()
    for row in rows:
        if origin_layer(row.get("source_layer")) not in LINKAGE_LAYERS:
            anchor_codons |= variant_codons(row.get("variant"), gene)
    if not anchor_codons:
        return rows, 0
    kept: list[dict] = []
    dropped = 0
    for row in rows:
        if origin_layer(row.get("source_layer")) in LINKAGE_LAYERS:
            codons = variant_codons(row.get("variant"), gene)
            if (
                codons
                and codons & anchor_codons
                and not grounded_in_source(row.get("variant"), gene, text_upper)
            ):
                dropped += 1
                continue
        kept.append(row)
    return kept, dropped


def trusted_count(
    value: int | None,
    *,
    field: str,
    trust_tier: str | None,
    field_trust: str | None,
) -> int | None:
    """Project one raw count through the persisted field-level trust mask."""
    if field_trust:
        try:
            field_state = json.loads(field_trust)
        except (TypeError, json.JSONDecodeError):
            field_state = None
        if isinstance(field_state, dict):
            return None if field_state.get(field) == "quarantine" else value
    return None if trust_tier == "quarantine" else value


def rows_for_gene(
    db: Path,
    pmids: set[str],
    exclude_layers: set[str],
    *,
    trust_mode: str = "all",
    identity_mode: str = "all",
):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    if trust_mode not in {"all", "trusted"}:
        raise ValueError("trust_mode must be 'all' or 'trusted'")
    if identity_mode not in {"all", "trusted"}:
        raise ValueError("identity_mode must be 'all' or 'trusted'")
    out: dict[str, list[dict]] = defaultdict(list)
    dropped = defaultdict(int)
    penetrance_columns = {
        str(row[1]) for row in con.execute("PRAGMA table_info(penetrance_data)")
    }
    variant_columns = {
        str(row[1]) for row in con.execute("PRAGMA table_info(variants)")
    }
    trust_tier_expr = "pd.trust_tier" if "trust_tier" in penetrance_columns else "NULL"
    field_trust_expr = (
        "pd.field_trust" if "field_trust" in penetrance_columns else "NULL"
    )
    if trust_mode == "trusted" and "trust_tier" not in penetrance_columns:
        con.close()
        raise ValueError(f"{db}: trusted count projection requires trust_tier")
    tables = {
        str(row[0])
        for row in con.execute("SELECT name FROM sqlite_master WHERE type='table'")
    }
    has_vf = "vf_enrichment" in tables
    if identity_mode == "trusted" and not has_vf:
        con.close()
        raise ValueError(f"{db}: trusted identity projection requires vf_enrichment")
    vf_select = (
        "vfe.variant_id AS vf_row, vfe.matched AS vf_matched, vfe.fp_class AS vf_class"
        if has_vf
        else "NULL AS vf_row, NULL AS vf_matched, NULL AS vf_class"
    )
    vf_join = (
        "LEFT JOIN vf_enrichment vfe ON vfe.variant_id = v.variant_id" if has_vf else ""
    )
    variant_class_expr = (
        "v.variant_class" if "variant_class" in variant_columns else "NULL"
    )
    structural_expr = (
        "v.structural_description"
        if "structural_description" in variant_columns
        else "NULL"
    )
    q = f"""
        SELECT vp.pmid              AS pmid,
               vp.variant_id        AS variant_id,
               vp.source_location   AS source_location,
               vp.key_quotes        AS key_quotes,
               vp.source_layer      AS source_layer,
               v.protein_notation   AS protein_notation,
               v.cdna_notation      AS cdna_notation,
               {variant_class_expr} AS variant_class,
               {structural_expr}    AS structural_description,
               pd.total_carriers_observed AS carriers,
               pd.affected_count    AS affected,
               pd.unaffected_count  AS unaffected,
               {trust_tier_expr}    AS trust_tier,
               {field_trust_expr}   AS field_trust,
               {vf_select}
        FROM variant_papers vp
        JOIN variants v ON v.variant_id = vp.variant_id
        {vf_join}
        LEFT JOIN penetrance_data pd
               ON pd.variant_id = vp.variant_id AND pd.pmid = vp.pmid
    """
    seen: set[tuple[str, str]] = set()
    for r in con.execute(q):
        pmid = str(r["pmid"])
        if pmid not in pmids:
            continue
        if identity_mode == "trusted":
            if r["vf_row"] is None:
                con.close()
                raise ValueError(
                    f"{db}: variant {r['variant_id']} lacks VariantFeatures clearance"
                )
            structural_is_trusted = bool(
                str(r["variant_class"] or "").strip().lower()
                in STRUCTURAL_ONLY_VARIANT_CLASSES
                and str(r["vf_class"] or "") in TRUSTED_STRUCTURAL_VF_CLASSES
                and canonical_structural_identity(r["structural_description"])
            )
            if (
                not bool(r["vf_matched"])
                and str(r["vf_class"] or "") not in TRUSTED_UNMATCHED_VF_CLASSES
                and not structural_is_trusted
            ):
                dropped[pmid] += 1
                continue
        layers = {t.strip() for t in (r["source_layer"] or "").split(",") if t.strip()}
        if exclude_layers and layers and layers <= exclude_layers:
            dropped[pmid] += 1
            continue
        var = notation(
            r["protein_notation"],
            r["cdna_notation"],
            r["structural_description"],
            r["variant_class"],
        )
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
        counts = {
            "carriers": r["carriers"],
            "affected": r["affected"],
            "unaffected": r["unaffected"],
        }
        if trust_mode == "trusted":
            counts = {
                field: trusted_count(
                    value,
                    field={
                        "carriers": "total_carriers",
                        "affected": "affected",
                        "unaffected": "unaffected",
                    }[field],
                    trust_tier=r["trust_tier"],
                    field_trust=r["field_trust"],
                )
                for field, value in counts.items()
            }
        out[pmid].append(
            {
                "variant": var,
                **counts,
                "evidence": quote[:2000],
                "source_location": (r["source_location"] or "unspecified")[:500],
                "source_layer": r["source_layer"],
            }
        )
    con.close()
    return out, dropped


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--production-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--exclude-layers",
        default="",
        help="comma list, e.g. clinvar,pubtator -- drops rows whose layers are all excluded",
    )
    ap.add_argument(
        "--trust-mode",
        choices=("all", "trusted"),
        default="all",
        help="all keeps raw counts; trusted masks quarantined count fields",
    )
    ap.add_argument(
        "--identity-mode",
        choices=("all", "trusted"),
        default="trusted",
        help=(
            "trusted holds ambiguous VariantFeatures identity classes out of the "
            "projection while retaining matched, novel-in-range, and cDNA-only variants"
        ),
    )
    ap.add_argument(
        "--keep-linkage-shadows",
        action="store_true",
        help=(
            "keep ungrounded clinvar/pubtator rows that sit at the exact codon "
            "of an independently-extracted row (default: excluded from the "
            "scored projection; DB rows are never touched)"
        ),
    )
    args = ap.parse_args()

    selection = json.loads((args.run_dir / "selection.json").read_text())
    excl = {t.strip() for t in args.exclude_layers.split(",") if t.strip()}

    wanted: dict[str, set[str]] = defaultdict(set)
    source_by_key: dict[tuple[str, str], str | None] = {}
    for p in selection["papers"]:
        wanted[p["gene"]].add(str(p["pmid"]))
        source_by_key[(p["gene"], str(p["pmid"]))] = p.get("source")
    linkage_shadows_excluded = 0

    per_gene = {}
    dbs = {}
    production_trace_manifests = []
    production_run_timing = []
    production_gold_access = []
    run_usage = empty_usage()
    usage_by_gene_pmid: dict[tuple[str, str], dict] = {}
    merged_away = 0
    for gene, pmids in sorted(wanted.items()):
        db = find_db(args.production_root, gene)
        dbs[gene] = str(db.relative_to(REPO)) if db.is_relative_to(REPO) else str(db)
        trace_manifest_path = db.parent / "llm_traces" / "trace_manifest.json"
        if not trace_manifest_path.is_file():
            raise SystemExit(
                f"no finalized LLM trace manifest for {gene}: {trace_manifest_path}"
            )
        trace_manifest = json.loads(trace_manifest_path.read_text())
        run_manifest_path = db.parent / "run_manifest.json"
        run_status_path = db.parent / "RUN_STATUS.json"
        run_manifest = (
            json.loads(run_manifest_path.read_text())
            if run_manifest_path.is_file()
            else {}
        )
        run_status = (
            json.loads(run_status_path.read_text()) if run_status_path.is_file() else {}
        )
        if (
            run_status.get("status") != "completed"
            or run_status.get("exit_code") != 0
            or run_status.get("stage_failures")
        ):
            raise SystemExit(f"{gene}: active production run did not complete cleanly")
        gold_access = run_status.get("gold_access")
        gold_disabled = bool(
            isinstance(gold_access, dict) and gold_access.get("disabled") is True
        )
        production_gold_access.append(
            {
                "gene": gene,
                "disabled": gold_disabled,
                "mode": (
                    gold_access.get("mode")
                    if isinstance(gold_access, dict)
                    else "legacy_status_missing_provenance"
                ),
            }
        )
        active_raw = str(run_status.get("active_db") or "").strip()
        active_db = Path(active_raw).expanduser()
        if not active_db.is_absolute():
            active_db = db.parent / active_db
        if active_db.resolve() != db.resolve():
            raise SystemExit(f"{gene}: projection DB is not RUN_STATUS active_db")
        verification = trace_manifest.get("verification") or {}
        if verification.get("level") != "write_time_verified":
            raise SystemExit(f"{gene}: production trace is not write-time verified")
        if verification.get("errors"):
            raise SystemExit(f"{gene}: production trace integrity errors are present")
        if int(trace_manifest.get("llm_call_count") or 0) <= 0:
            raise SystemExit(f"{gene}: production trace contains no model calls")
        start_timestamp = run_manifest.get("start_timestamp")
        duration_seconds = run_status.get("duration_seconds")
        timing_source = "RUN_STATUS.json"
        if duration_seconds is None:
            # An operator may interrupt only the optional post-run housekeeping
            # after the final DB/source-QC artifacts are complete.  Preserve the
            # useful end-to-end timing without using a later manifest rebuild's
            # mtime, which would inflate wall time.
            completion_candidates = [
                db,
                db.parent / "source_qc" / "source_acquisition_summary.json",
                db.parent / "layers" / "figure_reads" / "summary.json",
            ]
            completion_times = [
                path.stat().st_mtime for path in completion_candidates if path.is_file()
            ]
            if start_timestamp and completion_times:
                # gvf-run currently records a naive local timestamp; timestamp()
                # applies the host timezone, matching the filesystem mtime.
                start = datetime.fromisoformat(start_timestamp)
                duration_seconds = max(completion_times) - start.timestamp()
                timing_source = "final_artifact_mtime_minus_run_start"
        production_run_timing.append(
            {
                "gene": gene,
                "start_timestamp": start_timestamp,
                "duration_seconds": duration_seconds,
                "source": timing_source,
            }
        )
        production_trace_manifests.append(
            {
                "gene": gene,
                "manifest": (
                    str(trace_manifest_path.relative_to(REPO))
                    if trace_manifest_path.is_relative_to(REPO)
                    else str(trace_manifest_path)
                ),
                "sha256": file_digest(trace_manifest_path),
                "run_id": trace_manifest.get("run_id"),
                "llm_call_count": trace_manifest.get("llm_call_count"),
                "decision_event_count": trace_manifest.get("decision_event_count"),
                "integrity_level": (trace_manifest.get("verification") or {}).get(
                    "level"
                ),
            }
        )
        gene_usage, by_pmid = trace_usage(trace_manifest_path.parent)
        for field in (
            "input_tokens",
            "output_tokens",
            "total_tokens",
            "provider_seconds",
            "llm_calls",
            "successful_calls",
        ):
            run_usage[field] += gene_usage[field]
        for model, model_usage in gene_usage["models"].items():
            bucket = run_usage["models"].setdefault(model, empty_usage())
            bucket.pop("models", None)
            for field in (
                "input_tokens",
                "output_tokens",
                "total_tokens",
                "provider_seconds",
                "llm_calls",
                "successful_calls",
            ):
                bucket[field] += model_usage[field]
        usage_by_gene_pmid.update(
            {(gene, pmid): usage for pmid, usage in by_pmid.items()}
        )
        raw, _ = rows_for_gene(
            db,
            pmids,
            excl,
            trust_mode=args.trust_mode,
            identity_mode=args.identity_mode,
        )
        collapsed = {}
        for pmid, rows in raw.items():
            kept = merge_same_variant(rows, gene)
            merged_away += len(rows) - len(kept)
            if not args.keep_linkage_shadows:
                kept, shadows = drop_linkage_shadows(
                    kept, gene, source_by_key.get((gene, pmid))
                )
                linkage_shadows_excluded += shadows
            collapsed[pmid] = kept
        per_gene[gene] = collapsed

    papers = []
    for p in selection["papers"]:
        gene, pmid = p["gene"], str(p["pmid"])
        variants = per_gene[gene].get(pmid, [])
        paper_usage = usage_by_gene_pmid.get((gene, pmid), empty_usage())
        papers.append(
            {
                "gene": gene,
                "pmid": pmid,
                "tool": "text",
                "tool_rationale": TOOL_RATIONALE,
                "source_completeness": "corpus_as_locked",
                "elapsed_seconds": paper_usage["provider_seconds"],
                "token_usage": {
                    "telemetry_available": True,
                    **paper_usage,
                    "note": (
                        "Trace-derived. output_tokens is total minus input and "
                        "therefore conservatively includes billed reasoning tokens."
                    ),
                },
                "variants": variants,
            }
        )

    predictions = {
        "schema_version": 1,
        "run_id": selection["run_id"],
        "strategy": "production_gvf_run",
        "count_projection": args.trust_mode,
        "identity_projection": args.identity_mode,
        "excluded_source_layers": sorted(excl),
        "linkage_shadows_excluded": linkage_shadows_excluded,
        "source_databases": dbs,
        "production_trace_manifests": production_trace_manifests,
        "production_run_timing": production_run_timing,
        "production_gold_access": production_gold_access,
        "prelock_gold_usage": {
            "read_only_layer_scoring_possible": not all(
                row["disabled"] for row in production_gold_access
            ),
            "scores_feed_back_into_predictions": False,
            "provenance": (
                "all_production_run_statuses_record_gold_access_disabled"
                if all(row["disabled"] for row in production_gold_access)
                else "one_or_more_run_statuses_allow_or_do_not_prove_gold_disabled"
            ),
        },
        "started_at": None,
        "extraction_started_at": min(
            (
                row["start_timestamp"]
                for row in production_run_timing
                if row["start_timestamp"]
            ),
            default=None,
        ),
        "extraction_elapsed_seconds": sum(
            float(row["duration_seconds"] or 0) for row in production_run_timing
        ),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "token_usage": {
            "telemetry_available": True,
            **run_usage,
            "note": (
                "Trace-derived. output_tokens is total minus input and therefore "
                "conservatively includes billed reasoning tokens; provider_seconds "
                "is summed call latency, not end-to-end wall clock."
            ),
        },
        "papers": papers,
    }
    args.out.write_text(json.dumps(predictions, indent=2, sort_keys=True) + "\n")
    total = sum(len(p["variants"]) for p in papers)
    empty = sum(1 for p in papers if not p["variants"])
    print(
        f"{args.out}: {len(papers)} papers, {total} predicted variants, {empty} empty, "
        f"{merged_away} duplicate-notation rows merged, "
        f"{linkage_shadows_excluded} linkage codon-shadows excluded"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
