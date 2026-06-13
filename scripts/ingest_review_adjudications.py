#!/usr/bin/env python3
"""Ingest Variant_Browser collaborator adjudications into GVF's gold standard.

This is the round-trip half of the Variant_Browser integration. The publish
half is ``gvf-run --publish-review`` (which shells out to
``Variant_Browser/scripts/gvf_publish.sh``). This script reads the CSV produced
by the browser's export command:

    cd ~/GitRepos/Variant_Browser && set -a && source .env && set +a
    python manage.py export_adjudications [--gene GENE] --out adjudications.csv

and folds each reviewer verdict into a durable CORRECTION OVERLAY under
``gene_variant_fetcher_gold_standard/adjudications/``. The overlay keeps BOTH the
pipeline's extracted carrier counts AND the adjudicated corrections, so recall /
precision / MAE can be recomputed without ever mutating the raw extracted rows.

Matching — the round-trip key is ``(gene, source_notation, pmid)`` where
``source_notation`` is the GVF ``variants.protein_notation``. Each export row is
matched against a run DB's extracted (variant, paper) rows using the *same*
canonical notation bridging the recall scorer uses (``compare_variants``), so a
match here means a match there. Verdicts map to overlay actions:

  confirm        -> gold_confirmed     (extracted record is correct)
  correct_counts -> count_override     (apply corrected affected/unaffected/total)
  wrong_variant  -> false_positive     (extracted variant identity is wrong)
  wrong_paper    -> excluded           (variant/paper association is wrong)
  missing        -> followup_missing   (GVF missed it; queue for a human)
  other          -> followup_other     (free-text comment; queue for a human)

Outputs (idempotent — re-running replaces every gene present in the export):
  <out-dir>/<GENE>_review_adjudications.csv    per-gene correction overlay
  <out-dir>/review_followup_queue.csv          missing/other + unresolved matches
  <out-dir>/review_adjudications_summary.json  per-gene counts + net count deltas

Usage:
    python scripts/ingest_review_adjudications.py --export-csv adjudications.csv
    # attach extracted values from explicit DBs (otherwise auto-found under results/):
    python scripts/ingest_review_adjudications.py --export-csv adjudications.csv \\
        --db KCNH2=results/KCNH2/<ts>/KCNH2.db --db SCN5A=results/SCN5A/<ts>/SCN5A.db
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from cli.compare_variants import (  # noqa: E402
    aggregate_sqlite_data,
    extract_sqlite_data,
    introspect_sqlite,
    normalize_pmid,
    normalize_variant,
    to_canonical_form,
)

# Verdicts the browser emits, mapped to the action this overlay records.
VERDICT_TO_ACTION = {
    "confirm": "gold_confirmed",
    "correct_counts": "count_override",
    "wrong_variant": "false_positive",
    "wrong_paper": "excluded",
    "missing": "followup_missing",
    "other": "followup_other",
}

# Overlay schema. Extracted_* are the pipeline's values (kept), corrected_* are
# the adjudicated values (the overlay). Both survive so metrics can be recomputed.
OVERLAY_COLUMNS = [
    "gene",
    "pmid",
    "source_notation",  # round-trip key == GVF protein_notation
    "variant_label",
    "verdict",
    "action",
    "match_status",  # matched | unmatched | no_db
    "matched_protein_notation",
    "extracted_carriers",
    "extracted_affected",
    "extracted_unaffected",
    "corrected_affected",
    "corrected_unaffected",
    "corrected_total",
    "corrected_classification",
    "comment",
    "adjudicator",
    "updated_at",
]

FOLLOWUP_COLUMNS = [
    "gene",
    "pmid",
    "source_notation",
    "variant_label",
    "verdict",
    "action",
    "match_status",
    "reason",
    "comment",
    "adjudicator",
    "updated_at",
]


def _opt_int(value: Any) -> Optional[int]:
    """Parse an optional integer; blank/garbage -> None."""
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return None


def _variant_key(notation: str) -> str:
    """Canonicalize a protein notation the same way the scorer aggregates the DB.

    Mirrors ``aggregate_sqlite_data`` (canonical form, falling back to the
    normalized form) so an export row keys to the same bucket as the extracted
    (variant, paper) row would.
    """
    canon = to_canonical_form(notation)
    return canon if canon else normalize_variant(notation)


def parse_db_overrides(items: Optional[list[str]]) -> dict[str, Path]:
    """Parse repeatable ``--db GENE=path`` overrides into {GENE: Path}."""
    overrides: dict[str, Path] = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"--db expects GENE=path, got: {item!r}")
        gene, _, path = item.partition("=")
        overrides[gene.strip().upper()] = Path(path.strip()).expanduser()
    return overrides


def find_latest_db(gene: str, results_dir: Path) -> Optional[Path]:
    """Newest ``<GENE>.db`` (or any ``*.db``) under results/<GENE>/."""
    gene_root = results_dir / gene
    if not gene_root.is_dir():
        return None
    named = sorted(gene_root.glob(f"**/{gene}.db"), key=lambda p: p.stat().st_mtime)
    if named:
        return named[-1]
    anydb = sorted(gene_root.glob("**/*.db"), key=lambda p: p.stat().st_mtime)
    return anydb[-1] if anydb else None


def load_db_aggregate(db_path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    """Load a run DB and aggregate its extracted (pmid, variant) carrier rows.

    Reuses the recall scorer's own extraction + aggregation so matching is
    identical to scoring. Returns {} on any read failure (best-effort).
    """
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error:
        return {}
    try:
        table_info = introspect_sqlite(conn)
        df = extract_sqlite_data(conn, table_info)
    except Exception:  # noqa: BLE001 - a malformed DB shouldn't abort ingest
        return {}
    finally:
        conn.close()
    return aggregate_sqlite_data(df)


def _read_export(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def build_overlay_rows(
    export_rows: list[dict[str, str]],
    db_aggregates: dict[str, dict[tuple[str, str], dict[str, Any]]],
    have_db: set[str],
) -> list[dict[str, Any]]:
    """Match each export row to extracted values and classify its action."""
    out: list[dict[str, Any]] = []
    for raw in export_rows:
        gene = (raw.get("gene") or "").strip().upper()
        pmid = normalize_pmid(raw.get("pmid") or "")
        source_notation = (raw.get("source_notation") or "").strip()
        verdict = (raw.get("verdict") or "").strip().lower()
        action = VERDICT_TO_ACTION.get(verdict, "followup_other")

        agg = db_aggregates.get(gene)
        entry: Optional[dict[str, Any]] = None
        if agg is not None and source_notation:
            entry = agg.get((pmid, _variant_key(source_notation)))

        if gene not in have_db:
            match_status = "no_db"
        elif entry is not None:
            match_status = "matched"
        else:
            match_status = "unmatched"

        out.append(
            {
                "gene": gene,
                "pmid": pmid,
                "source_notation": source_notation,
                "variant_label": (raw.get("variant_label") or "").strip(),
                "verdict": verdict,
                "action": action,
                "match_status": match_status,
                "matched_protein_notation": (
                    entry.get("protein_notation") if entry else ""
                )
                or "",
                "extracted_carriers": entry.get("carriers_total") if entry else "",
                "extracted_affected": entry.get("affected_count") if entry else "",
                "extracted_unaffected": entry.get("unaffected_count") if entry else "",
                "corrected_affected": _blank(raw.get("corrected_affected")),
                "corrected_unaffected": _blank(raw.get("corrected_unaffected")),
                "corrected_total": _blank(raw.get("corrected_total")),
                "corrected_classification": _blank(raw.get("corrected_classification")),
                "comment": (raw.get("comment") or "").strip(),
                "adjudicator": (raw.get("adjudicator") or "").strip(),
                "updated_at": (raw.get("updated_at") or "").strip(),
            }
        )
    return out


def _blank(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _needs_followup(row: dict[str, Any]) -> Optional[str]:
    """Return a follow-up reason, or None if the row is cleanly resolved.

    ``missing``/``other`` always queue. An actionable verdict (confirm /
    correct_counts / wrong_variant / wrong_paper) that did not resolve to an
    extracted row also queues — the round-trip key didn't land on anything.
    """
    action = row["action"]
    if action in ("followup_missing", "followup_other"):
        return row["verdict"]
    if row["match_status"] == "no_db":
        return "no run DB available to verify the extracted record"
    if row["match_status"] == "unmatched":
        return "no extracted (variant, paper) row matched the round-trip key"
    return None


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Per-gene counts + net adjudicated count deltas (for metric recompute)."""
    per_gene: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "total": 0,
            "actions": defaultdict(int),
            "match_status": defaultdict(int),
            "count_override_matched": 0,
            "net_affected_delta": 0,
            "net_unaffected_delta": 0,
        }
    )
    for row in rows:
        g = per_gene[row["gene"]]
        g["total"] += 1
        g["actions"][row["action"]] += 1
        g["match_status"][row["match_status"]] += 1
        if row["action"] == "count_override" and row["match_status"] == "matched":
            ext_aff = _opt_int(row["extracted_affected"])
            ext_unaff = _opt_int(row["extracted_unaffected"])
            cor_aff = _opt_int(row["corrected_affected"])
            cor_unaff = _opt_int(row["corrected_unaffected"])
            if ext_aff is not None and cor_aff is not None:
                g["net_affected_delta"] += cor_aff - ext_aff
            if ext_unaff is not None and cor_unaff is not None:
                g["net_unaffected_delta"] += cor_unaff - ext_unaff
            g["count_override_matched"] += 1
    # Convert defaultdicts to plain dicts for JSON.
    return {
        gene: {
            "total": data["total"],
            "actions": dict(data["actions"]),
            "match_status": dict(data["match_status"]),
            "count_override_matched": data["count_override_matched"],
            "net_affected_delta": data["net_affected_delta"],
            "net_unaffected_delta": data["net_unaffected_delta"],
        }
        for gene, data in sorted(per_gene.items())
    }


def write_overlays(rows: list[dict[str, Any]], out_dir: Path) -> dict[str, Path]:
    """Write one correction-overlay CSV per gene (full replace per gene)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_gene[row["gene"]].append(row)
    written: dict[str, Path] = {}
    for gene, gene_rows in sorted(by_gene.items()):
        gene_rows.sort(key=lambda r: (r["pmid"], r["source_notation"]))
        path = out_dir / f"{gene}_review_adjudications.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OVERLAY_COLUMNS, extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(gene_rows)
        written[gene] = path
    return written


def write_followup_queue(rows: list[dict[str, Any]], out_dir: Path) -> Path:
    """Write the combined human follow-up queue."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "review_followup_queue.csv"
    queued = []
    for row in rows:
        reason = _needs_followup(row)
        if reason is None:
            continue
        queued.append(
            {
                "gene": row["gene"],
                "pmid": row["pmid"],
                "source_notation": row["source_notation"],
                "variant_label": row["variant_label"],
                "verdict": row["verdict"],
                "action": row["action"],
                "match_status": row["match_status"],
                "reason": reason,
                "comment": row["comment"],
                "adjudicator": row["adjudicator"],
                "updated_at": row["updated_at"],
            }
        )
    queued.sort(key=lambda r: (r["gene"], r["pmid"], r["source_notation"]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FOLLOWUP_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(queued)
    return path


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--export-csv",
        required=True,
        type=Path,
        help="CSV from Variant_Browser `manage.py export_adjudications`.",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO / "gene_variant_fetcher_gold_standard" / "adjudications",
        help="Where to write the overlay/queue/summary (default: gold adjudications/).",
    )
    ap.add_argument(
        "--db",
        action="append",
        metavar="GENE=PATH",
        help="Explicit run DB for a gene (repeatable). Overrides auto-discovery.",
    )
    ap.add_argument(
        "--results-dir",
        type=Path,
        default=REPO / "results",
        help="Root to auto-discover <GENE>.db when --db is not given for a gene.",
    )
    ap.add_argument(
        "--no-db",
        action="store_true",
        help="Skip DB matching entirely (overlay-only; all rows match_status=no_db).",
    )
    args = ap.parse_args()

    export_path = args.export_csv.expanduser()
    if not export_path.exists():
        raise SystemExit(f"export CSV not found: {export_path}")
    export_rows = _read_export(export_path)
    if not export_rows:
        print(f"No adjudication rows in {export_path}; nothing to ingest.")
        return 0

    genes = sorted({(r.get("gene") or "").strip().upper() for r in export_rows} - {""})
    overrides = parse_db_overrides(args.db)

    db_aggregates: dict[str, dict[tuple[str, str], dict[str, Any]]] = {}
    have_db: set[str] = set()
    if not args.no_db:
        for gene in genes:
            db_path = overrides.get(gene) or find_latest_db(gene, args.results_dir)
            if db_path and db_path.exists():
                agg = load_db_aggregate(db_path)
                db_aggregates[gene] = agg
                have_db.add(gene)
                print(
                    f"  {gene}: matched against {db_path} ({len(agg)} extracted pairs)"
                )
            else:
                print(
                    f"  {gene}: no run DB found — overlay rows will be match_status=no_db"
                )

    rows = build_overlay_rows(export_rows, db_aggregates, have_db)
    written = write_overlays(rows, args.out_dir)
    followup_path = write_followup_queue(rows, args.out_dir)
    summary = summarize(rows)
    summary_path = args.out_dir / "review_adjudications_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    n_followup = sum(1 for r in rows if _needs_followup(r) is not None)
    print(f"\nIngested {len(rows)} adjudications across {len(written)} gene(s):")
    for gene, path in written.items():
        s = summary[gene]
        actions = ", ".join(f"{k}={v}" for k, v in sorted(s["actions"].items()))
        print(f"  {gene}: {s['total']} rows ({actions}) -> {path}")
    print(f"Follow-up queue: {n_followup} rows -> {followup_path}")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
