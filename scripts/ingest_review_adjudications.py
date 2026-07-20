#!/usr/bin/env python3
"""Ingest Variant_Browser collaborator adjudications into GVF's gold standard.

This is the round-trip half of the Variant_Browser integration. The normal path
pulls a bearer-authenticated JSON snapshot directly from the live Azure-backed
review database and atomically replaces GVF's local ``review_gold.sqlite3``
cache. A lead-approved CSV remains supported only as an offline compatibility
input. The cache keeps BOTH the pipeline's extracted carrier counts AND the
adjudicated corrections, so recall / precision / MAE / RMSE can be recomputed
without mutating raw extracted rows.

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

Outputs (idempotent current view + immutable source history):
  <out-dir>/review_gold.sqlite3                live gold + metric overlays + sync audit
  <out-dir>/review_followup_queue.csv          missing/other + unresolved matches
  <out-dir>/review_adjudications_summary.json  per-gene counts + net count deltas

Usage:
    GVF_REVIEW_GOLD_TOKEN=... python scripts/ingest_review_adjudications.py \
        --source-url https://variantbrowser.org/review/api/gold-standard/
    # Offline compatibility input:
    python scripts/ingest_review_adjudications.py --export-csv adjudications.csv
    # attach extracted values from explicit DBs (otherwise auto-found under results/):
    python scripts/ingest_review_adjudications.py --export-csv adjudications.csv \\
        --db KCNH2=results/KCNH2/<ts>/KCNH2.db --db SCN5A=results/SCN5A/<ts>/SCN5A.db
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sqlite3
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
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
    "record_key",
    "gene",
    "pmid",
    "source_notation",  # round-trip key == GVF protein_notation
    "variant_label",
    "status",
    "revision",
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
    "source_reviewer_user_id",
    "source_reviewer",
    "decided_by_user_id",
    "decided_by",
    "adjudicator",
    "updated_at",
]

FOLLOWUP_COLUMNS = [
    "record_key",
    "gene",
    "pmid",
    "source_notation",
    "variant_label",
    "status",
    "revision",
    "verdict",
    "action",
    "match_status",
    "reason",
    "comment",
    "source_reviewer_user_id",
    "source_reviewer",
    "decided_by_user_id",
    "decided_by",
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
        # Adjudication ingest must see preserved raw evidence; trust is a later,
        # reversible scoring projection and must not alter source-key matching.
        df = extract_sqlite_data(conn, table_info, trust_mode="all")
    except Exception:  # noqa: BLE001 - a malformed DB shouldn't abort ingest
        return {}
    finally:
        conn.close()
    return aggregate_sqlite_data(df)


GOLD_EXPORT_MARKERS = {
    "record_key",
    "status",
    "revision",
    "source_reviewer_user_id",
    "source_reviewer",
    "decided_by_user_id",
    "decided_by",
}
ACCEPTED_GOLD_STATUSES = {"gold_standard", "adjudicated"}
LIVE_GOLD_SCHEMA_VERSION = 1
DEFAULT_LIVE_GOLD_URL = "https://variantbrowser.org/review/api/gold-standard/"
DEFAULT_CACHE_DB = (
    REPO
    / "gene_variant_fetcher_gold_standard"
    / "adjudications"
    / "review_gold.sqlite3"
)
CARDIAC_GENES = frozenset({"KCNH2", "KCNQ1", "RYR2", "SCN5A"})
BUILTIN_GOLD_TIERS: dict[str, dict[str, Any]] = {
    "all": {
        "description": "All current lead-approved records, including non-cardiac genes.",
        "gene_mode": "all",
        "genes": [],
    },
    "cardiac": {
        "description": "Canonical cardiac channelopathy genes only.",
        "gene_mode": "include",
        "genes": sorted(CARDIAC_GENES),
    },
    "noncardiac": {
        "description": "Current records outside the canonical cardiac gene set.",
        "gene_mode": "exclude",
        "genes": sorted(CARDIAC_GENES),
    },
}


class GoldSyncError(RuntimeError):
    """The live gold contract could not be authenticated or validated."""


def _validate_gold_rows(
    rows: list[dict[str, Any]],
    fields: set[str],
    *,
    source_label: str,
) -> list[dict[str, Any]]:
    missing = sorted(GOLD_EXPORT_MARKERS - fields)
    if missing:
        raise GoldSyncError(
            f"Expected Variant_Browser lead-approved gold from {source_label}; "
            f"missing approval columns: {', '.join(missing)}. Do not ingest the "
            "multi-reviewer export_adjudications audit file."
        )
    unsafe = sorted(
        {
            status
            for row in rows
            if (status := str(row.get("status") or "").strip().lower())
            not in ACCEPTED_GOLD_STATUSES
        }
    )
    if unsafe:
        raise GoldSyncError(
            "Gold source contains non-accepted status(es): "
            f"{', '.join(status or '<blank>' for status in unsafe)}. "
            "Disputed/withheld audit rows must not enter metric overlays."
        )
    record_keys = [str(row.get("record_key") or "").strip() for row in rows]
    if any(not key for key in record_keys):
        raise GoldSyncError("Gold source contains a blank record_key.")
    if len(set(record_keys)) != len(record_keys):
        raise GoldSyncError("Gold source contains duplicate current record_key values.")
    for row in rows:
        if not str(row.get("gene") or "").strip():
            raise GoldSyncError("Gold source contains a blank gene.")
        try:
            revision = int(row.get("revision") or 0)
        except (TypeError, ValueError) as exc:
            raise GoldSyncError("Gold source contains a non-integer revision.") from exc
        if revision < 1:
            raise GoldSyncError("Gold source contains a non-positive revision.")
        if not str(row.get("source_reviewer_user_id") or "").strip():
            raise GoldSyncError(
                "Gold source is missing the source reviewer's stable account ID."
            )
        if not str(row.get("decided_by_user_id") or "").strip():
            raise GoldSyncError(
                "Gold source is missing the approving lead's stable account ID."
            )
    return rows


def _read_export(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        rows = list(reader)
    try:
        return _validate_gold_rows(
            rows,
            fields,
            source_label="export_gold_standard CSV",
        )
    except GoldSyncError as exc:
        raise SystemExit(str(exc)) from exc


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Never forward the bearer token to a redirected host."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ANN001
        raise urllib.error.HTTPError(
            req.full_url,
            code,
            f"Redirect to {newurl} blocked for security",
            headers,
            fp,
        )


def _validated_source_url(source_url: str) -> str:
    parsed = urllib.parse.urlsplit(source_url)
    local_host = parsed.hostname in {"localhost", "127.0.0.1", "::1"}
    if parsed.scheme != "https" and not (parsed.scheme == "http" and local_host):
        raise GoldSyncError("Live gold sync requires HTTPS (except localhost tests).")
    if parsed.username or parsed.password or parsed.fragment:
        raise GoldSyncError("Live gold URL must not contain credentials or a fragment.")
    if not parsed.hostname:
        raise GoldSyncError("Live gold URL has no hostname.")
    return urllib.parse.urlunsplit(parsed)


def fetch_live_gold(
    source_url: str,
    token: str,
    *,
    timeout_s: int = 30,
    opener: Any = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Fetch and independently validate one complete Azure gold snapshot."""
    source_url = _validated_source_url(source_url)
    token = (token or "").strip()
    if len(token) < 32:
        raise GoldSyncError("GVF_REVIEW_GOLD_TOKEN is missing or too short.")
    request = urllib.request.Request(
        source_url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
            "User-Agent": "gvf-gold-sync/1",
        },
        method="GET",
    )
    opener = opener or urllib.request.build_opener(_NoRedirect())
    try:
        with opener.open(request, timeout=timeout_s) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as exc:
        raise GoldSyncError(
            f"Variant Browser gold sync returned HTTP {exc.code}."
        ) from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise GoldSyncError(f"Variant Browser gold sync failed: {exc}") from exc

    if not isinstance(payload, dict) or payload.get("ok") is not True:
        raise GoldSyncError("Variant Browser returned an invalid gold response.")
    if payload.get("schema_version") != LIVE_GOLD_SCHEMA_VERSION:
        raise GoldSyncError(
            f"Unsupported live gold schema version: {payload.get('schema_version')!r}."
        )
    rows = payload.get("records")
    columns = payload.get("columns")
    if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
        raise GoldSyncError("Live gold response records must be a list of objects.")
    if not isinstance(columns, list) or not all(
        isinstance(col, str) for col in columns
    ):
        raise GoldSyncError("Live gold response is missing its columns contract.")
    if payload.get("record_count") != len(rows):
        raise GoldSyncError("Live gold record_count does not match the payload.")
    _validate_gold_rows(rows, set(columns), source_label="live Azure API")

    canonical = json.dumps(
        rows,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    computed_revision = hashlib.sha256(canonical).hexdigest()
    if payload.get("dataset_revision") != computed_revision:
        raise GoldSyncError("Live gold dataset revision checksum does not match.")
    metadata = {
        "source_url": source_url,
        "source": str(payload.get("source") or ""),
        "schema_version": payload["schema_version"],
        "dataset_revision": computed_revision,
        "record_count": len(rows),
        "generated_at": str(payload.get("generated_at") or ""),
    }
    return rows, metadata


def build_overlay_rows(
    export_rows: list[dict[str, Any]],
    db_aggregates: dict[str, dict[tuple[str, str], dict[str, Any]]],
    have_db: set[str],
) -> list[dict[str, Any]]:
    """Match each export row to extracted values and classify its action."""
    out: list[dict[str, Any]] = []
    for raw in export_rows:
        # CSV inputs arrive as strings, while the live JSON contract preserves
        # native scalar types (notably ``revision`` as an integer).  Normalize
        # every text-bound field here so both transports follow the same path.
        gene = _blank(raw.get("gene")).upper()
        pmid = normalize_pmid(raw.get("pmid") or "")
        source_notation = _blank(raw.get("source_notation"))
        verdict = _blank(raw.get("verdict")).lower()
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
                "record_key": _blank(raw.get("record_key")),
                "gene": gene,
                "pmid": pmid,
                "source_notation": source_notation,
                "variant_label": _blank(raw.get("variant_label")),
                "status": _blank(raw.get("status")),
                "revision": _blank(raw.get("revision")),
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
                "comment": _blank(raw.get("comment")),
                "source_reviewer_user_id": _blank(raw.get("source_reviewer_user_id")),
                "source_reviewer": _blank(raw.get("source_reviewer")),
                "decided_by_user_id": _blank(raw.get("decided_by_user_id")),
                "decided_by": _blank(raw.get("decided_by")),
                "adjudicator": _blank(raw.get("adjudicator")),
                "updated_at": _blank(raw.get("updated_at")),
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


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def ensure_live_gold_schema(conn: sqlite3.Connection) -> None:
    """Create/migrate the live cache, immutable history, tiers, and audit ledgers."""
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS review_gold_records (
            record_key TEXT PRIMARY KEY,
            gene TEXT NOT NULL,
            pmid TEXT NOT NULL DEFAULT '',
            source_notation TEXT NOT NULL DEFAULT '',
            status TEXT NOT NULL,
            revision INTEGER NOT NULL,
            source_reviewer_user_id TEXT NOT NULL DEFAULT '',
            decided_by_user_id TEXT NOT NULL DEFAULT '',
            payload_json TEXT NOT NULL,
            source_dataset_revision TEXT NOT NULL,
            synced_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_review_gold_records_gene
            ON review_gold_records(gene, pmid);
        CREATE TABLE IF NOT EXISTS review_gold_overlays (
            record_key TEXT PRIMARY KEY,
            gene TEXT NOT NULL,
            pmid TEXT NOT NULL DEFAULT '',
            source_notation TEXT NOT NULL DEFAULT '',
            action TEXT NOT NULL,
            match_status TEXT NOT NULL,
            payload_json TEXT NOT NULL,
            source_dataset_revision TEXT NOT NULL,
            synced_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_review_gold_overlays_gene
            ON review_gold_overlays(gene, pmid);
        CREATE TABLE IF NOT EXISTS review_gold_sync_state (
            singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
            source_url TEXT NOT NULL,
            source TEXT NOT NULL,
            schema_version INTEGER NOT NULL,
            dataset_revision TEXT NOT NULL,
            record_count INTEGER NOT NULL,
            generated_at TEXT NOT NULL,
            synced_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS review_gold_snapshots (
            dataset_revision TEXT PRIMARY KEY,
            source_url TEXT NOT NULL,
            source TEXT NOT NULL,
            schema_version INTEGER NOT NULL,
            record_count INTEGER NOT NULL,
            cardiac_record_count INTEGER NOT NULL,
            noncardiac_record_count INTEGER NOT NULL,
            generated_at TEXT NOT NULL,
            first_synced_at TEXT NOT NULL,
            last_synced_at TEXT NOT NULL,
            sync_count INTEGER NOT NULL DEFAULT 1
        );
        CREATE TABLE IF NOT EXISTS review_gold_snapshot_records (
            dataset_revision TEXT NOT NULL,
            record_key TEXT NOT NULL,
            gene TEXT NOT NULL,
            pmid TEXT NOT NULL DEFAULT '',
            source_notation TEXT NOT NULL DEFAULT '',
            status TEXT NOT NULL,
            revision INTEGER NOT NULL,
            source_reviewer_user_id TEXT NOT NULL DEFAULT '',
            decided_by_user_id TEXT NOT NULL DEFAULT '',
            payload_json TEXT NOT NULL,
            PRIMARY KEY (dataset_revision, record_key)
        );
        CREATE INDEX IF NOT EXISTS idx_review_gold_snapshot_records_gene
            ON review_gold_snapshot_records(dataset_revision, gene, pmid);
        CREATE TABLE IF NOT EXISTS review_gold_sync_runs (
            sync_id INTEGER PRIMARY KEY AUTOINCREMENT,
            dataset_revision TEXT NOT NULL,
            previous_dataset_revision TEXT NOT NULL DEFAULT '',
            source_url TEXT NOT NULL,
            record_count INTEGER NOT NULL,
            cardiac_record_count INTEGER NOT NULL,
            noncardiac_record_count INTEGER NOT NULL,
            added_count INTEGER NOT NULL,
            updated_count INTEGER NOT NULL,
            removed_count INTEGER NOT NULL,
            generated_at TEXT NOT NULL,
            synced_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_review_gold_sync_runs_revision
            ON review_gold_sync_runs(dataset_revision, synced_at);
        CREATE TABLE IF NOT EXISTS review_gold_record_changes (
            sync_id INTEGER NOT NULL,
            dataset_revision TEXT NOT NULL,
            previous_dataset_revision TEXT NOT NULL DEFAULT '',
            record_key TEXT NOT NULL,
            gene TEXT NOT NULL,
            change_type TEXT NOT NULL CHECK (
                change_type IN ('added', 'updated', 'removed')
            ),
            previous_payload_json TEXT NOT NULL DEFAULT '',
            current_payload_json TEXT NOT NULL DEFAULT '',
            detected_at TEXT NOT NULL,
            PRIMARY KEY (sync_id, record_key)
        );
        CREATE INDEX IF NOT EXISTS idx_review_gold_record_changes_key
            ON review_gold_record_changes(record_key, sync_id);
        CREATE TABLE IF NOT EXISTS review_gold_tiers (
            name TEXT PRIMARY KEY,
            description TEXT NOT NULL DEFAULT '',
            gene_mode TEXT NOT NULL CHECK (gene_mode IN ('all', 'include', 'exclude')),
            genes_json TEXT NOT NULL DEFAULT '[]',
            builtin INTEGER NOT NULL DEFAULT 0,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            updated_by TEXT NOT NULL DEFAULT ''
        );
        CREATE TABLE IF NOT EXISTS review_gold_exclusions (
            record_key TEXT PRIMARY KEY,
            active INTEGER NOT NULL CHECK (active IN (0, 1)),
            actor TEXT NOT NULL,
            reason TEXT NOT NULL,
            updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS review_gold_exclusion_events (
            event_id INTEGER PRIMARY KEY AUTOINCREMENT,
            record_key TEXT NOT NULL,
            action TEXT NOT NULL CHECK (action IN ('exclude', 'restore')),
            actor TEXT NOT NULL,
            reason TEXT NOT NULL,
            created_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_review_gold_exclusion_events_key
            ON review_gold_exclusion_events(record_key, event_id);
        """
    )
    now = datetime.now(timezone.utc).isoformat()
    for name, definition in BUILTIN_GOLD_TIERS.items():
        conn.execute(
            """
            INSERT OR IGNORE INTO review_gold_tiers (
                name, description, gene_mode, genes_json, builtin,
                created_at, updated_at, updated_by
            ) VALUES (?, ?, ?, ?, 1, ?, ?, 'system')
            """,
            (
                name,
                definition["description"],
                definition["gene_mode"],
                _canonical_json(definition["genes"]),
                now,
                now,
            ),
        )


def _tier_definition_conn(conn: sqlite3.Connection, tier: str) -> dict[str, Any]:
    name = str(tier or "").strip().lower()
    row = conn.execute(
        "SELECT name, description, gene_mode, genes_json, builtin "
        "FROM review_gold_tiers WHERE name = ?",
        (name,),
    ).fetchone()
    if row is None:
        raise GoldSyncError(f"Unknown review-gold tier: {tier!r}.")
    genes = json.loads(row[3])
    if not isinstance(genes, list):
        raise GoldSyncError(f"Review-gold tier {name!r} has invalid genes JSON.")
    return {
        "name": row[0],
        "description": row[1],
        "gene_mode": row[2],
        "genes": sorted({str(g).strip().upper() for g in genes if str(g).strip()}),
        "builtin": bool(row[4]),
    }


def _tier_includes_gene(definition: dict[str, Any], gene: str) -> bool:
    normalized = str(gene or "").strip().upper()
    genes = set(definition["genes"])
    if definition["gene_mode"] == "all":
        return True
    if definition["gene_mode"] == "include":
        return normalized in genes
    return normalized not in genes


def _count_tier_records_conn(
    conn: sqlite3.Connection,
    tier: str,
    *,
    include_excluded: bool = False,
) -> int:
    definition = _tier_definition_conn(conn, tier)
    excluded = set()
    if not include_excluded:
        excluded = {
            row[0]
            for row in conn.execute(
                "SELECT record_key FROM review_gold_exclusions WHERE active = 1"
            )
        }
    return sum(
        1
        for record_key, gene in conn.execute(
            "SELECT record_key, gene FROM review_gold_records"
        )
        if record_key not in excluded and _tier_includes_gene(definition, gene)
    )


def _archive_existing_current_snapshot(conn: sqlite3.Connection) -> None:
    """Backfill a pre-v2 cache before its current view is replaced."""
    state = conn.execute(
        "SELECT source_url, source, schema_version, dataset_revision, "
        "record_count, generated_at, synced_at "
        "FROM review_gold_sync_state WHERE singleton = 1"
    ).fetchone()
    if state is None or not state[3]:
        return
    exists = conn.execute(
        "SELECT 1 FROM review_gold_snapshots WHERE dataset_revision = ?",
        (state[3],),
    ).fetchone()
    if exists is not None:
        return
    records = conn.execute(
        """
        SELECT record_key, gene, pmid, source_notation, status, revision,
               source_reviewer_user_id, decided_by_user_id, payload_json
        FROM review_gold_records
        """
    ).fetchall()
    cardiac_count = sum(1 for row in records if row[1] in CARDIAC_GENES)
    conn.execute(
        """
        INSERT INTO review_gold_snapshots (
            dataset_revision, source_url, source, schema_version, record_count,
            cardiac_record_count, noncardiac_record_count, generated_at,
            first_synced_at, last_synced_at, sync_count
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)
        """,
        (
            state[3],
            state[0],
            state[1],
            state[2],
            len(records),
            cardiac_count,
            len(records) - cardiac_count,
            state[5],
            state[6],
            state[6],
        ),
    )
    conn.executemany(
        """
        INSERT INTO review_gold_snapshot_records (
            dataset_revision, record_key, gene, pmid, source_notation, status,
            revision, source_reviewer_user_id, decided_by_user_id, payload_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [(state[3], *row) for row in records],
    )


def write_live_gold_cache(
    path: Path,
    gold_rows: list[dict[str, Any]],
    overlay_rows: list[dict[str, Any]],
    metadata: dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    """Update the current view while preserving every source snapshot and change."""
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    synced_at = datetime.now(timezone.utc).isoformat()
    conn = sqlite3.connect(path)
    try:
        with conn:
            ensure_live_gold_schema(conn)
            _archive_existing_current_snapshot(conn)
            previous_state = conn.execute(
                "SELECT dataset_revision FROM review_gold_sync_state "
                "WHERE singleton = 1"
            ).fetchone()
            previous_revision = previous_state[0] if previous_state else ""
            previous_rows = {
                row[0]: {"gene": row[1], "payload_json": row[2]}
                for row in conn.execute(
                    "SELECT record_key, gene, payload_json FROM review_gold_records"
                )
            }
            prepared_gold = [
                {
                    "record_key": str(row.get("record_key") or "").strip(),
                    "gene": str(row.get("gene") or "").strip().upper(),
                    "pmid": normalize_pmid(row.get("pmid") or ""),
                    "source_notation": str(row.get("source_notation") or "").strip(),
                    "status": str(row.get("status") or "").strip(),
                    "revision": int(row.get("revision") or 0),
                    "source_reviewer_user_id": str(
                        row.get("source_reviewer_user_id") or ""
                    ).strip(),
                    "decided_by_user_id": str(
                        row.get("decided_by_user_id") or ""
                    ).strip(),
                    "payload_json": _canonical_json(row),
                }
                for row in gold_rows
            ]
            current_rows = {
                row["record_key"]: {
                    "gene": row["gene"],
                    "payload_json": row["payload_json"],
                }
                for row in prepared_gold
            }
            added = sorted(set(current_rows) - set(previous_rows))
            removed = sorted(set(previous_rows) - set(current_rows))
            updated = sorted(
                key
                for key in set(previous_rows) & set(current_rows)
                if previous_rows[key]["payload_json"]
                != current_rows[key]["payload_json"]
            )
            cardiac_count = sum(
                1 for row in prepared_gold if row["gene"] in CARDIAC_GENES
            )
            conn.execute(
                """
                INSERT INTO review_gold_snapshots (
                    dataset_revision, source_url, source, schema_version,
                    record_count, cardiac_record_count, noncardiac_record_count,
                    generated_at, first_synced_at, last_synced_at, sync_count
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)
                ON CONFLICT(dataset_revision) DO UPDATE SET
                    last_synced_at=excluded.last_synced_at,
                    sync_count=review_gold_snapshots.sync_count + 1
                """,
                (
                    metadata["dataset_revision"],
                    metadata["source_url"],
                    metadata["source"],
                    metadata["schema_version"],
                    len(prepared_gold),
                    cardiac_count,
                    len(prepared_gold) - cardiac_count,
                    metadata["generated_at"],
                    synced_at,
                    synced_at,
                ),
            )
            conn.executemany(
                """
                INSERT OR IGNORE INTO review_gold_snapshot_records (
                    dataset_revision, record_key, gene, pmid, source_notation,
                    status, revision, source_reviewer_user_id,
                    decided_by_user_id, payload_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        metadata["dataset_revision"],
                        row["record_key"],
                        row["gene"],
                        row["pmid"],
                        row["source_notation"],
                        row["status"],
                        row["revision"],
                        row["source_reviewer_user_id"],
                        row["decided_by_user_id"],
                        row["payload_json"],
                    )
                    for row in prepared_gold
                ],
            )
            cursor = conn.execute(
                """
                INSERT INTO review_gold_sync_runs (
                    dataset_revision, previous_dataset_revision, source_url,
                    record_count, cardiac_record_count, noncardiac_record_count,
                    added_count, updated_count, removed_count, generated_at, synced_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    metadata["dataset_revision"],
                    previous_revision,
                    metadata["source_url"],
                    len(prepared_gold),
                    cardiac_count,
                    len(prepared_gold) - cardiac_count,
                    len(added),
                    len(updated),
                    len(removed),
                    metadata["generated_at"],
                    synced_at,
                ),
            )
            sync_id = int(cursor.lastrowid)
            change_rows = []
            for key in added:
                change_rows.append(
                    (
                        sync_id,
                        metadata["dataset_revision"],
                        previous_revision,
                        key,
                        current_rows[key]["gene"],
                        "added",
                        "",
                        current_rows[key]["payload_json"],
                        synced_at,
                    )
                )
            for key in updated:
                change_rows.append(
                    (
                        sync_id,
                        metadata["dataset_revision"],
                        previous_revision,
                        key,
                        current_rows[key]["gene"],
                        "updated",
                        previous_rows[key]["payload_json"],
                        current_rows[key]["payload_json"],
                        synced_at,
                    )
                )
            for key in removed:
                change_rows.append(
                    (
                        sync_id,
                        metadata["dataset_revision"],
                        previous_revision,
                        key,
                        previous_rows[key]["gene"],
                        "removed",
                        previous_rows[key]["payload_json"],
                        "",
                        synced_at,
                    )
                )
            conn.executemany(
                """
                INSERT INTO review_gold_record_changes (
                    sync_id, dataset_revision, previous_dataset_revision,
                    record_key, gene, change_type, previous_payload_json,
                    current_payload_json, detected_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                change_rows,
            )
            conn.execute("DELETE FROM review_gold_records")
            conn.execute("DELETE FROM review_gold_overlays")
            conn.executemany(
                """
                INSERT INTO review_gold_records (
                    record_key, gene, pmid, source_notation, status, revision,
                    source_reviewer_user_id, decided_by_user_id, payload_json,
                    source_dataset_revision, synced_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        row["record_key"],
                        row["gene"],
                        row["pmid"],
                        row["source_notation"],
                        row["status"],
                        row["revision"],
                        row["source_reviewer_user_id"],
                        row["decided_by_user_id"],
                        row["payload_json"],
                        metadata["dataset_revision"],
                        synced_at,
                    )
                    for row in prepared_gold
                ],
            )
            conn.executemany(
                """
                INSERT INTO review_gold_overlays (
                    record_key, gene, pmid, source_notation, action,
                    match_status, payload_json, source_dataset_revision, synced_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        row["record_key"],
                        row["gene"],
                        row["pmid"],
                        row["source_notation"],
                        row["action"],
                        row["match_status"],
                        _canonical_json(row),
                        metadata["dataset_revision"],
                        synced_at,
                    )
                    for row in overlay_rows
                ],
            )
            conn.execute(
                """
                INSERT INTO review_gold_sync_state (
                    singleton, source_url, source, schema_version,
                    dataset_revision, record_count, generated_at, synced_at
                ) VALUES (1, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(singleton) DO UPDATE SET
                    source_url=excluded.source_url,
                    source=excluded.source,
                    schema_version=excluded.schema_version,
                    dataset_revision=excluded.dataset_revision,
                    record_count=excluded.record_count,
                    generated_at=excluded.generated_at,
                    synced_at=excluded.synced_at
                """,
                (
                    metadata["source_url"],
                    metadata["source"],
                    metadata["schema_version"],
                    metadata["dataset_revision"],
                    metadata["record_count"],
                    metadata["generated_at"],
                    synced_at,
                ),
            )
            tier_counts = {
                name: _count_tier_records_conn(conn, name)
                for (name,) in conn.execute(
                    "SELECT name FROM review_gold_tiers ORDER BY name"
                )
            }
            excluded_count = conn.execute(
                "SELECT COUNT(*) FROM review_gold_exclusions WHERE active = 1"
            ).fetchone()[0]
            conn.execute("PRAGMA user_version = 2")
    finally:
        conn.close()
    return path, {
        "sync_id": sync_id,
        "previous_dataset_revision": previous_revision,
        "added_count": len(added),
        "updated_count": len(updated),
        "removed_count": len(removed),
        "cardiac_record_count": cardiac_count,
        "noncardiac_record_count": len(prepared_gold) - cardiac_count,
        "tier_counts": tier_counts,
        "excluded_count": excluded_count,
    }


def read_live_sync_state(path: Path) -> dict[str, Any]:
    """Read the non-secret audit manifest from an existing live cache."""
    if not path.exists():
        return {}
    try:
        conn = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    except sqlite3.Error:
        return {}
    conn.row_factory = sqlite3.Row
    try:
        row = conn.execute(
            "SELECT * FROM review_gold_sync_state WHERE singleton = 1"
        ).fetchone()
    except sqlite3.Error:
        return {}
    finally:
        conn.close()
    return dict(row) if row is not None else {}


def read_gold_tier_definition(path: Optional[Path], tier: str) -> dict[str, Any]:
    """Resolve a built-in or cache-defined tier without mutating the cache."""
    name = str(tier or "").strip().lower()
    builtin = BUILTIN_GOLD_TIERS.get(name)
    if path is None or not path.exists():
        if builtin is None:
            raise GoldSyncError(
                f"Review-gold tier {name!r} needs an existing cache definition."
            )
        return {"name": name, **builtin, "builtin": True}
    try:
        conn = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        raise GoldSyncError(f"Could not open review-gold cache: {exc}") from exc
    try:
        try:
            return _tier_definition_conn(conn, name)
        except sqlite3.Error:
            if builtin is not None:
                return {"name": name, **builtin, "builtin": True}
            raise GoldSyncError(f"Unknown review-gold tier: {name!r}.") from None
    finally:
        conn.close()


def gold_tier_includes_gene(path: Optional[Path], tier: str, gene: str) -> bool:
    return _tier_includes_gene(read_gold_tier_definition(path, tier), gene)


def current_gold_tier_counts(path: Path) -> dict[str, int]:
    """Return active record counts for every tier in a version-2 cache."""
    try:
        conn = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        raise GoldSyncError(f"Could not open review-gold cache: {exc}") from exc
    try:
        names = [row[0] for row in conn.execute("SELECT name FROM review_gold_tiers")]
        return {name: _count_tier_records_conn(conn, name) for name in names}
    except sqlite3.Error as exc:
        raise GoldSyncError(f"Review-gold cache has no tier metadata: {exc}") from exc
    finally:
        conn.close()


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
                "record_key": row["record_key"],
                "gene": row["gene"],
                "pmid": row["pmid"],
                "source_notation": row["source_notation"],
                "variant_label": row["variant_label"],
                "status": row["status"],
                "revision": row["revision"],
                "verdict": row["verdict"],
                "action": row["action"],
                "match_status": row["match_status"],
                "reason": reason,
                "comment": row["comment"],
                "source_reviewer_user_id": row["source_reviewer_user_id"],
                "source_reviewer": row["source_reviewer"],
                "decided_by_user_id": row["decided_by_user_id"],
                "decided_by": row["decided_by"],
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


def _load_db_aggregates_for_rows(
    export_rows: list[dict[str, Any]],
    *,
    db_overrides: dict[str, Path],
    results_dir: Path,
    no_db: bool,
    announce: bool = True,
) -> tuple[dict[str, dict[tuple[str, str], dict[str, Any]]], set[str]]:
    genes = sorted(
        {str(r.get("gene") or "").strip().upper() for r in export_rows} - {""}
    )
    db_aggregates: dict[str, dict[tuple[str, str], dict[str, Any]]] = {}
    have_db: set[str] = set()
    if no_db:
        return db_aggregates, have_db
    for gene in genes:
        db_path = db_overrides.get(gene) or find_latest_db(gene, results_dir)
        if db_path and db_path.exists():
            agg = load_db_aggregate(db_path)
            db_aggregates[gene] = agg
            have_db.add(gene)
            if announce:
                print(
                    f"  {gene}: matched against {db_path} ({len(agg)} extracted pairs)"
                )
        elif announce:
            print(
                f"  {gene}: no run DB found — overlay rows will be match_status=no_db"
            )
    return db_aggregates, have_db


def sync_live_gold(
    *,
    source_url: str,
    token: str,
    cache_db: Path = DEFAULT_CACHE_DB,
    out_dir: Optional[Path] = None,
    db_overrides: Optional[dict[str, Path]] = None,
    results_dir: Optional[Path] = None,
    no_db: bool = False,
    timeout_s: int = 30,
    announce: bool = True,
    opener: Any = None,
) -> dict[str, Any]:
    """Pull, validate, match, and atomically cache current Azure review gold."""
    out_dir = (out_dir or cache_db.parent).expanduser()
    results_dir = (results_dir or (REPO / "results")).expanduser()
    export_rows, metadata = fetch_live_gold(
        source_url,
        token,
        timeout_s=timeout_s,
        opener=opener,
    )
    db_aggregates, have_db = _load_db_aggregates_for_rows(
        export_rows,
        db_overrides=db_overrides or {},
        results_dir=results_dir,
        no_db=no_db,
        announce=announce,
    )
    overlay_rows = build_overlay_rows(export_rows, db_aggregates, have_db)
    try:
        cache_path, audit = write_live_gold_cache(
            cache_db,
            export_rows,
            overlay_rows,
            metadata,
        )
        followup_path = write_followup_queue(overlay_rows, out_dir)
        summary = summarize(overlay_rows)
        summary_path = out_dir / "review_adjudications_summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2) + "\n",
            encoding="utf-8",
        )
    except (OSError, sqlite3.Error, TypeError, ValueError) as exc:
        raise GoldSyncError(f"Failed to update GVF's live gold cache: {exc}") from exc
    state = {
        **metadata,
        **audit,
        "cache_db": str(cache_path),
        "followup_queue": str(followup_path),
        "summary_path": str(summary_path),
        "summary": summary,
    }
    if announce:
        print(
            f"Synced {metadata['record_count']} lead-approved gold record(s) "
            f"from Azure revision {metadata['dataset_revision'][:12]} -> {cache_path}"
        )
        print(
            "  tiers: "
            + ", ".join(
                f"{name}={count}" for name, count in audit["tier_counts"].items()
            )
            + f"; excluded={audit['excluded_count']}"
        )
        print(
            "  changes from prior current snapshot: "
            f"added={audit['added_count']}, updated={audit['updated_count']}, "
            f"removed={audit['removed_count']}"
        )
    return state


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--export-csv",
        type=Path,
        help=(
            "Offline compatibility input: lead-approved CSV from Variant_Browser "
            "export_gold_standard. Normal operation should use --source-url."
        ),
    )
    source.add_argument(
        "--source-url",
        nargs="?",
        const=DEFAULT_LIVE_GOLD_URL,
        help=(
            "Pull lead-approved gold directly from the Azure-backed Variant Browser "
            f"API (default when flag has no value: {DEFAULT_LIVE_GOLD_URL})."
        ),
    )
    ap.add_argument(
        "--source-token-env",
        default="GVF_REVIEW_GOLD_TOKEN",
        help="Environment variable containing the bearer token (default: %(default)s).",
    )
    ap.add_argument(
        "--cache-db",
        type=Path,
        default=DEFAULT_CACHE_DB,
        help="Live gold SQLite cache used by scoring (default: %(default)s).",
    )
    ap.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Live API timeout in seconds (default: %(default)s).",
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

    overrides = parse_db_overrides(args.db)
    if args.source_url:
        token = os.environ.get(args.source_token_env, "")
        try:
            sync_live_gold(
                source_url=args.source_url,
                token=token,
                cache_db=args.cache_db,
                out_dir=args.out_dir,
                db_overrides=overrides,
                results_dir=args.results_dir,
                no_db=args.no_db,
                timeout_s=args.timeout,
            )
        except GoldSyncError as exc:
            raise SystemExit(str(exc)) from exc
        return 0

    export_path = args.export_csv.expanduser()
    if not export_path.exists():
        raise SystemExit(f"export CSV not found: {export_path}")
    try:
        export_rows = _read_export(export_path)
    except GoldSyncError as exc:
        raise SystemExit(str(exc)) from exc
    if not export_rows:
        print(f"No adjudication rows in {export_path}; nothing to ingest.")
        return 0

    db_aggregates, have_db = _load_db_aggregates_for_rows(
        export_rows,
        db_overrides=overrides,
        results_dir=args.results_dir,
        no_db=args.no_db,
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
