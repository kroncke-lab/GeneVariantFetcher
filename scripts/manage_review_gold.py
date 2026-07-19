#!/usr/bin/env python3
"""Inspect and manage versioned Azure review-gold snapshots and tiers.

This tool never deletes immutable source snapshots. ``exclude`` removes records
from active metric tiers through an audited, reversible local ledger; ``restore``
re-enables them. Upstream approval changes remain owned by Variant Browser.
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.ingest_review_adjudications import (  # noqa: E402
    BUILTIN_GOLD_TIERS,
    DEFAULT_CACHE_DB,
    GoldSyncError,
    _canonical_json,
    _count_tier_records_conn,
    _tier_definition_conn,
    _tier_includes_gene,
    ensure_live_gold_schema,
)

TIER_NAME_RE = re.compile(r"^[a-z][a-z0-9_-]{0,31}$")


def _connect(path: Path, *, write: bool = False) -> sqlite3.Connection:
    path = path.expanduser()
    if not path.exists():
        raise GoldSyncError(f"Review-gold cache not found: {path}")
    if write:
        conn = sqlite3.connect(path)
    else:
        conn = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def _print_rows(rows: list[dict[str, Any]], columns: list[str]) -> None:
    if not rows:
        print("No rows.")
        return
    widths = {
        column: max(len(column), *(len(str(row.get(column, ""))) for row in rows))
        for column in columns
    }
    print("  ".join(column.ljust(widths[column]) for column in columns))
    print("  ".join("-" * widths[column] for column in columns))
    for row in rows:
        print(
            "  ".join(
                str(row.get(column, "")).ljust(widths[column]) for column in columns
            )
        )


def summary(path: Path) -> int:
    conn = _connect(path)
    try:
        state = conn.execute(
            "SELECT * FROM review_gold_sync_state WHERE singleton = 1"
        ).fetchone()
        if state is None:
            raise GoldSyncError("Review-gold cache has no current sync state.")
        snapshot_count = conn.execute(
            "SELECT COUNT(*) FROM review_gold_snapshots"
        ).fetchone()[0]
        sync_count = conn.execute(
            "SELECT COUNT(*) FROM review_gold_sync_runs"
        ).fetchone()[0]
        excluded_count = conn.execute(
            "SELECT COUNT(*) FROM review_gold_exclusions WHERE active = 1"
        ).fetchone()[0]
        tier_rows = []
        for tier_name, description, builtin in conn.execute(
            "SELECT name, description, builtin FROM review_gold_tiers ORDER BY name"
        ):
            tier_rows.append(
                {
                    "tier": tier_name,
                    "active_records": _count_tier_records_conn(conn, tier_name),
                    "builtin": "yes" if builtin else "no",
                    "description": description,
                }
            )
    finally:
        conn.close()
    print(f"Current revision: {state['dataset_revision']}")
    print(f"Generated in Azure: {state['generated_at']}")
    print(f"Last synced locally: {state['synced_at']}")
    print(f"Current source records: {state['record_count']}")
    print(f"Immutable snapshots: {snapshot_count}; recorded sync runs: {sync_count}")
    print(f"Active manual exclusions: {excluded_count}")
    print()
    _print_rows(tier_rows, ["tier", "active_records", "builtin", "description"])
    return 0


def list_snapshots(path: Path) -> int:
    conn = _connect(path)
    try:
        rows = [
            {
                "revision": row["dataset_revision"][:12],
                "records": row["record_count"],
                "cardiac": row["cardiac_record_count"],
                "noncardiac": row["noncardiac_record_count"],
                "pulls": row["sync_count"],
                "generated_at": row["generated_at"],
                "last_synced_at": row["last_synced_at"],
            }
            for row in conn.execute(
                "SELECT * FROM review_gold_snapshots ORDER BY first_synced_at DESC"
            )
        ]
    finally:
        conn.close()
    _print_rows(
        rows,
        [
            "revision",
            "records",
            "cardiac",
            "noncardiac",
            "pulls",
            "generated_at",
            "last_synced_at",
        ],
    )
    return 0


def list_syncs(path: Path, limit: int) -> int:
    conn = _connect(path)
    try:
        rows = [
            {
                "sync_id": row["sync_id"],
                "revision": row["dataset_revision"][:12],
                "records": row["record_count"],
                "cardiac": row["cardiac_record_count"],
                "noncardiac": row["noncardiac_record_count"],
                "added": row["added_count"],
                "updated": row["updated_count"],
                "removed": row["removed_count"],
                "synced_at": row["synced_at"],
            }
            for row in conn.execute(
                "SELECT * FROM review_gold_sync_runs ORDER BY sync_id DESC LIMIT ?",
                (limit,),
            )
        ]
    finally:
        conn.close()
    _print_rows(
        rows,
        [
            "sync_id",
            "revision",
            "records",
            "cardiac",
            "noncardiac",
            "added",
            "updated",
            "removed",
            "synced_at",
        ],
    )
    return 0


def list_changes(path: Path, limit: int) -> int:
    conn = _connect(path)
    try:
        rows = [
            {
                "sync_id": row["sync_id"],
                "change": row["change_type"],
                "gene": row["gene"],
                "record_key": row["record_key"],
                "detected_at": row["detected_at"],
            }
            for row in conn.execute(
                "SELECT * FROM review_gold_record_changes "
                "ORDER BY sync_id DESC, record_key LIMIT ?",
                (limit,),
            )
        ]
    finally:
        conn.close()
    _print_rows(rows, ["sync_id", "change", "gene", "record_key", "detected_at"])
    return 0


def list_records(
    path: Path,
    *,
    tier: str,
    snapshot: str,
    include_excluded: bool,
) -> int:
    conn = _connect(path)
    try:
        definition = _tier_definition_conn(conn, tier)
        excluded = {
            row[0]
            for row in conn.execute(
                "SELECT record_key FROM review_gold_exclusions WHERE active = 1"
            )
        }
        if snapshot == "current":
            query = (
                "SELECT record_key, gene, pmid, source_notation, status, revision "
                "FROM review_gold_records ORDER BY gene, pmid, source_notation"
            )
            source_rows = conn.execute(query)
        else:
            matches = conn.execute(
                "SELECT dataset_revision FROM review_gold_snapshots "
                "WHERE dataset_revision LIKE ? ORDER BY first_synced_at DESC",
                (f"{snapshot}%",),
            ).fetchall()
            if len(matches) != 1:
                raise GoldSyncError(
                    f"Snapshot prefix {snapshot!r} matched {len(matches)} revisions."
                )
            source_rows = conn.execute(
                "SELECT record_key, gene, pmid, source_notation, status, revision "
                "FROM review_gold_snapshot_records WHERE dataset_revision = ? "
                "ORDER BY gene, pmid, source_notation",
                (matches[0][0],),
            )
        rows = []
        for row in source_rows:
            if not _tier_includes_gene(definition, row[1]):
                continue
            is_excluded = row[0] in excluded
            if is_excluded and not include_excluded:
                continue
            rows.append(
                {
                    "record_key": row[0],
                    "gene": row[1],
                    "pmid": row[2],
                    "source_notation": row[3],
                    "status": row[4],
                    "revision": row[5],
                    "excluded": "yes" if is_excluded else "no",
                }
            )
    finally:
        conn.close()
    _print_rows(
        rows,
        [
            "record_key",
            "gene",
            "pmid",
            "source_notation",
            "status",
            "revision",
            "excluded",
        ],
    )
    return 0


def _target_record_keys(
    conn: sqlite3.Connection,
    *,
    record_keys: list[str],
    gene: str,
) -> list[str]:
    targets = {key.strip() for key in record_keys if key.strip()}
    if gene:
        targets.update(
            row[0]
            for row in conn.execute(
                "SELECT record_key FROM review_gold_records WHERE gene = ?",
                (gene.strip().upper(),),
            )
        )
    if not targets:
        raise GoldSyncError("Provide at least one --record-key or --gene.")
    known = {
        row[0]
        for row in conn.execute(
            "SELECT DISTINCT record_key FROM review_gold_snapshot_records "
            f"WHERE record_key IN ({','.join('?' for _ in targets)})",
            tuple(sorted(targets)),
        )
    }
    unknown = sorted(targets - known)
    if unknown:
        raise GoldSyncError(f"Unknown review-gold record key(s): {', '.join(unknown)}")
    return sorted(targets)


def set_exclusion(
    path: Path,
    *,
    record_keys: list[str],
    gene: str,
    active: bool,
    actor: str,
    reason: str,
) -> int:
    actor = actor.strip()
    reason = reason.strip()
    if not actor or not reason:
        raise GoldSyncError("--actor and --reason are required for the audit ledger.")
    conn = _connect(path, write=True)
    try:
        with conn:
            ensure_live_gold_schema(conn)
            targets = _target_record_keys(conn, record_keys=record_keys, gene=gene)
            now = datetime.now(timezone.utc).isoformat()
            action = "exclude" if active else "restore"
            conn.executemany(
                """
                INSERT INTO review_gold_exclusions (
                    record_key, active, actor, reason, updated_at
                ) VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(record_key) DO UPDATE SET
                    active=excluded.active,
                    actor=excluded.actor,
                    reason=excluded.reason,
                    updated_at=excluded.updated_at
                """,
                [(key, int(active), actor, reason, now) for key in targets],
            )
            conn.executemany(
                "INSERT INTO review_gold_exclusion_events "
                "(record_key, action, actor, reason, created_at) "
                "VALUES (?, ?, ?, ?, ?)",
                [(key, action, actor, reason, now) for key in targets],
            )
    finally:
        conn.close()
    print(
        f"{action.title()}d {len(targets)} record(s); source snapshots were retained."
    )
    return 0


def _parse_genes(value: str) -> list[str]:
    return sorted(
        {
            item.strip().upper()
            for item in re.split(r"[\s,]+", value or "")
            if item.strip()
        }
    )


def define_tier(
    path: Path,
    *,
    name: str,
    description: str,
    include_genes: str,
    exclude_genes: str,
    all_genes: bool,
    actor: str,
) -> int:
    name = name.strip().lower()
    actor = actor.strip()
    if not TIER_NAME_RE.fullmatch(name):
        raise GoldSyncError(
            "Tier name must start with a letter and contain only a-z, 0-9, _ or -."
        )
    if name in BUILTIN_GOLD_TIERS:
        raise GoldSyncError(f"Built-in tier {name!r} cannot be redefined.")
    choices = sum(bool(value) for value in (include_genes, exclude_genes, all_genes))
    if choices != 1:
        raise GoldSyncError(
            "Choose exactly one of --include-genes, --exclude-genes, or --all-genes."
        )
    if not actor:
        raise GoldSyncError("--actor is required for tier-definition audit metadata.")
    if include_genes:
        mode, genes = "include", _parse_genes(include_genes)
    elif exclude_genes:
        mode, genes = "exclude", _parse_genes(exclude_genes)
    else:
        mode, genes = "all", []
    if mode != "all" and not genes:
        raise GoldSyncError("A gene-based tier must contain at least one gene.")
    conn = _connect(path, write=True)
    try:
        with conn:
            ensure_live_gold_schema(conn)
            now = datetime.now(timezone.utc).isoformat()
            conn.execute(
                """
                INSERT INTO review_gold_tiers (
                    name, description, gene_mode, genes_json, builtin,
                    created_at, updated_at, updated_by
                ) VALUES (?, ?, ?, ?, 0, ?, ?, ?)
                ON CONFLICT(name) DO UPDATE SET
                    description=excluded.description,
                    gene_mode=excluded.gene_mode,
                    genes_json=excluded.genes_json,
                    updated_at=excluded.updated_at,
                    updated_by=excluded.updated_by
                """,
                (
                    name,
                    description.strip(),
                    mode,
                    _canonical_json(genes),
                    now,
                    now,
                    actor,
                ),
            )
    finally:
        conn.close()
    print(f"Defined tier {name!r}: mode={mode}, genes={','.join(genes) or '<all>'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=DEFAULT_CACHE_DB,
        help="Versioned review-gold SQLite cache (default: %(default)s).",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("summary", help="Show current revision, history, and tier counts.")
    sub.add_parser("snapshots", help="List immutable source snapshots.")
    syncs = sub.add_parser("syncs", help="List individual authenticated pulls.")
    syncs.add_argument("--limit", type=int, default=50)
    changes = sub.add_parser("changes", help="List added/updated/removed records.")
    changes.add_argument("--limit", type=int, default=100)
    records = sub.add_parser("records", help="List records in a tier/snapshot.")
    records.add_argument("--tier", default="cardiac")
    records.add_argument("--snapshot", default="current")
    records.add_argument("--include-excluded", action="store_true")
    for command in ("exclude", "restore"):
        action = sub.add_parser(
            command, help=f"{command.title()} active metric records."
        )
        action.add_argument("--record-key", action="append", default=[])
        action.add_argument("--gene", default="")
        action.add_argument("--actor", required=True)
        action.add_argument("--reason", required=True)
    tier = sub.add_parser("define-tier", help="Create/update a custom gene tier.")
    tier.add_argument("name")
    tier.add_argument("--description", default="")
    tier.add_argument("--include-genes", default="")
    tier.add_argument("--exclude-genes", default="")
    tier.add_argument("--all-genes", action="store_true")
    tier.add_argument("--actor", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "summary":
            return summary(args.db)
        if args.command == "snapshots":
            return list_snapshots(args.db)
        if args.command == "syncs":
            return list_syncs(args.db, args.limit)
        if args.command == "changes":
            return list_changes(args.db, args.limit)
        if args.command == "records":
            return list_records(
                args.db,
                tier=args.tier,
                snapshot=args.snapshot,
                include_excluded=args.include_excluded,
            )
        if args.command in {"exclude", "restore"}:
            return set_exclusion(
                args.db,
                record_keys=args.record_key,
                gene=args.gene,
                active=args.command == "exclude",
                actor=args.actor,
                reason=args.reason,
            )
        if args.command == "define-tier":
            return define_tier(
                args.db,
                name=args.name,
                description=args.description,
                include_genes=args.include_genes,
                exclude_genes=args.exclude_genes,
                all_genes=args.all_genes,
                actor=args.actor,
            )
    except (GoldSyncError, sqlite3.Error, json.JSONDecodeError) as exc:
        raise SystemExit(str(exc)) from exc
    raise AssertionError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
