#!/usr/bin/env python3
"""Refresh an existing GVF run after extractor/source-handling fixes.

This is the safe alternative to patching SQLite rows in place:

1. Select PMIDs whose extraction JSON is stale or under-counted relative to
   the current source files and deterministic parsers.
2. Re-extract only those PMIDs from source markdown into canonical
   ``<GENE>_PMID_<PMID>.json`` files, backing up replaced JSONs.
3. Rebuild a fresh SQLite DB from the complete extraction directory.
4. Optionally rerun recovery layers against the fresh DB. Gold is optional; if
   no gold exists, recovery still runs in DB-PMID mode and scoring is skipped.

The script is gene-agnostic and does not use gold standards to choose PMIDs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from harvesting.migrate_to_sqlite import (  # noqa: E402
    create_database_schema,
    migrate_extraction_directory,
)
from pipeline.extraction import ExpertExtractor  # noqa: E402
from pipeline.source_quality import is_usable_fulltext_source  # noqa: E402
from utils.models import Paper  # noqa: E402

logger = logging.getLogger("refresh_run_db")


@dataclass
class ReplayCandidate:
    pmid: str
    source_file: Path
    output_file: Path
    current_variants: int
    deterministic_variants: int
    reasons: list[str] = field(default_factory=list)


def _json_load(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _variant_count(data: dict[str, Any]) -> int:
    variants = data.get("variants")
    if isinstance(variants, list):
        return len(variants)
    meta_count = (data.get("extraction_metadata") or {}).get("total_variants_found")
    try:
        return int(meta_count)
    except (TypeError, ValueError):
        return 0


def _metadata_mentions_abstract_only(metadata: dict[str, Any]) -> bool:
    fields: list[str] = []
    for key in ("source_type", "source_kind", "notes"):
        value = metadata.get(key)
        if value:
            fields.append(str(value))
    challenges = metadata.get("challenges")
    if isinstance(challenges, list):
        fields.extend(str(item) for item in challenges)
    elif challenges:
        fields.append(str(challenges))
    joined = " ".join(fields).lower()
    return (
        bool(metadata.get("abstract_only"))
        or "abstract-only" in joined
        or "abstract only" in joined
        or "full text not available" in joined
        or "full text could not be retrieved" in joined
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_metadata(source_file: Path) -> dict[str, Any]:
    stat = source_file.stat()
    return {
        "source_file": str(source_file),
        "source_size_bytes": stat.st_size,
        "source_sha256": _sha256(source_file),
        "source_type": "fulltext",
        "abstract_only": False,
    }


def discover_source_files(harvest_dir: Path) -> dict[str, Path]:
    """Return PMID -> preferred extraction source, matching pipeline priority."""
    data_zones = {
        f.name.replace("_DATA_ZONES.md", ""): f
        for f in harvest_dir.glob("*_DATA_ZONES.md")
        if is_usable_fulltext_source(f)
    }
    cleaned = {
        f.name.replace("_CLEANED.md", ""): f
        for f in harvest_dir.glob("*_CLEANED.md")
        if is_usable_fulltext_source(f)
    }
    full_context = {
        f.name.replace("_FULL_CONTEXT.md", ""): f
        for f in harvest_dir.glob("*_FULL_CONTEXT.md")
        if is_usable_fulltext_source(f)
    }

    sources: dict[str, Path] = {}
    for pmid in sorted(set(data_zones) | set(cleaned) | set(full_context)):
        if pmid in data_zones:
            sources[pmid] = data_zones[pmid]
        elif pmid in cleaned:
            sources[pmid] = cleaned[pmid]
        else:
            sources[pmid] = full_context[pmid]
    return sources


def _dedup_count(variants: list[dict[str, Any]]) -> int:
    keys = {
        (
            str(v.get("cdna_notation") or "").lower(),
            str(v.get("protein_notation") or "").lower(),
            str(v.get("genomic_position") or "").lower(),
        )
        for v in variants
        if v.get("cdna_notation")
        or v.get("protein_notation")
        or v.get("genomic_position")
    }
    return len(keys)


def deterministic_variant_count(
    extractor: ExpertExtractor,
    source_file: Path,
    gene: str,
) -> int:
    if not is_usable_fulltext_source(source_file):
        return 0
    text = source_file.read_text(encoding="utf-8", errors="replace")
    variants: list[dict[str, Any]] = []
    variants.extend(extractor._parse_markdown_table_variants(text, gene))
    variants.extend(extractor._parse_fixed_width_table_variants(text, gene))
    variants.extend(extractor._parse_vertical_gene_table_variants(text, gene))
    return _dedup_count(variants)


def select_replay_candidates(
    *,
    gene: str,
    harvest_dir: Path,
    extraction_dir: Path,
    min_deterministic_variants: int,
    min_deterministic_lift: int,
    deterministic_lift_ratio: float,
    include_source_newer: bool,
    replay_missing_fingerprint: bool,
) -> list[ReplayCandidate]:
    sources = discover_source_files(harvest_dir)
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    candidates: list[ReplayCandidate] = []

    for pmid, source_file in sources.items():
        output_file = extraction_dir / f"{gene}_PMID_{pmid}.json"
        data = _json_load(output_file) if output_file.exists() else {}
        metadata = data.get("extraction_metadata") or {}
        current_count = _variant_count(data)
        deterministic_count = deterministic_variant_count(extractor, source_file, gene)

        reasons: list[str] = []
        if not output_file.exists():
            reasons.append("missing_extraction")
        if metadata and _metadata_mentions_abstract_only(metadata):
            reasons.append("stale_abstract_only")

        existing_sha = metadata.get("source_sha256")
        if existing_sha and existing_sha != _sha256(source_file):
            reasons.append("source_fingerprint_mismatch")
        elif replay_missing_fingerprint and output_file.exists() and not existing_sha:
            reasons.append("missing_source_fingerprint")

        if include_source_newer and output_file.exists():
            if source_file.stat().st_mtime > output_file.stat().st_mtime + 1:
                reasons.append("source_newer_than_extraction")

        deterministic_lift = deterministic_count - current_count
        ratio_ok = (
            deterministic_count >= current_count * deterministic_lift_ratio
            if current_count
            else deterministic_count >= min_deterministic_variants
        )
        if (
            deterministic_count >= min_deterministic_variants
            and deterministic_lift >= min_deterministic_lift
            and ratio_ok
        ):
            reasons.append(
                f"deterministic_parser_lift:{deterministic_count}>{current_count}"
            )

        if reasons:
            candidates.append(
                ReplayCandidate(
                    pmid=pmid,
                    source_file=source_file,
                    output_file=output_file,
                    current_variants=current_count,
                    deterministic_variants=deterministic_count,
                    reasons=reasons,
                )
            )
    return candidates


def write_candidates_csv(candidates: list[ReplayCandidate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "pmid",
                "source_file",
                "output_file",
                "current_variants",
                "deterministic_variants",
                "reasons",
            ],
        )
        writer.writeheader()
        for c in candidates:
            writer.writerow(
                {
                    "pmid": c.pmid,
                    "source_file": str(c.source_file),
                    "output_file": str(c.output_file),
                    "current_variants": c.current_variants,
                    "deterministic_variants": c.deterministic_variants,
                    "reasons": ";".join(c.reasons),
                }
            )


def replay_candidates(
    *,
    candidates: list[ReplayCandidate],
    gene: str,
    harvest_dir: Path,
    backup_dir: Path,
    tier_threshold: int,
    dry_run: bool,
) -> dict[str, Any]:
    if dry_run:
        return {"attempted": 0, "successful": 0, "failed": 0, "errors": []}

    backup_dir.mkdir(parents=True, exist_ok=True)
    extractor = ExpertExtractor(
        tier_threshold=tier_threshold, fulltext_dir=str(harvest_dir)
    )
    attempted = 0
    successful = 0
    errors: list[dict[str, str]] = []

    for candidate in candidates:
        attempted += 1
        backup_file = backup_dir / candidate.output_file.name
        if candidate.output_file.exists() and not backup_file.exists():
            shutil.copy2(candidate.output_file, backup_file)

        try:
            if not is_usable_fulltext_source(candidate.source_file):
                raise RuntimeError("source is abstract-only fallback or empty")

            text = candidate.source_file.read_text(encoding="utf-8", errors="replace")
            result = extractor.extract(
                Paper(pmid=candidate.pmid, full_text=text, gene_symbol=gene)
            )
            if not result.success or not result.extracted_data:
                raise RuntimeError(result.error or "empty extraction result")
            metadata = result.extracted_data.setdefault("extraction_metadata", {})
            metadata.update(_source_metadata(candidate.source_file))
            metadata["model_used"] = result.model_used
            candidate.output_file.write_text(
                json.dumps(result.extracted_data, indent=2),
                encoding="utf-8",
            )
            successful += 1
            logger.info(
                "replayed PMID %s: %s variants via %s",
                candidate.pmid,
                len(result.extracted_data.get("variants", []) or []),
                result.model_used,
            )
        except Exception as exc:  # noqa: BLE001
            errors.append({"pmid": candidate.pmid, "error": str(exc)})
            logger.warning("PMID %s replay failed: %s", candidate.pmid, exc)
            if backup_file.exists():
                shutil.copy2(backup_file, candidate.output_file)

    return {
        "attempted": attempted,
        "successful": successful,
        "failed": len(errors),
        "errors": errors,
        "backup_dir": str(backup_dir),
    }


def rebuild_db(
    extraction_dir: Path, output_db: Path, *, dry_run: bool
) -> dict[str, Any]:
    if dry_run:
        return {"output_db": str(output_db), "skipped": True}
    if output_db.exists():
        raise FileExistsError(f"Output DB already exists: {output_db}")
    conn = create_database_schema(str(output_db))
    try:
        stats = migrate_extraction_directory(conn, extraction_dir)
    finally:
        conn.close()
    stats["output_db"] = str(output_db)
    return stats


def run_recovery_layers(
    *,
    gene: str,
    db: Path,
    run_dir: Path,
    gold: Optional[Path],
    outdir: Path,
    skip_layers: list[str],
    dry_run: bool,
) -> Optional[dict[str, Any]]:
    if dry_run:
        return {"skipped": True, "outdir": str(outdir)}
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "run_all_layers.py"),
        "--gene",
        gene,
        "--db",
        str(db),
        "--pmc-dir",
        str(run_dir / "pmc_fulltext"),
        "--outdir",
        str(outdir),
    ]
    if gold:
        cmd.extend(["--gold", str(gold)])
    for layer in skip_layers:
        cmd.extend(["--skip", layer])
    logger.info("recovery → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr[-1000:] or result.stdout[-1000:])
    summary_path = outdir / "progression.json"
    if summary_path.exists():
        return json.loads(summary_path.read_text(encoding="utf-8"))
    return {"outdir": str(outdir), "stdout_tail": result.stdout[-1000:]}


def _default_gold(gene: str) -> Optional[Path]:
    path = (
        REPO_ROOT
        / "gene_variant_fetcher_gold_standard"
        / "normalized"
        / f"{gene}_recall_input.csv"
    )
    return path if path.exists() else None


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene", required=True, help="Gene symbol")
    p.add_argument("--run-dir", required=True, type=Path, help="Existing GVF run dir")
    p.add_argument("--harvest-dir", type=Path, default=None)
    p.add_argument("--extraction-dir", type=Path, default=None)
    p.add_argument("--output-db", type=Path, default=None)
    p.add_argument(
        "--replace-db",
        action="store_true",
        help="After refresh/recovery, back up <run-dir>/<GENE>.db and replace it.",
    )
    p.add_argument("--gold", type=Path, default=None)
    p.add_argument(
        "--layers-outdir",
        type=Path,
        default=None,
        help=(
            "Directory for recovery-layer outputs. Default: "
            "<refresh_dir>/layers, so each refresh keeps independent metrics."
        ),
    )
    p.add_argument("--skip-recovery", action="store_true")
    p.add_argument(
        "--skip-layer",
        action="append",
        default=[],
        choices=("clinvar", "pubtator", "figures"),
        help="Recovery layer to skip. Repeatable.",
    )
    p.add_argument("--tier-threshold", type=int, default=1)
    p.add_argument("--min-deterministic-variants", type=int, default=20)
    p.add_argument("--min-deterministic-lift", type=int, default=5)
    p.add_argument("--deterministic-lift-ratio", type=float, default=1.2)
    p.add_argument("--include-source-newer", action="store_true")
    p.add_argument("--replay-missing-fingerprint", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    gene = args.gene.upper()
    run_dir = args.run_dir.expanduser().resolve()
    harvest_dir = (args.harvest_dir or run_dir / "pmc_fulltext").expanduser()
    extraction_dir = (args.extraction_dir or run_dir / "extractions").expanduser()
    gold = args.gold.expanduser() if args.gold else _default_gold(gene)

    if not run_dir.is_dir():
        sys.exit(f"Run dir not found: {run_dir}")
    if not harvest_dir.is_dir():
        sys.exit(f"Harvest dir not found: {harvest_dir}")
    if not extraction_dir.is_dir():
        sys.exit(f"Extraction dir not found: {extraction_dir}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    refresh_dir = run_dir / f"refresh_{timestamp}"
    refresh_dir.mkdir(parents=True, exist_ok=True)

    candidates = select_replay_candidates(
        gene=gene,
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=args.min_deterministic_variants,
        min_deterministic_lift=args.min_deterministic_lift,
        deterministic_lift_ratio=args.deterministic_lift_ratio,
        include_source_newer=args.include_source_newer,
        replay_missing_fingerprint=args.replay_missing_fingerprint,
    )
    candidates_csv = refresh_dir / "replay_candidates.csv"
    write_candidates_csv(candidates, candidates_csv)
    logger.info("Selected %d replay candidates → %s", len(candidates), candidates_csv)

    replay_stats = replay_candidates(
        candidates=candidates,
        gene=gene,
        harvest_dir=harvest_dir,
        backup_dir=refresh_dir / "extraction_json_backup",
        tier_threshold=args.tier_threshold,
        dry_run=args.dry_run,
    )

    output_db = args.output_db
    if output_db is None:
        output_db = run_dir / f"{gene}.refresh_{timestamp}.db"
    output_db = output_db.expanduser()
    db_stats = rebuild_db(extraction_dir, output_db, dry_run=args.dry_run)

    recovery_summary = None
    if not args.skip_recovery:
        layers_outdir = (args.layers_outdir or refresh_dir / "layers").expanduser()
        recovery_summary = run_recovery_layers(
            gene=gene,
            db=output_db,
            run_dir=run_dir,
            gold=gold,
            outdir=layers_outdir,
            skip_layers=args.skip_layer,
            dry_run=args.dry_run,
        )

    replace_info = None
    active_db = output_db
    current_db = run_dir / f"{gene}.db"
    if args.replace_db and not args.dry_run:
        backup_db = run_dir / f"{gene}.db.before_refresh_{timestamp}.db"
        if current_db.exists():
            shutil.copy2(current_db, backup_db)
        shutil.copy2(output_db, current_db)
        active_db = current_db
        replace_info = {"backup_db": str(backup_db), "active_db": str(current_db)}

    summary = {
        "gene": gene,
        "run_dir": str(run_dir),
        "harvest_dir": str(harvest_dir),
        "extraction_dir": str(extraction_dir),
        "gold": str(gold) if gold else None,
        "candidate_count": len(candidates),
        "candidates_csv": str(candidates_csv),
        "replay": replay_stats,
        "db_rebuild": db_stats,
        "recovery": recovery_summary,
        "replace": replace_info,
        "active_db": str(active_db),
    }
    summary_path = refresh_dir / "refresh_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
