#!/usr/bin/env python3
"""Audit and extend existing FULL_CONTEXT.md artifacts without extraction.

This script is intentionally acquisition-only:

1. Audit target PMIDs for missing, usable, or retry-needed full text.
2. Optionally consolidate usable FULL_CONTEXT.md files already present under
   prior validation/result roots into the selected canonical run directories.
3. Optionally run the harvester download phase for remaining missing or bad
   PMIDs, using the current .env credentials, including ELSEVIER_INSTTOKEN.

It never runs LLM extraction, scoring, recovery layers, or DB migration.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from harvesting.html_body_fetcher import needs_reharvest  # noqa: E402
from utils.bootstrap import initialize_runtime  # noqa: E402

LOG = logging.getLogger("fulltext_acquisition")


@dataclass(frozen=True)
class GeneRun:
    gene: str
    run_dir: Path

    @property
    def pmc_dir(self) -> Path:
        return self.run_dir / "pmc_fulltext"

    @property
    def discovered_pmids_file(self) -> Path:
        return self.run_dir / f"{self.gene}_pmids.txt"

    @property
    def gold_pmids_file(self) -> Path:
        return (
            REPO_ROOT
            / "gene_variant_fetcher_gold_standard"
            / "normalized"
            / f"{self.gene}_recall_input.csv"
        )


DEFAULT_RUNS: tuple[GeneRun, ...] = (
    GeneRun(
        "KCNH2",
        REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/KCNH2/20260518_213938",
    ),
    GeneRun(
        "KCNQ1",
        REPO_ROOT / "validation_runs/20260517_203904/results/KCNQ1/20260517_204424",
    ),
    GeneRun(
        "RYR2",
        REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938",
    ),
    GeneRun(
        "SCN5A",
        REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938",
    ),
)


@dataclass
class PmidStatus:
    pmid: str
    status: str
    reason: str
    context_path: str
    context_bytes: int | None


@dataclass
class CopyAction:
    gene: str
    pmid: str
    action: str
    reason: str
    source_path: str
    destination_path: str
    source_bytes: int
    previous_bytes: int | None
    backup_path: str


def _rel(path: Path | str | None) -> str:
    if not path:
        return ""
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(p)


def _pmid_from_context_path(path: Path) -> str:
    return path.name.removesuffix("_FULL_CONTEXT.md")


def _read_pmid_lines(path: Path) -> list[str]:
    pmids: list[str] = []
    if not path.exists():
        return pmids
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        token = line.split("#", 1)[0].strip()
        if token.isdigit() and token not in pmids:
            pmids.append(token)
    return pmids


def _read_gold_pmids(path: Path) -> list[str]:
    pmids: list[str] = []
    if not path.exists():
        return pmids
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pmid = (row.get("pmid") or row.get("PMID") or "").strip()
            if pmid.isdigit() and pmid not in pmids:
                pmids.append(pmid)
    return pmids


def target_pmids(run: GeneRun, target: str, pmid_file: Path | None) -> list[str]:
    if target == "gold":
        return _read_gold_pmids(run.gold_pmids_file)
    if target == "discovered":
        return _read_pmid_lines(run.discovered_pmids_file)
    if target == "existing":
        return sorted(
            _pmid_from_context_path(path)
            for path in run.pmc_dir.glob("*_FULL_CONTEXT.md")
        )
    if target == "pmid-file":
        if pmid_file is None:
            raise ValueError("--target pmid-file requires --pmid-file")
        return _read_pmid_lines(pmid_file)
    raise ValueError(f"Unknown target: {target}")


def status_for_pmid(pmc_dir: Path, pmid: str) -> PmidStatus:
    path = pmc_dir / f"{pmid}_FULL_CONTEXT.md"
    if not path.exists():
        return PmidStatus(pmid, "missing", "missing", "", None)
    need, reason = needs_reharvest(path)
    status = "needs_retry" if need else "usable"
    return PmidStatus(
        pmid=pmid,
        status=status,
        reason=reason,
        context_path=_rel(path),
        context_bytes=path.stat().st_size,
    )


def audit_run(run: GeneRun, pmids: Iterable[str]) -> dict:
    rows = [status_for_pmid(run.pmc_dir, pmid) for pmid in pmids]
    counts: dict[str, int] = {}
    reasons: dict[str, int] = {}
    for row in rows:
        counts[row.status] = counts.get(row.status, 0) + 1
        reasons[row.reason] = reasons.get(row.reason, 0) + 1
    return {
        "gene": run.gene,
        "run_dir": _rel(run.run_dir),
        "pmc_dir": _rel(run.pmc_dir),
        "target_pmids": len(rows),
        "counts": counts,
        "reasons": reasons,
        "rows": [asdict(row) for row in rows],
    }


def build_context_index(scan_roots: Iterable[Path]) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    for root in scan_roots:
        if not root.exists():
            LOG.warning("scan root does not exist: %s", root)
            continue
        for path in root.glob("**/*_FULL_CONTEXT.md"):
            pmid = _pmid_from_context_path(path)
            if pmid.isdigit():
                index.setdefault(pmid, []).append(path)
    for paths in index.values():
        paths.sort(key=lambda p: p.stat().st_size, reverse=True)
    return index


def choose_usable_source(
    pmid: str, destination: Path, index: dict[str, list[Path]]
) -> Path | None:
    for candidate in index.get(pmid, []):
        if candidate.resolve() == destination.resolve():
            continue
        need, _ = needs_reharvest(candidate)
        if not need:
            return candidate
    return None


def _copy_associated_dirs(source: Path, destination: Path, *, dry_run: bool) -> None:
    pmid = _pmid_from_context_path(source)
    for suffix in ("_supplements", "_figures"):
        src = source.parent / f"{pmid}{suffix}"
        dst = destination.parent / f"{pmid}{suffix}"
        if not src.is_dir() or dst.exists():
            continue
        if not dry_run:
            shutil.copytree(src, dst)


def copy_context(
    *,
    gene: str,
    pmid: str,
    source: Path,
    destination: Path,
    dry_run: bool,
) -> CopyAction:
    previous_bytes = destination.stat().st_size if destination.exists() else None
    backup_path = ""
    action = "would_copy" if dry_run else "copied"

    if not dry_run:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            backup = destination.with_suffix(
                destination.suffix + ".pre_fulltext_acquisition_bak"
            )
            if not backup.exists():
                shutil.copy2(destination, backup)
            backup_path = _rel(backup)
        shutil.copy2(source, destination)
        _copy_associated_dirs(source, destination, dry_run=False)

    return CopyAction(
        gene=gene,
        pmid=pmid,
        action=action,
        reason="usable_prior_context",
        source_path=_rel(source),
        destination_path=_rel(destination),
        source_bytes=source.stat().st_size,
        previous_bytes=previous_bytes,
        backup_path=backup_path,
    )


def consolidate_from_prior(
    run: GeneRun,
    pmids: list[str],
    index: dict[str, list[Path]],
    *,
    dry_run: bool,
) -> list[CopyAction]:
    actions: list[CopyAction] = []
    for pmid in pmids:
        status = status_for_pmid(run.pmc_dir, pmid)
        if status.status == "usable":
            continue
        dst = run.pmc_dir / f"{pmid}_FULL_CONTEXT.md"
        source = choose_usable_source(pmid, dst, index)
        if source is None:
            continue
        actions.append(
            copy_context(
                gene=run.gene,
                pmid=pmid,
                source=source,
                destination=dst,
                dry_run=dry_run,
            )
        )
    return actions


def backup_harvest_logs(pmc_dir: Path) -> list[str]:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backed_up: list[str] = []
    for name in (
        "paywalled_missing.csv",
        "successful_downloads.csv",
        "abstract_only_fallback.csv",
        "manifest.json",
    ):
        path = pmc_dir / name
        if not path.exists():
            continue
        backup = path.with_name(f"{path.name}.pre_fulltext_acquisition_{stamp}.bak")
        shutil.copy2(path, backup)
        backed_up.append(_rel(backup))
    return backed_up


def run_harvest_only(
    run: GeneRun,
    pmids: list[str],
    *,
    delay: float,
    dry_run: bool,
) -> dict:
    selected = [
        pmid for pmid in pmids if status_for_pmid(run.pmc_dir, pmid).status != "usable"
    ]
    if dry_run or not selected:
        return {
            "gene": run.gene,
            "attempted": len(selected),
            "dry_run": dry_run,
            "backups": [],
        }

    backups = backup_harvest_logs(run.pmc_dir)

    from harvesting import PMCHarvester

    harvester = PMCHarvester(output_dir=str(run.pmc_dir), gene_symbol=run.gene)
    started = time.time()
    harvester.harvest(selected, delay=delay, run_scout=False)
    elapsed = round(time.time() - started, 2)
    return {
        "gene": run.gene,
        "attempted": len(selected),
        "dry_run": False,
        "elapsed_s": elapsed,
        "backups": backups,
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genes",
        default="KCNH2,KCNQ1,RYR2,SCN5A",
        help="Comma-separated genes from the built-in recall run map.",
    )
    parser.add_argument(
        "--target",
        choices=["gold", "discovered", "existing", "pmid-file"],
        default="gold",
        help="Which PMIDs to audit/acquire for each gene.",
    )
    parser.add_argument(
        "--pmid-file",
        type=Path,
        help="PMID list used when --target pmid-file.",
    )
    parser.add_argument(
        "--scan-root",
        action="append",
        type=Path,
        default=[REPO_ROOT / "validation_runs", REPO_ROOT / "results"],
        help="Root to scan for reusable prior FULL_CONTEXT.md files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Report directory. Default: recall_metrics/fulltext_acquisition/<timestamp>.",
    )
    parser.add_argument(
        "--consolidate",
        action="store_true",
        help="Copy usable prior FULL_CONTEXT.md files into the selected run dirs.",
    )
    parser.add_argument(
        "--harvest",
        action="store_true",
        help="Run PMCHarvester download-only pass for remaining missing/bad PMIDs.",
    )
    parser.add_argument(
        "--max-pmids",
        type=int,
        default=0,
        help="Cap harvest attempts per gene after consolidation (0 = no cap).",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=2.0,
        help="Delay between PMIDs during --harvest.",
    )
    parser.add_argument(
        "--allow-large",
        action="store_true",
        help="Allow --harvest over more than 500 remaining PMIDs per gene.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report planned actions without copying or harvesting.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    initialize_runtime()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    selected = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
    runs = [run for run in DEFAULT_RUNS if run.gene in selected]
    if not runs:
        print(f"No matching built-in genes for {args.genes}", file=sys.stderr)
        return 2

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = args.outdir or (
        REPO_ROOT / "recall_metrics" / "fulltext_acquisition" / stamp
    )
    outdir.mkdir(parents=True, exist_ok=True)

    index = build_context_index(path.expanduser().resolve() for path in args.scan_root)

    before: list[dict] = []
    after: list[dict] = []
    copy_actions: list[CopyAction] = []
    harvest_reports: list[dict] = []

    for run in runs:
        pmids = target_pmids(run, args.target, args.pmid_file)
        before_audit = audit_run(run, pmids)
        before.append(before_audit)

        if args.consolidate:
            actions = consolidate_from_prior(
                run, pmids, index, dry_run=bool(args.dry_run)
            )
            copy_actions.extend(actions)
            LOG.info(
                "%s: %d prior context(s) %s",
                run.gene,
                len(actions),
                "would copy" if args.dry_run else "copied",
            )

        remaining = [
            row["pmid"]
            for row in audit_run(run, pmids)["rows"]
            if row["status"] != "usable"
        ]
        if args.max_pmids and len(remaining) > args.max_pmids:
            remaining = remaining[: args.max_pmids]

        if args.harvest:
            uncapped_remaining = [
                row["pmid"]
                for row in audit_run(run, pmids)["rows"]
                if row["status"] != "usable"
            ]
            if (
                not args.max_pmids
                and not args.allow_large
                and len(uncapped_remaining) > 500
            ):
                print(
                    f"Refusing large harvest for {run.gene}: "
                    f"{len(uncapped_remaining)} PMIDs remain. Use --max-pmids "
                    "or --allow-large.",
                    file=sys.stderr,
                )
                return 4
            harvest_reports.append(
                run_harvest_only(
                    run,
                    remaining,
                    delay=args.delay,
                    dry_run=bool(args.dry_run),
                )
            )

        after.append(audit_run(run, pmids))

    copy_rows = [asdict(action) for action in copy_actions]
    write_csv(
        outdir / "copy_actions.csv",
        copy_rows,
        [
            "gene",
            "pmid",
            "action",
            "reason",
            "source_path",
            "destination_path",
            "source_bytes",
            "previous_bytes",
            "backup_path",
        ],
    )
    status_rows = []
    for phase, audits in (("before", before), ("after", after)):
        for audit in audits:
            for row in audit["rows"]:
                status_rows.append({"phase": phase, "gene": audit["gene"], **row})
    write_csv(
        outdir / "status.csv",
        status_rows,
        [
            "phase",
            "gene",
            "pmid",
            "status",
            "reason",
            "context_path",
            "context_bytes",
        ],
    )
    summary = {
        "created_at": stamp,
        "target": args.target,
        "genes": [run.gene for run in runs],
        "dry_run": bool(args.dry_run),
        "consolidate": bool(args.consolidate),
        "harvest": bool(args.harvest),
        "max_pmids": args.max_pmids,
        "scan_roots": [_rel(path) for path in args.scan_root],
        "before": [{k: v for k, v in audit.items() if k != "rows"} for audit in before],
        "after": [{k: v for k, v in audit.items() if k != "rows"} for audit in after],
        "copy_actions": copy_rows,
        "harvest_reports": harvest_reports,
        "outputs": {
            "status_csv": _rel(outdir / "status.csv"),
            "copy_actions_csv": _rel(outdir / "copy_actions.csv"),
            "summary_json": _rel(outdir / "summary.json"),
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    print(
        json.dumps({k: v for k, v in summary.items() if k != "copy_actions"}, indent=2)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
