"""
`gvf reharvest` — re-harvest abstract-only / paywall-scrape PMIDs using
`harvesting.html_body_fetcher`.

Scans a `pmc_fulltext/` directory, identifies FULL_CONTEXT.md files that need
re-harvest (abstract-only fallback header, tiny files, paywall sentinels), and
tries a fresh set of OA-friendly sources (Europe PMC, NCBI ELink → PMC,
Unpaywall HTML, Unpaywall PDF). If any source clears the quality gate, the
FULL_CONTEXT.md is overwritten in place and the recovered PMID is logged so
downstream re-extraction knows what to re-run.

Outputs (in --report-dir):
    reharvest_log.csv      — every PMID attempted with per-source outcomes
    recovered_pmids.txt    — PMIDs that gained real full text (feed into extract-folder)
    still_failed_pmids.txt — PMIDs that remain abstract-only
"""

from __future__ import annotations

import csv
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import typer
from typing_extensions import Annotated

logger = logging.getLogger(__name__)


def _iter_pmids_in_dir(pmc_dir: Path) -> list[str]:
    """Return PMIDs (file stems before `_FULL_CONTEXT`)."""
    out = []
    for p in pmc_dir.glob("*_FULL_CONTEXT.md"):
        stem = p.stem  # e.g. "12345_FULL_CONTEXT"
        if stem.endswith("_FULL_CONTEXT"):
            pmid = stem[: -len("_FULL_CONTEXT")]
            if pmid.isdigit():
                out.append(pmid)
    return sorted(out)


def _read_explicit_pmids(pmid_file: Path) -> list[str]:
    pmids = []
    for line in pmid_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # strip trailing comments
        tok = line.split("#", 1)[0].strip()
        if tok.isdigit():
            pmids.append(tok)
    return pmids


def reharvest_command(
    pmc_dir: Annotated[
        Path,
        typer.Argument(
            help="Directory containing `<PMID>_FULL_CONTEXT.md` files (e.g. results/KCNH2/<run>/pmc_fulltext)"
        ),
    ],
    report_dir: Annotated[
        Path,
        typer.Option(
            "--report-dir",
            help="Where to write reharvest_log.csv, recovered_pmids.txt, still_failed_pmids.txt",
        ),
    ],
    email: Annotated[
        Optional[str],
        typer.Option(
            "--email",
            help="Email for NCBI / Unpaywall (defaults to $NCBI_EMAIL).",
        ),
    ] = None,
    pmid_file: Annotated[
        Optional[Path],
        typer.Option(
            "--pmid-file",
            help="Limit to PMIDs listed in this file (one per line). Defaults to every needs-reharvest file in pmc_dir.",
        ),
    ] = None,
    max_pmids: Annotated[
        int,
        typer.Option(
            "--max-pmids",
            help="Stop after this many PMIDs (0 = no cap).",
        ),
    ] = 0,
    workers: Annotated[
        int,
        typer.Option(
            "--workers",
            help="Parallel fetcher threads (be modest — Europe PMC + Unpaywall are shared).",
        ),
    ] = 4,
    overwrite: Annotated[
        bool,
        typer.Option(
            "--overwrite/--no-overwrite",
            help="When True (default), recovered markdown overwrites the FULL_CONTEXT.md in place.",
        ),
    ] = True,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
):
    """Re-harvest abstract-only / paywall-scrape PMIDs using a multi-source HTML fallback."""
    from harvesting.html_body_fetcher import HtmlBodyFetcher, needs_reharvest

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    pmc_dir = pmc_dir.expanduser().resolve()
    if not pmc_dir.is_dir():
        typer.echo(f"⚠️  Not a directory: {pmc_dir}", err=True)
        raise typer.Exit(2)

    report_dir = report_dir.expanduser().resolve()
    report_dir.mkdir(parents=True, exist_ok=True)

    email = email or os.environ.get("NCBI_EMAIL")
    if not email:
        typer.echo("⚠️  No email provided and NCBI_EMAIL not in env.", err=True)
        raise typer.Exit(2)

    ncbi_key = os.environ.get("NCBI_API_KEY")

    # Build the candidate list
    if pmid_file:
        all_pmids = _read_explicit_pmids(pmid_file.expanduser())
    else:
        all_pmids = _iter_pmids_in_dir(pmc_dir)

    # Classify each PMID and keep only those needing re-harvest
    candidates: list[tuple[str, str]] = []  # (pmid, reason)
    ok_already = 0
    for pmid in all_pmids:
        fc = pmc_dir / f"{pmid}_FULL_CONTEXT.md"
        need, reason = needs_reharvest(fc)
        if need:
            candidates.append((pmid, reason))
        else:
            ok_already += 1

    if max_pmids and len(candidates) > max_pmids:
        candidates = candidates[:max_pmids]

    typer.echo(
        f"Found {len(all_pmids)} PMIDs in {pmc_dir}.\n"
        f"  needs re-harvest: {len(candidates)}\n"
        f"  already acceptable: {ok_already}\n"
    )
    if not candidates:
        typer.echo("Nothing to do.")
        return

    fetcher = HtmlBodyFetcher(email=email, ncbi_api_key=ncbi_key)

    # CSV log + recovered/failed lists
    log_path = report_dir / "reharvest_log.csv"
    recovered_path = report_dir / "recovered_pmids.txt"
    failed_path = report_dir / "still_failed_pmids.txt"

    log_fields = [
        "pmid",
        "pre_reason",
        "success",
        "source",
        "n_body_lines",
        "n_chars",
        "attempts",
        "error",
        "elapsed_s",
    ]

    recovered: list[str] = []
    failed: list[str] = []

    def _one(pmid: str, pre_reason: str) -> dict:
        t0 = time.time()
        try:
            res = fetcher.fetch_body(pmid)
        except Exception as exc:
            logger.exception(f"PMID {pmid}: fetcher crashed")
            return {
                "pmid": pmid,
                "pre_reason": pre_reason,
                "success": False,
                "source": "",
                "n_body_lines": 0,
                "n_chars": 0,
                "attempts": "[]",
                "error": f"fetcher crashed: {exc}",
                "elapsed_s": round(time.time() - t0, 2),
            }

        elapsed = round(time.time() - t0, 2)
        row = {
            "pmid": pmid,
            "pre_reason": pre_reason,
            "success": res.success,
            "source": res.source or "",
            "n_body_lines": res.n_body_lines,
            "n_chars": res.n_chars,
            "attempts": "; ".join(f"{s}:{st}" for s, st, _ in (res.attempts or [])),
            "error": res.error or "",
            "elapsed_s": elapsed,
        }

        if res.success and overwrite and res.markdown:
            # Write atomically: tmp + rename
            fc = pmc_dir / f"{pmid}_FULL_CONTEXT.md"
            tmp = fc.with_suffix(".md.new")
            header = (
                f"# REHARVESTED FULL TEXT\n\n"
                f"> Source: {res.source}\n"
                f"> Body lines: {res.n_body_lines}, chars: {res.n_chars}\n"
                f"> Recovered: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n"
            )
            tmp.write_text(header + res.markdown)
            tmp.replace(fc)
            logger.info(
                f"PMID {pmid}: RECOVERED via {res.source} ({res.n_body_lines} lines, {res.n_chars} chars, {elapsed}s)"
            )
        elif res.success:
            logger.info(
                f"PMID {pmid}: recovered ({res.source}, {res.n_chars} chars) — overwrite disabled, not writing"
            )
        else:
            logger.info(f"PMID {pmid}: still failed — {row['attempts']}")

        return row

    rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        future_to_pmid = {
            pool.submit(_one, pmid, reason): pmid for pmid, reason in candidates
        }
        done_count = 0
        total = len(candidates)
        for fut in as_completed(future_to_pmid):
            row = fut.result()
            rows.append(row)
            done_count += 1
            if done_count % 10 == 0 or done_count == total:
                n_ok = sum(1 for r in rows if r["success"])
                typer.echo(f"  [{done_count}/{total}] recovered so far: {n_ok}")
            if row["success"]:
                recovered.append(row["pmid"])
            else:
                failed.append(row["pmid"])

    # Write outputs
    with log_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=log_fields)
        w.writeheader()
        w.writerows(rows)

    recovered_path.write_text(
        "\n".join(sorted(recovered, key=int)) + "\n" if recovered else ""
    )
    failed_path.write_text("\n".join(sorted(failed, key=int)) + "\n" if failed else "")

    typer.echo("")
    typer.echo("=" * 60)
    typer.echo(
        f"Done. Recovered: {len(recovered)} / {len(candidates)} ({len(recovered) / len(candidates) * 100:.1f}%)"
    )
    typer.echo(f"  recovered list: {recovered_path}")
    typer.echo(f"  failed list:    {failed_path}")
    typer.echo(f"  full log:       {log_path}")
