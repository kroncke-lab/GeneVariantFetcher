#!/usr/bin/env python3
"""Recover open-access full text for an explicit PMID batch.

This is a source-acquisition helper, not a gold-standard recovery layer. It
tries PMC, Europe PMC, and Unpaywall routes through ``OARecoveryClient`` and
writes one ``*_FULL_CONTEXT.md`` file per successful PMID.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import Iterable

from dotenv import load_dotenv

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from harvesting.oa_recovery import OARecoveryClient

LOG = logging.getLogger("recover_paywall_oa")


def _env_email() -> str | None:
    return os.environ.get("NCBI_EMAIL") or os.environ.get("ENTREZ_EMAIL")


def _read_rows(
    path: Path,
    *,
    pmid_column: str,
    doi_column: str | None,
    limit: int | None,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pmid = (row.get(pmid_column) or "").strip()
            if not pmid:
                continue
            item = {"pmid": pmid}
            if doi_column:
                doi = (row.get(doi_column) or "").strip()
                if doi:
                    item["doi"] = doi
            rows.append(item)
            if limit and len(rows) >= limit:
                break
    return rows


def _write_context(
    *,
    output_dir: Path,
    pmid: str,
    source: str,
    markdown: str,
    write_cleaned: bool,
) -> tuple[Path, Path | None]:
    output_dir.mkdir(parents=True, exist_ok=True)
    context_path = output_dir / f"{pmid}_FULL_CONTEXT.md"
    context_path.write_text(
        f"# RECOVERED VIA {source}\n\n"
        f"PMID: {pmid}\nSource: {source}\nChars: {len(markdown)}\n\n"
        f"{markdown}",
        encoding="utf-8",
    )
    cleaned_path = None
    if write_cleaned:
        cleaned_path = output_dir / f"{pmid}_CLEANED.md"
        cleaned_path.write_text(markdown, encoding="utf-8")
    return context_path, cleaned_path


def recover_batch(
    rows: Iterable[dict[str, str]],
    *,
    output_dir: Path,
    email: str,
    ncbi_api_key: str | None,
    write_cleaned: bool,
) -> dict:
    client = OARecoveryClient(email=email, ncbi_api_key=ncbi_api_key)
    successes: list[dict] = []
    failures: list[dict] = []
    rows = list(rows)

    LOG.info("Recovering %d PMIDs via OA routes", len(rows))
    for index, row in enumerate(rows, start=1):
        pmid = row["pmid"]
        doi = row.get("doi")
        LOG.info("[%d/%d] PMID %s", index, len(rows), pmid)
        try:
            result = client.recover(pmid, doi=doi)
        except Exception as exc:  # pragma: no cover - defensive CLI guard
            LOG.exception("PMID %s failed with exception", pmid)
            failures.append({"pmid": pmid, "error": str(exc), "attempts": []})
            continue

        if result.success and result.markdown:
            context_path, cleaned_path = _write_context(
                output_dir=output_dir,
                pmid=pmid,
                source=result.source or "unknown",
                markdown=result.markdown,
                write_cleaned=write_cleaned,
            )
            successes.append(
                {
                    "pmid": pmid,
                    "source": result.source,
                    "pmcid": result.pmcid,
                    "chars": result.n_chars,
                    "context_path": str(context_path),
                    "cleaned_path": str(cleaned_path) if cleaned_path else None,
                }
            )
            continue

        failures.append(
            {
                "pmid": pmid,
                "error": result.error or "no acceptable source",
                "attempts": result.attempts,
            }
        )

    return {
        "attempted": len(rows),
        "recovered": len(successes),
        "failed": len(failures),
        "successes": successes,
        "failures": failures,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="CSV with PMID rows")
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for recovered *_FULL_CONTEXT.md files",
    )
    parser.add_argument("--pmid-column", default="pmid")
    parser.add_argument("--doi-column", default=None)
    parser.add_argument("--email", default=None, help="NCBI/Unpaywall email")
    parser.add_argument("--ncbi-api-key", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--write-cleaned",
        action="store_true",
        help="Also write *_CLEANED.md containing only recovered markdown",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Summary JSON path. Default: <output-dir>/oa_recovery_summary.json",
    )
    parser.add_argument("--verbose", "-v", action="store_true")
    return parser


def main() -> int:
    load_dotenv()
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    email = args.email or _env_email()
    if not email:
        raise SystemExit("Set NCBI_EMAIL or ENTREZ_EMAIL, or pass --email.")

    rows = _read_rows(
        args.input,
        pmid_column=args.pmid_column,
        doi_column=args.doi_column,
        limit=args.limit,
    )
    summary = recover_batch(
        rows,
        output_dir=args.output_dir,
        email=email,
        ncbi_api_key=args.ncbi_api_key or os.environ.get("NCBI_API_KEY"),
        write_cleaned=args.write_cleaned,
    )
    summary_path = args.summary or (args.output_dir / "oa_recovery_summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("attempted", "recovered", "failed")}))
    print(f"summary: {summary_path}")
    return 0 if summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
