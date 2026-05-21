"""Unlock paywalled Elsevier articles via ElsevierAPIClient + insttoken.

Writes each unlocked full-text body to ``{output_dir}/{PMID}_FULL_CONTEXT.md``
so the standard extraction discovery path (``cli.extract.find_input_files``
with ``--full-text``) finds it without per-PMID-dir plumbing. Existing stub
files at that path are preserved as ``.pre_insttoken_bak`` before overwrite.

Usage::

    .venv/bin/python scripts/test_insttoken_unlock.py \\
        --input  validation_runs/.../<GENE>/<TS>/pmc_fulltext/paywalled_missing.csv \\
        --output validation_runs/.../<GENE>/<TS>/pmc_fulltext

The same directory holds the rest of the run's full-text artifacts, so the
new files live alongside them rather than in a separate subdirectory.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path

# Allow running as a standalone script from the repo root
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# Match a standard Elsevier DOI. We deliberately DO NOT include `/` in the
# trailing character class -- DOIs do not contain `/` after the registrant
# prefix in practice, and including `/` causes the regex to greedily swallow
# URL path segments like ``.../doi/10.1016/j.cjca.2012.12.002/attachment/...``
# and mis-treat the attachment path as part of the DOI.
ELSEVIER_DOI_RE = re.compile(r"10\.10(16|53|67|54|06)/[A-Za-z0-9._\-()<>:;]+")


def load_env(env_path: Path) -> None:
    for line in env_path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or "=" not in s:
            continue
        k, _, v = s.partition("=")
        os.environ[k.strip()] = v.strip()


def extract_elsevier_doi(row: dict[str, str]) -> str | None:
    url = (row.get("URL") or "").strip()
    if not url:
        return None
    m = ELSEVIER_DOI_RE.search(url)
    if not m:
        return None
    return m.group(0).rstrip(").,;")


def backup_existing(path: Path) -> Path | None:
    """If ``path`` exists, rename it to ``<path>.pre_insttoken_bak``.

    Returns the backup path, or ``None`` if the original did not exist.
    Existing backup files are not overwritten -- the first backup is kept.
    """
    if not path.exists():
        return None
    bak = path.with_suffix(path.suffix + ".pre_insttoken_bak")
    if not bak.exists():
        shutil.copy2(path, bak)
    return bak


def write_full_context(
    pmid: str, markdown: str, output_dir: Path
) -> tuple[Path, Path | None]:
    """Write the unlocked markdown body to ``{PMID}_FULL_CONTEXT.md``.

    Returns ``(written_path, backup_path_or_None)``.
    """
    target = output_dir / f"{pmid}_FULL_CONTEXT.md"
    backup = backup_existing(target)
    target.write_text(markdown, encoding="utf-8")
    return target, backup


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help=(
            "Directory to write {PMID}_FULL_CONTEXT.md files into. Should be"
            " the run's existing pmc_fulltext/ so the new bodies sit alongside"
            " the rest of the run's full-text artifacts."
        ),
    )
    parser.add_argument(
        "--env",
        type=Path,
        default=_REPO_ROOT / ".env",
        help="Path to .env (default: <repo>/.env)",
    )
    parser.add_argument(
        "--summary-name",
        default="insttoken_unlock_results.csv",
        help="Filename for the unlock-summary CSV written into --output.",
    )
    args = parser.parse_args()

    if args.env.exists():
        load_env(args.env)

    api_key = os.environ.get("ELSEVIER_API_KEY")
    insttoken = os.environ.get("ELSEVIER_INSTTOKEN")
    if not api_key:
        print("ERROR: ELSEVIER_API_KEY not set", file=sys.stderr)
        return 2
    if not insttoken:
        print("ERROR: ELSEVIER_INSTTOKEN not set", file=sys.stderr)
        return 2

    from harvesting.elsevier_api import ElsevierAPIClient

    client = ElsevierAPIClient(api_key=api_key, insttoken=insttoken)

    args.output.mkdir(parents=True, exist_ok=True)

    # Deduplicate by (PMID, DOI)
    seen: set[tuple[str, str]] = set()
    candidates: list[tuple[str, str]] = []
    with args.input.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pmid = (row.get("PMID") or "").strip()
            doi = extract_elsevier_doi(row)
            if not pmid or not doi:
                continue
            key = (pmid, doi)
            if key in seen:
                continue
            seen.add(key)
            candidates.append(key)

    print(f"Found {len(candidates)} unique (PMID, Elsevier-DOI) pairs to probe.")
    print(f"Writing into: {args.output}")

    results_path = args.output / args.summary_name
    success = 0
    metadata_only = 0
    error_count = 0
    backed_up = 0

    with results_path.open("w", newline="") as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(
            [
                "PMID",
                "DOI",
                "outcome",
                "xml_chars",
                "markdown_chars",
                "has_body",
                "full_context_path",
                "backup_path",
                "error",
            ]
        )

        for pmid, doi in candidates:
            xml, xml_err = client.get_fulltext_by_doi(doi)
            xml_chars = len(xml) if xml else 0
            has_body = False
            if xml:
                has_body = any(
                    tag in xml
                    for tag in (
                        "<originalText",
                        "<ce:originalText",
                        "<body",
                        "<ce:body",
                        "<xocs:rawtext",
                    )
                )

            md, md_err = client.fetch_fulltext(doi=doi) if xml else (None, xml_err)
            md_chars = len(md) if md else 0

            outcome: str
            full_context_path = ""
            backup_path = ""

            if md and md_chars > 2000:
                outcome = "FULLTEXT"
                success += 1
                target, backup = write_full_context(pmid, md, args.output)
                full_context_path = str(target)
                if backup is not None:
                    backup_path = str(backup)
                    backed_up += 1
            elif xml and xml_chars > 50_000 and has_body:
                outcome = "FULLTEXT_XML_ONLY"
                success += 1
                # Write XML as a sidecar; markdown conversion underperformed.
                target = args.output / f"{pmid}_FULL_CONTEXT.md"
                backup = backup_existing(target)
                # Persist as a markdown wrapper around the original XML body
                # so the file is discoverable by the existing pipeline.
                target.write_text(
                    f"# MAIN TEXT\n\n<!-- xml_only fallback; "
                    f"xml_chars={xml_chars} -->\n\n{xml}",
                    encoding="utf-8",
                )
                full_context_path = str(target)
                if backup is not None:
                    backup_path = str(backup)
                    backed_up += 1
            elif xml:
                outcome = "METADATA_ONLY"
                metadata_only += 1
            else:
                outcome = "ERROR"
                error_count += 1

            err_str = md_err or xml_err or ""
            writer.writerow(
                [
                    pmid,
                    doi,
                    outcome,
                    xml_chars,
                    md_chars,
                    str(has_body),
                    full_context_path,
                    backup_path,
                    err_str,
                ]
            )
            print(
                f"  {pmid:>10}  {doi:<48}  {outcome:<20}  "
                f"xml={xml_chars:>7,}  md={md_chars:>6,}"
                + (f"  [bak]" if backup_path else "")
            )

    print()
    print(f"Summary CSV: {results_path}")
    print(
        f"  FULLTEXT(+XML_ONLY): {success}/{len(candidates)}  "
        f"METADATA_ONLY: {metadata_only}/{len(candidates)}  "
        f"ERROR: {error_count}/{len(candidates)}"
    )
    print(f"  pre-existing stub files backed up to .pre_insttoken_bak: {backed_up}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
