"""
Recover full text for PMIDs currently stuck as ABSTRACT-ONLY FALLBACK.

Reads a directory of `*_FULL_CONTEXT.md` files, identifies the ones flagged
as abstract-only, and tries the OA recovery chain (PMC HTML → Europe PMC XML
→ Europe PMC PDF render → Unpaywall) for each. On success the abstract-only
FULL_CONTEXT.md is overwritten with the recovered body and a `_OA_RECOVERED`
sidecar file records which source supplied the body.

Usage:
    python -m scripts.recover_abstract_only \
        --harvest-dir results/KCNH2/20260506_102238/pmc_fulltext \
        --out-summary results/KCNH2/20260506_102238/oa_recovery_summary.json

    # Dry run — classify only, do not overwrite:
    python -m scripts.recover_abstract_only --harvest-dir ... --dry-run

    # Subset by PMID list:
    python -m scripts.recover_abstract_only --harvest-dir ... \
        --pmid-file results/KCNH2/20260506_102238/gold_missing_pmids.txt
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

# Allow `python scripts/recover_abstract_only.py` from repo root
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from harvesting.oa_recovery import OARecoveryClient, is_acceptable_body  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-7s %(name)s :: %(message)s",
)
logger = logging.getLogger("oa_recovery")


ABSTRACT_ONLY_MARKER = "abstract-only fallback"


def load_env(env_path: Path) -> None:
    if not env_path.exists():
        return
    for line in env_path.read_text().splitlines():
        if line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))


def find_abstract_only(harvest_dir: Path) -> list[str]:
    """Return PMIDs whose FULL_CONTEXT.md starts with the abstract-only marker."""
    out = []
    for f in sorted(harvest_dir.glob("*_FULL_CONTEXT.md")):
        pmid = f.name.split("_")[0]
        if not pmid.isdigit():
            continue
        head = f.read_text(errors="ignore")[:500].lower()
        if ABSTRACT_ONLY_MARKER in head:
            out.append(pmid)
    return out


def extract_doi(md_text: str) -> str | None:
    """Pull a DOI out of the abstract-only markdown if one is present."""
    import re

    # Prefer "doi: <X>" or "doi.org/<X>" anchors when present
    anchored = re.search(
        r"(?:doi\.org/|doi:\s*)(10\.\d{4,9}/\S+)",
        md_text,
        flags=re.IGNORECASE,
    )
    candidate = anchored.group(1) if anchored else None
    if not candidate:
        # Fall back to bare-DOI match; \S+ keeps parens so Wiley legacy DOIs
        # (10.1002/(SICI)...) are not truncated at the first ')'.
        m = re.search(r"\b(10\.\d{4,9}/\S+)", md_text)
        if not m:
            return None
        candidate = m.group(1)

    # Trim trailing markdown/punctuation noise.
    candidate = candidate.rstrip(".,;:>])'\"")
    # Some abstract dumps include the DOI followed by URL fragments
    # (e.g. ".../asset/...", ".../main.assets/..."). Strip those.
    for marker in (
        "/asset/",
        "/main.assets/",
        "/full",
        "/abstract",
        "/pdf",
        "/figures/",
        "/figure/",
    ):
        idx = candidate.lower().find(marker)
        if idx > 0:
            candidate = candidate[:idx]
            break
    return candidate.lower()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harvest-dir", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, default=None)
    parser.add_argument(
        "--pmid-file",
        type=Path,
        default=None,
        help="Optional file containing PMIDs (one per line) to limit recovery to.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Try recovery but do not overwrite any FULL_CONTEXT.md",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Stop after this many PMIDs (0 = no limit).",
    )
    parser.add_argument("--env", type=Path, default=ROOT / ".env")
    args = parser.parse_args()

    load_env(args.env)
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        print("ERROR: NCBI_EMAIL not set in environment or .env", file=sys.stderr)
        return 2

    if not args.harvest_dir.is_dir():
        print(f"ERROR: not a directory: {args.harvest_dir}", file=sys.stderr)
        return 2

    if args.pmid_file:
        abstract_pmids = [
            p.strip()
            for p in args.pmid_file.read_text().splitlines()
            if p.strip() and not p.startswith("#") and p.strip().isdigit()
        ]
        logger.info(
            "Loaded %d PMIDs from %s (skipping abstract-only marker filter)",
            len(abstract_pmids),
            args.pmid_file,
        )
    else:
        abstract_pmids = find_abstract_only(args.harvest_dir)
        logger.info(
            "Found %d abstract-only PMIDs in %s", len(abstract_pmids), args.harvest_dir
        )

    if args.limit:
        abstract_pmids = abstract_pmids[: args.limit]

    client = OARecoveryClient(
        email=email,
        ncbi_api_key=os.environ.get("NCBI_API_KEY"),
    )

    summary = {"by_pmid": {}, "by_source": {}, "by_status": {"success": 0, "fail": 0}}
    t0 = time.time()

    for i, pmid in enumerate(abstract_pmids, 1):
        full_context = args.harvest_dir / f"{pmid}_FULL_CONTEXT.md"
        original_text = (
            full_context.read_text(errors="ignore") if full_context.exists() else ""
        )
        doi = extract_doi(original_text)

        logger.info(
            "[%d/%d] PMID %s  (doi=%s)", i, len(abstract_pmids), pmid, doi or "?"
        )
        try:
            result = client.recover(pmid, doi=doi)
        except Exception as e:
            logger.exception("recover crashed for %s: %s", pmid, e)
            summary["by_status"]["fail"] += 1
            summary["by_pmid"][pmid] = {"success": False, "error": str(e)}
            continue

        rec = {
            "success": result.success,
            "source": result.source,
            "pmcid": result.pmcid,
            "n_chars": result.n_chars,
            "attempts": [
                {"source": s, "status": st, "detail": str(d)[:200]}
                for s, st, d in result.attempts
            ],
            "error": result.error,
        }
        summary["by_pmid"][pmid] = rec

        if result.success:
            summary["by_status"]["success"] += 1
            summary["by_source"][result.source] = (
                summary["by_source"].get(result.source, 0) + 1
            )
            logger.info("  RECOVERED via %s (%d chars)", result.source, result.n_chars)
            if not args.dry_run:
                # Sanity gate: never overwrite if recovered body looks worse than what we have
                ok, reason = is_acceptable_body(result.markdown)
                if not ok:
                    logger.warning(
                        "  Refusing overwrite — recovered body failed gate: %s", reason
                    )
                    summary["by_pmid"][pmid]["overwrite"] = f"refused: {reason}"
                    continue
                # Backup original, write recovered, drop a sidecar manifest
                backup = full_context.with_suffix(".abstract.bak.md")
                if not backup.exists() and full_context.exists():
                    backup.write_text(original_text)
                header = "<!-- oa_recovery: source={src}, pmcid={pmc}, chars={n} -->\n".format(
                    src=result.source, pmc=result.pmcid or "-", n=result.n_chars
                )
                full_context.write_text(header + result.markdown)
                marker = args.harvest_dir / f"{pmid}_OA_RECOVERED.json"
                marker.write_text(json.dumps(rec, indent=2))
                summary["by_pmid"][pmid]["overwrite"] = "ok"
        else:
            summary["by_status"]["fail"] += 1
            logger.info(
                "  no source produced acceptable body (%d attempts)",
                len(result.attempts),
            )

    dur = time.time() - t0
    logger.info(
        "Done in %.1fs.  Success %d / %d (%.1f%%).  By source: %s",
        dur,
        summary["by_status"]["success"],
        len(abstract_pmids),
        100 * summary["by_status"]["success"] / max(1, len(abstract_pmids)),
        summary["by_source"],
    )

    if args.out_summary:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        args.out_summary.write_text(json.dumps(summary, indent=2))
        logger.info("Wrote summary → %s", args.out_summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
