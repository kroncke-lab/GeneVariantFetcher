#!/usr/bin/env python3
"""Recover Elsevier mmc supplements that the full-text-API route never fetched.

The Elsevier full-text API returns body markdown only; a paper recovered via
that route lands without its ``mmc`` supplement files (where the mutation tables
usually live). This is the single largest supplement-recall bucket
(~31% of the gold-missing-variant gap; see docs/SUPPLEMENT_ACQUISITION_PLAN.md).

For each target paper this script: reads the DOI from the corpus
``{pmid}_artifacts.json``, fetches the authenticated full-text XML
(``harvesting/elsevier_api.ElsevierAPIClient``, API key + insttoken), extracts
the ``1-s2.0-<PII>-mmc<N>.<ext>`` references, and downloads them from the open
ScienceDirect CDN into ``corpus/<GENE>/<PMID>/<PMID>_supplements/``. The existing
re-fold (wired into ``scripts/refresh_run_db.py``) then makes them visible to
Tier-3 on the next refresh.

  # all Elsevier-DOI corpus papers for a gene that lack supplements on disk:
  python scripts/fetch_elsevier_supplements.py --gene SCN5A
  # a specific set:
  python scripts/fetch_elsevier_supplements.py --gene SCN5A --pmids 29325976,20129283
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from harvesting.elsevier_api import ElsevierAPIClient  # noqa: E402
from harvesting.supplement_fold import (  # noqa: E402
    fold_supplements_into_full_context,
)

MANIFEST_NAME = ".gvf_elsevier_supplement_manifest.json"


@dataclass(frozen=True)
class PaperTarget:
    pmid: str
    paper_dir: Path
    doi: str
    reuse_supplement_dirs: tuple[Path, ...] = ()


def _doi_for(pdir: Path, pmid: str) -> str:
    for art in (
        pdir / f"{pmid}_artifacts.json",
        pdir / "result.json",
        pdir / pmid / "result.json",
    ):
        if not art.exists():
            continue
        try:
            doi = (json.loads(art.read_text()).get("doi") or "").strip()
        except (json.JSONDecodeError, OSError):
            continue
        if doi:
            return doi
    full_context = pdir / f"{pmid}_FULL_CONTEXT.md"
    if full_context.exists():
        try:
            text = full_context.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return ""
        labelled = re.search(
            r"doi:\s*(10\.\d{4,9}/[^\s\"'<>]+)", text, flags=re.IGNORECASE
        )
        if labelled:
            return labelled.group(1).strip().strip(".,;:)]}>")
    return ""


def _load_input(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return {
            str(row.get("PMID") or row.get("pmid") or "").strip(): str(
                row.get("DOI") or row.get("doi") or ""
            ).strip()
            for row in reader
            if str(row.get("PMID") or row.get("pmid") or "").strip()
        }


def discover_targets(
    *,
    gene: str,
    corpus: Path,
    harvest_dir: Path | None,
    wanted_pmids: set[str],
    input_dois: dict[str, str],
) -> list[PaperTarget]:
    targets: list[PaperTarget] = []
    if harvest_dir is not None:
        if not harvest_dir.is_dir():
            raise SystemExit(f"No harvest dir: {harvest_dir}")
        pmids = wanted_pmids or {
            path.name.replace("_FULL_CONTEXT.md", "")
            for path in harvest_dir.glob("*_FULL_CONTEXT.md")
        }
        paper_dirs = ((pmid, harvest_dir) for pmid in sorted(pmids))
    else:
        gene_dir = corpus / gene.upper()
        if not gene_dir.is_dir():
            raise SystemExit(f"No corpus dir for gene: {gene_dir}")
        paper_dirs = (
            (pdir.name, pdir)
            for pdir in sorted(gene_dir.iterdir())
            if pdir.is_dir() and (not wanted_pmids or pdir.name in wanted_pmids)
        )

    for pmid, pdir in paper_dirs:
        doi = input_dois.get(pmid) or _doi_for(pdir, pmid)
        if doi and ElsevierAPIClient.is_elsevier_doi(doi):
            reuse_dirs: tuple[Path, ...] = ()
            if harvest_dir is None:
                reuse_dirs = tuple(
                    sibling / pmid / f"{pmid}_supplements"
                    for sibling in sorted(corpus.iterdir())
                    if sibling.is_dir() and sibling.name.upper() != gene.upper()
                )
            targets.append(
                PaperTarget(
                    pmid=pmid,
                    paper_dir=pdir,
                    doi=doi,
                    reuse_supplement_dirs=reuse_dirs,
                )
            )
    return targets


def _is_complete_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 1000


def _cached_complete_refs(
    client: ElsevierAPIClient, supp_dir: Path, doi: str
) -> list[str]:
    manifest = supp_dir / MANIFEST_NAME
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    refs = data.get("refs") or []
    if data.get("doi") != doi or not isinstance(refs, list) or not refs:
        return []
    clean_refs = [str(ref) for ref in refs if str(ref).strip()]
    if len(clean_refs) != len(refs):
        return []
    expected = [supp_dir / client.supplement_local_name(ref) for ref in clean_refs]
    return clean_refs if all(_is_complete_file(path) for path in expected) else []


def _write_complete_manifest(supp_dir: Path, doi: str, refs: list[str]) -> None:
    supp_dir.mkdir(parents=True, exist_ok=True)
    (supp_dir / MANIFEST_NAME).write_text(
        json.dumps({"doi": doi, "refs": refs}, indent=2) + "\n",
        encoding="utf-8",
    )


def augment_paper(client: ElsevierAPIClient, target: PaperTarget) -> dict:
    """Fetch only missing mmc files, then idempotently fold the complete set."""
    supp_dir = target.paper_dir / f"{target.pmid}_supplements"
    cached_refs = _cached_complete_refs(client, supp_dir, target.doi)
    if cached_refs:
        folded = (
            fold_supplements_into_full_context(
                target.pmid,
                target.paper_dir,
                supplements_dir=supp_dir,
            )
            is not None
        )
        return {
            "pmid": target.pmid,
            "doi": target.doi,
            "refs": len(cached_refs),
            "new_files": [],
            "complete": True,
            "folded": folded,
            "error": "",
        }

    xml, err = client.get_fulltext_by_doi(target.doi)
    if not xml:
        folded = (
            fold_supplements_into_full_context(
                target.pmid,
                target.paper_dir,
                supplements_dir=supp_dir,
            )
            is not None
        )
        return {
            "pmid": target.pmid,
            "doi": target.doi,
            "refs": 0,
            "new_files": [],
            "complete": False,
            "folded": folded,
            "error": err or "no XML",
        }

    refs = client.extract_supplement_refs(xml)
    expected = {ref: supp_dir / client.supplement_local_name(ref) for ref in refs}
    before = {
        path.name: path.stat().st_size for path in expected.values() if path.exists()
    }
    for path in expected.values():
        if _is_complete_file(path):
            continue
        for reuse_dir in target.reuse_supplement_dirs:
            reusable = reuse_dir / path.name
            if _is_complete_file(reusable):
                supp_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(reusable, path)
                break
    missing = [ref for ref, path in expected.items() if not _is_complete_file(path)]
    if missing:
        client.download_supplements(xml, supp_dir)

    after = {
        path.name: path.stat().st_size for path in expected.values() if path.exists()
    }
    new_files = sorted(
        name for name, size in after.items() if size > 1000 and before.get(name) != size
    )
    complete = bool(refs) and all(_is_complete_file(path) for path in expected.values())
    if complete:
        _write_complete_manifest(supp_dir, target.doi, refs)
    folded = (
        fold_supplements_into_full_context(
            target.pmid,
            target.paper_dir,
            supplements_dir=supp_dir,
        )
        is not None
    )
    return {
        "pmid": target.pmid,
        "doi": target.doi,
        "refs": len(refs),
        "new_files": new_files,
        "complete": complete,
        "folded": folded,
        "error": "" if refs else "no mmc references in XML",
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument(
        "--harvest-dir",
        type=Path,
        default=None,
        help="Flat run pmc_fulltext directory instead of nested corpus layout.",
    )
    ap.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Optional CSV with PMID and DOI columns (used by gvf-run supplement queue).",
    )
    ap.add_argument(
        "--source-override-out",
        type=Path,
        default=None,
        help="Optional refresh_run_db source-override CSV for sources changed by folding.",
    )
    ap.add_argument("--pmids", default="", help="Comma-separated PMIDs (default: all).")
    ap.add_argument(
        "--include-with-supplements",
        action="store_true",
        help=(
            "Compatibility flag; partial sets are now always rechecked and complete "
            "sets are skipped per file."
        ),
    )
    ap.add_argument("--limit", type=int, default=0, help="Max papers to process.")
    args = ap.parse_args()

    client = ElsevierAPIClient(
        api_key=os.getenv("ELSEVIER_API_KEY"),
        insttoken=os.getenv("ELSEVIER_INSTTOKEN"),
    )
    if not client.is_available:
        raise SystemExit("ELSEVIER_API_KEY not set; cannot fetch supplements.")

    input_dois = _load_input(args.input)
    want = {p.strip() for p in args.pmids.split(",") if p.strip()} | set(input_dois)
    targets = discover_targets(
        gene=args.gene,
        corpus=Path(args.corpus),
        harvest_dir=args.harvest_dir,
        wanted_pmids=want,
        input_dois=input_dois,
    )

    if args.limit:
        targets = targets[: args.limit]
    print(
        f"{args.gene}: {len(targets)} Elsevier-DOI paper(s) to reconcile "
        "(existing files are never re-downloaded)"
    )

    papers_with_supp = 0
    total_files = 0
    papers_complete = 0
    folded_papers = 0
    source_overrides: list[dict[str, object]] = []
    for target in targets:
        result = augment_paper(client, target)
        new_files = result["new_files"]
        if new_files:
            papers_with_supp += 1
            total_files += len(new_files)
            print(
                f"  {target.pmid} ({target.doi}): +{len(new_files)} "
                f"supplement(s) {new_files}"
            )
        else:
            suffix = f" ({result['error']})" if result["error"] else ""
            print(f"  {target.pmid} ({target.doi}): no new files{suffix}")
        papers_complete += int(result["complete"])
        folded_papers += int(result["folded"])
        full_context = target.paper_dir / f"{target.pmid}_FULL_CONTEXT.md"
        if result["folded"] and full_context.exists():
            source_overrides.append(
                {
                    "gene": args.gene.upper(),
                    "pmid": target.pmid,
                    "action": "refresh_replay",
                    "route": "fetch_elsevier_supplements_only",
                    "available_context_path": str(full_context),
                    "available_context_bytes": full_context.stat().st_size,
                    "notes": "supplement set changed and was folded into full text",
                }
            )

    if args.source_override_out is not None:
        args.source_override_out.parent.mkdir(parents=True, exist_ok=True)
        with args.source_override_out.open("w", newline="", encoding="utf-8") as handle:
            fieldnames = [
                "gene",
                "pmid",
                "action",
                "route",
                "available_context_path",
                "available_context_bytes",
                "notes",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(source_overrides)

    print(
        f"\nDONE: recovered {total_files} supplement file(s) across "
        f"{papers_with_supp}/{len(targets)} papers; {papers_complete} complete mmc "
        f"sets; {folded_papers} FULL_CONTEXT file(s) updated. Re-extract to realize "
        "new variants."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
