#!/usr/bin/env python3
"""Fetch the supplements a paper advertised but that we never downloaded.

``figure_extractor`` records every supplement link it finds in a paper's markup,
but those links were write-only for years: rendered into ``FULL_CONTEXT.md`` as
``_link_: <href>`` lines and read by nothing. The real download path discovers
supplements through the PMC and publisher APIs instead, so on any paper those
APIs missed, links we already had went unused.

Measured on the corpus 2026-07-24: **470 papers carry file-shaped supplement
links and have an empty supplements directory**, 430 of them with a PMCID on
record. PMC serves every supplement out of the article's ``bin/`` directory, so
a bare filename plus a PMCID is fully resolvable - 1,467 of 1,550 links, free
and with no authentication.

This script closes the loop using existing machinery only:
``supplement_link_resolver`` turns recorded links into download jobs,
``PMCHarvester.download_supplement`` fetches them (with its PMC URL-variant
retries and HTML-error-page validation), ``process_supplement_files`` converts
them, and ``fold_supplements_into_full_context`` makes them visible to
extraction. Nothing here reimplements any of that.

  # see what is recoverable, no network:
  python scripts/fetch_linked_supplements.py --dry-run
  # fetch one gene:
  python scripts/fetch_linked_supplements.py --gene KCNQ1
  # everything, with a source-override CSV for refresh_run_db:
  python scripts/fetch_linked_supplements.py --source-override-out /tmp/overrides.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from urllib.parse import unquote, urlparse

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from harvesting.supplement_fold import (  # noqa: E402
    fold_supplements_into_full_context,
)
from harvesting.supplement_link_resolver import (  # noqa: E402
    ZIP_MAGIC,
    ResolutionReport,
    SupplementLink,
    europepmc_archive_job,
    links_from_artifacts,
    looks_like_supplement_file,
    parse_links_from_markdown,
    pmcid_from_artifacts,
    resolve_links,
)

# Europe PMC asks for restraint rather than a hard rate; one archive request per
# paper per ~0.4s is well inside anything either service objects to.
NCBI_SLEEP_SECONDS = 0.4

# When the per-file route is refused it is refused for the whole article, not
# one file: NCBI gates the article's bin/ directory. Each failed file costs
# ~18 requests (3 attempts x 6 URL variants) plus backoff, so a paper with ten
# gated links burns minutes to learn what the first two already proved. Stop
# after this many consecutive failures on a paper and record the rest as
# skipped. Only applies to the per-file fallback; the archive route is one
# request either way.
PER_FILE_FAILURE_BUDGET = 2


@dataclass
class PaperTarget:
    """A corpus paper with recorded supplement links we have not fetched."""

    gene: str
    pmid: str
    paper_dir: Path
    pmcid: Optional[str] = None
    doi: Optional[str] = None
    links: list[SupplementLink] = field(default_factory=list)
    on_disk: set[str] = field(default_factory=set)

    @property
    def supplements_dir(self) -> Path:
        return self.paper_dir / f"{self.pmid}_supplements"


def _read_artifacts(paper_dir: Path, pmid: str) -> dict:
    for candidate in (
        paper_dir / f"{pmid}_artifacts.json",
        paper_dir / "artifacts.json",
    ):
        if not candidate.is_file():
            continue
        try:
            data = json.loads(candidate.read_text(encoding="utf-8", errors="replace"))
        except (json.JSONDecodeError, OSError):
            continue
        if isinstance(data, dict):
            return data
    return {}


def _href_filename(href: str) -> str:
    """The filename an href names, ignoring any query string, lowercased."""
    return Path(unquote(urlparse(href).path or href)).name.lower()


def _files_on_disk(paper_dir: Path, pmid: str) -> set[str]:
    """Filenames already present in this paper's supplement directories.

    Scoped by PMID rather than a bare ``*_supplements`` glob: in a flat run
    harvest directory every paper shares one parent, so an unscoped glob would
    report another paper's files as this one's.
    """
    names: set[str] = set()
    for supp_dir in paper_dir.glob(f"{pmid}_supplements"):
        for path in supp_dir.rglob("*"):
            if path.is_file() and not path.name.startswith("."):
                names.add(path.name.lower())
    return names


def _collect_links(paper_dir: Path, pmid: str, artifacts: dict) -> list[SupplementLink]:
    """Structured links when the harvest recorded them, else parse the markdown.

    Newer harvests write ``supplement_links`` into the artifacts log. Everything
    harvested before that only has the rendered ``_link_:`` lines, which is the
    whole corpus backlog - so the markdown path is not a fallback, it is the
    main road for now.
    """
    links = links_from_artifacts(artifacts)
    if links:
        return links

    seen: set[str] = set()
    out: list[SupplementLink] = []
    for md in sorted(paper_dir.glob(f"{pmid}*.md")):
        try:
            text = md.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for link in parse_links_from_markdown(text, source=f"markdown:{md.name}"):
            if link.href in seen:
                continue
            seen.add(link.href)
            out.append(link)
    return out


def _paper_dirs(
    corpus: Path,
    harvest_dir: Optional[Path],
    genes: Optional[set[str]],
    wanted_pmids: Optional[set[str]],
):
    """Yield ``(gene, pmid, paper_dir)`` for both corpus and flat run layouts.

    ``corpus/<GENE>/<PMID>/`` gives one directory per paper; a run's
    ``pmc_fulltext/`` is flat, with every paper's files sharing one directory
    and distinguished only by their ``<PMID>_`` prefix.
    """
    if harvest_dir is not None:
        pmids = wanted_pmids or {
            path.name.replace("_FULL_CONTEXT.md", "")
            for path in harvest_dir.glob("*_FULL_CONTEXT.md")
        }
        gene = next(iter(genes)) if genes and len(genes) == 1 else ""
        for pmid in sorted(pmids):
            yield gene, pmid, harvest_dir
        return

    if not corpus.is_dir():
        return
    for gene_dir in sorted(corpus.iterdir()):
        if not gene_dir.is_dir() or gene_dir.name.startswith("."):
            continue
        if genes and gene_dir.name.upper() not in genes:
            continue
        for paper_dir in sorted(gene_dir.iterdir()):
            if not paper_dir.is_dir() or paper_dir.name.startswith("."):
                continue
            if wanted_pmids and paper_dir.name not in wanted_pmids:
                continue
            yield gene_dir.name, paper_dir.name, paper_dir


def discover_targets(
    corpus: Path,
    *,
    harvest_dir: Optional[Path] = None,
    genes: Optional[set[str]] = None,
    wanted_pmids: Optional[set[str]] = None,
    include_partial: bool = False,
) -> list[PaperTarget]:
    """Find papers with recorded links whose files are not on disk.

    By default only papers whose supplement directories hold *none* of the
    linked files are considered. ``include_partial`` widens that to any paper
    missing at least one linked file.
    """
    targets: list[PaperTarget] = []
    for gene, pmid, paper_dir in _paper_dirs(corpus, harvest_dir, genes, wanted_pmids):
        artifacts = _read_artifacts(paper_dir, pmid)
        links = [
            link
            for link in _collect_links(paper_dir, pmid, artifacts)
            if looks_like_supplement_file(link.href)
        ]
        if not links:
            continue

        on_disk = _files_on_disk(paper_dir, pmid)
        missing = [link for link in links if _href_filename(link.href) not in on_disk]
        if not missing:
            continue
        if on_disk and not include_partial:
            # Something was fetched for this paper already; leave it alone
            # unless the caller explicitly wants partial sets topped up.
            continue

        doi = artifacts.get("doi")
        targets.append(
            PaperTarget(
                gene=gene,
                pmid=pmid,
                paper_dir=paper_dir,
                pmcid=pmcid_from_artifacts(artifacts),
                doi=doi.strip() if isinstance(doi, str) else None,
                links=missing,
                on_disk=on_disk,
            )
        )
    return targets


def plan(target: PaperTarget) -> ResolutionReport:
    """Resolve a paper's recorded links, one job per link.

    Deliberately free of route selection: ``fetch_paper`` decides whether to
    use the Europe PMC archive or these per-file jobs. Keeping resolution pure
    is what lets the archive genuinely fall back to per-file when Europe PMC
    will not serve the article.
    """
    return resolve_links(target.links, pmid=target.pmid, pmcid=target.pmcid)


def _is_zip(path: Path) -> bool:
    try:
        with path.open("rb") as fh:
            return fh.read(4) == ZIP_MAGIC
    except OSError:
        return False


def _download_callback(harvester, archive_names: set[str], *, failure_budget=None):
    """Wrap download_supplement with an archive check and a per-paper breaker.

    Two things the generic downloader cannot know:

    * Europe PMC answers a non-open-access article with an XML ``errorBean``
      under a 200, so a "successful" download can be 164 bytes of error text
      named ``.zip``. Verify the magic bytes or the converter records a bogus
      supplement.
    * A refused per-file route is refused for the whole article. Once
      *failure_budget* files in a row have failed, stop making requests for
      this paper instead of paying ~18 requests per remaining file.
    """
    state = {"consecutive_failures": 0}

    def _cb(url, file_path, pmid, filename, supp):
        if (
            failure_budget is not None
            and state["consecutive_failures"] >= failure_budget
        ):
            print(f"    (skipping {filename}: article's per-file route is refused)")
            return False

        ok = harvester.download_supplement(
            url,
            file_path,
            pmid,
            filename,
            supp.get("base_url"),
            supp.get("original_url"),
        )
        state["consecutive_failures"] = 0 if ok else state["consecutive_failures"] + 1
        if not ok or filename not in archive_names:
            return ok
        try:
            with file_path.open("rb") as fh:
                head = fh.read(4)
        except OSError:
            return False
        if head == ZIP_MAGIC:
            return True
        file_path.unlink(missing_ok=True)
        print(f"    (not an archive: {filename} — article is not open access)")
        return False

    return _cb


def _run_jobs(harvester, target: PaperTarget, jobs: list[dict], archive: bool):
    """Download and convert one batch of jobs into the paper's supplements dir."""
    from harvesting.supplement_processing_service import process_supplement_files

    archive_names = {job["name"] for job in jobs} if archive else set()
    return process_supplement_files(
        supp_files=jobs,
        supplements_dir=target.supplements_dir,
        pmid=target.pmid,
        converter=harvester.converter,
        download_callback=_download_callback(
            harvester,
            archive_names,
            failure_budget=None if archive else PER_FILE_FAILURE_BUDGET,
        ),
        extract_figures=False,
        figures_dir=None,
        sleep_seconds=NCBI_SLEEP_SECONDS,
    )


def fetch_paper(harvester, target: PaperTarget, report: ResolutionReport) -> dict:
    """Download, convert, and fold one paper's supplements.

    Two routes, in order of what actually works. For a paper with a PMCID the
    Europe PMC archive endpoint returns every supplement in one ZIP; NCBI's
    per-file ``bin/`` URLs 403 an unattended client on every variant the
    scraper generates (verified live 2026-07-24). The per-file jobs are the
    fallback for non-PMC papers and for articles Europe PMC will not serve.
    """
    row = {
        "gene": target.gene,
        "pmid": target.pmid,
        "pmcid": target.pmcid or "",
        "planned": len(report.fetchable),
        "downloaded": 0,
        "unresolved": len(report.unresolved),
        "route": "",
        "folded": False,
        "new_files": [],
        "files_on_disk": 0,
        # Recorded here rather than rebuilt from gene/pmid: in a flat run
        # harvest dir there is no <gene>/<pmid>/ path to rebuild from.
        "full_context": str(target.paper_dir / f"{target.pmid}_FULL_CONTEXT.md"),
        "error": "",
    }

    results = []
    if target.pmcid:
        archive = _run_jobs(
            harvester, target, [europepmc_archive_job(target.pmcid).to_dict()], True
        )
        if archive.downloaded_count:
            row["route"] = "europepmc_archive"
            results.append(archive)

    if not results:
        # Europe PMC declined (closed access, or no archive for this article).
        # Fall back to the per-link URLs — 127 papers in the 2026-07-24 sweep
        # landed this way, so the per-file route is not universally dead.
        jobs = [job.to_dict() for job in report.fetchable]
        if not jobs:
            row["error"] = (
                "europepmc archive unavailable and no per-file link resolved"
                if target.pmcid
                else "no resolvable links"
            )
            return row
        per_file = _run_jobs(harvester, target, jobs, False)
        if per_file.downloaded_count:
            row["route"] = "per_file"
        results.append(per_file)

    downloaded = sum(r.downloaded_count for r in results)
    row["downloaded"] = downloaded
    # An archive downloads as one file but yields many; report what landed on
    # disk, not the number of HTTP fetches.
    new_files: list[str] = []
    for r in results:
        for fr in r.file_results:
            if not fr.downloaded:
                continue
            nested = [Path(n).name for n in fr.nested_files]
            new_files.extend(nested or [fr.filename])
    row["new_files"] = new_files
    row["files_on_disk"] = len(new_files)
    errors = [str(fr.error) for r in results for fr in r.file_results if fr.error]
    row["error"] = "; ".join(errors[:2])

    if downloaded:
        row["folded"] = (
            fold_supplements_into_full_context(
                target.pmid,
                target.paper_dir,
                supplements_dir=target.supplements_dir,
                converter=harvester.converter,
            )
            is not None
        )
    return row


def _print_plan(targets: list[PaperTarget]) -> tuple[int, int]:
    by_gene: dict[str, list[int]] = {}
    total_fetchable = 0
    total_unresolved = 0
    for target in targets:
        report = plan(target)
        total_fetchable += len(report.fetchable)
        total_unresolved += len(report.unresolved)
        row = by_gene.setdefault(target.gene, [0, 0, 0])
        row[0] += 1
        row[1] += len(report.fetchable)
        row[2] += len(report.unresolved)

    print(f"{'gene':10s} {'papers':>7s} {'fetchable':>10s} {'unresolved':>11s}")
    for gene, (papers, fetchable, unresolved) in sorted(by_gene.items()):
        print(f"{gene:10s} {papers:7d} {fetchable:10d} {unresolved:11d}")
    print(
        f"{'TOTAL':10s} {len(targets):7d} {total_fetchable:10d} {total_unresolved:11d}"
    )
    return total_fetchable, total_unresolved


def _write_report(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "gene",
                "pmid",
                "pmcid",
                "planned",
                "downloaded",
                "unresolved",
                "route",
                "files_on_disk",
                "folded",
                "full_context",
                "new_files",
                "error",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "new_files": ";".join(row["new_files"])})


def _write_source_overrides(path: Path, rows: list[dict]) -> None:
    """A refresh_run_db override row per paper whose FULL_CONTEXT actually grew."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pmid", "source_path"])
        for row in rows:
            if row["folded"] and row["full_context"]:
                writer.writerow([row["pmid"], row["full_context"]])


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--corpus", type=Path, default=REPO / "corpus")
    ap.add_argument(
        "--harvest-dir",
        type=Path,
        default=None,
        help=(
            "Flat run pmc_fulltext directory instead of the nested corpus "
            "layout (used by gvf-run's source-recovery step)."
        ),
    )
    ap.add_argument(
        "--gene",
        action="append",
        default=[],
        help="Restrict to a gene (repeatable). Default: every gene in the corpus.",
    )
    ap.add_argument("--pmids", default="", help="Comma-separated PMIDs.")
    ap.add_argument(
        "--include-partial",
        action="store_true",
        help="Also top up papers that already have some supplements on disk.",
    )
    ap.add_argument("--limit", type=int, default=0, help="Max papers to process.")
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what is resolvable and exit without fetching anything.",
    )
    ap.add_argument(
        "--report-out",
        type=Path,
        default=None,
        help="Write a per-paper CSV of outcomes.",
    )
    ap.add_argument(
        "--source-override-out",
        type=Path,
        default=None,
        help="Write a refresh_run_db source-override CSV for refolded papers.",
    )
    ap.add_argument(
        "--scratch-dir",
        type=Path,
        default=None,
        help=(
            "Working directory for the harvester's own CSV logs (never the "
            "corpus). Defaults beside the run when --harvest-dir is given, "
            "otherwise to .gvf_linked_supplements/ in the repo."
        ),
    )
    args = ap.parse_args()

    # Keep the harvester's own bookkeeping with the run it belongs to rather
    # than dropping an untracked directory in the repo root.
    scratch_dir = args.scratch_dir or (
        args.harvest_dir.parent / "linked_supplement_fetch"
        if args.harvest_dir is not None
        else REPO / ".gvf_linked_supplements"
    )

    genes = {g.strip().upper() for g in args.gene if g.strip()} or None
    wanted = {p.strip() for p in args.pmids.split(",") if p.strip()} or None

    if args.harvest_dir is not None and not args.harvest_dir.is_dir():
        # A run that skipped extraction has no harvest dir; that is not an
        # error, there is simply nothing to recover.
        print(f"No harvest dir at {args.harvest_dir}; nothing to recover.")
        return 0

    targets = discover_targets(
        args.corpus,
        harvest_dir=args.harvest_dir,
        genes=genes,
        wanted_pmids=wanted,
        include_partial=args.include_partial,
    )
    if args.limit:
        targets = targets[: args.limit]

    if not targets:
        print("No papers with unfetched supplement links found.")
        return 0

    print(
        f"{len(targets)} paper(s) with recorded supplement links not on disk "
        f"(corpus: {args.corpus})\n"
    )
    fetchable, _unresolved = _print_plan(targets)

    if args.dry_run:
        print("\n--dry-run: nothing fetched.")
        return 0
    if not fetchable:
        print("\nNothing resolvable; not fetching.")
        return 0

    # Constructed once. Publisher API keys are optional here - every resolvable
    # link in the backlog is a PMC bin/ URL, which needs no credentials.
    from harvesting.orchestrator import PMCHarvester

    harvester = PMCHarvester(output_dir=str(scratch_dir))

    rows: list[dict] = []
    print()
    for idx, target in enumerate(targets, 1):
        report = plan(target)
        if not report.fetchable and not target.pmcid:
            continue
        print(
            f"[{idx}/{len(targets)}] {target.gene} {target.pmid} "
            f"({target.pmcid or 'no pmcid'}): {len(report.fetchable)} linked file(s)"
        )
        try:
            row = fetch_paper(harvester, target, report)
        except Exception as exc:  # keep the sweep going on one bad paper
            row = {
                "gene": target.gene,
                "pmid": target.pmid,
                "pmcid": target.pmcid or "",
                "planned": len(report.fetchable),
                "downloaded": 0,
                "unresolved": len(report.unresolved),
                "route": "",
                "folded": False,
                "new_files": [],
                "files_on_disk": 0,
                "full_context": str(
                    target.paper_dir / f"{target.pmid}_FULL_CONTEXT.md"
                ),
                "error": f"{type(exc).__name__}: {exc}",
            }
        rows.append(row)
        if row["downloaded"]:
            print(
                f"    OK {row['files_on_disk']} file(s) via {row['route']}"
                f"{', folded' if row['folded'] else ', NOT folded'}"
            )
        else:
            print(f"    -- nothing downloaded ({row['error'] or 'no reason recorded'})")
        time.sleep(NCBI_SLEEP_SECONDS)

    papers_with_files = sum(1 for r in rows if r["downloaded"])
    total_files = sum(r["files_on_disk"] for r in rows)
    folded = sum(1 for r in rows if r["folded"])
    print(
        f"\n{papers_with_files}/{len(rows)} paper(s) gained {total_files} "
        f"supplement file(s); {folded} FULL_CONTEXT refolded."
    )

    if args.report_out:
        _write_report(args.report_out, rows)
        print(f"Wrote {args.report_out}")
    if args.source_override_out:
        _write_source_overrides(args.source_override_out, rows)
        print(f"Wrote {args.source_override_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
