#!/usr/bin/env python3
"""Publish the curated extraction-eval papers into Variant_Browser staging.

This is intentionally a staging/review helper, not the full-gene publish path.
It publishes the fixed curated paper set using the per-gene
``results/<GENE>/review_staging_test/<GENE>.db`` SQLite artifacts and the
generated ``benchmarks/curated_extraction_eval/pmids/<GENE>.txt`` manifests.
"""

from __future__ import annotations

import argparse
import csv
import os
import sqlite3
import subprocess
import sys
from datetime import datetime
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
BENCH = REPO / "benchmarks" / "curated_extraction_eval"
MANIFEST = BENCH / "manifest.csv"
PMID_DIR = BENCH / "pmids"
DEFAULT_REVIEW_REPO = REPO.parent / "Variant_Browser"

DISEASE_DEFAULTS = {
    "APOE": "APOE curated literature review",
    "BRCA1": "Hereditary breast and ovarian cancer",
    "BRCA2": "Hereditary breast and ovarian cancer",
    "KCNH2": "Long QT type 2",
    "KCNQ1": "Long QT type 1",
    "MYBPC3": "Hypertrophic cardiomyopathy",
    "RYR2": "Catecholaminergic polymorphic ventricular tachycardia",
    "SCN5A": "Brugada syndrome / Long QT type 3",
}


def load_manifest() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with MANIFEST.open(newline="") as f:
        for row in csv.DictReader(f):
            gene = row["gene"].strip().upper()
            pmid = row["pmid"].strip()
            out.setdefault(gene, []).append(pmid)
    return {gene: list(dict.fromkeys(pmids)) for gene, pmids in out.items()}


def parse_db_overrides(raw: list[str] | None) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in raw or []:
        if "=" not in item:
            raise SystemExit(f"--db must be GENE=/path/to.db, got {item!r}")
        gene, path = item.split("=", 1)
        gene = gene.strip().upper()
        if not gene:
            raise SystemExit(f"--db has empty gene in {item!r}")
        out[gene] = Path(path).expanduser()
    return out


def latest_gene_db(root: Path, gene: str) -> Path | None:
    """Return the newest DB for a gene under a benchmark extraction root."""
    search_roots = [root / gene, root]
    candidates: list[Path] = []
    seen: set[Path] = set()
    for search_root in search_roots:
        if not search_root.exists():
            continue
        for path in search_root.glob(f"**/{gene}*.db"):
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                candidates.append(path)
    if not candidates:
        return None
    return sorted(candidates, key=lambda p: (p.stat().st_mtime, str(p)))[-1]


def default_db_path_for(gene: str) -> Path:
    return REPO / "results" / gene / "review_staging_test" / f"{gene}.db"


def db_path_for(
    gene: str,
    *,
    db_overrides: dict[str, Path],
    extract_root: Path | None,
) -> Path:
    if gene in db_overrides:
        return db_overrides[gene]
    if extract_root is not None:
        found = latest_gene_db(extract_root, gene)
        if found is not None:
            return found
    return default_db_path_for(gene)


def sqlite_pmids(db_path: Path) -> set[str]:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        return {str(row[0]).strip() for row in con.execute("SELECT pmid FROM papers")}
    finally:
        con.close()


def count_rows(db_path: Path) -> dict[str, int]:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        out = {}
        for table in (
            "papers",
            "variant_papers",
            "penetrance_data",
            "individual_records",
            "phenotypes",
        ):
            try:
                out[table] = int(
                    con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                )
            except sqlite3.Error:
                out[table] = 0
        return out
    finally:
        con.close()


def parse_genes(raw: str | None, available: list[str]) -> list[str]:
    if not raw:
        return available
    wanted = [g.strip().upper() for g in raw.split(",") if g.strip()]
    bad = sorted(set(wanted) - set(available))
    if bad:
        raise SystemExit(f"Unknown gene(s) for curated manifest: {', '.join(bad)}")
    return wanted


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--review-repo", type=Path, default=DEFAULT_REVIEW_REPO)
    ap.add_argument(
        "--genes",
        help="comma-separated gene list; default is every gene in manifest.csv",
    )
    ap.add_argument(
        "--db",
        action="append",
        help="Use this DB for publish: GENE=/path/to.db (repeatable).",
    )
    ap.add_argument(
        "--extract-root",
        type=Path,
        help="Benchmark extract root containing per-gene GVF outputs; newest **/GENE*.db is used.",
    )
    ap.add_argument(
        "--run-root",
        type=Path,
        help="Benchmark validation run root; equivalent to --extract-root RUN_ROOT/extract when present.",
    )
    ap.add_argument(
        "--dataset-label", help="staging dataset slug stored by Variant_Browser"
    )
    ap.add_argument(
        "--dataset-note", default="101-paper curated GVF benchmark staging review"
    )
    ap.add_argument(
        "--create-pairs",
        action="store_true",
        help="allow Variant_Browser to create missing staging pairs/snapshots",
    )
    ap.add_argument(
        "--allow-missing-pmids",
        action="store_true",
        help="warn instead of failing if a DB lacks a manifest PMID",
    )
    ap.add_argument("--continue-on-error", action="store_true")
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print publish commands without mutating staging",
    )
    args = ap.parse_args()

    manifest = load_manifest()
    genes = parse_genes(args.genes, sorted(manifest))
    total_pmids = sum(len(manifest[g]) for g in genes)
    dataset_label = (
        args.dataset_label
        or f"gvf_curated_{total_pmids}_{datetime.now().strftime('%Y%m%d')}"
    )
    db_overrides = parse_db_overrides(args.db)
    extract_root = args.extract_root
    if args.run_root:
        run_extract = args.run_root / "extract"
        extract_root = run_extract if run_extract.exists() else args.run_root

    publish_script = args.review_repo / "scripts" / "gvf_publish.sh"
    if not publish_script.exists():
        raise SystemExit(f"Variant_Browser publish script not found: {publish_script}")

    print(f"Dataset label: {dataset_label}")
    print(f"Genes: {', '.join(genes)} ({total_pmids} manifest PMIDs)")
    failures: list[str] = []

    for gene in genes:
        db_path = db_path_for(
            gene, db_overrides=db_overrides, extract_root=extract_root
        )
        pmid_file = PMID_DIR / f"{gene}.txt"
        if not db_path.exists():
            raise SystemExit(f"{gene}: DB not found: {db_path}")
        if not pmid_file.exists():
            raise SystemExit(f"{gene}: PMID file not found: {pmid_file}")

        db_pmids = sqlite_pmids(db_path)
        missing = [pmid for pmid in manifest[gene] if pmid not in db_pmids]
        if missing and not args.allow_missing_pmids:
            raise SystemExit(
                f"{gene}: DB is missing manifest PMID(s): {', '.join(missing)}"
            )
        rows = count_rows(db_path)
        print(
            f"{gene}: {len(manifest[gene])} PMIDs, "
            f"db papers={rows['papers']}, variant_papers={rows['variant_papers']}, "
            f"penetrance={rows['penetrance_data']}, individuals={rows['individual_records']}, "
            f"phenotypes={rows['phenotypes']}"
        )

        cmd = [
            str(publish_script),
            gene,
            str(db_path),
            DISEASE_DEFAULTS.get(gene, f"{gene} curated review"),
        ]
        env = os.environ.copy()
        env.update(
            {
                "GVF_PMID_FILE": str(pmid_file),
                "GVF_DATASET_LABEL": dataset_label,
                "GVF_DATASET_NOTE": args.dataset_note,
                "GVF_CREATE_PAIR": "1" if args.create_pairs else "",
                "GVF_ALLOW_MISSING_PMIDS": "1" if args.allow_missing_pmids else "",
            }
        )
        if args.dry_run:
            print("  DRY RUN:", " ".join(cmd))
            continue
        res = subprocess.run(cmd, cwd=str(args.review_repo), env=env)
        if res.returncode != 0:
            msg = f"{gene}: publish failed with exit {res.returncode}"
            if not args.continue_on_error:
                raise SystemExit(msg)
            failures.append(msg)

    if failures:
        print("Failures:", file=sys.stderr)
        for msg in failures:
            print(f"  {msg}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
