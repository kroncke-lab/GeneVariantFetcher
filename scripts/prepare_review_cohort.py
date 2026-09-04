#!/usr/bin/env python3
"""Prepare a reproducible, evidence-rich PMID cohort for Variant Browser review.

The selector operates on a completed GVF run database.  Every selected paper has
at least one extracted variant.  A configurable quota must also have at least
one *trusted* penetrance row with explicit carrier, affected, and unaffected
counts.  Source-origin and publication-decade seeding prevents a high-volume
publisher or era from consuming the entire review queue.

This command only writes local cohort artifacts.  It never connects to or
mutates Variant Browser.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils.pmid_utils import is_valid_pmid  # noqa: E402


@dataclass(frozen=True)
class PaperCandidate:
    pmid: str
    title: str
    journal: str
    publication_date: str
    year: int | None
    decade: str
    doi: str
    pmc_id: str
    source_origin: str
    has_full_text: bool
    variant_links: int
    unique_variants: int
    count_rows: int
    carrier_rows: int
    affected_rows: int
    unaffected_rows: int
    complete_rows: int
    trusted_complete_rows: int
    total_carriers: int
    total_affected: int
    total_unaffected: int
    extraction_selected: bool
    secondary_llm_adjudicated: bool
    claim_verifier_calls: int

    @property
    def complete_for_quota(self) -> bool:
        return self.trusted_complete_rows > 0


def _safe_int(value: object) -> int:
    return int(value or 0)


def _year(value: str | None) -> int | None:
    text = (value or "").strip()
    if len(text) < 4 or not text[:4].isdigit():
        return None
    return int(text[:4])


def _decade(year: int | None) -> str:
    return f"{year // 10 * 10}s" if year else "unknown"


def _artifact_metadata(run_dir: Path, pmid: str) -> tuple[str, bool]:
    path = run_dir / "pmc_fulltext" / f"{pmid}_artifacts.json"
    if not path.exists():
        return "abstract-only", False
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return "artifact-unreadable", False
    source = ((payload.get("main_text") or {}).get("source") or "").strip()
    if not source:
        return "abstract-only", False
    return source, True


def _pubmed_metadata(run_dir: Path, pmid: str) -> dict:
    path = run_dir / "abstract_json" / f"{pmid}.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text()).get("metadata") or {}
    except (OSError, json.JSONDecodeError):
        return {}


def _trace_status(run_dir: Path) -> tuple[set[str], set[str], Counter[str]]:
    selected: set[str] = set()
    adjudicated: set[str] = set()
    claim_calls: Counter[str] = Counter()
    path = run_dir / "llm_traces" / "trace_index.jsonl"
    if not path.exists():
        return selected, adjudicated, claim_calls
    for line in path.read_text().splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        context = record.get("context") or {}
        pmid = str(context.get("pmid") or "").strip()
        if not pmid:
            continue
        stage = context.get("stage")
        if (
            record.get("record_type") == "decision_event"
            and stage == "paper_variant_extraction"
        ):
            selected.add(pmid)
        elif (
            record.get("record_type") == "llm_call"
            and stage == "paper_extraction_adjudication"
        ):
            adjudicated.add(pmid)
        elif (
            record.get("record_type") == "llm_call"
            and stage == "variant_claim_verification"
        ):
            claim_calls[pmid] += 1
    return selected, adjudicated, claim_calls


def load_candidates(run_dir: Path, gene: str) -> list[PaperCandidate]:
    db = run_dir / f"{gene}.db"
    if not db.exists():
        raise FileNotFoundError(f"run database not found: {db}")
    selected, adjudicated, claim_calls = _trace_status(run_dir)
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            """
            SELECT p.pmid, p.title, p.journal, p.publication_date, p.doi, p.pmc_id,
                   COUNT(*) AS variant_links,
                   COUNT(DISTINCT vp.variant_id) AS unique_variants,
                   SUM(CASE WHEN pd.penetrance_id IS NOT NULL THEN 1 ELSE 0 END) AS count_rows,
                   SUM(CASE WHEN pd.total_carriers_observed IS NOT NULL THEN 1 ELSE 0 END) AS carrier_rows,
                   SUM(CASE WHEN pd.affected_count IS NOT NULL THEN 1 ELSE 0 END) AS affected_rows,
                   SUM(CASE WHEN pd.unaffected_count IS NOT NULL THEN 1 ELSE 0 END) AS unaffected_rows,
                   SUM(CASE WHEN pd.total_carriers_observed > 0
                                  AND pd.affected_count IS NOT NULL
                                  AND pd.unaffected_count IS NOT NULL
                                  AND pd.total_carriers_observed >=
                                      pd.affected_count + pd.unaffected_count
                            THEN 1 ELSE 0 END) AS complete_rows,
                   SUM(CASE WHEN pd.trust_tier = 'trusted'
                                  AND pd.total_carriers_observed > 0
                                  AND pd.affected_count IS NOT NULL
                                  AND pd.unaffected_count IS NOT NULL
                                  AND pd.total_carriers_observed >=
                                      pd.affected_count + pd.unaffected_count
                            THEN 1 ELSE 0 END) AS trusted_complete_rows,
                   SUM(COALESCE(pd.total_carriers_observed, 0)) AS total_carriers,
                   SUM(COALESCE(pd.affected_count, 0)) AS total_affected,
                   SUM(COALESCE(pd.unaffected_count, 0)) AS total_unaffected
              FROM papers p
              JOIN variant_papers vp ON vp.pmid = p.pmid
              JOIN variants v ON v.variant_id = vp.variant_id
         LEFT JOIN penetrance_data pd
                ON pd.variant_id = vp.variant_id AND pd.pmid = vp.pmid
             WHERE UPPER(v.gene_symbol) = ?
          GROUP BY p.pmid, p.title, p.journal, p.publication_date, p.doi, p.pmc_id
            """,
            (gene.upper(),),
        ).fetchall()
    finally:
        con.close()

    candidates: list[PaperCandidate] = []
    for row in rows:
        pmid = str(row["pmid"] or "").strip()
        metadata = _pubmed_metadata(run_dir, pmid)
        publication_date = str(row["publication_date"] or metadata.get("year") or "")
        year = _year(publication_date)
        source, has_full_text = _artifact_metadata(run_dir, pmid)
        candidates.append(
            PaperCandidate(
                pmid=pmid,
                title=str(metadata.get("title") or row["title"] or ""),
                journal=str(metadata.get("journal") or row["journal"] or ""),
                publication_date=publication_date,
                year=year,
                decade=_decade(year),
                doi=str(row["doi"] or ""),
                pmc_id=str(row["pmc_id"] or ""),
                source_origin=source,
                has_full_text=has_full_text,
                variant_links=_safe_int(row["variant_links"]),
                unique_variants=_safe_int(row["unique_variants"]),
                count_rows=_safe_int(row["count_rows"]),
                carrier_rows=_safe_int(row["carrier_rows"]),
                affected_rows=_safe_int(row["affected_rows"]),
                unaffected_rows=_safe_int(row["unaffected_rows"]),
                complete_rows=_safe_int(row["complete_rows"]),
                trusted_complete_rows=_safe_int(row["trusted_complete_rows"]),
                total_carriers=_safe_int(row["total_carriers"]),
                total_affected=_safe_int(row["total_affected"]),
                total_unaffected=_safe_int(row["total_unaffected"]),
                extraction_selected=pmid in selected,
                secondary_llm_adjudicated=pmid in adjudicated,
                claim_verifier_calls=claim_calls[pmid],
            )
        )
    return candidates


def _rank(candidate: PaperCandidate) -> tuple:
    """Evidence-rich deterministic rank; ascending PMID is the final tie-break."""

    pmid_num = int(candidate.pmid) if candidate.pmid.isdigit() else math.inf
    return (
        -int(candidate.secondary_llm_adjudicated),
        -candidate.trusted_complete_rows,
        -candidate.complete_rows,
        -candidate.carrier_rows,
        -candidate.total_carriers,
        -candidate.unique_variants,
        pmid_num,
        candidate.pmid,
    )


def _balanced_pick(
    pool: Iterable[PaperCandidate],
    count: int,
    *,
    already_selected: Iterable[PaperCandidate] = (),
) -> list[PaperCandidate]:
    available = sorted(pool, key=_rank)
    picked: list[PaperCandidate] = []
    used = {paper.pmid for paper in already_selected}

    def add(paper: PaperCandidate | None) -> None:
        if paper is not None and paper.pmid not in used and len(picked) < count:
            picked.append(paper)
            used.add(paper.pmid)

    # Preserve every applicable secondary-adjudicator case: these are scarce and
    # particularly useful for comparing model and human judgments.
    for paper in available:
        if paper.secondary_llm_adjudicated:
            add(paper)

    # Seed every available acquisition source and publication decade before
    # filling by evidence richness.
    for source in sorted({paper.source_origin for paper in available}):
        add(next((paper for paper in available if paper.source_origin == source), None))
    for decade in sorted({paper.decade for paper in available}):
        add(next((paper for paper in available if paper.decade == decade), None))
    for paper in available:
        add(paper)
    return picked


def select_cohort(
    candidates: list[PaperCandidate],
    *,
    size: int,
    complete_count: int,
    min_source_origins: int,
    min_decades: int,
) -> list[PaperCandidate]:
    if complete_count > size:
        raise ValueError("complete-count cannot exceed cohort size")
    complete_pool = [paper for paper in candidates if paper.complete_for_quota]
    if len(complete_pool) < complete_count:
        raise ValueError(
            f"need {complete_count} trusted complete papers, found {len(complete_pool)}"
        )
    complete = _balanced_pick(complete_pool, complete_count)
    remainder_pool = [
        paper for paper in candidates if paper.pmid not in {p.pmid for p in complete}
    ]
    remainder = _balanced_pick(
        remainder_pool,
        size - complete_count,
        already_selected=complete,
    )
    cohort = complete + remainder

    if len(cohort) != size or len({paper.pmid for paper in cohort}) != size:
        raise ValueError(f"expected {size} unique papers, selected {len(cohort)}")
    invalid = sorted({paper.pmid for paper in cohort if not is_valid_pmid(paper.pmid)})
    if invalid:
        raise ValueError(f"selected non-PMID identifiers: {invalid}")
    if any(paper.variant_links < 1 for paper in cohort):
        raise ValueError("every selected paper must have at least one variant link")
    actual_complete = sum(paper.complete_for_quota for paper in cohort)
    if actual_complete < complete_count:
        raise ValueError(
            f"complete quota regressed: {actual_complete} < {complete_count}"
        )
    sources = {paper.source_origin for paper in cohort}
    decades = {paper.decade for paper in cohort}
    if len(sources) < min_source_origins:
        raise ValueError(
            f"need {min_source_origins} source origins, selected {len(sources)}"
        )
    if len(decades) < min_decades:
        raise ValueError(f"need {min_decades} decades, selected {len(decades)}")
    return cohort


def write_outputs(
    output_dir: Path,
    gene: str,
    cohort: list[PaperCandidate],
    *,
    run_dir: Path,
    complete_count: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pmid_path = output_dir / f"{gene}.txt"
    csv_path = output_dir / f"{gene}_manifest.csv"
    summary_path = output_dir / "summary.json"
    readme_path = output_dir / "README.md"

    pmid_path.write_text("\n".join(paper.pmid for paper in cohort) + "\n")
    fields = [
        "review_order",
        "pmid",
        "year",
        "decade",
        "source_origin",
        "has_full_text",
        "variant_links",
        "unique_variants",
        "count_rows",
        "carrier_rows",
        "affected_rows",
        "unaffected_rows",
        "complete_rows",
        "trusted_complete_rows",
        "complete_for_quota",
        "total_carriers",
        "total_affected",
        "total_unaffected",
        "extraction_selected",
        "secondary_llm_adjudicated",
        "claim_verifier_calls",
        "publication_date",
        "journal",
        "title",
        "doi",
        "pmc_id",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, paper in enumerate(cohort, 1):
            record = asdict(paper)
            record["review_order"] = index
            record["complete_for_quota"] = paper.complete_for_quota
            writer.writerow({field: record.get(field, "") for field in fields})

    source_counts = Counter(paper.source_origin for paper in cohort)
    decade_counts = Counter(paper.decade for paper in cohort)
    actual_complete = sum(paper.complete_for_quota for paper in cohort)
    summary = {
        "gene": gene,
        "run_dir": str(run_dir),
        "cohort_size": len(cohort),
        "selection_contract": {
            "all_papers_have_extracted_variants": True,
            "complete_paper_definition": (
                "at least one trusted row with total_carriers_observed > 0, "
                "non-null affected_count and unaffected_count, and total >= affected + unaffected"
            ),
            "minimum_complete_papers": complete_count,
        },
        "observed": {
            "complete_papers": actual_complete,
            "complete_fraction": actual_complete / len(cohort),
            "papers_with_carrier_rows": sum(paper.carrier_rows > 0 for paper in cohort),
            "secondary_llm_adjudicated_papers": sum(
                paper.secondary_llm_adjudicated for paper in cohort
            ),
            "source_origins": dict(sorted(source_counts.items())),
            "publication_decades": dict(sorted(decade_counts.items())),
            "publication_year_min": min(paper.year for paper in cohort if paper.year),
            "publication_year_max": max(paper.year for paper in cohort if paper.year),
            "variant_links": sum(paper.variant_links for paper in cohort),
            "unique_variant_paper_counts": sum(
                paper.unique_variants for paper in cohort
            ),
            "trusted_complete_rows": sum(
                paper.trusted_complete_rows for paper in cohort
            ),
        },
        "important_distinction": (
            "This is a queue for Variant Browser human adjudication. The secondary_llm_adjudicated "
            "field records the narrower subset that triggered GVF's second-reader adjudicator."
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    readme_path.write_text(
        f"""# {gene} 50-paper Variant Browser review cohort

This pinned cohort was selected from `{run_dir}` for private human adjudication
in Variant Browser. It is not a gold standard.

## Contract

- All {len(cohort)} papers have at least one extracted `{gene}` variant.
- {actual_complete}/{len(cohort)} papers ({actual_complete / len(cohort):.0%}) have at least one
  trusted variant row with explicit carrier, affected, and unaffected counts.
- Publication years span {summary["observed"]["publication_year_min"]}–{summary["observed"]["publication_year_max"]}.
- Source origins: {", ".join(f"{key}={value}" for key, value in sorted(source_counts.items()))}.
- The manifest distinguishes pipeline extraction selection, secondary LLM adjudication,
  and the human adjudication that will occur after publishing.

## Prepared staging commands

Run this only after BMPR2 variantFeatures readiness is verified and the EZproxy-backed
source gaps have been handled or explicitly accepted:

```bash
cd /Users/kronckbm/GitRepos/Variant_Browser
set -a; source .env; set +a
venv/bin/python manage.py import_features \\
  --gene BMPR2 \\
  --disease 'pulmonary arterial hypertension' \\
  --source local \\
  --db /Users/kronckbm/GitRepos/variantFeatures/data/variants.db \\
  --database staging
```

That feature import creates the private gene–disease pair and its first snapshot.
Then attach the pinned GVF paper cohort:

```bash
GVF_PMID_FILE={pmid_path.resolve()} \\
GVF_DATASET_LABEL=bmpr2_review_50_20260808 \\
GVF_DATASET_NOTE='BMPR2/PAH initial 50-paper human-review cohort' \\
bash /Users/kronckbm/GitRepos/Variant_Browser/scripts/gvf_publish.sh \\
  BMPR2 {run_dir.resolve() / "BMPR2.db"} 'pulmonary arterial hypertension'
```

Both commands write only to the private Variant Browser staging/review database.
The gene–disease pair does not currently exist. If a feature import is deliberately
deferred, add `GVF_CREATE_PAIR=1` to the carrier command to create a carrier-only
snapshot; the preferred path is features first. Do not set `GVF_FULL_DB_REPLACE`
for this cohort.
"""
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--size", type=int, default=50)
    parser.add_argument(
        "--complete-fraction",
        type=float,
        default=0.75,
        help="minimum fraction with trusted carrier+affected+unaffected rows",
    )
    parser.add_argument(
        "--complete-count",
        type=int,
        default=None,
        help="explicit complete-paper target (overrides --complete-fraction)",
    )
    parser.add_argument("--min-source-origins", type=int, default=4)
    parser.add_argument("--min-decades", type=int, default=3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gene = args.gene.upper()
    if args.size < 1:
        raise SystemExit("--size must be positive")
    if not 0 < args.complete_fraction <= 1:
        raise SystemExit("--complete-fraction must be in (0, 1]")
    complete_count = (
        args.complete_count
        if args.complete_count is not None
        else math.ceil(args.size * args.complete_fraction)
    )
    candidates = load_candidates(args.run_dir, gene)
    cohort = select_cohort(
        candidates,
        size=args.size,
        complete_count=complete_count,
        min_source_origins=args.min_source_origins,
        min_decades=args.min_decades,
    )
    write_outputs(
        args.output_dir,
        gene,
        cohort,
        run_dir=args.run_dir,
        complete_count=complete_count,
    )
    print(
        f"{gene}: selected {len(cohort)} papers; "
        f"{sum(p.complete_for_quota for p in cohort)} complete; "
        f"{len({p.source_origin for p in cohort})} sources; "
        f"{len({p.decade for p in cohort})} decades"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
