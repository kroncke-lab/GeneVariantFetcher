#!/usr/bin/env python3
"""Build deterministic pilot subsets from the GVF gold-standard package.

The output mirrors the top-level gold-standard layout enough that
``scripts/run_recall_suite.py --gold-dir gene_variant_fetcher_gold_standard/pilots``
can score against the pilot rows. The pilot manifest keeps the full clinical
counts linked for audit/provenance while the recall inputs stay in the reduced
GVF scoring shape.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable


BASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_GOLD_DIR = BASE_DIR / "gene_variant_fetcher_gold_standard"
DEFAULT_OUT_DIR = DEFAULT_GOLD_DIR / "pilots"
DEFAULT_GENES = ("KCNH2", "KCNQ1", "SCN5A")
EDGE_CASE_GENES = ("RYR2",)

# Existing integration-test PMIDs in tests/test_all_integrations.py. If a gene
# has any of these in its gold input, include one so API-client smoke failures
# can be compared against the older hand-verified test set.
KNOWN_API_PMIDS = ("24667783", "19841300", "20173333", "30036649", "19716085")

RECALL_COLUMNS = ("variant", "pmid", "carriers", "affected", "unaffected")
PILOT_COLUMNS = (
    "gene",
    "pilot_case",
    "pmid",
    "variant_rows",
    "unique_variants",
    "carriers",
    "affected",
    "unaffected",
    "zero_carrier_rows",
    "clinical_rows",
    "example_variants",
)


@dataclass
class PmidSummary:
    """Aggregated gold-standard counts for one gene/PMID."""

    pmid: str
    variant_rows: int = 0
    variants: set[str] = field(default_factory=set)
    carriers: int = 0
    affected: int = 0
    unaffected: int = 0
    zero_carrier_rows: int = 0
    clinical_rows: int = 0

    @property
    def unique_variants(self) -> int:
        return len(self.variants)

    @property
    def example_variants(self) -> str:
        return ";".join(sorted(self.variants)[:5])


def parse_int(value: object) -> int:
    """Parse a clinical count cell, treating blank/non-numeric values as zero."""
    if value is None:
        return 0
    text = str(value).strip()
    if not text:
        return 0
    try:
        return int(float(text))
    except ValueError:
        return 0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, rows: Iterable[dict[str, object]], fieldnames: tuple[str, ...]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def summarize_recall_rows(rows: list[dict[str, str]]) -> dict[str, PmidSummary]:
    summaries: dict[str, PmidSummary] = {}
    for row in rows:
        pmid = str(row.get("pmid", "")).strip()
        if not pmid:
            continue
        summary = summaries.setdefault(pmid, PmidSummary(pmid=pmid))
        summary.variant_rows += 1
        variant = str(row.get("variant", "")).strip()
        if variant:
            summary.variants.add(variant)
        carriers = parse_int(row.get("carriers"))
        summary.carriers += carriers
        summary.affected += parse_int(row.get("affected"))
        summary.unaffected += parse_int(row.get("unaffected"))
        if carriers == 0:
            summary.zero_carrier_rows += 1
    return summaries


def add_clinical_row_counts(
    summaries: dict[str, PmidSummary], clinical_rows: list[dict[str, str]]
) -> None:
    for row in clinical_rows:
        pmid = str(row.get("pmid", "")).strip()
        if pmid in summaries:
            summaries[pmid].clinical_rows += 1


def first_matching(
    summaries: dict[str, PmidSummary],
    predicate: Callable[[PmidSummary], bool],
    key: Callable[[PmidSummary], tuple],
    selected: set[str],
) -> PmidSummary | None:
    candidates = [
        summary
        for summary in summaries.values()
        if summary.pmid not in selected and predicate(summary)
    ]
    if not candidates:
        return None
    return sorted(candidates, key=key)[0]


def select_pilot_pmids(
    summaries: dict[str, PmidSummary], pmids_per_gene: int
) -> list[tuple[str, PmidSummary]]:
    """Select reproducible, behaviorally diverse PMIDs for one gene."""
    selected: list[tuple[str, PmidSummary]] = []
    seen: set[str] = set()

    case_specs: list[
        tuple[str, Callable[[PmidSummary], bool], Callable[[PmidSummary], tuple]]
    ] = [
        (
            "heavy_table",
            lambda item: True,
            lambda item: (-item.variant_rows, -item.unique_variants, item.pmid),
        ),
        (
            "affected_and_unaffected",
            lambda item: item.affected > 0 and item.unaffected > 0,
            lambda item: (
                -(item.affected + item.unaffected),
                -item.variant_rows,
                item.pmid,
            ),
        ),
        (
            "unaffected_only",
            lambda item: item.affected == 0 and item.unaffected > 0,
            lambda item: (-item.unaffected, -item.variant_rows, item.pmid),
        ),
        (
            "single_variant",
            lambda item: item.variant_rows == 1 and item.carriers > 0,
            lambda item: (-item.carriers, item.pmid),
        ),
        (
            "zero_carrier_row",
            lambda item: item.zero_carrier_rows > 0,
            lambda item: (-item.zero_carrier_rows, -item.variant_rows, item.pmid),
        ),
        (
            "existing_api_gold",
            lambda item: item.pmid in KNOWN_API_PMIDS,
            lambda item: (KNOWN_API_PMIDS.index(item.pmid), item.pmid),
        ),
    ]

    for case_name, predicate, key in case_specs:
        if len(selected) >= pmids_per_gene:
            break
        match = first_matching(summaries, predicate, key, seen)
        if match is None:
            continue
        selected.append((case_name, match))
        seen.add(match.pmid)

    fillers = sorted(
        (item for item in summaries.values() if item.pmid not in seen),
        key=lambda item: (-item.variant_rows, -item.unique_variants, item.pmid),
    )
    for item in fillers:
        if len(selected) >= pmids_per_gene:
            break
        selected.append(("coverage_fill", item))
        seen.add(item.pmid)

    return selected


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open() as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def write_readme(path: Path, genes: list[str], pmids_per_gene: int) -> None:
    gene_list = ", ".join(genes)
    text = f"""# Gold-Standard Pilot Set

Deterministic pilot subsets generated from `gene_variant_fetcher_gold_standard`.
The default clean pilot genes are `{gene_list}` with up to `{pmids_per_gene}`
PMIDs per gene.

## What To Compare Against

- Use `normalized/<GENE>_recall_input.csv` here for scoring pilots. These files
  have the same reduced schema as the full recall inputs: `variant`, `pmid`,
  `carriers`, `affected`, `unaffected`.
- Use the parent package's `normalized/<GENE>_clinical_counts.csv` for audit and
  provenance. Those files retain `source_row_id`, `source_row_ordinal`,
  `source_type`, raw/normalized variants, disease-specific affected columns,
  unaffected counts, and QC flags.
- `pilot_pmids.csv` records why each PMID was selected and gives aggregate
  counts for quick API smoke-test triage.

## Full vs Partial

The parent KCNH2, KCNQ1, and SCN5A `clinical_counts` files are full flat
snapshots of their exported source tables. The recall inputs are intentionally
partial: PubMed-indexed scoring rows only, collapsed to the counts GVF compares.
RYR2 is excluded by default because it is xlsx-derived and mixes literature rows
with population-only gnomAD rows; add it with `--include-ryr2` when testing that
edge case specifically.

## Regenerate

```bash
python scripts/build_gold_standard_pilots.py
```

Run live API smoke checks only when the local Python environment has GVF
dependencies installed:

```bash
python scripts/run_gold_standard_api_pilots.py --live
```
"""
    path.write_text(text)


def build_pilots(
    *,
    gold_dir: Path,
    out_dir: Path,
    genes: list[str],
    pmids_per_gene: int,
) -> dict:
    normalized_dir = gold_dir / "normalized"
    out_normalized = out_dir / "normalized"
    out_normalized.mkdir(parents=True, exist_ok=True)

    manifest = load_json(gold_dir / "manifest.json")
    pilot_rows: list[dict[str, object]] = []
    manifest_genes: dict[str, dict[str, object]] = {}

    for gene in genes:
        recall_path = normalized_dir / f"{gene}_recall_input.csv"
        clinical_path = normalized_dir / f"{gene}_clinical_counts.csv"
        if not recall_path.exists():
            raise FileNotFoundError(f"Missing recall input for {gene}: {recall_path}")
        if not clinical_path.exists():
            raise FileNotFoundError(
                f"Missing clinical counts for {gene}: {clinical_path}"
            )

        recall_rows = read_csv(recall_path)
        clinical_rows = read_csv(clinical_path)
        summaries = summarize_recall_rows(recall_rows)
        add_clinical_row_counts(summaries, clinical_rows)
        selected = select_pilot_pmids(summaries, pmids_per_gene)
        selected_pmids = {summary.pmid for _, summary in selected}

        subset_rows = [
            {column: row.get(column, "") for column in RECALL_COLUMNS}
            for row in recall_rows
            if row.get("pmid") in selected_pmids
        ]
        write_csv(
            out_normalized / f"{gene}_recall_input.csv", subset_rows, RECALL_COLUMNS
        )

        for case_name, summary in selected:
            pilot_rows.append(
                {
                    "gene": gene,
                    "pilot_case": case_name,
                    "pmid": summary.pmid,
                    "variant_rows": summary.variant_rows,
                    "unique_variants": summary.unique_variants,
                    "carriers": summary.carriers,
                    "affected": summary.affected,
                    "unaffected": summary.unaffected,
                    "zero_carrier_rows": summary.zero_carrier_rows,
                    "clinical_rows": summary.clinical_rows,
                    "example_variants": summary.example_variants,
                }
            )

        manifest_genes[gene] = {
            "source_recall_input": str(recall_path.relative_to(gold_dir)),
            "source_clinical_counts": str(clinical_path.relative_to(gold_dir)),
            "selected_pmid_count": len(selected),
            "pilot_recall_row_count": len(subset_rows),
            "pilot_pmids": [summary.pmid for _, summary in selected],
        }

    write_csv(out_dir / "pilot_pmids.csv", pilot_rows, PILOT_COLUMNS)
    write_readme(out_dir / "README.md", genes, pmids_per_gene)

    pilot_manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator_script": "scripts/build_gold_standard_pilots.py",
        "source_gold_dir": str(gold_dir),
        "source_manifest_generated_at_utc": manifest.get("generated_at_utc"),
        "source_manifest_commit": manifest.get("gvf_git_commit"),
        "pmids_per_gene_requested": pmids_per_gene,
        "genes": manifest_genes,
        "notes": [
            "Recall pilot files are PubMed-only scoring subsets.",
            "Clinical-count source files remain in the parent gold-standard package.",
            "RYR2 is excluded by default because its source is xlsx-derived, not Variant_Browser per-PMID SQL.",
        ],
    }
    write_json(out_dir / "manifest.json", pilot_manifest)
    return pilot_manifest


def parse_genes(value: str | None, include_ryr2: bool) -> list[str]:
    if value:
        genes = [item.strip().upper() for item in value.split(",") if item.strip()]
    else:
        genes = list(DEFAULT_GENES)
    if include_ryr2:
        for gene in EDGE_CASE_GENES:
            if gene not in genes:
                genes.append(gene)
    return genes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gold-dir",
        type=Path,
        default=DEFAULT_GOLD_DIR,
        help="Gold-standard package root.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Pilot package output directory.",
    )
    parser.add_argument(
        "--genes",
        help="Comma-separated gene list. Default: KCNH2,KCNQ1,SCN5A.",
    )
    parser.add_argument(
        "--include-ryr2",
        action="store_true",
        help="Include the xlsx-derived RYR2 edge-case pilot.",
    )
    parser.add_argument(
        "--pmids-per-gene",
        type=int,
        default=10,
        help="Maximum PMIDs to select per gene.",
    )
    args = parser.parse_args()

    if args.pmids_per_gene < 1:
        parser.error("--pmids-per-gene must be at least 1")

    genes = parse_genes(args.genes, args.include_ryr2)
    manifest = build_pilots(
        gold_dir=args.gold_dir.expanduser(),
        out_dir=args.out_dir.expanduser(),
        genes=genes,
        pmids_per_gene=args.pmids_per_gene,
    )
    for gene, info in manifest["genes"].items():
        print(
            f"{gene}: {info['selected_pmid_count']} PMIDs, "
            f"{info['pilot_recall_row_count']} recall rows"
        )
    print(f"Wrote pilot package: {args.out_dir.expanduser()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
