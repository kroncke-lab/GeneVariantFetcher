#!/usr/bin/env python3
"""Add adjudicated gold-v2 columns for manually reviewed recall rows.

The original `carriers`, `affected`, and `unaffected` columns are preserved.
This script only populates additive `gold_v2_*` columns plus provenance fields
for the small adjudicated set from the 2026-05 paper-failure probe.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

try:
    from scripts.recall_audit.common import repo_path
except ModuleNotFoundError:  # pragma: no cover
    from common import repo_path


V2_COLUMNS = [
    "gold_v2_carriers",
    "gold_v2_affected",
    "gold_v2_unaffected",
    "gold_v2_status",
    "gold_v2_note",
    "gold_v2_source",
]

ADJUDICATION_SOURCE = (
    "validation_runs/paper_failure_probe_20260520/gold_standard_adjudication.md"
)

ADJUDICATIONS: dict[tuple[str, str, str], dict[str, str]] = {
    ("KCNH2", "10086971", "S818L"): {
        "gold_v2_carriers": "2",
        "gold_v2_affected": "2",
        "gold_v2_unaffected": "0",
        "gold_v2_status": "schema_policy_paper_defined_affected",
        "gold_v2_note": (
            "Daughter carries S818L and QTc 472 ms meets the paper's affected "
            "definition; original 2/1/1 is retained for symptomatic-only policy."
        ),
    },
    ("KCNH2", "24606995", "A561V"): {
        "gold_v2_carriers": "1",
        "gold_v2_affected": "1",
        "gold_v2_unaffected": "0",
        "gold_v2_status": "confirmed_original",
        "gold_v2_note": "A561V table row supports one LQTS proband; model read of total=2 was table-layout confusion.",
    },
    ("RYR2", "30403697", "c.3599-9delT"): {
        "gold_v2_carriers": "1",
        "gold_v2_affected": "1",
        "gold_v2_unaffected": "",
        "gold_v2_status": "adjudicated_null_unaffected",
        "gold_v2_note": (
            "Variant is listed for subject 14 only; no explicit unaffected "
            "carrier is supported for this intronic VUS."
        ),
    },
    ("RYR2", "33686871", "G3118R"): {
        "gold_v2_carriers": "17",
        "gold_v2_affected": "6",
        "gold_v2_unaffected": "11",
        "gold_v2_status": "confirmed_original",
        "gold_v2_note": "Gold includes four severe cases plus two additional phenotype-positive siblings and pedigree-derived carriers.",
    },
    ("SCN5A", "14961552", "E1225K"): {
        "gold_v2_carriers": "3",
        "gold_v2_affected": "2",
        "gold_v2_unaffected": "1",
        "gold_v2_status": "confirmed_original_table_derived",
        "gold_v2_note": "The 25-carrier sentence is aggregate across SCN5A-positive families; E1225K count is table/pedigree-derived.",
    },
    ("SCN5A", "15671429", "D1275N"): {
        "gold_v2_carriers": "22",
        "gold_v2_affected": "7",
        "gold_v2_unaffected": "15",
        "gold_v2_status": "confirmed_original_pedigree_derived",
        "gold_v2_note": "Pedigree-derived mutation-carrier/phenotype count; sentence-level 22+1 haplotype count is too coarse.",
    },
    ("SCN5A", "20470418", "S1103Y"): {
        "gold_v2_carriers": "26",
        "gold_v2_affected": "17",
        "gold_v2_unaffected": "9",
        "gold_v2_status": "adjudicated_variant_carrier_count",
        "gold_v2_note": "Original 85/39/46 is the sampled cohort; only 26 were heterozygous carriers, split 17 SIDS and 9 controls.",
    },
}


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def _write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _variant_key(row: dict[str, Any], *, clinical_counts: bool) -> str:
    return str(
        row.get("variant_normalized" if clinical_counts else "variant") or ""
    ).strip()


def update_csv(path: Path, *, gene: str, clinical_counts: bool) -> int:
    fieldnames, rows = _read_rows(path)
    for column in V2_COLUMNS:
        if column not in fieldnames:
            fieldnames.append(column)

    updated = 0
    for row in rows:
        pmid = str(row.get("pmid") or "").strip()
        variant = _variant_key(row, clinical_counts=clinical_counts)
        decision = ADJUDICATIONS.get((gene, pmid, variant))
        if not decision:
            continue
        for column in V2_COLUMNS:
            row[column] = ""
        row.update(decision)
        row["gold_v2_source"] = ADJUDICATION_SOURCE
        updated += 1

    _write_rows(path, fieldnames, rows)
    return updated


def main() -> int:
    root = repo_path("gene_variant_fetcher_gold_standard/normalized")
    total = 0
    for gene in ("KCNH2", "RYR2", "SCN5A"):
        total += update_csv(
            root / f"{gene}_recall_input.csv",
            gene=gene,
            clinical_counts=False,
        )
        total += update_csv(
            root / f"{gene}_clinical_counts.csv",
            gene=gene,
            clinical_counts=True,
        )
    print(f"updated_rows={total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
