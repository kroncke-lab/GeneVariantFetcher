#!/usr/bin/env python3
"""Apply the 2026-09-08 SCN5A gold adjudications through the gold_v2 columns.

Approved by Brett Kroncke on 2026-09-08 ("I think your proposals are good") on
the proposal in ../README.md section 1. Original values are never overwritten:
the repository's gold_v2_{carriers,affected,unaffected,status,note,source}
columns carry the adjudicated values and the scorers already prefer them when
gold_v2_status is populated (utils/gold_standard.py). Exactly one status,
excluded_duplicate_current_cohort, removes a row from every score.

Two decisions:

* SCN5A 20129283 (Kapplinger 2010): gold stores one row per testing centre
  listed in Table 4. The paper prints one pooled count per nucleotide change and
  no per-centre split. Each multi-row variant keeps its first row with
  carriers = affected = the printed count and unaffected = 0; the other rows are
  excluded as duplicates of the same cohort. Distinct nucleotide changes with the
  same protein effect stay separate. K1493X (gold 2) takes the printed 1.
* SCN5A 30059973 (Baruteau 2018): gold records every carrier as affected. Table
  14 splits the twelve commonest variants by phenotype; the ECG-negative column
  becomes unaffected and the rest affected. The remaining 173 variants have no
  per-variant phenotype in the paper: carriers stay, affected/unaffected become
  explicit nulls.

Usage: apply_adjudications.py [--apply]   (default is a dry run that only prints)
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
GOLD = REPO / "gene_variant_fetcher_gold_standard/normalized/SCN5A_recall_input.csv"
TABLE4 = REPO / "corpus/SCN5A/20129283/20129283_FULL_CONTEXT.md"
TABLE14 = Path(__file__).resolve().parent / "gold_30059973_table14_vs_gold.csv"
SOURCE = "docs/evidence/phenotype_supply_20260908/item1_gold_adjudication/ADJUDICATIONS_20260908.md"
APPROVAL = "approved by Brett Kroncke 2026-09-08"


def norm(value: str) -> str:
    """Join gold and table spellings: drop p./c., the novelty asterisk the
    paper appends to new variants, and spacing; write stop codons as X."""
    value = value.strip().rstrip("*").upper()
    value = re.sub(r"^P\.", "", value)
    value = re.sub(r"^C\.", "", value)
    value = value.replace("*", "X")
    return re.sub(r"[^0-9A-Z+>_]", "", value)


def table4_rows() -> dict[str, list[dict]]:
    lines = TABLE4.read_text().split("\n")
    start = next(
        i
        for i, line in enumerate(lines)
        if line.startswith("| Region | Nucleotide change | Coding effect")
    )
    by_key: dict[str, list[dict]] = defaultdict(list)
    for line in lines[start + 2 :]:
        if line.startswith("###"):
            break
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 7 or not cells[5].isdigit():
            continue
        row = {
            "cdna": cells[1],
            "protein": cells[2].rstrip("*"),
            "n": int(cells[5]),
            "centres": cells[6],
        }
        keys = {norm(row["protein"])} if row["protein"] else set()
        keys.add(norm(row["cdna"]))
        for key in keys:
            if key:
                by_key[key].append(row)
    return by_key


def plan_20129283(
    rows: list[dict], t4: dict[str, list[dict]]
) -> list[tuple[dict, dict]]:
    """Return (row, update) pairs for PMID 20129283."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["pmid"] == "20129283":
            groups[row["variant"]].append(row)
    updates: list[tuple[dict, dict]] = []
    for variant, group in groups.items():
        table = t4.get(norm(variant), [])
        # dedupe table rows (a protein key and its cDNA key point at the same row)
        seen = set()
        table = [r for r in table if not (id(r) in seen or seen.add(id(r)))]
        if not table:
            if len(group) > 1:
                print(
                    f"  ! {variant}: {len(group)} gold rows, no Table 4 row found; left unchanged"
                )
            continue
        total = sum(r["n"] for r in table)
        printed = "; ".join(
            f"{r['cdna']} {r['protein']}".strip()
            + f" = {r['n']} (centres {r['centres']})"
            for r in table
        )
        if len(group) == 1:
            row = group[0]
            if int(row["carriers"]) != total:
                if variant == "K1493X":
                    updates.append(
                        (
                            row,
                            {
                                "gold_v2_carriers": str(total),
                                "gold_v2_affected": str(total),
                                "gold_v2_unaffected": "0",
                                "gold_v2_status": "adjudicated_variant_carrier_count",
                                "gold_v2_note": f"Table 4 prints {printed}; gold 2 conflated K1493del (a separate row). {APPROVAL}.",
                                "gold_v2_source": SOURCE,
                            },
                        )
                    )
                else:
                    print(
                        f"  ? {variant}: single gold row {row['carriers']} vs Table 4 {total}; reported, not changed"
                    )
            continue
        if len(group) == len(table) and sorted(
            int(r["carriers"]) for r in group
        ) == sorted(r["n"] for r in table):
            print(
                f"  = {variant}: {len(group)} gold rows match {len(table)} distinct nucleotide rows; kept"
            )
            continue
        keep, *dupes = group
        updates.append(
            (
                keep,
                {
                    "gold_v2_carriers": str(total),
                    "gold_v2_affected": str(total),
                    "gold_v2_unaffected": "0",
                    "gold_v2_status": "adjudicated_variant_carrier_count",
                    "gold_v2_note": (
                        f"Table 4 prints {printed}; the paper reports no per-centre split, so the "
                        f"{len(group)} one-row-per-centre gold rows are consolidated here. {APPROVAL}."
                    ),
                    "gold_v2_source": SOURCE,
                },
            )
        )
        for dupe in dupes:
            updates.append(
                (
                    dupe,
                    {
                        "gold_v2_carriers": "",
                        "gold_v2_affected": "",
                        "gold_v2_unaffected": "",
                        "gold_v2_status": "excluded_duplicate_current_cohort",
                        "gold_v2_note": f"Per-testing-centre duplicate of the consolidated {variant} Table 4 row. {APPROVAL}.",
                        "gold_v2_source": SOURCE,
                    },
                )
            )
    return updates


def plan_30059973(rows: list[dict]) -> list[tuple[dict, dict]]:
    split = {}
    with TABLE14.open(newline="") as handle:
        for r in csv.DictReader(handle):
            split[norm(r["variant_gold"])] = r
    updates = []
    for row in rows:
        if row["pmid"] != "30059973":
            continue
        key = norm(row["variant"])
        r = split.get(key)
        if r is not None:
            total = int(r["row_total"])
            neg = int(r["negative_ecg"])
            if int(row["carriers"]) != total:
                print(
                    f"  ? 30059973 {row['variant']}: gold carriers {row['carriers']} vs Table 14 total {total}"
                )
            updates.append(
                (
                    row,
                    {
                        "gold_v2_carriers": str(total),
                        "gold_v2_affected": str(total - neg),
                        "gold_v2_unaffected": str(neg),
                        "gold_v2_status": "adjudicated_source_phenotype_partition",
                        "gold_v2_note": (
                            f"Table 14 (source line {r['source_line']}): ECG-negative {neg}, isolated LQT3 {r['isolated_lqt3']}, "
                            f"BrS {r['isolated_brs1']}, PCCD {r['isolated_pccd']}, SSS {r['isolated_sss']}, DCM {r['isolated_dcm']}, "
                            f"overlap {r['overlap']}; affected = every phenotype-positive carrier. {APPROVAL}."
                        ),
                        "gold_v2_source": SOURCE,
                    },
                )
            )
        else:
            updates.append(
                (
                    row,
                    {
                        "gold_v2_carriers": row["carriers"],
                        "gold_v2_affected": "",
                        "gold_v2_unaffected": "",
                        "gold_v2_status": "adjudicated_phenotype_not_reported",
                        "gold_v2_note": (
                            "Table 5 lists per-variant occurrences without phenotype; Table 1 reports 196 of 442 "
                            "carriers with a negative ECG phenotype, so affected = carriers is not supported. "
                            f"Affected/unaffected are explicit nulls. {APPROVAL}."
                        ),
                        "gold_v2_source": SOURCE,
                    },
                )
            )
    return updates


def main(apply: bool) -> None:
    with GOLD.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)
    already = [
        r
        for r in rows
        if r["pmid"] in ("20129283", "30059973") and r.get("gold_v2_status")
    ]
    if already:
        print(
            f"{len(already)} rows already carry a gold_v2_status for these PMIDs; refusing to re-apply"
        )
        return
    t4 = table4_rows()
    plan = plan_20129283(rows, t4) + plan_30059973(rows)
    from collections import Counter

    print(
        "planned updates:",
        len(plan),
        dict(Counter(u["gold_v2_status"] for _, u in plan)),
    )
    for row, update in plan[:6]:
        print(
            "  e.g.",
            row["variant"],
            row["pmid"],
            f"{row['carriers']}/{row['affected']}/{row['unaffected']}",
            "->",
            update["gold_v2_status"],
            f"{update['gold_v2_carriers']}/{update['gold_v2_affected']}/{update['gold_v2_unaffected']}",
        )
    if not apply:
        print("dry run; pass --apply to write")
        return
    for row, update in plan:
        row.update(update)
    with GOLD.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {GOLD} ({len(rows)} rows, {len(plan)} adjudicated)")


if __name__ == "__main__":
    main(apply="--apply" in sys.argv[1:])
