#!/usr/bin/env python3
"""Build the standalone R01 phenotype-count recovery figure and source table.

Every identity-matched variant-paper row with an authoritative gold count enters
the analysis.  A missing automated count is scored as zero for evaluation while
remaining explicitly marked as missing in the source table; the script never
changes the underlying prediction artifact or database NULL semantics.

No plotting package is required. The script writes an editable SVG, a tidy CSV,
and a JSON summary directly from the locked evaluation artifacts. Convert the
SVG to PDF/PNG with a vector-aware renderer after generation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from html import escape
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
RUN_EVAL_PATH = REPO / "benchmarks" / "codex_paper_eval" / "run_eval.py"
DEFAULT_RUN = (
    REPO
    / "benchmarks"
    / "codex_paper_eval"
    / "runs"
    / "20260824_aha_table_sourcebound_gold118"
)
DEFAULT_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
DEFAULT_SVG = REPO / "docs" / "figures" / "r01_phenotype_count_recovery.svg"
DEFAULT_PNG = REPO / "docs" / "figures" / "r01_phenotype_count_recovery.png"
DEFAULT_PDF = REPO / "docs" / "figures" / "r01_phenotype_count_recovery.pdf"
DEFAULT_CSV = REPO / "docs" / "figures" / "data" / "r01_phenotype_count_recovery.csv"
DEFAULT_JSON = REPO / "docs" / "figures" / "data" / "r01_phenotype_count_recovery.json"

FIELDS = ("affected", "unaffected")
COLORS = {
    "navy": "#17324D",
    "blue": "#3B78A8",
    "blue_pale": "#DDEAF3",
    "teal": "#2A9D8F",
    "teal_pale": "#DDF1EE",
    "coral": "#C85A63",
    "coral_pale": "#F7E3E6",
    "slate": "#687584",
    "pale": "#D9DFE4",
    "divider": "#DCE3E8",
    "off_white": "#F8FAFC",
    "white": "#FFFFFF",
}


def load_run_eval():
    """Import the benchmark scorer without making it an installed package."""
    spec = importlib.util.spec_from_file_location("gvf_run_eval", RUN_EVAL_PATH)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise SystemExit(f"cannot import {RUN_EVAL_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_lock(run_dir: Path) -> dict[str, Any]:
    """Fail if the selection or predictions changed after the evaluation lock."""
    lock = json.loads((run_dir / "LOCK.json").read_text())
    checked = {
        "selection_sha256": sha256(run_dir / "selection.json"),
        "predictions_sha256": sha256(run_dir / "predictions.json"),
    }
    for key, observed in checked.items():
        if observed != lock.get(key):
            raise SystemExit(f"lock mismatch for {key}: {observed} != {lock.get(key)}")
    return {**lock, "verified_hashes": checked}


def gold_rows_with_provenance(
    run_eval, gold_root: Path, gene: str, pmid: str
) -> list[dict[str, Any]]:
    """Load the scorer's authoritative counts plus their adjudication source."""
    path = run_eval.gold_csv_path(gold_root, gene)
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for source in csv.DictReader(handle):
            if str(source.get("pmid", "")).strip() != pmid:
                continue
            if run_eval.gold_row_excluded(source):
                continue
            v2 = bool(str(source.get("gold_v2_status") or "").strip())
            rows.append(
                {
                    "variant": str(source.get("variant", "")).strip(),
                    **{
                        field: run_eval.authoritative_gold_count(
                            source, field, parser=run_eval.to_int
                        )
                        for field in run_eval.COUNT_FIELDS
                    },
                    "gold_value_set": "gold_v2" if v2 else "baseline",
                    "gold_v2_status": str(source.get("gold_v2_status") or ""),
                    "gold_source": str(source.get("gold_v2_source") or path),
                }
            )
    return rows


def paired_rows(
    run_eval,
    selection: dict[str, Any],
    predictions: dict[str, Any],
    gold_root: Path,
    run_id: str,
    locked_at: str,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, int]]]:
    """Reproduce identity matching and score missing matched counts as zero."""
    predicted_by_paper = {
        (paper["gene"], str(paper["pmid"])): paper for paper in predictions["papers"]
    }
    rows: list[dict[str, Any]] = []
    missingness = {
        field: {
            "matched": 0,
            "count_extracted": 0,
            "zero_filled": 0,
            "identity_miss": 0,
        }
        for field in FIELDS
    }
    for selected in selection["papers"]:
        gene, pmid = selected["gene"], str(selected["pmid"])
        predicted = predicted_by_paper.get((gene, pmid))
        if predicted is None:
            raise SystemExit(f"prediction missing for {gene} PMID {pmid}")
        gold = gold_rows_with_provenance(run_eval, gold_root, gene, pmid)
        merged, _merged_twins = run_eval.merge_notation_twins(
            predicted.get("variants", []), gene
        )
        used: set[int] = set()
        pairs: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for candidate in merged:
            hit = next(
                (
                    index
                    for index, manual in enumerate(gold)
                    if index not in used
                    and run_eval.matches(candidate["variant"], manual["variant"], gene)
                ),
                None,
            )
            if hit is None:
                continue
            used.add(hit)
            manual = gold[hit]
            pairs.append((candidate, manual))
        for candidate, manual in pairs:
            for field in FIELDS:
                ai_count = candidate.get(field)
                manual_count = manual.get(field)
                if manual_count is None:
                    continue
                count_extracted = ai_count is not None
                evaluation_count = int(ai_count) if count_extracted else 0
                missingness[field]["matched"] += 1
                if count_extracted:
                    missingness[field]["count_extracted"] += 1
                else:
                    missingness[field]["zero_filled"] += 1
                residual = evaluation_count - int(manual_count)
                rows.append(
                    {
                        "analysis_run": run_id,
                        "analysis_locked_at": locked_at,
                        "grain": "matched_variant_paper_measure",
                        "gene": gene,
                        "pmid": pmid,
                        "predicted_variant": candidate["variant"],
                        "manual_variant": manual["variant"],
                        "measure": field,
                        "manual_count": int(manual_count),
                        "ai_count_raw": int(ai_count) if count_extracted else None,
                        "ai_count_evaluation": evaluation_count,
                        "evaluation_zero_filled": int(not count_extracted),
                        "residual_ai_minus_manual": residual,
                        "absolute_difference": abs(residual),
                        "exact_match": int(residual == 0),
                        "zero_nonzero_disagreement": int(
                            (manual_count == 0) != (evaluation_count == 0)
                        ),
                        "prediction_source_layer": candidate.get("source_layer") or "",
                        "prediction_source_location": candidate.get("source_location")
                        or "",
                        "gold_value_set": manual["gold_value_set"],
                        "gold_v2_status": manual["gold_v2_status"],
                        "gold_source": manual["gold_source"],
                    }
                )
        for index, manual in enumerate(gold):
            if index in used:
                continue
            for field in FIELDS:
                if manual.get(field) is not None:
                    missingness[field]["identity_miss"] += 1
    rows.sort(
        key=lambda row: (
            row["measure"],
            row["gene"],
            int(row["pmid"]),
            row["manual_variant"],
        )
    )
    keys = [
        (row["measure"], row["gene"], row["pmid"], row["manual_variant"])
        for row in rows
    ]
    if len(keys) != len(set(keys)):
        raise SystemExit("duplicate plotted grain detected")
    return rows, missingness


def summarize(
    rows: list[dict[str, Any]],
    missingness: dict[str, dict[str, int]],
    report: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    summary: dict[str, dict[str, Any]] = {}

    def fraction(numerator, denominator):
        return numerator / denominator if denominator else None

    for field in FIELDS:
        values = [row for row in rows if row["measure"] == field]
        absolute = [int(row["absolute_difference"]) for row in values]
        asserted = int(report["overall"]["count"][field]["gold_asserted"])
        matched = len(values)
        gold_nonzero = [row for row in values if int(row["manual_count"]) != 0]
        zero_baseline_exact = sum(int(row["manual_count"]) == 0 for row in values)
        nonzero_exact = sum(int(row["exact_match"]) for row in gold_nonzero)
        record = {
            "gold_asserted": asserted,
            "identity_matched": matched,
            "identity_coverage": fraction(matched, asserted),
            "count_extracted": missingness[field]["count_extracted"],
            "count_extracted_fraction": fraction(
                missingness[field]["count_extracted"], matched
            ),
            "zero_filled": missingness[field]["zero_filled"],
            "zero_filled_fraction": fraction(
                missingness[field]["zero_filled"], matched
            ),
            "exact": sum(error == 0 for error in absolute),
            "exact_fraction": fraction(sum(error == 0 for error in absolute), matched),
            "zero_baseline_exact": zero_baseline_exact,
            "zero_baseline_fraction": fraction(zero_baseline_exact, matched),
            "increment_over_zero_baseline": fraction(
                sum(error == 0 for error in absolute) - zero_baseline_exact, matched
            ),
            "gold_nonzero": len(gold_nonzero),
            "nonzero_exact": nonzero_exact,
            "nonzero_exact_fraction": fraction(nonzero_exact, len(gold_nonzero)),
            "mae": statistics.fmean(absolute) if absolute else None,
            "rmse": (
                math.sqrt(statistics.fmean(error**2 for error in absolute))
                if absolute
                else None
            ),
            "median_absolute_difference": statistics.median(absolute)
            if absolute
            else None,
            "difference_equal_one": sum(error == 1 for error in absolute),
            "within_one": sum(error <= 1 for error in absolute),
            "greater_than_one": sum(error > 1 for error in absolute),
            "identity_miss": missingness[field]["identity_miss"],
            "manual_max": max(
                (int(row["manual_count"]) for row in values), default=None
            ),
            "ai_max": max(
                (int(row["ai_count_evaluation"]) for row in values), default=None
            ),
        }
        if record["count_extracted"] != int(
            report["overall"]["count"][field]["predicted"]
        ):
            raise SystemExit(f"{field} extracted-count total does not match report")
        if record["identity_matched"] + record["identity_miss"] != asserted:
            raise SystemExit(f"{field} missingness does not sum to gold assertions")
        if record["count_extracted"] + record["zero_filled"] != matched:
            raise SystemExit(f"{field} evaluation rows do not sum to identity matches")
        summary[field] = record
    return summary


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        fields = (
            list(rows[0])
            if rows
            else [
                "analysis_run",
                "analysis_locked_at",
                "grain",
                "gene",
                "pmid",
                "predicted_variant",
                "manual_variant",
                "measure",
                "manual_count",
                "ai_count_raw",
                "ai_count_evaluation",
                "evaluation_zero_filled",
                "residual_ai_minus_manual",
                "absolute_difference",
                "exact_match",
                "zero_nonzero_disagreement",
                "prediction_source_layer",
                "prediction_source_location",
                "gold_value_set",
                "gold_v2_status",
                "gold_source",
            ]
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fmt_pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def text_element(
    x: float,
    y: float,
    value: str,
    class_name: str,
    *,
    anchor: str | None = None,
    transform: str | None = None,
) -> str:
    attributes = [f'x="{x:.1f}"', f'y="{y:.1f}"', f'class="{class_name}"']
    if anchor:
        attributes.append(f'text-anchor="{anchor}"')
    if transform:
        attributes.append(f'transform="{escape(transform)}"')
    return f"<text {' '.join(attributes)}>{escape(value)}</text>"


PLOT_MAX = 110
PLOT_INSET = 46.0


def plot_coordinate(value: float, start: float, size: float) -> float:
    """Map a nonnegative count to a shared square-root plotting scale."""
    usable = size - 2 * PLOT_INSET
    return start + PLOT_INSET + math.sqrt(value / PLOT_MAX) * usable


def plot_y(value: float, top: float, size: float) -> float:
    usable = size - 2 * PLOT_INSET
    return top + size - PLOT_INSET - math.sqrt(value / PLOT_MAX) * usable


def diamond_points(cx: float, cy: float, radius: float) -> str:
    return " ".join(
        f"{x:.1f},{y:.1f}"
        for x, y in (
            (cx, cy - radius),
            (cx + radius, cy),
            (cx, cy + radius),
            (cx - radius, cy),
        )
    )


def grouped_coordinates(rows: list[dict[str, Any]], field: str) -> list[dict[str, Any]]:
    groups: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["measure"] == field:
            groups[(row["manual_count"], row["ai_count_evaluation"])].append(row)
    result = []
    for (manual_count, ai_count), members in sorted(groups.items()):
        result.append(
            {
                "manual_count": manual_count,
                "ai_count": ai_count,
                "members": members,
                "n": len(members),
                "zero_filled_n": sum(
                    int(member["evaluation_zero_filled"]) for member in members
                ),
                "exact": manual_count == ai_count,
                "zero_nonzero": (manual_count == 0) != (ai_count == 0),
            }
        )
    return result


def draw_bubble(
    group: dict[str, Any], field: str, plot_x: float, plot_top: float, size: float
) -> list[str]:
    cx = plot_coordinate(group["manual_count"], plot_x, size)
    cy = plot_y(group["ai_count"], plot_top, size)
    radius = bubble_radius(group["n"])
    color = COLORS["blue"] if field == "affected" else COLORS["teal"]
    pale = COLORS["blue_pale"] if field == "affected" else COLORS["teal_pale"]
    if group["exact"]:
        if group["zero_filled_n"]:
            fill, stroke, stroke_width = pale, color, 4.0
        else:
            fill, stroke, stroke_width = color, COLORS["navy"], 2.5
    else:
        fill = COLORS["coral_pale"] if group["zero_nonzero"] else COLORS["white"]
        stroke, stroke_width = COLORS["coral"], 3.0
    member_names = "; ".join(
        f"{row['gene']} {row['manual_variant']} PMID {row['pmid']}"
        for row in group["members"]
    )
    title = escape(
        f"{field}: manual {group['manual_count']}, evaluated AI {group['ai_count']}; "
        f"{group['n']} row(s), {group['zero_filled_n']} zero-filled: {member_names}"
    )
    opacity = 0.82 if group["exact"] else 0.96
    if field == "affected":
        shape = (
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{radius:.1f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" '
            f'fill-opacity="{opacity}"><title>{title}</title></circle>'
        )
    else:
        shape = (
            f'<polygon points="{diamond_points(cx, cy, radius)}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" '
            f'fill-opacity="{opacity}"><title>{title}</title></polygon>'
        )
    return [shape]


def bubble_radius(pair_count: int) -> float:
    """Return the marker radius used for one grouped plot coordinate."""
    if pair_count < 1:
        raise ValueError("pair_count must be positive")
    return min(28.0, 8.0 + 2.0 * math.sqrt(pair_count))


EXAMPLE_ANNOTATIONS = {
    ("affected", "KCNH2", "22338672", "K897T"): (74.0, -390.0),
    ("affected", "RYR2", "15466642", "I4848V"): (135.0, -130.0),
    ("affected", "RYR2", "25814417", "G357S"): (-195.0, 75.0),
    ("affected", "KCNQ1", "14678125", "V254M"): (-140.0, -284.0),
    ("affected", "SCN5A", "28339995", "D1790G"): (-118.0, -84.0),
    ("unaffected", "KCNH2", "25819988", "H562R"): (65.0, -180.0),
    ("unaffected", "RYR2", "25814417", "G357S"): (-279.0, -117.0),
    ("unaffected", "RYR2", "33315912", "P2328S"): (-274.0, -245.0),
    ("unaffected", "KCNQ1", "18808722", "L187P"): (-83.0, -129.0),
    ("unaffected", "KCNH2", "22338672", "R744X"): (-149.0, -294.0),
}

EXAMPLE_DOCKS = {
    ("affected", "SCN5A", "28339995", "D1790G"): "right",
    ("unaffected", "RYR2", "25814417", "G357S"): "right",
    ("unaffected", "RYR2", "33315912", "P2328S"): "right",
    ("unaffected", "KCNH2", "22338672", "R744X"): "right",
}


def annotation_connector_anchor(
    *,
    point_x: float,
    point_y: float,
    label_x: float,
    label_y: float,
    gene: str,
    variant: str,
    pmid: str,
    dock: str = "left",
) -> tuple[float, float]:
    """Return a designed dock on the facing side of a two-line label."""
    title = f"{gene} {variant}"
    title_width = len(title) * 30.0 * 0.54
    docks = {
        "left": (label_x - 12.0, label_y + 17.0),
        "right": (label_x + title_width + 8.0, label_y + 17.0),
        "top_left": (label_x - 12.0, label_y - 31.0),
        "under": (label_x + 0.55 * title_width, label_y + 50.0),
    }
    if dock not in docks:
        raise ValueError(f"unknown annotation dock: {dock}")
    return docks[dock]


def draw_example_annotations(
    rows: list[dict[str, Any]],
    field: str,
    plot_x: float,
    plot_top: float,
    size: float,
    *,
    annotations: dict[tuple[str, str, str, str], tuple[float, float]] | None = None,
    docks: dict[tuple[str, str, str, str], str] | None = None,
) -> list[str]:
    """Label a few source rows without turning the plot into a lookup table."""
    if annotations is None:
        annotations = EXAMPLE_ANNOTATIONS
    if docks is None:
        docks = EXAMPLE_DOCKS
    output: list[str] = []
    for (measure, gene, pmid, variant), (dx, dy) in annotations.items():
        if measure != field:
            continue
        matches = [
            row
            for row in rows
            if row["measure"] == measure
            and row["gene"] == gene
            and row["pmid"] == pmid
            and row["manual_variant"] == variant
        ]
        if not matches:
            # A later run may miss one of the illustrative identities. The
            # figure must still build and honestly omit that unavailable label.
            continue
        if len(matches) != 1:
            raise SystemExit(
                f"expected one annotation row for {measure} {gene} {variant} {pmid}"
            )
        row = matches[0]
        point_x = plot_coordinate(row["manual_count"], plot_x, size)
        point_y = plot_y(row["ai_count_evaluation"], plot_top, size)
        label_x, label_y = point_x + dx, point_y + dy
        anchor_x, anchor_y = annotation_connector_anchor(
            point_x=point_x,
            point_y=point_y,
            label_x=label_x,
            label_y=label_y,
            gene=gene,
            variant=variant,
            pmid=pmid,
            dock=docks.get((measure, gene, pmid, variant), "left"),
        )
        output.extend(
            [
                (
                    f'<path d="M{point_x:.1f},{point_y:.1f} '
                    f'L{anchor_x:.1f},{anchor_y:.1f}" '
                    'class="example-line"/>'
                ),
                text_element(
                    label_x,
                    label_y,
                    f"{gene} {variant}",
                    "example-title",
                ),
                text_element(
                    label_x,
                    label_y + 34,
                    f"PMID {pmid}",
                    "example-pmid",
                ),
            ]
        )
    return output


MARKER_SIZE_KEY_COUNTS = (5, 10, 25)


def draw_marker_size_key(plot_x: float, plot_top: float) -> list[str]:
    """Draw a compact marker-size key in the upper-left plot whitespace."""
    box_x, box_y = plot_x + 28.0, plot_top + 28.0
    box_width, box_height = 372.0, 174.0
    centers = (box_x + 72.0, box_x + 184.0, box_x + 306.0)
    marker_y = box_y + 89.0
    output = [
        '<g role="group" aria-label="Marker size key for variant-paper pairs">',
        (
            f'<rect x="{box_x:.1f}" y="{box_y:.1f}" width="{box_width:.1f}" '
            f'height="{box_height:.1f}" rx="12" class="size-key-box"/>'
        ),
        text_element(
            box_x + 18.0,
            box_y + 38.0,
            "Marker size · variant–paper pairs",
            "size-key-title",
        ),
    ]
    for cx, pair_count in zip(centers, MARKER_SIZE_KEY_COUNTS, strict=True):
        output.extend(
            [
                (
                    f'<circle cx="{cx:.1f}" cy="{marker_y:.1f}" '
                    f'r="{bubble_radius(pair_count):.1f}" class="size-key-marker"/>'
                ),
                text_element(
                    cx,
                    box_y + 148.0,
                    str(pair_count),
                    "size-key-value",
                    anchor="middle",
                ),
            ]
        )
    output.append("</g>")
    return output


CALLOUTS = {
    ("affected", "KCNH2", "22338672", 0, 7): (278, 556),
    ("affected", "KCNH2", "25819988", 9, 6): (620, 648),
    ("affected", "RYR2", "15466642", 4, 6): (265, 753),
    ("unaffected", "KCNH2", "25819988", 4, 7): (1485, 558),
    ("unaffected", "KCNH2", "20181576", 3, 1): (1450, 823),
}


def draw_callouts(
    rows: list[dict[str, Any]], field: str, plot_x: float, plot_top: float, size: float
) -> list[str]:
    output: list[str] = []
    for row in rows:
        key = (
            field,
            row["gene"],
            row["pmid"],
            row["manual_count"],
            row["ai_count"],
        )
        if key not in CALLOUTS:
            continue
        label_x, label_y = CALLOUTS[key]
        point_x = plot_coordinate(row["manual_count"], plot_x, size)
        point_y = plot_y(row["ai_count"], plot_top, size)
        box_width = 315
        output.extend(
            [
                (
                    f'<path d="M{point_x:.1f},{point_y:.1f} '
                    f'L{label_x - 8:.1f},{label_y + 22:.1f}" '
                    f'class="callout-line"/>'
                ),
                (
                    f'<rect x="{label_x - 12:.1f}" y="{label_y - 28:.1f}" '
                    f'width="{box_width}" height="82" rx="8" '
                    f'class="callout-box"/>'
                ),
                text_element(
                    label_x,
                    label_y,
                    f"{row['manual_variant']} · PMID {row['pmid']}",
                    "callout-title",
                ),
                text_element(
                    label_x,
                    label_y + 38,
                    f"manual {row['manual_count']}; AI {row['ai_count']}",
                    "callout-value",
                ),
            ]
        )
    return output


def render_svg(
    rows: list[dict[str, Any]],
    summary: dict[str, dict[str, Any]],
    paired_overlap: int,
    run_id: str,
    locked_at: str,
) -> str:
    width, height = 2160, 1260
    plot_top, plot_size = 438.0, 500.0
    plot_x = {"affected": 210.0, "unaffected": 1240.0}
    facet_color = {"affected": COLORS["blue"], "unaffected": COLORS["teal"]}
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" '
            f'height="4.2in" viewBox="0 0 {width} {height}" role="img" '
            'aria-labelledby="figure-title figure-desc">'
        ),
        '<title id="figure-title">Carrier-count agreement among paired '
        "variant-paper rows</title>",
        '<desc id="figure-desc">Paired-count coverage is 56 of 632 '
        "affected rows and 28 of 631 unaffected rows. Discrete bubble plots "
        "compare frozen manual and AI-extracted counts on matched variant-paper "
        "rows where both values are numeric.</desc>",
        "<metadata>",
        f"Analysis run: {escape(run_id)}",
        f"Locked at: {escape(locked_at)}",
        "Grain: matched variant-paper row and count measure; not an explicit cohort",
        f"Rows present in both affected and unaffected facets: {paired_overlap}",
        "Bubble area is proportional to the number of rows at a coordinate.",
        "</metadata>",
        "<defs>",
        "<style><![CDATA[",
        (
            'text { font-family: "Avenir Next Condensed", "Source Sans 3", '
            '"Helvetica Neue", Arial, sans-serif; fill: #17324D; }'
        ),
        ".title { font-size: 72px; font-weight: 600; letter-spacing: -0.4px; }",
        ".coverage-title { font-size: 39px; font-weight: 600; }",
        ".coverage-label { font-size: 38px; font-weight: 500; }",
        ".coverage-value { font-size: 38px; font-weight: 600; }",
        ".coverage-detail { font-size: 36px; font-weight: 400; fill: #687584; }",
        ".facet-title { font-size: 48px; font-weight: 600; }",
        ".facet-note { font-size: 36px; font-weight: 400; fill: #687584; }",
        ".tick { font-size: 36px; font-weight: 400; fill: #687584; }",
        ".axis-label { font-size: 38px; font-weight: 500; }",
        ".bubble-count { font-size: 34px; font-weight: 600; }",
        ".metric { font-size: 38px; font-weight: 500; }",
        ".metric-muted { font-size: 36px; font-weight: 400; fill: #687584; }",
        ".callout-title { font-size: 36px; font-weight: 600; }",
        ".callout-value { font-size: 36px; font-weight: 400; fill: #C85A63; }",
        ".grid { stroke: #DCE3E8; stroke-width: 2.5; }",
        ".axis { stroke: #17324D; stroke-width: 4; }",
        ".identity { stroke: #17324D; stroke-width: 5; stroke-dasharray: 13 10; }",
        ".card { fill: #FFFFFF; stroke: #DCE3E8; stroke-width: 3; }",
        ".callout-box { fill: #F8FAFC; fill-opacity: 0.96; }",
        ".callout-line { stroke: #C85A63; stroke-width: 3.5; fill: none; }",
        "]]></style>",
        (
            f'<clipPath id="affected-clip"><rect x="{plot_x["affected"]}" '
            f'y="{plot_top}" width="{plot_size}" height="{plot_size}"/></clipPath>'
        ),
        (
            f'<clipPath id="unaffected-clip"><rect x="{plot_x["unaffected"]}" '
            f'y="{plot_top}" width="{plot_size}" height="{plot_size}"/></clipPath>'
        ),
        "</defs>",
        f'<rect width="{width}" height="{height}" fill="{COLORS["off_white"]}"/>',
        text_element(
            80,
            94,
            "Carrier-count agreement among paired variant-paper rows",
            "title",
        ),
        '<rect x="80" y="126" width="2000" height="205" rx="18" class="card"/>',
        text_element(
            108,
            169,
            "Rows entering the agreement plots",
            "coverage-title",
        ),
        text_element(
            2018,
            169,
            "Colored segment = paired rows",
            "facet-note",
            anchor="end",
        ),
    ]

    bar_x, bar_width = 390.0, 1100.0
    coverage_y = {"affected": 213.0, "unaffected": 284.0}
    for field in FIELDS:
        stats = summary[field]
        y = coverage_y[field]
        lines.extend(
            [
                text_element(110, y + 11, field.capitalize(), "coverage-label"),
                (
                    f'<rect x="{bar_x}" y="{y - 17}" width="{bar_width}" '
                    f'height="26" rx="13" fill="{COLORS["pale"]}"/>'
                ),
                (
                    f'<rect x="{bar_x}" y="{y - 17}" '
                    f'width="{bar_width * stats["coverage"]:.1f}" height="26" '
                    f'rx="13" fill="{facet_color[field]}"/>'
                ),
                text_element(
                    1540,
                    y + 11,
                    f"{stats['paired']}/{stats['gold_asserted']} rows ({fmt_pct(stats['coverage'])})",
                    "coverage-value",
                ),
                text_element(
                    1540,
                    y + 42,
                    (
                        f"{stats['matched_ai_null']} AI null; "
                        f"{stats['identity_miss']} identity misses"
                    ),
                    "coverage-detail",
                ),
            ]
        )

    lines.extend(
        [
            text_element(115, 374, "Affected carriers", "facet-title"),
            text_element(1145, 374, "Unaffected carriers", "facet-title"),
            text_element(
                115,
                414,
                "Circles; area represents rows; numerals mark repeats",
                "facet-note",
            ),
            text_element(
                1145,
                414,
                "Diamonds; coral outline marks disagreement",
                "facet-note",
            ),
        ]
    )

    ticks = (0, 2, 4, 6, 8, 10, 12)
    for field in FIELDS:
        x0 = plot_x[field]
        for tick in ticks:
            x = plot_coordinate(tick, x0, plot_size)
            y = plot_y(tick, plot_top, plot_size)
            lines.extend(
                [
                    f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" '
                    f'y2="{plot_top + plot_size}" class="grid"/>',
                    f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + plot_size}" '
                    f'y2="{y:.1f}" class="grid"/>',
                    text_element(x, 978, str(tick), "tick", anchor="middle"),
                    text_element(x0 - 25, y + 12, str(tick), "tick", anchor="end"),
                ]
            )
        zero_x = plot_coordinate(0, x0, plot_size)
        zero_y = plot_y(0, plot_top, plot_size)
        end_x = plot_coordinate(13, x0, plot_size)
        end_y = plot_y(13, plot_top, plot_size)
        lines.extend(
            [
                f'<line x1="{x0}" y1="{plot_top + plot_size}" '
                f'x2="{x0 + plot_size}" y2="{plot_top + plot_size}" class="axis"/>',
                f'<line x1="{x0}" y1="{plot_top}" x2="{x0}" '
                f'y2="{plot_top + plot_size}" class="axis"/>',
                (
                    f'<line x1="{zero_x:.1f}" y1="{zero_y:.1f}" '
                    f'x2="{end_x:.1f}" y2="{end_y:.1f}" class="identity"/>'
                ),
                text_element(
                    x0 + plot_size / 2,
                    1020,
                    "Frozen manual count",
                    "axis-label",
                    anchor="middle",
                ),
            ]
        )
        groups = grouped_coordinates(rows, field)
        lines.append(f'<g clip-path="url(#{field}-clip)">')
        for group in groups:
            lines.extend(draw_bubble(group, field, x0, plot_top, plot_size))
        lines.append("</g>")
        lines.extend(draw_callouts(rows, field, x0, plot_top, plot_size))

    lines.append(
        text_element(
            73,
            plot_top + plot_size / 2,
            "AI-extracted count",
            "axis-label",
            anchor="middle",
            transform=f"rotate(-90 73 {plot_top + plot_size / 2})",
        )
    )

    card_y, card_height, card_width = 1040.0, 168.0, 970.0
    for field, card_x in (("affected", 80.0), ("unaffected", 1110.0)):
        stats = summary[field]
        two_column_x = (card_x + card_width * 0.25, card_x + card_width * 0.75)
        three_column_x = (
            card_x + card_width / 6,
            card_x + card_width / 2,
            card_x + 5 * card_width / 6,
        )
        lines.extend(
            [
                f'<rect x="{card_x}" y="{card_y}" width="{card_width}" '
                f'height="{card_height}" '
                'rx="15" class="card"/>',
                text_element(
                    two_column_x[0],
                    card_y + 51,
                    f"Paired {stats['paired']}/{stats['gold_asserted']} ({fmt_pct(stats['coverage'])})",
                    "metric",
                    anchor="middle",
                ),
                text_element(
                    two_column_x[1],
                    card_y + 51,
                    f"Exact {stats['exact']}/{stats['paired']} ({fmt_pct(stats['exact_fraction'])})",
                    "metric",
                    anchor="middle",
                ),
                text_element(
                    three_column_x[0],
                    card_y + 124,
                    f"Mean abs. diff. {stats['mae']:.2f}",
                    "metric-muted",
                    anchor="middle",
                ),
                text_element(
                    three_column_x[1],
                    card_y + 124,
                    f"Abs. diff. = 1: {stats['difference_equal_one']}/{stats['paired']}",
                    "metric-muted",
                    anchor="middle",
                ),
                text_element(
                    three_column_x[2],
                    card_y + 124,
                    f"> 1: {stats['greater_than_one']}/{stats['paired']}",
                    "metric-muted",
                    anchor="middle",
                ),
            ]
        )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def render_zero_filled_svg(
    rows: list[dict[str, Any]],
    summary: dict[str, dict[str, Any]],
    paired_overlap: int,
    run_id: str,
    locked_at: str,
    *,
    title: str = "Automated versus reference phenotype counts",
    example_annotations: dict[tuple[str, str, str, str], tuple[float, float]]
    | None = None,
    example_docks: dict[tuple[str, str, str, str], str] | None = None,
) -> str:
    """Render the publication figure using evaluation-time zero filling."""

    width, height = 2160, 1260
    plot_top, plot_size = 190.0, 900.0
    plot_x = {"affected": 100.0, "unaffected": 1160.0}
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" '
            f'height="4.2in" viewBox="0 0 {width} {height}" role="img" '
            'aria-labelledby="figure-title figure-desc">'
        ),
        f'<title id="figure-title">{escape(title)}</title>',
        '<desc id="figure-desc">All identity-matched variant-paper rows are '
        "scored. Missing automated counts are evaluated as zero. Agreement "
        "plots use a square-root count scale.</desc>",
        "<metadata>",
        f"Analysis run: {escape(run_id)}",
        f"Locked at: {escape(locked_at)}",
        "Evaluation rule: matched missing automated count is scored as zero",
        "Storage rule: missing automated count remains NULL in source data",
        "Identity misses are excluded from the plotted/scored rows",
        f"Rows present in both count facets: {paired_overlap}",
        "Marker radius is 8 + 2*sqrt(pair count), capped at 28 viewBox units",
        "Marker-size key values: 5, 10, and 25 variant-paper pairs",
        "</metadata>",
        "<defs>",
        "<style><![CDATA[",
        (
            'text { font-family: "Avenir Next Condensed", "Source Sans 3", '
            '"Helvetica Neue", Arial, sans-serif; fill: #17324D; }'
        ),
        ".title { font-size: 56px; font-weight: 600; letter-spacing: -0.3px; }",
        ".facet-title { font-size: 44px; font-weight: 600; }",
        ".tick { font-size: 29px; font-weight: 400; fill: #687584; }",
        ".axis-label { font-size: 32px; font-weight: 500; }",
        ".example-title { font-size: 30px; font-weight: 600; }",
        ".example-pmid { font-size: 27px; font-weight: 500; fill: #687584; }",
        ".example-line { stroke: #687584; stroke-width: 2.4; fill: none; }",
        ".size-key-box { fill: #FFFFFF; fill-opacity: 0.94; stroke: #DCE3E8; stroke-width: 2.2; }",
        ".size-key-title { font-size: 25px; font-weight: 600; }",
        ".size-key-value { font-size: 25px; font-weight: 600; fill: #687584; }",
        ".size-key-marker { fill: #FFFFFF; stroke: #17324D; stroke-width: 2.5; }",
        ".grid { stroke: #DCE3E8; stroke-width: 2.2; }",
        ".axis { stroke: #17324D; stroke-width: 3.5; fill: none; }",
        ".identity { stroke: #17324D; stroke-width: 4; stroke-dasharray: 12 10; }",
        " ]]></style>",
        (
            f'<clipPath id="affected-clip"><rect x="{plot_x["affected"]}" '
            f'y="{plot_top}" width="{plot_size}" height="{plot_size}"/></clipPath>'
        ),
        (
            f'<clipPath id="unaffected-clip"><rect x="{plot_x["unaffected"]}" '
            f'y="{plot_top}" width="{plot_size}" height="{plot_size}"/></clipPath>'
        ),
        "</defs>",
        f'<rect width="{width}" height="{height}" fill="{COLORS["off_white"]}"/>',
        text_element(
            80,
            82,
            title,
            "title",
        ),
        text_element(100, 161, "Affected individuals", "facet-title"),
        text_element(1160, 161, "Unaffected individuals", "facet-title"),
    ]

    ticks = (0, 1, 2, 5, 10, 20, 50, 100)
    for field in FIELDS:
        x0 = plot_x[field]
        for tick in ticks:
            x = plot_coordinate(tick, x0, plot_size)
            y = plot_y(tick, plot_top, plot_size)
            lines.extend(
                [
                    f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" '
                    f'y2="{plot_top + plot_size}" class="grid"/>',
                    f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + plot_size}" '
                    f'y2="{y:.1f}" class="grid"/>',
                    text_element(
                        x, plot_top + plot_size + 37, str(tick), "tick", anchor="middle"
                    ),
                ]
            )
            if field == "affected":
                lines.append(
                    text_element(x0 - 20, y + 10, str(tick), "tick", anchor="end")
                )
        zero_x = plot_coordinate(0, x0, plot_size)
        zero_y = plot_y(0, plot_top, plot_size)
        end_x = plot_coordinate(PLOT_MAX, x0, plot_size)
        end_y = plot_y(PLOT_MAX, plot_top, plot_size)
        lines.extend(
            [
                f'<rect x="{x0}" y="{plot_top}" width="{plot_size}" '
                f'height="{plot_size}" class="axis"/>',
                (
                    f'<line x1="{zero_x:.1f}" y1="{zero_y:.1f}" '
                    f'x2="{end_x:.1f}" y2="{end_y:.1f}" class="identity"/>'
                ),
                text_element(
                    x0 + plot_size / 2,
                    1188,
                    "Manual reference count",
                    "axis-label",
                    anchor="middle",
                ),
            ]
        )
        groups = grouped_coordinates(rows, field)
        if not groups:
            lines.append(
                text_element(
                    x0 + plot_size / 2,
                    plot_top + plot_size / 2,
                    "No identity-matched count rows",
                    "facet-title",
                    anchor="middle",
                )
            )
            lines.append(
                text_element(
                    x0 + plot_size / 2,
                    plot_top + plot_size / 2 + 48,
                    "Agreement is undefined; see the all-gold difference figure.",
                    "tick",
                    anchor="middle",
                )
            )
        lines.append(f'<g clip-path="url(#{field}-clip)">')
        for group in sorted(groups, key=lambda item: item["n"], reverse=True):
            lines.extend(draw_bubble(group, field, x0, plot_top, plot_size))
        lines.append("</g>")
        lines.extend(
            draw_example_annotations(
                rows,
                field,
                x0,
                plot_top,
                plot_size,
                annotations=example_annotations,
                docks=example_docks,
            )
        )

    # One key applies to circles and diamonds because both use the same radius
    # function. It is placed over otherwise-empty upper-left plot space.
    lines.extend(draw_marker_size_key(plot_x["affected"], plot_top))

    lines.append(
        text_element(
            32,
            plot_top + plot_size / 2,
            "Automated count (missing scored 0)",
            "axis-label",
            anchor="middle",
            transform=f"rotate(-90 32 {plot_top + plot_size / 2})",
        )
    )

    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def export_svg_copies(
    svg_path: Path,
    png_path: Path,
    pdf_path: Path,
    *,
    width: int = 2160,
    height: int = 1260,
) -> list[Path]:
    """Best-effort vector/raster copies without adding a Python dependency.

    librsvg is preferred when installed. macOS ``sips`` is the built-in
    fallback used by the workstation that runs the production evaluations.
    The editable SVG remains the required, portable figure artifact.
    """

    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run(
            [
                rsvg,
                "--format=png",
                f"--width={width}",
                "--output",
                str(png_path),
                str(svg_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [rsvg, "--format=pdf", "--output", str(pdf_path), str(svg_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return [png_path, pdf_path]

    qlmanage = shutil.which("qlmanage")
    ffmpeg = shutil.which("ffmpeg")
    sips = shutil.which("sips")
    if qlmanage and ffmpeg and sips:
        # Quick Look rasterizes the SVG at the requested native width instead
        # of upsampling a 72-dpi preview. It emits a square canvas with the SVG
        # anchored at the top, so crop exactly to the requested viewBox.
        with tempfile.TemporaryDirectory(prefix="gvf-phenotype-figure-") as tmp:
            tmp_dir = Path(tmp)
            preview_size = max(width, height)
            subprocess.run(
                [
                    qlmanage,
                    "-t",
                    "-s",
                    str(preview_size),
                    "-o",
                    str(tmp_dir),
                    str(svg_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            rendered = tmp_dir / f"{svg_path.name}.png"
            subprocess.run(
                [
                    ffmpeg,
                    "-y",
                    "-loglevel",
                    "error",
                    "-i",
                    str(rendered),
                    "-vf",
                    f"crop={width}:{height}:0:0",
                    "-frames:v",
                    "1",
                    str(png_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sips, "-s", "format", "pdf", str(svg_path), "--out", str(pdf_path)],
                check=True,
                capture_output=True,
                text=True,
            )
        return [png_path, pdf_path]

    if sips:
        with tempfile.TemporaryDirectory(prefix="gvf-phenotype-figure-") as tmp:
            intermediate = Path(tmp) / "phenotype_count_recovery.png"
            subprocess.run(
                [
                    sips,
                    "-s",
                    "format",
                    "png",
                    str(svg_path),
                    "--out",
                    str(intermediate),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [
                    sips,
                    "--resampleWidth",
                    str(width),
                    str(intermediate),
                    "--out",
                    str(png_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sips, "-s", "format", "pdf", str(svg_path), "--out", str(pdf_path)],
                check=True,
                capture_output=True,
                text=True,
            )
        return [png_path, pdf_path]

    print(
        "warning: no rsvg-convert or sips executable; retained editable SVG only",
        file=sys.stderr,
    )
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD)
    parser.add_argument("--svg-out", type=Path, default=DEFAULT_SVG)
    parser.add_argument("--png-out", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--pdf-out", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--csv-out", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON)
    args = parser.parse_args()

    lock = verify_lock(args.run_dir)
    selection = json.loads((args.run_dir / "selection.json").read_text())
    predictions = json.loads((args.run_dir / "predictions.json").read_text())
    report = json.loads((args.run_dir / "report.json").read_text())
    run_eval = load_run_eval()
    rows, missingness = paired_rows(
        run_eval,
        selection,
        predictions,
        args.gold_root,
        report["run_id"],
        lock["locked_at"],
    )
    summary = summarize(rows, missingness, report)
    identity_matched = int(report["overall"]["tp"])
    identity_missed = int(report["overall"]["fn"])
    identity_summary = {
        "matched": identity_matched,
        "missed": identity_missed,
        "gold_rows": identity_matched + identity_missed,
        "coverage": identity_matched / (identity_matched + identity_missed),
    }
    paired_keys = {
        field: {
            (row["gene"], row["pmid"], row["manual_variant"])
            for row in rows
            if row["measure"] == field
        }
        for field in FIELDS
    }
    paired_overlap = len(paired_keys["affected"] & paired_keys["unaffected"])

    write_csv(rows, args.csv_out)
    args.svg_out.parent.mkdir(parents=True, exist_ok=True)
    args.svg_out.write_text(
        render_zero_filled_svg(
            rows,
            summary,
            paired_overlap,
            report["run_id"],
            lock["locked_at"],
        ),
        encoding="utf-8",
    )
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(
        json.dumps(
            {
                "analysis_run": report["run_id"],
                "locked_at": lock["locked_at"],
                "integrity": lock["verified_hashes"],
                "row_grain": "identity-matched variant-paper row and count measure",
                "not_available": "explicit cohort and phenotype identifiers",
                "paired_rows_in_both_measures": paired_overlap,
                "variant_identity": identity_summary,
                "metrics": summary,
                "figure_contract": {
                    "question": (
                        "Among identity-matched variant-paper rows, how closely do "
                        "affected and unaffected counts agree when a missing automated "
                        "count is evaluated as zero?"
                    ),
                    "chart_family": (
                        "faceted square-root-scale bubble agreement plots with a "
                        "single variant-identity context line"
                    ),
                    "axes": "shared square-root count scale; 0 to 110 display domain",
                    "evaluation_missing_count": 0,
                    "storage_missing_count": None,
                    "marker_radius": "min(28, 8 + 2*sqrt(variant-paper pairs))",
                    "marker_size_key_pairs": list(MARKER_SIZE_KEY_COUNTS),
                    "labeled_variant_paper_pairs": [
                        {
                            "measure": measure,
                            "gene": gene,
                            "pmid": pmid,
                            "variant": variant,
                        }
                        for measure, gene, pmid, variant in EXAMPLE_ANNOTATIONS
                        if any(
                            row["measure"] == measure
                            and row["gene"] == gene
                            and row["pmid"] == pmid
                            and row["manual_variant"] == variant
                            for row in rows
                        )
                    ],
                    "renderer": "editable SVG; vector PDF and 300-dpi PNG downstream",
                    "fallback": (
                        "identity misses remain excluded; every identity-matched row "
                        "with an authoritative count is plotted"
                    ),
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    exported = export_svg_copies(args.svg_out, args.png_out, args.pdf_out)
    print(f"wrote {len(rows)} plotted rows to {args.csv_out}")
    print(f"wrote summary to {args.json_out}")
    print(f"wrote editable figure to {args.svg_out}")
    for path in exported:
        print(f"wrote figure copy to {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
