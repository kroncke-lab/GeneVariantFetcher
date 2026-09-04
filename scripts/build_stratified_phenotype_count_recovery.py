#!/usr/bin/env python3
"""Build one non-pooled view of every phenotype-count cohort opened to date.

The historical cardiac run and the mixed-gold candidate arms use different
sampling frames and projections. They therefore share a canvas and axes but
never a denominator, bubble, or performance summary.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import build_combined_phenotype_count_recovery as combined
from scripts import build_phenotype_count_recovery as figure


RUNS = REPO / "benchmarks" / "codex_paper_eval" / "runs"
HISTORICAL_RUN = RUNS / "20260902_false_zero_recovery_gold118"
OUT_DIR = REPO / "docs" / "figures" / "evaluated_phenotype_counts"
TITLE = "Automated versus reference phenotype counts across evaluated cohorts"

HISTORICAL_ANNOTATIONS = {
    key: (dx * 0.88, dy * 0.88) for key, (dx, dy) in figure.EXAMPLE_ANNOTATIONS.items()
}
HISTORICAL_DOCKS = figure.EXAMPLE_DOCKS
MIXED_ANNOTATIONS = {
    key: (dx * 0.82, dy * 0.82) for key, (dx, dy) in combined.POOLED_ANNOTATIONS.items()
}
MIXED_DOCKS = combined.POOLED_DOCKS


def read_cohort_rows(path: Path, cohort: str, source_id: str) -> list[dict[str, Any]]:
    rows = combined.read_rows(path, source_id)
    for row in rows:
        row["analysis_cohort"] = cohort
    return rows


def draw_cohort_row(
    *,
    rows: list[dict[str, Any]],
    row_id: str,
    label: str,
    plot_top: float,
    annotations: dict[tuple[str, str, str, str], tuple[float, float]],
    docks: dict[tuple[str, str, str, str], str],
) -> list[str]:
    plot_size = 810.0
    plot_x = {"affected": 145.0, "unaffected": 1205.0}
    label_y = plot_top - 90.0
    facet_y = plot_top - 37.0
    output = [figure.text_element(80, label_y, label, "cohort-title")]
    ticks = (0, 1, 2, 5, 10, 20, 50, 100)
    for field in figure.FIELDS:
        x0 = plot_x[field]
        output.append(
            figure.text_element(
                x0,
                facet_y,
                f"{field.capitalize()} individuals",
                "facet-title",
            )
        )
        for tick in ticks:
            x = figure.plot_coordinate(tick, x0, plot_size)
            y = figure.plot_y(tick, plot_top, plot_size)
            output.extend(
                [
                    f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" '
                    f'y2="{plot_top + plot_size}" class="grid"/>',
                    f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + plot_size}" '
                    f'y2="{y:.1f}" class="grid"/>',
                    figure.text_element(
                        x,
                        plot_top + plot_size + 36,
                        str(tick),
                        "tick",
                        anchor="middle",
                    ),
                ]
            )
            if field == "affected":
                output.append(
                    figure.text_element(x0 - 19, y + 9, str(tick), "tick", anchor="end")
                )
        zero_x = figure.plot_coordinate(0, x0, plot_size)
        zero_y = figure.plot_y(0, plot_top, plot_size)
        end_x = figure.plot_coordinate(figure.PLOT_MAX, x0, plot_size)
        end_y = figure.plot_y(figure.PLOT_MAX, plot_top, plot_size)
        output.extend(
            [
                f'<rect x="{x0}" y="{plot_top}" width="{plot_size}" '
                f'height="{plot_size}" class="axis"/>',
                f'<line x1="{zero_x:.1f}" y1="{zero_y:.1f}" '
                f'x2="{end_x:.1f}" y2="{end_y:.1f}" class="identity"/>',
                figure.text_element(
                    x0 + plot_size / 2,
                    plot_top + plot_size + 82,
                    "Manual reference count",
                    "axis-label",
                    anchor="middle",
                ),
                f'<g clip-path="url(#{row_id}-{field}-clip)">',
            ]
        )
        for group in sorted(
            figure.grouped_coordinates(rows, field),
            key=lambda item: item["n"],
            reverse=True,
        ):
            output.extend(figure.draw_bubble(group, field, x0, plot_top, plot_size))
        output.append("</g>")
        output.extend(
            figure.draw_example_annotations(
                rows,
                field,
                x0,
                plot_top,
                plot_size,
                annotations=annotations,
                docks=docks,
            )
        )
    output.append(
        figure.text_element(
            36,
            plot_top + plot_size / 2,
            "Automated count (missing scored 0)",
            "axis-label",
            anchor="middle",
            transform=f"rotate(-90 36 {plot_top + plot_size / 2})",
        )
    )
    return output


def render_svg(
    historical: list[dict[str, Any]],
    mixed: list[dict[str, Any]],
    *,
    historical_attempt_count: int,
    mixed_attempt_count: int,
) -> str:
    width, height = 2160, 2250
    plot_size = 810.0
    plot_tops = {"historical": 235.0, "mixed": 1260.0}
    plot_x = {"affected": 145.0, "unaffected": 1205.0}
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" '
            f'height="7.5in" viewBox="0 0 {width} {height}" role="img" '
            'aria-labelledby="figure-title figure-desc">'
        ),
        f'<title id="figure-title">{figure.escape(TITLE)}</title>',
        (
            '<desc id="figure-desc">Two separately scored cohorts share '
            "square-root axes. The historical cardiac cohort and opened mixed-gold "
            "candidate runs are not pooled. Missing automated counts are evaluated "
            "as zero only for scoring and display.</desc>"
        ),
        "<metadata>",
        "Historical source: 20260902_false_zero_recovery_gold118",
        "Mixed sources: scored candidate arms listed in versioned mixed-gold consumption logs",
        "Overlap rule: cohorts remain separate; no deduplication or pooled metric",
        "</metadata>",
        "<defs>",
        "<style><![CDATA[",
        (
            'text { font-family: "Avenir Next Condensed", "Source Sans 3", '
            '"Helvetica Neue", Arial, sans-serif; fill: #17324D; }'
        ),
        ".title { font-size: 54px; font-weight: 600; letter-spacing: -0.3px; }",
        ".cohort-title { font-size: 38px; font-weight: 600; fill: #40566D; }",
        ".facet-title { font-size: 36px; font-weight: 600; }",
        ".tick { font-size: 27px; font-weight: 400; fill: #687584; }",
        ".axis-label { font-size: 30px; font-weight: 500; }",
        ".example-title { font-size: 27px; font-weight: 600; }",
        ".example-pmid { font-size: 24px; font-weight: 500; fill: #687584; }",
        ".example-line { stroke: #687584; stroke-width: 2.3; fill: none; }",
        ".size-key-box { fill: #FFFFFF; fill-opacity: 0.94; stroke: #DCE3E8; stroke-width: 2.2; }",
        ".size-key-title { font-size: 25px; font-weight: 600; }",
        ".size-key-value { font-size: 25px; font-weight: 600; fill: #687584; }",
        ".size-key-marker { fill: #FFFFFF; stroke: #17324D; stroke-width: 2.5; }",
        ".grid { stroke: #DCE3E8; stroke-width: 2.2; }",
        ".axis { stroke: #17324D; stroke-width: 3.5; fill: none; }",
        ".identity { stroke: #17324D; stroke-width: 4; stroke-dasharray: 12 10; }",
        "]]></style>",
    ]
    for row_id, top in plot_tops.items():
        for field in figure.FIELDS:
            lines.append(
                f'<clipPath id="{row_id}-{field}-clip"><rect '
                f'x="{plot_x[field]}" y="{top}" width="{plot_size}" '
                f'height="{plot_size}"/></clipPath>'
            )
    lines.extend(
        [
            "</defs>",
            f'<rect width="{width}" height="{height}" fill="{figure.COLORS["off_white"]}"/>',
            figure.text_element(80, 78, TITLE, "title"),
        ]
    )
    lines.extend(
        draw_cohort_row(
            rows=historical,
            row_id="historical",
            label=(
                "Historical cardiac evaluation · "
                f"{historical_attempt_count} gene–paper attempts"
            ),
            plot_top=plot_tops["historical"],
            annotations=HISTORICAL_ANNOTATIONS,
            docks=HISTORICAL_DOCKS,
        )
    )
    lines.extend(
        figure.draw_marker_size_key(plot_x["affected"], plot_tops["historical"])
    )
    lines.extend(
        draw_cohort_row(
            rows=mixed,
            row_id="mixed",
            label=(
                "Opened mixed-gold candidate arms · "
                f"{mixed_attempt_count} gene–paper attempts"
            ),
            plot_top=plot_tops["mixed"],
            annotations=MIXED_ANNOTATIONS,
            docks=MIXED_DOCKS,
        )
    )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def attempt_keys(run_dirs: list[Path]) -> set[tuple[str, str]]:
    keys: set[tuple[str, str]] = set()
    for run_dir in run_dirs:
        selection = json.loads((run_dir / "selection.json").read_text())
        keys.update(
            (paper["gene"], str(paper["pmid"])) for paper in selection["papers"]
        )
    return keys


def build(out_dir: Path = OUT_DIR) -> dict[str, Path]:
    historical_lock = figure.verify_lock(HISTORICAL_RUN)
    historical_data = HISTORICAL_RUN / "figures" / "data"
    historical_summary = json.loads(
        (historical_data / "phenotype_count_recovery.json").read_text()
    )
    if historical_summary["integrity"] != historical_lock["verified_hashes"]:
        raise SystemExit("historical figure integrity does not match its lock")
    historical = read_cohort_rows(
        historical_data / "phenotype_count_recovery.csv",
        "historical_cardiac_gold118",
        HISTORICAL_RUN.name,
    )

    candidate_entries = combined.candidate_runs_from_logs(combined.DEFAULT_LOGS, RUNS)
    mixed: list[dict[str, Any]] = []
    mixed_summaries: list[dict[str, Any]] = []
    for tier_id, run_dir in candidate_entries:
        lock = figure.verify_lock(run_dir)
        data_dir = run_dir / "figures" / "data"
        summary = json.loads((data_dir / "phenotype_count_recovery.json").read_text())
        if summary["integrity"] != lock["verified_hashes"]:
            raise SystemExit(f"figure integrity does not match lock for {run_dir.name}")
        mixed.extend(
            read_cohort_rows(
                data_dir / "phenotype_count_recovery.csv",
                "opened_mixed_gold_candidates",
                tier_id,
            )
        )
        mixed_summaries.append(summary)

    historical_keys = {
        (row["measure"], row["gene"], row["pmid"], row["manual_variant"]): row
        for row in historical
    }
    mixed_keys = {
        (row["measure"], row["gene"], row["pmid"], row["manual_variant"]): row
        for row in mixed
    }
    overlap = set(historical_keys) & set(mixed_keys)
    count_disagreements = sum(
        historical_keys[key]["ai_count_evaluation"]
        != mixed_keys[key]["ai_count_evaluation"]
        for key in overlap
    )
    historical_attempts = attempt_keys([HISTORICAL_RUN])
    mixed_run_dirs = [run_dir for _, run_dir in candidate_entries]
    mixed_attempts = attempt_keys(mixed_run_dirs)

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "svg": out_dir / "phenotype_count_recovery_stratified.svg",
        "png": out_dir / "phenotype_count_recovery_stratified.png",
        "pdf": out_dir / "phenotype_count_recovery_stratified.pdf",
        "csv": out_dir / "phenotype_count_recovery_stratified.csv",
        "json": out_dir / "phenotype_count_recovery_stratified.json",
    }
    with paths["csv"].open("w", newline="", encoding="utf-8") as handle:
        all_rows = historical + mixed
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    paths["svg"].write_text(
        render_svg(
            historical,
            mixed,
            historical_attempt_count=len(historical_attempts),
            mixed_attempt_count=len(mixed_attempts),
        ),
        encoding="utf-8",
    )
    paths["json"].write_text(
        json.dumps(
            {
                "title": TITLE,
                "cohorts": {
                    "historical_cardiac_gold118": {
                        "run_id": HISTORICAL_RUN.name,
                        "gene_paper_attempts": len(historical_attempts),
                        "unique_pmids": len({pmid for _, pmid in historical_attempts}),
                        "plotted_rows": len(historical),
                        "metrics": historical_summary["metrics"],
                    },
                    "opened_mixed_gold_candidates": {
                        "run_ids": [run_dir.name for run_dir in mixed_run_dirs],
                        "tier_ids": [tier_id for tier_id, _ in candidate_entries],
                        "gene_paper_attempts": len(mixed_attempts),
                        "unique_pmids": len({pmid for _, pmid in mixed_attempts}),
                        "plotted_rows": len(mixed),
                        "protocol_note": (
                            "cumulative progress/failure view; candidate protocol may "
                            "change across opened disjoint tranches"
                        ),
                        "source_metrics": [item["metrics"] for item in mixed_summaries],
                    },
                },
                "overlap": {
                    "gene_paper_attempts": len(historical_attempts & mixed_attempts),
                    "plotted_measure_variant_keys": len(overlap),
                    "automated_count_disagreements": count_disagreements,
                    "rule": "shown in separate strata; never deduplicated or pooled",
                },
                "evaluation_rule": {
                    "missing_automated_count": 0,
                    "stored_missing_value": None,
                    "identity_misses": "excluded from plots",
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    figure.export_svg_copies(
        paths["svg"], paths["png"], paths["pdf"], width=2160, height=2250
    )
    return paths


def main() -> int:
    paths = build()
    for label, path in paths.items():
        print(f"wrote {label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
