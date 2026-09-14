"""Check plotted values and rendered text bounds without overwriting figures."""

import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parents[4] / "tmp/matplotlib")
)

import matplotlib.figure
import numpy as np
import pandas as pd
from matplotlib.text import Text

from audit_refresh import GEOMETRY, REFRESH, read, sha

HERE = Path(__file__).resolve().parent


def main():
    spec = importlib.util.spec_from_file_location(
        "reviewed_plots", REFRESH / "plot_residues.py"
    )
    plots = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plots)
    saved_figures = []
    original_save = matplotlib.figure.Figure.savefig
    receipts = []
    try:
        matplotlib.figure.Figure.savefig = lambda fig, *args, **kwargs: (
            saved_figures.append(fig)
        )
        for gene in GEOMETRY:
            data = pd.read_csv(REFRESH / "tables" / f"{gene}_residue_density.csv")
            frame = read(REFRESH / "analysis" / gene, "density_h3")
            plots.plot_individual(gene, data, frame)
            figure = saved_figures[-1]
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            canvas = figure.bbox
            bad = []
            for item in figure.findobj(match=Text):
                if not item.get_visible() or not item.get_text().strip():
                    continue
                # Matplotlib keeps invisible ticks outside the axis limits in its object tree.
                if item in sum(
                    [
                        axis.get_xticklabels() + axis.get_yticklabels()
                        for axis in figure.axes
                    ],
                    [],
                ):
                    box = item.get_window_extent(renderer)
                    if (
                        box.x1 < 0
                        or box.x0 > canvas.x1
                        or box.y1 < 0
                        or box.y0 > canvas.y1
                    ):
                        continue
                box = item.get_window_extent(renderer)
                if (
                    box.x0 < -1
                    or box.y0 < -1
                    or box.x1 > canvas.x1 + 1
                    or box.y1 > canvas.y1 + 1
                ):
                    bad.append(dict(text=item.get_text(), bbox=list(box.bounds)))
            assert not bad, (gene, bad)
            log_axis = figure.axes[1]
            positive = np.concatenate(
                [
                    data.loc[data[column].gt(0), column].to_numpy()
                    for column in [
                        "density",
                        "raw_variant_fraction_density",
                        "raw_kernel_one_prior_posterior",
                    ]
                ]
            )
            low, high = log_axis.get_ylim()
            assert low <= positive.min() and high >= positive.max(), (
                gene,
                low,
                positive.min(),
                high,
                positive.max(),
            )
            assert log_axis.get_yscale() == "log"
            # Root's displayed primary array contains NaNs, so line paths break at missing residues.
            line = figure.axes[0].lines[0]
            np.testing.assert_allclose(line.get_xdata(), data.canonical_pos)
            np.testing.assert_allclose(line.get_ydata(), data.density, equal_nan=True)
            own_line = figure.axes[2].lines[0]
            np.testing.assert_allclose(
                own_line.get_ydata(), data.own_posterior_mean, equal_nan=True
            )
            receipts.append(
                dict(
                    gene=gene,
                    visible_text_within_canvas=True,
                    all_positive_diagnostic_values_inside_log_axis=True,
                    positive_min=float(positive.min()),
                    log_axis_min=float(low),
                    primary_line_matches_table_with_gaps=True,
                    own_posterior_panel_matches_table=True,
                )
            )
            saved_figures.clear()
    finally:
        matplotlib.figure.Figure.savefig = original_save
    result = dict(
        individual_figures=receipts,
        audit_script_sha256=sha(__file__),
        plot_script_sha256=sha(REFRESH / "plot_residues.py"),
        outputs_overwritten=False,
    )
    (HERE / "plot_layout_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
