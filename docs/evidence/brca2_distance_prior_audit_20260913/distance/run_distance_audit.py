"""Trace frozen BRCA2 geometry, positive kernels, and normalized donor mass."""

import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parents[1]
REPO = EVIDENCE.parents[1]
FROZEN = EVIDENCE / "structural_sanity_20260913"
sys.path.insert(0, str(REPO.parent / "ProteinProximityAnalysis/src"))
from alphafold_rin.empirical_density import (
    empirical_variant_density,
    log_sigmoid_kernel,
    sigmoid_kernel,
)

spec = importlib.util.spec_from_file_location("distance_input", FROZEN / "run_brca2.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def save(frame, name):
    content = frame.to_csv(index=False, lineterminator="\n").encode()
    if name.endswith(".gz"):
        content = gzip.compress(content, mtime=0)
    assert len(content) < 1_150_000, (name, len(content))
    (HERE / name).write_bytes(content)


def shards(frame, stem, rows=2000):
    for number, start in enumerate(range(0, len(frame), rows), 1):
        save(frame.iloc[start : start + rows], f"{stem}.part{number:03d}.csv.gz")


def weighted_quantile(distance, weight, probabilities):
    order = np.argsort(distance, kind="stable")
    cumulative = np.cumsum(weight[order])
    cumulative /= cumulative[-1]
    indices = np.minimum(np.searchsorted(cumulative, probabilities), len(order) - 1)
    return distance[order[indices]]


def main():
    variants, prior, input_paths = source.frozen.load_variants("BRCA2")
    geometry_path = FROZEN / "geometry/BRCA2/primary_geometry.csv.gz"
    geometry = pd.read_csv(geometry_path, low_memory=False).fillna({"idr_segment": ""})
    primary = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted(
                (FROZEN / "analysis/BRCA2").glob("primary_density.part*.csv.gz")
            )
        ],
        ignore_index=True,
    )
    ids = variants.variant_id.tolist()
    assert primary.variant_id.tolist() == ids
    y = variants.posterior_mean.to_numpy()
    affected = variants.affected.to_numpy() > 0
    fraction = variants.affected.to_numpy() / variants.n.to_numpy()
    target_rows, context_rows = [], []
    thresholds = [3, 6, 10, 20, 30]
    for start in range(0, len(ids), 128):
        targets = ids[start : start + 128]
        result = empirical_variant_density(
            variants,
            geometry,
            target_ids=targets,
            include_context_weights=False,
            include_context_model=True,
        )
        model = result.context_model
        assert list(model.donor_ids) == ids
        selected = geometry.loc[
            geometry.geometry_state.isin(["structured", "idr"])
            & geometry.canonical_pos.isin(
                variants.canonical_pos.iloc[start : start + len(targets)]
            )
        ].reset_index(drop=True)
        for ti, target in enumerate(targets):
            contexts = np.flatnonzero(model.context_target_rows == ti)
            distances, normalized, per_context, bandwidth = [], [], [], []
            used_affected = set()
            nearest_affected = np.inf
            own = start + ti
            for context in contexts:
                if not np.isfinite(model.context_log_sums[context]):
                    continue
                grow = model.context_geometry_rows[context]
                record = selected.iloc[grow]
                lk = model.position_log_kernels[
                    grow, model.donor_position_indices
                ].copy()
                lk[own] = -np.inf
                available = np.isfinite(lk)
                indices = np.flatnonzero(available)
                logs = lk[available]
                # Invert K(d)=2/(1+exp(log(3)*d/3)) without overflow.
                distance = (
                    (np.log(2.0) - logs + np.log1p(-np.exp(logs - np.log(2.0))))
                    * 3
                    / np.log(3.0)
                )
                distance = np.maximum(distance, 0.0)
                q = np.exp(logs - model.context_log_sums[context])
                np.testing.assert_allclose(q.sum(), 1.0, atol=1e-12)
                local_affected = affected[indices]
                if local_affected.any():
                    nearest_affected = min(
                        nearest_affected, float(distance[local_affected].min())
                    )
                    used_affected.update(indices[local_affected].tolist())
                row = {
                    "variant_id": target,
                    "context_id": json.dumps(
                        [record.frame_id, record.chain, int(record.canonical_pos)],
                        separators=(",", ":"),
                    ),
                    "frame_id": record.frame_id,
                    "chain": record.chain,
                    "canonical_pos": record.canonical_pos,
                    "geometry_state": record.geometry_state,
                    "idr_segment": record.idr_segment,
                    "raw_kernel_sum": float(np.exp(model.context_log_sums[context])),
                    "nearest_donor_distance": float(distance.min()),
                    "affected_weight_share": float(q[local_affected].sum()),
                    "zero_affected_weight_share": float(q[~local_affected].sum()),
                    "affected_density_component": float(
                        q[local_affected] @ y[indices[local_affected]]
                    ),
                    "zero_affected_density_component": float(
                        q[~local_affected] @ y[indices[~local_affected]]
                    ),
                    "raw_count_neighborhood": float(q @ fraction[indices]),
                    "density_h3": float(q @ y[indices]),
                    "donors": len(indices),
                }
                for threshold in thresholds:
                    row[f"weight_beyond_{threshold}"] = float(
                        q[distance > threshold].sum()
                    )
                    row[f"raw_kernel_beyond_{threshold}"] = float(
                        np.exp(logs[distance > threshold]).sum()
                    )
                median, p90 = weighted_quantile(distance, q, [0.5, 0.9])
                row["weighted_distance_median"], row["weighted_distance_p90"] = (
                    float(median),
                    float(p90),
                )
                alternatives = {}
                for h in [1, 2, 3, 5]:
                    alt_logs = log_sigmoid_kernel(distance, h)
                    alt_q = np.exp(alt_logs - logsumexp(alt_logs))
                    alternatives[f"density_h{h}"] = float(alt_q @ y[indices])
                bandwidth.append(alternatives)
                distances.append(distance)
                normalized.append(q)
                per_context.append(row)
                context_rows.append(row)
            row = {
                "variant_id": target,
                "canonical_pos": int(variants.canonical_pos.iloc[own]),
                "affected": variants.affected.iloc[own],
                "unaffected": variants.unaffected.iloc[own],
                "posterior_mean": y[own],
                "density_source": primary.density_source.iloc[own],
                "supported_contexts": len(per_context),
            }
            if per_context:
                table = pd.DataFrame(per_context)
                row.update(
                    table.select_dtypes(include="number")
                    .drop(columns=["canonical_pos"])
                    .mean()
                    .to_dict()
                )
                row.update(pd.DataFrame(bandwidth).mean().to_dict())
                row["mean_context_nearest_donor_distance"] = row.pop(
                    "nearest_donor_distance"
                )
                row["nearest_donor_distance"] = table.nearest_donor_distance.min()
                dd, qq = (
                    np.concatenate(distances),
                    np.concatenate(normalized) / len(per_context),
                )
                median, p90 = weighted_quantile(dd, qq, [0.5, 0.9])
                row["weighted_distance_median"], row["weighted_distance_p90"] = (
                    float(median),
                    float(p90),
                )
                (
                    row["eligible_pair_distance_median"],
                    row["eligible_pair_distance_p90"],
                ) = np.quantile(dd, [0.5, 0.9])
                row["minimum_context_kernel_sum"] = table.raw_kernel_sum.min()
                row["maximum_context_weight_beyond_20"] = table.weight_beyond_20.max()
                row["nearest_affected_distance"] = (
                    nearest_affected if np.isfinite(nearest_affected) else np.nan
                )
                row["affected_donors_any_context"] = len(used_affected)
            target_rows.append(row)
        if start % 1280 == 0:
            print(f"Distance audit {start + len(targets)}/{len(ids)}", flush=True)
    targets, contexts = pd.DataFrame(target_rows), pd.DataFrame(context_rows)
    np.testing.assert_allclose(
        targets.density_h3, primary.density, atol=2e-12, equal_nan=True
    )
    np.testing.assert_allclose(
        targets.weight_beyond_20,
        primary.weight_share_beyond_20,
        atol=2e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        targets.raw_kernel_sum.fillna(0), primary.sum_kernel_weight, atol=2e-12
    )
    old_scenarios = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted(
                (FROZEN / "analysis/BRCA2").glob("density_scenarios.part*.csv.gz")
            )
        ],
        ignore_index=True,
    )
    for h in [2, 3, 5]:
        saved = (
            old_scenarios.loc[old_scenarios.scenario.eq(f"primary_com_h{h}")]
            .set_index("variant_id")
            .reindex(ids)
        )
        np.testing.assert_allclose(
            targets[f"density_h{h}"], saved.density, atol=2e-12, equal_nan=True
        )
    shards(targets, "variant_distance")
    shards(contexts, "context_distance")
    segments = []
    for (frame, chain, segment), group in geometry.loc[
        geometry.geometry_state.eq("idr")
    ].groupby(["frame_id", "chain", "idr_segment"]):
        members = variants.loc[variants.canonical_pos.isin(group.canonical_pos)]
        selected_targets = targets.loc[targets.variant_id.isin(members.variant_id)]
        segments.append(
            {
                "frame_id": frame,
                "chain": chain,
                "idr_segment": segment,
                "start": group.canonical_pos.min(),
                "end": group.canonical_pos.max(),
                "variants": len(members),
                "affected_variants": int(members.affected.gt(0).sum()),
                "affected": members.affected.sum(),
                "unaffected": members.unaffected.sum(),
                "median_posterior": members.posterior_mean.median()
                if len(members)
                else np.nan,
                "supported_density_targets": selected_targets.density_h3.notna().sum(),
                "median_density": selected_targets.density_h3.median()
                if selected_targets.density_h3.notna().any()
                else np.nan,
            }
        )
    segments = pd.DataFrame(segments)
    save(segments, "idr_segment_counts.csv")
    measures = [
        "weight_beyond_3",
        "weight_beyond_6",
        "weight_beyond_10",
        "weight_beyond_20",
        "weight_beyond_30",
        "weighted_distance_median",
        "weighted_distance_p90",
        "eligible_pair_distance_median",
        "eligible_pair_distance_p90",
        "raw_kernel_sum",
        "minimum_context_kernel_sum",
        "affected_weight_share",
        "zero_affected_weight_share",
        "nearest_affected_distance",
        "affected_density_component",
        "zero_affected_density_component",
        "raw_count_neighborhood",
    ]
    summaries = []
    for group_name, group in targets.loc[targets.density_h3.notna()].groupby(
        "density_source"
    ):
        for measure in measures:
            value = group[measure].dropna()
            summaries.append(
                {
                    "density_source": group_name,
                    "measure": measure,
                    "targets": len(group),
                    "nonmissing": len(value),
                    "mean": value.mean(),
                    "median": value.median(),
                    "p90": value.quantile(0.9),
                    "maximum": value.max(),
                    "minimum": value.min(),
                }
            )
    save(pd.DataFrame(summaries), "distance_summary.csv")
    sensitivity = []
    for group_name, group in targets.loc[targets.density_h3.notna()].groupby(
        "density_source"
    ):
        for h in [1, 2, 3, 5]:
            value = group[f"density_h{h}"]
            sensitivity.append(
                {
                    "density_source": group_name,
                    "half_distance": h,
                    "targets": len(group),
                    "median_density": value.median(),
                    "q25": value.quantile(0.25),
                    "q75": value.quantile(0.75),
                    "mean_abs_change_from_h3": abs(value - group.density_h3).mean(),
                    "max_abs_change_from_h3": abs(value - group.density_h3).max(),
                    "below_0_001": int((value <= 0.001).sum()),
                }
            )
    sensitivity = pd.DataFrame(sensitivity)
    save(sensitivity, "bandwidth_sensitivity.csv")
    examples = pd.concat(
        [
            targets.nlargest(3, "weight_beyond_20"),
            targets.nsmallest(3, "minimum_context_kernel_sum"),
            targets.loc[
                targets.density_h3.notna() & targets.affected_donors_any_context.eq(0)
            ].head(5),
        ]
    ).drop_duplicates("variant_id")
    save(examples, "concrete_examples.csv")
    check_ids = examples.variant_id.head(6).tolist()
    usable = geometry.loc[geometry.geometry_state.isin(["structured", "idr"])]
    fresh = empirical_variant_density(
        variants,
        usable,
        target_ids=check_ids,
        include_context_weights=True,
        backend="reference",
    )
    checks = []
    for (target, context), group in fresh.context_weights.groupby(
        ["variant_id", "context_id"]
    ):
        recorded = contexts.loc[
            contexts.variant_id.eq(target) & contexts.context_id.eq(context)
        ].iloc[0]
        direct_raw = 2 / (1 + np.exp(np.log(3) * group.distance.to_numpy() / 3))
        normalized_direct = direct_raw / direct_raw.sum()
        np.testing.assert_allclose(
            normalized_direct, group.context_normalized_weight, atol=2e-12
        )
        for cutoff in [10, 20, 30]:
            np.testing.assert_allclose(
                normalized_direct[group.distance.to_numpy() > cutoff].sum(),
                recorded[f"weight_beyond_{cutoff}"],
                atol=2e-12,
            )
        checks.append(
            {
                "variant_id": target,
                "context_id": context,
                "donors": len(group),
                "raw_kernel_sum": direct_raw.sum(),
            }
        )
    save(pd.DataFrame(checks), "independent_context_checks.csv")
    top = (
        fresh.context_weights.sort_values("context_normalized_weight", ascending=False)
        .groupby(["variant_id", "context_id"], sort=False)
        .head(15)
    )
    top = top.merge(
        variants[["variant_id", "affected", "unaffected", "posterior_mean"]],
        left_on="donor_id",
        right_on="variant_id",
        suffixes=("", "_donor"),
    )
    shards(top, "example_top_donors", rows=1500)
    zero_cases = targets.loc[
        targets.density_h3.notna() & targets.affected_donors_any_context.eq(0)
    ]
    receipt = {
        "prior": prior,
        "variants": len(variants),
        "supported": int(targets.density_h3.notna().sum()),
        "supported_contexts": len(contexts),
        "contexts_mass_beyond20_over_half": int(
            contexts.weight_beyond_20.gt(0.5).sum()
        ),
        "maximum_context_mass_beyond20": contexts.weight_beyond_20.max(),
        "minimum_context_raw_kernel_sum": contexts.raw_kernel_sum.min(),
        "targets_with_no_affected_donors": len(zero_cases),
        "zero_affected_donor_density_median": zero_cases.density_h3.median(),
        "zero_case_idr_segments_with_variants": int(
            ((segments.affected == 0) & (segments.variants > 0)).sum()
        ),
        "zero_case_idr_segments_supported_targets": int(
            segments.loc[segments.affected.eq(0), "supported_density_targets"].sum()
        ),
        "independent_explicit_contexts_checked": len(checks),
        "cached_h2_h3_h5_reproduced": True,
        "inputs": {
            str(path.relative_to(EVIDENCE)): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in [*input_paths, geometry_path]
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE / "checks.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    )
    make_figure(targets, sensitivity)
    print(
        json.dumps(
            {
                key: value
                for key, value in receipt.items()
                if key not in ["inputs", "prior"]
            },
            indent=2,
        )
    )


def make_figure(targets, sensitivity):
    os.environ.setdefault("MPLCONFIGDIR", str(REPO / "tmp/matplotlib_distance_audit"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    colors = {"polymer": "#3979a8", "structured": "#b66732"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    d = np.linspace(0, 40, 500)
    for h in [1, 2, 3, 5]:
        axes[0].plot(
            d, sigmoid_kernel(d, h), label=f"h={h} Å", linewidth=2 if h == 3 else 1
        )
    axes[0].set(
        xlabel="Donor distance (Å)",
        ylabel="Raw kernel weight",
        yscale="log",
        ylim=(1e-8, 1),
        title="Raw weights decay steeply",
    )
    axes[0].legend(frameon=False)
    bins = ["≤3", "3–6", "6–10", "10–20", "20–30", ">30"]
    shades = ["#163e63", "#3979a8", "#7fa9c8", "#abc6d9", "#d4dfc7", "#e3a47c"]
    for i, name in enumerate(["polymer", "structured"]):
        group = targets.loc[targets.density_source.eq(name)]
        tails = np.r_[
            1.0, [group[f"weight_beyond_{x}"].mean() for x in [3, 6, 10, 20, 30]], 0.0
        ]
        left = 0.0
        for mass, label, color in zip(-np.diff(tails), bins, shades):
            axes[1].barh(
                i, mass, left=left, color=color, label=label if i == 0 else None
            )
            left += mass
    axes[1].set(
        yticks=[0, 1],
        yticklabels=["Polymer", "Structured"],
        xlim=(0, 1),
        xlabel="Mean normalized donor mass",
        title="Most normalized mass stays nearby",
    )
    axes[1].xaxis.set_major_formatter(PercentFormatter(1))
    axes[1].legend(
        title="Distance (Å)",
        fontsize=8,
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.28),
    )
    for name, group in sensitivity.groupby("density_source"):
        axes[2].plot(
            group.half_distance,
            group.median_density,
            marker="o",
            label=name.capitalize(),
            color=colors[name],
        )
        axes[2].fill_between(
            group.half_distance, group.q25, group.q75, alpha=0.15, color=colors[name]
        )
    axes[2].set(
        xlabel="Kernel half-distance h (Å)",
        ylabel="Median neighborhood density",
        xticks=[1, 2, 3, 5],
        ylim=(0, 0.3),
        title="Shorter kernels barely shift the median",
    )
    axes[2].yaxis.set_major_formatter(PercentFormatter(1))
    axes[2].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "BRCA2: distance decay and normalization with frozen missense posteriors",
        fontsize=14,
    )
    fig.text(
        0.01,
        0.01,
        "6,326 supported targets; equal variant votes and equal supported contexts. Shading: variant interquartile range. No kernel cutoff. h=1 is a new sensitivity.",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.29, top=0.8, wspace=0.35)
    fig.savefig(HERE / "DISTANCE_AUDIT.png", dpi=160)


if __name__ == "__main__":
    main()
