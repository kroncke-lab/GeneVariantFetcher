"""Independent geometry-to-weight and selected outer-LOO reconstruction.

Does not call the PPA density engine. Frozen fractional-logistic optimization is
reused only after independently reconstructing and aligning its input features.
"""

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
RUN = HERE.parent
spec = importlib.util.spec_from_file_location(
    "frozen_fitter", RUN.parent / "gck_structural_pilot_20260912/pilot_statistics.py"
)
fitter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitter)


def read_parts(folder, stem):
    return pd.concat(
        [pd.read_csv(p) for p in sorted(folder.glob(f"{stem}.part*.csv.gz"))],
        ignore_index=True,
    )


def kernel(distance):
    return 2 * np.exp(-np.logaddexp(0, np.log(3) * distance / 3))


def geometry_kernels(variants, geometry):
    """One unnormalized row per actual residue/chain context, no copy counts."""
    pos = variants.canonical_pos.to_numpy(dtype=int)
    gp = geometry.canonical_pos.to_numpy(dtype=int)
    states = geometry.geometry_state.to_numpy()
    frames = geometry.frame_id.to_numpy()
    chains = geometry.chain.to_numpy()
    segments = geometry.idr_segment.fillna("").to_numpy()
    xyz = geometry[[f"com_{a}" for a in "xyz"]].to_numpy(dtype=float)
    answer = np.zeros((len(geometry), len(variants)))
    for i in range(len(geometry)):
        if states[i] == "structured":
            candidates = (states == "structured") & (frames == frames[i])
            distances = np.linalg.norm(xyz[candidates] - xyz[i], axis=1)
        elif states[i] == "idr":
            candidates = (
                (states == "idr")
                & (frames == frames[i])
                & (chains == chains[i])
                & (segments == segments[i])
            )
            distances = 3.8 * np.sqrt(abs(gp[candidates] - gp[i]))
        else:
            continue
        minimum = np.full(int(gp.max()), np.inf)
        np.minimum.at(minimum, gp[candidates] - 1, distances)
        answer[i] = kernel(minimum[pos - 1])
    return answer


def normalized_weights(variants, geometry, raw, excluded=None):
    output = np.zeros((len(variants), len(variants)))
    positions = geometry.canonical_pos.to_numpy()
    for i, pos in enumerate(variants.canonical_pos):
        contexts = raw[positions == pos].copy()
        contexts[:, i] = 0
        if excluded is not None:
            contexts[:, excluded] = 0
        sums = contexts.sum(axis=1)
        keep = sums > 0
        if keep.any():
            output[i] = (contexts[keep] / sums[keep, None]).mean(axis=0)
    return output


def means(weights, y):
    value = weights @ y
    value[weights.sum(axis=1) == 0] = np.nan
    return value


def audit_gene(gene):
    folder = RUN / "analysis" / gene
    variants = read_parts(folder, "input_variants")
    primary = read_parts(folder, "primary_density").set_index("variant_id")
    predictions = read_parts(folder, "loo_predictions")
    saved_weights = read_parts(folder, "primary_weights")
    geometry_name = {"HNF1A": "8PI8", "KCNQ1": "9U7F"}[gene]
    geometry_path = RUN / "geometry" / gene / f"{geometry_name}_canonical_geometry.csv"
    geometry = pd.read_csv(geometry_path)
    assert variants.variant_id.equals(variants.unit_id)
    assert variants.variant_id.equals(variants.canonical_variant_id)
    assert not variants.variant_id.duplicated().any()
    alleles = variants.member_alleles.dropna().str.split(";").explode()
    assert not alleles.duplicated().any(), "A genomic allele occurs in multiple units"
    assert variants.variant_type.eq("missense").all()
    assert variants.canonical_wt_status.eq("match").all()
    assert variants.alpha_empirical.nunique() == variants.beta_empirical.nunique() == 1
    ids = variants.variant_id.to_numpy()
    positions = variants.canonical_pos.to_numpy()
    y = variants.posterior_mean.to_numpy()
    raw = geometry_kernels(variants, geometry)
    weights = normalized_weights(variants, geometry, raw)
    d0 = means(weights, y)
    saved = (
        saved_weights.pivot(
            index="target_id", columns="donor_id", values="normalized_weight"
        )
        .reindex(index=ids, columns=ids)
        .fillna(0)
        .to_numpy()
    )
    np.testing.assert_allclose(saved, weights, atol=2e-13)
    np.testing.assert_allclose(d0, primary.loc[ids].density, atol=2e-13, equal_nan=True)
    assert (np.diag(weights) == 0).all()
    sequence_raw = kernel(3.8 * np.sqrt(abs(positions[:, None] - positions[None, :])))
    np.fill_diagonal(sequence_raw, 0)
    sequence = sequence_raw / sequence_raw.sum(axis=1, keepdims=True)
    s0 = means(sequence, y)
    selected = []
    for source in ["structured", "polymer", "mixed"]:
        keys = primary.index[primary.density_source.eq(source)].tolist()
        if keys:
            selected.append(keys[0])
    for model, values in [("density_unfitted", d0), ("sequence_unfitted", s0)]:
        model_rows = predictions.loc[predictions.model.eq(model)]
        expected = pd.Series(values, index=ids).reindex(model_rows.variant_id)
        np.testing.assert_allclose(expected, model_rows.prediction, atol=2e-13)
    comparisons = []
    for target in selected:
        heldout = int(np.flatnonzero(ids == target)[0])
        d = means(normalized_weights(variants, geometry, raw, excluded=heldout), y)
        seq_without = sequence_raw.copy()
        seq_without[:, heldout] = 0
        s = means(seq_without / seq_without.sum(axis=1, keepdims=True), y)
        for cohort in predictions.cohort.unique():
            rows = predictions.loc[predictions.cohort.eq(cohort)]
            if target not in set(rows.variant_id):
                continue
            target_set = set(rows.variant_id)
            training = np.flatnonzero(
                np.isin(ids, list(target_set)) & (ids != target) & np.isfinite(d + s)
            )
            features = {
                "intercept_only": (np.empty((len(training), 0)), np.empty((1, 0))),
                "density_fit": (d[training, None], d0[[heldout], None]),
                "sequence_fit": (s[training, None], s0[[heldout], None]),
            }
            if cohort == "am_common":
                am = variants.am.to_numpy()
                features.update(
                    {
                        "am_fit": (am[training, None], am[[heldout], None]),
                        "am_plus_density": (
                            np.column_stack([am[training], d[training]]),
                            np.array([[am[heldout], d0[heldout]]]),
                        ),
                        "am_plus_sequence": (
                            np.column_stack([am[training], s[training]]),
                            np.array([[am[heldout], s0[heldout]]]),
                        ),
                    }
                )
            for model, (train_x, test_x) in features.items():
                actual = float(
                    fitter.fractional_logistic(train_x, y[training], test_x)[0]
                )
                record = rows.loc[
                    rows.variant_id.eq(target) & rows.model.eq(model)
                ].iloc[0]
                np.testing.assert_allclose(actual, record.prediction, atol=1e-10)
                assert len(training) == record.training_targets
                comparisons.append(
                    dict(
                        variant_id=target,
                        cohort=cohort,
                        model=model,
                        training_targets=len(training),
                        absolute_prediction_error=abs(actual - record.prediction),
                    )
                )
    duplicates = variants.loc[variants.duplicated("protein_key", keep=False)]
    mixed = primary.loc[primary.density_source.eq("mixed")]
    # Derive this artifact even when the user-facing tables correctly suppress
    # the misleading intercept-only ranking diagnostic.
    intercept_rho = [
        dict(
            cohort=cohort,
            mechanical_loo_spearman=float(
                spearmanr(group.prediction, group.observed_fraction).statistic
            ),
        )
        for cohort, group in predictions.loc[
            predictions.model.eq("intercept_only")
        ].groupby("cohort")
    ]
    return dict(
        gene=gene,
        all_passed=True,
        variants=len(variants),
        full_gene_missense_prior_fixed=True,
        genomic_members_disjoint=True,
        identity_unit="frozen mixed genomic alleles and unresolved clinical protein aggregates; not a leave-one-protein-substitution-out experiment",
        distinct_protein_substitutions_with_multiple_DNA_units=int(
            duplicates.protein_key.nunique()
        ),
        density_source_counts=primary.density_source.value_counts().to_dict(),
        maximum_saved_weight_error=float(np.max(abs(saved - weights))),
        all_self_weights_zero=True,
        unsupported_preserved_as_NaN=int(np.isnan(d0).sum()),
        positive_normalized_pairs=int((weights > 0).sum()),
        same_residue_distinct_variant_pairs=int(
            ((weights > 0) & (positions[:, None] == positions[None, :])).sum()
        ),
        context_policy="Nearest donor copy per target context, normalize separately in each supported context, then equal context average",
        independent_outer_LOO_comparisons=comparisons,
        mixed_variants=mixed.reset_index()[
            [
                "variant_id",
                "canonical_pos",
                "context_density_min",
                "context_density_max",
                "context_spread",
            ]
        ].to_dict("records"),
        intercept_only_spearman=intercept_rho,
        geometry_sha256=hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256(
            (RUN / "run_structure.py").read_bytes()
        ).hexdigest(),
        summary_sha256=hashlib.sha256((RUN / "summarize.py").read_bytes()).hexdigest(),
    )


if __name__ == "__main__":
    results = [audit_gene(gene) for gene in ["HNF1A", "KCNQ1"]]
    (HERE / "independent_runner_checks.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n"
    )
    for result in results:
        print(
            result["gene"],
            result["density_source_counts"],
            "all independent checks passed",
        )
