"""Independently audit held-out donor removal and common training eligibility.

Reads frozen outputs; writes only analysis/fold_eligibility_audit.json. PPA is
read-only. The reference oracle explicitly reconstructs every newly unsupported
target plus a deterministic set of retained targets for each affected fold.
"""

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PPA = HERE.parents[3] / "ProteinProximityAnalysis/src"
sys.path.insert(0, str(PPA))

from alphafold_rin.empirical_density import (  # noqa: E402
    empirical_variant_density,
    sigmoid_kernel,
)


def parts(folder, stem):
    return pd.concat(
        [pd.read_csv(p) for p in sorted(folder.glob(f"{stem}.part*.csv.gz"))],
        ignore_index=True,
    )


def audit_gene(gene, config):
    folder = HERE / "analysis" / gene
    variants = parts(folder, "input_variants")
    primary = parts(folder, "primary_density")
    predictions = parts(folder, "loo_predictions")
    supported = primary.loc[primary.density.notna()].reset_index(drop=True)
    ids = variants.variant_id.tolist()
    lookup = {key: i for i, key in enumerate(ids)}
    target_ids = supported.variant_id.tolist()
    target_lookup = {key: i for i, key in enumerate(target_ids)}
    y = variants.posterior_mean.to_numpy()
    positions = variants.canonical_pos.to_numpy()
    raw = sigmoid_kernel(
        3.8 * np.sqrt(abs(supported.canonical_pos.to_numpy()[:, None] - positions))
    )
    for row, key in enumerate(target_ids):
        raw[row, lookup[key]] = 0
    sequence = raw / raw.sum(axis=1, keepdims=True)
    base_sequence = sequence @ y
    new_schema = "training_targets_lost_support" in predictions
    singleton_ids = supported.loc[supported.donor_count.eq(1), "variant_id"].tolist()
    weights = None
    if singleton_ids:
        weights = parts(folder, "primary_weights")
    cohorts, affected = [], []
    for cohort, rows in predictions.groupby("cohort"):
        heldout_ids = rows.variant_id.drop_duplicates().tolist()
        heldout_set = set(heldout_ids)
        target_rows = np.array([target_lookup[key] for key in heldout_ids])
        donor_cols = np.array([lookup[key] for key in heldout_ids])
        sub = sequence[np.ix_(target_rows, donor_cols)]
        max_location = np.unravel_index(sub.argmax(), sub.shape)
        min_remaining = float(1 - sub[max_location])
        assert min_remaining > 0.01, "Dominant-deletion cancellation needs review"
        direct_row = sequence[target_rows[max_location[0]]].copy()
        direct_row[donor_cols[max_location[1]]] = 0
        min_direct = float(direct_row.sum())
        np.testing.assert_allclose(min_remaining, min_direct, atol=1e-14)
        proven_losses = {}
        for target in singleton_ids:
            donor = weights.loc[weights.target_id.eq(target), "donor_id"].tolist()
            assert len(donor) == 1 and donor[0] != target
            if target in heldout_set and donor[0] in heldout_set:
                proven_losses.setdefault(donor[0], []).append(target)
        for heldout, fold in rows.groupby("variant_id"):
            losses = proven_losses.get(heldout, [])
            assert fold.training_targets.nunique() == 1
            assert int(fold.training_targets.iloc[0]) == len(heldout_ids) - 1 - len(
                losses
            )
            if new_schema:
                assert fold.training_targets_lost_support.nunique() == 1
                assert int(fold.training_targets_lost_support.iloc[0]) == len(losses)
            if losses:
                affected.append((cohort, heldout, losses, heldout_ids, fold))
        cohorts.append(
            {
                "cohort": cohort,
                "heldout_variants": len(heldout_ids),
                "models": sorted(rows.model.unique()),
                "same_training_target_count_across_models": True,
                "folds_losing_training_support": len(proven_losses),
                "sequence_min_remaining_normalized_weight_all_folds": min_remaining,
                "sequence_direct_sum_at_minimum": min_direct,
                "sequence_minimum_target": heldout_ids[max_location[0]],
                "sequence_minimum_heldout_donor": heldout_ids[max_location[1]],
            }
        )
    details = []
    if affected:
        frame = next(
            f for f in config["frames"] if config["primary"] == f["name"] + "_com_h3"
        )
        geometry = pd.read_csv(HERE / frame["path"]).fillna({"idr_segment": ""})
        models = []
        for start in range(0, len(target_ids), 128):
            result = empirical_variant_density(
                variants,
                geometry,
                target_ids=target_ids[start : start + 128],
                include_context_weights=False,
                include_context_model=True,
            )
            models.append(result.context_model)
        rebuilt = pd.concat([model.density() for model in models]).reindex(target_ids)
        np.testing.assert_allclose(rebuilt, supported.density, atol=1e-12)
        rebuilt_exclusions = {}
        for heldout in sorted({item[1] for item in affected}):
            excluded = pd.concat(
                [model.density(excluded_variant_ids=[heldout]) for model in models]
            ).reindex(target_ids)
            newly_missing = excluded.index[excluded.isna()].tolist()
            check_ids = list(
                dict.fromkeys(
                    newly_missing
                    + [heldout]
                    + target_ids[:: max(1, len(target_ids) // 12)]
                )
            )
            reference = (
                empirical_variant_density(
                    variants,
                    geometry,
                    target_ids=check_ids,
                    excluded_variant_ids=[heldout],
                    include_context_weights=False,
                    backend="reference",
                )
                .summary.set_index("variant_id")
                .density.reindex(check_ids)
            )
            np.testing.assert_allclose(
                excluded.reindex(check_ids), reference, atol=1e-11, equal_nan=True
            )
            rebuilt_exclusions[heldout] = (excluded, check_ids)
        for cohort, heldout, losses, heldout_ids, fold in affected:
            excluded, check_ids = rebuilt_exclusions[heldout]
            training = [key for key in heldout_ids if key != heldout]
            actual_losses = [key for key in training if not np.isfinite(excluded[key])]
            assert sorted(actual_losses) == sorted(losses)
            train_indices = [target_lookup[key] for key in training]
            col = lookup[heldout]
            remaining = 1 - sequence[:, col]
            updated = (base_sequence - sequence[:, col] * y[col]) / remaining
            direct = sequence.copy()
            direct[:, col] = 0
            direct /= direct.sum(axis=1, keepdims=True)
            direct_mean = direct @ y
            max_error = float(
                np.max(abs(updated[train_indices] - direct_mean[train_indices]))
            )
            assert max_error < 1e-12
            retained = [key for key in training if key not in losses]
            assert len(retained) == int(fold.training_targets.iloc[0])
            lost_positions = variants.loc[
                variants.variant_id.isin(losses), "canonical_pos"
            ].unique()
            context = geometry.loc[
                geometry.canonical_pos.isin(lost_positions),
                [
                    "frame_id",
                    "chain",
                    "canonical_pos",
                    "geometry_state",
                    "idr_segment",
                    "plddt",
                ],
            ].to_dict("records")
            names = variants.set_index("variant_id").protein_key
            details.append(
                {
                    "cohort": cohort,
                    "heldout_id": heldout,
                    "heldout_protein_key": names[heldout],
                    "lost_training_ids": losses,
                    "lost_training_protein_keys": [names[key] for key in losses],
                    "lost_target_contexts": context,
                    "training_targets": len(retained),
                    "training_ids_sha256_sorted_newline": hashlib.sha256(
                        ("\n".join(sorted(retained)) + "\n").encode()
                    ).hexdigest(),
                    "training_target_selection_shared_before_all_model_fits": True,
                    "context_model_all_training_targets_reconstructed": True,
                    "reference_target_ids": check_ids,
                    "reference_oracle_agrees_including_lost_targets": True,
                    "sequence_min_remaining_weight_training_targets": float(
                        remaining[train_indices].min()
                    ),
                    "sequence_exclusion_max_error_vs_direct_positive_renormalization": max_error,
                }
            )
    return {
        "gene": gene,
        "loss_column_present": new_schema,
        "supported_singleton_donor_targets": singleton_ids,
        "zero_loss_proof_when_no_singletons": "At least two distinct eligible donors for every supported target: deleting one leaves support. Counts additionally equal cohort size minus one in every saved model row.",
        "cohorts": cohorts,
        "affected_folds": details,
    }


def main():
    config = json.loads((HERE / "geometry_config.json").read_text())
    report = {
        "scope": "Independent audit of frozen outputs; does not refit models or alter outcomes.",
        "reason": "LDLR A234T/A234V each have only the other as an eligible donor in a one-residue candidate IDR. Outer LOO removes that donor globally; both model comparisons must omit the newly unsupported training target.",
        "ppa_source_sha256": {
            name: hashlib.sha256(
                (PPA / "alphafold_rin" / name).read_bytes()
            ).hexdigest()
            for name in ["empirical_density.py", "empirical_context.py"]
        },
        "runner_sha256": hashlib.sha256(
            (HERE / "run_structure.py").read_bytes()
        ).hexdigest(),
        "genes": [
            audit_gene(gene, config[gene])
            for gene in ["HNF1A", "LDLR", "KCNQ1", "BRCA2"]
        ],
    }
    target = HERE / "analysis/fold_eligibility_audit.json"
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(target)


if __name__ == "__main__":
    main()
