"""All single-donor exclusions from a frozen PPA ContextDensityModel.

This helper accelerates the outer variant-LOO feature calculation only. It does
not fit priors, modify posterior values, or change donor eligibility or geometry.
Each column removes one canonical donor globally, in addition to the target's
own identity already removed by PPA. Other alleles at that residue remain.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def all_excluded_density(model, *, donor_batch_size: int = 512) -> pd.DataFrame:
    """Return target-by-heldout-donor densities with exact context normalization.

    Rows are ``model.target_ids`` and columns are ``model.donor_ids`` in their
    existing order. Donor aliases have already been consolidated by PPA; target
    aliases remain separate rows. Every supported target context is normalized
    after exclusion and remaining supported contexts receive equal weight.

    Donor blocks bound the temporary context-by-donor arrays; the returned dense
    matrix requires 8 * targets * donors bytes. Call on modest target batches
    (for example 128), then store the result in the caller's ordered memmap.
    A column with any context weight >= .99 uses PPA's existing stable single-
    exclusion method, including log-space recovery and lost-context handling.
    Neither the model nor its arrays are mutated.
    """
    if (
        isinstance(donor_batch_size, (bool, np.bool_))
        or not isinstance(donor_batch_size, (int, np.integer))
        or donor_batch_size < 1
    ):
        raise ValueError("donor_batch_size must be a positive integer")
    target_count, donor_count = len(model.target_ids), len(model.donor_ids)
    answer = np.full((target_count, donor_count), np.nan)
    valid = np.isfinite(model.context_log_sums) & np.isfinite(model.context_densities)
    context_rows = model.context_target_rows[valid]
    geometry_rows = model.context_geometry_rows[valid]
    log_sums = model.context_log_sums[valid]
    densities = model.context_densities[valid]
    supported = np.bincount(context_rows, minlength=target_count)
    own_ids = model.target_canonical_ids[context_rows]
    lookup = {identity: column for column, identity in enumerate(model.donor_ids)}
    own_columns = np.array(
        [lookup.get(identity, -1) for identity in own_ids], dtype=int
    )
    for start in range(0, donor_count, donor_batch_size):
        stop = min(start + donor_batch_size, donor_count)
        width = stop - start
        # Advanced indexing creates a new block, never a view of cached logs.
        q = model.position_log_kernels[
            geometry_rows[:, None],
            model.donor_position_indices[None, start:stop],
        ]
        q -= log_sums[:, None]
        own = (own_columns >= start) & (own_columns < stop)
        q[np.flatnonzero(own), own_columns[own] - start] = -np.inf
        np.exp(q, out=q)
        dominant = np.any(q >= 0.99, axis=0)
        # Dominant columns are replaced below by PPA's stable implementation.
        # Zero their provisional q to avoid undefined 0/0 intermediate values.
        q[:, dominant] = 0
        values = (densities[:, None] - q * model.donor_means[None, start:stop]) / (
            1 - q
        )
        sums = np.zeros((target_count, width))
        np.add.at(sums, context_rows, values)
        np.divide(
            sums,
            supported[:, None],
            out=answer[:, start:stop],
            where=supported[:, None] > 0,
        )
        for local_column in np.flatnonzero(dominant):
            column = start + local_column
            fallback = model.density(excluded_variant_ids=[model.donor_ids[column]])
            answer[:, column] = fallback.reindex(model.target_ids).to_numpy()
    return pd.DataFrame(
        answer,
        index=pd.Index(model.target_ids, name="variant_id"),
        columns=pd.Index(model.donor_ids, name="excluded_donor_id"),
    )
