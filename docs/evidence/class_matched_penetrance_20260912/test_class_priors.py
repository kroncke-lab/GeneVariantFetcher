"""Class separation and count-update invariants for empirical priors."""

import numpy as np
import pandas as pd
import pytest

from fit_class_priors import GENES, fit_types, select_type


def fixture():
    rows = []
    for gene in GENES:
        for kind, pairs in [
            ("missense", [(0, 1), (1, 1), (1, 2)]),
            ("nonsense", [(0, 1), (5, 1), (7, 1)]),
        ]:
            for i, (a, u) in enumerate(pairs):
                rows.append(
                    {
                        "gene": gene,
                        "unit_id": f"{gene}-{kind}-{i}",
                        "vclass": kind
                        if i != 1 or kind != "nonsense"
                        else "stop_gained",
                        "canonical_wt_status": "match",
                        "origin": "population_only"
                        if a == 0
                        else "literature_and_population",
                        "affected": a,
                        "unaffected_literature": 0,
                        "gnomad_carriers": u,
                        "unaffected": u,
                        "n": a + u,
                    }
                )
    return pd.DataFrame(rows)


def primary(frame):
    table = fit_types(frame)[0]
    return table.loc[
        table.scope.isin(["canonical_missense", "canonical_nonsense"])
    ].reset_index(drop=True)


def test_other_variant_classes_cannot_change_either_prior():
    base = fixture()
    extra = pd.concat(
        [
            base.assign(
                vclass=kind,
                unit_id=base.unit_id + "-" + kind,
                affected=10000,
                n=10000 + base.unaffected,
            )
            for kind in [
                "synonymous",
                "frameshift",
                "splice",
                "noncoding",
                "stop_lost",
                "start_lost",
            ]
        ],
        ignore_index=True,
    )
    pd.testing.assert_frame_equal(
        primary(base), primary(pd.concat([base, extra], ignore_index=True))
    )


def test_nonsense_counts_cannot_change_missense_prior():
    base = fixture()
    changed = base.copy()
    mask = changed.vclass.isin(["nonsense", "stop_gained"])
    changed.loc[mask, "affected"] += 5
    changed.loc[mask, "n"] += 5
    before, after = primary(base), primary(changed)
    pd.testing.assert_frame_equal(
        before.loc[before.variant_type.eq("missense")],
        after.loc[after.variant_type.eq("missense")],
    )
    assert not np.allclose(
        before.loc[before.variant_type.eq("nonsense"), "mean"],
        after.loc[after.variant_type.eq("nonsense"), "mean"],
    )


def test_equivalent_stop_gained_label_is_same_type():
    base = fixture()
    labels = base.copy()
    labels["vclass"] = labels.vclass.replace({"stop_gained": "nonsense"})
    pd.testing.assert_frame_equal(primary(base), primary(labels))


def test_unaffected_singleton_updates_beta_only_and_preserves_nonzero_mean():
    moments, post, _, _ = fit_types(fixture())
    singleton = post.loc[post.affected.eq(0) & post.unaffected.eq(1)]
    np.testing.assert_allclose(singleton.posterior_alpha, singleton.alpha_empirical)
    np.testing.assert_allclose(singleton.posterior_beta, singleton.beta_empirical + 1)
    mean = singleton.alpha_empirical / (
        singleton.alpha_empirical + singleton.beta_empirical
    )
    assert ((singleton.posterior_mean > 0) & (singleton.posterior_mean < mean)).all()
    for r in moments.itertuples():
        assert r.unaffected_singleton_posterior == pytest.approx(
            r.mean * r.strength / (r.strength + 1)
        )


def test_noncanonical_identity_never_enters_matching_prior():
    base = fixture()
    bad = base.assign(canonical_wt_status="mismatch", unit_id=base.unit_id + "-bad")
    pd.testing.assert_frame_equal(
        primary(base), primary(pd.concat([base, bad], ignore_index=True))
    )


def test_unknown_type_fails_without_fallback():
    with pytest.raises(ValueError, match="Unsupported prior type"):
        select_type(fixture(), "truncating")
