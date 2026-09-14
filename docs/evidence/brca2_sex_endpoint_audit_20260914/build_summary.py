"""Reproduce denominator sensitivities and record bounded source adjudications.

This produces audit companions. It does not refit priors or replace the frozen
residue maps, and the denominator sensitivities are not women-only estimates.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE = HERE.parent / "residue_density_refresh_20260914"


def main():
    hashes = {}

    def record(path):
        hashes[str(path.relative_to(REPO))] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        return path

    def units(gene):
        return pd.concat(
            [
                pd.read_csv(record(p))
                for p in sorted(
                    (BASE / f"analysis/{gene}").glob("missense_posteriors.part*.csv.gz")
                )
            ],
            ignore_index=True,
        )

    priors = pd.read_csv(record(BASE / "tables/prior_summary.csv"))
    priors = priors.loc[
        priors.scenario.eq("refreshed")
        & priors.input_scope.eq("primary")
        & priors.variant_type.eq("missense")
    ].set_index("gene", verify_integrity=True)
    brca = units("BRCA2").set_index("key", verify_integrity=True)
    a, b = priors.loc["BRCA2", ["alpha_empirical", "beta_empirical"]]
    probe = json.loads(record(HERE / "sex_probe.json").read_text())
    keys = {
        "13-32363389-G-T": "K2729N",
        "13-32363369-G-C": "D2723H",
        "13-32380043-C-T": "R3052W",
        "13-32362595-G-C": "W2626C",
        "13-32332592-A-C": "N372H",
        "13-32363225-A-G": "I2675V",
    }
    rows = []
    for receipt in probe:
        raw = record(REPO / receipt["source"])
        assert hashes[str(raw.relative_to(REPO))] == receipt["source_sha256"]
        key = keys[receipt["variant_id"]]
        row = brca.loc[key]
        assert row.member_alleles == receipt["variant_id"]
        joint = receipt["assays"]["joint"]
        carriers = joint["ac"] - joint["homozygote_count"]
        sexes = {s["id"]: s for s in joint["sex"]}
        xx = sexes["XX"]["ac"] - sexes["XX"]["homozygote_count"]
        xy = sexes["XY"]["ac"] - sexes["XY"]["homozygote_count"]
        assert xx + xy == carriers == row.gnomad_carriers
        expected = (a + row.affected) / (
            a + b + row.affected + row.unaffected_literature + carriers
        )
        np.testing.assert_allclose(expected, row.posterior_mean)
        rows.append(
            dict(
                key=key,
                allele=receipt["variant_id"],
                affected_all_sex=row.affected,
                literature_unaffected_all_sex=row.unaffected_literature,
                gnomad_all_sex_carriers=carriers,
                gnomad_XX_carriers=xx,
                gnomad_XY_carriers=xy,
                frozen_all_sex_own_posterior=expected,
                denominator_only_counterfactual_own_posterior=(a + row.affected)
                / (a + b + row.affected + row.unaffected_literature + xx),
                clinically_sex_matched=False,
                prior_refitted=False,
            )
        )
    pd.DataFrame(rows).to_csv(
        HERE / "sex_denominator_sensitivity.csv", index=False, lineterminator="\n"
    )

    diagnostics = []
    for gene, prior in priors.iterrows():
        x = units(gene)
        n = x.affected + x.unaffected
        assert n.ge(1).all()
        y = x.affected / n
        w = 1 - 1 / (n + 0.01)
        mu = np.average(y, weights=w)
        var = np.sum(w * (y - mu) ** 2) / len(x)
        strength = mu * (1 - mu) / var - 1
        normalized_var = np.average((y - mu) ** 2, weights=w)
        normalized_strength = mu * (1 - mu) / normalized_var - 1
        np.testing.assert_allclose([mu, strength], [prior["mean"], prior.strength])
        diagnostics.append(
            dict(
                gene=gene,
                units=len(x),
                singletons=int(n.eq(1).sum()),
                unaffected_singletons=int((n.eq(1) & x.affected.eq(0)).sum()),
                singleton_unit_share=float(n.eq(1).mean()),
                singleton_fit_weight_share=float(w[n.eq(1)].sum() / w.sum()),
                historical_prior_mean=float(mu),
                historical_strength=float(strength),
                historical_unaffected_singleton_posterior=float(
                    mu * strength / (strength + 1)
                ),
                normalized_mse_strength=float(normalized_strength),
                normalized_mse_unaffected_singleton_posterior=float(
                    mu * normalized_strength / (normalized_strength + 1)
                ),
                normalized_mse_is_only_a_diagnostic=True,
            )
        )
    (HERE / "prior_weight_diagnostic.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n"
    )

    source = (
        REPO
        / "results/grant_e2e_20260909/shards/BRCA2_2/BRCA2/20260909_115043/pmc_fulltext/17100994_CLEANED.md"
    )
    content = record(source).read_text()
    assert "| 16 | 4862A>C | E1581D | MS | 48 | 3 | 1 |" in content
    assert "| 27 | 10462A>G | I3412V | MS | 26 | 5 | 110 |" in content
    male = (
        REPO
        / "results/grant_e2e_20260909/shards/BRCA2_1/BRCA2/20260909_115040/pmc_fulltext/20927582_CLEANED.md"
    )
    male_text = record(male).read_text()
    assert "D2723H in two unrelated Caucasians cases" in male_text
    observations = pd.read_csv(
        record(HERE / "retained_missense_observation_scope.csv.gz")
    )
    decisions = [
        dict(
            key="E1581D",
            pmid=17100994,
            observation_variant_id=40000141,
            A_to_remove=48,
            U_to_add=0,
            scope="all_BRCA2_endpoints",
            reason="Wrong gene: T1 is the BRCA1 table; its 48 cases and 3 controls cannot be BRCA2 observations.",
            source_path=str(source.relative_to(REPO)),
            source_lines="65-107; E1581D at 85",
        ),
        dict(
            key="I3412V",
            pmid=17100994,
            observation_variant_id=40000163,
            A_to_remove=0,
            U_to_add=5,
            scope="BRCA2_breast_cancer_source_table",
            reason="T2 contains 26 patient observations and 5 normal-control observations; 110 is BIC entries, not people.",
            source_path=str(source.relative_to(REPO)),
            source_lines="47; 112-146; I3412V at 145",
        ),
        dict(
            key="D2723H",
            pmid=20927582,
            observation_variant_id=20001234,
            A_to_remove=2,
            U_to_add=0,
            scope="female_breast_cancer_only",
            reason="Two confirmed affected male carriers; retain in male breast-cancer endpoint and exclude from a female numerator.",
            source_path=str(male.relative_to(REPO)),
            source_lines="31; 88",
        ),
    ]
    for d in decisions:
        match = observations.loc[
            observations.key.eq(d["key"])
            & observations.pmid.eq(d["pmid"])
            & observations.observation_variant_id.eq(d["observation_variant_id"])
        ]
        assert len(match) == 1 and match.affected.iloc[0] >= d["A_to_remove"]
        d["baseline_affected"] = float(match.affected.iloc[0])
        d["baseline_unaffected"] = float(match.unaffected.iloc[0])
        d["applied_to_frozen_model"] = False
    pd.DataFrame(decisions).to_csv(
        HERE / "reviewed_source_patch.csv", index=False, lineterminator="\n"
    )
    record(Path(__file__).resolve())
    (HERE / "summary_input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    print(pd.DataFrame(rows).to_string(index=False))
    print(
        "Validated six exact allele counts, five prior diagnostics, three source decisions."
    )


if __name__ == "__main__":
    main()
