"""Summarize the count, prior and distance diagnostics behind the residue plots."""

from pathlib import Path
import json

import pandas as pd


HERE = Path(__file__).resolve().parent
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]


def main():
    summaries, bandwidth = [], []
    for gene in GENES:
        folder = HERE / "analysis" / gene
        frame = pd.concat(
            [
                pd.read_csv(p, low_memory=False)
                for p in sorted(folder.glob("density_h3.part*.csv.gz"))
            ],
            ignore_index=True,
        )
        fit = pd.read_csv(folder / "prior_comparison.csv")
        prior = fit.loc[
            fit.scenario.eq("refreshed") & fit.variant_type.eq("missense")
        ].iloc[0]
        supported = frame.loc[frame.density.notna()]
        zero_donors = supported.loc[supported.affected_donor_count.eq(0)]
        prior_share = supported.prior_component / supported.density
        summaries.append(
            dict(
                gene=gene,
                variants=len(frame),
                supported=len(supported),
                zero_affected_units=int(frame.affected.eq(0).sum()),
                zero_affected_singletons=int(
                    (frame.affected.eq(0) & frame.n.eq(1)).sum()
                ),
                singleton_unaffected_posterior=prior.alpha_empirical
                / (prior.strength + 1),
                median_prior_retention=supported.neighborhood_prior_retention.median(),
                median_prior_component_share=prior_share.median(),
                median_weight_share_beyond_20=supported.weight_share_beyond_20.median(),
                max_weight_share_beyond_20=supported.weight_share_beyond_20.max(),
                no_affected_donor_variants=len(zero_donors),
                no_affected_donor_residues=zero_donors.canonical_pos.nunique(),
                zero_donor_primary_min=zero_donors.density.min(),
                zero_donor_primary_max=zero_donors.density.max(),
            )
        )
        bandwidth.append(pd.read_csv(folder / "bandwidth_sensitivity.csv"))
    pd.DataFrame(summaries).to_csv(
        HERE / "tables/sanity_summary.csv", index=False, lineterminator="\n"
    )
    pd.concat(bandwidth, ignore_index=True).to_csv(
        HERE / "tables/bandwidth_sensitivity.csv", index=False, lineterminator="\n"
    )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
