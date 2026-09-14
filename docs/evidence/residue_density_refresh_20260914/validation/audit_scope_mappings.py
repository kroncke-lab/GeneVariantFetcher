"""Trace source type labels versus canonical identities actually used by each fit."""

import json
from pathlib import Path

import pandas as pd

from audit_refresh import GEOMETRY, REFRESH, EVIDENCE, read, close, sha, clinical_check

HERE = Path(__file__).resolve().parent


def main():
    rows, summaries = [], []
    for gene in GEOMETRY:
        source, _ = clinical_check(gene)
        ledger = read(REFRESH / "analysis" / gene, "clinical_join_ledger").set_index(
            "key"
        )
        for kind in ["missense", "nonsense"]:
            prior_source = source.loc[
                source.vclass.isin(
                    [kind, "stop_gained"] if kind == "nonsense" else [kind]
                )
                & source.canonical_wt_status.eq("match")
                & (source.affected_literature + source.unaffected_literature).gt(0)
            ]
            used = read(REFRESH / "analysis" / gene, kind + "_posteriors")
            used = used.loc[used.literature_key.notna()].set_index("literature_key")
            omitted = set(prior_source.index) - set(used.index)
            resolved = set(used.index) - set(prior_source.index)
            a_difference = (
                source.loc[list(resolved), "affected_literature"].sum()
                - source.loc[list(omitted), "affected_literature"].sum()
            )
            u_difference = (
                source.loc[list(resolved), "unaffected_literature"].sum()
                - source.loc[list(omitted), "unaffected_literature"].sum()
            )
            close(
                [
                    prior_source.affected_literature.sum() + a_difference,
                    prior_source.unaffected_literature.sum() + u_difference,
                ],
                [used.affected.sum(), used.unaffected_literature.sum()],
            )
            summaries.append(
                dict(
                    gene=gene,
                    variant_type=kind,
                    source_literal_canonical_class_A=prior_source.affected_literature.sum(),
                    source_literal_canonical_class_U=prior_source.unaffected_literature.sum(),
                    union_fitted_clinical_A=used.affected.sum(),
                    union_fitted_clinical_U=used.unaffected_literature.sum(),
                    net_A_scope_change=a_difference,
                    net_U_scope_change=u_difference,
                )
            )
            for direction, keys in [
                ("source_label_not_in_fit", omitted),
                ("canonical_join_resolved_into_fit", resolved),
            ]:
                for key in sorted(keys):
                    s = source.loc[key]
                    reason = (
                        ledger.loc[key, "join_exclusion_reason"]
                        if key in ledger.index
                        else "no_remaining_clinical_carrier_evidence"
                    )
                    if direction == "canonical_join_resolved_into_fit":
                        reason = used.loc[key, "join_method"]
                    rows.append(
                        dict(
                            gene=gene,
                            variant_type=kind,
                            key=key,
                            direction=direction,
                            source_vclass=s.vclass,
                            source_canonical_wt_status=s.canonical_wt_status,
                            source_A=s.affected_literature,
                            source_U=s.unaffected_literature,
                            reason=reason,
                            resolved_protein_key=used.loc[key, "protein_key"]
                            if key in used.index
                            else "",
                        )
                    )
    # The KCNQ1 sensitivity keeps event splits but still removes total-carrier/count mistakes.
    original = pd.read_csv(
        EVIDENCE
        / "population_inclusive_penetrance_20260912/audit/literature_identity_flags.csv.gz",
        low_memory=False,
    )
    expected = original.loc[original.gene.eq("KCNQ1")].set_index("key")[
        ["affected_literature", "unaffected_literature"]
    ]
    decisions = pd.read_csv(REFRESH / "source/KCNQ1/observation_corrections.csv")
    for key, group in decisions.loc[decisions.remove_from_event_sensitivity].groupby(
        "key"
    ):
        expected.loc[key] += [
            group.affected_restored.sum() - group.affected_removed.sum(),
            group.unaffected_restored.sum() - group.unaffected_removed.sum(),
        ]
    events = pd.read_csv(
        REFRESH / "source/KCNQ1/clinical_input_cardiac_events.csv.gz", low_memory=False
    ).set_index("key")
    close(expected.reindex(events.index), events[expected.columns])
    pd.DataFrame(rows).to_csv(HERE / "clinical_type_scope_changes.csv", index=False)
    pd.DataFrame(summaries).to_csv(
        HERE / "clinical_type_scope_summary.csv", index=False
    )
    receipt = dict(
        all_type_scope_deltas_reproduced=True,
        kcnq1_event_sensitivity_source_arithmetic_verified=True,
        source_to_union_transition_rows=len(rows),
        audit_script_sha256=sha(__file__),
    )
    (HERE / "scope_mapping_audit.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
