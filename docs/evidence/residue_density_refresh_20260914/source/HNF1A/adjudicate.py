#!/usr/bin/env python3
"""Apply bounded HNF1A endpoint/somatic adjudications to frozen clinical keys."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
OTHER = HERE.parent / "other_genes"
FLAGS = (
    REPO
    / "docs/evidence/population_inclusive_penetrance_20260912/audit/literature_identity_flags.csv.gz"
)
OBS = (
    REPO.parent
    / "BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results/HNF1A_protocol/observations.csv"
)

# Values are source-based counts for the specified endpoint at the reported
# follow-up; None means that the observation is quarantined, not observed zero.
PLAN = [
    (
        24915262,
        "E508K",
        None,
        None,
        "quarantine_T2_endpoint",
        "Case-control carriers with type 2 diabetes; no established MODY endpoint for affected cases or MODY-negative status of controls.",
        ["52 affected carriers", "type 2 diabetes"],
    ),
    (
        28116330,
        "I27L",
        None,
        None,
        "quarantine_T2_endpoint",
        "38 diabetic and 23 control carriers are explicitly from a type-2-diabetes study; not established MODY affected/unaffected counts.",
        ["38/33", "23/11", "T2D"],
    ),
    (
        20172480,
        "I27L",
        None,
        None,
        "quarantine_unresolved_variant_and_endpoint_assignment",
        "Abstract-only source describes 31 total subjects split into MODY3=10, MODY2=15 and type2=6 with I27L; ownership of all 31 as I27L carriers and HNF1A-MODY endpoint is unestablished. No unsupported 31-to-6 MODY replacement.",
        ["31 previously diagnosed", "6 subjects"],
    ),
    (
        30121369,
        "L214Q",
        None,
        None,
        "quarantine_somatic_tumor",
        "Somatic HNF1A alteration in hepatocellular neoplasm; not an affected germline MODY carrier.",
        ["Two somatic mutations", "L214Q"],
    ),
    (
        30121369,
        "E32X",
        None,
        None,
        "quarantine_somatic_tumor",
        "Somatic HNF1A alteration in hepatocellular neoplasm; not an affected germline MODY carrier.",
        ["Two somatic mutations", "E32*"],
    ),
    (
        29101032,
        "Q511L",
        None,
        None,
        "quarantine_somatic_tumor",
        "Variant detected only in tumor and not adjacent non-tumorous liver; not an affected germline MODY carrier.",
        ["Q511L", "not in the adjacent"],
    ),
    (
        38133737,
        "K205E",
        None,
        None,
        "quarantine_somatic_tumor",
        "HNF1A is wild type in healthy liver and no pathogenic germline HNF1A variant identified; tumor variant does not establish a germline MODY carrier.",
        ["K205E", "healthy liver tissue was wild type"],
    ),
    (
        26631547,
        "S247T",
        None,
        None,
        "quarantine_somatic_tumor",
        "Somatic HCC variant studied at DNA/RNA/protein levels; an expressed tumor allele is not a germline MODY carrier.",
        ["S247T", "somatic mutations"],
    ),
    (
        31483937,
        "R272H",
        4,
        1,
        "replace_adenoma_endpoint_with_reported_diabetes",
        "Five confirmed germline carriers: index, mother, one sister and brother developed diabetes; the other sister has no diabetes onset at follow-up. Frozen 5/0 corresponds to adenoma burden, not reported diabetes 4/1.",
        ["DM was diagnosed 8 years later", "she does not have symptoms of DM"],
    ),
    (
        14598263,
        "R229X",
        0,
        3,
        "replace_adenoma_endpoint_with_reported_diabetes",
        "Family A has three confirmed germline carriers (A1,A2,A3), all with normal fasting glucose and no reported diabetes. Two had liver adenoma; frozen 2/1 is the wrong phenotype partition.",
        ["R229X", "fasting glucose levels"],
    ),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame, name):
    kwargs = dict(index=False, lineterminator="\n", float_format="%.12g")
    if name.endswith(".gz"):
        kwargs["compression"] = dict(method="gzip", mtime=0)
    frame.to_csv(HERE / name, **kwargs)
    assert (HERE / name).stat().st_size < 1_200_000


def main():
    baseline_all = pd.read_csv(FLAGS)
    baseline = baseline_all.loc[baseline_all.gene.eq("HNF1A")].copy()
    obs = pd.read_csv(OBS)
    metadata = pd.read_csv(OTHER / "source_metadata_screen.csv.gz")
    metadata = metadata.loc[metadata.gene.eq("HNF1A")].copy()
    metadata["pmid"] = pd.to_numeric(metadata.pmid)
    prior_checks = json.loads(
        (
            REPO
            / "docs/evidence/brca2_distance_prior_audit_20260913/source/cross_gene_checks.json"
        ).read_text()
    )
    assert sha(OBS) == prior_checks["original_observation_hashes"]["HNF1A"]
    counts = obs.groupby("key")[["affected", "unaffected"]].sum()
    np.testing.assert_allclose(
        baseline.affected_literature, baseline.key.map(counts.affected).fillna(0)
    )
    np.testing.assert_allclose(
        baseline.unaffected_literature, baseline.key.map(counts.unaffected).fillna(0)
    )
    decisions, source_hashes = [], {}
    for pmid, key, a, u, action, reason, required in PLAN:
        selected = obs.loc[obs.pmid.eq(pmid) & obs.key.eq(key)]
        assert len(selected) == 1
        r = selected.iloc[0]
        meta = metadata.loc[
            metadata.pmid.eq(pmid)
            & metadata.key.eq(key)
            & metadata.variant_id.eq(r.variant_id)
        ].iloc[0]
        source = Path(meta.db_path).parent / "pmc_fulltext" / f"{pmid}_CLEANED.md"
        content = source.read_text()
        for term in required:
            assert term.lower() in content.lower(), (pmid, key, term)
        source_hashes[str(source.relative_to(REPO))] = sha(source)
        lines = content.splitlines()
        matched_lines = [
            str(i)
            for i, line in enumerate(lines, 1)
            if any(term.lower() in line.lower() for term in required)
        ]
        decisions.append(
            dict(
                gene="HNF1A",
                key=key,
                pmid=pmid,
                variant_id=int(r.variant_id),
                vclass=r.vclass,
                original_affected=float(r.affected),
                original_unaffected=float(r.unaffected),
                A_to_remove=float(r.affected),
                U_to_remove=float(r.unaffected),
                A_to_restore=0 if a is None else a,
                U_to_restore=0 if u is None else u,
                adjudicated_affected=a,
                adjudicated_unaffected=u,
                action=action,
                reason=reason,
                source_status="abstract_only"
                if "ABSTRACT-ONLY" in content
                else "full_primary_text",
                source_path=str(source.relative_to(REPO)),
                source_lines=";".join(matched_lines),
                source_location=meta.source_location,
                original_key_quotes=meta.key_quotes,
                primary_url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                unaffected_semantics="No diabetes at reported follow-up; not lifetime absence"
                if u is not None
                else "unknown_not_zero",
            )
        )
    decisions = pd.DataFrame(decisions)
    save(decisions, "observation_decisions.csv")
    selected_ids = set(zip(decisions.key, decisions.pmid, decisions.variant_id))
    original_ids = list(zip(obs.key, obs.pmid, obs.variant_id))
    retained = obs.loc[[i not in selected_ids for i in original_ids]].copy()
    replacements = decisions.loc[decisions.adjudicated_affected.notna()].rename(
        columns={
            "adjudicated_affected": "affected",
            "adjudicated_unaffected": "unaffected",
        }
    )
    complete = pd.concat(
        [
            retained[["key", "pmid", "variant_id", "affected", "unaffected"]],
            replacements[["key", "pmid", "variant_id", "affected", "unaffected"]],
        ],
        ignore_index=True,
    )
    new_counts = complete.groupby("key")[["affected", "unaffected"]].sum()
    updated = baseline.copy()
    updated["refresh_baseline_affected_literature"] = updated.affected_literature
    updated["refresh_baseline_unaffected_literature"] = updated.unaffected_literature
    updated["affected_literature"] = updated.key.map(new_counts.affected).fillna(0)
    updated["unaffected_literature"] = updated.key.map(new_counts.unaffected).fillna(0)
    updated["n_literature"] = (
        updated.affected_literature + updated.unaffected_literature
    )
    updated["drop_clinical_key_zero_clinical_evidence"] = updated.n_literature.eq(0)
    pmids = complete.groupby("key").pmid.agg(
        lambda x: ";".join(map(str, sorted(set(x))))
    )
    updated["pmids"] = updated.key.map(pmids).fillna("")
    delta = decisions.groupby("key")[
        ["A_to_remove", "U_to_remove", "A_to_restore", "U_to_restore"]
    ].sum()
    for col in delta:
        updated[f"refresh_{col}"] = updated.key.map(delta[col]).fillna(0)
    np.testing.assert_allclose(
        updated.affected_literature,
        updated.refresh_baseline_affected_literature
        - updated.refresh_A_to_remove
        + updated.refresh_A_to_restore,
    )
    np.testing.assert_allclose(
        updated.unaffected_literature,
        updated.refresh_baseline_unaffected_literature
        - updated.refresh_U_to_remove
        + updated.refresh_U_to_restore,
    )
    assert (
        updated.affected_literature.ge(0).all()
        and updated.unaffected_literature.ge(0).all()
    )
    assert (
        updated.affected_literature.sum() == 633
        and updated.unaffected_literature.sum() == 83
    )
    save(updated, "clinical_input.csv.gz")
    save(updated.loc[updated.key.isin(delta.index)], "clinical_key_deltas.csv")
    save(complete, "corrected_observations.csv.gz")
    # Optional broader-diabetes sensitivity restores the two established T2
    # case/control observations, while retaining all other source repairs.
    broad = updated[list(baseline.columns)].copy()
    for pmid, key in [(24915262, "E508K"), (28116330, "I27L")]:
        r = decisions.loc[decisions.pmid.eq(pmid) & decisions.key.eq(key)].iloc[0]
        idx = broad.index[broad.key.eq(key)][0]
        broad.loc[idx, "affected_literature"] += r.original_affected
        broad.loc[idx, "unaffected_literature"] += r.original_unaffected
        broad.loc[idx, "pmids"] = ";".join(
            sorted(set(str(broad.loc[idx, "pmids"]).split(";")) | {str(pmid)})
        ).strip(";")
    broad["n_literature"] = broad.affected_literature + broad.unaffected_literature
    broad["drop_clinical_key_zero_clinical_evidence"] = broad.n_literature.eq(0)
    save(broad, "broader_diabetes_sensitivity_clinical_input.csv.gz")
    checks = dict(
        passed=True,
        original_observations=len(obs),
        original_A=float(obs.affected.sum()),
        original_U=float(obs.unaffected.sum()),
        adjudicated_observations=len(decisions),
        old_adjudicated_A=float(decisions.original_affected.sum()),
        old_adjudicated_U=float(decisions.original_unaffected.sum()),
        primary_A=float(updated.affected_literature.sum()),
        primary_U=float(updated.unaffected_literature.sum()),
        restored_source_A=int(decisions.A_to_restore.sum()),
        restored_source_U=int(decisions.U_to_restore.sum()),
        clinical_keys=len(updated),
        retained_clinical_keys=int(updated.n_literature.gt(0).sum()),
        baseline_hashes={str(FLAGS.relative_to(REPO)): sha(FLAGS), str(OBS): sha(OBS)},
        source_hashes=source_hashes,
        output_ledger_sha256=sha(HERE / "clinical_input.csv.gz"),
        scope="Top three original affected sources plus already-flagged somatic and adenoma endpoint examples; no claim all remaining sources are MODY-specific. Frameshift/splice flags outside requested missense/nonsense fits remain unchanged.",
        population_rule="Drop zero-clinical-evidence keys and rebuild/rejoin population union; retain all observed gnomAD alleles.",
    )
    (HERE / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in checks.items()
                if k not in ["source_hashes", "baseline_hashes"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
