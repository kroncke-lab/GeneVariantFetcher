#!/usr/bin/env python3
"""Reproduce bounded LDLR FH endpoint corrections from exact frozen source rows."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
FLAGS = (
    REPO
    / "docs/evidence/population_inclusive_penetrance_20260912/audit/literature_identity_flags.csv.gz"
)
OBS = (
    REPO.parent
    / "BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results/LDLR_protocol/observations.csv"
)
META = HERE.parent / "other_genes/source_metadata_screen.csv.gz"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(df, name):
    kwargs = dict(index=False, lineterminator="\n", float_format="%.12g")
    if name.endswith(".gz"):
        kwargs["compression"] = dict(method="gzip", mtime=0)
    df.to_csv(HERE / name, **kwargs)
    assert (HERE / name).stat().st_size < 1_200_000


def compact(s):
    return re.sub(r"\s+", "", s).replace("\\", "")


def cells(s):
    return [x.strip() for x in s.strip().strip("|").split("|")]


def source_path(meta, pmid):
    return Path(meta.db_path).parent / "pmc_fulltext" / f"{pmid}_CLEANED.md"


def main():
    baseline = pd.read_csv(FLAGS).query('gene == "LDLR"').copy()
    obs = pd.read_csv(OBS)
    meta = pd.read_csv(META).query('gene == "LDLR"').copy()
    meta["pmid"] = pd.to_numeric(meta.pmid)
    prior_checks = json.loads(
        (
            REPO
            / "docs/evidence/brca2_distance_prior_audit_20260913/source/cross_gene_checks.json"
        ).read_text()
    )
    assert sha(OBS) == prior_checks["original_observation_hashes"]["LDLR"]
    old_counts = obs.groupby("key")[["affected", "unaffected"]].sum()
    np.testing.assert_allclose(
        baseline.affected_literature, baseline.key.map(old_counts.affected).fillna(0)
    )
    np.testing.assert_allclose(
        baseline.unaffected_literature,
        baseline.key.map(old_counts.unaffected).fillna(0),
    )
    assert (
        len(obs) == 1554 and obs.affected.sum() == 4595 and obs.unaffected.sum() == 216
    )
    hashes, sources = {}, {}
    for pmid in [29802317, 25463123, 31491741]:
        path = source_path(meta.loc[meta.pmid.eq(pmid)].iloc[0], pmid)
        sources[pmid] = path.read_text().splitlines()
        hashes[str(path.relative_to(REPO))] = sha(path)
    mi_text = "\n".join(sources[29802317])
    assert "9,956 cases and 8,373 controls" in mi_text
    assert "All MI cases including discovery and replication stages" in mi_text
    assert "we did not check and exclude subjects with FH" in mi_text
    greek_text = "\n".join(sources[25463123])
    assert (
        "subjects participating in the present study were clinically diagnosed with heFH"
        in greek_text
    )
    assert (
        "identified 140 pediatric index cases and 141 relatives with heFH" in greek_text
    )
    greek = {}
    in_table = False
    for lineno, line in enumerate(sources[25463123], 1):
        if line.startswith("Table 2 List of detected LDLR mutations"):
            in_table = True
        elif in_table and line.startswith("|"):
            c = cells(line)
            if len(c) != 8 or not c[1].startswith("c."):
                continue
            index = 0 if c[4].startswith("-") else int(c[4].split()[0])
            relative = int(c[5])
            key = compact(c[1])
            assert key not in greek
            greek[key] = dict(
                source_line=lineno,
                source_row=line,
                source_cdna=key,
                source_preprotein=c[2],
                source_mature_protein=c[3],
                source_index_cases=index,
                source_relatives=relative,
            )
        elif in_table and greek:
            break
    assert len(greek) == 26
    assert sum(x["source_index_cases"] for x in greek.values()) == 140
    assert sum(x["source_relatives"] for x in greek.values()) == 141
    save(pd.DataFrame(greek.values()), "Greek_Table2_exact_rows.csv")
    decisions = []
    selected = obs.loc[obs.pmid.isin([29802317, 25463123])].copy()
    assert len(selected) == 80
    for r in selected.itertuples():
        mr = meta.loc[
            meta.pmid.eq(r.pmid) & meta.key.eq(r.key) & meta.variant_id.eq(r.variant_id)
        ]
        assert len(mr) == 1
        mr = mr.iloc[0]
        provenance = dict(
            source_path=str(source_path(mr, r.pmid).relative_to(REPO)),
            source_location=mr.source_location,
            original_key_quotes=mr.key_quotes,
        )
        if r.pmid == 29802317:
            q = "\n".join(json.loads(mr.key_quotes))
            targets = re.findall(r"Target row: (.*)", q)
            assert len(targets) == 1
            matches = [
                i
                for i, line in enumerate(sources[r.pmid], 1)
                if compact(line) == compact(targets[0])
            ]
            assert len(matches) == 1
            c = cells(targets[0])
            assert c[0] == "19" and c[8] == r.key and float(c[5]) == r.affected
            assert float(c[6]) == r.unaffected == 0
            provenance.update(
                source_line=matches[0],
                source_row=targets[0],
                source_chrom=c[0],
                source_pos=int(c[1]),
                source_ref=c[2],
                source_alt=c[3],
                source_cases=int(c[5]),
                source_controls=int(c[6]),
            )
            a, u = None, None
            action = "quarantine_MI_endpoint_FH_unknown"
            reason = "Exact LDLR chr19 allele row reports myocardial infarction cases/controls. Study did not check/exclude FH; neither FH affected nor FH unaffected counts are established by MI status or an FH database annotation."
            endpoint = "myocardial_infarction"
        else:
            row = greek[compact(r.cdna)]
            # Exact cDNA + original preprotein position/reference avoids mature-number offsets.
            expected = re.search(r"p\.\(([A-Za-z]{3})(\d+)", row["source_preprotein"])
            assert expected and int(expected.group(2)) == r.aa_pos
            aa3 = {
                "Cys": "C",
                "Asp": "D",
                "Glu": "E",
                "Gly": "G",
                "Lys": "K",
                "Gln": "Q",
                "Arg": "R",
                "Ser": "S",
                "Val": "V",
                "Tyr": "Y",
            }
            assert aa3[expected.group(1)] == r.aa_ref
            provenance.update(row)
            a, u = row["source_index_cases"] + row["source_relatives"], 0
            assert a == r.affected + r.unaffected
            action = (
                "retain_correct_FH_total"
                if (a == r.affected and u == r.unaffected)
                else "replace_index_relative_partition_with_FH_carrier_total"
            )
            reason = "Table2 detected-variant row gives mutation-positive index cases and relatives, all drawn from the clinically heFH cohort. Full table sums 140 index/141 relatives, matching molecularly diagnosed totals. Relatives are affected carriers, not unaffected; do not add them twice when frozen A already equals total."
            endpoint = "clinically_diagnosed_heterozygous_FH"
        decisions.append(
            dict(
                gene="LDLR",
                key=r.key,
                pmid=int(r.pmid),
                variant_id=int(r.variant_id),
                vclass=r.vclass,
                original_affected=r.affected,
                original_unaffected=r.unaffected,
                A_to_remove=r.affected,
                U_to_remove=r.unaffected,
                A_to_restore=0 if a is None else a,
                U_to_restore=0 if u is None else u,
                adjudicated_affected=a,
                adjudicated_unaffected=u,
                action=action,
                reason=reason,
                source_endpoint=endpoint,
                primary_url=f"https://pubmed.ncbi.nlm.nih.gov/{r.pmid}/",
                unaffected_semantics="unknown_not_zero"
                if u is None
                else "all enumerated variant carriers in this cohort clinically diagnosed heFH; no unaffected carriers enumerated",
                **provenance,
            )
        )
    decisions = pd.DataFrame(decisions)
    assert not decisions.duplicated(["key", "pmid", "variant_id"]).any()
    save(decisions, "observation_decisions.csv.gz")
    selected_ids = set(zip(decisions.key, decisions.pmid, decisions.variant_id))
    retained = obs.loc[
        [i not in selected_ids for i in zip(obs.key, obs.pmid, obs.variant_id)]
    ].copy()
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
    for direction, label in [("affected", "A"), ("unaffected", "U")]:
        np.testing.assert_allclose(
            updated[f"{direction}_literature"],
            updated[f"refresh_baseline_{direction}_literature"]
            - updated[f"refresh_{label}_to_remove"]
            + updated[f"refresh_{label}_to_restore"],
        )
        assert updated[f"{direction}_literature"].ge(0).all()
    assert (
        updated.affected_literature.sum() == 4535
        and updated.unaffected_literature.sum() == 193
    )
    save(updated, "clinical_input.csv.gz")
    save(updated.loc[updated.key.isin(delta.index)], "clinical_key_deltas.csv")
    save(complete, "corrected_observations.csv.gz")
    # Coverage uses frozen canonical missense keys, before any count-driven unit regrouping.
    files = sorted(
        (
            REPO
            / "docs/evidence/class_matched_penetrance_20260912/analysis/empirical_posteriors"
        ).glob("LDLR.part*.csv.gz")
    )
    frozen = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
    keys = set(
        frozen.loc[frozen.variant_type.eq("missense"), "literature_key"].dropna()
    )
    eligible = obs.loc[obs.key.isin(keys)]
    top = (
        eligible.groupby("pmid")
        .agg(
            observation_rows=("key", "size"),
            affected=("affected", "sum"),
            unaffected=("unaffected", "sum"),
        )
        .sort_values("affected", ascending=False)
        .head(3)
        .reset_index()
    )
    assert set(top.pmid) == {29802317, 25463123, 31491741}
    assert (
        top.observation_rows.sum() == 71
        and top.affected.sum() == 197
        and top.unaffected.sum() == 19
    )
    top["review_status"] = top.pmid.map(
        {
            29802317: "quarantine_MI_endpoint_FH_unknown",
            25463123: "correct_index_relative_partition",
            31491741: "retained_source_consistent_clinical_FH_counts",
        }
    )
    save(top, "canonical_missense_top3_coverage.csv")
    save(
        meta.loc[meta.key.isin(keys) & meta.pmid.isin(top.pmid)],
        "canonical_missense_top3_original_observations.csv.gz",
    )
    corrected_eligible = complete.loc[complete.key.isin(keys)]
    assert (
        corrected_eligible.affected.sum() == 1346
        and corrected_eligible.unaffected.sum() == 95
    )
    retained_reviews = []
    third_text = "\n".join(sources[31491741])
    assert (
        "801 clinically diagnosed HeFH patients, 650 of whom were unrelated"
        in third_text
    )
    for key, token, count in [
        ("C338S", "p.(Cys338Ser) (n = 23", 23),
        ("D433H", "p.(Asp433His) (n = 20", 20),
        ("L568V", "p.(Leu568Val) (n = 19", 19),
    ]:
        r = eligible.loc[eligible.pmid.eq(31491741) & eligible.key.eq(key)]
        assert len(r) == 1 and r.affected.iloc[0] == count and r.unaffected.iloc[0] == 0
        matched = [
            i
            for i, line in enumerate(sources[31491741], 1)
            if compact(token) in compact(line)
        ]
        assert matched
        retained_reviews.append(
            dict(
                gene="LDLR",
                key=key,
                pmid=31491741,
                variant_id=int(r.variant_id.iloc[0]),
                original_affected=count,
                original_unaffected=0,
                decision="retain_source_consistent_FH_count",
                source_lines=";".join(map(str, matched)),
                reason="Explicit per-variant count among unrelated clinically diagnosed HeFH patients; no symptomatic-CAD partition inferred.",
                primary_url="https://pubmed.ncbi.nlm.nih.gov/31491741/",
            )
        )
    save(pd.DataFrame(retained_reviews), "retained_source_review.csv")
    checks = dict(
        passed=True,
        original_observations=len(obs),
        original_A=float(obs.affected.sum()),
        original_U=float(obs.unaffected.sum()),
        reviewed_decision_rows=len(decisions),
        quarantine_MI_rows=64,
        quarantine_MI_A=83,
        greek_rows=16,
        greek_original_A=67,
        greek_original_U=23,
        greek_new_A=90,
        greek_new_U=0,
        greek_already_correct_rows=int(
            decisions.action.eq("retain_correct_FH_total").sum()
        ),
        primary_A=float(updated.affected_literature.sum()),
        primary_U=float(updated.unaffected_literature.sum()),
        clinical_keys=len(updated),
        retained_clinical_keys=int(updated.n_literature.gt(0).sum()),
        canonical_missense_original_rows=len(eligible),
        canonical_missense_original_A=float(eligible.affected.sum()),
        canonical_missense_original_U=float(eligible.unaffected.sum()),
        canonical_missense_new_A=float(corrected_eligible.affected.sum()),
        canonical_missense_new_U=float(corrected_eligible.unaffected.sum()),
        baseline_hashes={
            str(FLAGS.relative_to(REPO)): sha(FLAGS),
            str(OBS): sha(OBS),
            str(META.relative_to(REPO)): sha(META),
        },
        source_hashes=hashes,
        output_ledger_sha256=sha(HERE / "clinical_input.csv.gz"),
        scope="Three highest affected-count sources among frozen canonical missense keys:71 rows/197A/19U; exact two-paper proofs also applied to selected noncanonical and other-type rows with same table semantics. No exhaustive validation of remaining sources.",
        population_rule="Drop zero-clinical-evidence keys before rejoining complete observed gnomAD population alleles. Quarantine means unknown clinical endpoint, never inferred unaffected.",
    )
    (HERE / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in checks.items()
                if k not in ["baseline_hashes", "source_hashes"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
