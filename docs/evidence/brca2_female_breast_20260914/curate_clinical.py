"""Apply source repairs and build a strict, source-verified female BC subset.

Archived databases stay immutable. Every retained or excluded observation has a
ledger entry; an unknown clinical sex/endpoint never becomes an unaffected count.
"""

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OLD = HERE.parent / "brca2_distance_prior_audit_20260913/source"
BASE = HERE.parent / "residue_density_refresh_20260914"
HASHES = {}


def pin(path):
    HASHES[str(path.relative_to(REPO))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read(path):
    return pd.read_csv(pin(path), low_memory=False)


def save(frame, name):
    blob = frame.to_csv(index=False, lineterminator="\n").encode()
    (HERE / name).write_bytes(
        gzip.compress(blob, mtime=0) if name.endswith(".gz") else blob
    )


def source(pmid):
    paths = list(
        (REPO / "results/grant_e2e_20260909/shards").glob(
            f"BRCA2_*/BRCA2/*/pmc_fulltext/{pmid}_CLEANED.md"
        )
    )
    assert len(paths) == 1, (pmid, paths)
    return pin(paths[0])


def main():
    clinical = read(BASE / "source/BRCA2/BRCA2_literature_identity_corrected.csv.gz")
    obs = read(OLD / "original_observations.csv.gz")
    obs["source_row_id"] = [f"original:{i}" for i in range(len(obs))]
    for path in [
        OLD / "BRCA2_catalog_quarantine_observations.csv",
        BASE / "source/BRCA2/BRCA2_observation_decisions.csv.gz",
    ]:
        for row in read(path).itertuples():
            mask = (
                obs.key.eq(row.original_key)
                & obs.pmid.eq(row.pmid)
                & obs.variant_id.eq(row.original_observation_variant_id)
            )
            assert mask.sum() == 1
            obs.loc[mask, "affected"] -= row.A_to_remove
            obs.loc[mask, "unaffected"] -= row.U_to_remove
    restored = read(OLD / "BRCA2_clinical_restoration_observations.csv")
    extra = pd.DataFrame(
        dict(
            key=restored.protein_key,
            pmid=restored.pmid,
            affected=restored.A_to_restore,
            unaffected=restored.U_to_restore,
            variant_id=-np.arange(1, len(restored) + 1),
            vclass="missense",
            source_row_id=[f"restored:{i}" for i in range(len(restored))],
        )
    )
    obs = pd.concat([obs, extra], ignore_index=True)
    assert obs[["affected", "unaffected"]].ge(0).all().all()
    obs = obs.loc[(obs.affected + obs.unaffected).gt(0)].copy()
    totals = (
        obs.groupby("key")[["affected", "unaffected"]]
        .sum()
        .reindex(clinical.key)
        .fillna(0)
    )
    np.testing.assert_allclose(
        totals.to_numpy(),
        clinical[["affected_literature", "unaffected_literature"]].to_numpy(),
    )

    # Reviewable repairs independent of the new sex restriction.
    patches = [
        (
            "E1581D",
            17100994,
            40000141,
            48,
            0,
            0,
            "BRCA1 table, not BRCA2",
            "65-107; row 85",
        ),
        (
            "I3412V",
            17100994,
            40000163,
            0,
            0,
            5,
            "Recover explicit normal-control column; BIC entries are not people",
            "47; 112-145",
        ),
        (
            "V109G",
            16672066,
            50017508,
            142,
            137,
            0,
            "P27 genotype block, not BRCA2; GG+GT = 17+125 cases and 15+122 controls",
            "152; 272-274",
        ),
    ]
    repair_rows = []
    for key, pmid, vid, remove_a, remove_u, add_u, reason, lines in patches:
        text = source(pmid).read_text()
        assert key == "I3412V" or (
            "E1581D" in text if key == "E1581D" else "| P27 | TT |" in text
        )
        mask = obs.key.eq(key) & obs.pmid.eq(pmid) & obs.variant_id.eq(vid)
        assert mask.sum() == 1
        repair_rows.append(
            dict(
                key=key,
                pmid=pmid,
                variant_id=vid,
                A_removed=remove_a,
                U_removed=remove_u,
                U_restored=add_u,
                reason=reason,
                source_lines=lines,
                source_path=str(source(pmid).relative_to(REPO)),
            )
        )
        obs.loc[mask, "affected"] -= remove_a
        obs.loc[mask, "unaffected"] += add_u - remove_u
    assert obs[["affected", "unaffected"]].ge(0).all().all()
    save(pd.DataFrame(repair_rows), "source_repairs.csv")

    def aggregate(frame, name):
        # Carry identity evidence, not stale counters from earlier corrections.
        result = clinical[
            [
                c
                for c in clinical
                if not c.startswith(
                    ("original_", "catalog_", "refresh_", "clinical_A_", "clinical_U_")
                )
            ]
        ].copy()
        t = frame.groupby("key")[["affected", "unaffected"]].sum()
        result["affected_literature"] = result.key.map(t.affected).fillna(0)
        result["unaffected_literature"] = result.key.map(t.unaffected).fillna(0)
        result["n_literature"] = (
            result.affected_literature + result.unaffected_literature
        )
        result["pmids"] = result.key.map(
            frame.groupby("key").pmid.agg(lambda s: ";".join(map(str, sorted(set(s)))))
        ).fillna("")
        result["source_key_status"] = name
        result["drop_clinical_key_zero_clinical_evidence"] = result.n_literature.eq(0)
        result = result.loc[result.n_literature.gt(0)].copy()
        save(result, name + ".csv.gz")
        return result

    aggregate(obs, "clinical_source_corrected_all_sex")
    obs["female_affected"] = 0.0
    obs["female_unaffected"] = 0.0
    obs["eligibility_reason"] = "unresolved_source_sex_germline_or_breast_endpoint"
    # Cohort-level sex evidence is valid for its rows; individual-level tags are
    # not required when the primary cohort is explicitly all women.
    cohorts = {
        28283652: (
            "31; 7; 19",
            "Analyses were restricted to either Caucasian or Asian women",
            "female_BCAC_germline_case_control",
        ),
        28222693: (
            "21-35; 45-47; 78",
            "1,530 women with no personal history of breast cancer",
            "female_MyBrCa_blood_case_control",
        ),
        29755871: (
            "31; 37-41; 71",
            "156 Kazakhstan women",
            "female_Kazakhstan_blood_case_control",
        ),
        38863777: (
            "7; 25-29; 70-74; 117-131",
            "Seventy women",
            "female_Kurdish_blood_breast_cases",
        ),
        33278427: (
            "25-31; 152-156",
            "103 unrelated Egyptian female breast cancer patients",
            "female_Egypt_blood_breast_cases",
        ),
        37719058: (
            "7; 17-25; 111",
            "selected women with TNBC",
            "female_Kenya_saliva_breast_cases",
        ),
        37816281: (
            "17; 23",
            "54-year-old white female",
            "female_breast_patient_PBMC_not_cell_count",
        ),
    }
    cohort_rows = []
    for pmid, (lines, anchor, reason) in cohorts.items():
        path = source(pmid)
        assert anchor in path.read_text(), (pmid, anchor)
        mask = obs.pmid.eq(pmid)
        obs.loc[mask, ["female_affected", "female_unaffected"]] = obs.loc[
            mask, ["affected", "unaffected"]
        ].to_numpy()
        obs.loc[mask, "eligibility_reason"] = reason
        cohort_rows.append(
            dict(
                pmid=pmid,
                source_path=str(path.relative_to(REPO)),
                source_lines=lines,
                decision=reason,
            )
        )
    # Keep the already verified normal women in Han et al.; sex of its 793
    # patient cohort is not explicitly enumerated in the cached Methods.
    mask = obs.key.eq("I3412V") & obs.pmid.eq(17100994)
    obs.loc[mask, "female_unaffected"] = 5
    obs.loc[mask, "eligibility_reason"] = (
        "five_confirmed_female_controls_case_sex_unresolved"
    )
    cohort_rows.append(
        dict(
            pmid=17100994,
            source_path=str(source(17100994).relative_to(REPO)),
            source_lines="47; 145",
            decision="retain_5_female_controls_only",
        )
    )
    # Do not sum a shared allele across a regional report and BCAC while
    # participant overlap is unresolved. This is not proof of overlap.
    overlap = obs.pmid.eq(28222693) & obs.key.eq("K2729N")
    obs.loc[overlap, ["female_affected", "female_unaffected"]] = 0
    obs.loc[overlap, "eligibility_reason"] = (
        "possible_BCAC_cohort_overlap_shared_allele_not_added"
    )

    # The largest original study has separate female and male supplements.
    # Rebuild from female SD1 instead of guessing sex from duplicate DB keys.
    pmid = 30287823
    text_path = source(pmid)
    assert "We analyzed women and men separately" in text_path.read_text()
    workbook = next(text_path.parent.glob("30287823_supplements/**/*MOESM4*.xlsx"))
    pin(workbook)
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    sheet = book["SD1"]
    rows = list(sheet.values)
    book.close()
    assert "in women" in rows[0][0]
    assert "7,051 cases" in rows[1][11] and "11,241 controls" in rows[1][12]
    cdna_keys = (
        clinical.dropna(subset=["cdna"]).groupby("cdna").key.agg(lambda s: set(s))
    )
    recovered, unmapped = [], []
    for excel_row, r in enumerate(rows[2:], 3):
        if r[5] != "BRCA2" or not any(
            t in str(r[6]).split("&") for t in ["missense_variant", "stop_gained"]
        ):
            continue
        if r[7] not in cdna_keys or len(cdna_keys[r[7]]) != 1:
            unmapped.append(
                dict(
                    excel_row=excel_row,
                    cdna=r[7],
                    protein=r[9],
                    reason="no_unique_existing_clinical_identity",
                )
            )
            continue
        key = next(iter(cdna_keys[r[7]]))
        counts = []
        for freq, total in [(r[11], 7051), (r[12], 11241)]:
            # >=98% variant call rate does not give exact per-allele denominators.
            # Allow 620 missing calls in either stratum, conservatively above
            # 2% of all 30,926 initially described samples. Keep only an integer
            # identifiable throughout that range and five-decimal rounding.
            lo = max(0, float(freq) - 0.000005) * (total - 620)
            hi = min(1, float(freq) + 0.000005) * total
            if np.ceil(lo - 1e-9) != np.floor(hi + 1e-9):
                counts = []
                break
            counts.append(int(np.ceil(lo - 1e-9)))
        if len(counts) != 2:
            unmapped.append(
                dict(
                    excel_row=excel_row,
                    cdna=r[7],
                    protein=r[9],
                    reason="carrier_count_not_identified_with_missing_call_and_rounding_bounds",
                )
            )
            continue
        if not sum(counts):
            continue
        recovered.append(
            dict(
                key=key,
                pmid=pmid,
                variant_id=-100000 - excel_row,
                affected=counts[0],
                unaffected=counts[1],
                female_affected=counts[0],
                female_unaffected=counts[1],
                source_row_id=f"female_SD1:{excel_row}",
                eligibility_reason="female_SD1_unique_integer_from_carrier_frequency",
                source_cdna=r[7],
                source_protein=r[9],
                source_excel_row=excel_row,
            )
        )
    obs.loc[obs.pmid.eq(pmid), "eligibility_reason"] = (
        "superseded_by_female_SD1_source_rows"
    )
    obs = pd.concat([obs, pd.DataFrame(recovered)], ignore_index=True)
    obs["status"] = np.where(
        (obs.female_affected + obs.female_unaffected).gt(0),
        "retained_or_partitioned",
        "excluded",
    )
    save(
        obs[
            [
                "source_row_id",
                "key",
                "pmid",
                "variant_id",
                "affected",
                "unaffected",
                "female_affected",
                "female_unaffected",
                "status",
                "eligibility_reason",
                "source_cdna",
                "source_protein",
                "source_excel_row",
            ]
        ],
        "observation_eligibility.csv.gz",
    )
    save(pd.DataFrame(unmapped), "female_SD1_identity_queue.csv")
    cohort_rows.append(
        dict(
            pmid=pmid,
            source_path=str(workbook.relative_to(REPO)),
            source_lines="SD1; Methods lines 120,128",
            decision="female_source_table_rebuilt_unique_integer_counts",
        )
    )
    save(pd.DataFrame(cohort_rows), "cohort_decisions.csv")
    female = obs.loc[obs.status.eq("retained_or_partitioned")].copy()
    female[["affected", "unaffected"]] = female[
        ["female_affected", "female_unaffected"]
    ].to_numpy()
    aggregate(female, "clinical_female_breast")
    attrition = (
        obs.groupby(["pmid", "eligibility_reason"])[
            ["affected", "unaffected", "female_affected", "female_unaffected"]
        ]
        .sum()
        .reset_index()
    )
    save(attrition, "clinical_attrition.csv")
    assert not female.key.eq("E1581D").any()
    assert not (female.key.eq("D2723H") & female.pmid.eq(20927582)).any()
    assert not (female.key.eq("V109G") & female.pmid.eq(16672066)).any()
    pin(Path(__file__).resolve())
    (HERE / "clinical_input_hashes.json").write_text(
        json.dumps(HASHES, indent=2) + "\n"
    )
    print(
        f"Female source rows: {len(female)}; studies: {female.pmid.nunique()}; A={female.affected.sum():g}; U={female.unaffected.sum():g}; SD1 identity queue={len(unmapped)}"
    )


if __name__ == "__main__":
    main()
