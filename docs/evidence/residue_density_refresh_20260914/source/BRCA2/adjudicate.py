#!/usr/bin/env python3
"""Adjudicate two frozen BRCA2 source queues without changing original evidence.

The original affected fact is traced to its exact public supplement row. A
tumor observation without established germline status, or an unpartitioned
pooled carrier total, cannot contribute an established hereditary-disease A/U.
Unknown is quarantined; it is never converted to an unaffected observation.
"""

from __future__ import annotations

import hashlib
import json
import re
import sqlite3
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
OLD = REPO / "docs/evidence/brca2_distance_prior_audit_20260913/source"
RUN = REPO / "results/grant_e2e_20260909/shards/BRCA2_3/BRCA2/20260909_115046"
BASELINE_HASH = "8202b4021bb0466f00c31fbcdc1a50bbd6474c1831eb1f06d4ee80fb13722cf6"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame, name):
    kwargs = dict(index=False, lineterminator="\n", float_format="%.12g")
    if name.endswith(".gz"):
        kwargs["compression"] = dict(method="gzip", mtime=0)
    frame.to_csv(HERE / name, **kwargs)
    assert (HERE / name).stat().st_size < 1_200_000


def text(value):
    return str(value or "").strip()


def norm(value):
    return re.sub(r"\s+", " ", text(value))


def workbook(path, sheet):
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        return list(book[sheet].values)
    finally:
        book.close()


def main():
    baseline_path = OLD / "BRCA2_literature_identity_corrected.csv.gz"
    assert sha(baseline_path) == BASELINE_HASH
    baseline = pd.read_csv(baseline_path)
    observations = pd.read_csv(OLD / "original_observations.csv.gz")
    queue = pd.read_csv(OLD / "remaining_BRCA2_source_queue.csv.gz")
    selected = observations.loc[
        observations.pmid.isin([36385461, 33054725])
        & observations.vclass.isin(["missense", "nonsense"])
    ].copy()
    assert len(selected) == 508 and selected.affected.eq(1).all()
    assert selected.unaffected.eq(0).all()
    assert (selected.variant_id // 10_000_000).eq(5).all()

    source_md = {p: RUN / f"pmc_fulltext/{p}_CLEANED.md" for p in [36385461, 33054725]}
    md_lines = {p: path.read_text().splitlines() for p, path in source_md.items()}
    source_line_index = {}
    for pmid, lines in md_lines.items():
        source_line_index[pmid] = {}
        for i, line in enumerate(lines, 1):
            source_line_index[pmid].setdefault(norm(line), []).append(i)
    asian_path = next(
        (RUN / "pmc_fulltext/36385461_supplements").rglob("IJC-152-1159-s011.xlsx")
    )
    tcga_path = next(
        (RUN / "pmc_fulltext/33054725_supplements").rglob(
            "12885_2020_7481_MOESM3_ESM.xlsx"
        )
    )
    sysucc_path = next(
        (RUN / "pmc_fulltext/33054725_supplements").rglob(
            "12885_2020_7481_MOESM4_ESM.xlsx"
        )
    )
    asian = workbook(asian_path, "BRCA2")
    tcga = workbook(tcga_path, "Sheet1")
    assert asian[1][64] == "Carrier number***"
    assert norm(asian[2][65]) == "Cancer" and norm(asian[2][66]) == "Non-cancer"
    assert tcga[1][0] == "ID" and tcga[1][3] == "Muts"
    assert any("data on germline BRCA mutations" in line for line in md_lines[33054725])
    assert any(
        '"1" is given for the variants without carrier number' in line
        for line in md_lines[36385461]
    )

    db_path = RUN / "BRCA2.db"
    frozen_manifest = pd.read_csv(OLD / "cross_gene_db_provenance.csv")
    expected_db_hash = frozen_manifest.loc[
        frozen_manifest.gene.eq("BRCA2") & frozen_manifest.db_index.eq(5), "sha256"
    ].item()
    assert sha(db_path) == expected_db_hash
    decisions, snapshots = [], []
    with sqlite3.connect(f"file:{db_path}?mode=ro", uri=True) as db:
        for row in selected.itertuples():
            local_id = int(row.variant_id) % 10_000_000
            args = (local_id, int(row.pmid))
            vp = db.execute(
                "SELECT source_location,count_provenance,source_layer FROM variant_papers WHERE variant_id=? AND pmid=?",
                args,
            ).fetchone()
            facts = db.execute(
                "SELECT fact_value,source_table,source_row,evidence_quote,provenance_kind FROM fact_provenance "
                "WHERE variant_id=? AND pmid=? AND fact_type='affected_count' AND evidence_quote IS NOT NULL",
                args,
            ).fetchall()
            # The selected A=1 always has an exact evidence row; other facts may
            # report different pooled totals and must not be interpreted as A.
            candidates = [f for f in facts if float(f[0]) == 1]
            unique = {f[3]: f for f in candidates}
            assert len(unique) == 1, (row.key, row.pmid, unique)
            quote, fact = next(iter(unique.items()))
            cells = [v.strip() for v in quote.strip().strip("|").split("|")]
            matches = source_line_index[row.pmid].get(norm(quote), [])
            assert matches, (
                row.key,
                row.pmid,
                "source row absent from full primary source",
            )
            record = dict(
                gene="BRCA2",
                original_key=row.key,
                pmid=int(row.pmid),
                original_observation_variant_id=int(row.variant_id),
                db_variant_id=local_id,
                vclass=row.vclass,
                original_affected=float(row.affected),
                original_unaffected=float(row.unaffected),
                A_to_remove=float(row.affected),
                U_to_remove=float(row.unaffected),
                A_to_restore=0,
                U_to_restore=0,
                source_db_location=vp[0],
                affected_fact_table=fact[1],
                affected_fact_row=fact[2],
                affected_fact_provenance_kind=fact[4],
                source_markdown_lines=";".join(map(str, matches)),
                source_germline_status="unknown",
                source_affected_carriers=None,
                source_unaffected_carriers=None,
                action="quarantine_from_hereditary_endpoint_fit",
                adjudication_status="resolved_for_current_fit_endpoint_unknown",
                primary_url=f"https://pmc.ncbi.nlm.nih.gov/articles/{'PMC10098510' if row.pmid == 36385461 else 'PMC7556962'}/",
            )
            if row.pmid == 36385461:
                assert fact[1] == "Table T24"
                assert cells[9] == "BRCA2" and cells[7].startswith("c.")
                matched = [
                    (i + 1, r)
                    for i, r in enumerate(asian)
                    if len(r) > 14 and norm(r[7]) == cells[7] and norm(r[8]) == cells[8]
                ]
                assert len(matched) == 1, (row.key, cells[:15], len(matched))
                xls_row, r = matched[0]
                assert all(norm(r[i]) == cells[i] for i in range(15)), (
                    row.key,
                    "identity mismatch",
                )
                record.update(
                    source_file=asian_path.name,
                    source_sheet="BRCA2",
                    source_excel_row=xls_row,
                    source_table="Supplementary Table S7B",
                    source_cdna=r[7],
                    source_protein=r[8],
                    source_consequence=r[11],
                    source_carrier_total=r[64],
                    source_tested_cancer=r[65],
                    source_tested_noncancer=r[66],
                    source_tested_total=r[67],
                    source_references=r[69],
                    source_countries=r[68],
                    source_endpoint="Pooled cancer and non-cancer tested cohorts; carrier endpoint unpartitioned",
                    source_classification=r[13],
                    reason="Clinical classification was used as phenotype and A=1 inferred from an inventory row. Pooled carrier totals and tested-cohort denominators do not establish variant-specific hereditary-endpoint A/U; source may substitute carrier total 1 when unavailable; source overlap unadjudicated.",
                )
            else:
                assert fact[1] == "Table T4" and cells[1] == "BRCA2"
                matched = [
                    (i + 1, r)
                    for i, r in enumerate(tcga)
                    if norm(r[0]) == cells[0] and norm(r[3]) == cells[3]
                ]
                assert len(matched) == 1, (row.key, cells, len(matched))
                xls_row, r = matched[0]
                assert all(norm(r[i]) == cells[i] for i in range(6)), (
                    row.key,
                    "identity mismatch",
                )
                record.update(
                    source_file=tcga_path.name,
                    source_sheet="Sheet1",
                    source_excel_row=xls_row,
                    source_table="Supplementary Table 1, TCGA",
                    source_subject=r[0],
                    source_protein=r[3],
                    source_consequence=r[5],
                    source_endpoint=r[2],
                    source_carrier_total=1,
                    source_germline_status="unavailable_for_TCGA_per_primary_paper",
                    reason="Actual tumor sample row; germline mutation data explicitly unavailable for this cohort. The observed tumor cannot be counted as an established hereditary variant carrier. Do not label all rows proven somatic or infer unaffected status.",
                )
            decisions.append(record)
            snapshots.append(
                dict(
                    original_key=row.key,
                    pmid=int(row.pmid),
                    original_observation_variant_id=int(row.variant_id),
                    source_location=vp[0],
                    count_provenance=vp[1],
                    source_layer=vp[2],
                    affected_fact_quote=quote,
                )
            )

    decisions = pd.DataFrame(decisions)
    assert (
        decisions[["original_key", "pmid", "original_observation_variant_id"]]
        .duplicated()
        .sum()
        == 0
    )
    expected_queue_ids = set(zip(queue.key, queue.pmid, queue.variant_id))
    decision_ids = set(
        zip(
            decisions.original_key,
            decisions.pmid,
            decisions.original_observation_variant_id,
        )
    )
    assert expected_queue_ids <= decision_ids and len(expected_queue_ids) == 227
    decisions["in_previous_canonical_missense_queue"] = [
        identity in expected_queue_ids
        for identity in zip(
            decisions.original_key,
            decisions.pmid,
            decisions.original_observation_variant_id,
        )
    ]
    save(decisions, "BRCA2_observation_decisions.csv.gz")
    save(pd.DataFrame(snapshots), "selected_source_records.csv.gz")

    updated = baseline.copy()
    delta = decisions.groupby("original_key")[
        ["A_to_remove", "U_to_remove", "A_to_restore", "U_to_restore"]
    ].sum()
    assert set(delta.index) <= set(updated.key)
    updated["refresh_baseline_affected_literature"] = updated.affected_literature
    updated["refresh_baseline_unaffected_literature"] = updated.unaffected_literature
    for col in delta.columns:
        updated[f"refresh_{col}"] = updated.key.map(delta[col]).fillna(0)
    updated["affected_literature"] += (
        updated.refresh_A_to_restore - updated.refresh_A_to_remove
    )
    updated["unaffected_literature"] += (
        updated.refresh_U_to_restore - updated.refresh_U_to_remove
    )
    assert (
        updated.affected_literature.ge(0).all()
        and updated.unaffected_literature.ge(0).all()
    )
    updated["n_literature"] = (
        updated.affected_literature + updated.unaffected_literature
    )
    updated["drop_clinical_key_zero_clinical_evidence"] = updated.n_literature.eq(0)
    # Keep other papers and all five previous clinical restorations. The exact
    # selected identities, not the PMID globally, are removed from membership.
    excluded = observations.apply(
        lambda r: (r.key, r.pmid, r.variant_id) in decision_ids, axis=1
    )
    old_quarantine = pd.read_csv(OLD / "BRCA2_catalog_quarantine_observations.csv")
    old_ids = set(
        zip(
            old_quarantine.original_key,
            old_quarantine.pmid,
            old_quarantine.original_observation_variant_id,
        )
    )
    excluded |= observations.apply(
        lambda r: (r.key, r.pmid, r.variant_id) in old_ids, axis=1
    )
    retained = observations.loc[~excluded, ["key", "pmid", "affected", "unaffected"]]
    restored = pd.read_csv(OLD / "BRCA2_clinical_restoration_observations.csv")
    restorations = pd.DataFrame(
        dict(
            key=restored.protein_key,
            pmid=40664060,
            affected=restored.A_to_restore,
            unaffected=restored.U_to_restore,
        )
    )
    complete_counts = pd.concat([retained, restorations])
    counts = complete_counts.groupby("key")[["affected", "unaffected"]].sum()
    np.testing.assert_allclose(
        updated.affected_literature, updated.key.map(counts.affected).fillna(0)
    )
    np.testing.assert_allclose(
        updated.unaffected_literature, updated.key.map(counts.unaffected).fillna(0)
    )
    pmids = complete_counts.groupby("key").pmid.agg(
        lambda s: ";".join(map(str, sorted(set(s))))
    )
    updated["pmids"] = updated.key.map(pmids).fillna("")
    for key in restored.protein_key:
        b, u = (
            baseline.loc[baseline.key.eq(key)].iloc[0],
            updated.loc[updated.key.eq(key)].iloc[0],
        )
        assert u.affected_literature == b.affected_literature and "40664060" in u.pmids
    save(updated, "BRCA2_literature_identity_corrected.csv.gz")
    save(updated.loc[updated.key.isin(delta.index)], "BRCA2_clinical_key_deltas.csv.gz")
    save(restored, "preserved_40664060_restorations.csv")

    # The other cohort is inspected for explicit germline/endpoint restoration;
    # it is not assumed equivalent to the selected TCGA observations.
    sysucc = workbook(sysucc_path, "Sheet1")
    sysucc_rows = []
    for i, r in enumerate(sysucc[2:], 3):
        if r[1] == "BRCA2" and norm(r[5]).lower() == "germline":
            sysucc_rows.append(
                dict(
                    source_excel_row=i,
                    sample_id=r[0],
                    gene=r[1],
                    extraneous_repeated_cancer_column=r[2],
                    observed_endpoint=r[3],
                    variant=r[4],
                    germline_status=r[5],
                    consequence=r[7],
                    restoration_A=0,
                    restoration_U=0,
                    decision="No restoration to breast/ovarian endpoint: actual SYSUCC endpoint is not breast/ovarian; not assumed unaffected.",
                )
            )
    assert not any(
        "Breast" in r["observed_endpoint"] or "Ovarian" in r["observed_endpoint"]
        for r in sysucc_rows
    )
    save(pd.DataFrame(sysucc_rows), "SYSUCC_germline_review.csv")
    summary = decisions.groupby(["pmid", "vclass"], as_index=False).agg(
        observations=("original_key", "size"),
        A_removed=("A_to_remove", "sum"),
        U_removed=("U_to_remove", "sum"),
        A_restored=("A_to_restore", "sum"),
        canonical_missense_queue=("in_previous_canonical_missense_queue", "sum"),
    )
    save(summary, "decision_summary.csv")
    receipt = dict(
        passed=True,
        baseline_commit="c8e34c96",
        baseline_clinical_ledger_sha256=BASELINE_HASH,
        frozen_db_hash_matches=True,
        exact_source_row_matches=len(decisions),
        queued_canonical_missense_observations_resolved=227,
        all_typed_observations_adjudicated=len(decisions),
        clinical_keys=len(updated),
        retained_clinical_keys=int(updated.n_literature.gt(0).sum()),
        preserved_previous_real_clinical_restorations=5,
        all_original_affected_delta=float(
            updated.affected_literature.sum() - baseline.affected_literature.sum()
        ),
        all_original_unaffected_delta=float(
            updated.unaffected_literature.sum() - baseline.unaffected_literature.sum()
        ),
        unknown_policy="Quarantine from hereditary endpoint A/U fit; retain source totals and unresolved endpoints; never assign U from missing germline or endpoint data.",
        scope="Two known BRCA2 source queues; missense and nonsense source observations only. Other paper/type observations remain unchanged; this is not full source validation.",
        population_rule="Root must rebuild the union after dropping zero-clinical-evidence keys; all observed gnomAD alleles retained and rejoined by exact identity.",
        hashes={
            str(p.relative_to(REPO)): sha(p)
            for p in [
                baseline_path,
                OLD / "original_observations.csv.gz",
                OLD / "remaining_BRCA2_source_queue.csv.gz",
                db_path,
                asian_path,
                tcga_path,
                sysucc_path,
                *source_md.values(),
            ]
        },
        output_ledger_sha256=sha(HERE / "BRCA2_literature_identity_corrected.csv.gz"),
    )
    (HERE / "checks.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(json.dumps({k: v for k, v in receipt.items() if k != "hashes"}, indent=2))


if __name__ == "__main__":
    main()
