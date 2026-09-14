"""Bounded audit of KCNQ1's three largest frozen affected-count papers."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
IDENTITY = (
    ROOT / "docs/evidence/population_inclusive_penetrance_20260912/audit/"
    "literature_identity_flags.csv.gz"
)
OBS = (
    ROOT.parent / "BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results/"
    "KCNQ1_protocol/observations.csv"
)
DB = ROOT / "validation_runs/canonical_baseline/KCNQ1.db"
DB_SHA = "9d974e08fd0510192abbe027bb5dee6c298ee28e45d53323c4a5830985ae0398"
XLSX = ROOT / "corpus/KCNQ1/32893267/32893267_supplements/mmc2.xlsx"
AA = dict(
    zip(
        "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter".split(),
        "A R N D C Q E G H I L K M F P S T W Y V X".split(),
    )
)
PMIDS = ["32893267", "23856471", "21244686"]
PMC = {
    "21244686": "PMC3032654",
    "23856471": "PMC3864834",
    "32893267": "PMC7790744",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame: pd.DataFrame, name: str) -> None:
    payload = frame.to_csv(index=False, lineterminator="\n").encode()
    if name.endswith(".gz"):
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, mode="wb", filename="", mtime=0) as f:
            f.write(payload)
        payload = buffer.getvalue()
    assert len(payload) < 1_200_000, name
    (OUT / name).write_bytes(payload)


def protein_key(value: str) -> str:
    value = str(value).replace("p.", "").replace("(", "").replace(")", "")
    value = value.replace("*", "X")
    for three, one in AA.items():
        value = value.replace(three, one)
    return value


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    assert (ROOT / "corpus").is_symlink() and (ROOT / "corpus").is_dir()
    assert sha(DB) == DB_SHA
    observations = pd.read_csv(OBS, dtype={"pmid": str})
    observations["frozen_csv_line"] = observations.index + 2
    literature = pd.read_csv(IDENTITY)
    literature = literature[literature.gene.eq("KCNQ1")].copy()
    assert not literature.key.duplicated().any()
    counts = observations.groupby("key")[["affected", "unaffected"]].sum()
    for row in literature.itertuples():
        assert counts.loc[row.key, "affected"] == row.affected_literature
        assert counts.loc[row.key, "unaffected"] == row.unaffected_literature
    ranking = (
        observations.groupby("pmid")
        .agg(
            rows=("key", "size"),
            affected=("affected", "sum"),
            unaffected=("unaffected", "sum"),
        )
        .sort_values("affected", ascending=False)
        .reset_index()
    )
    assert ranking.head(3).pmid.tolist() == PMIDS
    save(ranking, "frozen_paper_count_ranking.csv")

    with sqlite3.connect(f"file:{DB}?mode=ro", uri=True) as connection:
        for pmid in PMIDS:
            recorded = connection.execute(
                "SELECT pmc_id FROM papers WHERE pmid=?", (pmid,)
            ).fetchone()[0]
            assert recorded is None or recorded == PMC[pmid]
        db_rows = pd.read_sql_query(
            "SELECT p.*, v.cdna_notation, v.protein_notation, "
            "vp.source_location, vp.additional_notes, vp.key_quotes, "
            "vp.count_provenance, vp.source_layer FROM penetrance_data p "
            "JOIN variants v USING (variant_id) "
            "LEFT JOIN variant_papers vp USING (variant_id, pmid) "
            "WHERE p.pmid IN ('32893267','23856471','21244686') "
            "ORDER BY p.pmid, p.variant_id, p.penetrance_id",
            connection,
        )
    save(db_rows, "selected_source_db_rows.csv.gz")

    table = pd.read_excel(XLSX, sheet_name="S4", header=4)
    table["xlsx_row"] = table.index + 6
    table = table[table.Gene.eq("KCNQ1")].copy()
    table["protein_key"] = table.Protein.map(protein_key)
    table["source_case_count"] = table["Cases-Europe"] + table["Cases-Japan"]
    checks = []
    corrections = []
    for row in observations[observations.pmid.eq("32893267")].itertuples():
        match = (
            table[table.CDS.eq(row.cdna)]
            if pd.notna(row.cdna)
            else table[table.protein_key.eq(row.key)]
        )
        assert len(match) == 1, (row.key, row.cdna)
        source = match.iloc[0]
        exact = int(source.source_case_count) == row.affected
        checks.append(
            {
                "key": row.key,
                "pmid": row.pmid,
                "frozen_csv_line": row.frozen_csv_line,
                "frozen_A": row.affected,
                "frozen_U": row.unaffected,
                "source_CDS": source.CDS,
                "source_protein": source.Protein,
                "source_variant_type": source.Vartype,
                "source_cases_Europe": source["Cases-Europe"],
                "source_cases_Japan": source["Cases-Japan"],
                "source_case_count": source.source_case_count,
                "source_sheet": "S4",
                "source_xlsx_row": source.xlsx_row,
                "exact_frozen_count_match": exact,
                "endpoint": "clinically_diagnosed_ECG_assessed_LQTS",
                "decision": "retain_explicit_case_count"
                if exact
                else "remove_duplicate_A1",
            }
        )
        if not exact:
            assert row.key in ["A344=", "c.477+5G>A"]
            assert row.affected == 1 and row.unaffected == 0
            retained = observations[
                observations.key.eq(row.key)
                & observations.pmid.eq(row.pmid)
                & observations.affected.eq(source.source_case_count)
            ]
            assert len(retained) == 1
            corrections.append(
                decision(
                    row,
                    "remove_duplicate_A1",
                    f"Supplementary Table S4, XLSX row {int(source.xlsx_row)}",
                    "An additional A1 repeats the same variant/paper and the larger explicit S4 case count. Retain the verified Europe+Japan count once.",
                    True,
                )
            )
    assert len(checks) == 270
    assert sum(not x["exact_frozen_count_match"] for x in checks) == 2
    save(pd.DataFrame(checks), "32893267_case_count_checks.csv")

    for row in observations[observations.pmid.eq("21244686")].itertuples():
        assert row.key in ["G589D", "Y171X"]
        corrections.append(
            decision(
                row,
                "quarantine_unpartitioned_genotype_carrier_total",
                "Methods, Patients; Table 1/Table 2; archived additional_notes",
                (
                    "492 is all available G589D mutation carriers, not a count with manifest LQTS under an explicit ECG/clinical threshold. QTc is continuous, while syncope history exists for only 488. No broad LQTS A/U partition can be reconstructed."
                    if row.key == "G589D"
                    else "The single Y171X carrier is described only as also carrying G589D. No variant-specific clinical/ECG phenotype partition is provided; compound genotype retained as provenance."
                ),
                True,
            )
        )
    for row in observations[observations.pmid.eq("23856471")].itertuples():
        if row.key == "c.386+18089C>T":
            corrections.append(
                decision(
                    row,
                    "quarantine_whole_study_size_assigned_to_modifier_allele",
                    "Abstract; Methods discovery/replication cohorts; archived additional_notes",
                    "560 equals 224 discovery participants +152 South African +184 Finnish participants. This is not the number of rs2074238 T-allele carriers. The source quotes do not support A560/U0 for this allele.",
                    True,
                )
            )
        else:
            assert row.key in ["A341V", "G589D"]
            corrections.append(
                decision(
                    row,
                    "separate_cardiac_events_endpoint",
                    "Methods, Study population/inclusion criteria and Replication populations",
                    "The A/U split is cardiac events before age35 versus event-free older untreated carriers. Event-free LQTS patients are not established LQTS-unaffected. Preserve the valid event split only in the explicitly mixed endpoint sensitivity.",
                    False,
                )
            )
    delta = pd.DataFrame(corrections)
    assert len(delta) == 7
    save(delta, "observation_corrections.csv")
    save(
        delta[delta.decision.eq("separate_cardiac_events_endpoint")],
        "preserved_cardiac_event_observations.csv",
    )

    ledger = literature.copy()
    ledger["old_affected_literature"] = ledger.affected_literature
    ledger["old_unaffected_literature"] = ledger.unaffected_literature
    for sensitivity, filename in [
        (False, "clinical_input.csv.gz"),
        (True, "clinical_input_cardiac_events.csv.gz"),
    ]:
        arm = ledger.copy()
        selected = delta[delta.remove_from_event_sensitivity] if sensitivity else delta
        summed = selected.groupby("key")[
            ["affected_removed", "unaffected_removed"]
        ].sum()
        for name in ["affected_removed", "unaffected_removed"]:
            arm[name] = arm.key.map(summed[name]).fillna(0)
        arm["affected_restored"] = 0
        arm["unaffected_restored"] = 0
        arm["affected_literature"] -= arm.affected_removed
        arm["unaffected_literature"] -= arm.unaffected_removed
        arm["new_affected_literature"] = arm.affected_literature
        arm["new_unaffected_literature"] = arm.unaffected_literature
        arm["n_literature"] = arm.affected_literature + arm.unaffected_literature
        arm["zero_evidence_after_removal"] = arm.n_literature.eq(0)
        arm["review_status"] = "retained_not_fully_source_adjudicated"
        changed = (arm.affected_removed + arm.unaffected_removed).gt(0)
        arm.loc[changed, "review_status"] = "reviewed_source_counts_corrected"
        arm["cardiac_event_split_retained"] = sensitivity & arm.key.isin(
            ["A341V", "G589D"]
        )
        arm["endpoint_scope"] = (
            "mixed_clinical_LQTS_plus_explicit_cardiac_events_sensitivity"
            if sensitivity
            else "clinical_ECG_LQTS_bounded_source_correction"
        )
        assert arm[["affected_literature", "unaffected_literature"]].ge(0).all().all()
        save(arm, filename)

    sources = [IDENTITY, OBS, DB, XLSX, Path(__file__)]
    sources += [ROOT / f"corpus/KCNQ1/{p}/{p}_CLEANED.md" for p in PMIDS]
    primary = pd.read_csv(OUT / "clinical_input.csv.gz")
    sensitivity = pd.read_csv(OUT / "clinical_input_cardiac_events.csv.gz")
    eligible = primary.vclass.eq("missense") & primary.canonical_wt_status.eq("match")
    result = {
        "scope": "largest_three_frozen_affected_count_papers_only",
        "paper_ranking": PMIDS,
        "clinical_keys": len(primary),
        "adjudicated_correction_rows": len(delta),
        "primary_affected_removed_all_types": float(delta.affected_removed.sum()),
        "primary_unaffected_removed_all_types": float(delta.unaffected_removed.sum()),
        "primary_missense_affected": float(
            primary.loc[eligible, "affected_literature"].sum()
        ),
        "primary_missense_unaffected_literature": float(
            primary.loc[eligible, "unaffected_literature"].sum()
        ),
        "sensitivity_missense_affected": float(
            sensitivity.loc[eligible, "affected_literature"].sum()
        ),
        "sensitivity_missense_unaffected_literature": float(
            sensitivity.loc[eligible, "unaffected_literature"].sum()
        ),
        "gnomad_counts_edited": False,
        "unknowns_converted_to_unaffected": 0,
        "clinical_ECG_LQTS_cases_removed_for_no_symptoms": 0,
        "full_source_or_person_ownership_adjudication": False,
        "source_databases_edited": False,
        "source_sha256": {str(p): sha(p) for p in sources},
    }
    (OUT / "source_checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in result.items() if k != "source_sha256"}, indent=2)
    )


def decision(
    row, action: str, location: str, rationale: str, remove_event: bool
) -> dict:
    record = row._asdict()
    record.pop("Index", None)
    source = ROOT / f"corpus/KCNQ1/{row.pmid}/{row.pmid}_CLEANED.md"
    record.update(
        {
            "gene": "KCNQ1",
            "decision": action,
            "affected_removed": row.affected,
            "unaffected_removed": row.unaffected,
            "affected_restored": 0,
            "unaffected_restored": 0,
            "remove_from_event_sensitivity": remove_event,
            "source_location": location,
            "source_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{PMC[row.pmid]}/",
            "source_path": str(source.relative_to(ROOT)),
            "source_sha256": sha(source),
            "rationale": rationale,
            "unknown_converted_to_unaffected": False,
        }
    )
    return record


if __name__ == "__main__":
    main()
