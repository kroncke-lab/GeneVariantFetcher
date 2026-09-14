"""Apply bounded, source-adjudicated GCK endpoint corrections to frozen keys.

This is an analysis ledger, not a mutation of the source extraction databases.
The primary excludes PMID 36208030's diabetes-only A/U split. A broader endpoint
sensitivity retains those explicitly flagged diabetes-diagnosis observations.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
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
    "GCK_protocol/observations.csv"
)

# Each key is an exact frozen variant x paper observation, not a gene-wide or
# mechanism-only blacklist. Clinical MODY observations for activating enzymes
# such as V62M and N180D are intentionally retained.
DECISIONS = {
    ("R397L", "15644838"): (
        "Case report, infant IV:1 and parents III:1/III:2; Abstract",
        "Replace frozen A1/U2: exclude the homozygous neonatal-diabetes infant; restore the two heterozygous parents as MODY A2/U0. Their measured fasting glucose is 6.1 and 5.8 mmol/L and source explicitly says GCK-MODY.",
    ),
    ("N254H", "26587058"): (
        "Results, subjects III-1/III-3, their mothers and Figure 1",
        "Replace frozen A4/U0: two compound-heterozygous neonatal-diabetes infants leave the MODY endpoint; restore the two heterozygous mothers with measured mild fasting hyperglycemia as A2/U0.",
    ),
    ("H424Y", "39610513"): (
        "Abstract, Case presentation",
        "Quarantine A1/U2: the proband is homozygous with neonatal hyperglycemia. Both parents are heterozygous, but their individual mild-hyperglycemia status is not resolved by the available source; do not call them unaffected.",
    ),
    ("T209M", "25555642"): (
        "Results, two siblings and mosaic father; Figure 1",
        "Retain two constitutionally heterozygous siblings with GCK-MODY as A2/U0. Exclude U1 for their father because his variant is low-level mosaic in blood, not a constitutional heterozygous carrier; source gives fasting glucose 69 mg/dL.",
    ),
    ("V455L", "32928245"): (
        "Case presentation G.V.; Table 2 patient 1",
        "The V455L carrier has measured hypoglycemia and hyperinsulinism treated with diazoxide, not MODY.",
    ),
    ("V455L", "34532767"): (
        "Results, family III-4/II-2/II-4/II-5; Figure 1",
        "Four HI-ascertained carriers with persistent low fasting glucose. One later develops type 2 diabetes; the four affected observations are not MODY cases and no MODY-unaffected status is inferred.",
    ),
    ("A456V", "11916951"): (
        "Abstract, proband and mother",
        "Both carriers have biochemical hypoglycemia; mother lacks symptoms, not hypoglycemia. Remove frozen A1 and U1 from the MODY endpoint.",
    ),
    ("Y214C", "15277402"): (
        "Abstract, baby girl with de novo Y214C",
        "Hypoglycemic seizures and activating glucokinase, not a MODY case.",
    ),
    ("V452L", "19053014"): (
        "Clinical presentation and mutation analysis",
        "The index child has glucokinase hyperinsulinism and an activating V452L variant.",
    ),
    ("M197I", "19336674"): (
        "Results, child 3; Abstract",
        "Child 3 has glucokinase hyperinsulinism with hypoglycemia.",
    ),
    ("W99L", "19336674"): (
        "Results, child 2; Abstract",
        "Child 2 has glucokinase hyperinsulinism with hypoglycemia.",
    ),
    ("T103S", "21454522"): (
        "Results, Novel GCK Mutations in Two Families with HH; Figure 1B",
        "HH family 2; source includes a grandmother with late-onset type 2 diabetes. The frozen A3 is not a source-resolved MODY count. Do not infer three unique HH individuals or add a MODY case from that grandmother.",
    ),
    ("V389L", "21454522"): (
        "Results, Novel GCK Mutations in Two Families with HH; Figure 1A",
        "Proband and affected father in HH family 1; two hypoglycemia observations.",
    ),
    ("I211F", "23274908"): (
        "Results, GCK mutation in case 6; Figure 7",
        "Somatic I211F in pathological pancreatic islets, absent in blood; congenital hyperinsulinism, not germline MODY.",
    ),
    ("V389L", "24890200"): (
        "Patient descriptions I.1, II.2, III.2, III.7 and proband; Figure 2",
        "Five frozen affected observations comprise four biochemical HH observations and one unassessed carrier. Remove all from MODY; unknown is not unaffected.",
    ),
    ("V455M", "26399329"): (
        "Case report",
        "Drug-resistant congenital hyperinsulinism in the child. Later follow-up may repeat this person; no unique-person claim.",
    ),
    ("V452L", "27802864"): (
        "Case presentation",
        "Adult case with documented recurrent hypoglycemia and activating Val452Leu.",
    ),
    ("M197T", "28163940"): (
        "Case reports, grandmother, son and grandchild",
        "Three biochemically hypoglycemic carriers, including asymptomatic adults.",
    ),
    ("W99C", "28247534"): (
        "Abstract, Results",
        "The novel Trp99Cys variant occurs in a hyperinsulinism case and is activating.",
    ),
    ("A456V", "28458896"): (
        "Case presentation, 34-year-old woman",
        "Familial hyperinsulinemic hypoglycemia in the mother; newborn without the allele is not a carrier denominator.",
    ),
    ("W99R", "30352420"): (
        "Results and Table 2, case 29",
        "Congenital hyperinsulinemia; c.295T>C in Results/Table 2 conflicts with c.295C>T in Discussion and frozen annotation. No nucleotide identity repair here.",
    ),
    ("T65I", "34184638"): (
        "Case presentation and genetic investigation",
        "51-year-old woman has GCK-related hyperinsulinemic hypoglycemia, not MODY.",
    ),
    ("V455M", "34635134"): (
        "Table 1, patient 1",
        "Drug-unresponsive hyperinsulinemic hypoglycemia. Possible follow-up overlap with 26399329 remains a person-ownership issue.",
    ),
    ("W99R", "34635134"): (
        "Table 1, patients 2 and 3",
        "Both patients have drug-unresponsive hyperinsulinemic hypoglycemia.",
    ),
    ("W99R", "34680961"): (
        "Abstract and family results, Table 1",
        "Four relatives have hyperinsulinemic hypoglycemia; literature-review cases are not added.",
    ),
    ("E67V", "37725835"): (
        "Table 1, patients 1 and 2",
        "Both carriers are GCK hyperinsulinism cases.",
    ),
    ("S64P", "37725835"): (
        "Table 1, patient 6; Results, Genetics",
        "Mosaic pancreatic variant in a hyperinsulinism case.",
    ),
    ("V91L", "37725835"): (
        "Table 1, patient 3; Figure 2 and Results, Genetics",
        "Proband with HI inherited V91L from an affected father previously treated for HI.",
    ),
    ("W99C", "37725835"): (
        "Table 1, patient 8",
        "Hyperinsulinism with persistent hypoglycemia after near-total pancreatectomy.",
    ),
    ("Y215C", "37725835"): (
        "Table 1, patients 4 and 5",
        "Both carriers have hyperinsulinism with low plasma glucose.",
    ),
    ("I211F", "38033998"): (
        "Table 1, patients 6 and 12; GCK variant results",
        "Two CHI patients, one heterozygous and one mosaic, not MODY cases.",
    ),
    ("A456V", "38963811"): (
        "Table 3 and Results, GCK",
        "One CHI proband with GCK A456V and an ABCC8 VUS; hypoglycemia began on day 1.",
    ),
    ("M197V", "40302972"): (
        "Patient genotype/phenotype table, patient 118; CHI inclusion criteria",
        "The child belongs to a monogenic congenital hyperinsulinism cohort.",
    ),
    ("R447L", "42184599"): (
        "Supplementary Table 5, HI GCK Arg447Leu row",
        "Low-level mosaic call in the CHI arm. ddPCR confirmation caveats do not supply a germline MODY observation.",
    ),
    ("Y214C", "42184599"): (
        "Supplementary Table 5, HI GCK Tyr214Cys row",
        "Low-level mosaic call in the CHI arm, not germline MODY.",
    ),
}

RESTORED = {
    ("R397L", "15644838"): (2, 0),
    ("N254H", "26587058"): (2, 0),
    ("T209M", "25555642"): (2, 0),
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame: pd.DataFrame, name: str) -> None:
    path = OUT / name
    payload = frame.to_csv(index=False, lineterminator="\n").encode()
    if path.suffix == ".gz":
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, mode="wb", filename="", mtime=0) as f:
            f.write(payload)
        payload = buffer.getvalue()
    path.write_bytes(payload)
    assert len(payload) < 1_200_000, path


def source_path(pmid: str) -> Path:
    paths = sorted(
        (ROOT / "results/grant_e2e_20260909/shards").glob(
            f"GCK_*/GCK/*/pmc_fulltext/{pmid}_CLEANED.md"
        )
    )
    if paths:
        return paths[0]
    path = ROOT / f"corpus/GCK/{pmid}/{pmid}_CLEANED.md"
    assert path.is_file(), pmid
    return path


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    literature = pd.read_csv(IDENTITY)
    literature = literature[literature.gene.eq("GCK")].copy()
    observations = pd.read_csv(OBS, dtype={"pmid": str})
    assert not observations.duplicated(["key", "pmid"]).any()
    assert not literature.key.duplicated().any()
    counts = observations.groupby("key")[["affected", "unaffected"]].sum()
    for row in literature.itertuples():
        assert counts.loc[row.key, "affected"] == row.affected_literature
        assert counts.loc[row.key, "unaffected"] == row.unaffected_literature

    decisions = []
    pins = {str(IDENTITY): sha(IDENTITY), str(OBS): sha(OBS)}
    for (key, pmid), (location, rationale) in DECISIONS.items():
        selected = observations[observations.key.eq(key) & observations.pmid.eq(pmid)]
        assert len(selected) == 1, (key, pmid)
        row = selected.iloc[0].to_dict()
        source = source_path(pmid)
        pins[str(source)] = sha(source)
        row.update(
            {
                "gene": "GCK",
                "decision": "quarantine_non_MODY_endpoint",
                "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "source_path": str(source.relative_to(ROOT)),
                "source_sha256": sha(source),
                "source_location": location,
                "rationale": rationale,
                "affected_removed": row["affected"],
                "unaffected_removed": row["unaffected"],
                "affected_restored": RESTORED.get((key, pmid), (0, 0))[0],
                "unaffected_restored": RESTORED.get((key, pmid), (0, 0))[1],
                "unknown_converted_to_unaffected": False,
            }
        )
        decisions.append(row)
    delta = pd.DataFrame(decisions)
    save(delta, "GCK_observation_corrections.csv")
    sums = delta.groupby("key")[
        [
            "affected_removed",
            "unaffected_removed",
            "affected_restored",
            "unaffected_restored",
        ]
    ].sum()
    literature["old_affected_literature"] = literature.affected_literature
    literature["old_unaffected_literature"] = literature.unaffected_literature
    for col in [
        "affected_removed",
        "unaffected_removed",
        "affected_restored",
        "unaffected_restored",
    ]:
        literature[col] = literature.key.map(sums[col]).fillna(0)
    literature["affected_literature"] -= literature.affected_removed
    literature["unaffected_literature"] -= literature.unaffected_removed
    literature["affected_literature"] += literature.affected_restored
    literature["unaffected_literature"] += literature.unaffected_restored
    literature["n_literature"] = (
        literature.affected_literature + literature.unaffected_literature
    )
    literature["new_affected_literature"] = literature.affected_literature
    literature["new_unaffected_literature"] = literature.unaffected_literature
    literature["zero_evidence_after_removal"] = literature.n_literature.eq(0)
    literature["review_status"] = "retained_not_fully_source_adjudicated"
    changed = (literature.affected_removed + literature.unaffected_removed).gt(0)
    literature.loc[changed, "review_status"] = "wrong_endpoint_observations_removed"
    diagnosis = observations[observations.pmid.eq("36208030")].copy()
    assert len(diagnosis) == 15
    literature["diagnosis_only_affected"] = literature.key.map(
        diagnosis.set_index("key").affected
    ).fillna(0)
    literature["diagnosis_only_unaffected"] = literature.key.map(
        diagnosis.set_index("key").unaffected
    ).fillna(0)
    literature["has_diabetes_only_endpoint_pending"] = literature.key.isin(
        diagnosis.key
    )
    literature["strict_endpoint_affected_literature"] = (
        literature.affected_literature - literature.diagnosis_only_affected
    )
    literature["strict_endpoint_unaffected_literature"] = (
        literature.unaffected_literature - literature.diagnosis_only_unaffected
    )
    literature["canonical_identity_repaired"] = False
    assert (
        literature[
            [
                "affected_literature",
                "unaffected_literature",
                "strict_endpoint_affected_literature",
                "strict_endpoint_unaffected_literature",
            ]
        ]
        .ge(0)
        .all()
        .all()
    )
    save(literature, "clinical_input_diabetes_proxy.csv.gz")
    save(literature[changed], "GCK_changed_clinical_keys.csv")
    save(
        literature[
            literature.vclass.eq("missense")
            & ~literature.canonical_wt_status.eq("match")
        ],
        "GCK_canonical_identity_queue.csv",
    )
    diagnosis["decision"] = "quarantine_A_and_U_as_unresolved_MODY"
    diagnosis["broad_endpoint_sensitivity"] = "retain_flagged_diabetes_proxy"
    diagnosis["source_location"] = (
        "Methods diabetes phenotypes and Supplementary variant table"
    )
    diagnosis["source_url"] = "https://pmc.ncbi.nlm.nih.gov/articles/PMC9659663/"
    diagnosis["rationale"] = (
        "Diagnosis-code A/U split is not a measured mild-hyperglycemia split. "
        "Do not interpret no diagnosed diabetes as absence of GCK hyperglycemia."
    )
    save(diagnosis, "GCK_diabetes_endpoint_pending.csv")
    primary = literature.copy()
    primary["affected_literature"] = primary.strict_endpoint_affected_literature
    primary["unaffected_literature"] = primary.strict_endpoint_unaffected_literature
    primary["new_affected_literature"] = primary.affected_literature
    primary["new_unaffected_literature"] = primary.unaffected_literature
    primary["affected_removed"] += primary.diagnosis_only_affected
    primary["unaffected_removed"] += primary.diagnosis_only_unaffected
    primary["n_literature"] = (
        primary.affected_literature + primary.unaffected_literature
    )
    primary["zero_evidence_after_removal"] = primary.n_literature.eq(0)
    primary.loc[primary.has_diabetes_only_endpoint_pending, "review_status"] = (
        "diagnosis_proxy_quarantined_from_MODY_primary"
    )
    save(primary, "clinical_input.csv.gz")
    for pmid in ["36208030", "36257325", "15677479", "34496959", "28331372"]:
        source = source_path(pmid)
        pins[str(source)] = sha(source)
    eligible = literature.vclass.eq("missense") & literature.canonical_wt_status.eq(
        "match"
    )
    summary = {
        "scope": "bounded source-specific GCK MODY endpoint correction",
        "clinical_keys": len(literature),
        "observation_rows_adjudicated": len(delta),
        "affected_observations_removed": float(delta.affected_removed.sum()),
        "unaffected_observations_removed": float(delta.unaffected_removed.sum()),
        "affected_observations_restored": float(delta.affected_restored.sum()),
        "unaffected_observations_restored": float(delta.unaffected_restored.sum()),
        "changed_keys": int(changed.sum()),
        "keys_without_remaining_clinical_evidence": int(
            literature.zero_evidence_after_removal.sum()
        ),
        "eligible_missense_old_A": float(
            literature.loc[eligible, "old_affected_literature"].sum()
        ),
        "eligible_missense_corrected_A": float(
            literature.loc[eligible, "affected_literature"].sum()
        ),
        "eligible_missense_old_U": float(
            literature.loc[eligible, "old_unaffected_literature"].sum()
        ),
        "eligible_missense_corrected_U": float(
            literature.loc[eligible, "unaffected_literature"].sum()
        ),
        "primary_eligible_missense_A": float(
            primary.loc[eligible, "affected_literature"].sum()
        ),
        "primary_eligible_missense_U": float(
            primary.loc[eligible, "unaffected_literature"].sum()
        ),
        "primary_keys_without_clinical_evidence": int(
            primary.zero_evidence_after_removal.sum()
        ),
        "diagnosis_only_pending_A": float(diagnosis.affected.sum()),
        "diagnosis_only_pending_U": float(diagnosis.unaffected.sum()),
        "pmid_36257325_usable_frozen_observations": int(
            observations.pmid.eq("36257325").sum()
        ),
        "unknowns_converted_to_unaffected": 0,
        "canonical_identity_repairs": 0,
        "full_source_or_person_ownership_adjudication": False,
        "gnomad_counts_edited": False,
        "source_databases_edited": False,
        "script_sha256": sha(Path(__file__)),
        "source_sha256": pins,
    }
    (OUT / "source_checks.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in summary.items() if k != "source_sha256"}, indent=2)
    )


if __name__ == "__main__":
    main()
