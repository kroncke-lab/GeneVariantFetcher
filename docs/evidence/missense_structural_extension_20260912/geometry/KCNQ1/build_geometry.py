"""Acquire human KCNQ1 cardiac biological assemblies and canonical geometry.

Run with ProteinProximityAnalysis/.venv/bin/python. Public downloads are cached
under GVF results/; only compact geometry and provenance are durable artifacts.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
RAW = REPO / "results/missense_structural_extension_20260912/raw/KCNQ1"
PDB_IDS = ("9U7F", "9UC8")
TARGET_CHAINS = ("A", "D", "G", "J")
PRIMARY_STUDY = "https://www.nature.com/articles/s41422-025-01182-9"
ATOMIC_MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
SIDECHAIN_ATOMS = {
    "A": "CB",
    "R": "CB CG CD NE CZ NH1 NH2",
    "N": "CB CG OD1 ND2",
    "D": "CB CG OD1 OD2",
    "C": "CB SG",
    "Q": "CB CG CD OE1 NE2",
    "E": "CB CG CD OE1 OE2",
    "G": "",
    "H": "CB CG ND1 CD2 CE1 NE2",
    "I": "CB CG1 CG2 CD1",
    "L": "CB CG CD1 CD2",
    "K": "CB CG CD CE NZ",
    "M": "CB CG SD CE",
    "F": "CB CG CD1 CD2 CE1 CE2 CZ",
    "P": "CB CG CD",
    "S": "CB OG",
    "T": "CB OG1 CG2",
    "W": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "Y": "CB CG CD1 CD2 CE1 CE2 CZ OH",
    "V": "CB CG1 CG2",
}
SOURCES = {
    "P51787.json": "https://rest.uniprot.org/uniprotkb/P51787.json",
    "alphafold_api.json": "https://alphafold.ebi.ac.uk/api/prediction/P51787",
}
for _pdb in PDB_IDS:
    SOURCES[f"{_pdb}-assembly1.cif.gz"] = (
        f"https://files.rcsb.org/download/{_pdb.lower()}-assembly1.cif.gz"
    )
    SOURCES[f"{_pdb}_sifts.json"] = (
        f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{_pdb.lower()}"
    )
    SOURCES[f"{_pdb}_entry.json"] = f"https://data.rcsb.org/rest/v1/core/entry/{_pdb}"
    SOURCES[f"{_pdb}_assembly.json"] = (
        f"https://data.rcsb.org/rest/v1/core/assembly/{_pdb}/1"
    )
    for _entity in (1, 2, 3):
        SOURCES[f"{_pdb}_entity{_entity}.json"] = (
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{_pdb}/{_entity}"
        )


def get_one(item):
    name, url = item
    path = RAW / name
    if not path.exists():
        response = requests.get(url, timeout=90)
        response.raise_for_status()
        path.write_bytes(response.content)
    return name


def acquire():
    RAW.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=6) as pool:
        for name in pool.map(get_one, SOURCES.items()):
            print(name, flush=True)
    af = json.loads((RAW / "alphafold_api.json").read_text())
    candidates = [x for x in af if x["uniprotAccession"] == "P51787"]
    if len(candidates) != 1:
        raise ValueError("Expected one canonical KCNQ1 AlphaFold model")
    SOURCES["AF_P51787.cif"] = candidates[0]["cifUrl"]
    get_one(("AF_P51787.cif", SOURCES["AF_P51787.cif"]))


def read_cif(name):
    path = RAW / name
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        table = MMCIF2Dict(handle)
    with opener(path, "rt") as handle:
        structure = MMCIFParser(
            QUIET=True, auth_chains=False, auth_residues=False
        ).get_structure(name, handle)
    if len(structure) != 1:
        raise ValueError("Expected one model of the actual biological assembly")
    return table, structure[0]


def identity_mapping(table):
    mapping = {}
    fields = (
        "label_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_ins_code",
    )
    for chain, label, author, author_chain, insertion in zip(
        *(table[f"_atom_site.{field}"] for field in fields), strict=True
    ):
        if label in (".", "?"):
            continue
        value = (author, author_chain, "" if insertion in (".", "?") else insertion)
        key = chain, int(label)
        if key in mapping and mapping[key] != value:
            raise ValueError("Conflicting author-to-label residue mapping")
        mapping[key] = value
    return mapping


def coordinates(residue):
    aa = seq1(residue.resname)
    ca = np.asarray(residue["CA"].coord, dtype=float)
    sidechain = [a for a in residue if a.name not in BACKBONE and a.element != "H"]
    expected = set(SIDECHAIN_ATOMS[aa].split())
    missing = sorted(expected - {a.name for a in sidechain})
    extra = sorted({a.name for a in sidechain} - expected)
    if extra:
        raise ValueError(f"Unexpected sidechain atoms: {extra}")
    if missing:
        return None, ca, missing
    com = (
        ca
        if aa == "G"
        else np.average(
            np.asarray([a.coord for a in sidechain], dtype=float),
            axis=0,
            weights=[ATOMIC_MASS[a.element] for a in sidechain],
        )
    )
    if not np.isfinite(com).all() or not np.isfinite(ca).all():
        raise ValueError("Nonfinite coordinates")
    return com, ca, []


def write_rows(name, rows):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def classify_segments(rows, pdb, chain):
    start = None
    for index in range(len(rows) + 1):
        eligible = index < len(rows) and rows[index]["geometry_state"] == "idr"
        if eligible and start is None:
            start = index
        elif not eligible and start is not None:
            segment = f"{pdb}_{chain}_AF_lt50_{start + 1}_{index}"
            for row in rows[start:index]:
                row["idr_segment"] = segment
            start = None


def build():
    old_report = HERE / "geometry_identity_report.json"
    if old_report.exists():
        for source in json.loads(old_report.read_text())["sources"]:
            if (
                hashlib.sha256(
                    (REPO / source["local_raw_file"]).read_bytes()
                ).hexdigest()
                != source["sha256"]
            ):
                raise ValueError("Cached source bytes changed after geometry freeze")
    uniprot = json.loads((RAW / "P51787.json").read_text())
    canonical = uniprot["sequence"]["value"]
    if len(canonical) != 676:
        raise ValueError("Expected canonical 676-aa P51787 sequence")
    (HERE / "P51787-1.fasta").write_text(
        ">P51787-1 canonical KCNQ1; UniProt displayed sequence\n" + canonical + "\n"
    )
    afmeta = json.loads((RAW / "alphafold_api.json").read_text())[0]
    if afmeta["sequence"] != canonical or afmeta["latestVersion"] != 6:
        raise ValueError("Dated analysis requires exact canonical AlphaFold v6")
    _, af_model = read_cif("AF_P51787.cif")
    af_chains = list(af_model.get_chains())
    if len(af_chains) != 1:
        raise ValueError("AlphaFold confidence reference must be a single monomer")
    af_residues = {r.id[1]: r for r in af_chains[0] if r.id[0] == " "}
    if (
        set(af_residues) != set(range(1, 677))
        or "".join(seq1(af_residues[p].resname) for p in range(1, 677)) != canonical
    ):
        raise ValueError("AlphaFold residue numbers/WT sequence are not canonical")
    confidence = {p: float(r["CA"].bfactor) for p, r in af_residues.items()}
    frames, qc_rows = {}, []
    incomplete_rows = []
    for pdb in PDB_IDS:
        assembly = json.loads((RAW / f"{pdb}_assembly.json").read_text())
        entry = json.loads((RAW / f"{pdb}_entry.json").read_text())
        entities = {
            e: json.loads((RAW / f"{pdb}_entity{e}.json").read_text())
            for e in (1, 2, 3)
        }
        if (
            assembly["rcsb_assembly_info"]["polymer_entity_instance_count_protein"]
            != 12
        ):
            raise ValueError("Expected four KCNQ1, four KCNE1 and four CaM chains")
        target = entities[1]
        if (
            target["entity_poly"]["pdbx_seq_one_letter_code_can"] != canonical
            or target["entity_poly"]["rcsb_mutation_count"] != 0
        ):
            raise ValueError("KCNQ1 deposited reference sequence has changed")
        alignments = [
            a
            for a in target["rcsb_polymer_entity_align"]
            if a["reference_database_accession"] == "P51787"
            and a["provenance_source"] == "SIFTS"
        ]
        if len(alignments) != 1 or alignments[0]["aligned_regions"] != [
            {"entity_beg_seq_id": 1, "length": 676, "ref_beg_seq_id": 1}
        ]:
            raise ValueError("Expected complete 1:1 SIFTS sequence mapping")
        pdbe = json.loads((RAW / f"{pdb}_sifts.json").read_text())[pdb.lower()][
            "UniProt"
        ]["P51787"]["mappings"]
        if {m["struct_asym_id"] for m in pdbe} != set(TARGET_CHAINS):
            raise ValueError(
                "PDBe target chain identity differs from expected tetramer"
            )
        for mapping in pdbe:
            if (
                mapping["unp_start"],
                mapping["unp_end"],
                mapping["start"]["residue_number"],
                mapping["end"]["residue_number"],
                mapping["identity"],
            ) != (1, 676, 1, 676, 1.0):
                raise ValueError("PDBe numbering or WT mapping changed")
        table, model = read_cif(f"{pdb}-assembly1.cif.gz")
        asym_entity = dict(
            zip(table["_struct_asym.id"], table["_struct_asym.entity_id"], strict=True)
        )
        if {c for c, e in asym_entity.items() if e == "1"} != set(TARGET_CHAINS):
            raise ValueError(
                "Assembly file target chains differ from entity/SIFTS records"
            )
        identity = identity_mapping(table)
        partner_metadata = []
        for entity_id, accession in ((1, "P51787"), (2, "P0DP23"), (3, "P15382")):
            entity = entities[entity_id]
            ids = entity["rcsb_polymer_entity_container_identifiers"]
            if ids["uniprot_ids"] != [accession] or len(ids["asym_ids"]) != 4:
                raise ValueError(
                    "Biological assembly partner identity/stoichiometry mismatch"
                )
            source_taxa = {
                s["ncbi_taxonomy_id"] for s in entity["rcsb_entity_source_organism"]
            }
            if source_taxa != {9606}:
                raise ValueError("Expected human KCNQ1, KCNE1 and CaM proteins")
            seq = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
            observed_counts = {}
            for chain in ids["asym_ids"]:
                residues = [r for r in model[chain] if r.id[0] == " "]
                observed_counts[chain] = len(residues)
                if any(seq1(r.resname) != seq[r.id[1] - 1] for r in residues):
                    raise ValueError(
                        "Assembly coordinate WT disagrees with entity sequence"
                    )
            partner_metadata.append(
                {
                    "entity_id": entity_id,
                    "uniprot_accession": accession,
                    "name": entity["rcsb_polymer_entity"]["pdbx_description"],
                    "chains": ids["asym_ids"],
                    "deposited_reference_length": len(seq),
                    "observed_residues_by_chain": observed_counts,
                    "variant_donor_target": entity_id == 1,
                    "rcsb_sifts_alignment": entity["rcsb_polymer_entity_align"],
                }
            )
        all_rows = []
        for chain in TARGET_CHAINS:
            residues = {r.id[1]: r for r in model[chain] if r.id[0] == " "}
            rows = []
            for pos, aa in enumerate(canonical, 1):
                row = {
                    "frame_id": pdb,
                    "chain": chain,
                    "canonical_pos": pos,
                    "aa": aa,
                    "aa_ref": aa,
                    "structure_aa": "",
                    "author_pos": "",
                    "author_chain": "",
                    "insertion_code": "",
                    "label_seq_id": pos,
                    "observed": pos in residues,
                    "geometry_state": "missing",
                    "idr_segment": "",
                    "geometry_state_source": "unresolved_or_outside_construct",
                    "com_x": "",
                    "com_y": "",
                    "com_z": "",
                    "ca_x": "",
                    "ca_y": "",
                    "ca_z": "",
                    "plddt": confidence[pos],
                    "plddt_source": "AF-P51787-F1-model_v6",
                    "experimental_b_factor": "",
                    "coordinate_source": "none",
                    "mapping_method": "SIFTS_RCSB_PDBe_1to1_exact_sequence_WT",
                    "coordinate_note": "no_AF_coordinate_hybrid",
                    "in_experimental_construct": 76 <= pos <= 620,
                    "ca_geometry_state": "missing",
                }
                if pos in residues:
                    residue = residues[pos]
                    if (
                        seq1(residue.resname) != aa
                        or not row["in_experimental_construct"]
                    ):
                        raise ValueError("Observed WT or construct numbering mismatch")
                    row["structure_aa"] = seq1(residue.resname)
                    row["author_pos"], row["author_chain"], row["insertion_code"] = (
                        identity[(chain, pos)]
                    )
                    row["experimental_b_factor"] = float(residue["CA"].bfactor)
                    com, ca, missing = coordinates(residue)
                    row["ca_geometry_state"] = "structured"
                    for axis, value in zip("xyz", ca, strict=True):
                        row[f"ca_{axis}"] = round(float(value), 8)
                    row["coordinate_source"] = pdb
                    if com is None:
                        row["geometry_state_source"] = (
                            "incomplete_experimental_sidechain_COM_unavailable"
                        )
                        row["coordinate_note"] = (
                            "incomplete_sidechain;CA_retained_for_explicit_metric_sensitivity"
                        )
                        incomplete_rows.append(
                            {
                                "frame_id": pdb,
                                "chain": chain,
                                "canonical_pos": pos,
                                "aa_ref": aa,
                                "missing_atoms": " ".join(missing),
                            }
                        )
                    else:
                        row["geometry_state"] = "structured"
                        row["geometry_state_source"] = (
                            "observed_experimental_coordinates"
                        )
                        row["coordinate_note"] = (
                            "glycine_CA_fallback"
                            if aa == "G"
                            else "complete_sidechain_heavy_atom_COM"
                        )
                        for axis, value in zip("xyz", com, strict=True):
                            row[f"com_{axis}"] = round(float(value), 8)
                elif confidence[pos] < 50:
                    row["geometry_state"] = "idr"
                    row["ca_geometry_state"] = "idr"
                    row["geometry_state_source"] = (
                        "independent_AF_plddt_lt50_candidate_IDR"
                    )
                    row["coordinate_note"] = (
                        "polymer_only;no_experimental_or_AF_coordinates"
                    )
                elif confidence[pos] < 70:
                    row["geometry_state"] = "ambiguous"
                    row["ca_geometry_state"] = "ambiguous"
                    row["geometry_state_source"] = (
                        "independent_AF_plddt50_to70;geometry_unavailable"
                    )
                rows.append(row)
            classify_segments(rows, pdb, chain)
            counts = dict(Counter(r["geometry_state"] for r in rows))
            qc_rows.append(
                {
                    "frame_id": pdb,
                    "chain": chain,
                    "canonical_rows": len(rows),
                    "observed_CA": len(residues),
                    **{
                        state: counts.get(state, 0)
                        for state in ("structured", "idr", "ambiguous", "missing")
                    },
                }
            )
            all_rows.extend(rows)
        if len({(r["chain"], r["canonical_pos"]) for r in all_rows}) != 2704:
            raise ValueError("Canonical row grain/coverage failed")
        write_rows(f"{pdb}_canonical_geometry.csv", all_rows)
        frames[pdb] = {
            "assembly_id": "1",
            "source_type": "experimental_biological_assembly",
            "status": "actual_KCNQ1_tetramer_with_partial_resolved_coverage",
            "target_accession": "P51787",
            "target_chains": list(TARGET_CHAINS),
            "assembly_protein_stoichiometry": "KCNQ1_4:KCNE1_4:CaM_4",
            "state": "apo_condition_closed" if pdb == "9U7F" else "PIP2_added_open",
            "resolution_angstrom": entry["rcsb_entry_info"]["resolution_combined"],
            "canonical_rows": 2704,
            "full_length_reference_sequence_identity": 1.0,
            "observed_coordinate_WT_mismatches": 0,
            "entities": partner_metadata,
            "assembly_record": assembly["pdbx_struct_assembly"],
            "construct": {
                "retained_KCNQ1_range": [76, 620],
                "type": "KCNE1_linker_fusion_to_N_C_truncated_KCNQ1",
                "source": PRIMARY_STUDY,
                "note": "PDB full canonical SEQRES does not establish full-length experimental construct",
            },
            "per_chain_geometry": [q for q in qc_rows if q["frame_id"] == pdb],
            "candidate_IDR_segments": sorted(
                {r["idr_segment"] for r in all_rows if r["idr_segment"]}
            ),
        }
    write_rows("geometry_summary.csv", qc_rows)
    write_rows("incomplete_sidechains.csv", incomplete_rows)
    report = {
        "gene": "KCNQ1",
        "canonical_accession": "P51787-1",
        "canonical_length": 676,
        "canonical_sequence_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "primary_frame": "9U7F",
        "sensitivity_frame": "9UC8",
        "frames": frames,
        "metric": "complete_sidechain_heavy_atom_mass_weighted_COM;glycine_CA",
        "source_scope": "Actual assembly coordinates only; AF supplies confidence without coordinates",
        "disorder_rule": "Absent experimental position and independently mapped AF pLDDT<50; contiguous per-chain segments only; candidate not experimental proof",
        "assembly_partner_scope": "Full partner coordinates retained in raw assemblies; canonical geometry rows and density donors restricted to KCNQ1 chains",
        "primary_study": PRIMARY_STUDY,
        "sources": [
            {
                "name": name,
                "url": url,
                "local_raw_file": str((RAW / name).relative_to(REPO)),
                "sha256": hashlib.sha256((RAW / name).read_bytes()).hexdigest(),
                "bytes": (RAW / name).stat().st_size,
            }
            for name, url in sorted(SOURCES.items())
        ],
        "checks": {
            "exact_full_canonical_sequence": True,
            "independent_RCSB_PDBe_SIFTS_agreement": True,
            "all_target_copy_WT_checks": True,
            "partner_identity_and_stoichiometry": True,
            "full_676_rows_per_target_copy": True,
            "no_AF_coordinates_spliced_into_assembly": True,
            "missing_sidechain_COM_withheld": True,
        },
    }
    (HERE / "geometry_identity_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(qc_rows, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-only", action="store_true")
    args = parser.parse_args()
    acquire()
    if args.download_only:
        return
    build()


if __name__ == "__main__":
    main()
