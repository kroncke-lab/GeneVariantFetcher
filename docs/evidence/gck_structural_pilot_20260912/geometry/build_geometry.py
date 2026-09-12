"""Fetch public GCK structures and build canonical, full-length geometry manifests.

Run with the ProteinProximityAnalysis virtualenv; raw downloads are cached under
GVF results/. No experimental/AlphaFold coordinate frames are combined.
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

import requests
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RAW = REPO / "results/gck_structural_pilot_20260912/raw"
SOURCES = {
    "P35557.json": "https://rest.uniprot.org/uniprotkb/P35557.json",
    "P35557-2.fasta": "https://rest.uniprot.org/uniprotkb/P35557-2.fasta",
    "alphafold_api.json": "https://alphafold.ebi.ac.uk/api/prediction/P35557",
}
ATOMIC_MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
EXPECTED_SIDECHAIN_ATOMS = {
    "A": 1,
    "R": 7,
    "N": 4,
    "D": 4,
    "C": 2,
    "Q": 5,
    "E": 5,
    "G": 0,
    "H": 6,
    "I": 4,
    "L": 4,
    "K": 5,
    "M": 4,
    "F": 7,
    "P": 3,
    "S": 2,
    "T": 3,
    "W": 10,
    "Y": 8,
    "V": 3,
}
LOOP_DISORDER_SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3531626/"
for _pdb in ("1V4S", "1V4T"):
    SOURCES[f"{_pdb}-assembly1.cif.gz"] = (
        f"https://files.rcsb.org/download/{_pdb.lower()}-assembly1.cif.gz"
    )
    for _kind, _suffix in (("entry", ""), ("assembly", "/1"), ("polymer_entity", "/1")):
        SOURCES[f"{_pdb}_{_kind}.json"] = (
            f"https://data.rcsb.org/rest/v1/core/{_kind}/{_pdb}{_suffix}"
        )


def download_one(item):
    name, url = item
    target = RAW / name
    if not target.exists():
        result = requests.get(url, timeout=90)
        result.raise_for_status()
        target.write_bytes(result.content)
    return name, url


def download():
    RAW.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=6) as pool:
        for name, _ in pool.map(download_one, SOURCES.items()):
            print(name, flush=True)
    af = json.loads((RAW / "alphafold_api.json").read_text())
    candidates = [x for x in af if x["uniprotAccession"] == "P35557"]
    if len(candidates) != 1:
        raise ValueError(
            f"Expected one full GCK AlphaFold model, found {len(candidates)}"
        )
    for suffix, key in (("cif", "cifUrl"), ("pae.json", "paeDocUrl")):
        name = f"AF_P35557.{suffix}"
        SOURCES[name] = candidates[0][key]
        download_one((name, SOURCES[name]))


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
        raise ValueError("Expected a single structural model")
    chains = list(structure[0].get_chains())
    if len(chains) != 1 or chains[0].id != "A":
        raise ValueError(
            f"Expected one protein chain A, found {[c.id for c in chains]}"
        )
    return table, {r.id[1]: r for r in chains[0] if r.id[0] == " "}


def coordinates(residue):
    aa = seq1(residue.resname)
    atoms = [a for a in residue if a.name not in BACKBONE and a.element != "H"]
    if len(atoms) != EXPECTED_SIDECHAIN_ATOMS[aa]:
        raise ValueError(f"Incomplete sidechain: {residue.id} {aa}")
    ca = np.asarray(residue["CA"].coord, dtype=float)
    if aa == "G":
        com = ca
    else:
        com = np.average(
            np.asarray([a.coord for a in atoms], dtype=float),
            axis=0,
            weights=[ATOMIC_MASS[a.element] for a in atoms],
        )
    return com, ca


def atom_identity(table):
    fields = [
        "label_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_ins_code",
    ]
    mapping = {}
    for chain, label, auth, auth_chain, insertion in zip(
        *(table[f"_atom_site.{field}"] for field in fields), strict=True
    ):
        if label in (".", "?"):
            continue
        value = (auth, auth_chain, "" if insertion in (".", "?") else insertion)
        key = (chain, int(label))
        if key in mapping and mapping[key] != value:
            raise ValueError("Ambiguous author/label residue identity")
        mapping[key] = value
    return mapping


def idr_segments(plddt):
    segments = {}
    start = None
    for pos in range(1, len(plddt) + 2):
        is_idr = pos in plddt and plddt[pos] < 50
        if is_idr and start is None:
            start = pos
        if not is_idr and start is not None:
            for j in range(start, pos):
                segments[j] = f"AF_lt50_{start}_{pos - 1}"
            start = None
    return segments


def empty_row(frame, pos, aa, confidence, segments):
    return {
        "frame_id": frame,
        "chain": "A",
        "canonical_pos": pos,
        "aa": aa,
        "aa_ref": aa,
        "structure_aa": "",
        "author_pos": "",
        "author_chain": "",
        "insertion_code": "",
        "label_seq_id": "",
        "observed": False,
        "geometry_state": "missing",
        "idr_segment": "",
        "geometry_state_source": "no_usable_coordinates",
        "com_x": "",
        "com_y": "",
        "com_z": "",
        "ca_x": "",
        "ca_y": "",
        "ca_z": "",
        "plddt": confidence[pos],
        "plddt_source": "AF-P35557-F1-model_v6",
        "experimental_b_factor": "",
        "coordinate_source": "none",
        "mapping_method": "canonical_reference_only",
        "coordinate_note": "",
    }


def set_coords(row, residue, identity, *, experimental):
    com, ca = coordinates(residue)
    row.update(
        {
            "structure_aa": seq1(residue.resname),
            "author_pos": identity[0],
            "author_chain": identity[1],
            "insertion_code": identity[2],
            "label_seq_id": residue.id[1],
            "observed": True,
        }
    )
    if row["structure_aa"] != row["aa_ref"]:
        raise ValueError(f"Coordinate WT differs from canonical: {row}")
    if row["geometry_state"] == "structured":
        for prefix, xyz in (("com", com), ("ca", ca)):
            for axis, value in zip("xyz", xyz, strict=True):
                row[f"{prefix}_{axis}"] = round(float(value), 8)
        row["coordinate_source"] = row["frame_id"]
        row["coordinate_note"] = (
            "glycine_CA_fallback"
            if row["aa_ref"] == "G"
            else "complete_sidechain_heavy_atom_mass_weighted_COM"
        )
    else:
        row["coordinate_note"] = "predicted_coordinates_withheld_below_plddt70"
    if experimental:
        row["experimental_b_factor"] = float(residue["CA"].bfactor)


def write_rows(name, rows):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build():
    frozen_report = HERE / "geometry_identity_report.json"
    if frozen_report.exists():
        for source in json.loads(frozen_report.read_text())["sources"]:
            path = REPO / source["local_raw_file"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
                raise ValueError(f"Frozen public-source bytes changed: {path.name}")
    uniprot = json.loads((RAW / "P35557.json").read_text())
    canonical = uniprot["sequence"]["value"]
    if len(canonical) != 465:
        raise ValueError("Canonical GCK sequence version/length changed")
    iso2 = "".join((RAW / "P35557-2.fasta").read_text().splitlines()[1:])
    alternative = next(
        f for f in uniprot["features"] if f.get("featureId") == "VSP_002074"
    )
    alt = alternative["alternativeSequence"]
    if alternative["location"] != {
        "start": {"value": 1, "modifier": "EXACT"},
        "end": {"value": 15, "modifier": "EXACT"},
    }:
        raise ValueError("Isoform-2 alternative-sequence mapping changed")
    if (
        canonical[:15] != alt["originalSequence"]
        or iso2 != alt["alternativeSequences"][0] + canonical[15:]
    ):
        raise ValueError("Isoform-2 sequence fails UniProt feature reconstruction")
    (HERE / "P35557-1.fasta").write_text(
        ">P35557-1 canonical GCK; UniProt displayed sequence\n" + canonical + "\n"
    )
    afmeta = next(
        a
        for a in json.loads((RAW / "alphafold_api.json").read_text())
        if a["uniprotAccession"] == "P35557"
    )
    if afmeta["latestVersion"] != 6:
        raise ValueError("This dated pilot is frozen to AlphaFold model v6")
    if afmeta["sequence"] != canonical or (
        afmeta["sequenceStart"],
        afmeta["sequenceEnd"],
    ) != (1, 465):
        raise ValueError("AlphaFold model does not cover the full canonical sequence")
    af_table, af_residues = read_cif("AF_P35557.cif")
    if set(af_residues) != set(range(1, 466)):
        raise ValueError("AlphaFold residue numbering/coverage is not canonical")
    if "".join(seq1(af_residues[p].resname) for p in range(1, 466)) != canonical:
        raise ValueError("AlphaFold coordinates fail full WT identity check")
    confidence = {
        pos: float(residue["CA"].bfactor) for pos, residue in af_residues.items()
    }
    segments = idr_segments(confidence)
    af_identity = atom_identity(af_table)
    reports = {}
    af_rows = []
    for pos, aa in enumerate(canonical, 1):
        row = empty_row("AF_P35557", pos, aa, confidence, segments)
        row["geometry_state"] = (
            "structured"
            if confidence[pos] >= 70
            else "idr"
            if confidence[pos] < 50
            else "ambiguous"
        )
        row["geometry_state_source"] = "AlphaFold_plddt_threshold"
        row["idr_segment"] = segments.get(pos, "")
        row["mapping_method"] = (
            "full_canonical_sequence_and_all_coordinate_WT_exact_match"
        )
        set_coords(row, af_residues[pos], af_identity[("A", pos)], experimental=False)
        af_rows.append(row)
    write_rows("AF_P35557_canonical_geometry.csv", af_rows)
    reports["AF_P35557"] = {
        "source_type": "predicted_full_monomer",
        "model_id": afmeta["modelEntityId"],
        "model_version": afmeta["latestVersion"],
        "residue_count": 465,
        "geometry_state_counts": dict(Counter(r["geometry_state"] for r in af_rows)),
        "candidate_idr_segments": sorted(set(segments.values())),
        "biological_unit": "monomer; consistent with UniProt SUBUNIT",
        "full_coordinate_sequence_identity": 1.0,
    }
    excluded_rows = []
    for pdb in ("1V4S", "1V4T"):
        assembly = json.loads((RAW / f"{pdb}_assembly.json").read_text())
        entity = json.loads((RAW / f"{pdb}_polymer_entity.json").read_text())
        entry = json.loads((RAW / f"{pdb}_entry.json").read_text())
        info = assembly["rcsb_assembly_info"]
        if (
            assembly["pdbx_struct_assembly"]["oligomeric_count"] != 1
            or info["polymer_entity_instance_count_protein"] != 1
        ):
            raise ValueError(f"{pdb} assembly is not an actual biological monomer")
        if entity["entity_poly"]["rcsb_mutation_count"] != 0:
            raise ValueError(f"{pdb} is a mutant structure")
        table, residues = read_cif(f"{pdb}-assembly1.cif.gz")
        protein_seq = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
        if (
            "".join(table["_entity_poly.pdbx_seq_one_letter_code_can"][0].split())
            != protein_seq
        ):
            raise ValueError("Assembly and entry polymer sequences differ")
        references = [
            a
            for a in entity["rcsb_polymer_entity_align"]
            if a["reference_database_accession"] == "P35557"
            and a["provenance_source"] == "SIFTS"
        ]
        if len(references) != 1:
            raise ValueError("Ambiguous or absent canonical SIFTS mapping")
        mapping = {}
        for region in references[0]["aligned_regions"]:
            e, c, length = (
                region[k] for k in ("entity_beg_seq_id", "ref_beg_seq_id", "length")
            )
            segment = protein_seq[e - 1 : e - 1 + length]
            if segment != canonical[c - 1 : c - 1 + length]:
                raise ValueError("SIFTS mapping fails sequence identity check")
            # Independent exact-sequence alignment must identify the same unique start.
            if canonical.find(segment) != c - 1 or canonical.count(segment) != 1:
                raise ValueError("Canonical sequence alignment is ambiguous")
            for offset in range(length):
                mapping[e + offset] = c + offset
        if set(mapping.values()) != set(range(16, 466)):
            raise ValueError("Common isoform-1/isoform-2 region changed")
        identity = atom_identity(table)
        rows = [
            empty_row(pdb, pos, aa, confidence, segments)
            for pos, aa in enumerate(canonical, 1)
        ]
        inverse_mapping = {c: e for e, c in mapping.items()}
        for row in rows:
            pos = row["canonical_pos"]
            label = inverse_mapping.get(pos)
            if label in residues:
                residue = residues[label]
                if seq1(residue.resname) != protein_seq[label - 1]:
                    raise ValueError(
                        "Observed residue differs from deposited polymer sequence"
                    )
                row["geometry_state"] = "structured"
                row["geometry_state_source"] = "observed_experimental_coordinates"
                row["mapping_method"] = (
                    "SIFTS_segment_plus_unique_exact_sequence_alignment_and_WT_check"
                )
                set_coords(row, residue, identity[("A", label)], experimental=True)
            else:
                row["label_seq_id"] = label or ""
                row["mapping_method"] = (
                    "SIFTS_sequence_mapped_unobserved"
                    if label
                    else "unmapped_isoform_specific_canonical_N_terminus"
                )
                if confidence[pos] < 50:
                    row["geometry_state"] = "idr"
                    row["geometry_state_source"] = "AlphaFold_plddt_lt50_candidate"
                    row["idr_segment"] = segments[pos]
                    row["coordinate_note"] = (
                        "missing_experimental;candidate_IDR_from_AF_plddt_lt50;polymer_only"
                    )
                elif confidence[pos] < 70:
                    row["geometry_state"] = "ambiguous"
                    row["geometry_state_source"] = "AlphaFold_plddt50_to70"
                    row["coordinate_note"] = (
                        "missing_experimental;AF_plddt50_to70;no_coordinates"
                    )
                else:
                    row["coordinate_note"] = (
                        "missing_experimental;no_AF_coordinate_hybrid"
                    )
        if pdb == "1V4T":
            # Preserve a strict missing-loop control. Primary state-specific
            # disorder follows experimental evidence rather than a different
            # conformation's high AlphaFold confidence.
            control = [dict(row, frame_id="1V4T_missing_loop_control") for row in rows]
            write_rows("1V4T_missing_loop_control.csv", control)
            if canonical[156] != "E" or canonical[178] != "N":
                raise ValueError(
                    "Published Glu157-Asn179 loop endpoints do not match canonical"
                )
            scheme = {
                int(label): (int(pdb_pos), seq1(mon))
                for label, pdb_pos, mon, chain in zip(
                    table["_pdbx_poly_seq_scheme.seq_id"],
                    table["_pdbx_poly_seq_scheme.pdb_seq_num"],
                    table["_pdbx_poly_seq_scheme.mon_id"],
                    table["_pdbx_poly_seq_scheme.asym_id"],
                    strict=True,
                )
                if chain == "A"
            }
            for pos in range(157, 180):
                row = rows[pos - 1]
                label = inverse_mapping[pos]
                if row["observed"] or scheme[label] != (pos, canonical[pos - 1]):
                    raise ValueError(
                        "Published missing-loop coordinates/numbering check fails"
                    )
                row["geometry_state"] = "idr"
                row["idr_segment"] = "1V4T_experimental_disordered_157_179"
                row["geometry_state_source"] = (
                    "experimental_state_specific_disorder_PMC3531626"
                )
                row["coordinate_note"] = (
                    "experimentally_disordered_active_site_loop;polymer_only;no_AF_coordinate_hybrid"
                )
        for label, residue in residues.items():
            if label not in mapping:
                auth = identity[("A", label)]
                excluded_rows.append(
                    {
                        "frame_id": pdb,
                        "chain": "A",
                        "label_seq_id": label,
                        "author_pos": auth[0],
                        "structure_aa": seq1(residue.resname),
                        "reason": "outside_canonical_SIFTS_region;isoform_specific_or_construct_residue",
                    }
                )
        write_rows(f"{pdb}_canonical_geometry.csv", rows)
        reports[pdb] = {
            "source_type": "experimental_biological_assembly",
            "assembly_id": "1",
            "biological_unit": assembly["pdbx_struct_assembly"]["oligomeric_details"],
            "assembly_assignment": assembly["pdbx_struct_assembly"]["details"],
            "protein_chain_count": info["polymer_entity_instance_count_protein"],
            "deposited_polymer_residues": len(protein_seq),
            "observed_polymer_residues": len(residues),
            "mapped_sequence_residues": len(mapping),
            "canonical_rows": len(rows),
            "geometry_state_counts": dict(Counter(r["geometry_state"] for r in rows)),
            "unmapped_observed_residues_excluded": sum(
                r["frame_id"] == pdb for r in excluded_rows
            ),
            "sifts_segments": references[0]["aligned_regions"],
            "mapped_WT_identity": 1.0,
            "resolution_angstrom": entry["rcsb_entry_info"]["resolution_combined"],
            "ligand_state": "glucose_and_activator_bound"
            if pdb == "1V4S"
            else "glucose_free",
            "conformation_interpretation": "active_closed"
            if pdb == "1V4S"
            else "inactive_super_open",
            "state_source": LOOP_DISORDER_SOURCE,
        }
        if pdb == "1V4T":
            reports[pdb]["experimental_disorder_segment"] = {
                "canonical_start": 157,
                "canonical_end": 179,
                "endpoint_identity": "E157-N179",
                "residues": 23,
                "source": LOOP_DISORDER_SOURCE,
                "source_sections": [
                    "In silico dynamic and conformational effects of ATP binding",
                    "Modelling of the hGK apoenzyme in the super-open conformation",
                ],
                "numbering_validation": "Each missing label maps through SIFTS to the identically numbered mmCIF pdb_seq_num and exact canonical AA; no isoform offset assumed",
                "geometry": "3.8*sqrt(canonical sequence separation) within this same-chain segment only",
                "AF_confidence_does_not_override_experimental_state_disorder": True,
                "control_file": "1V4T_missing_loop_control.csv",
            }
    write_rows("excluded_noncanonical_coordinates.csv", excluded_rows)
    report = {
        "canonical_accession": "P35557-1",
        "canonical_length": len(canonical),
        "canonical_sequence_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "isoform_2_length": len(iso2),
        "isoform_2_alternative_feature": alternative,
        "uniprot_subunit_annotation": [
            c for c in uniprot["comments"] if c["commentType"] == "SUBUNIT"
        ],
        "coordinate_metric": "mass-weighted complete sidechain heavy-atom COM; glycine CA fallback",
        "sensitivity_coordinate_metric": "CA",
        "coordinate_completeness_policy": "fail if any mapped sidechain lacks expected heavy atoms",
        "experimental_B_factors_are_not_plddt": True,
        "plddt_policy": "AF>=70 structured; AF<50 candidate IDR; 50<=AF<70 ambiguous; experimental present coordinates trusted; source-confirmed state-specific experimental disorder overrides AF",
        "polymer_eligibility_policy": "both endpoints same contiguous candidate-IDR segment and same chain; no coordinate hybrid",
        "af_coordinate_hybrid_used": False,
        "structures": reports,
        "sources": [
            {
                "local_raw_file": str((RAW / name).relative_to(REPO)),
                "url": url,
                "sha256": hashlib.sha256((RAW / name).read_bytes()).hexdigest(),
                "bytes": (RAW / name).stat().st_size,
            }
            for name, url in sorted(SOURCES.items())
        ],
    }
    (HERE / "geometry_identity_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(reports, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-only", action="store_true")
    args = parser.parse_args()
    download()
    if args.download_only:
        return
    build()


if __name__ == "__main__":
    main()
