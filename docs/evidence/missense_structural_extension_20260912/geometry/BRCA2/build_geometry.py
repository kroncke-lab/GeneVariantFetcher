"""Build domain-local BRCA2 geometry from official biological assemblies.

These are small bound peptides, not full-length BRCA2 or a merged global frame.
All absent canonical residues remain missing; no polymer model is introduced.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
RAW = REPO / "results/missense_structural_extension_20260912/raw/BRCA2"
PDB_IDS = ("7LDG", "8PBC")
SOURCES = {
    "P51587.json": "https://rest.uniprot.org/uniprotkb/P51587.json",
    "best_structures.json": "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/P51587",
    "MSE_chem_comp.json": "https://data.rcsb.org/rest/v1/core/chemcomp/MSE",
}
for _pdb in PDB_IDS:
    SOURCES[f"{_pdb}-assembly1.cif.gz"] = (
        f"https://files.rcsb.org/download/{_pdb.lower()}-assembly1.cif.gz"
    )
    SOURCES[f"{_pdb}_sifts.json"] = (
        f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{_pdb.lower()}"
    )
    for _kind, _suffix in (("entry", ""), ("assembly", "/1"), ("polymer_entity", "/2")):
        SOURCES[f"{_pdb}_{_kind}.json"] = (
            f"https://data.rcsb.org/rest/v1/core/{_kind}/{_pdb}{_suffix}"
        )


def acquire_one(item):
    name, url = item
    target = RAW / name
    if not target.exists():
        response = requests.get(url, timeout=90)
        response.raise_for_status()
        target.write_bytes(response.content)
    return name


def acquire():
    RAW.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=6) as pool:
        for name in pool.map(acquire_one, SOURCES.items()):
            print(name, flush=True)


def load_coordinate_helpers():
    path = HERE.parent / "KCNQ1/build_geometry.py"
    spec = importlib.util.spec_from_file_location("cardiac_geometry_helpers", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_cif(name):
    path = RAW / name
    with gzip.open(path, "rt") as handle:
        table = MMCIF2Dict(handle)
    with gzip.open(path, "rt") as handle:
        structure = MMCIFParser(
            QUIET=True, auth_chains=False, auth_residues=False
        ).get_structure(name, handle)
    if len(structure) != 1:
        raise ValueError("Expected one actual assembly coordinate model")
    return table, structure[0]


def parent_aa(residue):
    return seq1(residue.resname, custom_map={"MSE": "M"})


def coordinate_data(residue, helpers):
    if residue.resname != "MSE":
        return helpers.coordinates(residue)
    ca = np.asarray(residue["CA"].coord, dtype=float)
    return None, ca, ["modified_MSE_not_native_MET_COM"]


def write_rows(name, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = stream.getvalue().encode()
    path = HERE / name
    if name.endswith(".gz"):
        data = gzip.compress(data, mtime=0)
    path.write_bytes(data)


def build():
    prior_report = HERE / "geometry_identity_report.json"
    if prior_report.exists():
        for source in json.loads(prior_report.read_text())["sources"]:
            if (
                hashlib.sha256(
                    (REPO / source["local_raw_file"]).read_bytes()
                ).hexdigest()
                != source["sha256"]
            ):
                raise ValueError("Frozen source bytes changed")
    helpers = load_coordinate_helpers()
    canonical = json.loads((RAW / "P51587.json").read_text())["sequence"]["value"]
    if len(canonical) != 3418:
        raise ValueError("Expected canonical 3418-aa BRCA2")
    (HERE / "P51587-1.fasta").write_text(
        ">P51587-1 canonical BRCA2; UniProt displayed sequence\n" + canonical + "\n"
    )
    reports, summaries, combined, excluded, incomplete = {}, [], [], [], []
    coverage = set()
    for pdb in PDB_IDS:
        table, model = read_cif(f"{pdb}-assembly1.cif.gz")
        entity = json.loads((RAW / f"{pdb}_polymer_entity.json").read_text())
        assembly = json.loads((RAW / f"{pdb}_assembly.json").read_text())
        entry = json.loads((RAW / f"{pdb}_entry.json").read_text())
        sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
        if entity["entity_poly"]["rcsb_mutation_count"] != 0:
            raise ValueError("BRCA2 construct has amino-acid substitutions")
        if {s["ncbi_taxonomy_id"] for s in entity["rcsb_entity_source_organism"]} != {
            9606
        }:
            raise ValueError("Expected human BRCA2")
        alignments = [
            a
            for a in entity["rcsb_polymer_entity_align"]
            if a["reference_database_accession"] == "P51587"
            and a["provenance_source"] == "SIFTS"
        ]
        if len(alignments) != 1:
            raise ValueError("Ambiguous or missing SIFTS canonical mapping")
        mapping = {}
        for region in alignments[0]["aligned_regions"]:
            e, c, length = (
                region[k] for k in ("entity_beg_seq_id", "ref_beg_seq_id", "length")
            )
            fragment = sequence[e - 1 : e - 1 + length]
            if (
                fragment != canonical[c - 1 : c - 1 + length]
                or canonical.count(fragment) != 1
            ):
                raise ValueError(
                    "SIFTS mapping fails independent unique exact sequence alignment"
                )
            for offset in range(length):
                mapping[e + offset] = c + offset
        expected_mapping = (
            {p: p + 2269 for p in range(2, 67)}
            if pdb == "7LDG"
            else {p: p + 3259 for p in range(1, 50)}
        )
        if mapping != expected_mapping:
            raise ValueError("Dated mapped peptide interval changed")
        pdbe = json.loads((RAW / f"{pdb}_sifts.json").read_text())[pdb.lower()][
            "UniProt"
        ]["P51587"]["mappings"]
        pdbe_chains = {m["struct_asym_id"] for m in pdbe}
        for m in pdbe:
            if m["identity"] != 1.0 or (
                m["start"]["residue_number"],
                m["end"]["residue_number"],
                m["unp_start"],
                m["unp_end"],
            ) != (
                min(mapping),
                max(mapping),
                min(mapping.values()),
                max(mapping.values()),
            ):
                raise ValueError("PDBe and RCSB SIFTS mappings disagree")
        asym_entity = dict(
            zip(table["_struct_asym.id"], table["_struct_asym.entity_id"], strict=True)
        )
        target_chains = [chain for chain, e in asym_entity.items() if e == "2"]
        expected_chains = (
            {"B", "D", "B-2", "D-2"} if pdb == "7LDG" else set("LMNOPQRSTU")
        )
        if (
            set(target_chains) != expected_chains
            or {c.split("-", 1)[0] for c in target_chains} != pdbe_chains
        ):
            raise ValueError("Assembly expansion or target chain mapping changed")
        chem_parent = json.loads((RAW / "MSE_chem_comp.json").read_text())["chem_comp"][
            "mon_nstd_parent_comp_id"
        ]
        if pdb == "7LDG" and chem_parent != ["MET"]:
            raise ValueError("Selenomethionine parent mapping not verified")
        identity = helpers.identity_mapping(table)
        inverse = {v: k for k, v in mapping.items()}
        all_rows = []
        for chain in target_chains:
            residues = {
                r.id[1]: r for r in model[chain] if (chain, r.id[1]) in identity
            }
            if len(residues) != len(list(model[chain])):
                raise ValueError("Unexpected unlabelled atoms in target peptide chain")
            for label, residue in residues.items():
                if label not in mapping:
                    excluded.append(
                        {
                            "frame_id": pdb,
                            "chain": chain,
                            "label_seq_id": label,
                            "resname": residue.resname,
                            "reason": "outside_canonical_SIFTS_mapping",
                        }
                    )
                    continue
                pos = mapping[label]
                if parent_aa(residue) != canonical[pos - 1]:
                    raise ValueError("Observed WT differs from canonical sequence")
            rows = []
            for pos, aa in enumerate(canonical, 1):
                label = inverse.get(pos)
                residue = residues.get(label)
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
                    "label_seq_id": label or "",
                    "observed": residue is not None,
                    "geometry_state": "missing",
                    "idr_segment": "",
                    "geometry_state_source": "outside_resolved_bound_peptide",
                    "com_x": "",
                    "com_y": "",
                    "com_z": "",
                    "ca_x": "",
                    "ca_y": "",
                    "ca_z": "",
                    "plddt": "",
                    "plddt_source": "none",
                    "experimental_b_factor": "",
                    "coordinate_source": "none",
                    "mapping_method": "SIFTS_RCSB_PDBe_unique_exact_sequence_WT"
                    if label
                    else "canonical_reference_only",
                    "coordinate_note": "no_polymer_no_AF_no_cross_domain_frame",
                    "in_experimental_construct": label is not None,
                    "ca_geometry_state": "missing",
                    "modified_residue": "",
                }
                if residue is not None:
                    com, ca, missing = coordinate_data(residue, helpers)
                    row["structure_aa"] = parent_aa(residue)
                    row["author_pos"], row["author_chain"], row["insertion_code"] = (
                        identity[(chain, label)]
                    )
                    row["experimental_b_factor"] = float(residue["CA"].bfactor)
                    row["coordinate_source"] = pdb
                    row["ca_geometry_state"] = "structured"
                    row["modified_residue"] = (
                        "MSE_parent_MET" if residue.resname == "MSE" else ""
                    )
                    for axis, value in zip("xyz", ca, strict=True):
                        row[f"ca_{axis}"] = round(float(value), 8)
                    if missing:
                        incomplete.append(
                            {
                                "frame_id": pdb,
                                "chain": chain,
                                "canonical_pos": pos,
                                "missing_atoms": " ".join(missing),
                            }
                        )
                        row["geometry_state_source"] = (
                            "modified_MSE_native_MET_COM_unavailable"
                            if residue.resname == "MSE"
                            else "incomplete_experimental_sidechain_COM_unavailable"
                        )
                        row["coordinate_note"] = (
                            "observed_CA_retained;native_MET_COM_withheld"
                            if residue.resname == "MSE"
                            else "incomplete_sidechain_COM_withheld"
                        )
                    else:
                        row["geometry_state"] = "structured"
                        row["geometry_state_source"] = (
                            "observed_bound_peptide_experimental_coordinates"
                        )
                        row["coordinate_note"] = (
                            "glycine_CA_fallback"
                            if aa == "G"
                            else "complete_sidechain_heavy_atom_COM"
                        )
                        for axis, value in zip("xyz", com, strict=True):
                            row[f"com_{axis}"] = round(float(value), 8)
                        coverage.add(pos)
                rows.append(row)
            counts = Counter(r["geometry_state"] for r in rows)
            summaries.append(
                {
                    "frame_id": pdb,
                    "chain": chain,
                    "canonical_rows": len(rows),
                    "observed_CA": sum(r["observed"] for r in rows),
                    "structured": counts["structured"],
                    "idr": 0,
                    "ambiguous": 0,
                    "missing": counts["missing"],
                }
            )
            all_rows.extend(rows)
        if len(
            {(r["frame_id"], r["chain"], r["canonical_pos"]) for r in all_rows}
        ) != 3418 * len(target_chains):
            raise ValueError("Canonical full-length row grain failed")
        if any(r["geometry_state"] == "idr" or r["idr_segment"] for r in all_rows):
            raise ValueError("Domain-local geometry must not invent polymer support")
        write_rows(f"{pdb}_canonical_geometry.csv.gz", all_rows)
        write_rows(
            f"{pdb}_ca_canonical_geometry.csv.gz",
            [dict(row, geometry_state=row["ca_geometry_state"]) for row in all_rows],
        )
        combined.extend(all_rows)
        descriptions = dict(
            zip(table["_entity.id"], table["_entity.pdbx_description"], strict=True)
        )
        entity_summary = []
        for entity_id, description in descriptions.items():
            chains = [c for c, e in asym_entity.items() if e == entity_id]
            entity_summary.append(
                {
                    "entity_id": entity_id,
                    "description": description,
                    "chains": chains,
                    "variant_donor_target": entity_id == "2",
                }
            )
        positions = sorted(
            {
                r["canonical_pos"]
                for r in all_rows
                if r["geometry_state"] == "structured"
            }
        )
        reports[pdb] = {
            "assembly_id": "1",
            "status": "domain_local_exploratory_only_not_full_gene",
            "canonical_construct_range": [min(mapping.values()), max(mapping.values())],
            "target_chain_ids_in_downloaded_assembly": target_chains,
            "target_sequence_chains_before_assembly_expansion": sorted(pdbe_chains),
            "structured_unique_canonical_positions": positions,
            "unique_residue_coverage_fraction": len(positions) / 3418,
            "resolution_angstrom": entry["rcsb_entry_info"]["resolution_combined"],
            "assembly_record": assembly["pdbx_struct_assembly"],
            "partner_context": entity_summary,
            "source_study_pubmed": entry["rcsb_primary_citation"].get(
                "pdbx_database_id_PubMed", ""
            ),
            "observed_WT_mismatches": 0,
            "SIFTS_sequence_alignment": alignments[0],
            "copy_context_note": "Deposited assembly has B/D/B-2/D-2 BRCA2 fragments; paper's 4:2 MEILB2:BRCA2 description is not four complete BRCA2 polypeptides. Preserve chain fragments without inventing B-D connections."
            if pdb == "7LDG"
            else "Ten BRCA2 peptide copies on a finite eleven-RAD51 filament model; end contexts are not assumed identical.",
        }
    write_rows("BRCA2_canonical_geometry.csv.gz", combined)
    write_rows(
        "BRCA2_ca_canonical_geometry.csv.gz",
        [dict(row, geometry_state=row["ca_geometry_state"]) for row in combined],
    )
    write_rows("geometry_summary.csv", summaries)
    (HERE / "coordinate_exclusions.json").write_text(
        json.dumps(
            {
                "unmapped_observed_residues": excluded,
                "incomplete_sidechains": incomplete,
            },
            indent=2,
        )
        + "\n"
    )
    report = {
        "gene": "BRCA2",
        "canonical_accession": "P51587-1",
        "canonical_length": 3418,
        "canonical_sequence_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "status": "two_separate_domain_local_frames;not_full_gene_structural_coverage",
        "combined_unique_structured_positions": len(coverage),
        "combined_canonical_coverage_fraction": len(coverage) / 3418,
        "frames": reports,
        "polymer_policy": "none;all_unresolved_or_outside_construct_positions_remain_missing",
        "frame_policy": "Separate frame_id values are mandatory; no cross-domain Euclidean distances",
        "modified_residue_policy": "MSE maps to canonicalM only after PDB parentMET verification; native methionine COM withheld; observed CAlpha retained for explicit sensitivity",
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
        "coordinate_helper": {
            "file": str((HERE.parent / "KCNQ1/build_geometry.py").relative_to(REPO)),
            "sha256": hashlib.sha256(
                (HERE.parent / "KCNQ1/build_geometry.py").read_bytes()
            ).hexdigest(),
        },
        "checks": {
            "strict_parent_WT_match": True,
            "SIFTS_unique_sequence_alignment": True,
            "no_incomplete_COM_used": True,
            "full_canonical_missing_rows_materialized": True,
            "no_polymers_or_AF_splices": True,
            "separate_frames_retained": True,
        },
    }
    (HERE / "geometry_identity_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(
        json.dumps(
            {"combined_unique_residues": len(coverage), "summaries": summaries},
            indent=2,
        )
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-only", action="store_true")
    args = parser.parse_args()
    acquire()
    if not args.download_only:
        build()


if __name__ == "__main__":
    main()
