"""Build frozen HNF1A/LDLR canonical geometry from public biological sources.

Run with ProteinProximityAnalysis/.venv/bin/python. --download fetches public
sources into ignored results/; a normal run rebuilds from those source bytes.
This shared builder writes only the adjacent HNF1A and LDLR evidence folders.
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
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
GEOMETRY = HERE.parent
RAW = ROOT / "results/missense_structural_extension_20260912/raw"
CONFIG = {
    "HNF1A": {
        "accession": "P20823",
        "length": 631,
        "transcript": "ENST00000257555.11",
        "selected": ["8PI8", "1IC8"],
        "candidates": ["8PI8", "1IC8", "2GYP"],
        "copies": 2,
    },
    "LDLR": {
        "accession": "P01130",
        "length": 860,
        "transcript": "ENST00000558518.6",
        "selected": ["1N7D"],
        "candidates": ["1N7D", "9BD8", "9BDE", "9COO"],
        "copies": 1,
    },
}
MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
SIDECHAIN_COUNTS = dict(
    zip(
        "ARNDCQEGHILKMFPSTWYV",
        [1, 7, 4, 4, 2, 5, 5, 0, 6, 4, 4, 5, 4, 7, 3, 2, 3, 10, 8, 3],
        strict=True,
    )
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_urls(gene):
    cfg = CONFIG[gene]
    acc = cfg["accession"]
    urls = {
        f"{acc}.json": f"https://rest.uniprot.org/uniprotkb/{acc}.json",
        "alphafold_api.json": f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",
    }
    for pdb in cfg["candidates"]:
        urls[f"{pdb}_entry.json"] = f"https://data.rcsb.org/rest/v1/core/entry/{pdb}"
        urls[f"{pdb}_assembly.json"] = (
            f"https://data.rcsb.org/rest/v1/core/assembly/{pdb}/1"
        )
        entry_path = RAW / gene / f"{pdb}_entry.json"
        if entry_path.exists():
            for entity in json.loads(entry_path.read_text())[
                "rcsb_entry_container_identifiers"
            ]["polymer_entity_ids"]:
                urls[f"{pdb}_polymer_entity_{entity}.json"] = (
                    f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity}"
                )
    for pdb in cfg["selected"]:
        urls[f"{pdb}-assembly1.cif.gz"] = (
            f"https://files.rcsb.org/download/{pdb.lower()}-assembly1.cif.gz"
        )
    api = RAW / gene / "alphafold_api.json"
    if api.exists():
        models = [
            m for m in json.loads(api.read_text()) if m["uniprotAccession"] == acc
        ]
        if len(models) != 1:
            raise ValueError("Expected exactly one full canonical AlphaFold model")
        urls[f"AF_{acc}.cif"] = models[0]["cifUrl"]
        if gene == "LDLR":
            urls[f"AF_{acc}_pae.json"] = models[0]["paeDocUrl"]
    return urls


def download(gene):
    (RAW / gene).mkdir(parents=True, exist_ok=True)

    def fetch(item):
        name, url = item
        path = RAW / gene / name
        if not path.exists():
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            path.write_bytes(response.content)
        return name

    # Second pass discovers entity IDs and model URLs from the first pass.
    for _ in range(2):
        with ThreadPoolExecutor(max_workers=6) as pool:
            list(pool.map(fetch, source_urls(gene).items()))


def read_cif(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        table = MMCIF2Dict(handle)
    with opener(path, "rt") as handle:
        structure = MMCIFParser(
            QUIET=True, auth_chains=False, auth_residues=False
        ).get_structure(path.stem, handle)
    if len(structure) != 1:
        raise ValueError("Expected one coordinate model")
    identities = {}
    fields = [
        "label_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_ins_code",
    ]
    for chain, label, auth, author_chain, insertion in zip(
        *(table[f"_atom_site.{f}"] for f in fields), strict=True
    ):
        if label in (".", "?"):
            continue
        key = (chain, int(label))
        val = (auth, author_chain, "" if insertion in (".", "?") else insertion)
        if key in identities and identities[key] != val:
            raise ValueError("Ambiguous atom author/label numbering")
        identities[key] = val
    return table, structure[0], identities


def coordinate_values(residue):
    aa = seq1(residue.resname)
    if aa not in SIDECHAIN_COUNTS or "CA" not in residue:
        return None, None, "noncanonical_or_no_CA"
    ca = np.asarray(residue["CA"].coord, dtype=float)
    atoms = [a for a in residue if a.name not in BACKBONE and a.element != "H"]
    if len(atoms) != SIDECHAIN_COUNTS[aa]:
        return None, ca, "incomplete_sidechain"
    com = (
        ca
        if aa == "G"
        else np.average(
            np.asarray([a.coord for a in atoms]),
            axis=0,
            weights=[MASS[a.element] for a in atoms],
        )
    )
    return (
        com,
        ca,
        "glycine_CA_fallback" if aa == "G" else "complete_sidechain_mass_weighted_COM",
    )


def ranges(positions):
    positions = sorted(set(positions))
    if not positions:
        return []
    answer, first, last = [], positions[0], positions[0]
    for pos in positions[1:]:
        if pos != last + 1:
            answer.append([first, last])
            first = pos
        last = pos
    return answer + [[first, last]]


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def initialize_row(gene, frame, chain, pos, aa, confidence, annotated_idr):
    row = dict(
        frame_id=frame,
        chain=chain,
        canonical_pos=pos,
        aa_ref=aa,
        structure_aa="",
        deposited_aa="",
        author_pos="",
        author_chain="",
        insertion_code="",
        label_seq_id="",
        deposited_in_construct=False,
        observed=False,
        geometry_state="missing",
        idr_segment="",
        geometry_state_source="no_usable_coordinates",
        com_x="",
        com_y="",
        com_z="",
        ca_x="",
        ca_y="",
        ca_z="",
        plddt=confidence,
        plddt_source=f"AF-{CONFIG[gene]['accession']}-F1-model_v6",
        uniprot_disorder=pos in annotated_idr,
        experimental_b_factor="",
        coordinate_source="none",
        mapping_method="canonical_reference_only",
        coordinate_note="",
    )
    # These regions cannot be called disordered solely because a free monomer
    # model is uncertain: a known oligomerization domain, a cleaved signal
    # peptide, or a membrane helix has a different biological explanation.
    if gene == "HNF1A" and pos <= 32:
        row["geometry_state_source"] = (
            "known_dimerization_domain_not_present_in_DNA_binding_construct"
        )
    elif gene == "LDLR" and pos <= 21:
        row["geometry_state_source"] = (
            "cleaved_signal_peptide_not_mature_receptor_geometry"
        )
    elif gene == "LDLR" and 789 <= pos <= 810:
        row["geometry_state_source"] = "membrane_helix_not_IDR"
    elif pos in annotated_idr or confidence < 50:
        row["geometry_state"] = "idr"
        row["geometry_state_source"] = (
            "UniProt_disorder_annotation"
            if pos in annotated_idr
            else "AF_plddt_lt50_candidate_IDR"
        )
    elif confidence < 70:
        row["geometry_state"] = "ambiguous"
        row["geometry_state_source"] = "AF_plddt50_to70_no_experimental_coordinates"
    return row


def set_observed(row, residue, identity, *, experimental, allow_coordinates):
    deposited = seq1(residue.resname)
    row.update(
        deposited_aa=deposited,
        observed=True,
        author_pos=identity[0],
        author_chain=identity[1],
        insertion_code=identity[2],
        label_seq_id=residue.id[1],
    )
    if experimental and "CA" in residue:
        row["experimental_b_factor"] = float(residue["CA"].bfactor)
    if deposited != row["aa_ref"]:
        row.update(
            geometry_state="ambiguous",
            geometry_state_source="engineered_nonWT_residue_excluded",
            coordinate_note="WT_mismatch_no_geometry",
        )
        return
    row["structure_aa"] = deposited
    com, ca, note = coordinate_values(residue)
    if not allow_coordinates:
        row["coordinate_note"] = (
            "predicted_coordinates_withheld_confidence_or_biology_gate"
        )
        return
    row["coordinate_note"] = note
    if ca is not None:
        for axis, value in zip("xyz", ca, strict=True):
            row[f"ca_{axis}"] = round(float(value), 8)
    if com is None:
        row.update(
            geometry_state="ambiguous",
            geometry_state_source="incomplete_sidechain_COM_unavailable",
        )
        return
    for axis, value in zip("xyz", com, strict=True):
        row[f"com_{axis}"] = round(float(value), 8)
    row.update(
        geometry_state="structured",
        geometry_state_source="observed_experimental_coordinates"
        if experimental
        else "AF_plddt_ge70",
        coordinate_source=row["frame_id"],
    )


def assign_idr_segments(rows):
    for row in rows:
        row["idr_segment"] = ""
    for chain in sorted({r["chain"] for r in rows}):
        chain_rows = [r for r in rows if r["chain"] == chain]
        idr_ranges = ranges(
            r["canonical_pos"] for r in chain_rows if r["geometry_state"] == "idr"
        )
        for first, last in idr_ranges:
            for row in chain_rows:
                if first <= row["canonical_pos"] <= last:
                    row["idr_segment"] = f"candidate_IDR_{first}_{last}"
                    row["coordinate_note"] = (
                        "sequence_polymer_only_no_experimental_AF_coordinate_hybrid"
                    )


def validate_rows(rows, canonical, copies, metric="com"):
    assert len(rows) == copies * len(canonical)
    keys = [(r["frame_id"], r["chain"], r["canonical_pos"]) for r in rows]
    assert len(keys) == len(set(keys))
    for row in rows:
        assert row["aa_ref"] == canonical[row["canonical_pos"] - 1]
        if row["geometry_state"] == "structured":
            assert row["structure_aa"] == row["aa_ref"]
            assert np.isfinite([row[f"{metric}_{a}"] for a in "xyz"]).all()
        if row["geometry_state"] == "idr":
            assert row["idr_segment"]
            assert all(row[f"{m}_{a}"] == "" for m in ("com", "ca") for a in "xyz")
    for chain in {r["chain"] for r in rows}:
        assert {r["canonical_pos"] for r in rows if r["chain"] == chain} == set(
            range(1, len(canonical) + 1)
        )


def mapped_positions(entity, canonical, acc):
    sequence = "".join(entity["entity_poly"]["pdbx_seq_one_letter_code_can"].split())
    refs = [
        a
        for a in entity.get("rcsb_polymer_entity_align", [])
        if a["reference_database_accession"] == acc
        and a["provenance_source"] == "SIFTS"
    ]
    if len(refs) != 1:
        raise ValueError("Missing or ambiguous canonical SIFTS mapping")
    mapping, mismatches = {}, []
    for segment in refs[0]["aligned_regions"]:
        e, c, length = (
            segment[k] for k in ("entity_beg_seq_id", "ref_beg_seq_id", "length")
        )
        query = sequence[e - 1 : e - 1 + length]
        # Independent ungapped sequence placement must have one optimal start.
        # SIFTS segments are ungapped; gaps occur between segments, never hidden
        # by an author-number offset. Expected engineered mismatches are logged.
        scores = [
            sum(x != y for x, y in zip(query, canonical[s : s + length], strict=True))
            for s in range(len(canonical) - length + 1)
        ]
        if scores.count(min(scores)) != 1 or scores.index(min(scores)) != c - 1:
            raise ValueError("SIFTS canonical placement is not independently unique")
        for j in range(length):
            if e + j in mapping or c + j in mapping.values():
                raise ValueError("Duplicate canonical SIFTS mapping")
            mapping[e + j] = c + j
            if sequence[e + j - 1] != canonical[c + j - 1]:
                mismatches.append(
                    dict(
                        label_seq_id=e + j,
                        canonical_pos=c + j,
                        aa_ref=canonical[c + j - 1],
                        deposited_aa=sequence[e + j - 1],
                    )
                )
    return sequence, mapping, mismatches, refs[0]["aligned_regions"]


def experimental(gene, pdb, canonical, confidence, disorder):
    raw = RAW / gene
    cfg = CONFIG[gene]
    assembly = json.loads((raw / f"{pdb}_assembly.json").read_text())
    entry = json.loads((raw / f"{pdb}_entry.json").read_text())
    entities = [
        json.loads(p.read_text())
        for p in sorted(raw.glob(f"{pdb}_polymer_entity_*.json"))
    ]
    targets = [
        e
        for e in entities
        if cfg["accession"]
        in e["rcsb_polymer_entity_container_identifiers"].get("uniprot_ids", [])
    ]
    if len(targets) != 1:
        raise ValueError("Expected one target protein entity")
    entity = targets[0]
    sequence, mapping, mismatches, aligned = mapped_positions(
        entity, canonical, cfg["accession"]
    )
    expected_mutations = [(515, "N", "Q"), (657, "N", "Q")] if pdb == "1N7D" else []
    if [
        (m["canonical_pos"], m["aa_ref"], m["deposited_aa"]) for m in mismatches
    ] != expected_mutations:
        raise ValueError("Unexpected construct mutations")
    if entity["entity_poly"]["rcsb_mutation_count"] != len(expected_mutations):
        raise ValueError("Mutation metadata disagrees with sequence")
    table, model, identity = read_cif(raw / f"{pdb}-assembly1.cif.gz")
    target_id = entity["rcsb_polymer_entity_container_identifiers"]["entity_id"]
    chains = [
        c
        for c, e in zip(
            table["_struct_asym.id"], table["_struct_asym.entity_id"], strict=True
        )
        if e == target_id and c in model
    ]
    if len(chains) != cfg["copies"]:
        raise ValueError("Target copy count in actual assembly differs")
    if (
        assembly["rcsb_assembly_info"]["polymer_entity_instance_count_protein"]
        != cfg["copies"]
    ):
        raise ValueError("Unexpected biological assembly protein copies")
    if (
        gene == "HNF1A"
        and assembly["rcsb_assembly_info"]["polymer_entity_instance_count_DNA"] != 2
    ):
        raise ValueError("HNF1A frame must retain its DNA-bound biological context")
    cif_seq = dict(
        zip(
            table["_entity_poly.entity_id"],
            table["_entity_poly.pdbx_seq_one_letter_code_can"],
            strict=True,
        )
    )[target_id]
    if "".join(cif_seq.split()) != sequence:
        raise ValueError("Assembly and metadata entity sequences disagree")
    inverse = {c: e for e, c in mapping.items()}
    rows, excluded = [], []
    for chain in chains:
        residues = {r.id[1]: r for r in model[chain] if r.id[0] == " "}
        for label in sorted(set(range(1, len(sequence) + 1)) - set(mapping)):
            excluded.append(
                dict(
                    frame_id=pdb,
                    chain=chain,
                    label_seq_id=label,
                    deposited_aa=sequence[label - 1],
                    observed=label in residues,
                    reason="unmapped_construct_sequence_not_canonical",
                )
            )
        for pos, aa in enumerate(canonical, 1):
            row = initialize_row(gene, pdb, chain, pos, aa, confidence[pos], disorder)
            label = inverse.get(pos)
            row.update(
                deposited_in_construct=label is not None, label_seq_id=label or ""
            )
            if label is not None:
                row["mapping_method"] = (
                    "SIFTS_plus_unique_sequence_placement_per_residue_WT_check"
                )
                row["deposited_aa"] = sequence[label - 1]
            if label in residues:
                residue = residues[label]
                if seq1(residue.resname) != sequence[label - 1]:
                    raise ValueError(
                        "Observed residue identity differs from construct sequence"
                    )
                set_observed(
                    row,
                    residue,
                    identity[(chain, label)],
                    experimental=True,
                    allow_coordinates=True,
                )
            rows.append(row)
    assign_idr_segments(rows)
    validate_rows(rows, canonical, cfg["copies"])
    report = dict(
        source_type="partial_experimental_biological_assembly",
        title=entry["struct"]["title"],
        assembly_id="1",
        target_label_chains=chains,
        target_author_chains=entity["rcsb_polymer_entity_container_identifiers"][
            "auth_asym_ids"
        ],
        full_length_functional_unit=False,
        resolution_angstrom=entry["rcsb_entry_info"]["resolution_combined"],
        construct_canonical_ranges=ranges(mapping.values()),
        sifts_aligned_regions=aligned,
        independently_unique_sequence_placement=True,
        engineered_WT_mismatches_excluded=mismatches,
        geometry_state_counts=dict(Counter(r["geometry_state"] for r in rows)),
        per_chain={
            c: dict(
                observed=sum(r["observed"] for r in rows if r["chain"] == c),
                structured=sum(
                    r["geometry_state"] == "structured" for r in rows if r["chain"] == c
                ),
                observed_canonical_ranges=ranges(
                    r["canonical_pos"]
                    for r in rows
                    if r["chain"] == c and r["observed"]
                ),
            )
            for c in chains
        },
        assembly_definition=assembly["pdbx_struct_assembly"],
        assembly_counts=assembly["rcsb_assembly_info"],
        polymer_partners=[
            dict(
                description=e["rcsb_polymer_entity"].get("pdbx_description"),
                type=e["entity_poly"]["type"],
                **e["rcsb_polymer_entity_container_identifiers"],
            )
            for e in entities
            if e is not entity
        ],
    )
    return rows, report, excluded


def ldlr_domain_sensitivity(rows, uniprot):
    labels = {}
    for feature in uniprot["features"]:
        if feature["type"] in ("Domain", "Transmembrane"):
            a, b = (feature["location"][k]["value"] for k in ("start", "end"))
            for pos in range(a, b + 1):
                labels[pos] = f"domain_{a}_{b}"
    # Six class-B repeats are blades of one beta propeller, not independent
    # folding domains. Retain their intrapropeller geometry in this sensitivity.
    for pos in range(397, 659):
        labels[pos] = "beta_propeller_397_658"
    unassigned = [
        r["canonical_pos"]
        for r in rows
        if r["geometry_state"] == "structured" and r["canonical_pos"] not in labels
    ]
    for a, b in ranges(unassigned):
        for pos in range(a, b + 1):
            labels[pos] = f"unassigned_local_block_{a}_{b}"
    out = []
    for row in rows:
        copy = dict(row)
        block = (
            labels[row["canonical_pos"]]
            if row["geometry_state"] == "structured"
            else "IDR"
            if row["geometry_state"] == "idr"
            else "unavailable"
        )
        copy["frame_id"] = f"AF_P01130_{block}"
        out.append(copy)
    return out


def pae_audit(rows, raw):
    matrix = np.asarray(
        json.loads((raw / "AF_P01130_pae.json").read_text())[0][
            "predicted_aligned_error"
        ]
    )
    if matrix.shape != (860, 860) or not np.isfinite(matrix).all():
        raise ValueError("Invalid full canonical LDLR PAE matrix")
    structured = [r for r in rows if r["geometry_state"] == "structured"]
    positions = np.array([r["canonical_pos"] for r in structured])
    sym = np.maximum(matrix, matrix.T)[np.ix_(positions - 1, positions - 1)]
    xyz = np.array([[r[f"ca_{axis}"] for axis in "xyz"] for r in structured])
    distances = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
    mask = (
        np.triu(np.ones_like(distances, dtype=bool), 1)
        & (np.abs(positions[:, None] - positions[None, :]) > 20)
        & (distances <= 20)
    )
    pairs = []
    for i, j in zip(*np.where(mask & (sym > 10)), strict=True):
        pairs.append(
            dict(
                canonical_pos_1=int(positions[i]),
                canonical_pos_2=int(positions[j]),
                ca_distance_angstrom=round(float(distances[i, j]), 6),
                max_directional_pae_angstrom=float(sym[i, j]),
            )
        )
    return dict(
        criterion="structured CA<=20A; sequence separation>20; max(PAE_ij,PAE_ji)",
        pairs=int(mask.sum()),
        pae_gt10_pairs=int((mask & (sym > 10)).sum()),
        pae_gt20_pairs=int((mask & (sym > 20)).sum()),
        quantiles=dict(
            zip(
                ["min", "q25", "median", "q75", "max"],
                np.quantile(sym[mask], [0, 0.25, 0.5, 0.75, 1]).tolist(),
                strict=True,
            )
        ),
        meaning="pLDDT gates local residue confidence only; interdomain distances remain uncertain; domain-separated geometry provided as conservative sensitivity",
    ), pairs


def build(gene):
    cfg, out, raw = CONFIG[gene], GEOMETRY / gene, RAW / gene
    provenance_path = out / "geometry_provenance.json"
    if provenance_path.exists():
        for source in json.loads(provenance_path.read_text())["sources"]:
            if sha(ROOT / source["path"]) != source["sha256"]:
                raise ValueError("Frozen public source hash changed")
    acc = cfg["accession"]
    uniprot = json.loads((raw / f"{acc}.json").read_text())
    canonical = uniprot["sequence"]["value"]
    if (
        uniprot["primaryAccession"] != acc
        or len(canonical) != cfg["length"]
        or uniprot["organism"]["taxonId"] != 9606
        or gene not in {g["geneName"]["value"] for g in uniprot["genes"]}
    ):
        raise ValueError("Canonical human gene assignment failed")
    transcript = [
        x
        for x in uniprot["uniProtKBCrossReferences"]
        if x["database"] == "Ensembl"
        and x["id"] == cfg["transcript"]
        and x.get("isoformId") == acc + "-1"
    ]
    if len(transcript) != 1:
        raise ValueError("Frozen missense transcript does not map to canonical protein")
    previous = (
        ROOT
        / f"docs/evidence/population_inclusive_penetrance_20260912/population/{gene}_canonical.fasta"
    )
    if canonical != "".join(previous.read_text().splitlines()[1:]):
        raise ValueError("Canonical sequence changed since missense annotation freeze")
    (out / f"{acc}-1.fasta").write_text(
        f">{acc}-1 {gene} canonical; {cfg['transcript']}\n{canonical}\n"
    )
    af = next(
        x
        for x in json.loads((raw / "alphafold_api.json").read_text())
        if x["uniprotAccession"] == acc
    )
    if (
        af["latestVersion"] != 6
        or af["sequence"] != canonical
        or (af["sequenceStart"], af["sequenceEnd"]) != (1, len(canonical))
    ):
        raise ValueError("AF v6 complete canonical model verification failed")
    _, model, identity = read_cif(raw / f"AF_{acc}.cif")
    if [c.id for c in model] != ["A"]:
        raise ValueError("Expected one full AF monomer")
    residues = {r.id[1]: r for r in model["A"] if r.id[0] == " "}
    if (
        set(residues) != set(range(1, len(canonical) + 1))
        or "".join(seq1(residues[p].resname) for p in sorted(residues)) != canonical
    ):
        raise ValueError("AF coordinate sequence WT identity failed")
    confidence = {p: float(r["CA"].bfactor) for p, r in residues.items()}
    disorder_features = [
        f
        for f in uniprot["features"]
        if f["type"] == "Region" and f.get("description") == "Disordered"
    ]
    disorder = {
        p
        for f in disorder_features
        for p in range(
            f["location"]["start"]["value"], f["location"]["end"]["value"] + 1
        )
    }
    write_csv(
        out / f"AF_{acc}_residue_confidence.csv",
        [
            dict(
                canonical_pos=p,
                aa_ref=canonical[p - 1],
                plddt=confidence[p],
                uniprot_disorder=p in disorder,
            )
            for p in confidence
        ],
    )
    reports, outputs, excluded = {}, [], []
    if gene == "LDLR":
        rows = []
        for pos, aa in enumerate(canonical, 1):
            row = initialize_row(
                gene, "AF_P01130", "A", pos, aa, confidence[pos], disorder
            )
            row.update(
                deposited_in_construct=True,
                mapping_method="full_canonical_sequence_and_coordinate_WT_exact_match",
            )
            set_observed(
                row,
                residues[pos],
                identity[("A", pos)],
                experimental=False,
                allow_coordinates=confidence[pos] >= 70 and pos > 21,
            )
            rows.append(row)
        assign_idr_segments(rows)
        validate_rows(rows, canonical, 1)
        write_csv(out / "AF_P01130_canonical_geometry.csv", rows)
        outputs.append("AF_P01130_canonical_geometry.csv")
        split = ldlr_domain_sensitivity(rows, uniprot)
        validate_rows(split, canonical, 1)
        write_csv(out / "AF_P01130_domain_separated_geometry.csv", split)
        outputs.append("AF_P01130_domain_separated_geometry.csv")
        pae, pairs = pae_audit(rows, raw)
        write_csv(out / "AF_P01130_uncertain_close_pairs.csv", pairs)
        reports["AF_P01130"] = dict(
            source_type="predicted_full_canonical_monomer",
            biological_policy="PPA P01130 explicit full-length AF monomer override",
            full_length_functional_unit=False,
            model_id=af["modelEntityId"],
            model_version=6,
            sequence_identity=1.0,
            geometry_state_counts=dict(Counter(r["geometry_state"] for r in rows)),
            excluded_signal_peptide=[1, 21],
            pae_audit=pae,
        )
    for pdb in cfg["selected"]:
        rows, report, exclusions = experimental(
            gene, pdb, canonical, confidence, disorder
        )
        name = f"{pdb}_canonical_geometry.csv"
        write_csv(out / name, rows)
        outputs.append(name)
        reports[pdb] = report
        excluded.extend(exclusions)
        if pdb == "1N7D":
            ca_rows = [dict(r, frame_id="1N7D_CA") for r in rows]
            for row in ca_rows:
                if (
                    row["geometry_state_source"]
                    == "incomplete_sidechain_COM_unavailable"
                    and row["ca_x"] != ""
                ):
                    row.update(
                        geometry_state="structured",
                        geometry_state_source="WT_observed_CA_only_incomplete_sidechain",
                        coordinate_source="1N7D",
                    )
            assign_idr_segments(ca_rows)
            validate_rows(ca_rows, canonical, 1, metric="ca")
            write_csv(out / "1N7D_CA_canonical_geometry.csv", ca_rows)
            outputs.append("1N7D_CA_canonical_geometry.csv")
            reports["1N7D_CA"] = dict(
                source_type="CA_only_sensitivity_partial_mutant_ectodomain",
                allowed_metric="ca",
                geometry_state_counts=dict(
                    Counter(r["geometry_state"] for r in ca_rows)
                ),
            )
    if excluded:
        write_csv(out / "excluded_construct_residues.csv", excluded)
    sources = [
        dict(
            path=str((raw / name).relative_to(ROOT)),
            url=url,
            sha256=sha(raw / name),
            bytes=(raw / name).stat().st_size,
        )
        for name, url in sorted(source_urls(gene).items())
    ]
    candidates = []
    for pdb in cfg["candidates"]:
        entry = json.loads((raw / f"{pdb}_entry.json").read_text())
        assembly = json.loads((raw / f"{pdb}_assembly.json").read_text())
        targets = [
            json.loads(p.read_text())
            for p in sorted(raw.glob(f"{pdb}_polymer_entity_*.json"))
            if acc
            in json.loads(p.read_text())[
                "rcsb_polymer_entity_container_identifiers"
            ].get("uniprot_ids", [])
        ]
        candidates.append(
            dict(
                pdb=pdb,
                title=entry["struct"]["title"],
                selected=pdb in cfg["selected"],
                target_entities=[
                    dict(
                        entity=e["rcsb_polymer_entity_container_identifiers"],
                        mutation_count=e["entity_poly"]["rcsb_mutation_count"],
                        polymer_sequence_length=e["entity_poly"][
                            "rcsb_sample_sequence_length"
                        ],
                    )
                    for e in targets
                ],
                assembly=assembly["pdbx_struct_assembly"],
                assembly_total_modeled_residues=assembly["rcsb_assembly_info"][
                    "modeled_polymer_monomer_count"
                ],
                selection_note="selected; see geometry report"
                if pdb in cfg["selected"]
                else "D-Ala20 analog; not WT dimerization geometry"
                if pdb == "2GYP"
                else "ApoB complex candidate inventoried; full target sequence deposition is not full observed coverage; target coordinates not acquired in this bounded run",
            )
        )
    provenance = dict(
        gene=gene,
        date="2026-09-12",
        accession=acc,
        canonical_length=len(canonical),
        canonical_sequence_sha256=hashlib.sha256(canonical.encode()).hexdigest(),
        canonical_transcript=transcript[0],
        exact_match_prior_canonical_sequence=True,
        full_gene_monomer_coordinates_used=gene == "LDLR",
        uniprot_disorder_features=disorder_features,
        geometry_policy=dict(
            primary_coordinate="complete sidechain heavy-atom mass-weighted COM; glycine CA fallback",
            missingness="Experimental absence alone is not IDR. AF<50 candidate or UniProt disorder annotation, with explicit signal/TM/known-dimerization exclusions; observed experimental coordinates take precedence.",
            alpha_fold="pLDDT>=70 for Cartesian coordinates; no AF coordinates added to experimental frame",
            polymer="Only contiguous final-state IDR, same chain; root engine supplies 3.8*sqrt(sequence separation)",
            copies="Only target-gene chains are geometry donors; partners retained in assembly provenance",
            WT="Exclude every observed nonWT construct residue, without relabeling as canonical WT",
        ),
        candidates=candidates,
        frames=reports,
        sources=sources,
    )
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    checks = dict(
        all_passed=True,
        geometry_files=[
            dict(
                file=name,
                sha256=sha(out / name),
                bytes=(out / name).stat().st_size,
                LF_only=b"\r" not in (out / name).read_bytes(),
            )
            for name in outputs
        ],
        source_count=len(sources),
        canonical_WT_checks="every structured row passed",
        exact_canonical_positions_per_chain=True,
        no_duplicate_context_positions=True,
        no_non_gene_donor_chains=True,
        engineered_mismatch_geometry_excluded=True,
        no_IDR_coordinates=True,
        no_cross_frame_hybrid=True,
    )
    (out / "geometry_checks.json").write_text(
        json.dumps(checks, indent=2, sort_keys=True) + "\n"
    )
    print(
        gene,
        json.dumps({k: v["geometry_state_counts"] for k, v in reports.items()}),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", choices=list(CONFIG))
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    for target_gene in [args.gene] if args.gene else CONFIG:
        if args.download:
            download(target_gene)
        build(target_gene)
