"""Build BRCA2 local ordered frames plus canonical polymer neighborhoods.

Reads pinned earlier GVF evidence and locally cached AFv6 fragments. Never merges
fragment coordinate systems, downloads sources, or edits prior evidence.
"""

import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parents[2]
REPO = HERE.parents[4]
OLD = EVIDENCE / "missense_structural_extension_20260912"
FRAGMENTS = REPO.parent / "ProteinProximityAnalysis/output/BRCA2/structure/fragments"
UNIPROT = REPO / "results/missense_structural_extension_20260912/raw/BRCA2/P51587.json"
ARCHIVE = (
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar"
)
FIELDS = [
    "frame_id",
    "chain",
    "canonical_pos",
    "aa_ref",
    "geometry_state",
    "ca_geometry_state",
    "idr_segment",
    "com_x",
    "com_y",
    "com_z",
    "ca_x",
    "ca_y",
    "ca_z",
    "plddt",
    "geometry_source",
    "source_residue",
    "source_note",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, rows):
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    payload = text.getvalue().encode()
    assert b"\r" not in payload
    blob = gzip.compress(payload, mtime=0) if name.endswith(".gz") else payload
    if len(blob) >= 1_150_000:
        raise ValueError(f"Artifact size exceeds bound: {name} {len(blob)}")
    path = HERE / name
    path.write_bytes(blob)
    return {"file": name, "rows": len(rows), "bytes": len(blob), "sha256": sha(path)}


def empty_row(frame, chain, position, aa):
    row = dict.fromkeys(FIELDS, "")
    row.update(
        frame_id=frame,
        chain=chain,
        canonical_pos=position,
        aa_ref=aa,
        geometry_state="missing",
        ca_geometry_state="missing",
    )
    return row


def segment_labels(positions, prefix):
    positions = set(positions)
    segments = []
    for position in sorted(positions):
        if not segments or position != segments[-1][-1] + 1:
            segments.append([position])
        else:
            segments[-1].append(position)
    return {
        position: f"{prefix}_{block[0]}_{block[-1]}"
        for block in segments
        for position in block
    }


def polymer_rows(sequence, positions, policy, confidence):
    labels = segment_labels(positions, "polymer")
    rows = []
    for position, aa in enumerate(sequence, 1):
        row = empty_row("BRCA2_canonical_polymer", "A", position, aa)
        row["plddt"] = confidence[position]["max"]
        row["geometry_source"] = policy
        if position in labels:
            row.update(
                geometry_state="idr",
                ca_geometry_state="idr",
                idr_segment=labels[position],
                source_residue=position,
            )
        rows.append(row)
    return rows


def build():
    HERE.mkdir(parents=True, exist_ok=True)
    uniprot = json.loads(UNIPROT.read_text())
    sequence = uniprot["sequence"]["value"]
    assert len(sequence) == 3418 and uniprot["primaryAccession"] == "P51587"
    helper_path = OLD / "geometry/KCNQ1/build_geometry.py"
    spec = importlib.util.spec_from_file_location(
        "frozen_coordinate_helper", helper_path
    )
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    uniprot_disorder = set()
    uniprot_regions = []
    for feature in uniprot["features"]:
        if feature.get("description") == "Disordered":
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            uniprot_disorder.update(range(start, end + 1))
            uniprot_regions.append(feature)
    experiment_path = OLD / "geometry/BRCA2/BRCA2_canonical_geometry.csv.gz"
    experimental = pd.read_csv(experiment_path, low_memory=False).fillna("")
    assert all(
        sequence[int(r.canonical_pos) - 1] == r.aa_ref
        for r in experimental.itertuples()
    )
    experimental_positions = set(
        experimental.loc[
            experimental.ca_geometry_state.eq("structured"), "canonical_pos"
        ]
    )
    assert len(experimental_positions) == 73
    experimental_rows = []
    for old in experimental.to_dict("records"):
        row = empty_row(
            old["frame_id"], old["chain"], int(old["canonical_pos"]), old["aa_ref"]
        )
        for column in FIELDS:
            if column in old:
                row[column] = old[column]
        row["geometry_source"] = "experimental_biological_assembly"
        row["source_residue"] = old.get("author_pos", "")
        row["source_note"] = old.get("coordinate_note", "")
        experimental_rows.append(row)
    by_position = {position: [] for position in range(1, len(sequence) + 1)}
    fragments = []
    confidence_rows = []
    source_files = []
    for path in sorted(
        FRAGMENTS.glob("AF-P51587-F*-model_v6.pdb.gz"),
        key=lambda p: int(re.search(r"-F(\d+)-", p.name)[1]),
    ):
        number = int(re.search(r"-F(\d+)-", path.name)[1])
        text = gzip.decompress(path.read_bytes()).decode()
        structure = PDBParser(QUIET=True).get_structure(path.name, io.StringIO(text))
        assert len(structure) == 1 and [chain.id for chain in structure[0]] == ["A"]
        residues = list(structure[0]["A"])
        observed = "".join(seq1(residue.resname) for residue in residues)
        matches = [
            start
            for start in range(len(sequence))
            if sequence.startswith(observed, start)
        ]
        assert len(matches) == 1, (path, matches)
        start = matches[0] + 1
        assert start == 1 + 200 * (number - 1)
        assert len(residues) == min(1400, len(sequence) - start + 1)
        frame = f"AF_P51587_F{number}_v6"
        records = []
        for offset, residue in enumerate(residues):
            position = start + offset
            assert residue.id == (" ", offset + 1, " ")
            assert seq1(residue.resname) == sequence[position - 1]
            com, ca, missing = helper.coordinates(residue)
            assert not missing and com is not None
            plddt = float(residue["CA"].bfactor)
            assert 0 <= plddt <= 100
            by_position[position].append((frame, plddt))
            records.append((position, com, ca, plddt, offset + 1))
            confidence_rows.append(
                {
                    "frame_id": frame,
                    "canonical_pos": position,
                    "local_pos": offset + 1,
                    "aa_ref": sequence[position - 1],
                    "plddt": plddt,
                }
            )
        fragments.append((frame, records))
        source_files.append(
            {
                "file": str(path.relative_to(REPO.parent)),
                "sha256": sha(path),
                "archive_url": ARCHIVE,
                "archive_member": path.name,
                "canonical_start": start,
                "canonical_end": start + len(residues) - 1,
                "unique_exact_sequence_mapping": True,
                "observed_WT_mismatches": 0,
            }
        )
    assert len(fragments) == 12 and all(by_position.values())
    confidence = {
        pos: {
            "min": min(x[1] for x in values),
            "max": max(x[1] for x in values),
            "median": float(np.median([x[1] for x in values])),
            "count": len(values),
        }
        for pos, values in by_position.items()
    }
    # Uniform exclusive routing includes donor eligibility, not just targets.
    idr = {
        pos for pos, c in confidence.items() if c["max"] < 50
    } - experimental_positions
    af_ordered = {
        pos for pos, c in confidence.items() if c["max"] >= 70
    } - experimental_positions
    unresolved = set(by_position) - experimental_positions - idr - af_ordered
    af_rows = []
    for frame, records in fragments:
        present = {record[0]: record for record in records}
        for position, aa in enumerate(sequence, 1):
            row = empty_row(frame, "A", position, aa)
            row["geometry_source"] = "AF_fragment_local"
            if position in present:
                _, com, ca, plddt, local = present[position]
                row.update(plddt=plddt, source_residue=local)
                if position in af_ordered and plddt >= 70:
                    row.update(
                        geometry_state="structured", ca_geometry_state="structured"
                    )
                    for axis, center, alpha in zip("xyz", com, ca, strict=True):
                        row[f"com_{axis}"] = float(center)
                        row[f"ca_{axis}"] = float(alpha)
                elif position in experimental_positions:
                    row["source_note"] = "withheld_experimental_precedence"
                elif position in idr:
                    row["source_note"] = "withheld_global_polymer_assignment"
                else:
                    row["source_note"] = "withheld_local_confidence_below70"
            af_rows.append(row)
    primary = (
        experimental_rows
        + af_rows
        + polymer_rows(
            sequence, idr, "all_covering_AF_fragments_pLDDT_below50", confidence
        )
    )
    conservative = experimental_rows + polymer_rows(
        sequence, idr, "all_covering_AF_fragments_pLDDT_below50", confidence
    )
    broad = (
        experimental_rows
        + af_rows
        + polymer_rows(
            sequence,
            idr | unresolved,
            "all_unresolved_as_polymer_sensitivity",
            confidence,
        )
    )
    outputs = []
    summaries = []
    for name, rows in [
        ("primary_geometry.csv.gz", primary),
        ("experimental_idr_geometry.csv.gz", conservative),
        ("all_unresolved_polymer_geometry.csv.gz", broad),
    ]:
        identities = {(r["frame_id"], r["chain"], r["canonical_pos"]) for r in rows}
        assert len(identities) == len(rows)
        structured = {
            r["canonical_pos"] for r in rows if r["geometry_state"] == "structured"
        }
        polymers = {r["canonical_pos"] for r in rows if r["geometry_state"] == "idr"}
        assert not structured & polymers
        for row in rows:
            assert sequence[row["canonical_pos"] - 1] == row["aa_ref"]
            if row["geometry_state"] != "structured":
                assert all(row[f"com_{axis}"] == "" for axis in "xyz")
            if row["geometry_state"] == "idr":
                assert (
                    row["frame_id"] == "BRCA2_canonical_polymer" and row["idr_segment"]
                )
        outputs.append(save(name, rows))
        summaries.append(
            {
                "policy": name.removesuffix("_geometry.csv.gz"),
                "rows": len(rows),
                "structured_positions": len(structured),
                "polymer_positions": len(polymers),
                "usable_COM_positions": len(structured | polymers),
                "unavailable_COM_positions": len(sequence) - len(structured | polymers),
                "frames": len({r["frame_id"] for r in rows}),
            }
        )
    ledger = []
    for position, aa in enumerate(sequence, 1):
        c = confidence[position]
        route = (
            "experimental"
            if position in experimental_positions
            else "AF_ordered"
            if position in af_ordered
            else "polymer"
            if position in idr
            else "ambiguous_unresolved"
        )
        ledger.append(
            {
                "canonical_pos": position,
                "aa_ref": aa,
                "overlapping_AF_fragments": c["count"],
                "min_plddt": c["min"],
                "max_plddt": c["max"],
                "median_plddt": c["median"],
                "fragment_disagreement_at50": c["min"] < 50 <= c["max"],
                "uniprot_disorder": position in uniprot_disorder,
                "uniprot_disorder_conflicts_max50": position in uniprot_disorder
                and c["max"] >= 50,
                "exclusive_primary_route": route,
                "local_COM_available": position != 2322
                and route != "ambiguous_unresolved",
            }
        )
    outputs.append(save("canonical_confidence_map.csv.gz", ledger))
    outputs.append(save("fragment_confidence.csv.gz", confidence_rows))
    outputs.append(save("geometry_summary.csv", summaries))
    report = {
        "canonical_accession": "P51587-1",
        "canonical_length": len(sequence),
        "canonical_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "prior_omission": "Earlier run only used experimental7LDG/8PBC and AF API404; local AF proteome fragments were not consulted. It explicitly withheld polymer outside constructs, excluding full-gene IDRs.",
        "primary": "primary_geometry.csv.gz",
        "summaries": summaries,
        "exclusive_routing": "Experimental observed CAlpha position first, including modified2322; else AF max-pLDDT>=70 and individual fragment pLDDT>=70; else all-fragment max<50 to polymer; remaining ambiguous unavailable. The same exclusive source applies to targets and donors.",
        "experimental_precedence_positions": sorted(experimental_positions),
        "modified_COM_exclusion": {
            "position": 2322,
            "parent": "MET",
            "observed": "MSE",
            "COM": "withheld",
            "CA": "retained",
        },
        "AF_fragment_local_policy": "Never merge coordinates or compute across frame IDs; all individually>=70 contexts average equally after normalizing donors. Frames are model alternatives, not independent observations or actual protein copies. No PAE files available; relative domain placements inside a fragment remain uncertain.",
        "polymer_policy": "One canonical chain frame across full sequence, outside experimental constructs allowed; only same contiguous segment, same chain; b=3.8 Angstrom, exponent=.5; kernel h3 normalizedsigmoid with positive tails.",
        "disorder_policy": "Conservative consensus: maximum across all covering fragment pLDDTs<50. UniProt MobiDB-lite annotation recorded independently; disagreements not overridden. Low confidence indicates candidate disorder, not experimental proof.",
        "broad_policy": "Retain experimental and AF ordered routes; all remaining unavailable positions become explicitly assumed-polymer, except observed modified2322 COM remains withheld. Adjacent assumed/consensus polymer residues merge into contiguous segments.",
        "confidence_counts": {
            "max_below50": sum(c["max"] < 50 for c in confidence.values()),
            "max_50_to70": sum(50 <= c["max"] < 70 for c in confidence.values()),
            "max_atleast70": sum(c["max"] >= 70 for c in confidence.values()),
            "disagree_at50": sum(
                c["min"] < 50 <= c["max"] for c in confidence.values()
            ),
        },
        "uniprot_disorder_regions": uniprot_regions,
        "limitation": "Full canonical mapping with partial experimental and predicted fragment-local geometry plus polymer; not a full native BRCA2 biological unit, not a global tertiary structure.",
        "sources": source_files
        + [
            {
                "file": str(UNIPROT.relative_to(REPO)),
                "url": "https://rest.uniprot.org/uniprotkb/P51587.json",
                "sha256": sha(UNIPROT),
            },
            {
                "file": str(experiment_path.relative_to(REPO)),
                "sha256": sha(experiment_path),
            },
            {
                "file": str(
                    (OLD / "geometry/BRCA2/geometry_identity_report.json").relative_to(
                        REPO
                    )
                ),
                "sha256": sha(OLD / "geometry/BRCA2/geometry_identity_report.json"),
            },
            {"file": str(helper_path.relative_to(REPO)), "sha256": sha(helper_path)},
        ],
        "outputs": outputs,
        "checks": {
            "all_3418_positions_mapped": True,
            "all_overlapping_sequences_exact": True,
            "zero_WT_mismatches": True,
            "no_incomplete_native_COM": True,
            "exclusive_sources": True,
            "no_cross_frame_geometry": True,
            "one_polymer_context_per_position": True,
            "LF_payloads": True,
        },
    }
    (HERE / "geometry_manifest.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    (HERE / "P51587-1.fasta").write_text(
        ">P51587-1 BRCA2 canonical\n" + sequence + "\n"
    )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    build()
