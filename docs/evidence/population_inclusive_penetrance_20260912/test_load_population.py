"""Source-loading checks for the full-locus/canonical-footprint union."""

import hashlib
import json

import pandas as pd
import pytest

import rebuild_union


def region_row(identity="7-110-A-T", **changes):
    chrom, pos, ref, alt = identity.split("-")
    row = {
        "gene": "GCK",
        "variant_id": identity,
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "gnomad_release": "4.1.1",
        "joint_ac": 8,
        "joint_an": 1000,
        "joint_hom": 1,
        "gnomad_carriers": 7,
        "qc_pass": True,
        "population_eligible": True,
        "exome_ac": 5,
        "exome_an": 800,
        "exome_hom": 1,
        "exome_filters": "",
        "genome_ac": 3,
        "genome_an": 200,
        "genome_hom": 0,
        "genome_filters": "",
        "joint_filters": "",
        "reconstructed_joint_filter": "PASS",
    }
    row.update(changes)
    return row


def annotation_row(identity="7-110-A-T", **changes):
    row = region_row(identity)
    row.update(
        canonical_transcript_id="ENST00000403799.8",
        transcript_id="ENST00000403799.8",
        hgvsc="c.29C>T",
        hgvsp="p.Ala10Val",
        protein_key="A10V",
        aa_ref="A",
        aa_pos=10,
        aa_alt="V",
        canonical_wt_status="match",
        canonical_annotation_status="canonical_transcript_endpoint",
        coding_or_splice=True,
        vclass="missense",
    )
    row.update(changes)
    return row


@pytest.fixture
def source_files(tmp_path, monkeypatch):
    monkeypatch.setattr(rebuild_union, "HERE", tmp_path)
    monkeypatch.setattr(rebuild_union, "GENES", ["GCK"])
    parent = tmp_path / "population"
    full = parent / "full_locus"
    full.mkdir(parents=True)
    footprint = pd.DataFrame(
        [
            annotation_row(),
            annotation_row(
                "7-95-C-G",
                hgvsc="c.29-5C>G",
                hgvsp="",
                protein_key="",
                aa_ref="",
                aa_pos=None,
                aa_alt="",
                vclass="splice",
            ),
        ]
    )
    footprint.to_csv(parent / "population_variants.csv.gz", index=False)

    def write_region(rows):
        path = full / "GCK_population_variants.csv.gz"
        pd.DataFrame(rows).to_csv(path, index=False)
        manifest = {
            "gene": "GCK",
            "interval": {"chrom": "7", "start": 100, "stop": 200},
            "source_rows": len(rows),
            "complete_interval_coverage_verified": True,
            "output_files": [
                {
                    "file": path.name,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "rows": len(rows),
                    "bytes": path.stat().st_size,
                }
            ],
        }
        (full / "GCK_provenance.json").write_text(json.dumps(manifest))
        return path

    write_region([region_row(), region_row("7-150-C-T")])
    return write_region


def test_exact_annotation_merge_preserves_region_extras_and_boundary_footprint(
    source_files,
):
    frame, sources = rebuild_union.load_population()
    assert len(frame) == 3 and frame.variant_id.is_unique
    frame = frame.set_index("variant_id")
    assert frame.loc["7-110-A-T", "protein_key"] == "A10V"
    assert bool(frame.loc["7-110-A-T", "in_gene_span"])
    assert bool(frame.loc["7-110-A-T", "in_coding_footprint"])
    extra = frame.loc["7-150-C-T"]
    assert extra.vclass == "region_only_unannotated"
    assert extra.canonical_wt_status == "unavailable"
    assert not bool(extra.coding_or_splice)
    assert pd.isna(extra.protein_key) and pd.isna(extra.hgvsp)
    assert pd.isna(extra.aa_ref) and pd.isna(extra.aa_pos)
    assert not bool(extra.in_coding_footprint) and bool(extra.in_gene_span)
    boundary = frame.loc["7-95-C-G"]
    assert bool(boundary.in_coding_footprint) and not bool(boundary.in_gene_span)
    assert len(sources) == 3


@pytest.mark.parametrize(
    "changes",
    [
        {"joint_ac": 9},
        {"qc_pass": False},
        {"exome_ac": 6},
        {"genome_filters": "AS_VQSR"},
    ],
)
def test_overlapping_count_and_assay_qc_differences_are_rejected(source_files, changes):
    source_files([region_row(**changes), region_row("7-150-C-T")])
    with pytest.raises(ValueError, match="Population overlap changed"):
        rebuild_union.load_population()


def test_missing_interior_footprint_allele_is_not_relabelled_as_boundary(source_files):
    source_files([region_row("7-150-C-T")])
    with pytest.raises(
        ValueError,
        match="(?i)missing.*(span|interval|interior)|(?i:interior).*(missing|absent)",
    ):
        rebuild_union.load_population()


def test_file_digest_is_checked_before_consuming_csv(source_files):
    path = source_files([region_row(), region_row("7-150-C-T")])
    path.write_bytes(b"corrupted CSV content")
    with pytest.raises(ValueError, match="Population source hash mismatch"):
        rebuild_union.load_population()


def test_duplicate_alleles_across_shards_are_rejected(source_files):
    source_files([region_row(), region_row(), region_row("7-150-C-T")])
    with pytest.raises(ValueError, match="Duplicate full-locus allele"):
        rebuild_union.load_population()


def test_absent_assay_values_compare_equal_without_invented_zeroes(source_files):
    parent = rebuild_union.HERE / "population"
    footprint = pd.read_csv(parent / "population_variants.csv.gz")
    columns = ["genome_ac", "genome_an", "genome_hom"]
    footprint[columns] = float("nan")
    footprint["genome_filters"] = "ABSENT"
    footprint["joint_ac"] = 5
    footprint["gnomad_carriers"] = 4
    footprint.to_csv(parent / "population_variants.csv.gz", index=False)
    source_files(
        [
            region_row(
                genome_ac=None,
                genome_an=None,
                genome_hom=None,
                genome_filters="ABSENT",
                joint_ac=5,
                gnomad_carriers=4,
            ),
            region_row("7-150-C-T"),
        ]
    )
    frame, _ = rebuild_union.load_population()
    row = frame.loc[frame.variant_id.eq("7-110-A-T")].iloc[0]
    assert pd.isna(row.genome_ac) and row.genome_filters == "ABSENT"
