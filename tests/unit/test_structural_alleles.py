"""Generic structural-allele matching against RefSeq exon maps."""

from pipeline.reference_validation import load_reference_protein
from utils.structural_alleles import (
    expand_structural_keys,
    parse_delta_peptide,
    parse_exon_event,
    structural_keys_match,
    _unique_peptide_span,
)


def test_exon_event_parser_is_notation_agnostic():
    assert parse_exon_event("EXON 3 DELETION") == ("del", 3, 3)
    assert parse_exon_event("deletion of exons 15-28") == ("del", 15, 28)
    assert parse_exon_event("exon 9–10 dup") == ("dup", 9, 10)
    assert parse_exon_event("del:exon3") == ("del", 3, 3)


def test_in_frame_exon_deletion_equals_protein_range():
    assert structural_keys_match("EXON 3 DELETION", "p.Asn57_Gly91del", "RYR2")
    keys = expand_structural_keys("EXON 3 DELETION", "RYR2")
    assert "N57_G91del" in keys
    assert expand_structural_keys("p.Asn57_Gly91del", "RYR2") & keys


def test_exon_deletion_does_not_cross_genes():
    assert not structural_keys_match("EXON 3 DELETION", "p.Asn57_Gly91del", "KCNH2")
    assert not structural_keys_match("EXON 3 DELETION", "p.Asn57_Gly91del", "SCN5A")


def test_delta_peptide_uses_unique_reference_hit():
    assert parse_delta_peptide("ΔKPQ") == "KPQ"
    assert parse_delta_peptide("delta KPQ") == "KPQ"
    assert structural_keys_match("ΔKPQ", "p.K1505_Q1507del", "SCN5A")
    assert not structural_keys_match("ΔKPQ", "p.K1505_Q1507del", "KCNH2")


def test_ambiguous_peptide_fails_closed():
    seq = load_reference_protein("SCN5A")
    assert seq is not None
    # Alanine is not unique in SCN5A, so a 1-residue scan would be unsafe;
    # the matcher only accepts a peptide with exactly one occurrence.
    assert _unique_peptide_span(seq, "A") is None
    assert _unique_peptide_span(seq, "KPQ") == (1505, 1507)


def test_ins_at_position_uses_reference_residue():
    seq = load_reference_protein("KCNH2")
    assert seq is not None
    keys = expand_structural_keys("INS GCT 628", "KCNH2")
    assert f"{seq[627]}628ins" in keys
