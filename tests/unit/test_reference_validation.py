"""Unit tests for the reference-transcript validation gate (B1)."""

from pipeline.reference_validation import (
    is_high_confidence_reject,
    is_simple_substitution,
    load_reference_protein,
    parse_reference_residue,
    validate_protein_variant,
)

# Synthetic reference protein: positions 1..7 = M A R G S T C
_SEQ = "MARGSTC"


def test_parse_reference_residue_forms():
    assert parse_reference_residue("R190W") == ("R", 190)
    assert parse_reference_residue("p.Arg190Trp") == ("R", 190)
    assert parse_reference_residue("Arg190Trp") == ("R", 190)
    assert parse_reference_residue("R190X") == ("R", 190)  # nonsense
    assert parse_reference_residue("p.Arg190fs") == ("R", 190)  # frameshift start
    # Non-substitution / non-protein notations have nothing to residue-check.
    assert parse_reference_residue("c.570A>G") is None
    assert parse_reference_residue("190") is None
    assert parse_reference_residue("Xyz12Ala") is None  # not an amino acid
    assert parse_reference_residue("") is None


def test_validate_ok_when_reference_matches():
    assert validate_protein_variant("X", "A2V", sequence=_SEQ).status == "ok"
    assert validate_protein_variant("X", "R3W", sequence=_SEQ).status == "ok"
    assert validate_protein_variant("X", "p.Ala2Val", sequence=_SEQ).status == "ok"


def test_validate_rejects_residue_mismatch():
    v = validate_protein_variant("X", "A3V", sequence=_SEQ)  # pos 3 is R, not A
    assert v.is_reject
    assert v.reason == "reference_residue_mismatch"
    assert v.expected == "R" and v.observed == "A"


def test_validate_rejects_out_of_range():
    v = validate_protein_variant("X", "C99R", sequence=_SEQ)
    assert v.is_reject
    assert v.reason == "position_out_of_range"


def test_validate_unknown_when_no_residue_or_no_sequence():
    # cDNA notation: nothing to residue-check.
    assert validate_protein_variant("X", "c.5A>G", sequence=_SEQ).status == "unknown"
    # No sequence available for the gene: never reject (turnkey-safe).
    v = validate_protein_variant("UNCACHEDGENE", "R3W", sequence=None)
    assert v.status == "unknown"
    assert v.reason == "no_reference_sequence"


def test_load_reference_protein_from_cache_dir(tmp_path):
    (tmp_path / "FAKE.fasta").write_text(">FAKE ref\nMARG\nSTC\n", encoding="utf-8")
    assert load_reference_protein("FAKE", cache_dir=tmp_path) == "MARGSTC"
    assert load_reference_protein("MISSING", cache_dir=tmp_path) is None


def test_validate_uses_cache_dir(tmp_path):
    (tmp_path / "FAKE.fasta").write_text(">FAKE\nMARGSTC\n", encoding="utf-8")
    # pos 3 is R; citing A there must reject using the cached sequence.
    v = validate_protein_variant("FAKE", "A3V", cache_dir=tmp_path)
    assert v.is_reject and v.reason == "reference_residue_mismatch"
    assert validate_protein_variant("FAKE", "R3W", cache_dir=tmp_path).status == "ok"


def test_is_simple_substitution():
    assert is_simple_substitution("A561V")
    assert is_simple_substitution("p.Arg190Trp")
    assert is_simple_substitution("R190X")  # nonsense (single-letter stop)
    assert is_simple_substitution("Y616*")  # nonsense (* stop)
    assert is_simple_substitution("p.Arg190Ter")
    # Frameshifts / indels / extensions are NOT simple substitutions.
    assert not is_simple_substitution("p.Gln1122fs*147")
    assert not is_simple_substitution("p.Lys674fs")
    assert not is_simple_substitution("p.Glu90del")
    assert not is_simple_substitution("p.Ala100_Glu101insLys")
    assert not is_simple_substitution("c.570A>G")
    assert not is_simple_substitution("")


def test_is_high_confidence_reject():
    # pos 1 = M, pos 2..61 = A.
    seq = "M" + "A" * 60
    # Simple substitution mismatch beyond the N-terminal window -> droppable.
    v = validate_protein_variant("X", "R40W", sequence=seq)
    assert v.is_reject and is_high_confidence_reject("R40W", v)
    # Same mismatch INSIDE the N-terminal isoform window -> not droppable.
    v_nt = validate_protein_variant("X", "R5W", sequence=seq)
    assert v_nt.is_reject and not is_high_confidence_reject("R5W", v_nt)
    # Frameshift reject -> not droppable even though far from the N-terminus.
    v_fs = validate_protein_variant("X", "R40fs", sequence=seq)
    assert v_fs.is_reject and not is_high_confidence_reject("R40fs", v_fs)
    # A non-reject verdict is never droppable.
    v_ok = validate_protein_variant("X", "A40V", sequence=seq)
    assert not v_ok.is_reject and not is_high_confidence_reject("A40V", v_ok)
