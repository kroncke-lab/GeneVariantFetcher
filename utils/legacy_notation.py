"""Honest identity handling for strict source-only legacy indel notation.

Older papers can report a BIC-style nucleotide label such as ``4321delAC``
without a transcript or HGVS prefix.  That string is useful curator evidence,
but prepending ``c.`` invents a coordinate system the source did not provide.
This module recognizes only the narrow, high-signal grammar we can preserve as
an identity without claiming it is HGVS.
"""

from __future__ import annotations

import re
from typing import Any, MutableMapping


LEGACY_NOTATION_GENES = frozenset({"BRCA1", "BRCA2"})


_LEGACY_INDEL_RE = re.compile(
    r"^(?P<position>\d{3,5}(?:(?:_|-)\d{3,5}(?:[+-]\d+)?)?)"
    r"(?P<operation>(?i:del|dup|ins))"
    # BIC writes either the affected bases (``4321delAC``) or their count
    # (``4184del4``, ``1294del40``). Both are the same label class; accepting
    # only the base form silently discarded the count form, and the extraction
    # prompt now routes bare labels here rather than into ``cdna_notation``.
    r"(?P<payload>[ACGT]{1,20}|\d{1,3})"
    r"(?P<insert>(?i:ins)(?:[ACGT]{1,20}|\d{1,3}))?$"
)


def normalize_legacy_notation(value: Any) -> str | None:
    """Return the normalized strict legacy indel label, or ``None``.

    Whitespace introduced by PDF/table conversion is ignored.  Bases must be
    uppercase in the source-shaped value; this keeps ordinary prose such as
    ``120delta`` outside the grammar.  The operation is normalized lowercase.
    """

    compact = re.sub(r"\s+", "", str(value or "").strip())
    match = _LEGACY_INDEL_RE.fullmatch(compact)
    if not match:
        return None
    return (
        f"{match.group('position')}"
        f"{match.group('operation').lower()}"
        f"{match.group('payload')}"
        f"{(match.group('insert') or '').lower()}"
    )


# Protein-level legacy spellings older channelopathy papers use for a variant
# they never write in HGVS. ``L860fsX89`` / ``L860fsx89`` is a frameshift at
# Leu860 with the stop 89 codons downstream; ``1795insD`` / ``1570insI`` is an
# amino-acid insertion after residue 1795 / 1570. The scorer already bridges
# both to gold (``matches("L860fsX", "L860fsx89")`` and
# ``matches("P.Y1795_E1796INSD", "1795insD")`` are True), but the extraction
# filter dropped them as nameless rows because the model left every identity
# field null and kept only ``source_notation`` (tranche 01, SCN5A 16764707:
# three gold rows lost this way).
_PROTEIN_FRAMESHIFT_RE = re.compile(
    r"^(?P<ref>[ACDEFGHIKLMNPQRSTVWY])(?P<pos>\d{1,5})fs(?:[xX*]\d{0,4})?$"
)
_PROTEIN_LEGACY_INDEL_RE = re.compile(
    r"^(?P<pos>\d{1,5})(?P<op>(?i:ins|del))(?P<aa>[A-Z]{1,6})$"
)
# N is the universal nucleotide wildcard; a payload spelled only in these
# letters is a nucleotide indel (or ambiguous), never a protein-level one.
_NUCLEOTIDES = frozenset("ACGTN")


def promote_source_only_protein_identity(
    variant: MutableMapping[str, Any],
    gene_symbol: Any = None,
) -> str | None:
    """Give a source-only legacy *protein* spelling an identity field.

    Only fires when ``cdna_notation``, ``protein_notation``, and
    ``legacy_notation`` are all empty and ``source_notation`` is one of two
    narrow shapes:

    * a single-letter frameshift ``<Ref><pos>fs[X<n>]`` → ``protein_notation``
      is set to the canonical ``<Ref><pos>fsX`` the rest of the pipeline uses;
    * a residue-position amino-acid indel ``<pos>(ins|del)<AA>`` whose payload
      is not spelled only in A/C/G/T → ``legacy_notation`` keeps the verbatim
      token (a payload of only A/C/G/T letters is ambiguous with a nucleotide
      indel and is left to the cDNA rules).

    Returns the identity written, or ``None`` when nothing changed. Never
    invents coordinates or a transcript.
    """

    for key in ("cdna_notation", "protein_notation", "legacy_notation"):
        if str(variant.get(key) or "").strip():
            return None
    token = re.sub(r"\s+", "", str(variant.get("source_notation") or "").strip())
    if not token:
        return None
    fs = _PROTEIN_FRAMESHIFT_RE.fullmatch(token)
    if fs:
        identity = f"{fs.group('ref')}{fs.group('pos')}fsX"
        variant["protein_notation"] = identity
        return identity
    indel = _PROTEIN_LEGACY_INDEL_RE.fullmatch(token)
    if indel:
        payload = indel.group("aa")
        if set(payload.upper()) <= _NUCLEOTIDES:
            return None
        identity = f"{indel.group('pos')}{indel.group('op').lower()}{payload}"
        variant["legacy_notation"] = identity
        return identity
    return None


def gene_supports_legacy_notation(gene_symbol: Any) -> bool:
    """Return whether bare strict indels are a known legacy identity for a gene."""

    return str(gene_symbol or "").strip().upper() in LEGACY_NOTATION_GENES


def _prefixed_count_form(value: Any) -> str | None:
    """Return the legacy body of an invalid c.-prefixed deleted-count form."""

    compact = re.sub(r"\s+", "", str(value or "").strip())
    if not compact.lower().startswith("c."):
        return None
    body = normalize_legacy_notation(compact[2:])
    if body and re.search(
        r"(?:del|dup|ins)\d{1,3}(?:ins\d{1,3})?$", body, re.IGNORECASE
    ):
        return body
    return None


def preserve_source_only_legacy_identity(
    variant: MutableMapping[str, Any],
    gene_symbol: Any = None,
) -> str | None:
    """Move a strict source-only legacy identity out of ``cdna_notation``.

    The helper mutates an extraction-shaped mapping.  A legacy identity is
    If an upstream layer fabricated ``c.`` while preserving the verbatim bare
    string in ``source_notation``, the fabricated cDNA value is cleared even
    when the row also has a real protein/genomic identity. Bare strict indels
    are treated as legacy only for genes whose historical nomenclature uses
    them (currently BRCA1/BRCA2); for other genes they remain cDNA.
    """

    effective_gene = gene_symbol or variant.get("gene_symbol")
    legacy_gene = gene_supports_legacy_notation(effective_gene)
    protein = str(variant.get("protein_notation") or "").strip()
    genomic = str(variant.get("genomic_position") or "").strip()
    structural = str(variant.get("structural_description") or "").strip()
    cdna = re.sub(r"\s+", "", str(variant.get("cdna_notation") or "").strip())
    source_legacy = normalize_legacy_notation(variant.get("source_notation"))
    explicit_legacy = normalize_legacy_notation(variant.get("legacy_notation"))
    bare_cdna_legacy = normalize_legacy_notation(cdna)
    prefixed_count = _prefixed_count_form(variant.get("source_notation")) or (
        _prefixed_count_form(cdna)
    )

    if not legacy_gene and effective_gene:
        # Prefixless simple indels are common omitted-prefix HGVS in non-BRCA
        # literature. Keeping them in a BIC-only field would hide real BMPR2
        # variants from publication. Promote only the same strict token; never
        # infer a coordinate from free text.
        source_identity = explicit_legacy or source_legacy or bare_cdna_legacy
        if source_identity and (
            not cdna
            or bare_cdna_legacy
            or cdna.casefold() == f"c.{source_identity}".casefold()
        ):
            variant["cdna_notation"] = f"c.{source_identity}"
        variant["legacy_notation"] = None
        return None

    if bare_cdna_legacy:
        variant["cdna_notation"] = None
        variant["legacy_notation"] = bare_cdna_legacy
        return bare_cdna_legacy
    if prefixed_count and cdna.casefold() == f"c.{prefixed_count}".casefold():
        variant["cdna_notation"] = None
        variant["legacy_notation"] = prefixed_count
        return prefixed_count
    if source_legacy and cdna.casefold() == f"c.{source_legacy}".casefold():
        variant["cdna_notation"] = None
        variant["legacy_notation"] = source_legacy
        return source_legacy
    if explicit_legacy and cdna.casefold() == f"c.{explicit_legacy}".casefold():
        # Conflicting identity classes must resolve conservatively. An
        # explicit source-only field is enough to reject the prefixed form as
        # asserted HGVS even when a real protein identity is also present.
        variant["cdna_notation"] = None
        variant["legacy_notation"] = explicit_legacy
        return explicit_legacy

    source_backed = explicit_legacy and explicit_legacy == source_legacy
    if (protein or cdna or genomic or structural) and not source_backed:
        # ``legacy_notation`` is an identity class, not another alias field.
        # Verbatim aliases remain available in variant_papers.source_notation.
        variant["legacy_notation"] = None
        return None
    if source_backed and (protein or cdna or genomic or structural):
        # The paper itself printed this bare label, so it is a real identity
        # rather than a duplicated alias. Keep it beside the other class: the
        # co-present value may be junk that a later guard nulls (a transcript
        # accession in the cDNA slot is a common model slip), and discarding
        # the source-proven label first would delete the whole variant.
        variant["legacy_notation"] = explicit_legacy
        return explicit_legacy

    legacy = explicit_legacy or source_legacy
    variant["legacy_notation"] = legacy
    return legacy
