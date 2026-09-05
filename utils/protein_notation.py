"""Strict validation for the protein-notation subset accepted by GVF.

The extraction and migration paths must make the same keep/drop decision.  Keep
the grammar here so a migration replay cannot silently reject notation that the
extractor accepted (or vice versa).  This intentionally validates a practical
HGVS subset; normalization remains the responsibility of
``utils.variant_normalizer``.
"""

from __future__ import annotations

import re


AA3_PATTERN = (
    r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|"
    r"Ser|Thr|Trp|Tyr|Val|Ter|Stop|Xaa)"
)
AA1_PATTERN = r"[ACDEFGHIKLMNPQRSTVWY]"
AMINO_ACID_PATTERN = rf"(?:{AA3_PATTERN}|{AA1_PATTERN})"
AMINO_ACID_SEQUENCE_PATTERN = rf"(?:{AMINO_ACID_PATTERN})+"

# Some legacy clinical tables use a residue-position label for a splice event
# instead of HGVS cDNA notation.  Keep this as a separate, exact alternative:
# ``sp``/``splice`` is only valid immediately after one residue position, not
# after a range and not as the prefix of an arbitrary trailing annotation.
LEGACY_SPLICE_LABEL_PATTERN = rf"(?:p\.)?{AMINO_ACID_PATTERN}\d{{1,4}}(?:splice|sp)"

# Supported changes include substitutions/termination, frameshifts with an
# optional changed residue and stop distance, and simple/range indels.  Examples
# that exercise the less obvious branches:
#   p.Gly262Alafs*98       changed residue before ``fs``
#   p.Gly24fsTer58         three-letter termination marker
#   p.Ala100_Glu101insLys  three-letter insertion sequence
#   p.Arg176delinsLysGly   multi-residue replacement sequence
PROTEIN_NOTATION_RE = re.compile(
    rf"(?:{LEGACY_SPLICE_LABEL_PATTERN}"
    rf"|(?:p\.)?{AMINO_ACID_PATTERN}"
    r"\d{1,4}"
    rf"(?:[_-]{AMINO_ACID_PATTERN}\d{{1,4}})?"
    rf"(?:{AMINO_ACID_PATTERN}|[X*?=]"
    rf"|(?:{AMINO_ACID_PATTERN})?fs(?:Ter|Stop|X|\*)?\d*"
    rf"|delins{AMINO_ACID_SEQUENCE_PATTERN}"
    rf"|del(?:{AMINO_ACID_SEQUENCE_PATTERN})?"
    rf"|dup(?:{AMINO_ACID_SEQUENCE_PATTERN})?"
    rf"|ins(?:{AMINO_ACID_SEQUENCE_PATTERN})?))",
    re.IGNORECASE,
)

_AA3_TO_1 = dict(
    zip(
        "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".split(),
        "ARNDCQEGHILKMFPSTWYV",
    )
)
_RESIDUE = rf"(?:{'|'.join(_AA3_TO_1)}|{AA1_PATTERN})"
_PARTIAL_INSERTION = re.compile(
    rf"^p\.(?P<start>[1-9]\d{{0,3}})[_-](?P<right>{_RESIDUE})"
    rf"(?P<end>[1-9]\d{{0,3}})ins(?P<insert>(?:{_RESIDUE}){{1,10}})$",
    re.IGNORECASE,
)


def normalize_partial_protein_insertion(value: str) -> str | None:
    """Preserve a source insertion with an unspecified left reference residue.

    ``p.1570-Phe1571insIle`` explicitly identifies adjacent flanking positions
    and an inserted residue. Convert only amino-acid spelling and the separator;
    never infer the missing amino acid or a transcript/cDNA allele. Nonadjacent
    positions, unknown residues, stop codons, and trailing prose are refused.
    """
    match = _PARTIAL_INSERTION.fullmatch(value)
    if not match or int(match["end"]) != int(match["start"]) + 1:
        return None

    def one_letter(token: str) -> str:
        return _AA3_TO_1[token.capitalize()] if len(token) == 3 else token.upper()

    inserted = "".join(
        one_letter(token) for token in re.findall(_RESIDUE, match["insert"], re.I)
    )
    return f"p.{match['start']}_{one_letter(match['right'])}{match['end']}ins{inserted}"


def is_supported_protein_notation(value: str) -> bool:
    """Accept the strict grammar or the explicitly bounded source shorthand."""
    return bool(
        PROTEIN_NOTATION_RE.fullmatch(value)
        or normalize_partial_protein_insertion(value)
    )
