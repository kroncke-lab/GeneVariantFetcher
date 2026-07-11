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

# Supported changes include substitutions/termination, frameshifts with an
# optional changed residue and stop distance, and simple/range indels.  Examples
# that exercise the less obvious branches:
#   p.Gly262Alafs*98       changed residue before ``fs``
#   p.Gly24fsTer58         three-letter termination marker
#   p.Ala100_Glu101insLys  three-letter insertion sequence
#   p.Arg176delinsLysGly   multi-residue replacement sequence
PROTEIN_NOTATION_RE = re.compile(
    rf"(?:p\.)?{AMINO_ACID_PATTERN}"
    r"\d{1,4}"
    rf"(?:[_-]{AMINO_ACID_PATTERN}\d{{1,4}})?"
    rf"(?:{AMINO_ACID_PATTERN}|[X*?=]"
    rf"|(?:{AMINO_ACID_PATTERN})?fs(?:Ter|Stop|X|\*)?\d*"
    rf"|delins{AMINO_ACID_SEQUENCE_PATTERN}"
    rf"|del(?:{AMINO_ACID_SEQUENCE_PATTERN})?"
    rf"|dup(?:{AMINO_ACID_SEQUENCE_PATTERN})?"
    rf"|ins(?:{AMINO_ACID_SEQUENCE_PATTERN})?)",
    re.IGNORECASE,
)


def is_valid_protein_notation(value: str) -> bool:
    """Return whether *value* is entirely covered by the GVF protein grammar."""

    return bool(PROTEIN_NOTATION_RE.fullmatch(value))
