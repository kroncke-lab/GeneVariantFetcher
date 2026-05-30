"""Reference-transcript validation gate (B1).

A gold-free precision gate: reject a protein variant whose cited reference amino
acid does not match the residue at that position in the gene's reference protein,
or whose position is out of range. It is a **no-op for any gene without a cached
reference sequence**, so it is turnkey-safe on uncached / unknown genes — it can
only ever reject when it has ground truth to reject against.

Reference protein sequences live under ``data/reference_sequences/<GENE>.fasta``
and are populated on demand by ``scripts/fetch_reference_sequences.py`` (NCBI
Entrez; network). The expected length for the curated genes is cross-checked
against ``utils.variant_normalizer.PROTEIN_LENGTHS`` when a sequence is loaded, so
a wrong accession is caught rather than silently trusted.

This module deliberately has no third-party dependency (it only reads a FASTA
text file); only the fetch script needs ``Bio.Entrez``.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_REF_DIR = Path(__file__).resolve().parents[1] / "data" / "reference_sequences"

# Standard amino-acid three-letter -> one-letter map (plus Sec/Pyl).
_AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}
_AA1 = set(_AA3_TO_1.values())

# Leading "<ref-AA><position>" of a protein change, with optional p. prefix.
# Captures the reference residue + position regardless of what follows
# (substitution, nonsense, frameshift, indel start), so we can residue-check any
# notation that names a concrete reference residue at a concrete position.
_LEADING_REF_RE = re.compile(r"^\s*(?:p\.?\s*)?([A-Za-z]{3}|[A-Za-z])(\d+)")

# A clean single-residue substitution or nonsense: <ref><pos><alt> where alt is a
# single residue, a stop (X / * / Ter), or a 3-letter residue. Frameshifts and
# indels are excluded by the keyword guard below — for those the leading residue
# is the first *affected* codon, whose identity is notation-convention dependent
# and a poor target for a hard reference-residue check.
_SIMPLE_SUB_RE = re.compile(
    r"^\s*(?:p\.?\s*)?(?:[A-Za-z]{3}|[A-Za-z])\d+(?:[A-Za-z]{3}|[A-Za-z]|\*)\s*$"
)
_COMPLEX_NOTATION_TOKENS = ("fs", "del", "ins", "dup", "ext")

# Residues at/below this position are in the N-terminal window where alternative
# transcript isoforms most commonly diverge from the canonical RefSeq protein
# (e.g. KCNH2 isoform 1a vs 1b). A reference-residue mismatch there is far more
# likely an isoform/numbering artifact than a true false positive, so such a
# reject is reported but never acted on by a destructive (drop) policy.
NTERM_ISOFORM_GUARD = 30

# in-process cache: gene (upper) -> sequence string, or None if known-absent.
_SEQ_CACHE: dict[str, Optional[str]] = {}


@dataclass(frozen=True)
class ReferenceVerdict:
    """Outcome of validating one variant against a reference protein.

    ``status`` is one of:
      * ``"ok"``      — reference residue matches at the cited position;
      * ``"reject"``  — position out of range or reference residue mismatch;
      * ``"unknown"`` — cannot check (no cached sequence, or notation has no
                        concrete reference residue+position, e.g. cDNA-only).
    """

    status: str
    reason: str
    expected: Optional[str] = None
    observed: Optional[str] = None

    @property
    def is_reject(self) -> bool:
        return self.status == "reject"


def _ref_one_letter(token: str) -> Optional[str]:
    """Map a 1- or 3-letter reference-AA token to a single letter, or None."""
    token = token.strip()
    if len(token) == 1:
        up = token.upper()
        return up if up in _AA1 else None
    return _AA3_TO_1.get(token.upper())


def parse_reference_residue(protein_notation: str) -> Optional[tuple[str, int]]:
    """Return ``(ref_one_letter, position)`` for a protein change, else None.

    None means the notation has no concrete leading reference residue+position to
    check (cDNA notation, bare positions, malformed strings) — the gate then
    returns ``unknown`` rather than rejecting.
    """
    if not protein_notation:
        return None
    m = _LEADING_REF_RE.match(str(protein_notation))
    if not m:
        return None
    ref = _ref_one_letter(m.group(1))
    if ref is None:
        return None
    try:
        pos = int(m.group(2))
    except ValueError:
        return None
    if pos < 1:
        return None
    return ref, pos


def _read_fasta_sequence(path: Path) -> Optional[str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    seq = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    return seq or None


def load_reference_protein(
    gene: str, *, cache_dir: Optional[Path] = None
) -> Optional[str]:
    """Return the cached reference protein sequence for ``gene``, or None.

    Reads ``<cache_dir>/<GENE>.fasta`` (default ``data/reference_sequences/``).
    Results are memoized per process. When a curated expected length is known
    (``utils.variant_normalizer.PROTEIN_LENGTHS``) and the loaded sequence does
    not match it, a warning is logged and the sequence is still returned (the
    caller decides; length-mismatch usually means a stale/wrong accession).
    """
    key = (gene or "").upper()
    if cache_dir is None and key in _SEQ_CACHE:
        return _SEQ_CACHE[key]

    base = cache_dir or _REF_DIR
    seq = _read_fasta_sequence(base / f"{key}.fasta")

    if seq is not None:
        try:
            from utils.variant_normalizer import PROTEIN_LENGTHS

            expected = PROTEIN_LENGTHS.get(key)
        except Exception:
            expected = None
        if expected and len(seq) != expected:
            logger.warning(
                "reference sequence for %s has length %d but expected %d "
                "(stale/wrong accession?)",
                key,
                len(seq),
                expected,
            )
    if cache_dir is None:
        _SEQ_CACHE[key] = seq
    return seq


def validate_protein_variant(
    gene: str,
    protein_notation: str,
    *,
    sequence: Optional[str] = None,
    cache_dir: Optional[Path] = None,
) -> ReferenceVerdict:
    """Validate a protein variant's reference residue against the transcript.

    ``sequence`` may be supplied directly (e.g. in tests); otherwise it is loaded
    from the per-gene cache. Returns ``unknown`` (never ``reject``) whenever there
    is no concrete residue to check or no reference sequence is available, so the
    gate is safe to run on every gene.
    """
    parsed = parse_reference_residue(protein_notation)
    if parsed is None:
        return ReferenceVerdict("unknown", "no_reference_residue_in_notation")
    ref_aa, pos = parsed

    if sequence is None:
        sequence = load_reference_protein(gene, cache_dir=cache_dir)
    if not sequence:
        return ReferenceVerdict("unknown", "no_reference_sequence")

    if pos > len(sequence):
        return ReferenceVerdict(
            "reject",
            "position_out_of_range",
            expected=str(len(sequence)),
            observed=str(pos),
        )
    actual = sequence[pos - 1].upper()
    if actual != ref_aa:
        return ReferenceVerdict(
            "reject", "reference_residue_mismatch", expected=actual, observed=ref_aa
        )
    return ReferenceVerdict("ok", "match", expected=actual, observed=ref_aa)


def is_simple_substitution(protein_notation: str) -> bool:
    """True for a clean single-residue substitution or nonsense change.

    e.g. ``A561V``, ``p.Arg190Trp``, ``R190X``, ``Y616*``, ``p.Arg190Ter``.
    False for frameshifts (``fs``), indels (``del``/``ins``/``dup``), and
    extensions (``ext``) — for those the leading residue is the first *affected*
    codon, whose identity is notation-convention dependent, so a hard
    reference-residue reject on them is low-confidence.
    """
    if not protein_notation:
        return False
    s = str(protein_notation)
    low = s.lower()
    if any(tok in low for tok in _COMPLEX_NOTATION_TOKENS):
        return False
    return bool(_SIMPLE_SUB_RE.match(s))


def is_high_confidence_reject(
    protein_notation: str,
    verdict: "ReferenceVerdict",
    *,
    nterm_guard: int = NTERM_ISOFORM_GUARD,
) -> bool:
    """True when a ``reject`` verdict is safe for a destructive (drop) policy.

    A reject is high-confidence — and therefore droppable — only when it is a
    clean single-residue substitution/nonsense (:func:`is_simple_substitution`)
    at a position beyond the N-terminal isoform-divergence window
    (``nterm_guard``). Frameshift/indel rejects and N-terminal rejects are
    real but ambiguous (notation convention / alternative isoform numbering),
    so they should be *flagged but kept*, never dropped. Non-reject verdicts
    are never droppable.
    """
    if not verdict.is_reject:
        return False
    parsed = parse_reference_residue(protein_notation)
    if parsed is None:
        return False
    _ref, pos = parsed
    if pos <= nterm_guard:
        return False
    return is_simple_substitution(protein_notation)


def clear_cache() -> None:
    """Drop the in-process sequence cache (tests / after populating the cache)."""
    _SEQ_CACHE.clear()
