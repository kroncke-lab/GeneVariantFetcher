"""
Variant notation utilities for normalization and comparison.

Handles conversion between different variant notation formats:
- HGVS protein: p.Ala561Val, p.A561V
- Short protein: A561V, Ala561Val
- HGVS cDNA: c.1682C>T
- Legacy notations: A561V, 561A>V
"""

import re
from typing import Optional, Tuple

# Three-letter to one-letter amino acid mapping
AA_THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "X",  # Stop codon
    "Sec": "U",  # Selenocysteine
    "Pyl": "O",  # Pyrrolysine
}

# One-letter to three-letter (reverse mapping)
AA_ONE_TO_THREE = {v: k for k, v in AA_THREE_TO_ONE.items()}


def normalize_protein_variant(variant: str) -> str:
    """
    Normalize a protein variant to short format (e.g., A561V).

    Args:
        variant: Protein variant in any format (p.Ala561Val, A561V, etc.)

    Returns:
        Normalized variant in short format (A561V)
    """
    if not variant:
        return ""

    variant = variant.strip()

    # Remove p. prefix
    variant = re.sub(r"^p\.", "", variant, flags=re.IGNORECASE)

    # Convert three-letter amino acids to one-letter
    # Handle format like Ala561Val
    for three, one in AA_THREE_TO_ONE.items():
        # Case-insensitive replacement
        variant = re.sub(three, one, variant, flags=re.IGNORECASE)

    # Handle stop codons: * -> X, Ter -> X
    variant = variant.replace("*", "X")
    variant = re.sub(r"Ter", "X", variant, flags=re.IGNORECASE)

    # Uppercase the result
    variant = variant.upper()

    # Handle frameshift notation - normalize all fs variants to FSX format
    # A100fsX10 -> A100FS, A100fs -> A100FS, A1017fsX -> A1017FS
    variant = re.sub(r"FS[X*]?\d*$", "FS", variant)
    variant = re.sub(r"FSX$", "FS", variant)

    return variant


def normalize_cdna_variant(variant: str) -> str:
    """
    Normalize a cDNA variant notation.

    Args:
        variant: cDNA variant (c.1682C>T, 1682C>T, etc.)

    Returns:
        Normalized cDNA variant (c.1682C>T)
    """
    if not variant:
        return ""

    variant = variant.strip()

    # Ensure c. prefix (lowercase)
    if variant.lower().startswith("c."):
        variant = "c." + variant[2:]
    elif re.match(r"^\d+[ACGT]", variant, re.IGNORECASE):
        # Looks like a cDNA variant without prefix
        variant = "c." + variant
    elif re.match(r"^-?\d+[+-]\d+", variant):  # Intronic variants
        variant = "c." + variant

    # Uppercase nucleotides but keep c. lowercase
    if variant.startswith("c."):
        rest = variant[2:]
        rest = re.sub(r"([acgt])", lambda m: m.group(1).upper(), rest)
        variant = "c." + rest

    return variant


def normalize_variant(variant: str) -> str:
    """
    Normalize any variant notation to a comparable format.

    Protein variants are normalized to short format (A561V).
    cDNA variants are normalized with c. prefix.

    Args:
        variant: Variant in any format

    Returns:
        Normalized variant string
    """
    if not variant:
        return ""

    variant = variant.strip()

    # Detect variant type and normalize accordingly
    if variant.lower().startswith("p.") or re.match(r"^[A-Z][a-z]{2}\d+", variant):
        # Protein variant with p. prefix or 3-letter AA
        return normalize_protein_variant(variant)
    elif variant.lower().startswith("c.") or re.match(r"^\d+[ACGT]", variant):
        # cDNA variant
        return normalize_cdna_variant(variant)
    elif re.match(r"^[A-Z]\d+[A-Z*]$", variant):
        # Short protein format (A561V)
        return variant.upper().replace("*", "X")
    elif re.match(r"^[A-Z]\d+", variant, re.IGNORECASE):
        # Protein variant that needs normalization (frameshifts, deletions, etc.)
        # e.g., A83fsX, A100del, G1006dup
        return normalize_protein_variant(variant)
    else:
        # Unknown format - return uppercase
        return variant.upper()


def variants_match(v1: str, v2: str) -> bool:
    """
    Check if two variant notations refer to the same variant.

    Args:
        v1: First variant notation
        v2: Second variant notation

    Returns:
        True if variants match after normalization
    """
    norm1 = normalize_variant(v1)
    norm2 = normalize_variant(v2)

    if not norm1 or not norm2:
        return False

    return norm1 == norm2


def extract_variant_components(variant: str) -> Optional[Tuple[str, int, str]]:
    """
    Extract components from a protein variant.

    Args:
        variant: Protein variant notation

    Returns:
        Tuple of (ref_aa, position, alt_aa) or None if parsing fails
    """
    normalized = normalize_protein_variant(variant)

    if not normalized:
        return None

    # Match pattern like A561V
    match = re.match(r"^([A-Z])(\d+)([A-Z]|FS|DEL|DUP|INS)$", normalized)
    if match:
        return (match.group(1), int(match.group(2)), match.group(3))

    return None


def is_protein_variant(variant: str) -> bool:
    """Check if variant appears to be a protein-level notation."""
    if not variant:
        return False

    variant = variant.strip().lower()

    # Check for protein indicators
    if variant.startswith("p."):
        return True

    # Check for amino acid patterns
    if re.match(r"^[a-z]{3}\d+", variant):  # Ala561...
        return True

    if re.match(r"^[a-z]\d+[a-z*]$", variant):  # A561V
        return True

    return False


def is_cdna_variant(variant: str) -> bool:
    """Check if variant appears to be a cDNA-level notation."""
    if not variant:
        return False

    variant = variant.strip().lower()

    # Check for cDNA indicators
    if variant.startswith("c."):
        return True

    # Check for nucleotide patterns without prefix
    if re.match(r"^\d+[acgt]>[acgt]", variant):
        return True

    return False


# Convenience function for comparing variant lists
def find_matching_variants(extracted: list, expected: list) -> Tuple[list, list, list]:
    """
    Compare extracted variants against expected variants.

    Args:
        extracted: List of extracted variant strings
        expected: List of expected variant strings

    Returns:
        Tuple of (matched, missed, extra) where:
        - matched: List of (extracted, expected) tuples that match
        - missed: List of expected variants not found in extracted
        - extra: List of extracted variants not found in expected
    """
    # Normalize all variants
    norm_extracted = {normalize_variant(v): v for v in extracted if v}
    norm_expected = {normalize_variant(v): v for v in expected if v}

    matched = []
    missed = []
    extra = []

    # Find matches
    for norm_exp, orig_exp in norm_expected.items():
        if norm_exp in norm_extracted:
            matched.append((norm_extracted[norm_exp], orig_exp))
        else:
            missed.append(orig_exp)

    # Find extras
    for norm_ext, orig_ext in norm_extracted.items():
        if norm_ext not in norm_expected:
            extra.append(orig_ext)

    return matched, missed, extra
