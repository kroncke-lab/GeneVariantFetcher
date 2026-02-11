"""
DEPRECATED: Use utils.variant_normalizer instead.

This module is kept for backward compatibility only.
All functions delegate to the canonical variant_normalizer module.

Migration guide:
    # Old:
    from utils.variant_utils import normalize_variant

    # New:
    from utils.variant_normalizer import normalize_variant

The new module provides additional features:
- Gene-specific validation (protein length checks)
- Non-target variant detection (TP53/KRAS hotspots)
- KCNH2 variant aliases from ClinVar
- Splice variant handling (IVS notation)
- create_variant_key() for aggregation
"""

import warnings
from typing import Optional, Tuple, List

# Import everything from canonical module
from utils.variant_normalizer import (
    normalize_variant,
    normalize_protein_variant,
    normalize_cdna_variant,
    variants_match,
    find_matching_variants,
    VariantNormalizer,
    AA_MAP as AA_THREE_TO_ONE_NEW,
    AA_MAP_REVERSE as AA_ONE_TO_THREE_NEW,
)

# Preserve old constant names for backward compatibility
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
    "Ter": "X",
    "Sec": "U",
    "Pyl": "O",
}
AA_ONE_TO_THREE = {v: k for k, v in AA_THREE_TO_ONE.items()}


def extract_variant_components(variant: str) -> Optional[Tuple[str, int, str]]:
    """
    Extract components from a protein variant.

    Args:
        variant: Protein variant notation

    Returns:
        Tuple of (ref_aa, position, alt_aa) or None if parsing fails
    """
    import re

    normalized = normalize_protein_variant(variant)
    if not normalized:
        return None

    # Match pattern like A561V
    match = re.match(r"^([A-Z])(\d+)([A-Z]|FS|DEL|DUP|INS|X)$", normalized.upper())
    if match:
        return (match.group(1), int(match.group(2)), match.group(3))
    return None


def is_protein_variant(variant: str) -> bool:
    """Check if variant appears to be a protein-level notation."""
    import re

    if not variant:
        return False
    variant = variant.strip().lower()
    if variant.startswith("p."):
        return True
    if re.match(r"^[a-z]{3}\d+", variant):  # Ala561...
        return True
    if re.match(r"^[a-z]\d+[a-z*]$", variant):  # A561V
        return True
    return False


def is_cdna_variant(variant: str) -> bool:
    """Check if variant appears to be a cDNA-level notation."""
    import re

    if not variant:
        return False
    variant = variant.strip().lower()
    if variant.startswith("c."):
        return True
    if re.match(r"^\d+[acgt]>[acgt]", variant):
        return True
    return False


# Export all public names
__all__ = [
    "normalize_variant",
    "normalize_protein_variant",
    "normalize_cdna_variant",
    "variants_match",
    "find_matching_variants",
    "extract_variant_components",
    "is_protein_variant",
    "is_cdna_variant",
    "AA_THREE_TO_ONE",
    "AA_ONE_TO_THREE",
]
