"""
Variant Normalizer for GeneVariantFetcher

CANONICAL NORMALIZER - Use this module for all variant normalization.
(Consolidated from variant_utils.py and improved_variant_normalizer.py)

Normalizes variant nomenclature to standard HGVS-style forms and validates
positions against known gene/protein lengths.

Features:
- Protein normalization: p.Ala561Val ↔ A561V (both directions)
- cDNA normalization with prefix handling and intronic variants
- Splice variant handling (IVS notation)
- Frameshift normalization (fs, fsX, fs*, fsX10 → fs*)
- Non-target variant detection (TP53/KRAS/BRAF/PIK3CA hotspots)
- Fuzzy position matching for off-by-one errors
- KCNH2 variant aliases from ClinVar
- Gene-specific validation (protein length, position bounds)

Usage:
    # Standalone function (drop-in replacement for variant_utils)
    from utils.variant_normalizer import normalize_variant
    normalized = normalize_variant("p.Ala561Val")  # Returns "A561V"

    # Class-based for gene-specific validation
    from utils.variant_normalizer import VariantNormalizer
    norm = VariantNormalizer("KCNH2")
    result = norm.get_all_forms("A561V")
"""

import re
import logging
from typing import Optional, Dict, Any, Tuple, Set, List

logger = logging.getLogger(__name__)

# Amino acid single-letter to three-letter code mapping
AA_MAP = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
    "*": "Ter",
    "X": "Xaa",
}

# Reverse mapping (three-letter to single-letter)
AA_MAP_REVERSE = {v: k for k, v in AA_MAP.items()}
AA_MAP_REVERSE["Ter"] = "*"
AA_MAP_REVERSE["Stop"] = "*"
AA_MAP_REVERSE["Xaa"] = "X"  # Unknown amino acid

# Protein lengths for common cardiac genes
PROTEIN_LENGTHS = {
    "KCNH2": 1159,
    "KCNQ1": 676,
    "SCN5A": 2016,
    "KCNE1": 129,
    "KCNE2": 123,
    "KCNJ2": 427,
    "CACNA1C": 2221,
    "RYR2": 4967,
}

# Non-target gene hotspots (TP53, KRAS, BRAF, PIK3CA)
# These are commonly extracted by mistake when searching for cardiac gene variants
NON_TARGET_HOTSPOTS = {
    # TP53 (most common cancer gene mutations)
    "R175H",
    "R248W",
    "R248Q",
    "R249S",
    "R273H",
    "R273C",
    "R282W",
    "G245S",
    "G245D",
    "Y220C",
    "C176F",
    "C176Y",
    "C242F",
    "C242S",
    "H179R",
    "H179Y",
    "M237I",
    "S241F",
    "S241C",
    "C277F",
    "Y163C",
    # KRAS
    "G12D",
    "G12V",
    "G12C",
    "G12R",
    "G12A",
    "G12S",
    "G13D",
    "Q61H",
    "Q61L",
    "Q61R",
    "Q61G",
    # BRAF
    "V600E",
    "V600K",
    # PIK3CA
    "E545K",
    "H1047R",
    "H1047L",
}

# KCNH2 variant aliases from ClinVar - maps alternative names to canonical form
# Key: canonical single-letter form, Values: alternative representations
KCNH2_VARIANT_ALIASES = {
    # Well-known founder mutations with multiple names
    "A561V": ["A561V", "p.Ala561Val", "Ala561Val", "c.1682C>T"],
    "R534C": ["R534C", "p.Arg534Cys", "Arg534Cys", "c.1600C>T"],
    "A614V": ["A614V", "p.Ala614Val", "Ala614Val", "c.1841C>T"],
    "G628S": ["G628S", "p.Gly628Ser", "Gly628Ser", "c.1882G>A"],
    "N470D": ["N470D", "p.Asn470Asp", "Asn470Asp", "c.1408A>G"],
    "R176W": ["R176W", "p.Arg176Trp", "Arg176Trp", "c.526C>T"],
    "T613M": ["T613M", "p.Thr613Met", "Thr613Met", "c.1838C>T"],
    "A422T": ["A422T", "p.Ala422Thr", "Ala422Thr", "c.1264G>A"],
    "R784W": ["R784W", "p.Arg784Trp", "Arg784Trp", "c.2350C>T"],
    "Y611H": ["Y611H", "p.Tyr611His", "Tyr611His", "c.1831T>C"],
    # Common splice variants
    "c.2398+1G>A": ["c.2398+1G>A", "IVS10+1G>A"],
    "c.1129-1G>A": ["c.1129-1G>A", "IVS5-1G>A"],
}

# Build reverse lookup for aliases from hardcoded dict
_KCNH2_ALIAS_LOOKUP = {}
for canonical, aliases in KCNH2_VARIANT_ALIASES.items():
    for alias in aliases:
        _KCNH2_ALIAS_LOOKUP[alias.upper()] = canonical

# Load comprehensive alias dictionary from JSON if available
_KCNH2_COMPREHENSIVE_ALIASES = {}
try:
    import json
    from pathlib import Path

    _alias_path = Path(__file__).parent / "kcnh2_variant_aliases.json"
    if _alias_path.exists():
        with open(_alias_path) as f:
            _alias_data = json.load(f)
            _KCNH2_COMPREHENSIVE_ALIASES = _alias_data.get("aliases", {})
            logger.debug(
                f"Loaded {len(_KCNH2_COMPREHENSIVE_ALIASES)} KCNH2 aliases from JSON"
            )
except Exception as e:
    logger.warning(f"Could not load KCNH2 alias dictionary: {e}")


def _lookup_alias(variant: str, gene: str = "KCNH2") -> Optional[str]:
    """
    Look up a variant in the comprehensive alias dictionary.
    Returns the canonical form if found, None otherwise.
    """
    if gene.upper() != "KCNH2":
        return None

    v_upper = variant.upper().strip()

    # Try comprehensive aliases first (from JSON)
    if v_upper in _KCNH2_COMPREHENSIVE_ALIASES:
        return _KCNH2_COMPREHENSIVE_ALIASES[v_upper]

    # Fall back to hardcoded aliases
    if v_upper in _KCNH2_ALIAS_LOOKUP:
        return _KCNH2_ALIAS_LOOKUP[v_upper]

    return None


# IVS to cDNA splice notation mapping for KCNH2
# This is a partial list - in practice would need full gene coordinates
KCNH2_IVS_MAP = {
    "IVS1": "c.175",
    "IVS2": "c.453",
    "IVS3": "c.575",
    "IVS4": "c.787",
    "IVS5": "c.1129",
    "IVS6": "c.1351",
    "IVS7": "c.1592",
    "IVS8": "c.1759",
    "IVS9": "c.2006",
    "IVS10": "c.2398",
    "IVS11": "c.2599",
    "IVS12": "c.2769",
    "IVS13": "c.3040",
}

# Patterns for variant parsing
PROTEIN_PATTERNS = [
    # Already normalized three-letter: p.Ala561Val
    re.compile(
        r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|[fs\*]|[delins\*])$", re.IGNORECASE
    ),
    # Three-letter without prefix: Ala561Val
    re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", re.IGNORECASE),
    # Single-letter conversion: A561V -> p.Ala561Val
    re.compile(r"^(?:p\.)?([A-Z])(\d+)([A-Z])$", re.IGNORECASE),
    # Frame-shift with single-letter: A193fsX, A193fsX10, A193fs*10
    re.compile(r"^(?:p\.)?([A-Z])(\d+)(fsX\d*|fs\*\d*|fs)$", re.IGNORECASE),
    # Frame-shift with three-letter: Ala193fsX, Ala193fs*10
    re.compile(r"^(?:p\.)?([A-Z][a-z]{2})(\d+)(fsX\d*|fs\*\d*|fs)$", re.IGNORECASE),
    # Insertion/deletion: G184Del, G184Ins, G184_G185del
    re.compile(r"^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(Del|Ins|Dup)$", re.IGNORECASE),
    # Stop/truncation: R864stop, D864sp, R864*, R864X, R864Ter
    re.compile(
        r"^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(sp|stop|trunc|\*|Ter|X)$", re.IGNORECASE
    ),
    # Extended deletion/dup with range: p.Gly184_Leu186del
    re.compile(
        r"^(?:p\.)?([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)(del|dup|ins)$",
        re.IGNORECASE,
    ),
]

CDNA_PATTERNS = [
    # c.1234A>G
    re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$"),
    # c.1234delA, c.1234dupG
    re.compile(r"^c\.(\d+)(del|dup|ins)([ACGT]*)$", re.IGNORECASE),
    # c.1234_1235delAG
    re.compile(r"^c\.(\d+)_(\d+)(del|dup|ins)([ACGT]*)$", re.IGNORECASE),
    # c.1234+1G>A (intronic - donor)
    re.compile(r"^c\.(\d+[\+]\d+)([ACGT])>([ACGT])$"),
    # c.1234-1G>A (intronic - acceptor)
    re.compile(r"^c\.(\d+[\-]\d+)([ACGT])>([ACGT])$"),
    # c.1234+1del or c.1234-1del (intronic indels)
    re.compile(r"^c\.(\d+[\+\-]\d+)(del|dup|ins)([ACGT]*)$", re.IGNORECASE),
]

# IVS (splice) patterns - legacy notation
SPLICE_PATTERNS = [
    # IVS10+1G>A, IVS5-2A>G
    re.compile(r"^IVS(\d+)([\+\-]\d+)([ACGT])>([ACGT])$", re.IGNORECASE),
    # IVS10+1del
    re.compile(r"^IVS(\d+)([\+\-]\d+)(del|dup|ins)([ACGT]*)$", re.IGNORECASE),
]


class VariantNormalizer:
    """Normalize variant nomenclature to standard forms."""

    def __init__(self, gene_symbol: str):
        """
        Initialize normalizer for a specific gene.

        Args:
            gene_symbol: Gene symbol (e.g., 'KCNH2')
        """
        self.gene_symbol = gene_symbol.upper()
        self.protein_length = PROTEIN_LENGTHS.get(self.gene_symbol)

    def normalize_protein(self, variant: str) -> Optional[str]:
        """
        Normalize protein variant to p.Xxx###Yyy format.

        Args:
            variant: Variant string in various formats including fs, del, ins patterns

        Returns:
            Normalized variant string or None if unparseable

        Examples:
            'A561V' -> 'p.Ala561Val'
            'p.A561V' -> 'p.Ala561Val'
            'Ala561Val' -> 'p.Ala561Val'
            'A193fsX' -> 'p.Ala193fs*'
            'G184Del' -> 'p.Gly184del'
            'R864sp' -> 'p.Arg864*'
        """
        if not variant:
            return None

        variant = variant.strip()
        original_variant = variant

        for pattern in PROTEIN_PATTERNS:
            match = pattern.match(variant)
            if match:
                groups = match.groups()
                if len(groups) != 3:
                    continue

                ref, pos, suffix = groups

                # Convert reference to three-letter code
                ref_3 = self._to_three_letter(ref.upper())
                if not ref_3:
                    continue

                suffix_upper = suffix.upper()

                # Handle different variant types
                if suffix_upper.startswith("FS"):
                    # Frame-shift mutations
                    normalized_suffix = "fs*"
                elif suffix_upper in ["X", "TER"]:
                    # Stop codon
                    normalized_suffix = "*"
                elif suffix_upper in ["SP", "STOP", "TRUNC"]:
                    # Alternative stop annotations
                    normalized_suffix = "*"
                elif suffix_upper in ["DEL", "INS"]:
                    # Insertion/deletion
                    normalized_suffix = suffix_upper.lower()
                elif len(suffix) == 1 and suffix.isalpha():
                    # Single amino acid change
                    alt_3 = self._to_three_letter(suffix.upper())
                    if alt_3:
                        normalized_suffix = alt_3
                    else:
                        continue
                elif len(suffix) == 3 and suffix[0].isupper():
                    # Three-letter amino acid change
                    normalized_suffix = suffix.capitalize()
                elif "*" in suffix:
                    # Stop codon variants
                    normalized_suffix = "*"
                else:
                    continue

                return f"p.{ref_3}{pos}{normalized_suffix}"

        # If no pattern matched, return None
        logger.debug(f"Could not normalize protein variant: {original_variant}")
        return None

    def normalize_to_single_letter(self, variant: str) -> Optional[str]:
        """
        Normalize protein variant to single-letter format: A561V

        Handles: A561V, p.Ala561Val, Ala561Val, p.A561V

        Normalization rules:
        - Strip p. prefix
        - Remove trailing asterisks (DMS artifacts): A561V* → A561V
        - Remove trailing annotations: R176W(het) → R176W
        - Case normalization: a561v → A561V
        - 3-letter → 1-letter: Ala561Val → A561V
        - Frameshift standardization: fs*46 → fsX, fsXNN → fsX, Alafs → AfsX

        Args:
            variant: Variant string in various formats

        Returns:
            Single-letter format (e.g., 'A561V') or None if unparseable
        """
        if not variant:
            return None

        v = variant.strip()

        # Remove p. prefix
        if v.lower().startswith("p."):
            v = v[2:]

        # Remove trailing DMS asterisks (A561V* → A561V, but keep R864* as stop)
        if v.endswith("*") and len(v) >= 2 and v[-2].isalpha():
            v = v[:-1]

        # Remove trailing annotations
        v = re.sub(r"\([^)]*\)$", "", v)  # (het), (hom)
        v = re.sub(r"/[Ww][Tt]$", "", v)  # /WT
        v = re.sub(r"/[\+\-]$", "", v)  # /+, /-
        v = v.strip()

        # Single-letter format (handles both lowercase and uppercase): A561V, a561v, R864*, R864X
        m = re.match(r"^([A-Za-z])(\d+)([A-Za-z\*X])$", v, re.IGNORECASE)
        if m:
            ref = m.group(1).upper()
            pos = m.group(2)
            alt = m.group(3).upper()
            # Normalize * to X for stop codons
            if alt == "*":
                alt = "X"
            return f"{ref}{pos}{alt}"

        # Three-letter: Ala561Val, Arg864Ter, Arg864*
        m = re.match(
            r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Ter|X)$", v, re.IGNORECASE
        )
        if m:
            ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
            alt_raw = m.group(3)
            if alt_raw.upper() in ["*", "TER", "X"]:
                alt = "X"
            else:
                alt = AA_MAP_REVERSE.get(alt_raw.capitalize())
            if ref and alt:
                return f"{ref}{m.group(2)}{alt}"

        # === FRAMESHIFT STANDARDIZATION ===
        # All frameshift variants → REFposfsX (drop position after fs)
        # Patterns: fs*46, fsX10, fs*10, fs, fsXNN, fs+NX, fsTerN

        # Pattern 1: Three-letter with frameshift: Ala193fs*46, Gly262Alafs*98, Arg518fs*
        m = re.match(
            r"^([A-Z][a-z]{2})(\d+)(?:[A-Z][a-z]{2})?fs[\*X]?\d*$", v, re.IGNORECASE
        )
        if m:
            ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
            if ref:
                return f"{ref}{m.group(2)}fsX"

        # Pattern 2: Single-letter with frameshift: A193fsX, A193fs*46, A193fs
        m = re.match(r"^([A-Za-z])(\d+)fs[\*X]?\d*$", v, re.IGNORECASE)
        if m:
            return f"{m.group(1).upper()}{m.group(2)}fsX"

        # Pattern 3: fs+NX format: A1077fs+X*, A193fs+34X
        m = re.match(r"^([A-Za-z])(\d+)fs\+?\d*X?\*?$", v, re.IGNORECASE)
        if m:
            return f"{m.group(1).upper()}{m.group(2)}fsX"

        # Pattern 4: fsTer format: Gly24fsTer58
        m = re.match(
            r"^([A-Z][a-z]{2})(\d+)(?:[A-Z][a-z]{2})?fsTer\d*$", v, re.IGNORECASE
        )
        if m:
            ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
            if ref:
                return f"{ref}{m.group(2)}fsX"

        # Pattern 5: Three-letter stop: Arg518*
        m = re.match(r"^([A-Z][a-z]{2})(\d+)\*$", v, re.IGNORECASE)
        if m:
            ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
            if ref:
                return f"{ref}{m.group(2)}X"

        # Stop variants: R864sp, R864stop, R864Ter
        m = re.match(r"^([A-Za-z])(\d+)(sp|stop|Ter)$", v, re.IGNORECASE)
        if m:
            return f"{m.group(1).upper()}{m.group(2)}X"

        return None

    def normalize_cdna(self, variant: str) -> Optional[str]:
        """
        Normalize cDNA variant to c.###X>Y format.

        Enhanced to handle missing c. prefix and intronic variants.

        Args:
            variant: Variant string

        Returns:
            Normalized variant string or None if unparseable
        """
        if not variant:
            return None

        variant = variant.strip()

        # Already has c. prefix - check if valid pattern
        if variant.startswith("c."):
            for pattern in CDNA_PATTERNS:
                if pattern.match(variant):
                    return variant
            # Return as-is if it starts with c. (may be valid but complex)
            return variant

        # Try adding c. prefix if missing
        prefixed = f"c.{variant}"
        for pattern in CDNA_PATTERNS:
            if pattern.match(prefixed):
                return prefixed

        # Check if it looks like cDNA (numbers and nucleotides) - more flexible
        if re.match(r"^\d+[\+\-\d]*[ACGT_>delinsdup]", variant):
            return f"c.{variant}"

        logger.debug(f"Could not normalize cDNA variant: {variant}")
        return None

    def validate_position(self, position: int) -> bool:
        """
        Check if position is valid for this gene's protein.

        Args:
            position: Amino acid position

        Returns:
            True if valid, False otherwise
        """
        if self.protein_length is None:
            # Unknown gene, can't validate
            return True

        return 1 <= position <= self.protein_length

    def extract_position(self, variant: str) -> Optional[int]:
        """
        Extract amino acid position from protein variant.

        Args:
            variant: Variant string

        Returns:
            Position as integer, or None if not extractable
        """
        single = self.normalize_to_single_letter(variant)
        if single:
            m = re.match(r"^[A-Z](\d+)", single)
            if m:
                return int(m.group(1))

        # Fallback: try to extract from three-letter or p. notation
        m = re.search(r"(\d+)", variant)
        if m:
            return int(m.group(1))

        return None

    def is_non_target_variant(self, variant: str) -> Tuple[bool, Optional[str]]:
        """
        Check if variant is likely from a non-target gene.

        Detects common TP53/KRAS/BRAF/PIK3CA hotspots and variants
        with positions exceeding the target gene's protein length.

        Args:
            variant: Variant string

        Returns:
            Tuple of (is_non_target, reason)
        """
        single = self.normalize_to_single_letter(variant)
        if not single:
            return False, None

        # Check hotspots (remove fsX suffix for comparison)
        check_var = re.sub(r"fsX$", "", single)
        if check_var in NON_TARGET_HOTSPOTS:
            return True, f"Known hotspot in TP53/KRAS/BRAF/PIK3CA"

        # Check position against protein length
        m = re.match(r"^[A-Z](\d+)", single)
        if m:
            pos = int(m.group(1))
            if self.protein_length and pos > self.protein_length:
                return (
                    True,
                    f"Position {pos} exceeds {self.gene_symbol} length ({self.protein_length})",
                )

        return False, None

    def get_all_forms(self, variant: str) -> Dict[str, str]:
        """
        Get all normalized forms of a variant.

        Args:
            variant: Variant string

        Returns:
            Dict with keys: original, single, three, cdna (as available)
        """
        result = {"original": variant}

        single = self.normalize_to_single_letter(variant)
        if single:
            result["single"] = single

        three = self.normalize_protein(variant)
        if three:
            result["three"] = three

        cdna = self.normalize_cdna(variant)
        if cdna:
            result["cdna"] = cdna

        return result

    def get_canonical_form(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Return normalized variant with all fields populated.

        Args:
            variant: Variant dict with protein_notation, cdna_notation, etc.

        Returns:
            Enhanced variant dict with normalized fields
        """
        result = variant.copy()

        # Normalize protein notation
        protein = variant.get("protein_notation") or variant.get("protein_change")
        if protein:
            normalized = self.normalize_protein(protein)
            if normalized:
                result["protein_notation"] = normalized
                result["protein_notation_normalized"] = True

                # Extract position and validate
                match = re.search(r"(\d+)", normalized)
                if match:
                    pos = int(match.group(1))
                    result["position"] = pos
                    result["position_valid"] = self.validate_position(pos)

                # Check if non-target
                is_non_target, reason = self.is_non_target_variant(protein)
                if is_non_target:
                    result["is_non_target"] = True
                    result["non_target_reason"] = reason
            else:
                result["protein_notation_normalized"] = False

        # Normalize cDNA notation
        cdna = variant.get("cdna_notation") or variant.get("cdna_change")
        if cdna:
            normalized = self.normalize_cdna(cdna)
            if normalized:
                result["cdna_notation"] = normalized
                result["cdna_notation_normalized"] = True
            else:
                result["cdna_notation_normalized"] = False

        # Create canonical key for deduplication
        result["canonical_key"] = self._create_canonical_key(result)

        return result

    def _to_three_letter(self, aa: str) -> Optional[str]:
        """Convert amino acid to three-letter code."""
        if len(aa) == 1:
            return AA_MAP.get(aa.upper())
        elif len(aa) == 3:
            # Already three-letter, just capitalize properly
            return aa.capitalize()
        return None

    def _to_single_letter(self, aa: str) -> Optional[str]:
        """Convert amino acid to single-letter code."""
        if len(aa) == 1:
            return aa.upper()
        elif len(aa) == 3:
            return AA_MAP_REVERSE.get(aa.capitalize())
        return None

    def _create_canonical_key(self, variant: Dict[str, Any]) -> str:
        """Create a canonical key for variant deduplication."""
        gene = self.gene_symbol
        protein = variant.get("protein_notation", "")
        cdna = variant.get("cdna_notation", "")

        # Prefer protein notation for key
        if protein:
            return f"{gene}:{protein}"
        elif cdna:
            return f"{gene}:{cdna}"
        else:
            return f"{gene}:unknown"


def normalize_variant_list(variants: list, gene_symbol: str) -> Tuple[list, dict]:
    """
    Normalize a list of variants and deduplicate.

    Args:
        variants: List of variant dicts
        gene_symbol: Gene symbol for validation

    Returns:
        Tuple of (normalized_variants, stats_dict)
    """
    normalizer = VariantNormalizer(gene_symbol)

    seen_keys = set()
    normalized = []
    stats = {
        "total_input": len(variants),
        "normalized": 0,
        "duplicates_removed": 0,
        "normalization_failed": 0,
        "non_target_filtered": 0,
    }

    for v in variants:
        result = normalizer.get_canonical_form(v)
        key = result.get("canonical_key", "")

        if key in seen_keys:
            stats["duplicates_removed"] += 1
            continue

        # Track non-target variants
        if result.get("is_non_target"):
            stats["non_target_filtered"] += 1

        seen_keys.add(key)
        normalized.append(result)

        if result.get("protein_notation_normalized") or result.get(
            "cdna_notation_normalized"
        ):
            stats["normalized"] += 1
        else:
            stats["normalization_failed"] += 1

    stats["total_output"] = len(normalized)

    return normalized, stats


def normalize_frameshift(variant: str) -> Optional[str]:
    """
    Normalize all frameshift naming conventions to a standard form: REFposfsX

    Handles:
    - L987fs, L987fsX, L987fsX10, L987fs*10 → L987fsX
    - p.Leu987fs, p.Leu987fsTer10, p.Leu987Profs*10 → L987fsX
    - 987fs, fs987 → position extracted if possible

    Args:
        variant: Frameshift variant in any notation

    Returns:
        Normalized form (L987fsX) or None if not a frameshift
    """
    if not variant:
        return None

    v = variant.strip()

    # Remove p. prefix
    if v.lower().startswith("p."):
        v = v[2:]

    # Pattern 1: Single letter + position + fs variants
    # L987fs, L987fsX, L987fsX10, L987fs*10, L987fs*
    m = re.match(r"^([A-Za-z])(\d+)fs[\*X]?\d*$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}fsX"

    # Pattern 2: Three-letter + position + optional second AA + fs
    # Leu987fs, Leu987Profs, Leu987Profs*10, Leu987fsTer10
    m = re.match(
        r"^([A-Z][a-z]{2})(\d+)(?:[A-Z][a-z]{2})?fs[\*X]?\d*$", v, re.IGNORECASE
    )
    if m:
        ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
        if ref:
            return f"{ref}{m.group(2)}fsX"

    # Pattern 3: fsTer format (Leu987fsTer10, L987fsTer)
    m = re.match(
        r"^([A-Za-z]|[A-Z][a-z]{2})(\d+)(?:[A-Za-z]{0,3})?fsTer\d*$", v, re.IGNORECASE
    )
    if m:
        ref_raw = m.group(1)
        if len(ref_raw) == 1:
            ref = ref_raw.upper()
        else:
            ref = AA_MAP_REVERSE.get(ref_raw.capitalize())
        if ref:
            return f"{ref}{m.group(2)}fsX"

    # Pattern 4: Position-only fs (987fs, fs987) - less reliable
    m = re.match(r"^(\d+)fs[\*X]?\d*$", v, re.IGNORECASE)
    if m:
        return f"?{m.group(1)}fsX"  # Mark as ambiguous

    m = re.match(r"^fs(\d+)$", v, re.IGNORECASE)
    if m:
        return f"?{m.group(1)}fsX"  # Mark as ambiguous

    return None


def normalize_nonsense(variant: str) -> Optional[str]:
    """
    Normalize nonsense/stop variants to standard form: REFposX

    Handles:
    - W1001X, W1001*, p.Trp1001Ter, p.Trp1001* → W1001X
    - R864stop, R864sp → R864X

    Args:
        variant: Stop/nonsense variant in any notation

    Returns:
        Normalized form (W1001X) or None if not a nonsense variant
    """
    if not variant:
        return None

    v = variant.strip()

    # Remove p. prefix
    if v.lower().startswith("p."):
        v = v[2:]

    # Pattern 1: Single letter + position + stop indicators
    # W1001X, W1001*, R864X
    m = re.match(r"^([A-Za-z])(\d+)[\*X]$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}X"

    # Pattern 2: Three-letter + position + Ter/*
    # Trp1001Ter, Trp1001*, Arg864Ter
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(Ter|\*)$", v, re.IGNORECASE)
    if m:
        ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
        if ref:
            return f"{ref}{m.group(2)}X"

    # Pattern 3: stop/sp suffix
    # R864stop, R864sp
    m = re.match(r"^([A-Za-z])(\d+)(stop|sp)$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}X"

    # Pattern 4: Three-letter with stop
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(stop|sp)$", v, re.IGNORECASE)
    if m:
        ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
        if ref:
            return f"{ref}{m.group(2)}X"

    return None


def normalize_deletion(variant: str) -> Optional[str]:
    """
    Normalize deletion variants to standard form: REFposdel

    Handles:
    - del552, p.Leu552del, L552del → L552del
    - L552_L555del → L552_L555del (range preserved)

    Args:
        variant: Deletion variant in any notation

    Returns:
        Normalized form (L552del) or None if not a deletion
    """
    if not variant:
        return None

    v = variant.strip()

    # Remove p. prefix
    if v.lower().startswith("p."):
        v = v[2:]

    # Pattern 1: Single letter + position + del
    # L552del
    m = re.match(r"^([A-Za-z])(\d+)del$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}del"

    # Pattern 2: Three-letter + position + del
    # Leu552del
    m = re.match(r"^([A-Z][a-z]{2})(\d+)del$", v, re.IGNORECASE)
    if m:
        ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
        if ref:
            return f"{ref}{m.group(2)}del"

    # Pattern 3: Position-only del
    # del552
    m = re.match(r"^del(\d+)$", v, re.IGNORECASE)
    if m:
        return f"?{m.group(1)}del"  # Mark as ambiguous

    # Pattern 4: Range deletion (preserve range)
    # L552_L555del, Leu552_Leu555del
    m = re.match(
        r"^([A-Za-z]|[A-Z][a-z]{2})(\d+)_([A-Za-z]|[A-Z][a-z]{2})(\d+)del$",
        v,
        re.IGNORECASE,
    )
    if m:
        ref1_raw, pos1, ref2_raw, pos2 = m.groups()
        if len(ref1_raw) == 1:
            ref1 = ref1_raw.upper()
        else:
            ref1 = AA_MAP_REVERSE.get(ref1_raw.capitalize(), ref1_raw[0].upper())
        if len(ref2_raw) == 1:
            ref2 = ref2_raw.upper()
        else:
            ref2 = AA_MAP_REVERSE.get(ref2_raw.capitalize(), ref2_raw[0].upper())
        return f"{ref1}{pos1}_{ref2}{pos2}del"

    return None


def normalize_duplication(variant: str) -> Optional[str]:
    """
    Normalize duplication variants to standard form: REFposdup

    Handles:
    - dup552, p.Leu552dup, L552dup → L552dup

    Args:
        variant: Duplication variant in any notation

    Returns:
        Normalized form (L552dup) or None if not a duplication
    """
    if not variant:
        return None

    v = variant.strip()

    # Remove p. prefix
    if v.lower().startswith("p."):
        v = v[2:]

    # Pattern 1: Single letter + position + dup
    m = re.match(r"^([A-Za-z])(\d+)dup$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}dup"

    # Pattern 2: Three-letter + position + dup
    m = re.match(r"^([A-Z][a-z]{2})(\d+)dup$", v, re.IGNORECASE)
    if m:
        ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
        if ref:
            return f"{ref}{m.group(2)}dup"

    # Pattern 3: Position-only dup
    m = re.match(r"^dup(\d+)$", v, re.IGNORECASE)
    if m:
        return f"?{m.group(1)}dup"  # Mark as ambiguous

    # Pattern 4: Range duplication
    m = re.match(
        r"^([A-Za-z]|[A-Z][a-z]{2})(\d+)_([A-Za-z]|[A-Z][a-z]{2})(\d+)dup$",
        v,
        re.IGNORECASE,
    )
    if m:
        ref1_raw, pos1, ref2_raw, pos2 = m.groups()
        if len(ref1_raw) == 1:
            ref1 = ref1_raw.upper()
        else:
            ref1 = AA_MAP_REVERSE.get(ref1_raw.capitalize(), ref1_raw[0].upper())
        if len(ref2_raw) == 1:
            ref2 = ref2_raw.upper()
        else:
            ref2 = AA_MAP_REVERSE.get(ref2_raw.capitalize(), ref2_raw[0].upper())
        return f"{ref1}{pos1}_{ref2}{pos2}dup"

    return None


def get_variant_type(variant: str) -> str:
    """
    Determine the type of a protein variant.

    Returns one of: 'missense', 'nonsense', 'frameshift', 'deletion',
    'duplication', 'insertion', 'splice', 'unknown'
    """
    if not variant:
        return "unknown"

    v = variant.upper()

    if "FS" in v:
        return "frameshift"
    if v.endswith("X") or v.endswith("*") or "TER" in v.upper() or "STOP" in v.upper():
        return "nonsense"
    if "DEL" in v:
        return "deletion"
    if "DUP" in v:
        return "duplication"
    if "INS" in v:
        return "insertion"
    if "IVS" in v or "+" in v or "-" in v:
        return "splice"

    # Missense: single AA change
    if re.match(r"^[A-Z]\d+[A-Z]$", v):
        return "missense"

    return "unknown"


def match_variants_fuzzy(
    variant1: str,
    variant2: str,
    gene_symbol: str = "KCNH2",
    position_tolerance: int = 1,
) -> Tuple[bool, str]:
    """
    Check if two variants match using fuzzy matching rules.

    Tries in order:
    1. Exact match after normalization
    2. Type-specific matching (frameshift, nonsense, deletion, etc.)
    3. Position-based matching with tolerance for same variant type

    Args:
        variant1: First variant string
        variant2: Second variant string
        gene_symbol: Gene context for normalization
        position_tolerance: Position offset to allow (default ±1)

    Returns:
        Tuple of (is_match, match_type) where match_type describes how they matched
    """
    if not variant1 or not variant2:
        return False, "no_input"

    normalizer = VariantNormalizer(gene_symbol)

    # Step 1: Try exact match after standard normalization
    norm1 = normalize_variant(variant1, gene_symbol)
    norm2 = normalize_variant(variant2, gene_symbol)

    if norm1 and norm2 and norm1 == norm2:
        return True, "exact_normalized"

    # Step 2: Type-specific normalization
    type1 = get_variant_type(variant1)
    type2 = get_variant_type(variant2)

    # Frameshift matching
    if type1 == "frameshift" or type2 == "frameshift":
        fs1 = normalize_frameshift(variant1)
        fs2 = normalize_frameshift(variant2)
        if fs1 and fs2:
            # Remove ambiguous marker for comparison
            fs1_clean = fs1.replace("?", "")
            fs2_clean = fs2.replace("?", "")
            if fs1_clean == fs2_clean:
                return True, "frameshift_normalized"

    # Nonsense matching
    if type1 == "nonsense" or type2 == "nonsense":
        ns1 = normalize_nonsense(variant1)
        ns2 = normalize_nonsense(variant2)
        if ns1 and ns2 and ns1 == ns2:
            return True, "nonsense_normalized"

    # Deletion matching
    if type1 == "deletion" or type2 == "deletion":
        del1 = normalize_deletion(variant1)
        del2 = normalize_deletion(variant2)
        if del1 and del2:
            del1_clean = del1.replace("?", "")
            del2_clean = del2.replace("?", "")
            if del1_clean == del2_clean:
                return True, "deletion_normalized"

    # Duplication matching
    if type1 == "duplication" or type2 == "duplication":
        dup1 = normalize_duplication(variant1)
        dup2 = normalize_duplication(variant2)
        if dup1 and dup2:
            dup1_clean = dup1.replace("?", "")
            dup2_clean = dup2.replace("?", "")
            if dup1_clean == dup2_clean:
                return True, "duplication_normalized"

    # Step 3: Position-based fuzzy matching (for same type)
    if type1 == type2 and type1 != "unknown":
        pos1 = normalizer.extract_position(variant1)
        pos2 = normalizer.extract_position(variant2)

        if pos1 and pos2 and abs(pos1 - pos2) <= position_tolerance:
            # For missense, also check ref and alt amino acids
            if type1 == "missense":
                single1 = normalizer.normalize_to_single_letter(variant1)
                single2 = normalizer.normalize_to_single_letter(variant2)
                if single1 and single2:
                    m1 = re.match(r"^([A-Z])\d+([A-Z])$", single1)
                    m2 = re.match(r"^([A-Z])\d+([A-Z])$", single2)
                    if m1 and m2:
                        ref1, alt1 = m1.groups()
                        ref2, alt2 = m2.groups()
                        if ref1 == ref2 and alt1 == alt2:
                            return True, f"position_fuzzy_{pos1 - pos2:+d}"
            else:
                # For non-missense types, position match is sufficient
                return True, f"position_fuzzy_{pos1 - pos2:+d}"

    return False, "no_match"


def match_variants_to_baseline(
    extracted: List[str],
    baseline: Set[str],
    gene_symbol: str = "KCNH2",
    fuzzy_position: bool = True,
    position_tolerance: int = 1,
) -> Dict[str, Any]:
    """
    Match extracted variants to a baseline set with improved normalization.

    Supports exact matching, normalized matching, cDNA prefix matching,
    fuzzy position matching (±1 for off-by-one errors), and type-specific
    normalization for frameshifts, nonsense, deletions, and duplications.

    ENHANCED (2026-02-10): Now uses match_variants_fuzzy() for comprehensive
    variant matching including all frameshift conventions, nonsense variants,
    deletions, and duplications.

    Args:
        extracted: List of extracted variant strings
        baseline: Set of baseline variant strings
        gene_symbol: Target gene symbol
        fuzzy_position: Allow ±1 position matching
        position_tolerance: Position offset to allow (default 1)

    Returns:
        Dict with matches, unmatched, filtered_non_target, and stats
    """
    normalizer = VariantNormalizer(gene_symbol)

    # Build baseline position index for fuzzy matching
    baseline_by_pos = {}
    baseline_singles = set()
    baseline_normalized = {}  # Map normalized form -> original

    for v in baseline:
        single = normalizer.normalize_to_single_letter(v)
        if single:
            baseline_singles.add(single)
            baseline_normalized[single] = v
            pos = normalizer.extract_position(v)
            if pos:
                if pos not in baseline_by_pos:
                    baseline_by_pos[pos] = []
                baseline_by_pos[pos].append(v)

        # Also store type-specific normalized forms
        fs_norm = normalize_frameshift(v)
        if fs_norm:
            baseline_normalized[fs_norm.replace("?", "")] = v

        ns_norm = normalize_nonsense(v)
        if ns_norm:
            baseline_normalized[ns_norm] = v

        del_norm = normalize_deletion(v)
        if del_norm:
            baseline_normalized[del_norm.replace("?", "")] = v

        dup_norm = normalize_duplication(v)
        if dup_norm:
            baseline_normalized[dup_norm.replace("?", "")] = v

    # Also index cDNA forms
    baseline_cdna = set()
    for v in baseline:
        if v.startswith("c.") or v.startswith("IVS"):
            baseline_cdna.add(v)

    results = {
        "matches": [],
        "unmatched": [],
        "filtered_non_target": [],
        "stats": {
            "total_input": len(extracted),
            "exact_matches": 0,
            "normalized_matches": 0,
            "fuzzy_matches": 0,
            "cdna_matches": 0,
            "frameshift_matches": 0,
            "nonsense_matches": 0,
            "deletion_matches": 0,
            "duplication_matches": 0,
            "filtered": 0,
            "unmatched": 0,
        },
    }

    for v in extracted:
        # Check if non-target gene variant
        is_non_target, reason = normalizer.is_non_target_variant(v)
        if is_non_target:
            results["filtered_non_target"].append((v, reason))
            results["stats"]["filtered"] += 1
            continue

        matched_to = None
        match_type = None

        # Try exact match
        if v in baseline:
            matched_to = v
            match_type = "exact"
            results["stats"]["exact_matches"] += 1

        # Try normalized single-letter match
        if not matched_to:
            single = normalizer.normalize_to_single_letter(v)
            if single and single in baseline_singles:
                matched_to = baseline_normalized.get(single, single)
                match_type = "normalized"
                results["stats"]["normalized_matches"] += 1

        # Try type-specific normalized match
        if not matched_to:
            # Frameshift
            fs_norm = normalize_frameshift(v)
            if fs_norm:
                fs_key = fs_norm.replace("?", "")
                if fs_key in baseline_normalized:
                    matched_to = baseline_normalized[fs_key]
                    match_type = "frameshift_normalized"
                    results["stats"]["frameshift_matches"] += 1

        if not matched_to:
            # Nonsense
            ns_norm = normalize_nonsense(v)
            if ns_norm and ns_norm in baseline_normalized:
                matched_to = baseline_normalized[ns_norm]
                match_type = "nonsense_normalized"
                results["stats"]["nonsense_matches"] += 1

        if not matched_to:
            # Deletion
            del_norm = normalize_deletion(v)
            if del_norm:
                del_key = del_norm.replace("?", "")
                if del_key in baseline_normalized:
                    matched_to = baseline_normalized[del_key]
                    match_type = "deletion_normalized"
                    results["stats"]["deletion_matches"] += 1

        if not matched_to:
            # Duplication
            dup_norm = normalize_duplication(v)
            if dup_norm:
                dup_key = dup_norm.replace("?", "")
                if dup_key in baseline_normalized:
                    matched_to = baseline_normalized[dup_key]
                    match_type = "duplication_normalized"
                    results["stats"]["duplication_matches"] += 1

        # Try cDNA normalization
        if not matched_to:
            cdna = normalizer.normalize_cdna(v)
            if cdna and cdna in baseline_cdna:
                matched_to = cdna
                match_type = "cdna_prefix"
                results["stats"]["cdna_matches"] += 1

        # Try fuzzy position matching (±tolerance)
        if not matched_to and fuzzy_position:
            pos = normalizer.extract_position(v)
            single = normalizer.normalize_to_single_letter(v)
            var_type = get_variant_type(v)

            if pos and single:
                # Parse ref and alt from single
                m = re.match(r"^([A-Z])(\d+)(.+)$", single)
                if m:
                    ref, _, alt = m.groups()
                    for delta in range(-position_tolerance, position_tolerance + 1):
                        if delta == 0:
                            continue
                        test_pos = pos + delta
                        if test_pos in baseline_by_pos:
                            for bv in baseline_by_pos[test_pos]:
                                # Check type compatibility
                                b_type = get_variant_type(bv)
                                if var_type != b_type:
                                    continue

                                bsingle = normalizer.normalize_to_single_letter(bv)
                                if bsingle:
                                    bm = re.match(r"^([A-Z])(\d+)(.+)$", bsingle)
                                    if bm:
                                        b_ref, _, b_alt = bm.groups()
                                        if ref == b_ref and alt == b_alt:
                                            matched_to = bv
                                            match_type = f"fuzzy_pos_{delta:+d}"
                                            results["stats"]["fuzzy_matches"] += 1
                                            break
                        if matched_to:
                            break

        if matched_to:
            results["matches"].append(
                {"extracted": v, "matched_to": matched_to, "match_type": match_type}
            )
        else:
            results["unmatched"].append(v)
            results["stats"]["unmatched"] += 1

    return results


# =============================================================================
# STANDALONE FUNCTIONS (drop-in replacement for variant_utils.py)
# =============================================================================


def normalize_variant(variant: str, gene_symbol: str = "KCNH2") -> str:
    """
    Normalize any variant notation to a comparable format.

    This is the primary normalization function - use this for deduplication
    and comparison. Normalizes to single-letter format for protein variants
    (A561V) and c. prefix format for cDNA.

    DROP-IN REPLACEMENT for utils.variant_utils.normalize_variant()

    Normalization rules applied (in order):
    1. Strip whitespace and common prefixes
    2. Remove trailing asterisks from DMS studies: A561V* → A561V
    3. Remove trailing annotations: R176W(het) → R176W, A561V/WT → A561V
    4. Case normalization: a561v → A561V
    5. 3-letter → 1-letter amino acid: p.Ala561Val → A561V
    6. Frameshift standardization: fs*46 → fsX, fsXNN → fsX
    7. p. prefix removal: p.R176W → R176W

    Args:
        variant: Variant in any format (p.Ala561Val, A561V, c.1682C>T, etc.)
        gene_symbol: Gene symbol for context (default: KCNH2)

    Returns:
        Normalized variant string (uppercase, consistent format)

    Examples:
        >>> normalize_variant('p.Ala561Val')
        'A561V'
        >>> normalize_variant('A561V*')
        'A561V'
        >>> normalize_variant('R176W(het)')
        'R176W'
        >>> normalize_variant('a561v')
        'A561V'
        >>> normalize_variant('c.1682C>T')
        'c.1682C>T'
        >>> normalize_variant('IVS10+1G>A')
        'c.2398+1G>A'
    """
    if not variant:
        return ""

    variant = variant.strip()

    # === RULE 1: Remove trailing asterisks (DMS study artifacts) ===
    # A561V* → A561V, p.Arg231Cys* → p.Arg231Cys
    # But preserve * when it's a stop codon: R864* should stay R864* (later normalized to R864X)
    # Pattern: if * is at the very end AFTER a letter (not a number), remove it
    if variant.endswith("*"):
        # Check if this is likely a DMS asterisk vs stop codon
        # Stop codon: ends in digit + * (R864*, Arg864*)
        # DMS artifact: ends in letter + * (A561V*, Cys*)
        if len(variant) >= 2 and variant[-2].isalpha():
            variant = variant[:-1]

    # === RULE 2: Remove trailing annotations ===
    # R176W(het) → R176W, A561V/WT → A561V, G584S(hom) → G584S
    # Remove parenthetical annotations
    variant = re.sub(r"\([^)]*\)$", "", variant)
    # Remove /WT or /wt suffix
    variant = re.sub(r"/[Ww][Tt]$", "", variant)
    # Remove other slash annotations like /+ or /-
    variant = re.sub(r"/[\+\-]$", "", variant)

    variant = variant.strip()

    # === RULE 3: Early case normalization for simple variants ===
    # Handle lowercase variants like a561v → A561V before other processing
    # This catches variants that might not match patterns otherwise
    simple_lower_match = re.match(r"^([a-z])(\d+)([a-z\*])$", variant)
    if simple_lower_match:
        variant = f"{simple_lower_match.group(1).upper()}{simple_lower_match.group(2)}{simple_lower_match.group(3).upper()}"

    # Check comprehensive KCNH2 alias lookup first
    if gene_symbol.upper() == "KCNH2":
        # Try comprehensive aliases (from JSON) first
        canonical = _lookup_alias(variant, gene_symbol)
        if canonical:
            # Return canonical form (single-letter for protein, as-is for cDNA)
            if canonical.startswith("c."):
                return canonical
            return canonical

    # Create normalizer for the gene
    normalizer = VariantNormalizer(gene_symbol)

    # Try to detect variant type and normalize
    v_upper = variant.upper()
    v_lower = variant.lower()

    # Handle IVS notation (splice variants)
    ivs_match = re.match(
        r"^IVS(\d+)([\+\-]\d+)([ACGT])>([ACGT])$", variant, re.IGNORECASE
    )
    if ivs_match:
        ivs_num, offset, ref, alt = ivs_match.groups()
        # Try to convert to c. notation using gene-specific map
        if gene_symbol.upper() == "KCNH2":
            base_pos = KCNH2_IVS_MAP.get(f"IVS{ivs_num}")
            if base_pos:
                return f"{base_pos}{offset}{ref.upper()}>{alt.upper()}"
        # If no mapping, return cleaned up IVS notation
        return f"IVS{ivs_num}{offset}{ref.upper()}>{alt.upper()}"

    # === RULE 4: Protein variant detection ===
    # Try normalize_to_single_letter first - it handles most cases now

    # Check for protein-like patterns: p.XXX, Ala123Val, A123V, frameshifts
    is_protein_like = (
        v_lower.startswith("p.")
        or re.match(r"^[A-Z][a-z]{2}\d+", variant, re.IGNORECASE)
        or re.match(r"^[A-Za-z]\d+[A-Za-z\*X]$", variant, re.IGNORECASE)
        or re.match(r"^[A-Za-z]\d+fs", variant, re.IGNORECASE)
        or re.match(r"^[A-Za-z]\d+(stop|sp|ter)$", variant, re.IGNORECASE)
    )

    if is_protein_like:
        single = normalizer.normalize_to_single_letter(variant)
        if single:
            return single

    # cDNA variant detection
    if v_lower.startswith("c.") or re.match(r"^\d+[\+\-]?\d*[ACGT]", variant):
        cdna = normalizer.normalize_cdna(variant)
        if cdna:
            return cdna

    # Short protein format (A561V) - normalize case and stop codon
    # This is a fallback if normalize_to_single_letter didn't catch it
    m = re.match(r"^([A-Za-z])(\d+)([A-Za-z\*X])$", variant, re.IGNORECASE)
    if m:
        ref = m.group(1).upper()
        pos = m.group(2)
        alt = m.group(3).upper()
        if alt == "*":
            alt = "X"
        return f"{ref}{pos}{alt}"

    # Frameshift variants - ensure normalization
    if re.match(r"^[A-Za-z]\d+fs", variant, re.IGNORECASE):
        single = normalizer.normalize_to_single_letter(variant)
        if single:
            return single

    # Stop variants
    if re.match(r"^[A-Za-z]\d+(stop|sp|ter|\*|X)$", variant, re.IGNORECASE):
        single = normalizer.normalize_to_single_letter(variant)
        if single:
            return single

    # Deletion/insertion variants
    if re.match(r"^[A-Za-z]\d+(del|ins|dup)", variant, re.IGNORECASE):
        # Normalize to uppercase with standard suffix
        m = re.match(r"^([A-Za-z])(\d+)(del|ins|dup)(.*)$", variant, re.IGNORECASE)
        if m:
            return f"{m.group(1).upper()}{m.group(2)}{m.group(3).lower()}{m.group(4).upper()}"

    # Unknown format - return uppercase
    return variant.upper()


def normalize_protein_variant(variant: str) -> str:
    """
    Normalize a protein variant to short format (e.g., A561V).

    DROP-IN REPLACEMENT for utils.variant_utils.normalize_protein_variant()

    Args:
        variant: Protein variant in any format (p.Ala561Val, A561V, etc.)

    Returns:
        Normalized variant in short format (A561V)
    """
    if not variant:
        return ""

    normalizer = VariantNormalizer("UNKNOWN")
    single = normalizer.normalize_to_single_letter(variant)
    return single if single else variant.strip().upper()


def normalize_cdna_variant(variant: str) -> str:
    """
    Normalize a cDNA variant notation.

    DROP-IN REPLACEMENT for utils.variant_utils.normalize_cdna_variant()

    Args:
        variant: cDNA variant (c.1682C>T, 1682C>T, etc.)

    Returns:
        Normalized cDNA variant (c.1682C>T)
    """
    if not variant:
        return ""

    normalizer = VariantNormalizer("UNKNOWN")
    cdna = normalizer.normalize_cdna(variant)
    return cdna if cdna else variant.strip()


def variants_match(v1: str, v2: str, gene_symbol: str = "KCNH2") -> bool:
    """
    Check if two variant notations refer to the same variant.

    DROP-IN REPLACEMENT for utils.variant_utils.variants_match()

    Args:
        v1: First variant notation
        v2: Second variant notation
        gene_symbol: Gene symbol for context

    Returns:
        True if variants match after normalization
    """
    norm1 = normalize_variant(v1, gene_symbol)
    norm2 = normalize_variant(v2, gene_symbol)

    if not norm1 or not norm2:
        return False

    return norm1 == norm2


def find_matching_variants(
    extracted: List[str], expected: List[str], gene_symbol: str = "KCNH2"
) -> Tuple[List[Tuple[str, str]], List[str], List[str]]:
    """
    Compare extracted variants against expected variants.

    DROP-IN REPLACEMENT for utils.variant_utils.find_matching_variants()

    Args:
        extracted: List of extracted variant strings
        expected: List of expected variant strings
        gene_symbol: Gene symbol for context

    Returns:
        Tuple of (matched, missed, extra) where:
        - matched: List of (extracted, expected) tuples that match
        - missed: List of expected variants not found in extracted
        - extra: List of extracted variants not found in expected
    """
    # Normalize all variants
    norm_extracted = {normalize_variant(v, gene_symbol): v for v in extracted if v}
    norm_expected = {normalize_variant(v, gene_symbol): v for v in expected if v}

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


def create_variant_key(variant: Dict[str, Any], gene_symbol: str = "KCNH2") -> str:
    """
    Create a normalized key for variant grouping/deduplication.

    This is the key function for aggregation - ensures variants in different
    notations (p.Arg534Cys vs R534C) are grouped together.

    Args:
        variant: Variant dict with protein_notation, cdna_notation, etc.
        gene_symbol: Gene symbol for context

    Returns:
        Normalized variant key for grouping
    """
    # Try protein notation first (most common)
    protein = variant.get("protein_notation") or variant.get("protein_change")
    if protein:
        normalized = normalize_variant(protein, gene_symbol)
        if normalized:
            return f"{gene_symbol}:{normalized}"

    # Fall back to cDNA
    cdna = variant.get("cdna_notation") or variant.get("cdna_change")
    if cdna:
        normalized = normalize_variant(cdna, gene_symbol)
        if normalized:
            return f"{gene_symbol}:{normalized}"

    # Last resort: genomic position
    genomic = variant.get("genomic_position")
    if genomic:
        return f"{gene_symbol}:{genomic.strip()}"

    return f"{gene_symbol}:unknown_variant"


if __name__ == "__main__":
    # Test the normalizer
    print("=" * 60)
    print("Variant Normalizer Test Suite")
    print("=" * 60)

    norm = VariantNormalizer("KCNH2")

    # Test protein normalization
    test_variants = [
        "A561V",
        "p.Ala561Val",
        "Ala561Val",
        "p.A561V",
        "A193fsX",
        "G584S",
        "R864stop",
        "G184Del",
        "p.Arg534Cys",
        "R534C",  # Common variant in multiple notations
    ]

    print("\n1. Protein Normalization (Class-based):")
    for v in test_variants:
        forms = norm.get_all_forms(v)
        print(
            f"  {v:15} -> single: {forms.get('single', 'N/A'):12} three: {forms.get('three', 'N/A')}"
        )

    # Test standalone normalize_variant function
    print("\n2. Standalone normalize_variant() Function:")
    standalone_tests = [
        "p.Ala561Val",  # Three-letter protein
        "A561V",  # Already normalized
        "c.1682C>T",  # cDNA
        "1682C>T",  # cDNA without prefix
        "IVS10+1G>A",  # Splice variant (IVS notation)
        "A193fsX10",  # Frameshift with position
        "R864*",  # Stop codon
        "p.Gly628Ser",  # KCNH2 variant
    ]
    for v in standalone_tests:
        normalized = normalize_variant(v)
        print(f"  {v:18} -> {normalized}")

    # === NEW TESTS: Normalization gap fixes (2026-02-10) ===
    print("\n8. NEW Normalization Gap Fixes (2026-02-10):")

    # Test 8a: Trailing asterisk removal (DMS artifacts)
    print("\n  8a. Trailing Asterisk Removal (DMS artifacts):")
    asterisk_tests = [
        ("A561V*", "A561V"),  # DMS artifact
        ("p.Arg231Cys*", "R231C"),  # Three-letter with asterisk
        ("M291T*", "M291T"),  # DMS study variant
        ("R864*", "R864X"),  # Stop codon - should convert * to X
    ]
    for v, expected in asterisk_tests:
        result = normalize_variant(v)
        status = "✓" if result == expected else f"✗ (got {result})"
        print(f"    {v:18} -> {result:10} expected {expected:10} {status}")

    # Test 8b: Trailing annotation removal
    print("\n  8b. Trailing Annotation Removal:")
    annotation_tests = [
        ("R176W(het)", "R176W"),
        ("A561V/WT", "A561V"),
        ("G584S(hom)", "G584S"),
        ("R534C/+", "R534C"),
    ]
    for v, expected in annotation_tests:
        result = normalize_variant(v)
        status = "✓" if result == expected else f"✗ (got {result})"
        print(f"    {v:18} -> {result:10} expected {expected:10} {status}")

    # Test 8c: Case normalization
    print("\n  8c. Case Normalization:")
    case_tests = [
        ("a561v", "A561V"),
        ("r176w", "R176W"),
        ("g584s", "G584S"),
    ]
    for v, expected in case_tests:
        result = normalize_variant(v)
        status = "✓" if result == expected else f"✗ (got {result})"
        print(f"    {v:18} -> {result:10} expected {expected:10} {status}")

    # Test 8d: Frameshift standardization
    print("\n  8d. Frameshift Standardization:")
    fs_tests = [
        ("A193fs*46", "A193fsX"),
        ("p.Gly262Alafs*98", "G262fsX"),
        ("G262fsX10", "G262fsX"),
        ("A1077fs+X*", "A1077fsX"),
        ("p.Arg518fs*", "R518fsX"),
        ("Gly24fsTer58", "G24fsX"),
    ]
    for v, expected in fs_tests:
        result = normalize_variant(v)
        status = "✓" if result == expected else f"✗ (got {result})"
        print(f"    {v:18} -> {result:10} expected {expected:10} {status}")

    # Test 8e: p. prefix removal (should already work)
    print("\n  8e. p. Prefix Removal:")
    prefix_tests = [
        ("p.R176W", "R176W"),
        ("p.A561V", "A561V"),
        ("p.Ala561Val", "A561V"),
    ]
    for v, expected in prefix_tests:
        result = normalize_variant(v)
        status = "✓" if result == expected else f"✗ (got {result})"
        print(f"    {v:18} -> {result:10} expected {expected:10} {status}")

    # Test non-target detection
    non_target_tests = [
        "R248W",  # TP53 hotspot
        "G12D",  # KRAS hotspot
        "V600E",  # BRAF hotspot
        "P2006A",  # Position > KCNH2 length
        "A561V",  # Valid KCNH2 variant
    ]

    print("\n3. Non-Target Variant Detection:")
    for v in non_target_tests:
        is_non, reason = norm.is_non_target_variant(v)
        status = f"⚠️  NON-TARGET: {reason}" if is_non else "✓ Valid"
        print(f"  {v:12} -> {status}")

    # Test cDNA normalization
    cdna_tests = [
        "c.1234A>G",
        "1234A>G",  # Missing c. prefix
        "c.3152+1G>A",  # Intronic
        "3152+1G>A",  # Missing prefix, intronic
    ]

    print("\n4. cDNA Normalization:")
    for v in cdna_tests:
        normalized = norm.normalize_cdna(v)
        print(f"  {v:18} -> {normalized}")

    # Test variant matching
    print("\n5. variants_match() Function:")
    match_tests = [
        ("p.Ala561Val", "A561V"),
        ("p.Arg534Cys", "R534C"),
        ("c.1682C>T", "1682C>T"),
        ("A561V", "G584S"),  # Should NOT match
    ]
    for v1, v2 in match_tests:
        result = variants_match(v1, v2)
        status = "✓ MATCH" if result else "✗ no match"
        print(f"  {v1:15} vs {v2:10} -> {status}")

    # Test create_variant_key for aggregation
    print("\n6. create_variant_key() for Aggregation:")
    key_tests = [
        {"protein_notation": "p.Ala561Val", "cdna_notation": "c.1682C>T"},
        {"protein_notation": "A561V"},
        {"protein_notation": "p.Arg534Cys"},
        {"cdna_notation": "c.1600C>T"},
    ]
    for vdict in key_tests:
        key = create_variant_key(vdict)
        print(f"  {vdict} -> {key}")

    # Test baseline matching
    print("\n7. Baseline Matching with Fuzzy Position:")
    baseline = {"A561V", "G584S", "R534C", "c.1234A>G"}
    extracted = ["p.Ala561Val", "G585S", "R248W", "1234A>G", "unknown"]

    results = match_variants_to_baseline(
        extracted, baseline, "KCNH2", fuzzy_position=True
    )

    print(f"  Matches: {len(results['matches'])}")
    for m in results["matches"]:
        print(f"    {m['extracted']} -> {m['matched_to']} ({m['match_type']})")

    print(f"  Filtered (non-target): {len(results['filtered_non_target'])}")
    for v, reason in results["filtered_non_target"]:
        print(f"    {v}: {reason}")

    print(f"  Unmatched: {results['unmatched']}")
    print(f"  Stats: {results['stats']}")

    print("\n" + "=" * 60)
    print("All tests complete!")
