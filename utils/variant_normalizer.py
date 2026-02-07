"""
Variant Normalizer for GeneVariantFetcher

Normalizes variant nomenclature to standard HGVS-style forms and validates
positions against known gene/protein lengths.
"""

import re
import logging
from typing import Optional, Dict, Any, Tuple

logger = logging.getLogger(__name__)

# Amino acid single-letter to three-letter code mapping
AA_MAP = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Xaa'
}

# Reverse mapping (three-letter to single-letter)
AA_MAP_REVERSE = {v: k for k, v in AA_MAP.items()}
AA_MAP_REVERSE['Ter'] = '*'
AA_MAP_REVERSE['Stop'] = '*'

# Protein lengths for common cardiac genes
PROTEIN_LENGTHS = {
    'KCNH2': 1159,
    'KCNQ1': 676,
    'SCN5A': 2016,
    'KCNE1': 129,
    'KCNE2': 123,
    'KCNJ2': 427,
    'CACNA1C': 2221,
    'RYR2': 4967,
}

# Patterns for variant parsing
PROTEIN_PATTERNS = [
    # Already normalized three-letter: p.Ala561Val
    re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|[fs\*]|[delins\*])$', re.IGNORECASE),
    
    # Three-letter without prefix: Ala561Val
    re.compile(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', re.IGNORECASE),
    
    # Single-letter conversion: A561V -> p.Ala561Val
    re.compile(r'^(?:p\.)?([A-Z])(\d+)([A-Z])$', re.IGNORECASE),
    
    # Frame-shift with single-letter: A193fsX
    re.compile(r'^(?:p\.)?([A-Z])(\d+)(fsX|fs\*|fs)$', re.IGNORECASE),
    
    # Frame-shift with three-letter: Ala193fsX
    re.compile(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)(fsX|fs\*|fs)$', re.IGNORECASE),
    
    # Insertion/deletion: G184Del, G184Ins
    re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(Del|Ins)$', re.IGNORECASE),
    
    # Stop/truncation: R864stop, D864sp
    re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(sp|stop|trunc|\*|Ter)$', re.IGNORECASE),
]

CDNA_PATTERNS = [
    # c.1234A>G
    re.compile(r'^c\.(\d+)([ACGT])>([ACGT])$'),
    # c.1234delA, c.1234dupG
    re.compile(r'^c\.(\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
    # c.1234_1235delAG
    re.compile(r'^c\.(\d+)_(\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
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
                if suffix_upper.startswith('FS'):
                    # Frame-shift mutations
                    normalized_suffix = 'fs*'
                elif suffix_upper in ['X', 'TER']:
                    # Stop codon
                    normalized_suffix = '*'
                elif suffix_upper in ['SP', 'STOP', 'TRUNC']:
                    # Alternative stop annotations
                    normalized_suffix = '*'
                elif suffix_upper in ['DEL', 'INS']:
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
                elif '*' in suffix:
                    # Stop codon variants
                    normalized_suffix = '*'
                else:
                    continue
                
                return f"p.{ref_3}{pos}{normalized_suffix}"
        
        # If no pattern matched, return None
        logger.debug(f"Could not normalize protein variant: {original_variant}")
        return None
    
    def normalize_cdna(self, variant: str) -> Optional[str]:
        """
        Normalize cDNA variant to c.###X>Y format.
        
        Args:
            variant: Variant string
            
        Returns:
            Normalized variant string or None if unparseable
        """
        if not variant:
            return None
            
        variant = variant.strip()
        
        # Already in good format
        for pattern in CDNA_PATTERNS:
            if pattern.match(variant):
                return variant
        
        # Try adding c. prefix if missing
        if not variant.startswith('c.'):
            prefixed = f"c.{variant}"
            for pattern in CDNA_PATTERNS:
                if pattern.match(prefixed):
                    return prefixed
        
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
        protein = variant.get('protein_notation') or variant.get('protein_change')
        if protein:
            normalized = self.normalize_protein(protein)
            if normalized:
                result['protein_notation'] = normalized
                result['protein_notation_normalized'] = True
                
                # Extract position and validate
                match = re.search(r'(\d+)', normalized)
                if match:
                    pos = int(match.group(1))
                    result['position'] = pos
                    result['position_valid'] = self.validate_position(pos)
            else:
                result['protein_notation_normalized'] = False
        
        # Normalize cDNA notation
        cdna = variant.get('cdna_notation') or variant.get('cdna_change')
        if cdna:
            normalized = self.normalize_cdna(cdna)
            if normalized:
                result['cdna_notation'] = normalized
                result['cdna_notation_normalized'] = True
            else:
                result['cdna_notation_normalized'] = False
        
        # Create canonical key for deduplication
        result['canonical_key'] = self._create_canonical_key(result)
        
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
        protein = variant.get('protein_notation', '')
        cdna = variant.get('cdna_notation', '')
        
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
        'total_input': len(variants),
        'normalized': 0,
        'duplicates_removed': 0,
        'normalization_failed': 0,
    }
    
    for v in variants:
        result = normalizer.get_canonical_form(v)
        key = result.get('canonical_key', '')
        
        if key in seen_keys:
            stats['duplicates_removed'] += 1
            continue
            
        seen_keys.add(key)
        normalized.append(result)
        
        if result.get('protein_notation_normalized') or result.get('cdna_notation_normalized'):
            stats['normalized'] += 1
        else:
            stats['normalization_failed'] += 1
    
    stats['total_output'] = len(normalized)
    
    return normalized, stats
