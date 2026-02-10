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
AA_MAP_REVERSE['Xaa'] = 'X'  # Unknown amino acid

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

# Non-target gene hotspots (TP53, KRAS, BRAF, PIK3CA)
# These are commonly extracted by mistake when searching for cardiac gene variants
NON_TARGET_HOTSPOTS = {
    # TP53 (most common cancer gene mutations)
    'R175H', 'R248W', 'R248Q', 'R249S', 'R273H', 'R273C', 'R282W',
    'G245S', 'G245D', 'Y220C', 'C176F', 'C176Y', 'C242F', 'C242S',
    'H179R', 'H179Y', 'M237I', 'S241F', 'S241C', 'C277F', 'Y163C',
    # KRAS
    'G12D', 'G12V', 'G12C', 'G12R', 'G12A', 'G12S', 'G13D',
    'Q61H', 'Q61L', 'Q61R', 'Q61G',
    # BRAF
    'V600E', 'V600K',
    # PIK3CA
    'E545K', 'H1047R', 'H1047L',
}

# KCNH2 variant aliases from ClinVar - maps alternative names to canonical form
# Key: canonical single-letter form, Values: alternative representations
KCNH2_VARIANT_ALIASES = {
    # Well-known founder mutations with multiple names
    'A561V': ['A561V', 'p.Ala561Val', 'Ala561Val', 'c.1682C>T'],
    'R534C': ['R534C', 'p.Arg534Cys', 'Arg534Cys', 'c.1600C>T'],
    'A614V': ['A614V', 'p.Ala614Val', 'Ala614Val', 'c.1841C>T'],
    'G628S': ['G628S', 'p.Gly628Ser', 'Gly628Ser', 'c.1882G>A'],
    'N470D': ['N470D', 'p.Asn470Asp', 'Asn470Asp', 'c.1408A>G'],
    'R176W': ['R176W', 'p.Arg176Trp', 'Arg176Trp', 'c.526C>T'],
    'T613M': ['T613M', 'p.Thr613Met', 'Thr613Met', 'c.1838C>T'],
    'A422T': ['A422T', 'p.Ala422Thr', 'Ala422Thr', 'c.1264G>A'],
    'R784W': ['R784W', 'p.Arg784Trp', 'Arg784Trp', 'c.2350C>T'],
    'Y611H': ['Y611H', 'p.Tyr611His', 'Tyr611His', 'c.1831T>C'],
    # Common splice variants
    'c.2398+1G>A': ['c.2398+1G>A', 'IVS10+1G>A'],
    'c.1129-1G>A': ['c.1129-1G>A', 'IVS5-1G>A'],
}

# Build reverse lookup for aliases
_KCNH2_ALIAS_LOOKUP = {}
for canonical, aliases in KCNH2_VARIANT_ALIASES.items():
    for alias in aliases:
        _KCNH2_ALIAS_LOOKUP[alias.upper()] = canonical

# IVS to cDNA splice notation mapping for KCNH2
# This is a partial list - in practice would need full gene coordinates
KCNH2_IVS_MAP = {
    'IVS1': 'c.175',
    'IVS2': 'c.453', 
    'IVS3': 'c.575',
    'IVS4': 'c.787',
    'IVS5': 'c.1129',
    'IVS6': 'c.1351',
    'IVS7': 'c.1592',
    'IVS8': 'c.1759',
    'IVS9': 'c.2006',
    'IVS10': 'c.2398',
    'IVS11': 'c.2599',
    'IVS12': 'c.2769',
    'IVS13': 'c.3040',
}

# Patterns for variant parsing
PROTEIN_PATTERNS = [
    # Already normalized three-letter: p.Ala561Val
    re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|[fs\*]|[delins\*])$', re.IGNORECASE),
    
    # Three-letter without prefix: Ala561Val
    re.compile(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', re.IGNORECASE),
    
    # Single-letter conversion: A561V -> p.Ala561Val
    re.compile(r'^(?:p\.)?([A-Z])(\d+)([A-Z])$', re.IGNORECASE),
    
    # Frame-shift with single-letter: A193fsX, A193fsX10, A193fs*10
    re.compile(r'^(?:p\.)?([A-Z])(\d+)(fsX\d*|fs\*\d*|fs)$', re.IGNORECASE),
    
    # Frame-shift with three-letter: Ala193fsX, Ala193fs*10
    re.compile(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)(fsX\d*|fs\*\d*|fs)$', re.IGNORECASE),
    
    # Insertion/deletion: G184Del, G184Ins, G184_G185del
    re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(Del|Ins|Dup)$', re.IGNORECASE),
    
    # Stop/truncation: R864stop, D864sp, R864*, R864X, R864Ter
    re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(sp|stop|trunc|\*|Ter|X)$', re.IGNORECASE),
    
    # Extended deletion/dup with range: p.Gly184_Leu186del
    re.compile(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)(del|dup|ins)$', re.IGNORECASE),
]

CDNA_PATTERNS = [
    # c.1234A>G
    re.compile(r'^c\.(\d+)([ACGT])>([ACGT])$'),
    # c.1234delA, c.1234dupG
    re.compile(r'^c\.(\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
    # c.1234_1235delAG
    re.compile(r'^c\.(\d+)_(\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
    # c.1234+1G>A (intronic - donor)
    re.compile(r'^c\.(\d+[\+]\d+)([ACGT])>([ACGT])$'),
    # c.1234-1G>A (intronic - acceptor)
    re.compile(r'^c\.(\d+[\-]\d+)([ACGT])>([ACGT])$'),
    # c.1234+1del or c.1234-1del (intronic indels)
    re.compile(r'^c\.(\d+[\+\-]\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
]

# IVS (splice) patterns - legacy notation
SPLICE_PATTERNS = [
    # IVS10+1G>A, IVS5-2A>G
    re.compile(r'^IVS(\d+)([\+\-]\d+)([ACGT])>([ACGT])$', re.IGNORECASE),
    # IVS10+1del
    re.compile(r'^IVS(\d+)([\+\-]\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
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
    
    def normalize_to_single_letter(self, variant: str) -> Optional[str]:
        """
        Normalize protein variant to single-letter format: A561V
        
        Handles: A561V, p.Ala561Val, Ala561Val, p.A561V
        
        Args:
            variant: Variant string in various formats
            
        Returns:
            Single-letter format (e.g., 'A561V') or None if unparseable
        """
        if not variant:
            return None
        
        v = variant.strip()
        if v.startswith('p.'):
            v = v[2:]
        
        # Already single-letter
        m = re.match(r'^([A-Z])(\d+)([A-Z\*X])$', v)
        if m:
            return f"{m.group(1)}{m.group(2)}{m.group(3)}"
        
        # Three-letter: Ala561Val
        m = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Ter|X)$', v)
        if m:
            ref = AA_MAP_REVERSE.get(m.group(1).capitalize())
            alt_raw = m.group(3)
            if alt_raw in ['*', 'Ter', 'X']:
                alt = 'X'
            else:
                alt = AA_MAP_REVERSE.get(alt_raw.capitalize())
            if ref and alt:
                return f"{ref}{m.group(2)}{alt}"
        
        # Frameshift: A193fsX
        m = re.match(r'^([A-Z]|[A-Z][a-z]{2})(\d+)(fs[X\*]?|fsX\d*)$', v, re.IGNORECASE)
        if m:
            ref_raw = m.group(1)
            if len(ref_raw) == 1:
                ref = ref_raw.upper()
            else:
                ref = AA_MAP_REVERSE.get(ref_raw.capitalize())
            if ref:
                return f"{ref}{m.group(2)}fsX"
        
        # Stop variants
        m = re.match(r'^([A-Z])(\d+)(sp|stop|X|\*)$', v, re.IGNORECASE)
        if m:
            return f"{m.group(1)}{m.group(2)}X"
        
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
        if variant.startswith('c.'):
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
        if re.match(r'^\d+[\+\-\d]*[ACGT_>delinsdup]', variant):
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
            m = re.match(r'^[A-Z](\d+)', single)
            if m:
                return int(m.group(1))
        
        # Fallback: try to extract from three-letter or p. notation
        m = re.search(r'(\d+)', variant)
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
        check_var = re.sub(r'fsX$', '', single)
        if check_var in NON_TARGET_HOTSPOTS:
            return True, f"Known hotspot in TP53/KRAS/BRAF/PIK3CA"
        
        # Check position against protein length
        m = re.match(r'^[A-Z](\d+)', single)
        if m:
            pos = int(m.group(1))
            if self.protein_length and pos > self.protein_length:
                return True, f"Position {pos} exceeds {self.gene_symbol} length ({self.protein_length})"
        
        return False, None
    
    def get_all_forms(self, variant: str) -> Dict[str, str]:
        """
        Get all normalized forms of a variant.
        
        Args:
            variant: Variant string
            
        Returns:
            Dict with keys: original, single, three, cdna (as available)
        """
        result = {'original': variant}
        
        single = self.normalize_to_single_letter(variant)
        if single:
            result['single'] = single
        
        three = self.normalize_protein(variant)
        if three:
            result['three'] = three
        
        cdna = self.normalize_cdna(variant)
        if cdna:
            result['cdna'] = cdna
        
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
                
                # Check if non-target
                is_non_target, reason = self.is_non_target_variant(protein)
                if is_non_target:
                    result['is_non_target'] = True
                    result['non_target_reason'] = reason
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
        'non_target_filtered': 0,
    }
    
    for v in variants:
        result = normalizer.get_canonical_form(v)
        key = result.get('canonical_key', '')
        
        if key in seen_keys:
            stats['duplicates_removed'] += 1
            continue
        
        # Track non-target variants
        if result.get('is_non_target'):
            stats['non_target_filtered'] += 1
            
        seen_keys.add(key)
        normalized.append(result)
        
        if result.get('protein_notation_normalized') or result.get('cdna_notation_normalized'):
            stats['normalized'] += 1
        else:
            stats['normalization_failed'] += 1
    
    stats['total_output'] = len(normalized)
    
    return normalized, stats


def match_variants_to_baseline(
    extracted: List[str],
    baseline: Set[str],
    gene_symbol: str = 'KCNH2',
    fuzzy_position: bool = True
) -> Dict[str, Any]:
    """
    Match extracted variants to a baseline set with improved normalization.
    
    Supports exact matching, normalized matching, cDNA prefix matching,
    and optional fuzzy position matching (±1 for off-by-one errors).
    
    Args:
        extracted: List of extracted variant strings
        baseline: Set of baseline variant strings
        gene_symbol: Target gene symbol
        fuzzy_position: Allow ±1 position matching
    
    Returns:
        Dict with matches, unmatched, filtered_non_target, and stats
    """
    normalizer = VariantNormalizer(gene_symbol)
    
    # Build baseline position index for fuzzy matching
    baseline_by_pos = {}
    baseline_singles = set()
    for v in baseline:
        single = normalizer.normalize_to_single_letter(v)
        if single:
            baseline_singles.add(single)
            pos = normalizer.extract_position(v)
            if pos:
                if pos not in baseline_by_pos:
                    baseline_by_pos[pos] = []
                baseline_by_pos[pos].append(v)
    
    # Also index cDNA forms
    baseline_cdna = set()
    for v in baseline:
        if v.startswith('c.') or v.startswith('IVS'):
            baseline_cdna.add(v)
    
    results = {
        'matches': [],
        'unmatched': [],
        'filtered_non_target': [],
        'stats': {
            'total_input': len(extracted),
            'exact_matches': 0,
            'normalized_matches': 0,
            'fuzzy_matches': 0,
            'cdna_matches': 0,
            'filtered': 0,
            'unmatched': 0,
        }
    }
    
    for v in extracted:
        # Check if non-target gene variant
        is_non_target, reason = normalizer.is_non_target_variant(v)
        if is_non_target:
            results['filtered_non_target'].append((v, reason))
            results['stats']['filtered'] += 1
            continue
        
        matched_to = None
        match_type = None
        
        # Try exact match
        if v in baseline:
            matched_to = v
            match_type = 'exact'
            results['stats']['exact_matches'] += 1
        
        # Try normalized single-letter match
        if not matched_to:
            single = normalizer.normalize_to_single_letter(v)
            if single and single in baseline_singles:
                # Find the original baseline form
                for bv in baseline:
                    if normalizer.normalize_to_single_letter(bv) == single:
                        matched_to = bv
                        match_type = 'normalized'
                        results['stats']['normalized_matches'] += 1
                        break
        
        # Try cDNA normalization
        if not matched_to:
            cdna = normalizer.normalize_cdna(v)
            if cdna and cdna in baseline_cdna:
                matched_to = cdna
                match_type = 'cdna_prefix'
                results['stats']['cdna_matches'] += 1
        
        # Try fuzzy position matching (±1)
        if not matched_to and fuzzy_position:
            pos = normalizer.extract_position(v)
            single = normalizer.normalize_to_single_letter(v)
            if pos and single:
                # Parse ref and alt from single
                m = re.match(r'^([A-Z])(\d+)(.+)$', single)
                if m:
                    ref, _, alt = m.groups()
                    for delta in [-1, 1]:
                        test_pos = pos + delta
                        if test_pos in baseline_by_pos:
                            for bv in baseline_by_pos[test_pos]:
                                bsingle = normalizer.normalize_to_single_letter(bv)
                                if bsingle:
                                    bm = re.match(r'^([A-Z])(\d+)(.+)$', bsingle)
                                    if bm:
                                        b_ref, _, b_alt = bm.groups()
                                        if ref == b_ref and alt == b_alt:
                                            matched_to = bv
                                            match_type = f'fuzzy_pos_{delta:+d}'
                                            results['stats']['fuzzy_matches'] += 1
                                            break
                        if matched_to:
                            break
        
        if matched_to:
            results['matches'].append({
                'extracted': v,
                'matched_to': matched_to,
                'match_type': match_type
            })
        else:
            results['unmatched'].append(v)
            results['stats']['unmatched'] += 1
    
    return results


# =============================================================================
# STANDALONE FUNCTIONS (drop-in replacement for variant_utils.py)
# =============================================================================

def normalize_variant(variant: str, gene_symbol: str = 'KCNH2') -> str:
    """
    Normalize any variant notation to a comparable format.
    
    This is the primary normalization function - use this for deduplication
    and comparison. Normalizes to single-letter format for protein variants
    (A561V) and c. prefix format for cDNA.
    
    DROP-IN REPLACEMENT for utils.variant_utils.normalize_variant()
    
    Args:
        variant: Variant in any format (p.Ala561Val, A561V, c.1682C>T, etc.)
        gene_symbol: Gene symbol for context (default: KCNH2)
        
    Returns:
        Normalized variant string (uppercase, consistent format)
        
    Examples:
        >>> normalize_variant('p.Ala561Val')
        'A561V'
        >>> normalize_variant('c.1682C>T')
        'c.1682C>T'
        >>> normalize_variant('IVS10+1G>A')
        'c.2398+1G>A'
    """
    if not variant:
        return ""
    
    variant = variant.strip()
    
    # Check KCNH2 alias lookup first
    if gene_symbol.upper() == 'KCNH2':
        canonical = _KCNH2_ALIAS_LOOKUP.get(variant.upper())
        if canonical:
            # Return canonical form (single-letter for protein, as-is for cDNA)
            if canonical.startswith('c.'):
                return canonical
            return canonical
    
    # Create normalizer for the gene
    normalizer = VariantNormalizer(gene_symbol)
    
    # Try to detect variant type and normalize
    v_upper = variant.upper()
    v_lower = variant.lower()
    
    # Handle IVS notation (splice variants)
    ivs_match = re.match(r'^IVS(\d+)([\+\-]\d+)([ACGT])>([ACGT])$', variant, re.IGNORECASE)
    if ivs_match:
        ivs_num, offset, ref, alt = ivs_match.groups()
        # Try to convert to c. notation using gene-specific map
        if gene_symbol.upper() == 'KCNH2':
            base_pos = KCNH2_IVS_MAP.get(f'IVS{ivs_num}')
            if base_pos:
                return f"{base_pos}{offset}{ref.upper()}>{alt.upper()}"
        # If no mapping, return cleaned up IVS notation
        return f"IVS{ivs_num}{offset}{ref.upper()}>{alt.upper()}"
    
    # Protein variant detection
    if v_lower.startswith('p.') or re.match(r'^[A-Z][a-z]{2}\d+', variant):
        # Protein variant with p. prefix or 3-letter AA
        single = normalizer.normalize_to_single_letter(variant)
        if single:
            return single
    
    # cDNA variant detection  
    if v_lower.startswith('c.') or re.match(r'^\d+[\+\-]?\d*[ACGT]', variant):
        cdna = normalizer.normalize_cdna(variant)
        if cdna:
            return cdna
    
    # Short protein format (A561V) - already normalized
    if re.match(r'^[A-Z]\d+[A-Z*X]$', v_upper):
        return v_upper.replace('*', 'X')
    
    # Frameshift variants
    if re.match(r'^[A-Z]\d+fs', variant, re.IGNORECASE):
        single = normalizer.normalize_to_single_letter(variant)
        if single:
            return single
    
    # Stop variants
    if re.match(r'^[A-Z]\d+(stop|sp|ter|\*|X)$', variant, re.IGNORECASE):
        single = normalizer.normalize_to_single_letter(variant)
        if single:
            return single
    
    # Deletion/insertion variants
    if re.match(r'^[A-Z]\d+(del|ins|dup)', variant, re.IGNORECASE):
        # Normalize to uppercase with standard suffix
        m = re.match(r'^([A-Z])(\d+)(del|ins|dup)(.*)$', variant, re.IGNORECASE)
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
    
    normalizer = VariantNormalizer('UNKNOWN')
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
    
    normalizer = VariantNormalizer('UNKNOWN')
    cdna = normalizer.normalize_cdna(variant)
    return cdna if cdna else variant.strip()


def variants_match(v1: str, v2: str, gene_symbol: str = 'KCNH2') -> bool:
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
    extracted: List[str], 
    expected: List[str],
    gene_symbol: str = 'KCNH2'
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


def create_variant_key(variant: Dict[str, Any], gene_symbol: str = 'KCNH2') -> str:
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
    protein = variant.get('protein_notation') or variant.get('protein_change')
    if protein:
        normalized = normalize_variant(protein, gene_symbol)
        if normalized:
            return f"{gene_symbol}:{normalized}"
    
    # Fall back to cDNA
    cdna = variant.get('cdna_notation') or variant.get('cdna_change')
    if cdna:
        normalized = normalize_variant(cdna, gene_symbol)
        if normalized:
            return f"{gene_symbol}:{normalized}"
    
    # Last resort: genomic position
    genomic = variant.get('genomic_position')
    if genomic:
        return f"{gene_symbol}:{genomic.strip()}"
    
    return f"{gene_symbol}:unknown_variant"


if __name__ == '__main__':
    # Test the normalizer
    print("=" * 60)
    print("Variant Normalizer Test Suite")
    print("=" * 60)
    
    norm = VariantNormalizer('KCNH2')
    
    # Test protein normalization
    test_variants = [
        'A561V', 'p.Ala561Val', 'Ala561Val', 'p.A561V',
        'A193fsX', 'G584S', 'R864stop', 'G184Del',
        'p.Arg534Cys', 'R534C',  # Common variant in multiple notations
    ]
    
    print("\n1. Protein Normalization (Class-based):")
    for v in test_variants:
        forms = norm.get_all_forms(v)
        print(f"  {v:15} -> single: {forms.get('single', 'N/A'):12} three: {forms.get('three', 'N/A')}")
    
    # Test standalone normalize_variant function
    print("\n2. Standalone normalize_variant() Function:")
    standalone_tests = [
        'p.Ala561Val',   # Three-letter protein
        'A561V',         # Already normalized
        'c.1682C>T',     # cDNA
        '1682C>T',       # cDNA without prefix
        'IVS10+1G>A',    # Splice variant (IVS notation)
        'A193fsX10',     # Frameshift with position
        'R864*',         # Stop codon
        'p.Gly628Ser',   # KCNH2 variant
    ]
    for v in standalone_tests:
        normalized = normalize_variant(v)
        print(f"  {v:18} -> {normalized}")
    
    # Test non-target detection
    non_target_tests = [
        'R248W',    # TP53 hotspot
        'G12D',     # KRAS hotspot
        'V600E',    # BRAF hotspot
        'P2006A',   # Position > KCNH2 length
        'A561V',    # Valid KCNH2 variant
    ]
    
    print("\n3. Non-Target Variant Detection:")
    for v in non_target_tests:
        is_non, reason = norm.is_non_target_variant(v)
        status = f"⚠️  NON-TARGET: {reason}" if is_non else "✓ Valid"
        print(f"  {v:12} -> {status}")
    
    # Test cDNA normalization  
    cdna_tests = [
        'c.1234A>G',
        '1234A>G',       # Missing c. prefix
        'c.3152+1G>A',   # Intronic
        '3152+1G>A',     # Missing prefix, intronic
    ]
    
    print("\n4. cDNA Normalization:")
    for v in cdna_tests:
        normalized = norm.normalize_cdna(v)
        print(f"  {v:18} -> {normalized}")
    
    # Test variant matching
    print("\n5. variants_match() Function:")
    match_tests = [
        ('p.Ala561Val', 'A561V'),
        ('p.Arg534Cys', 'R534C'),
        ('c.1682C>T', '1682C>T'),
        ('A561V', 'G584S'),  # Should NOT match
    ]
    for v1, v2 in match_tests:
        result = variants_match(v1, v2)
        status = "✓ MATCH" if result else "✗ no match"
        print(f"  {v1:15} vs {v2:10} -> {status}")
    
    # Test create_variant_key for aggregation
    print("\n6. create_variant_key() for Aggregation:")
    key_tests = [
        {'protein_notation': 'p.Ala561Val', 'cdna_notation': 'c.1682C>T'},
        {'protein_notation': 'A561V'},
        {'protein_notation': 'p.Arg534Cys'},
        {'cdna_notation': 'c.1600C>T'},
    ]
    for vdict in key_tests:
        key = create_variant_key(vdict)
        print(f"  {vdict} -> {key}")
    
    # Test baseline matching
    print("\n7. Baseline Matching with Fuzzy Position:")
    baseline = {'A561V', 'G584S', 'R534C', 'c.1234A>G'}
    extracted = ['p.Ala561Val', 'G585S', 'R248W', '1234A>G', 'unknown']
    
    results = match_variants_to_baseline(extracted, baseline, 'KCNH2', fuzzy_position=True)
    
    print(f"  Matches: {len(results['matches'])}")
    for m in results['matches']:
        print(f"    {m['extracted']} -> {m['matched_to']} ({m['match_type']})")
    
    print(f"  Filtered (non-target): {len(results['filtered_non_target'])}")
    for v, reason in results['filtered_non_target']:
        print(f"    {v}: {reason}")
    
    print(f"  Unmatched: {results['unmatched']}")
    print(f"  Stats: {results['stats']}")
    
    print("\n" + "=" * 60)
    print("All tests complete!")
