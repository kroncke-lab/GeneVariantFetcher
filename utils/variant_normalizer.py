"""
Variant Normalizer for GeneVariantFetcher

Normalizes variant nomenclature to standard HGVS-style forms and validates
positions against known gene/protein lengths.

Merged features from improved_variant_normalizer.py:
- Non-target variant detection (TP53/KRAS/BRAF/PIK3CA hotspots)
- Fuzzy position matching for off-by-one errors  
- Enhanced cDNA normalization with prefix handling
- normalize_to_single_letter() and get_all_forms() methods
- match_variants_to_baseline() function
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
    # c.1234+1G>A (intronic)
    re.compile(r'^c\.(\d+[\+\-]\d+)([ACGT])>([ACGT])$'),
    # c.1234-1G>A (intronic)
    re.compile(r'^c\.(\d+[\+\-]\d+)(del|dup|ins)([ACGT]*)$', re.IGNORECASE),
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


if __name__ == '__main__':
    # Test the normalizer
    print("=" * 60)
    print("Variant Normalizer Test Suite")
    print("=" * 60)
    
    norm = VariantNormalizer('KCNH2')
    
    # Test protein normalization
    test_variants = [
        'A561V', 'p.Ala561Val', 'Ala561Val', 'p.A561V',
        'A193fsX', 'G584S', 'R864stop', 'G184Del'
    ]
    
    print("\n1. Protein Normalization:")
    for v in test_variants:
        forms = norm.get_all_forms(v)
        print(f"  {v:15} -> single: {forms.get('single', 'N/A'):12} three: {forms.get('three', 'N/A')}")
    
    # Test non-target detection
    non_target_tests = [
        'R248W',    # TP53 hotspot
        'G12D',     # KRAS hotspot
        'V600E',    # BRAF hotspot
        'P2006A',   # Position > KCNH2 length
        'A561V',    # Valid KCNH2 variant
    ]
    
    print("\n2. Non-Target Variant Detection:")
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
    
    print("\n3. cDNA Normalization:")
    for v in cdna_tests:
        normalized = norm.normalize_cdna(v)
        print(f"  {v:18} -> {normalized}")
    
    # Test baseline matching
    print("\n4. Baseline Matching with Fuzzy Position:")
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
