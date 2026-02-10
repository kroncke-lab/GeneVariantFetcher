#!/usr/bin/env python3
"""
Improved Variant Normalizer for GeneVariantFetcher

Enhancements:
1. Better handling of all variant nomenclature forms
2. Non-KCNH2 variant detection (TP53/KRAS hotspots, position > 1159)
3. Fuzzy position matching for off-by-one errors
4. cDNA prefix normalization
"""

import re
import logging
from typing import Optional, Dict, Any, Tuple, Set, List

logger = logging.getLogger(__name__)

# Amino acid mappings
AA_MAP = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Xaa'
}

AA_MAP_REVERSE = {v: k for k, v in AA_MAP.items()}
AA_MAP_REVERSE['Ter'] = '*'
AA_MAP_REVERSE['Stop'] = '*'

# Protein lengths
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


class ImprovedVariantNormalizer:
    """Enhanced variant normalizer with multi-format matching."""
    
    def __init__(self, gene_symbol: str):
        self.gene_symbol = gene_symbol.upper()
        self.protein_length = PROTEIN_LENGTHS.get(self.gene_symbol)
    
    def normalize_to_single_letter(self, variant: str) -> Optional[str]:
        """
        Normalize protein variant to single-letter format: A561V
        
        Handles: A561V, p.Ala561Val, Ala561Val, p.A561V
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
    
    def normalize_to_three_letter(self, variant: str) -> Optional[str]:
        """
        Normalize to three-letter HGVS format: p.Ala561Val
        """
        single = self.normalize_to_single_letter(variant)
        if not single:
            return None
        
        # Parse single letter
        m = re.match(r'^([A-Z])(\d+)([A-Z\*X]|fsX)$', single)
        if m:
            ref = AA_MAP.get(m.group(1))
            pos = m.group(2)
            alt_raw = m.group(3)
            
            if alt_raw == 'fsX':
                alt = 'fs*'
            elif alt_raw in ['X', '*']:
                alt = '*'
            else:
                alt = AA_MAP.get(alt_raw)
            
            if ref and alt:
                return f"p.{ref}{pos}{alt}"
        
        return None
    
    def normalize_cdna(self, variant: str) -> Optional[str]:
        """
        Normalize cDNA variant, ensuring c. prefix.
        """
        if not variant:
            return None
        
        v = variant.strip()
        if v.startswith('c.'):
            return v
        
        # Check if it looks like cDNA (numbers and nucleotides)
        if re.match(r'^\d+[\+\-\d]*[ACGT_>delinsdup]', v):
            return f"c.{v}"
        
        return None
    
    def is_non_target_variant(self, variant: str) -> Tuple[bool, Optional[str]]:
        """
        Check if variant is likely from a non-target gene.
        
        Returns: (is_non_target, reason)
        """
        single = self.normalize_to_single_letter(variant)
        if not single:
            return False, None
        
        # Check hotspots
        # Remove fsX suffix for hotspot check
        check_var = re.sub(r'fsX$', '', single)
        if check_var in NON_TARGET_HOTSPOTS:
            return True, f"Known hotspot in TP53/KRAS/BRAF/PIK3CA"
        
        # Check position
        m = re.match(r'^[A-Z](\d+)', single)
        if m:
            pos = int(m.group(1))
            if self.protein_length and pos > self.protein_length:
                return True, f"Position {pos} exceeds {self.gene_symbol} length ({self.protein_length})"
        
        return False, None
    
    def extract_position(self, variant: str) -> Optional[int]:
        """Extract amino acid position from protein variant."""
        single = self.normalize_to_single_letter(variant)
        if single:
            m = re.match(r'^[A-Z](\d+)', single)
            if m:
                return int(m.group(1))
        return None
    
    def get_all_forms(self, variant: str) -> Dict[str, str]:
        """
        Get all normalized forms of a variant.
        
        Returns dict with keys: single, three, original
        """
        result = {'original': variant}
        
        single = self.normalize_to_single_letter(variant)
        if single:
            result['single'] = single
        
        three = self.normalize_to_three_letter(variant)
        if three:
            result['three'] = three
        
        cdna = self.normalize_cdna(variant)
        if cdna:
            result['cdna'] = cdna
        
        return result


def match_variants_to_baseline(
    extracted: List[str],
    baseline: Set[str],
    gene_symbol: str = 'KCNH2',
    fuzzy_position: bool = True
) -> Dict[str, Any]:
    """
    Match extracted variants to a baseline set with improved normalization.
    
    Args:
        extracted: List of extracted variant strings
        baseline: Set of baseline variant strings
        gene_symbol: Target gene symbol
        fuzzy_position: Allow ±1 position matching
    
    Returns:
        Dict with matches, unmatched, filtered, and stats
    """
    normalizer = ImprovedVariantNormalizer(gene_symbol)
    
    # Build baseline position index
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
    
    # Also add cDNA forms
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
        # Check if non-target gene
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
    norm = ImprovedVariantNormalizer('KCNH2')
    
    test_variants = [
        'A561V', 'p.Ala561Val', 'Ala561Val', 'p.A561V',
        'A193fsX', 'G584S', 'R248W', 'P2006A',
        '1558-1G>C', 'c.3152+1G>A'
    ]
    
    print("Test variant normalization:")
    for v in test_variants:
        forms = norm.get_all_forms(v)
        is_non, reason = norm.is_non_target_variant(v)
        print(f"  {v}:")
        print(f"    Forms: {forms}")
        if is_non:
            print(f"    ⚠️  Non-target: {reason}")
