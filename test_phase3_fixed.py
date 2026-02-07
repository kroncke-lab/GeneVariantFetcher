#!/usr/bin/env python3
"""
Standalone test for variant normalizer improvements.
"""

import re
import logging
from typing import Optional

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

class FixedVariantNormalizer:
    """Fixed version of the variant normalizer."""
    
    def __init__(self, gene_symbol: str):
        self.gene_symbol = gene_symbol.upper()
        
        # Patterns for variant parsing
        self.protein_patterns = [
            # Standard three-letter: p.Ala561Val or Ala561Val
            re.compile(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', re.IGNORECASE),
            # Single-letter: p.A561V or A561V
            re.compile(r'^(?:p\.)?([A-Z])(\d+)([A-Z])$', re.IGNORECASE),
            # Stop codon: p.A561* or A561*, p.Ter561*
            re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)([\*]|Ter)$', re.IGNORECASE),
            # Frame-shift: A193fsX, A193fs*, p.Ala193fsX
            re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2}|\*)(\d+)(fsX|fs\*|fs)$', re.IGNORECASE),
            # Insertion/deletion: G184Del, G184Ins
            re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(Del|Ins|del|ins)$', re.IGNORECASE),
            # Special termination: sp, stop, trunc suffixes
            re.compile(r'^(?:p\.)?([A-Z]|[A-Z][a-z]{2})(\d+)(sp|stop|trunc)$', re.IGNORECASE),
        ]
    
    def _to_three_letter(self, aa: str) -> Optional[str]:
        """Convert amino acid to three-letter code."""
        if not aa:
            return None
        
        aa = aa.strip().upper()
        if len(aa) == 1:
            return AA_MAP.get(aa)
        elif len(aa) == 3:
            return aa.capitalize()
        return None
    
    def normalize_protein(self, variant: str) -> Optional[str]:
        """Normalize protein variant to p.Xxx###Yyy format."""
        if not variant:
            return None
            
        variant = variant.strip()
        original_variant = variant
        
        for pattern in self.protein_patterns:
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
                normalized_suffix = None
                
                # Handle different variant types based on suffix
                if 'FS' in suffix_upper:
                    normalized_suffix = 'fs*'
                elif suffix_upper in ['X', 'TER']:
                    normalized_suffix = '*'
                elif suffix_upper in ['SP', 'STOP', 'TRUNC']:
                    normalized_suffix = '*'
                elif suffix_upper in ['DEL', 'INS']:
                    normalized_suffix = suffix_upper.lower()
                elif len(suffix) == 1 and suffix.isalpha():
                    alt_3 = self._to_three_letter(suffix.upper())
                    if alt_3:
                        normalized_suffix = alt_3
                    else:
                        continue
                elif len(suffix) == 3:
                    # Already three-letter
                    normalized_suffix = suffix.capitalize()
                else:
                    continue
                
                return f"p.{ref_3}{pos}{normalized_suffix}"
        
        return None

def test_improvements():
    """Test the normalization improvements."""
    normalizer = FixedVariantNormalizer('KCNH2')
    
    # Major test cases based on analysis
    test_cases = [
        # KCNH2 notation issues
        ('A561V', 'p.Ala561Val'),
        ('A561T', 'p.Ala561Thr'),
        ('R752W', 'p.Arg752Trp'),
        ('A193fsX', 'p.Ala193fs*'),
        ('A1017fsX', 'p.Ala1017fs*'),
        ('D864sp', 'p.Asp864*'),
        ('G184Del', 'p.Gly184del'),
        ('G189Ins', 'p.Gly189ins'),
        ('R100X', 'p.Arg100*'),
        ('R100stop', 'p.Arg100*'),
        
        # Different input formats
        ('p.A561V', 'p.Ala561Val'),
        ('Ala561Val', 'p.Ala561Val'),
        ('p.Ala561Val', 'p.Ala561Val'),
        ('A561*', 'p.Ala561*'),
        ('ALA561VAL', 'p.Ala561Val'),
        
        # Edge cases
        ('A1V', 'p.Ala1Val'),
        ('M1T', 'p.Met1Thr'),
    ]
    
    total = len(test_cases)
    successful = 0
    correct = 0
    
    print("Testing Variant Normalizer Improvements")
    print("=" * 60)
    
    for raw_variant, expected in test_cases:
        result = normalizer.normalize_protein(raw_variant)
        
        if result is not None:
            successful += 1
            if result == expected:
                correct += 1
                status = "✓"
            else:
                status = "⚠"
        else:
            status = "✗"
            result = "FAILED"
        
        print(f"{status} {raw_variant:15} -> {result:20} (expected: {expected})")
    
    print("\n" + "=" * 60)
    print(f"Normalization Rate: {successful}/{total} = {100*successful/total:.1f}%")
    print(f"Accuracy: {correct}/{total} = {100*correct/total:.1f}%")
    
    # Test against some real data samples
    print("\n" + "Testing actual Excel data samples...")
    
    real_samples = [
        'A561V', 'A561T', 'R100Q', 'R181Q', 'R752W', 'R854W', 'A1038E', 'A1058E', 'C566R',
        'G628S', 'G745A', 'G803Y', 'N588D', 'N629S', 'D609G', 'T613M', 'S818L', 'L552S'
    ]
    
    real_success = 0
    for variant in real_samples[:10]:  # Test first 10
        normalized = normalizer.normalize_protein(variant)
        if normalized:
            real_success += 1
            print(f"  {variant} -> {normalized}")
        else:
            print(f"  {variant} -> FAILED")
    
    print(f"\nReal data coverage: {real_success}/{len(real_samples[:10])} = {100*real_success/10:.1f}%")
    
    return {
        'total_test_cases': total,
        'normalized': successful,
        'correct': correct,
        'accuracy': 100*correct/total
    }

if __name__ == "__main__":
    results = test_improvements()
    print("\n" + "=" * 60)
    print("PHASE 3 IMPROVEMENT SUMMARY")
    print(f"Fixed variant notation coverage: {results['accuracy']:.1f}%")
    print("Addressed single-letter→three-letter and missing 'p.' prefix issues")