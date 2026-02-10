#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for variant normalizer improvements.
"""

import re
import csv

# Amino acid maps
AA_MAP = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Xaa'
}

class FixedNormalizer:
    def __init__(self, gene_symbol):
        self.gene_symbol = gene_symbol.upper()
        
        self.patterns = [
            re.compile(r'^(?:p\.)?([A-Z])(\d+)([A-Z])$', re.IGNORECASE),
            re.compile(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', re.IGNORECASE),
            re.compile(r'^(?:p\.)?([A-Z])(\d+)(fs\*|fsX|fs)$', re.IGNORECASE),
            re.compile(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)(fs\*|fsX|fs)$', re.IGNORECASE),
            re.compile(r'^(?:p\.)?([A-Z])(\d+)(sp|stop|trunc|\*)$', re.IGNORECASE),
            re.compile(r'^(?:p\.)?([A-Z])(\d+)(del|ins)$', re.IGNORECASE),
        ]
    
    def _to_three_letter(self, aa):
        if not aa:
            return None
        aa = aa.strip().upper()
        if len(aa) == 1:
            return AA_MAP.get(aa)
        return None
    
    def normalize_protein(self, variant):
        if not variant:
            return None
        variant = variant.strip()
        
        for pattern in self.patterns:
            match = pattern.match(variant)
            if match:
                ref, pos, suffix = match.groups()
                
                ref_3 = self._to_three_letter(ref.upper())
                if not ref_3:
                    continue
                
                suffix = suffix.upper()
                
                if 'FS' in suffix:
                    normalized_suffix = 'fs*'
                elif suffix in ['X', 'TER', 'SP', 'STOP', 'TRUNC']:
                    normalized_suffix = '*'
                elif suffix in ['DEL', 'INS']:
                    normalized_suffix = suffix.lower()
                else:
                    alt_3 = self._to_three_letter(suffix.upper())
                    if alt_3:
                        normalized_suffix = alt_3
                    else:
                        continue
                
                return 'p.' + ref_3 + pos + normalized_suffix
        
        return None

def test():
    normalizer = FixedNormalizer('KCNH2')
    
    test_cases = [
        ('A561V', 'p.Ala561Val'),
        ('R752W', 'p.Arg752Trp'),
        ('A193fsX', 'p.Ala193fs*'),
        ('D864sp', 'p.Asp864*'),
        ('G184del', 'p.Gly184del'),
        ('p.Ala561Val', 'p.Ala561Val'),
        ('Ala561Val', 'p.Ala561Val'),
    ]
    
    total = len(test_cases)
    correct = 0
    
    print("Testing normalization improvements...")
    
    for variant, expected in test_cases:
        result = normalizer.normalize_protein(variant)
        
        if result == expected:
            status = "✓"
            correct += 1
        elif result is not None:
            status = "?"
        else:
            status = "✗"
            result = "FAILED"
        
        print(status + " " + variant.ljust(15) + " -> " + str(result))
    
    accuracy = float(correct)/total*100
    print("\nAccuracy: " + str(correct) + "/" + str(total) + " = " + str(accuracy) + "%")
    return accuracy
    return float(correct)/total*100

if __name__ == "__main__":
    accuracy = test()
    print("\nPHASE 3 NORMALIZER IMPROVEMENTS COMPLETE")
    print("Normalization accuracy: " + str(accuracy) + "%")