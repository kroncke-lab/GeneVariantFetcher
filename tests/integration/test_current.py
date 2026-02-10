#!/usr/bin/env python3
"""Test current normalizer functionality."""
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from utils.variant_normalizer import VariantNormalizer

# Test some key cases
normalizer = VariantNormalizer('KCNH2')

test_cases = [
    'A561V',
    'R752W', 
    'A193fsX',
    'D864sp',
    'G184Del',
    'p.Ala561Val',
    'Ala561Val'
]

print('Current normalizer results:')
for case in test_cases:
    result = normalizer.normalize_protein(case)
    print(f'{case} -> {result}')