#!/usr/bin/env python3
"""Test the updated normalizer with comprehensive aliases."""

import sys
from pathlib import Path

# Add the repo root to path to import local utils
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

from utils.variant_normalizer import (
    normalize_variant, 
    _lookup_alias, 
    _KCNH2_COMPREHENSIVE_ALIASES
)

print('=== ALIAS DICTIONARY LOADED ===')
print(f'Aliases loaded: {len(_KCNH2_COMPREHENSIVE_ALIASES)}')

print()
print('=== ALIAS LOOKUP TESTS ===')
test_cases = [
    'p.Ala561Val',
    'P.ALA561VAL',
    'ALA561VAL',
    'a561v',
    'A561V',
    'p.Gly262Alafs*98',
    'G262FSX',
    'p.Arg176Trp',
    'W1001X',
    'W1001*',
    'P72Q',
    'p.Pro72Gln',
    'PRO72GLN',
    'S660L',
    'p.Ser660Leu',
    'G628del',
    'GLY628DEL',
]

print()
print(f"{'Input':<25} {'Alias Lookup':<20} {'Normalized'}")
print("-" * 60)
for v in test_cases:
    lookup = _lookup_alias(v, 'KCNH2')
    normalized = normalize_variant(v, 'KCNH2')
    print(f'{v:<25} {str(lookup):<20} {normalized}')

# Now check if key missed variants can be normalized
print()
print('=== KEY MISSED VARIANTS ===')
missed_test = [
    ('p.Gly262Alafs*98', 'G262fsX'),
    ('p.Arg176Trp', 'R176W'),
    ('p.Val644Phe', 'V644F'),
    ('p.Glu58Lys', 'E58K'),
    ('p.Gly903Arg', 'G903R'),
    ('p.Ser660Leu', 'S660L'),
    ('p.Trp1001Ter', 'W1001X'),
    ('P72Q', 'P72Q'),  # Should already match
]

print(f"\n{'Extracted Form':<25} {'Expected Gold':<15} {'Normalized':<15} {'Match?'}")
print("-" * 70)
for extracted, expected in missed_test:
    normalized = normalize_variant(extracted, 'KCNH2')
    match = "✓" if normalized == expected else f"✗ (got {normalized})"
    print(f'{extracted:<25} {expected:<15} {normalized:<15} {match}')
