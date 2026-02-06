#!/usr/bin/env python3
import json
from pathlib import Path

extractions_dir = Path('/mnt/temp2/kronckbm/gvf_output/KCNH2_clean')
total_variants = 0
unique_variants = set()
paper_count = 0

for f in extractions_dir.glob('*.json'):
    try:
        data = json.loads(f.read_text())
        variants = data.get('variants', [])
        paper_count += 1
        total_variants += len(variants)
        for v in variants:
            protein = v.get('protein_notation', '') or ''
            cdna = v.get('cdna_notation', '') or ''
            key = (protein, cdna)
            if protein or cdna:
                unique_variants.add(key)
    except Exception as e:
        pass

print(f'CLEAN (regex disabled):')
print(f'  Papers: {paper_count}')
print(f'  Total variants: {total_variants}')
print(f'  Unique variants: {len(unique_variants)}')
