#!/usr/bin/env python3
"""Run variant comparison."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher')

import subprocess
result = subprocess.run([
    '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/venv/bin/python',
    '-m', 'cli.compare_variants',
    '--excel', '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls',
    '--sqlite', '/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db',
    '--outdir', '/mnt/temp2/kronckbm/gvf_output/KCNH2_comparison_20260204',
    '--variant_match_mode', 'fuzzy',
    '--fuzzy_threshold', '0.85'
], cwd='/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher', capture_output=True, text=True)

print(result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr)
