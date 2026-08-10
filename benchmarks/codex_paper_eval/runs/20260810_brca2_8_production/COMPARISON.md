# Production strategy on the BRCA2 arm — first non-cardiac baseline

Run 2026-08-10. The 8 BRCA2 gold papers (`brca2_8_papers_20260810.tsv`) through
the same calibrated replay as `20260726_fixed48_production`: full `gvf-run`
pipeline (`--pmid-file --no-source-recovery`, all layers, azure_ai/grok-4.3
Tier 3 + gpt-5.6-sol verification), DB projected by `db_to_predictions.py`,
locked before scoring. Scored against the curator-adjudicated
`benchmarks/curated_extraction_eval/gold_overrides/BRCA2_recall_input.csv`
(Azure lead-approved adjudications — NOT the manual cardiac gold standard;
report separately). Wall clock 34.7 min; all 8 papers resolved full text,
0 abstract-only stubs; 60 traced LLM calls (report under
`production_runs/BRCA2/20260810_133029/`, untracked).

## Overall, vs the cardiac production baseline

| run | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| BRCA2 production (all layers) | 89 | 400 | 21 | 18.2% | **80.9%** | 0.297 |
| BRCA2 production (paper-derived layers only) | 81 | 384 | 29 | 17.4% | 73.6% | 0.282 |
| cardiac production (all layers, 48 papers) | 789 | 985 | 212 | 44.5% | 78.8% | 0.569 |

Variant recall on usable-source BRCA2 papers lands at the cardiac production
level (80.9% vs 78.8%) — the missense-centric-grammar concern does not show up
as a recall collapse on this arm, though these 8 papers were curated to have
usable source, so acquisition failure modes are out of frame by construction.

## Precision is dominated by known subset-gold papers

| PMID | TP | FP | FN | note |
|---|---:|---:|---:|---|
| 26848529 | 3 | 118 | 0 | gold is 3 Azure-confirmed variants; subset follow-up OPEN |
| 26833046 | 4 | 80 | 0 | gold restricted to 4 founder variants; count-semantics follow-up OPEN; registry notes 77 wrong-variant FPs historically |
| 15365993 | 19 | 106 | 8 | mutation spectrum with polymorphism/VUS mix — extras are mostly identity rows without gold entries |
| 18489799 | 29 | 56 | 12 | largest gold table (41 rows) |
| others (4 papers) | 34 | 40 | 1 | 10398279, 21356067, 22655046, 25802882 — all at 0.94–1.00 recall |

The two subset-gold papers alone contribute 198 of the 400 FPs (49.5%). Until
the 26833046/26848529 curator follow-ups close, BRCA2 precision (and to a
lesser degree count MAE) is a property of the answer key's deliberate scope,
not of the pipeline. Recall and per-paper FN counts are meaningful now.

## Counts

| field | supplied / gold assertions | count recall | MAE |
|---|---:|---:|---:|
| carriers | 45 / 110 | 41% | 0.31 |
| affected | 37 / 110 | 34% | 1.38 |
| unaffected | 37 / 110 | 34% | 0.32 |

Same shape as cardiac production: when the pipeline commits a count it is
accurate (carriers MAE 0.31); the deficit is omission (counts supplied on only
~34–41% of gold assertions). Count metrics are identical between the
all-layers and paper-only views because linkage rows carry no counts.
