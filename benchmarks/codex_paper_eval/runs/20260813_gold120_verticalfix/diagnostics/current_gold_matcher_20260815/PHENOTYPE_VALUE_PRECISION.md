# Gold-120 affected/unaffected exact-match (2026-08-16)

Same locked predictions as `PRECISION_AND_COST.md`. No new extraction.
Locked parent `report.json` not rewritten.

Agent briefing: `docs/AFFECTED_UNAFFECTED_PRECISION.md`.

## Metric

Exact / supplied on matched rows where both sides have a value.

| | affected | unaffected |
| --- | ---: | ---: |
| Baseline | 51/74 (68.9%), MAE 0.905, 23 errors | 37/56 (66.1%), MAE 0.964, 19 errors |
| After `phenotype_count_guard` | 40/55 (72.7%), MAE 0.600, 15 errors | 27/40 (67.5%), MAE 0.750, 13 errors |

The guard removes family/figure copies (`affected == carriers`, no split).
It does not invent replacements.

## Errors the guard removes

16029385 M124R 16/0 vs 9/7; 18808722 A490T/K897T 7/0 vs 2/5; 22338672
R744X 4/0 vs 1/3 and K897T affected=7 with no carriers; 29016797 S631A
3 vs 0; 24667783 A300T 2 vs 0; 28491547 una 0 vs 3; 25463374 una 0 vs 1;
25171853 P1048Sfs 3 vs 1.

## Errors that remain (need a better read, not a wider NULL)

See `docs/AFFECTED_UNAFFECTED_PRECISION.md` §4. Largest leftover MAE:
SCN5A 15671429 figure split 23/2 vs 7/15; KCNH2 25819988 symptoms 6/7 vs
9/4.
