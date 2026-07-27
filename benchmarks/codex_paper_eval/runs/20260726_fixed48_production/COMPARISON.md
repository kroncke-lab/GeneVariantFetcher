# Production strategy vs single-model protocols — locked 48-paper cardiac set

Run 2026-07-26. Same `selection.json` as `20260724_fixed48_sol`, asserted identical
on all 48 `(gene, pmid, source_sha256)` triples. 12 papers each for KCNH2, KCNQ1,
RYR2, SCN5A. Scored with `run_eval.py score`, same matcher and same gold
(`gene_variant_fetcher_gold_standard/normalized/<GENE>_recall_input.csv`).

## Overall

| run | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| production (all layers) | 789 | 985 | 212 | 44.5% | **78.8%** | 0.569 |
| production (paper-derived layers only) | 670 | 588 | 331 | 53.3% | 66.9% | 0.593 |
| gpt-5.6-sol, single model | 724 | 102 | 277 | 87.7% | 72.3% | **0.793** |
| grok-4.3, single model | 599 | 61 | 402 | 90.8% | 59.8% | 0.721 |

Production finds the most gold (recall 78.8%, +6.5 pts over sol) and is by far the
noisiest (985 false positives vs sol's 102).

## The false positives are almost all count-free

Only **56 of 985** production FPs (5.7%) assert any count. The rest are variant
identities with no clinical claim attached — mostly ClinVar linkage and regex
sweep hits.

| source_layer | FP | FP carrying a count | TP from same layer |
|---|---:|---:|---:|
| regex_table | 429 | 28 | 381 |
| clinvar | 395 | 0 | 118 |
| regex_text | 81 | 0 | 25 |
| llm_table | 55 | 27 | 181 |
| figure | 13 | 0 | 65 |
| llm_text | 10 | 1 | 18 |
| pubtator | 2 | 0 | 1 |

Restricted to predictions that actually assert a count — the pipeline's real
product — the precision gap nearly closes:

| run | count-bearing TP | count-bearing FP | precision |
|---|---:|---:|---:|
| production | 330 | 56 | 85.5% |
| gpt-5.6-sol | 642 | 61 | 91.3% |
| grok-4.3 | 521 | 57 | 90.1% |

## Counts are where the single model wins

| run | carriers recall / MAE | affected recall / MAE | unaffected recall / MAE |
|---|---|---|---|
| production | 33% / 1.42 | 29% / 0.14 | 23% / 0.36 |
| gpt-5.6-sol | 64% / 0.48 | 33% / 0.69 | 17% / 0.85 |
| grok-4.3 | 52% / 0.61 | 21% / 0.71 | 7% / 0.99 |

Production emits roughly 2.4× as many variant identities as sol (1,774 vs 724)
but attaches a carrier count to half as many of the matched ones. Its affected-count
MAE is much better than sol's (0.14 vs 0.69) at similar recall, so when production
does commit a count it is accurate; the deficit is omission, not arithmetic.

Count metrics are identical between the all-layers and paper-only views because
ClinVar/PubTator linkage rows carry no counts at all.

## Per gene, recall / precision

| gene | production (all) | production (paper) | sol | grok |
|---|---|---|---|---|
| KCNH2 | 88.0% / 34.3% | 78.4% / 49.5% | 70.5% / 84.1% | 48.3% / 97.9% |
| KCNQ1 | 73.8% / 46.9% | 42.8% / 40.0% | 63.8% / 83.2% | 62.7% / 81.7% |
| RYR2 | 67.9% / 56.1% | 66.5% / 58.0% | 76.1% / 93.3% | 72.9% / 93.5% |
| SCN5A | 83.6% / 54.9% | 81.8% / 70.6% | 81.4% / 91.8% | 58.6% / 93.5% |

KCNQ1 recall drops 31 points when ClinVar/PubTator rows are removed — that gene's
production recall leans heavily on database linkage rather than on reading the paper.

## Cost and time

| run | wall clock | tokens |
|---|---|---|
| production | 7,139 s (1 h 59 m) — KCNH2 20.9m, KCNQ1 20.4m, RYR2 38.7m, SCN5A 39.0m | not aggregated by `gvf-run` |
| gpt-5.6-sol | 4,555 s | 881,425 (644k in / 237k out) |
| grok-4.3 | 1,807 s wall, 2,704 s summed | 945,426 (687k in / 258k out) |

Production spend is not comparable: `gvf-run` does not aggregate per-run token
usage, so its cost column is empty by construction rather than by measurement.

## Caveats — read before quoting these numbers

1. **Source parity does not fully hold.** 11 of 48 papers were read from different
   bytes than the sol run read, because supplement folding runs inside the extract
   step and is not suppressed by `--no-source-recovery`. Most are re-folds
   (KCNH2 24667783: 52,591 → 79,975 bytes). One is a known regression:
   **KCNQ1 17470695** — sol read the richer `_CLEANED.md` (37,690 B) via the
   harness's `choose_source` Pareto selection, production took `_FULL_CONTEXT.md`
   (25,597 B) because `usable_sources` still prefers it. PR #174 fixed this in the
   harness; production never got the fix.
2. **Production's recall is understated by two skipped papers.** The extraction
   circuit breaker discarded KCNQ1 25087618 and RYR2 19926015 outright
   (`Contains HTML/markup garbage (3 patterns)`) over 1.70% and 0.54% markup
   density respectively. Those papers hold 1 and 40 gold rows; all 41 count as FN.
   Sol read the same files and extracted from both.
3. **Duplicate notations were merged before scoring.** Production stores one variant
   several times across layers (e.g. `p.Leu552Ser c.1655T>C` from llm_text with
   counts plus `p.L552S` from pubtator without). 360 such rows were collapsed using
   the harness's own `matches()`, keeping counts from the paper-derived layer, so
   production is not charged a false positive for a variant it got right. Without
   this merge its FP count would be higher and meaninglessly so.
4. **KCNH2 11844290 is a shared zero,** not a production failure: its corpus source
   is 9,580 chars of abstract plus a 190-char supplement fold, the regex scanner
   found 0 candidates, and both production and sol returned 0 variants.
5. **`LOCK.json` was written by the conversion flow, not `run_eval.py lock`,**
   because `command_lock` requires exact per-paper token telemetry that `gvf-run`
   does not emit. The guarantee it encodes still holds: predictions were finalized
   before any gold value was read.

## Reproduce

```bash
benchmarks/codex_paper_eval/runs/20260726_fixed48_production/run_production.sh
.venv/bin/python benchmarks/codex_paper_eval/runs/20260726_fixed48_production/db_to_predictions.py \
  --out benchmarks/codex_paper_eval/runs/20260726_fixed48_production/predictions.json
.venv/bin/python benchmarks/codex_paper_eval/run_eval.py score \
  --run-dir benchmarks/codex_paper_eval/runs/20260726_fixed48_production
```
