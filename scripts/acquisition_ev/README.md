# Acquisition Expected-Value (EV) prototype

Predict, from a paper's **abstract alone**, whether it is worth acquiring — i.e.
how many phenotyped carriers and distinct genetic variants a full extraction would
likely yield — so the un-downloaded tail can be ranked and manual acquisition
effort spent on the highest-payoff papers first.

This is the deterministic **v1 + its honest evaluation**. It establishes that
abstract-level signal ranks true full-paper yield well above chance. It is the
evidence gate before wiring an acquisition score into the worklist and the Variant
Browser paper list.

## What it does

`predict_yield.py`:

1. Loads the **true full-paper gold counts** (unique variants, carriers) per PMID
   from `gene_variant_fetcher_gold_standard/normalized/<GENE>_recall_input.csv`.
2. Fetches each paper's **abstract + title + metadata** from PubMed (cached under
   `.cache/`, gitignored). No full text, no LLM.
3. Computes deterministic features (regex signals mirroring
   `pipeline/extraction_priority.py`) and a **fixed, unfit** score with three
   exposed pieces:
   - `p_relevant` — probability the paper has phenotyped individuals *with* variant
     identities;
   - `est_carriers`, `est_variants` — volume floors;
   - `ev_score` — the combined ranking value.
4. Evaluates the score against the true counts: Spearman rank correlation and a
   **yield-capture curve** (acquire the top-K% by score → what fraction of all gold
   carriers/variants is captured), vs a random baseline and an oracle ceiling.

Because the score is a fixed formula with no fitting, a win here cannot be
overfitting.

```bash
python scripts/acquisition_ev/predict_yield.py --email you@example.com
# --genes KCNH2,SCN5A to subset;  --refresh to ignore the abstract cache
```

Artifacts: `eval_output/metrics.json`, `eval_output/per_paper_scores.csv`. Exit code
is 0 when the score beats the noise thresholds, 2 otherwise.

## Result on the cardiac ion-channel gold standards (n=1502)

| Set | Spearman(EV, carriers) | Spearman(EV, variants) | carriers capture@20% | variants capture@20% |
| --- | ---: | ---: | ---: | ---: |
| **Pooled** | **+0.43** | **+0.35** | **0.53** (random 0.20, oracle 0.85) | **0.51** (random 0.20, oracle 0.76) |
| KCNH2 | +0.30 | +0.43 | 0.39 | 0.41 |
| KCNQ1 | +0.39 | +0.23 | 0.53 | 0.53 |
| SCN5A | +0.46 | +0.31 | 0.61 | 0.53 |
| RYR2  | +0.58 | +0.55 | 0.48 | 0.52 |

All correlations are positive and significant (p ≤ 5e-5 per gene, p≈0 pooled).
Acquiring the top 20% of papers by score recovers roughly **half of all carriers
and variants** — a **~2.5–2.7× lift** over random, and it beats a naive
abstract-length baseline (~0.28 capture) by ~2×.

**Verdict: signal, not noise.** Abstract-level ranking is good enough to prioritize
the acquisition tail.

## What carries the signal (and what doesn't)

Raw per-signal Spearman vs true carriers (pooled):

- `carrier_mentions` (people-noun density) **+0.51** — strongest single signal
- `cohort_size` (largest integer next to a people/family noun, or `n=`) **+0.37**
- `original_data` (mutation/genotype/cohort vocabulary) **+0.34**
- `variant_mentions` (c./p./rs notation in the abstract) **−0.16** — *anti-correlated*

The negative `variant_mentions` correlation is the key finding: abstracts that spell
out specific variants tend to be **single-variant case reports / functional
studies** (low carrier count), while high-yield cohort papers say "we screened N
patients" and list the variants only in a **supplement the abstract never shows**.
This confirms the anticipated risk — abstracts under-represent exactly the
supplement-heavy papers that matter most — and it means the useful predictor is
**"large phenotyped cohort" (people/cohort signals), not "counts variants in the
abstract."** The `est_variants = variant_mentions` proxy is therefore weak; treat it
as a floor only.

## Known limitations (→ v2 levers)

- **Population/GWAS/biobank studies are over-ranked.** Papers with a huge N but few
  per-gene curated carriers (e.g. PMID 23465283, cohort 6500 → 0–3 gold carriers)
  float to the top on `cohort_size`. Next: penalize association/population/GWAS
  context, or down-weight papers mentioning many genes.
- **No-abstract papers can't be scored.** 116/1502 (7.7%) have no PubMed abstract
  and default to the bottom — but they hold only 1.7% of gold carriers, so the cost
  is small. A title/journal-only fallback would recover them.
- **Relevance (`p_relevant`) is not yet tested against negatives.** The gold set is
  all-positive, so this eval only validates *volume ranking*. Testing whether
  `p_relevant` separates data papers from reviews/irrelevant needs a negative set
  (`benchmarks/curated_extraction_eval/negative_cases/` or filter-dropped papers).
- **`est_carriers`/`est_variants` are ordinal, not calibrated.** To publish real
  expected counts (not just a rank), fit/calibrate them against these gold labels.

## Next steps once this is trusted

1. Add the population/GWAS penalty and re-measure (cheap, should lift precision@K).
2. Optional LLM abstract-scorer for a sharper `p_relevant` + count estimate on the
   promising tail (Tier1→Tier2 style two-stage).
3. Run the scorer over the **full PubMed candidate universe** (via `cli/discover.py`),
   not just gold PMIDs, to surface the ranked manual-acquisition list.
4. Surface `p_relevant` / `est_carriers` / `ev_score` per paper in Variant Browser
   (see `docs/VARIANT_BROWSER_INTEGRATION.md`) so the to-acquire list is sortable
   there.
