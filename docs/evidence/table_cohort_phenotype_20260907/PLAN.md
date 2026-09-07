# Prospective plan: table-cohort phenotype projection on tranche 04 (2026-09-07)

Written before the tranche is opened. Brett's request: raise affected-count
exactness "way up", adversarially, with agy and grok consulted, and a $100
API budget for testing.

## What is being tested

`pipeline/table_cohort_phenotype.py` (new, default on) projects a
deterministic table row's per-variant people count onto `affected` when the
table's own caption or count-column header names a disease case series, and
onto `unaffected` (with a closed `affected = 0`) when the table names a
control cohort counted in people. It is a post-extraction step over
code-parsed rows only; it never touches identity, carriers, model-authored
rows, or rows that already carry a phenotype value. The paper-ascertainment
tier is implemented but OFF (`GVF_TABLE_COHORT_PAPER_ASCERTAINMENT`), following
both reviewers.

Calibration on the two opened locks (zero-LLM replay, `replay_tier1/`):
positive-gold affected exact recovery 203 → 567 of 1,118 (18.2% → 50.7%) on
the four cardiac genes; 366 newly supplied values exact, 46 wrong, all 46 in
SCN5A 20129283 where gold stores one row of 1 per testing centre while Table 4
prints the pooled count (the pipeline's carriers are already scored wrong on
those rows for the same reason); 0 new affected values on rows where gold
reports a real split; identity, counted extras, carriers and unaffected
unchanged.

## Design fixed before scoring

- Open registered continuation tranche 04 (120 attempts, 110 PMIDs) in
  consumption order. One paid arm, created as the registry `baseline` arm
  because the registry requires a baseline before a candidate. **This arm runs
  current `main` including the projection (default on).** It is not the
  campaign's historical nine-file `506a949c` protocol, and the tranche's
  registered paired identity decision is not being made here; a later candidate
  arm may still pair against it. Recorded here so the label cannot mislead.
- The count endpoint is a **within-arm ablation**: after lock and score, the
  `strip` mode of `scripts/replay_table_cohort_phenotype.py` nulls every field
  whose provenance source is `table_cohort_phenotype_v1` in the locked
  extraction output and rescores. The two arms differ only by this
  deterministic step, so provider run-to-run variance cannot explain the delta
  (the failure mode of the previous paired tranches).
- Preregistered rules, evaluated on the four cardiac genes and reported for all
  genes:
  1. **Exactness of the newly supplied set** (grok's rule): among affected
     values supplied only in the `on` arm and matched to a gold row, exact /
     (exact + wrong) ≥ 90% after excluding rows whose gold has several rows for
     the same variant in the same PMID (the per-centre convention above); the
     unexcluded figure is reported alongside.
  2. **No manufactured split**: zero new affected values on rows where gold
     `affected != carriers`.
  3. **Nothing else moves**: identity TP/FP/FN, counted-extra rows, carrier
     supply/exactness and unaffected supply/exactness identical between arms
     (by construction; verified by the scorer).
  4. **Recovery forecast**: positive-gold affected exact recovery rises by at
     least half the calibration effect relative to the addressable pool, i.e.
     the module supplies affected on ≥ 30% of matched rows that have a
     deterministic carrier count and no affected value (calibration: 62%). If
     tranche 04 contains no table the module classifies, the test is
     uninformative, not a failure, and is reported as such.
  5. Every newly supplied wrong value is listed with its caption, column and
     gold row; any wrong value outside the per-centre gold convention is
     adjudicated individually before promotion.
- Gold is not read before the arm's production completion and lock. Frozen
  source bytes only; no cached predictions, prior databases, corpus sync or
  review publication. Failed/missing-source attempts stay in the denominator.
- Budget: one arm at the registry estimate (about $12.5, headroom $15.7). No
  Anthropic usage. A second arm (tranche 05 confirmation) is spent only if rules
  1–4 pass on tranche 04.

## What would make me reject the change even if recovery jumps

A 28237968-style proband/relatives table receiving `affected = carriers`; an
allele count turned into people; newly supplied exactness below 90% after the
per-centre exclusion; any counted-extra increase; any PMID special-casing.
