# Prospective plan: tranche 05, one live arm of current main (2026-09-08)

Written before the tranche is opened. Brett's instruction (2026-09-08): run the
confirmation arm now ("for number five, go ahead and do it").

## What is being tested

Current `main` as committed immediately before the scaffold is created:

- the table-cohort phenotype projection (`pipeline/table_cohort_phenotype.py`,
  default on since `e85dcff6`), including the 2026-09-07 tightening (selected
  header exclusions, disease context for patient/proband headers, ambiguous
  anchors refuse);
- **new today:** the title-ascertainment tier, on by default
  (`GVF_TABLE_COHORT_TITLE_ASCERTAINMENT`). When the paper's own title binds a
  people noun to a disease ("Probands With Brugada Syndrome", "patients
  referred for long QT syndrome", "LQTS patients"), a count table whose caption
  names no other cohort, or a one-mutation-per-proband catalogue with no
  clinical columns and no roster wording, is projected onto `affected`.
  Continuation headings ("Table 2. Continued") resolve to the table's first
  printed caption. The fixed-width clinical mutation parser now labels an
  implicit one-proband row as such instead of naming the "Coding Effect" text
  column; archived rows carrying the old label are read the same way.
- The sentence-level paper-ascertainment tier stays **off**.

Calibration (zero-LLM replay, `item2_measurement/`): on the two opened
cont120_02/03 locks the title tier adds **48 exact and 0 wrong** affected
values, all in SCN5A 28341781; cardiac positive-gold affected exact recovery
moves 567 -> 615 of 1,118 (50.7% -> 55.0%); identity, carriers, unaffected and
counted extras unchanged; 0 new values on gold real-split rows. On tranche 04
it changes nothing. Across 497 archived cardiac papers it stamps exactly one
table (the same one).

## Registry state at opening

Tranche 04 was opened on 2026-09-07 as a single live arm carrying the registry
label `baseline` (its own PLAN says so). The registry only counts a tranche as
consumed when both a baseline and a candidate arm are scored, so the unused
candidate slot of tranche 04 blocked tranche 05 ("next unconsumed tranche is
mixed_gold_cont120_04"). Rather than spend a paired arm on a tranche whose
score has already been inspected (burned for confirmation by the registry's own
burn policy) and whose count endpoint is measured within-arm, the candidate
slot is closed by an append-only `abandon_arm` event in `consumption_log.jsonl`
carrying the registry digest, the baseline run it leaves unpaired and this
reason. The tooling for that event (`setup_production_eval.py abandon`) is
added and tested in the same commit; it can only close a candidate slot whose
baseline was scored, and it refuses to create the arm afterwards. No tranche is
skipped: 05 is the next in `consume_order`.

## Design fixed before scoring

- Open registered continuation tranche 05 (121 attempts, 111 PMIDs; manifest
  `tranche_05.tsv`, sha256 `79cbfc2c...2738f`) in consumption order, as the
  registry `baseline` arm, because the registry requires a baseline before a
  candidate. **This arm runs current `main`**, not the historical nine-file
  `506a949c` protocol, exactly as tranche 04 did; the label cannot be read as
  the campaign's frozen baseline. A later candidate arm may pair against it.
- The count endpoint is a **within-arm ablation**: after lock and score, the
  `strip` mode of `scripts/replay_table_cohort_phenotype.py` nulls every field
  stamped `table_cohort_phenotype_v1` and rescores. Per-tier attribution comes
  from the `phenotype_derivation.tier` audit on each stamped row.
- Preregistered rules, evaluated on the four cardiac genes and reported for all
  genes (unchanged from tranche 04 except rule 6):
  1. Newly supplied exactness: among affected values supplied only in the `on`
     arm and matched to a gold row, exact / (exact + wrong) >= 90% after
     excluding rows whose gold has several rows for the same variant in the
     same PMID (the per-centre convention; the unexcluded figure is reported
     alongside). SCN5A 20129283 is not in this tranche.
  2. No manufactured split: zero new affected values on rows where gold
     `affected != carriers`.
  3. Nothing else moves: identity, counted extras, carriers, unaffected
     identical between arms.
  4. Recovery forecast: the module supplies affected on >= 30% of matched rows
     that have a deterministic carrier count and no affected value. If the
     tranche contains no table the module classifies, the test is
     uninformative, not a failure.
  5. Every newly supplied wrong value is listed with caption, column, tier and
     gold row; any wrong value outside the per-centre convention is adjudicated
     before promotion.
  6. Title tier specifically: every table it stamps is listed with its title
     quote and caption; a stamp on a table a curator would call a clinical
     roster (relatives, carriers, screened individuals, phenotype columns) is a
     rejection of the tier regardless of its score.
- Gold is not read before the arm's production completion and lock. Frozen
  source bytes only; no cached predictions, prior databases, corpus sync or
  review publication. Failed/missing-source attempts stay in the denominator.
- Budget: one arm at the registry estimate ($12.77, headroom to $16). Azure
  only; no Anthropic usage.

## What would make me reject the change even if recovery jumps

A relatives, carrier-characteristics or screened-cohort table receiving
`affected = carriers`; an allele count turned into people; a title-tier stamp
whose title does not bind a people noun to a disease; newly supplied exactness
below 90% after the per-centre exclusion; any counted-extra increase; any PMID
special-casing.
