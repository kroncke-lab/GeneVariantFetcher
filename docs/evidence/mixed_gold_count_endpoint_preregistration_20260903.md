# Preregistration: secondary count endpoint for the mixed-gold tranches

Written 2026-09-03 **before** the tranche 02 candidate arm
(`runs/20260903_protocol_mixed02_candidate`, extraction in flight) is locked
or scored, and before any tranche 03 arm exists. Committed so the ordering is
verifiable from git history against `consumption_log.jsonl`.

## Why a second endpoint, and why now

The registry's only scored endpoint is paper-derived identity recall (+1.0 pp
observed, PMID-cluster bootstrap lower bound ≥ −1 pp) with an identity
precision guard (lower bound ≥ −2 pp). Two facts established on 2026-09-03
make that bar unreachable for a reading-protocol change on these tranches:

- The gold-blind source-presence sweep puts 15.8% of runnable gold rows behind
  a hard acquisition ceiling and 28.7% behind the wide one
  (`gold_source_presence_sweep_20260903.md`). Inside a paired comparison the
  corpus is frozen, so no candidate can move those rows.
- Stage-by-stage root cause of the two opened baselines found **5 reachable
  rows of 87 misses** on tranche 01 and **≤ 5 of 35** on tranche 02: the
  +1.0 pp bar is 3 rows on tranche 02 and about 2.4 on tranche 01, i.e. the
  candidate must recover essentially every reachable row to pass.

The campaign's stated goal ("get variant recall up and get MAE better") names
counts explicitly, and the tranche 01 candidate moved them by an order of
magnitude (carriers supplied 48 → 125 of 242 rows, conditional MAE 0.81 →
0.10) while the identity rule saw +0.83 pp. That result **motivated adding a
count endpoint**; that is disclosed here rather than hidden. The thresholds
below are chosen from first principles, not from any tranche's numbers.

## What does not change

- The identity rule stays **primary**. Its verdicts on tranche 01
  (`reject_or_revise_candidate`) and on any later tranche stand as recorded;
  this endpoint never re-adjudicates them.
- Identity **non-inferiority is a hard guard** on the count endpoint: a count
  pass requires the same one-sided 95% PMID-cluster bootstrap lower bounds the
  identity rule uses (recall ≥ −1 pp, precision ≥ −2 pp). Counts may not
  improve by dropping identities.
- Confirmation semantics are unchanged: a count pass on one tranche advances
  only to the next unopened tranche, where the **identical** candidate runtime
  must pass the same rule.
- Nothing here loosens the identity rule or its thresholds.

## Endpoint definitions (computable from the locked `report.json` alone)

All quantities are over the paper-derived lane, per gene-paper attempt,
pairing candidate against baseline on the same tranche.

1. **End-to-end carrier error.** For every gold row with an asserted
   carrier value, the error is `|gold − predicted|`, where `predicted` is the
   matched prediction's carrier value if the row was matched and the value
   was supplied, and **0 otherwise** (identity miss or abstention both count
   as the full error). This is the repository's existing end-to-end
   definition (`cli/compare_variants.py` end-to-end MAE), which
   `docs/RECALL_STATUS.md` already requires to be quoted alongside the
   conditional MAE. Aggregate: mean over the tranche's asserted gold rows.
2. **Carrier coverage on matched rows.** Among gold rows matched to a
   prediction, the share whose prediction supplies a carrier value.
3. **Conditional carrier MAE** (diagnostic only, not a pass criterion, because
   abstention improves it).

Same three quantities for affected and unaffected are **reported**, not
scored: their gold uses an explicit-zero convention the pipeline deliberately
does not infer, so a pooled rule would reward guessing zeros.

## Paired decision rule (secondary, from tranche 02 candidate onward)

A candidate **passes the count endpoint** on a tranche when all hold:

- identity guard: recall LB ≥ −1 pp and precision LB ≥ −2 pp (as above);
- end-to-end carrier MAE: `candidate − baseline` observed **≤ −0.05
  carriers per gold row**, and its one-sided 95% PMID-cluster bootstrap
  **upper** bound **< 0**;
- carrier coverage on matched rows: observed `candidate − baseline ≥ 0` and
  its one-sided 95% lower bound ≥ −5 pp (coverage may not be traded away).

Bootstrap: the registry's existing resampler (10,000 PMID-cluster resamples,
seed 2026090301), applied to the per-attempt sums so the statistic is exact
for the whole tranche.

A pass is a **discovery**; acceptance requires the identical candidate to
pass on the next unopened tranche. A count pass with an identity fail of the
primary rule is reported as "count pass / identity not passed", never as a
pass of the campaign's primary endpoint.

## Application

- Tranche 02 candidate (v2, being extracted when this was written): scored
  under both rules; because the count endpoint was fixed after tranche 01 was
  inspected and before tranche 02's candidate was scored, a tranche 02 count
  pass is a **discovery**, not a confirmation.
- Tranche 03: identical v2 runtime against a fresh frozen baseline is the
  confirmation for whichever rule(s) tranche 02 passed.
- Implementation: `run_eval.py compare` gains a `secondary_count_endpoint`
  block computed from `papers[].matched_variants`, `papers[].count`, and the
  locked per-row values; the registry gains the rule text above under
  `evaluation_design.secondary_endpoints`. The code change lands only after
  the tranche 02 candidate is locked (the harness file is inside the runtime
  fingerprint) and before tranche 03's baseline is scaffolded.
