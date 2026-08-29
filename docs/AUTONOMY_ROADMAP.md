# Autonomy-at-Scale Roadmap

Last reconciled: 2026-08-12.

This document defines the target operating model and design boundaries. It is
not a second task list: executable follow-ups live only in
[`../TASKS.md`](../TASKS.md), current measurements live in
[`RECALL_STATUS.md`](RECALL_STATUS.md), and landed changes are recorded in
[`PROTOCOL_CHANGELOG.md`](PROTOCOL_CHANGELOG.md).

## Goal

Run GVF unattended across hundreds of genes and tens of thousands of papers,
with automated quality gates as the primary control. Human adjudication is an
exception path for calibrated ambiguity, not a required per-paper or per-run
step. The system must generalize beyond cardiac missense evidence to truncating,
structural, segregation, and case-control evidence without per-gene heuristics.

## Landed foundation

- `gvf-run` writes machine-readable run status and returns non-zero when a
  completeness stage fails.
- Run manifests record code/configuration provenance and resolved model routes.
- SQLite migration uses per-file rollback boundaries.
- `pipeline/trust_gate.py` assigns count facts to `trusted` or `quarantine`
  without rewriting the raw observation.
- The comparison scorer defaults to the trusted count projection while retaining
  raw rows for identity/recall audit.
- End-to-end count error distinguishes an asserted zero from an unobserved value.
- Deterministic extraction leaves phenotype partitions NULL unless the source
  asserts them; it does not manufacture unaffected zero by complement.
- Variant Browser adjudications are lead-gated, versioned, auditable, and
  excluded from the cardiac headline unless their scope permits otherwise.

## Target trust record

Every scientific observation should retain:

```text
raw value                 immutable source observation
source location/quote     evidence pointer
count role                patient | carrier | cohort total | control | population
trust tier                trusted | quarantine
trust reasons/version     reproducible rule decision
```

Soft quarantine is preferred to deletion or NULL mutation. Re-tiering after a
rule correction must not require reconstructing the source observation.

## Transferable decision rules

Prioritize rules that remain meaningful across gene classes:

| Rule class | Example |
| --- | --- |
| Internal arithmetic | asserted partitions cannot exceed an asserted carrier total |
| Scope consistency | a screened cohort total is not a per-variant carrier count |
| Population contamination | gnomAD allele counts are not clinical carriers |
| Evidence strength | a non-null count without a usable locator is lower trust |
| Within-paper consistency | a single extreme count needs stronger local evidence |

Avoid cardiac-shaped priors, phenotype-specific vocabularies, or fixed magnitude
ceilings as universal truth. Unknown-gene position validation must fail closed
for the trusted tier until a transcript/length authority is available.

## Fleet acceptance boundary

A release is not accepted merely because aggregate cardiac recall is unchanged.
The gate must report, by stratum:

- trusted-tier precision from an adjudicated sample;
- recall of known quarantine classes;
- end-to-end count error, including missed observations;
- assertion coverage and NULL denominators;
- quarantine rate, so a gate cannot appear accurate by hiding everything; and
- cardiac, BRCA, and at least one cold-gene result separately.

Variant Browser should receive quarantine diffs and calibration samples rather
than every paper. The canonical rollout tiers are defined in
`benchmarks/evaluation_tiers/` (its `README.md` and `registry.json` are the
authority for how many there are); larger cohorts run only after smaller gates
pass.

## Generalization boundary

BRCA readiness requires notation-by-class support for legacy/BIC indels, IVS
forms, and exon-level CNVs; registered transcript/length metadata; a separately
labeled silver standard; and a blind-spot instrument for variant-like tokens the
parser cannot normalize. Do not copy a large KCNH2-specific alias file for each
new gene.

Count recovery remains carrier-only and default off until measured. The parked
per-paper final check remains an optional source-grounded reviewer, not a default
trust substitute. The retired partial `unassessed_count` experiment is not part
of the schema; unobserved phenotype partitions remain NULL as specified in
`TASKS.md` and `EXTRACTION_CONTRACT.md`.
