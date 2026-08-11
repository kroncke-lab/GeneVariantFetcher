# Collaborator-gold count-semantics projection — 50 papers

## Outcome

The active strategy-comparison and scored count cohort is now the fixed cardiac 48
plus only the two BRCA2 papers with lead-approved Variant Browser adjudications
by Nate: PMIDs 26833046 and 26848529. Six internally derived BRCA2 papers were
removed from active membership. Historical 56-paper artifacts remain unchanged.

No extraction was rerun. Filtering the locked A1 predictions and reference rows
to the active 50-paper cohort gives:

| Slice | Supplied carrier observations | Absolute error after repair | Carrier MAE |
| --- | ---: | ---: | ---: |
| Cardiac 48 | 326 | 16 | 0.0491 |
| BRCA2 collaborator 2 | 3 | 4 | 1.3333 |
| Combined active 50 | 329 | 20 | 0.0608 |

The combined carrier MAE changes from 0.9058 (298/329) under the legacy answer
key to 0.0608 (20/329) after the cardiac count-scope/scorer repair. Carrier
count recall is 329/1007 (32.67%). The combined number is dominated by the 326
cardiac observations and must not obscure the much weaker BRCA2 two-paper
stratum.

## Provenance and limitations

- Both active BRCA2 papers appear in the lead-approved Variant Browser gold
  snapshot with source reviewer `nate` and approving lead `kronckbm`.
- PMID 26833046 still has family-versus-individual carrier-scope uncertainty.
- PMID 26848529 is a reviewed-positive subset, not exhaustive paper-level gold.
- The six removed papers remain in the historical curated benchmark and its
  dated prediction artifacts, but are excluded from active review publishing
  and new strategy-comparison manifests.

Machine-readable arithmetic and source-lock digests are in `metrics.json`.
