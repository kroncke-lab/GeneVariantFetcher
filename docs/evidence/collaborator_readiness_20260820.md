# Collaborator extraction readiness audit

Structural QC is not an empirical precision estimate. Publication remains on hold until a source-grounded sample is adjudicated and the Variant Browser imports trusted fields/quarantines.

| Gene | Gate | Papers | Full text | Variants | Nameless | Species-scope links | VF matched | VF held | VF quarantined |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | FAIL | 50/50 | 50/50 | 5348 | 0 | 0 | 2186/5348 (40.9%) | 346 | 324 |
| BRCA2 | FAIL | 45/45 | 44/45 | 1917 | 0 | 0 | 958/1917 (50.0%) | 129 | 137 |
| BMPR2 | FAIL | 50/50 | 45/50 | 554 | 0 | 0 | 303/554 (54.7%) | 56 | 32 |

## Gate failures

- **BRCA1:** mandatory VariantFeatures enrichment/quarantine gate leaves 346 ambiguous live variant identities
- **BRCA2:** mandatory VariantFeatures enrichment/quarantine gate leaves 129 ambiguous live variant identities
- **BMPR2:** mandatory VariantFeatures enrichment/quarantine gate leaves 56 ambiguous live variant identities

## Publication decision

**HOLD.** Passing these invariants is necessary but cannot establish paper-level precision. Keep public publication disabled pending a source-grounded sample and a trusted-field-only Variant Browser import.

## 2026-08-21 staging follow-up

The trusted-field-only private staging imports are complete and database-level
invariants pass. Public publication remains HOLD pending human adjudication; see
`docs/evidence/collaborator_staging_verification_20260821.md`.
